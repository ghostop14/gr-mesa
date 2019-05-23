/* -*- c++ -*- */
/* 
 * Copyright 2019 ghostop14.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "LongTermIntegrator_impl.h"
#include <volk/volk.h>

namespace gr {
  namespace mesa {

    LongTermIntegrator::sptr
    LongTermIntegrator::make(int fftsize, bool normalize)
    {
      return gnuradio::get_initial_sptr
        (new LongTermIntegrator_impl(fftsize,normalize));
    }

    /*
     * The private constructor
     */
    LongTermIntegrator_impl::LongTermIntegrator_impl(int fftsize, bool normalize)
      : gr::sync_block("LongTermIntegrator",
              gr::io_signature::make(1, 1, sizeof(float)*fftsize),
              gr::io_signature::make(1, 1, sizeof(float)*fftsize))
    {
    	d_fftsize = fftsize;
    	d_normalize = normalize;

		size_t memAlignment = volk_get_alignment();
    	aggBuffer = (float *)volk_malloc(fftsize*sizeof(float),memAlignment);
    	for (int i=0;i<d_fftsize;i++)
    		aggBuffer[i] = 0.0;

    	startTime = std::chrono::steady_clock::now();

		threadRunning = false;
		stopThread = false;
		readThread = new boost::thread(boost::bind(&LongTermIntegrator_impl::runThread, this));

        message_port_register_out(pmt::mp("runtime"));
    }

	void LongTermIntegrator_impl::runThread() {
		threadRunning = true;
		int loopCount = 0;
		int maxLoop = 30 / 0.01;  // 30 seconds

		usleep(10000); // let primary thread start.
    	startTime = std::chrono::steady_clock::now();

		while (!stopThread) {
			std::chrono::time_point<std::chrono::steady_clock> curTimestamp = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed_seconds = curTimestamp - startTime;
			float sec = elapsed_seconds.count();
			stringstream stream;
			string timeStr;
			if (sec < 60.0) {
				stream << std::fixed << std::setprecision(2) << sec;
				timeStr = stream.str() + " seconds";
			}
			else {
				stream << std::fixed << std::setprecision(2) << (sec/60.0);
				timeStr =stream.str() + " minutes";
			}

			// std::cout << "Sending " << timeStr << std::endl;

			message_port_pub(pmt::mp("runtime"),pmt::intern(timeStr));

			loopCount = 0;

			while (!stopThread && loopCount < maxLoop) {
				// check increment
				usleep(10000); // sleep 10 millisec
				loopCount++;
			}
		}

		threadRunning = false;
	}

    bool LongTermIntegrator_impl::stop() {
    	if (readThread) {
    		stopThread = true;

    		while (threadRunning) {
    			usleep(10000); // sleep 10 millisec
    		}

    		delete readThread;
    		readThread = NULL;
    	}

    	if (aggBuffer) {
    		volk_free(aggBuffer);
    		aggBuffer = NULL;
    	}

    	return true;
    }

    /*
     * Our virtual destructor.
     */
    LongTermIntegrator_impl::~LongTermIntegrator_impl()
    {
    	bool retval = stop();
    }

    int
    LongTermIntegrator_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
        const float *in = (const float *) input_items[0];
        float *out = (float *) output_items[0];
        int noi = noutput_items * d_fftsize;
        int baseIndex;
        float max=in[0];
        uint32_t maxIndex;


        // Vectors have to be done individually to map into aggBuffer;
        for (int curVector=0;curVector<noutput_items;curVector++) {
        	// aggBuffer[i] = aggBuffer[i] + in[curVector*d_fftsize + i];
    		volk_32f_x2_add_32f(aggBuffer,aggBuffer,&in[curVector*d_fftsize],d_fftsize);

    		// out[curVector*d_fftsize] = aggBuffer[i]
        	memcpy(&out[baseIndex],aggBuffer,d_fftsize*sizeof(float));

        	if (d_normalize) {
            	// find max
            	volk_32f_index_max_32u(&maxIndex,aggBuffer,d_fftsize);

            	if (aggBuffer[maxIndex] > max)
            		max = aggBuffer[maxIndex];
        	}
        }

        if (d_normalize) {
            // now normalize
            // out[i] = out[i] / max * 100.0;

            volk_32f_s32f_multiply_32f(out,out,100.0/max,noi);
        }

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace mesa */
} /* namespace gr */

