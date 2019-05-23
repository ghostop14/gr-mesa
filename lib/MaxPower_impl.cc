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
#include "MaxPower_impl.h"

namespace gr {
  namespace mesa {

    MaxPower::sptr
    MaxPower::make(float sampleRate, int fft_size, float squelchThreshold, float framesToAvg, bool produceOut, float stateThreshold, float holdUpSec)
    {
      return gnuradio::get_initial_sptr
        (new MaxPower_impl(sampleRate, fft_size, squelchThreshold, framesToAvg,produceOut,stateThreshold, holdUpSec));
    }

    /*
     * The private constructor
     */
    MaxPower_impl::MaxPower_impl(float sampleRate, int fft_size, float squelchThreshold, float framesToAvg, bool produceOut, float stateThreshold, float holdUpSec)
      : gr::sync_block("MaxPower",
              gr::io_signature::make(0, 1, sizeof(gr_complex)),
              gr::io_signature::make(0, 0, 0))
    {
    	d_startInitialized = false;
    	d_holdUpSec = holdUpSec;
    	curState = false;
    	d_stateThreshold = stateThreshold;

    	// std::cout << "[debug] << hold time: " << d_holdUpSec << " sec, threshold: " << stateThreshold << std::endl;

    	d_sampleRate = sampleRate;
    	d_framesToAvg = framesToAvg;
    	d_fftSize = fft_size;
    	d_produceOut = produceOut;
    	d_squelchThreshold = squelchThreshold; // This is also available in the energy analyzer, but this saves us a fn call in work

    	// Create energy analyzer
    	pEnergyAnalyzer = new EnergyAnalyzer(d_fftSize,squelchThreshold, 0.0);

    	// buffer capacity is for n seconds.  framestoavg * d_fftSize is the samples / block.  sample rate / that gets you
    	// blocks / sec.  Times seconds to avg gets you how many calls you need to average.
    	// float secondsToAvg = 0.2;
    	// float fBufferCapacity = secondsToAvg * d_sampleRate / (d_framesToAvg * d_fftSize);

    	// iBufferCapacity = (int)fBufferCapacity;
    	iBufferCapacity = 6;

    	if (iBufferCapacity == 0) {
    		iBufferCapacity = 1;
    	}

    	// std::cout << "Starting maxpower with a " << iBufferCapacity << " averaging buffer." << std::endl;

    	maxBuffer = new boost::circular_buffer<float>(iBufferCapacity);

		message_port_register_in(pmt::mp("msgin"));
        set_msg_handler(pmt::mp("msgin"), boost::bind(&MaxPower_impl::handleMsgIn, this, _1) );

        message_port_register_out(pmt::mp("out"));
        message_port_register_out(pmt::mp("maxpower"));
        message_port_register_out(pmt::mp("state"));

        gr::block::set_output_multiple(d_fftSize*d_framesToAvg);

    }

    /*
     * Our virtual destructor.
     */
    MaxPower_impl::~MaxPower_impl()
    {
    	bool retVal = stop();
    }

    bool MaxPower_impl::stop() {
    	if (pEnergyAnalyzer) {
    		delete pEnergyAnalyzer;
    		pEnergyAnalyzer = NULL;
    	}

    	if (maxBuffer) {
    		delete maxBuffer;
    		maxBuffer = NULL;
    	}

    	return true;
    }

    void MaxPower_impl::handleMsgIn(pmt::pmt_t msg) {
		pmt::pmt_t inputMetadata = pmt::car(msg);
		pmt::pmt_t data = pmt::cdr(msg);
		size_t noutput_items = pmt::length(data);
		const gr_complex *cc_samples;

		cc_samples = pmt::c32vector_elements(data,noutput_items);

		int retVal = processData(noutput_items,cc_samples);
    }

    int MaxPower_impl::processData(int noutput_items,const gr_complex *in) {
      	gr::thread::scoped_lock guard(d_mutex);

      FloatVector maxSpectrum;

      // last boolean param indicates to use the squelch for values below the configured squelch threshold.
      long samplesProcessed = pEnergyAnalyzer->maxHold(in, noutput_items, maxSpectrum, true);

      float maxPower = pEnergyAnalyzer->maxPower(maxSpectrum);

      maxBuffer->push_back(maxPower);
      float maxTotal = 0.0;

      // Wait till we have enough to do an average before we start sending data.
      if (maxBuffer->size() >= iBufferCapacity) {
          // Note: circular buffer size() returns the number of items in the buffer.
          // If the number of items < capacity, the number of items in the buffer are returned,
          // so the following loop should always be good.
          for (int i=0;i<maxBuffer->size();i++) {
        	  maxTotal += (*maxBuffer)[i];
          }

          float maxAvg = maxTotal / maxBuffer->size();

          // Send maxpower message
    	  pmt::pmt_t meta = pmt::make_dict();
    	  meta = pmt::dict_add(meta, pmt::mp("decisionvalue"), pmt::from_float(maxAvg));
    	  meta = pmt::dict_add(meta, pmt::mp("maxpower"), pmt::from_float(maxAvg));
    	  meta = pmt::dict_add(meta, pmt::mp("squelch"), pmt::from_float(d_squelchThreshold));

    	  pmt::pmt_t pdu = pmt::cons( meta, pmt::PMT_NIL );

    	  message_port_pub(pmt::mp("maxpower"),pdu);

    	  // Test our state conditions
    	  if (maxAvg >= d_stateThreshold) {
    		  // We're over our threshold.  Let's see if we need to notify.
    		  holdTime = std::chrono::steady_clock::now();
    		  d_startInitialized = true;
    		  // std::cout << "[Debug] Power above threshold" << std::endl;
    		  if (!curState) {
    			  // std::cout << "[Debug] Sending state true" << std::endl;
    			  sendState(true);
    			  curState = true;
    		  }
    	  }
    	  else {
    		  // We're below our threshold
    		  if (curState && d_startInitialized) {
    			  // We only need to worry about this if we were high.
    			  std::chrono::time_point<std::chrono::steady_clock> curTimestamp = std::chrono::steady_clock::now();
    			  // if (std::chrono::duration_cast<std::chrono::milliseconds>(curTimestamp-holdTime).count() > (double)d_holdUpSec) {
    			  if (std::chrono::duration<double>(curTimestamp-holdTime).count() > (double)d_holdUpSec) {
        			  // std::cout << "[Debug] Sending state false" << std::endl;
    				  sendState(false);
    				  curState = false;
    			  }
    		  }
    	  }

    	  // Send data message
    	  if (d_produceOut) {
    		  pmt::pmt_t data_out(pmt::init_c32vector(noutput_items, in));
    		  pmt::pmt_t datapdu = pmt::cons( meta, data_out );

    		  message_port_pub(pmt::mp("out"),datapdu);
    	  }
      }

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }


    void MaxPower_impl::sendState(bool state) {
        pmt::pmt_t meta = pmt::make_dict();

        if (state)
        	meta = pmt::dict_add(meta, pmt::mp("state"), pmt::mp(1));
        else
        	meta = pmt::dict_add(meta, pmt::mp("state"), pmt::mp(0));

        pmt::pmt_t pdu = pmt::cons( meta, pmt::PMT_NIL );

		message_port_pub(pmt::mp("state"),pdu);

    }

    int
    MaxPower_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
        const gr_complex *in = (const gr_complex *) input_items[0];

        return processData(noutput_items,in);
    }

  } /* namespace mesa */
} /* namespace gr */

