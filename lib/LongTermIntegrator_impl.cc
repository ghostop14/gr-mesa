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

#include "LongTermIntegrator_impl.h"
#include <gnuradio/io_signature.h>
#include <volk/volk.h>
#include <iomanip>
namespace gr {
namespace mesa {

LongTermIntegrator::sptr LongTermIntegrator::make(int fftsize, bool normalize)
{
    return gnuradio::make_block_sptr<LongTermIntegrator_impl>(fftsize, normalize);
}

/*
 * The private constructor
 */
LongTermIntegrator_impl::LongTermIntegrator_impl(int fftsize, bool normalize)
: gr::sync_block("LongTermIntegrator",
		gr::io_signature::make(1, 1, sizeof(float) * fftsize),
		gr::io_signature::make(1, 1, sizeof(float) * fftsize)) {
	d_fftsize = fftsize;
	d_normalize = normalize;

	size_t memAlignment = volk_get_alignment();
	aggBuffer = (double *)volk_malloc(fftsize * sizeof(double), memAlignment);
	for (int i = 0; i < d_fftsize; i++)
		aggBuffer[i] = 0.0;

	startTime = std::chrono::steady_clock::now();

	threadRunning = false;
	stopThread = false;
	readThread =
			new boost::thread(boost::bind(&LongTermIntegrator_impl::runThread, this));

	message_port_register_out(pmt::mp("runtime"));
}

void LongTermIntegrator_impl::runThread() {
	threadRunning = true;
	int loopCount = 0;
	int maxLoop = 5 / 0.01; // 5 seconds

	usleep(10000); // let primary thread start.
	startTime = std::chrono::steady_clock::now();

	while (!stopThread) {
		std::chrono::time_point<std::chrono::steady_clock> curTimestamp =
				std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds;

		elapsed_seconds = curTimestamp - startTime;

		float sec = elapsed_seconds.count();
		stringstream stream;
		string timeStr;
		if (sec < 60.0) {
			stream << std::fixed << std::setprecision(2) << sec;
			timeStr = stream.str() + " seconds";
		} else {
			stream << std::fixed << std::setprecision(2) << (sec / 60.0);
			timeStr = stream.str() + " minutes";
		}

		// std::cout << "Sending " << timeStr << std::endl;

		message_port_pub(pmt::mp("runtime"), pmt::intern(timeStr));

		loopCount = 0;

		while (!stopThread && loopCount < maxLoop) {
			// check increment
			usleep(10000); // sleep 10 millisec
			loopCount++;
		}
	}

	threadRunning = false;
}

void LongTermIntegrator_impl::reset(bool bReset) {
	if (bReset) {
		gr::thread::scoped_lock guard(d_mutex);

		// Zero out all aggregation buckets
		for (int i = 0; i < d_fftsize; i++)
			aggBuffer[i] = 0.0;

		// Reset integration time
		startTime = std::chrono::steady_clock::now();
	}
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
LongTermIntegrator_impl::~LongTermIntegrator_impl() { bool retval = stop(); }

int LongTermIntegrator_impl::work(int noutput_items,
		gr_vector_const_void_star &input_items,
		gr_vector_void_star &output_items) {
	const float *in = (const float *)input_items[0];
	float *out = (float *)output_items[0];
	int noi = noutput_items * d_fftsize;
	uint32_t maxIndex;

	gr::thread::scoped_lock guard(d_mutex);

	// Vectors have to be done individually to map into aggBuffer;
	for (int curVector = 0; curVector < noutput_items; curVector++) {
		// Switched to double precision to avoid accumulation errors with floats.
		for (int i=0;i<d_fftsize;i++) {
			// aggBuffer[i] += in[curVector*d_fftsize + i];

			aggBuffer[i] += in[curVector*d_fftsize + i];
			out[curVector*d_fftsize+i] = aggBuffer[i];
		}
		// volk_32f_x2_add_32f(aggBuffer, aggBuffer, &in[curVector * d_fftsize], d_fftsize);
		// out[curVector*d_fftsize+i] = aggBuffer[i]
		// memcpy(&out[curVector*d_fftsize], aggBuffer, d_fftsize * sizeof(float));

		if (d_normalize) {
			// now normalize
			// find max in the agg buffer.
			// aggBuffer is double so can't use volk
			float max = aggBuffer[0];
			for (int i=1;i<d_fftsize;i++) {
				if (aggBuffer[i] > max) {
					max = aggBuffer[i];
				}
			}

			// Can use volk here since out is float
			if (max != 0.0)
				volk_32f_s32f_multiply_32f(&out[curVector*d_fftsize], &out[curVector*d_fftsize], 100.0 / max, d_fftsize);
		}
	}

	// Tell runtime system how many output items we produced.
	return noutput_items;
}

void LongTermIntegrator_impl::setup_rpc() {
#ifdef GR_CTRLPORT
	// Setters
	add_rpc_variable(
			rpcbasic_sptr(new rpcbasic_register_set<LongTermIntegrator_impl, bool>(
					alias(), "reset", &LongTermIntegrator_impl::reset, pmt::mp(false),
					pmt::mp(true), pmt::mp(false), "bool", "reset", RPC_PRIVLVL_MIN,
					DISPTIME | DISPOPTSTRIP)));

#endif /* GR_CTRLPORT */
}
} /* namespace mesa */
} /* namespace gr */
