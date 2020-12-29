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

#include "AvgToMsg_impl.h"
#include <gnuradio/io_signature.h>

#include <volk/volk.h>

namespace gr {
namespace mesa {

AvgToMsg::sptr AvgToMsg::make(int veclen) {
    return gnuradio::make_block_sptr<AvgToMsg_impl>(veclen);
}

/*
 * The private constructor
 */
AvgToMsg_impl::AvgToMsg_impl(int veclen)
    : gr::sync_block("AvgToMsg",
                     gr::io_signature::make(1, 1, sizeof(float) * veclen),
                     gr::io_signature::make(0, 0, 0)) {
  d_veclen = veclen;
  b_hold = false;
  b_useStdDev = true;

  message_port_register_out(pmt::mp("avg"));
}

/*
 * Our virtual destructor.
 */
AvgToMsg_impl::~AvgToMsg_impl() {}

int AvgToMsg_impl::work(int noutput_items,
                        gr_vector_const_void_star &input_items,
                        gr_vector_void_star &output_items) {
  const float *in = (const float *)input_items[0];

  float sum = 0.0;
  int vecstart;
  long nitems = d_veclen * noutput_items;

  volk_32f_accumulator_s32f(&sum, in, nitems);

  float avg = sum / (d_veclen * noutput_items);

  /*
  if (b_useStdDev) {
      // include in result only if value is within 2 std_dev (98% of all samples
  should fall here) float standardDeviation = 0.0; for(long i = 0; i < nitems;
  ++i) { standardDeviation += pow(in[i] - avg, 2);
      }
      standardDeviation = sqrt(standardDeviation / (float)nitems);

      float TwoStdDev = 2.0*standardDeviation;

      long itemsUsed = 0;
      float newSum = 0.0;

      for(long i = 0; i < nitems; ++i) {
              if (fabs(in[i] - avg) <= TwoStdDev) {
                      itemsUsed++;
                      newSum += in[i];
              }
      }

      if (itemsUsed > 0) {
              avg = newSum / (float)itemsUsed;
      }

  }
      */

  if (!b_hold)
    message_port_pub(pmt::mp("avg"), pmt::from_float(avg));

  // Tell runtime system how many output items we produced.
  return noutput_items;
}

void AvgToMsg_impl::setHold(bool newValue) { b_hold = newValue; }

void AvgToMsg_impl::setup_rpc() {
#ifdef GR_CTRLPORT
  // Setters
  add_rpc_variable(rpcbasic_sptr(new rpcbasic_register_set<AvgToMsg_impl, bool>(
      alias(), "hold", &AvgToMsg_impl::setHold, pmt::mp(false), pmt::mp(true),
      pmt::mp(false), "bool", "hold", RPC_PRIVLVL_MIN,
      DISPTIME | DISPOPTSTRIP)));

#endif /* GR_CTRLPORT */
}
} /* namespace mesa */
} /* namespace gr */
