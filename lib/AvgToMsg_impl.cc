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
#include "AvgToMsg_impl.h"

#include <volk/volk.h>

namespace gr {
  namespace mesa {

  AvgToMsg::sptr
  AvgToMsg::make(int veclen)
  {
    return gnuradio::get_initial_sptr
      (new AvgToMsg_impl(veclen));
  }

  /*
   * The private constructor
   */
  AvgToMsg_impl::AvgToMsg_impl(int veclen)
    : gr::sync_block("AvgToMsg",
            gr::io_signature::make(1, 1, sizeof(float)*veclen),
            gr::io_signature::make(0, 0, 0))
  {
  	d_veclen = veclen;
      message_port_register_out(pmt::mp("avg"));
  }

  /*
   * Our virtual destructor.
   */
  AvgToMsg_impl::~AvgToMsg_impl()
  {
  }

  int
  AvgToMsg_impl::work(int noutput_items,
      gr_vector_const_void_star &input_items,
      gr_vector_void_star &output_items)
  {
    const float *in = (const float *) input_items[0];

    float sum=0.0;
    int vecstart;

    /*
    for (int nvecs=0;nvecs<noutput_items;nvecs++) {
  	  vecstart = nvecs*d_veclen;
        for (int i=0;i<d_veclen;i++) {
      	sum += in[vecstart + i];
        }
    }
	*/
    volk_32f_accumulator_s32f(&sum,in,(d_veclen*noutput_items));

    float avg = sum / (d_veclen*noutput_items);

    message_port_pub(pmt::mp("avg"),pmt::from_float(avg));

    // Tell runtime system how many output items we produced.
    return noutput_items;
  }

  } /* namespace mesa */
} /* namespace gr */

