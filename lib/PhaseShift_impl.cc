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
#include "PhaseShift_impl.h"

#include <volk/volk.h>

namespace gr {
  namespace mesa {

  PhaseShift::sptr
  PhaseShift::make(float shift_in_radians)
  {
    return gnuradio::get_initial_sptr
      (new PhaseShift_impl(shift_in_radians));
  }

  /*
   * The private constructor
   */
  PhaseShift_impl::PhaseShift_impl(float shift_in_radians)
    : gr::sync_block("PhaseShift",
            gr::io_signature::make(1, 1, sizeof(gr_complex)),
            gr::io_signature::make(1, 1, sizeof(gr_complex)))
  {
  	d_shift_in_radians = shift_in_radians;

		message_port_register_in(pmt::mp("shift_rad"));
      set_msg_handler(pmt::mp("shift_rad"), boost::bind(&PhaseShift_impl::handleMsgIn, this, _1) );
  }

  /*
   * Our virtual destructor.
   */
  PhaseShift_impl::~PhaseShift_impl()
  {
  }

  void PhaseShift_impl::handleMsgIn(pmt::pmt_t msg) {
		setShift(pmt::to_float(msg));
  }

  float PhaseShift_impl::getShift() const {
  	return d_shift_in_radians;
  }

  void PhaseShift_impl::setShift(float newValue) {
      gr::thread::scoped_lock guard(d_mutex);
  	d_shift_in_radians = newValue;
  }

  int
  PhaseShift_impl::work(int noutput_items,
      gr_vector_const_void_star &input_items,
      gr_vector_void_star &output_items)
  {
    const gr_complex *in = (const gr_complex *) input_items[0];
    gr_complex *out = (gr_complex *) output_items[0];

    gr::thread::scoped_lock guard(d_mutex);

    if (d_shift_in_radians != 0.0) {
        gr_complex shift = gr_complex(cos(d_shift_in_radians),sin(d_shift_in_radians));
        /*
        for (int i=0;i<noutput_items;i++) {
      	  out[i] = in[i] * shift;
        }
        */
        volk_32fc_s32fc_multiply_32fc(out,in,shift,noutput_items);
    }
    else {
  	  memcpy(out,in,sizeof(gr_complex)*noutput_items);
    }

    // Tell runtime system how many output items we produced.
    return noutput_items;
  }

  void
	PhaseShift_impl::setup_rpc()
  {
#ifdef GR_CTRLPORT
  	// Getters
    add_rpc_variable(
      rpcbasic_sptr(new rpcbasic_register_get<PhaseShift_impl, float>(
	  alias(), "shift",
	  &AdvFileSink_impl::getShift,
    pmt::mp(0.0), pmt::mp(4.0), pmt::mp(0.0),
    "rad", "shift", RPC_PRIVLVL_MIN,
    DISPTIME | DISPOPTSTRIP)));

    // Setters
    add_rpc_variable(
      rpcbasic_sptr(new rpcbasic_register_set<PhaseShift_impl, float>(
	  alias(), "shift",
	  &AdvFileSink_impl::setShift,
    pmt::mp(0.0), pmt::mp(4.0), pmt::mp(0.0),
    "rad", "shift", RPC_PRIVLVL_MIN,
    DISPTIME | DISPOPTSTRIP)));

#endif /* GR_CTRLPORT */
  }

  } /* namespace mesa */
} /* namespace gr */

