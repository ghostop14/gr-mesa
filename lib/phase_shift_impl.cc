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

#include "phase_shift_impl.h"
#include <gnuradio/io_signature.h>

#include <volk/volk.h>

namespace gr {
namespace mesa {

phase_shift::sptr phase_shift::make(float shift_in_radians)
{
    return gnuradio::make_block_sptr<phase_shift_impl>(shift_in_radians);
}

/*
 * The private constructor
 */
phase_shift_impl::phase_shift_impl(float shift_in_radians)
    : gr::sync_block("phase_shift",
                     gr::io_signature::make(1, 1, sizeof(gr_complex)),
                     gr::io_signature::make(1, 1, sizeof(gr_complex))) {
  d_shift_in_radians = shift_in_radians;
  d_shift_cc = gr_complex(cos(d_shift_in_radians), sin(d_shift_in_radians));

  message_port_register_in(pmt::mp("shift_rad"));
  set_msg_handler(pmt::mp("shift_rad"),
                  [this](pmt::pmt_t msg) { this->handle_msg_in(msg); });
}

/*
 * Our virtual destructor.
 */
phase_shift_impl::~phase_shift_impl() {}

void phase_shift_impl::handle_msg_in(pmt::pmt_t msg) {
  set_shift(pmt::to_float(msg));
}

float phase_shift_impl::get_shift() const { return d_shift_in_radians; }

void phase_shift_impl::set_shift(float newValue) {
  gr::thread::scoped_lock guard(d_setlock);
  d_shift_in_radians = newValue;
  d_shift_cc = gr_complex(cos(d_shift_in_radians), sin(d_shift_in_radians));
}

int phase_shift_impl::work(int noutput_items,
                           gr_vector_const_void_star &input_items,
                           gr_vector_void_star &output_items) {
  const gr_complex *in = (const gr_complex *)input_items[0];
  gr_complex *out = (gr_complex *)output_items[0];

  gr::thread::scoped_lock guard(d_setlock);

  if (d_shift_in_radians != 0.0) {
    volk_32fc_s32fc_multiply_32fc(out, in, d_shift_cc, noutput_items);
  } else {
    memcpy(out, in, sizeof(gr_complex) * noutput_items);
  }

  // Tell runtime system how many output items we produced.
  return noutput_items;
}

} /* namespace mesa */
} /* namespace gr */
