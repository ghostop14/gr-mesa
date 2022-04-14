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

#include "ioselector_impl.h"
#include <gnuradio/io_signature.h>
#include <functional>
namespace gr {
namespace mesa {

ioselector::sptr ioselector::make(
    int numinputs, int numoutputs, int inputport, int outputport, int itemsize)
{
    return gnuradio::make_block_sptr<ioselector_impl>(
        numinputs, numoutputs, inputport, outputport, itemsize);
}

/*
 * The private constructor
 */
ioselector_impl::ioselector_impl(int numinputs, int numoutputs, int inputport,
                                 int outputport, int itemsize)
    : gr::sync_block("ioselector",
                     gr::io_signature::make(1, numinputs, itemsize),
                     gr::io_signature::make(1, numoutputs, itemsize)) {
  d_numinputs = numinputs;
  d_numoutputs = numoutputs;
  d_curInput = inputport;
  d_curOutput = outputport;
  d_itemsize = itemsize;

  message_port_register_in(pmt::mp("inputindex"));
  set_msg_handler(pmt::mp("inputindex"),
                  std::bind(&ioselector_impl::handleMsgInputIndex, this, std::placeholders::_1));
  message_port_register_in(pmt::mp("outputindex"));
  set_msg_handler(
      pmt::mp("outputindex"),
      std::bind(&ioselector_impl::handleMsgOutputIndex, this, std::placeholders::_1));
}

/*
 * Our virtual destructor.
 */
ioselector_impl::~ioselector_impl() {}

void ioselector_impl::handleMsgInputIndex(pmt::pmt_t msg) {
  pmt::pmt_t data = pmt::cdr(msg);

  if (pmt::is_integer(data)) {
    int newPort = pmt::to_long(data);

    if (newPort < (d_numinputs)) {
      set_input_index(newPort);
    }
  }
}

void ioselector_impl::handleMsgOutputIndex(pmt::pmt_t msg) {
  pmt::pmt_t data = pmt::cdr(msg);

  if (pmt::is_integer(data)) {
    int newPort = pmt::to_long(data);
    if (newPort < (d_numoutputs))
      set_output_index(newPort);
  }
}

void ioselector_impl::set_input_index(int newValue) { d_curInput = newValue; }

void ioselector_impl::set_output_index(int newValue) { d_curOutput = newValue; }

int ioselector_impl::work(int noutput_items,
                          gr_vector_const_void_star &input_items,
                          gr_vector_void_star &output_items) {
  const char *in = (const char *)input_items[d_curInput];
  char *out = (char *)output_items[d_curOutput];

  int noi = noutput_items * d_itemsize;
  memcpy(out, in, noi);

  // Tell runtime system how many output items we produced.
  return noutput_items;
}

} /* namespace mesa */
} /* namespace gr */
