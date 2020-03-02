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

#ifndef INCLUDED_MESA_PHASE_SHIFT_IMPL_H
#define INCLUDED_MESA_PHASE_SHIFT_IMPL_H

#include <mesa/phase_shift.h>

namespace gr {
namespace mesa {

class phase_shift_impl : public phase_shift {
private:
  float d_shift_in_radians;
  gr_complex d_shift_cc;

public:
  phase_shift_impl(float shift_in_radians);
  ~phase_shift_impl();

  virtual float get_shift() const;
  virtual void set_shift(float newValue);

  void handle_msg_in(pmt::pmt_t msg);

  // Where all the action really happens
  int work(int noutput_items, gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
};

} // namespace mesa
} // namespace gr

#endif /* INCLUDED_MESA_PHASE_SHIFT_IMPL_H */
