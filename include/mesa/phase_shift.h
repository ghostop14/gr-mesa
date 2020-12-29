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

#ifndef INCLUDED_MESA_PHASESHIFT_H
#define INCLUDED_MESA_PHASESHIFT_H

#include <gnuradio/sync_block.h>
#include <mesa/api.h>

namespace gr {
namespace mesa {

/*!
 * \brief This block will shift the incoming signal by the specified radians.
 * \ingroup mesa
 *
 */
class MESA_API phase_shift : virtual public gr::sync_block {
public:
  typedef std::shared_ptr<phase_shift> sptr;

  /*!
   * \brief Return a shared_ptr to a new instance of mesa::phase_shift.
   *
   * To avoid accidental use of raw pointers, mesa::phase_shift's
   * constructor is in a private implementation
   * class. mesa::phase_shift::make is the public interface for
   * creating new instances.
   */
  static sptr make(float shift_in_radians);
  virtual float get_shift() const = 0;
  virtual void set_shift(float newValue) = 0;
};

} // namespace mesa
} // namespace gr

#endif /* INCLUDED_MESA_PHASESHIFT_H */
