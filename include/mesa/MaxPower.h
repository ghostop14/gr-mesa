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

#ifndef INCLUDED_MESA_MAXPOWER_H
#define INCLUDED_MESA_MAXPOWER_H

#include <gnuradio/sync_block.h>
#include <mesa/api.h>

namespace gr {
namespace mesa {

/*!
 * \brief <+description of block+>
 * \ingroup mesa
 *
 */
class MESA_API MaxPower : virtual public gr::sync_block {
public:
  typedef std::shared_ptr<MaxPower> sptr;

  /*!
   * \brief Return a shared_ptr to a new instance of mesa::MaxPower.
   *
   * To avoid accidental use of raw pointers, mesa::MaxPower's
   * constructor is in a private implementation
   * class. mesa::MaxPower::make is the public interface for
   * creating new instances.
   */
  static sptr make(double sampleRate, int fft_size, float squelchThreshold,
                   float framesToAvg, bool produceOut, float stateThreshold,
                   float holdUpSec);

  virtual float getSquelchThreshold() const = 0;
  virtual void setSquelchThreshold(float newValue) = 0;
  virtual float getStateThreshold() const = 0;
  virtual void setStateThreshold(float newValue) = 0;
  virtual float getHoldTime() const = 0;
  virtual void setHoldTime(float newValue) = 0;
};

} // namespace mesa
} // namespace gr

#endif /* INCLUDED_MESA_MAXPOWER_H */
