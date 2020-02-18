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

#ifndef INCLUDED_MESA_MAXPOWER_IMPL_H
#define INCLUDED_MESA_MAXPOWER_IMPL_H

#include "signals_mesa.h"
#include <boost/circular_buffer.hpp>
#include <chrono>
#include <ctime>
#include <mesa/MaxPower.h>

using namespace MesaSignals;

namespace gr {
namespace mesa {

class MaxPower_impl : public MaxPower {
private:
  // Nothing to declare in this block.
  boost::mutex d_mutex;
  EnergyAnalyzer *pEnergyAnalyzer;
  double d_sampleRate;
  int d_framesToAvg;
  float d_squelchThreshold;

  int d_fftSize;

  bool d_produceOut;

  int iBufferCapacity;

  boost::circular_buffer<float> *maxBuffer;

  bool d_startInitialized;
  float d_holdUpSec;
  bool curState;
  float d_stateThreshold;
  std::chrono::time_point<std::chrono::steady_clock> holdTime;

  virtual void handleMsgIn(pmt::pmt_t msg);

  virtual int processData(int noutput_items, const gr_complex *in);
  virtual void sendState(bool state);

  virtual float calcAverage();

public:
  MaxPower_impl(double sampleRate, int fft_size, float squelchThreshold,
                float framesToAvg, bool produceOut, float stateThreshold,
                float holdUpSec);
  ~MaxPower_impl();

  void setup_rpc();
  virtual float getSquelchThreshold() const;
  virtual void setSquelchThreshold(float newValue);
  virtual float getStateThreshold() const;
  virtual void setStateThreshold(float newValue);
  virtual float getHoldTime() const;
  virtual void setHoldTime(float newValue);

  virtual bool stop();

  // Where all the action really happens
  int work(int noutput_items, gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
};

} // namespace mesa
} // namespace gr

#endif /* INCLUDED_MESA_MAXPOWER_IMPL_H */
