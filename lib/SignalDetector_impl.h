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

#ifndef INCLUDED_MESA_SIGNALDETECTOR_IMPL_H
#define INCLUDED_MESA_SIGNALDETECTOR_IMPL_H

#include "signals_mesa.h"
#include <chrono>
#include <ctime>
#include <mesa/SignalDetector.h>

using namespace MesaSignals;

#define SIGDETECTOR_METHOD_SEPARATESIGNALS 1
#define SIGDETECTOR_METHOD_BOXOUTSIDEIN 2

namespace gr {
namespace mesa {

class SignalDetector_impl : public SignalDetector {
protected:
  boost::mutex d_mutex;
  EnergyAnalyzer *pEnergyAnalyzer;
  int d_detectionMethod;

  gr_complex *pMsgOutBuff;
  int msgBufferSize;

  double d_sampleRate;
  double d_centerFreq;
  double d_minWidthHz;
  double d_maxWidthHz;
  int d_framesToAvg;

  int d_fftSize;
  bool d_enableDebug;

  bool d_genSignalPDUs;

  std::chrono::time_point<std::chrono::steady_clock> startup, endup;
  bool d_startInitialized;
  float d_holdUpSec;

  // Methods
  float calcMinDutyCycle();
  virtual int processData(int noutput_items, const gr_complex *in,
                          gr_complex *out, pmt::pmt_t *pMetadata);
  void sendState(bool state);

public:
  SignalDetector_impl(int fftsize, float squelchThreshold, double minWidthHz,
                      double maxWidthHz, double radioCenterFreq,
                      double sampleRate, float holdUpSec, int framesToAvg,
                      bool genSignalPDUs, bool enableDebug,
                      int detectionMethod);
  virtual ~SignalDetector_impl();

  virtual bool stop();

  void handleMsgIn(pmt::pmt_t msg);

  virtual float getSquelch() const;
  virtual void setSquelch(float newValue);

  virtual double getCenterFrequency() const;
  virtual void setCenterFrequency(double newValue);

  virtual double getMinWidthHz() const;
  virtual void setMinWidthHz(double newValue);

  virtual double getMaxWidthHz() const;
  virtual void setMaxWidthHz(double newValue);

  // Where all the action really happens
  int work(int noutput_items, gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
};

} // namespace mesa
} // namespace gr

#endif /* INCLUDED_MESA_SIGNALDETECTOR_IMPL_H */
