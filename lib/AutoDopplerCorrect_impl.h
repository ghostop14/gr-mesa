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

#ifndef INCLUDED_MESA_AUTODOPPLERCORRECT_IMPL_H
#define INCLUDED_MESA_AUTODOPPLERCORRECT_IMPL_H

#include "signals_mesa.h"
#include <chrono>
#include <ctime>
#include <gnuradio/fxpt_nco.h>
#include <mesa/AutoDopplerCorrect.h>

using namespace MesaSignals;

#define AUTODOPPLER_METHOD_CLOSESTSIGNAL 1
#define AUTODOPPLER_METHOD_BOXOUTSIDEIN 2

namespace gr {
namespace mesa {

class AutoDopplerCorrect_impl : public AutoDopplerCorrect {
protected:
  boost::mutex d_mutex;

  EnergyAnalyzer *pEnergyAnalyzer;
  int d_detectionMethod;

  gr::fxpt_nco d_nco;

  gr_complex *pMsgOutBuff;
  int msgBufferSize;

  double d_sampleRate;
  double d_centerFreq;
  double d_maxDrift;
  double d_expectedWidth;
  int d_shiftHolddownMS;
  bool d_processMessages;

  int d_framesToAvg;
  int d_fftSize;
  double d_minWidthHz;
  double d_maxWidthHz;
  bool d_startInitialized;
  float d_holdUpSec;

  double d_currentFreqShiftDelta;

  std::chrono::time_point<std::chrono::steady_clock> lastSeen, lastShifted;

  virtual void sendMessageData(gr_complex *data, long datasize,
                               double signalCenterFreq, double signalWidth,
                               float maxPower, pmt::pmt_t *pMetadata);
  void sendState(bool state);

public:
  AutoDopplerCorrect_impl(double freq, double sampleRate, double maxDrift,
                          double minWidth, double expectedWidth,
                          int shiftHolddownMS, int fft_size,
                          float squelchThreshold, int framesToAvg,
                          float holdUpSec, bool processMessages,
                          int detectionMethod);
  ~AutoDopplerCorrect_impl();

  virtual bool stop();

  // Needed to be public for debug testing
  virtual int processData(int noutput_items, const gr_complex *in,
                          gr_complex *out, pmt::pmt_t *pMetadata,
                          bool testMode = false);

  void handleMsgIn(pmt::pmt_t msg);

  // Where all the action really happens
  int work(int noutput_items, gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);

  virtual float getSquelch() const;
  virtual void setSquelch(float newValue);

  virtual double getCenterFrequency() const;
  virtual void setCenterFrequency(double newValue);

  virtual double getMinWidthHz() const;
  virtual void setMinWidthHz(double newValue);

  virtual double getExpectedWidth() const;
  virtual void setExpectedWidth(double newValue);

  virtual double getMaxDrift() const;
  virtual void setMaxDrift(double newValue);
};

} // namespace mesa
} // namespace gr

#endif /* INCLUDED_MESA_AUTODOPPLERCORRECT_IMPL_H */
