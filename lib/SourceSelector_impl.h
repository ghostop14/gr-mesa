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

#ifndef INCLUDED_MESA_SOURCESELECTOR_IMPL_H
#define INCLUDED_MESA_SOURCESELECTOR_IMPL_H

#include <chrono>
#include <ctime>
#include <mesa/SourceSelector.h>

#include <queue>

using namespace std;

namespace gr {
namespace mesa {
class SourceSelector_impl : public SourceSelector {
protected:
  boost::mutex d_mutex;
  boost::mutex d_queuemutex;
  float d_holdTime;
  int d_numInputs;
  int d_defaultInput;
  int d_inputBlockSize;
  int d_currentInput;

  // Data queue management
  queue<gr_complex> dataQueue;
  // long maxQueueSize;
  bool limitQueue;

  long minQueueLength;
  long initialDataQueueRequirement;
  bool initialQueueSizeMet;

  // Max Power for each input
  float maxPower[4];

  bool d_startInitialized;
  std::chrono::time_point<std::chrono::steady_clock> lastShifted;

  int maxPowerIndex();
  void queueData(pmt::pmt_t msg);

  int getDataAvailable();

  void sendNewPortMsg(int port);

  virtual void handleMsg(pmt::pmt_t msg, int port);

public:
  SourceSelector_impl(float holdTime, int numInputs, int defaultInput,
                      int inputBlockSize);
  ~SourceSelector_impl();
  virtual bool stop();

  void handleMsgIn1(pmt::pmt_t msg);
  void handleMsgIn2(pmt::pmt_t msg);
  void handleMsgIn3(pmt::pmt_t msg);
  void handleMsgIn4(pmt::pmt_t msg);

  // Where all the action really happens
  int work(int noutput_items, gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
};

} // namespace mesa
} // namespace gr

#endif /* INCLUDED_MESA_SOURCESELECTOR_IMPL_H */
