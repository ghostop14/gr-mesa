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

#ifndef INCLUDED_MESA_LONGTERMINTEGRATOR_IMPL_H
#define INCLUDED_MESA_LONGTERMINTEGRATOR_IMPL_H

#include "signals_mesa.h"
#include <mesa/LongTermIntegrator.h>

using namespace std;
using namespace MesaSignals;
#include <boost/thread/thread.hpp>
#include <chrono>
#include <ctime>

namespace gr {
namespace mesa {

class LongTermIntegrator_impl : public LongTermIntegrator {
private:
  int d_fftsize;
  bool d_normalize;
  double *aggBuffer;
  boost::mutex d_mutex;

  std::chrono::time_point<std::chrono::steady_clock> startTime;

  boost::thread *readThread = NULL;
  bool threadRunning;
  bool stopThread;

  void runThread();

public:
  LongTermIntegrator_impl(int fftsize, bool normalize);
  ~LongTermIntegrator_impl();

  virtual void reset(bool bReset);
  void setup_rpc();

  bool stop();

  // Where all the action really happens
  int work(int noutput_items, gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
};

} // namespace mesa
} // namespace gr

#endif /* INCLUDED_MESA_LONGTERMINTEGRATOR_IMPL_H */
