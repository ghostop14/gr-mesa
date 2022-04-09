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

#include "SourceSelector_impl.h"
#include <gnuradio/io_signature.h>
#include <functional>
namespace gr {
namespace mesa {

SourceSelector::sptr
SourceSelector::make(float holdTime, int numInputs, int defaultInput, int inputBlockSize)
{
    return gnuradio::make_block_sptr<SourceSelector_impl>(
        holdTime, numInputs, defaultInput, inputBlockSize);
}

/*
 * The private constructor
 */
SourceSelector_impl::SourceSelector_impl(float holdTime, int numInputs,
                                         int defaultInput, int inputBlockSize)
    : gr::sync_block("SourceSelector", gr::io_signature::make(0, 0, 0),
                     gr::io_signature::make(0, 1, sizeof(gr_complex))) {
  d_holdTime = holdTime;
  d_numInputs = numInputs;
  d_defaultInput = defaultInput;
  d_inputBlockSize = inputBlockSize;

  if (defaultInput <= 0)
    defaultInput = 1;
  d_currentInput = defaultInput;

  d_startInitialized = false;

  limitQueue = false;

  // Initial anti-jitter buffer
  minQueueLength = d_inputBlockSize * 2;
  initialDataQueueRequirement = d_inputBlockSize * 6;
  initialQueueSizeMet = false;

  std::cout << "[Source Selector] Buffering initial frames..." << std::endl;

  for (int i = 0; i < 4; i++)
    maxPower[i] = -999.0;

  message_port_register_in(pmt::mp("in1"));
  set_msg_handler(pmt::mp("in1"),
                  std::bind(&SourceSelector_impl::handleMsgIn1, this, std::placeholders::_1));
  message_port_register_in(pmt::mp("in2"));
  set_msg_handler(pmt::mp("in2"),
                  std::bind(&SourceSelector_impl::handleMsgIn2, this, std::placeholders::_1));
  message_port_register_in(pmt::mp("in3"));
  set_msg_handler(pmt::mp("in3"),
                  std::bind(&SourceSelector_impl::handleMsgIn3, this, std::placeholders::_1));
  message_port_register_in(pmt::mp("in4"));
  set_msg_handler(pmt::mp("in4"),
                  std::bind(&SourceSelector_impl::handleMsgIn4, this, std::placeholders::_1));

  message_port_register_out(pmt::mp("inputport"));

  if (inputBlockSize > 0)
    gr::block::set_output_multiple(inputBlockSize);
}

bool SourceSelector_impl::stop() { return true; }

/*
 * Our virtual destructor.
 */
SourceSelector_impl::~SourceSelector_impl() { bool retVal = stop(); }

int SourceSelector_impl::maxPowerIndex() {
  int maxIndex = 0;
  float curMax = maxPower[0];

  for (int i = 1; i < 4; i++) {
    if (maxPower[i] > curMax) {
      maxIndex = i;
    }
  }

  return maxIndex;
}

int SourceSelector_impl::getDataAvailable() {
  gr::thread::scoped_lock guard(d_queuemutex);
  return dataQueue.size();
}

void SourceSelector_impl::queueData(pmt::pmt_t msg) {
  // Since we removed the output, just return.
  /*
  pmt::pmt_t data = pmt::cdr(msg);

  if (data == pmt::PMT_NIL)
          return;

  size_t vecSize = pmt::length(data);
  const gr_complex *cc_samples;
  cc_samples = pmt::c32vector_elements(data,vecSize);

  // queue the data
gr::thread::scoped_lock guard(d_queuemutex);
  for (long i=0;i<vecSize;i++)
          dataQueue.push(cc_samples[i]);
  */
}

void SourceSelector_impl::sendNewPortMsg(int port) {
  // send index out
  pmt::pmt_t pdu =
      pmt::cons(pmt::intern("inputport"), pmt::from_long(port - 1));
  message_port_pub(pmt::mp("inputport"), pdu);
}

void SourceSelector_impl::handleMsg(pmt::pmt_t msg, int port) {
  pmt::pmt_t meta = pmt::car(msg);

  // Take a look at max power to see what we want to do.
  float maxVal = pmt::to_float(
      pmt::dict_ref(meta, pmt::mp("decisionvalue"), pmt::mp(-999.0)));

  // Set the currently sent power as the power for the specified port
  maxPower[port - 1] = maxVal;

  // Get the port with the current max power
  int iMaxPowerPort = maxPowerIndex() + 1;

  if ((d_currentInput == port) && (iMaxPowerPort == port)) {
    // If the max power port is the current port, then we'll queue the data.
    queueData(msg);
  } else {
    // Need to see what didn't match.  If we're not max, we don't care and we
    // can just drop through The missing "else" to the if statement below.

    if (iMaxPowerPort == port) {
      // We're here because the max power port is not d_currentInput, it's this
      // port.  Which means the max power port has changed.

      // Need to check our hold-down timers
      if (!d_startInitialized) {
        // First time the port has changed.  Go ahead and move it as we may just
        // be initializing.
        d_startInitialized = true;
        d_currentInput = port;
        lastShifted =
            std::chrono::steady_clock::now(); // Initialize the shifted timer.

        // We haven't initialized prior to this, so this is locking on to the
        // first max power.  It may Hop a bit as the engine starts here.
        queueData(msg);
        sendNewPortMsg(port);
      } else {
        float powerDiff =
            fabs(maxPower[port - 1] - maxPower[d_currentInput - 1]);
        float pctDiff = fabs(powerDiff / maxPower[d_currentInput - 1]) *
                        100.0; // Convert to percent of maxPower

        // we're initialized so let's see if we're within our holddown period.
        std::chrono::time_point<std::chrono::steady_clock> curTimestamp =
            std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds =
            curTimestamp - lastShifted;
        if (elapsed_seconds.count() > (double)d_holdTime) {
          d_currentInput = port;
          lastShifted = curTimestamp; // Reset the shifted timer.
          queueData(msg);
          sendNewPortMsg(port);
        } // elapsed_seconds
          /*
           * The else to this that drops through is that we're not the max port
           * and we're within our hold-down   timer, so we're not allowed to shift.
           * Which means we have to just drop the data.  So there's no queueing.
           */
      }   // else initialized

    } // iMaxPower == port
  }
}

void SourceSelector_impl::handleMsgIn1(pmt::pmt_t msg) { handleMsg(msg, 1); }

void SourceSelector_impl::handleMsgIn2(pmt::pmt_t msg) { handleMsg(msg, 2); }

void SourceSelector_impl::handleMsgIn3(pmt::pmt_t msg) { handleMsg(msg, 3); }

void SourceSelector_impl::handleMsgIn4(pmt::pmt_t msg) { handleMsg(msg, 4); }

int SourceSelector_impl::work(int noutput_items,
                              gr_vector_const_void_star &input_items,
                              gr_vector_void_star &output_items) {
  int curQueueSize = getDataAvailable();

  gr_complex *out = (gr_complex *)output_items[0];

  if ((!initialQueueSizeMet && (curQueueSize < initialDataQueueRequirement)) ||
      (curQueueSize == 0)) {
    // std::cout << "iU";
    // Return zeros to keep the flowgraph running
    memset((void *)out, 0x00, noutput_items * sizeof(gr_complex));
    return noutput_items;
  }

  initialQueueSizeMet = true;

  int itemsProduced = noutput_items;

  if (curQueueSize < noutput_items)
    itemsProduced = curQueueSize;

  gr::thread::scoped_lock guard(d_queuemutex);
  for (long i = 0; i < itemsProduced; i++) {
    out[i] = dataQueue.front();
    dataQueue.pop();
  }

  return itemsProduced;
}

} /* namespace mesa */
} /* namespace gr */
