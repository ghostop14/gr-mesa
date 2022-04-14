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

#include "AutoDopplerCorrect_impl.h"
#include <gnuradio/io_signature.h>
#include <volk/volk.h>
#include <functional>
#ifndef NDEBUG
#define PRINTDEBUG
//#else
// #define PRINTDEBUG
#endif

namespace gr {
namespace mesa {

AutoDopplerCorrect::sptr AutoDopplerCorrect::make(double freq,
                                                  double sampleRate,
                                                  double maxDrift,
                                                  double minWidth,
                                                  double expectedWidth,
                                                  int shiftHolddownMS,
                                                  int fft_size,
                                                  float squelchThreshold,
                                                  int framesToAvg,
                                                  float holdUpSec,
                                                  bool processMessages,
                                                  int detectionMethod)
{
    return gnuradio::make_block_sptr<AutoDopplerCorrect_impl>(freq,
                                                              sampleRate,
                                                              maxDrift,
                                                              minWidth,
                                                              expectedWidth,
                                                              shiftHolddownMS,
                                                              fft_size,
                                                              squelchThreshold,
                                                              framesToAvg,
                                                              holdUpSec,
                                                              processMessages,
                                                              detectionMethod);
}

/*
 * The private constructor
 */
AutoDopplerCorrect_impl::AutoDopplerCorrect_impl(
    double freq, double sampleRate, double maxDrift, double minWidth,
    double expectedWidth, int shiftHolddownMS, int fft_size,
    float squelchThreshold, int framesToAvg, float holdUpSec,
    bool processMessages, int detectionMethod)
    : gr::sync_block("AutoDopplerCorrect",
                     gr::io_signature::make(1, 1, sizeof(gr_complex)),
                     gr::io_signature::make(1, 1, sizeof(gr_complex))) {
  d_detectionMethod = detectionMethod;

  d_sampleRate = sampleRate;
  d_centerFreq = freq;

  d_maxDrift = maxDrift;
  d_expectedWidth = expectedWidth;
  d_shiftHolddownMS = shiftHolddownMS;
  d_processMessages = processMessages;

  d_framesToAvg = framesToAvg;
  d_fftSize = fft_size;

  d_startInitialized = false;
  d_holdUpSec = holdUpSec;

  pMsgOutBuff = NULL;
  msgBufferSize = 0;

  // Set up input filter
  // -------------------
  /*
  float filterFactor = 1.1;  // How much to assume filter is off.  1.1 = 10%
  float driftFactor = 1.5;  // How much the drift estimate may be off.

  // Decided to move this outside the block.  Integrating a FIR filter here with
  an FFT block
  // Was just easier to move it out.
  // d_detectionFilterWidth = (d_expectedWidth / 2.0 * filterFactor +
  d_maxDrift*driftFactor);
          */

  d_currentFreqShiftDelta = 0.0;
  d_nco.set_freq(0.0);

  // Set up energy detector
  // -------------------
  float hzPerBucket = d_sampleRate / d_fftSize;
  d_minWidthHz = minWidth;
  d_maxWidthHz = d_expectedWidth * 1.4;

  float binsForMinHz = d_minWidthHz / hzPerBucket;
  float minDutyCycle = binsForMinHz / d_fftSize;

  // Create energy analyzer
  pEnergyAnalyzer =
      new EnergyAnalyzer(d_fftSize, squelchThreshold, minDutyCycle);
  //    	std::cout << "min duty cycle: " << minDutyCycle << std::endl;

  // Make sure we have a multiple of fftsize coming in
  // -------------------
  gr::block::set_output_multiple(fft_size * d_framesToAvg);

  // Set up PDUs
  // -------------------
  message_port_register_in(pmt::mp("msgin"));
  set_msg_handler(pmt::mp("msgin"),
                  std::bind(&AutoDopplerCorrect_impl::handleMsgIn, this, std::placeholders::_1));

  message_port_register_out(pmt::mp("signaldetect"));
  message_port_register_out(pmt::mp("msgout"));
  message_port_register_out(pmt::mp("freq_info"));
  message_port_register_out(pmt::mp("freq_shift"));
  message_port_register_out(pmt::mp("state"));
}

bool AutoDopplerCorrect_impl::stop() {
  if (pEnergyAnalyzer) {
    delete pEnergyAnalyzer;
    pEnergyAnalyzer = NULL;
  }

  if (pMsgOutBuff) {
    volk_free(pMsgOutBuff);
    msgBufferSize = 0;
    pMsgOutBuff = NULL;
  }

  return true;
}

/*
 * Our virtual destructor.
 */
AutoDopplerCorrect_impl::~AutoDopplerCorrect_impl() { bool retVal = stop(); }

float AutoDopplerCorrect_impl::getSquelch() const {
  return pEnergyAnalyzer->getThreshold();
}

void AutoDopplerCorrect_impl::setSquelch(float newValue) {
  pEnergyAnalyzer->setThreshold(newValue);
}

double AutoDopplerCorrect_impl::getCenterFrequency() const {
  return d_centerFreq;
}

void AutoDopplerCorrect_impl::setCenterFrequency(double newValue) {
  // Since we moved the center frequency, if we're in a signal, let's send an
  // end message, then reset center freq and shift.
  gr::thread::scoped_lock guard(d_mutex);

  if (d_startInitialized) {
    pmt::pmt_t meta = pmt::make_dict();

    meta = pmt::dict_add(meta, pmt::mp("state"), pmt::mp(0));

    pmt::pmt_t pdu = pmt::cons(meta, pmt::PMT_NIL);
    message_port_pub(pmt::mp("signaldetect"), pdu);

    d_startInitialized = 0.0;
  }

  d_currentFreqShiftDelta = 0.0;

  d_centerFreq = newValue;
}

double AutoDopplerCorrect_impl::getMinWidthHz() const { return d_minWidthHz; }

void AutoDopplerCorrect_impl::setMinWidthHz(double newValue) {
  double hzPerBucket = d_sampleRate / d_fftSize;
  double binsForMinHz = d_minWidthHz / hzPerBucket;
  double minDutyCycle = binsForMinHz / d_fftSize;

  // Create energy analyzer
  pEnergyAnalyzer->setDutyCycle(minDutyCycle);

  d_minWidthHz = newValue;
}

double AutoDopplerCorrect_impl::getExpectedWidth() const {
  return d_expectedWidth;
}

void AutoDopplerCorrect_impl::setExpectedWidth(double newValue) {
  d_expectedWidth = newValue;
  d_maxWidthHz = d_expectedWidth * 1.4;
}

double AutoDopplerCorrect_impl::getMaxDrift() const { return d_maxDrift; }

void AutoDopplerCorrect_impl::setMaxDrift(double newValue) {
  d_maxDrift = newValue;
}

void AutoDopplerCorrect_impl::sendState(bool state) {
  int newState;
  if (state) {
    newState = 1;
  } else {
    newState = 0;
  }

  pmt::pmt_t pdu = pmt::cons(pmt::intern("state"), pmt::from_long(newState));

  message_port_pub(pmt::mp("state"), pdu);
}
void AutoDopplerCorrect_impl::handleMsgIn(pmt::pmt_t msg) {
  if (!d_processMessages)
    return;

  pmt::pmt_t inputMetadata = pmt::car(msg);
  pmt::pmt_t data = pmt::cdr(msg);
  size_t noutput_items = pmt::length(data);
  const gr_complex *cc_samples;

  cc_samples = pmt::c32vector_elements(data, noutput_items);

  if (noutput_items > msgBufferSize) {
    if (pMsgOutBuff)
      volk_free(pMsgOutBuff);

    size_t memAlignment = volk_get_alignment();
    pMsgOutBuff =
        (SComplex *)volk_malloc(noutput_items * sizeof(SComplex), memAlignment);
  }

  int result =
      processData(noutput_items, cc_samples, pMsgOutBuff, &inputMetadata);
}

void AutoDopplerCorrect_impl::sendMessageData(gr_complex *data, long datasize,
                                              double signalCenterFreq,
                                              double signalWidth,
                                              float maxPower,
                                              pmt::pmt_t *pMetadata) {
  if (!d_processMessages) {
    return;
  }

  pmt::pmt_t data_out(pmt::init_c32vector(datasize, data));

  if (!pMetadata) {
    pmt::pmt_t meta = pmt::make_dict();

    meta =
        pmt::dict_add(meta, pmt::mp("radioCenterFreq"), pmt::mp(d_centerFreq));
    meta = pmt::dict_add(meta, pmt::mp("sampleRate"), pmt::mp(d_sampleRate));
    meta = pmt::dict_add(meta, pmt::mp("signalCenterFreq"),
                         pmt::mp(signalCenterFreq));
    meta = pmt::dict_add(meta, pmt::mp("widthHz"), pmt::mp(signalWidth));
    meta = pmt::dict_add(meta, pmt::mp("maxPower"), pmt::mp(maxPower));

    pmt::pmt_t pdu = pmt::cons(meta, data_out);
    message_port_pub(pmt::mp("msgout"), pdu);
  } else {
    if (!pmt::dict_has_key(*pMetadata, pmt::mp("radioCenterFreq")))
      *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("radioCenterFreq"),
                                 pmt::mp(d_centerFreq));
    if (!pmt::dict_has_key(*pMetadata, pmt::mp("sampleRate")))
      *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("sampleRate"),
                                 pmt::mp(d_sampleRate));
    if (!pmt::dict_has_key(*pMetadata, pmt::mp("signalCenterFreq")))
      *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("signalCenterFreq"),
                                 pmt::mp(signalCenterFreq));
    if (!pmt::dict_has_key(*pMetadata, pmt::mp("widthHz")))
      *pMetadata =
          pmt::dict_add(*pMetadata, pmt::mp("widthHz"), pmt::mp(signalWidth));
    if (!pmt::dict_has_key(*pMetadata, pmt::mp("maxPower")))
      *pMetadata =
          pmt::dict_add(*pMetadata, pmt::mp("maxPower"), pmt::mp(maxPower));

    pmt::pmt_t pdu = pmt::cons(*pMetadata, data_out);
    message_port_pub(pmt::mp("msgout"), pdu);
  }
}

int AutoDopplerCorrect_impl::processData(int noutput_items,
                                         const gr_complex *in, gr_complex *out,
                                         pmt::pmt_t *pMetadata, bool testMode) {
  gr::thread::scoped_lock guard(d_mutex);

  FloatVector maxSpectrum;

  // last boolean param indicates to use the squelch for values below the
  // configured squelch threshold.
  long samplesProcessed =
      pEnergyAnalyzer->maxHold(in, noutput_items, maxSpectrum, true);

#ifdef PRINTDEBUG
  static int debugPrint = 300;
  if (debugPrint == 300) {
    printArray(maxSpectrum, "Max Hold");
    debugPrint = 0;
  }

  debugPrint++;
#endif

  // Now look if we have signals
  int numSignals = 0;
  SignalOverviewVector signalVector;

  // Last bool param says stop looking on the first detected signal or not.
  if (d_detectionMethod == AUTODOPPLER_METHOD_CLOSESTSIGNAL) {
    // Look for the closest signal
    numSignals = pEnergyAnalyzer->findSignals(
        (const float *)&maxSpectrum[0], d_sampleRate, d_centerFreq,
        d_minWidthHz, d_maxWidthHz, signalVector, false);
  } else {
    // This uses a boxing method, outside-in looking for a signal.
    // If you have a channelized signal, this approach will work better.

    SignalOverview signalOverview;
    numSignals = pEnergyAnalyzer->findSingleSignal(
        (const float *)&maxSpectrum[0], d_sampleRate, d_centerFreq,
        d_minWidthHz, signalOverview);

    if (numSignals > 0) {
      signalVector.push_back(signalOverview);
    }
  }

  // start initialized tracks if we've picked up a signal and we're in a "high"
  // / have signal state
  bool justDetectedSignal = false; // first detection
  bool lostSignal = false;
  bool signalPresent = false;
  bool inHoldDown = false;
  int closestIndex = 0;
  int closestDelta = 1.0e6;
  double curDelta;

  // Before processing, need to see if the detected signal is within our
  // expected range or a neighboring signal
  if (numSignals > 0) {
    for (int i = 0; i < signalVector.size(); i++) {
      curDelta = fabs(signalVector[i].centerFreqHz - d_centerFreq);
      if (curDelta < closestDelta) {
        closestIndex = i;
      }
    }

    curDelta = fabs(signalVector[closestIndex].centerFreqHz - d_centerFreq);
    if (curDelta > d_maxDrift) {
      // Still not the right signal.  Zero out signal counter.
      numSignals = 0;
    }
  }

  bool foundGoodSignal = false;

  for (int i = 0; i < signalVector.size(); i++) {
    double curDelta = fabs(signalVector[i].centerFreqHz - d_centerFreq);
    if (curDelta <= d_maxDrift)
      foundGoodSignal = true;

    if (curDelta < closestDelta) {
      closestIndex = i;
    }
  }

  if (!foundGoodSignal)
    numSignals = 0;

  if (numSignals > 0) { // We have a detection
    signalPresent = true;

    if (!d_startInitialized) {
      // Haven't seen a signal in a while, start the rising edge.
      lastSeen = std::chrono::steady_clock::now();
      lastShifted = lastSeen;
      d_startInitialized = true;
      justDetectedSignal = true;

      // set initial offset
      d_currentFreqShiftDelta =
          d_centerFreq - signalVector[closestIndex].centerFreqHz;

      d_nco.set_freq(2 * M_PI * d_currentFreqShiftDelta / d_sampleRate);

      // send msg notification
      if (!pMetadata) {
        pmt::pmt_t meta = pmt::make_dict();

        // NOTE: freq matches the key looked for by the signal source block if
        // needed.
        meta = pmt::dict_add(meta, pmt::mp("numsignals"),
                             pmt::mp((int)signalVector.size()));
        meta = pmt::dict_add(meta, pmt::mp("decisionvalue"),
                             pmt::mp((int)signalVector.size()));
        meta = pmt::dict_add(meta, pmt::mp("closestsignalnum"),
                             pmt::mp(closestIndex + 1));
        meta = pmt::dict_add(meta, pmt::mp("freq"),
                             pmt::mp(d_currentFreqShiftDelta));
        meta = pmt::dict_add(meta, pmt::mp("freqoffset"),
                             pmt::mp(d_currentFreqShiftDelta));
        meta = pmt::dict_add(meta, pmt::mp("signalcenterfreq"),
                             pmt::mp(signalVector[closestIndex].centerFreqHz));
        meta = pmt::dict_add(meta, pmt::mp("trackingcenterfreq"),
                             pmt::mp(signalVector[closestIndex].centerFreqHz));
        meta = pmt::dict_add(meta, pmt::mp("widthHz"),
                             pmt::mp(signalVector[closestIndex].widthHz));
        meta = pmt::dict_add(meta, pmt::mp("signalpower"),
                             pmt::mp(signalVector[closestIndex].maxPower));
        pmt::pmt_t pdu = pmt::cons(meta, pmt::PMT_NIL);

        message_port_pub(pmt::mp("freq_shift"),
                         pmt::from_double(d_currentFreqShiftDelta));

        if (!testMode)
          message_port_pub(pmt::mp("freq_info"), pdu);
      } else {
        *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("numsignals"),
                                   pmt::mp((int)signalVector.size()));
        *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("decisionvalue"),
                                   pmt::mp((int)signalVector.size()));
        *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("closestsignalnum"),
                                   pmt::mp(closestIndex + 1));
        *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("widthHz"),
                                   pmt::mp(signalVector[closestIndex].widthHz));
        *pMetadata =
            pmt::dict_add(*pMetadata, pmt::mp("signalpower"),
                          pmt::mp(signalVector[closestIndex].maxPower));
        if (!pmt::dict_has_key(*pMetadata, pmt::mp("freq")))
          *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("freq"),
                                     pmt::mp(d_currentFreqShiftDelta));
        if (!pmt::dict_has_key(*pMetadata, pmt::mp("freqoffset")))
          *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("freqoffset"),
                                     pmt::mp(d_currentFreqShiftDelta));
        if (!pmt::dict_has_key(*pMetadata, pmt::mp("signalcenterfreq")))
          *pMetadata =
              pmt::dict_add(*pMetadata, pmt::mp("signalcenterfreq"),
                            pmt::mp(signalVector[closestIndex].centerFreqHz));
        if (!pmt::dict_has_key(*pMetadata, pmt::mp("trackingcenterfreq")))
          *pMetadata =
              pmt::dict_add(*pMetadata, pmt::mp("trackingcenterfreq"),
                            pmt::mp(signalVector[closestIndex].centerFreqHz));

        message_port_pub(pmt::mp("freq_shift"),
                         pmt::from_double(d_currentFreqShiftDelta));

        pmt::pmt_t pdu = pmt::cons(*pMetadata, pmt::PMT_NIL);
        if (!testMode)
          message_port_pub(pmt::mp("freq_info"), pdu);
      }
    } else {
      // We're continuing to see a signal.  Move the end indicator
      lastSeen = std::chrono::steady_clock::now();
    }
  } else {                    // No Detection
    if (d_startInitialized) { // We had a signal so we can track losing it.
      // Before we say we've lost it, let's see if we're within our hold timer
      std::chrono::time_point<std::chrono::steady_clock> curTimestamp =
          std::chrono::steady_clock::now();
      std::chrono::duration<double> elapsed_seconds = curTimestamp - lastSeen;
      if (elapsed_seconds.count() > (double)d_holdUpSec) {
        // No detection and we've exceeded our hold window.  Reset start
        // tracker.
        d_startInitialized = false;
        lostSignal = true;

        d_currentFreqShiftDelta = 0.0;
        d_nco.set_freq(0.0);

        message_port_pub(pmt::mp("freq_shift"),
                         pmt::from_double(d_currentFreqShiftDelta));
      } else {
        inHoldDown = true;
      }
    }
  }

  // ----------------  Signal Processing Here ----------------------------------
  // For now, just copy in to out.  This is where we'll have to do the shift and
  // calcs eventually.
  /*
   * Steps:
   * 1. Determine if we're in a hold-down period
   * 2. If we're out of a hold-down period, Look at detected signal center
   * frequency and adjust
   * 3. Multiply to offset the signal into outbuff. (requires a generated sine
   * wave)
   */

  if (signalPresent) {
    // We have a signal, lets see if we're allowed to update our shift
    std::chrono::time_point<std::chrono::steady_clock> curTimestamp =
        std::chrono::steady_clock::now();

    if (std::chrono::duration_cast<std::chrono::milliseconds>(curTimestamp -
                                                              lastShifted)
            .count() > (double)d_shiftHolddownMS) {
      // We can update our shift
      lastShifted = curTimestamp;

      double tmpShift = d_centerFreq - signalVector[0].centerFreqHz;

      if (tmpShift != d_currentFreqShiftDelta &&
          (fabs(tmpShift) <= d_maxDrift)) {
        d_currentFreqShiftDelta =
            d_centerFreq - signalVector[closestIndex].centerFreqHz;
        d_nco.set_freq(2 * M_PI * d_currentFreqShiftDelta / d_sampleRate);

        if (!pMetadata) {
          pmt::pmt_t meta = pmt::make_dict();

          // NOTE: freq matches the key looked for by the signal source block if
          // needed.
          meta = pmt::dict_add(meta, pmt::mp("numsignals"),
                               pmt::mp((int)signalVector.size()));
          meta = pmt::dict_add(meta, pmt::mp("decisionvalue"),
                               pmt::mp((int)signalVector.size()));
          meta = pmt::dict_add(meta, pmt::mp("closestsignal"),
                               pmt::mp(closestIndex + 1));
          meta = pmt::dict_add(meta, pmt::mp("freq"),
                               pmt::mp(d_currentFreqShiftDelta));
          meta = pmt::dict_add(meta, pmt::mp("freqoffset"),
                               pmt::mp(d_currentFreqShiftDelta));
          meta =
              pmt::dict_add(meta, pmt::mp("signalcenterfreq"),
                            pmt::mp(signalVector[closestIndex].centerFreqHz));
          meta =
              pmt::dict_add(meta, pmt::mp("trackingcenterfreq"),
                            pmt::mp(signalVector[closestIndex].centerFreqHz));
          meta = pmt::dict_add(meta, pmt::mp("widthHz"),
                               pmt::mp(signalVector[closestIndex].widthHz));
          meta = pmt::dict_add(meta, pmt::mp("signalpower"),
                               pmt::mp(signalVector[closestIndex].maxPower));

          message_port_pub(pmt::mp("freq_shift"),
                           pmt::from_double(d_currentFreqShiftDelta));

          pmt::pmt_t pdu = pmt::cons(meta, pmt::PMT_NIL);
          if (!testMode)
            message_port_pub(pmt::mp("freq_info"), pdu);
        } else {
          *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("numsignals"),
                                     pmt::mp((int)signalVector.size()));
          *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("decisionvalue"),
                                     pmt::mp((int)signalVector.size()));
          *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("closestsignal"),
                                     pmt::mp(closestIndex + 1));
          *pMetadata =
              pmt::dict_add(*pMetadata, pmt::mp("widthHz"),
                            pmt::mp(signalVector[closestIndex].widthHz));
          *pMetadata =
              pmt::dict_add(*pMetadata, pmt::mp("signalpower"),
                            pmt::mp(signalVector[closestIndex].maxPower));
          if (!pmt::dict_has_key(*pMetadata, pmt::mp("freq")))
            *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("freq"),
                                       pmt::mp(d_currentFreqShiftDelta));
          if (!pmt::dict_has_key(*pMetadata, pmt::mp("freqoffset")))
            *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("freqoffset"),
                                       pmt::mp(d_currentFreqShiftDelta));
          if (!pmt::dict_has_key(*pMetadata, pmt::mp("signalcenterfreq")))
            *pMetadata =
                pmt::dict_add(*pMetadata, pmt::mp("signalcenterfreq"),
                              pmt::mp(signalVector[closestIndex].centerFreqHz));
          if (!pmt::dict_has_key(*pMetadata, pmt::mp("trackingcenterfreq")))
            *pMetadata =
                pmt::dict_add(*pMetadata, pmt::mp("trackingcenterfreq"),
                              pmt::mp(signalVector[closestIndex].centerFreqHz));

          message_port_pub(pmt::mp("freq_shift"),
                           pmt::from_double(d_currentFreqShiftDelta));

          pmt::pmt_t pdu = pmt::cons(*pMetadata, pmt::PMT_NIL);
          if (!testMode)
            message_port_pub(pmt::mp("freq_info"), pdu);
        }
      }
    }

    if (d_currentFreqShiftDelta == 0.0) {
      // Not shifting so just memcpy
      // std::cout << "[Debug] Signal present but shift = 0.0" << std::endl;
      memcpy((void *)out, (void *)in, noutput_items * sizeof(gr_complex));
    } else {
      // shift signal
      // std::cout << "[Debug] Shifting signal at " <<
      // signalVector[closestIndex].centerFreqHz << " by " <<
      // d_currentFreqShiftDelta << "Hz" << std::endl;
      // Generate next block of signals direct into output
      d_nco.sincos(out, noutput_items, 1.0);
      // multiply times input.
      volk_32fc_x2_multiply_32fc(out, out, in, noutput_items);
    }
  } else {
    if (!inHoldDown) {
      // Not shifting so just memcpy
      // std::cout << "[Debug] No signal.  Just copying. (NumSignals = " <<
      // numSignals << ")" << std::endl;
      memcpy((void *)out, (void *)in, noutput_items * sizeof(gr_complex));
    } else {
      // shift signal
      // std::cout << "[Debug] numsignals=0 but in hold-down: " <<
      // signalVector[closestIndex].centerFreqHz << " by " <<
      // d_currentFreqShiftDelta << "Hz" << std::endl; Generate next block of
      // signals direct into output
      d_nco.sincos(out, noutput_items, 1.0);
      // multiply times input.
      volk_32fc_x2_multiply_32fc(out, out, in, noutput_items);
    }
  }

  // Have to do extra copies if we're sending the PDU so only send it if we're
  // configured to.
  if (d_processMessages)
    sendMessageData(out, noutput_items, signalVector[closestIndex].centerFreqHz,
                    signalVector[closestIndex].widthHz,
                    signalVector[closestIndex].maxPower, pMetadata);
  // ---------------------------------------------------------------------------

  // PDU Output:
  // signalState:
  // 1 - Signal just acquired
  // 0 - Signal just lost
  //

  // If just detected signal, send new PDU
  if (justDetectedSignal) {
    pmt::pmt_t meta = pmt::make_dict();

    meta = pmt::dict_add(meta, pmt::mp("state"), pmt::mp(1));
    meta = pmt::dict_add(meta, pmt::mp("decisionvalue"),
                         pmt::mp((int)signalVector.size()));
    meta = pmt::dict_add(meta, pmt::mp("numsignals"),
                         pmt::mp((int)signalVector.size()));
    meta = pmt::dict_add(meta, pmt::mp("closestsignal"),
                         pmt::mp(closestIndex + 1));

    pmt::pmt_t pdu = pmt::cons(meta, pmt::PMT_NIL);
    if (!testMode)
      message_port_pub(pmt::mp("signaldetect"), pdu);

    sendState(true);
  }
  // if Just lost signal, send PDU
  if (lostSignal) {
    pmt::pmt_t meta = pmt::make_dict();
    meta = pmt::dict_add(meta, pmt::mp("state"), pmt::mp(0));
    meta = pmt::dict_add(meta, pmt::mp("decisionvalue"), pmt::mp(0));
    meta = pmt::dict_add(meta, pmt::mp("numsignals"),
                         pmt::mp((int)signalVector.size()));

    pmt::pmt_t pdu = pmt::cons(meta, pmt::PMT_NIL);
    if (!testMode)
      message_port_pub(pmt::mp("signaldetect"), pdu);

    sendState(false);
  }

  // Tell runtime system how many output items we produced.
  return noutput_items;
}

int AutoDopplerCorrect_impl::work(int noutput_items,
                                  gr_vector_const_void_star &input_items,
                                  gr_vector_void_star &output_items) {
  const gr_complex *in = (const gr_complex *)input_items[0];
  gr_complex *out = (gr_complex *)output_items[0];

  return processData(noutput_items, in, out, NULL);
}

} /* namespace mesa */
} /* namespace gr */
