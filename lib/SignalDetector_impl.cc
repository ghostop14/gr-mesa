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

#include <cstring> // memcpy

#include "SignalDetector_impl.h"
#include <gnuradio/io_signature.h>
#include <functional>
namespace gr {
namespace mesa {

SignalDetector::sptr SignalDetector::make(int fftsize,
                                          float squelchThreshold,
                                          double minWidthHz,
                                          double maxWidthHz,
                                          double radioCenterFreq,
                                          double sampleRate,
                                          float holdUpSec,
                                          int framesToAvg,
                                          bool genSignalPDUs,
                                          bool enableDebug,
                                          int detectionMethod)
{
    return gnuradio::make_block_sptr<SignalDetector_impl>(fftsize,
                                                          squelchThreshold,
                                                          minWidthHz,
                                                          maxWidthHz,
                                                          radioCenterFreq,
                                                          sampleRate,
                                                          holdUpSec,
                                                          framesToAvg,
                                                          genSignalPDUs,
                                                          enableDebug,
                                                          detectionMethod);
}

/*
 * The private constructor
 */
SignalDetector_impl::SignalDetector_impl(int fftsize, float squelchThreshold,
                                         double minWidthHz, double maxWidthHz,
                                         double radioCenterFreq,
                                         double sampleRate, float holdUpSec,
                                         int framesToAvg, bool genSignalPDUs,
                                         bool enableDebug, int detectionMethod)
    : gr::sync_block("SignalDetector",
                     gr::io_signature::make(1, 1, sizeof(gr_complex)),
                     gr::io_signature::make(1, 1, sizeof(gr_complex))) {
  pMsgOutBuff = NULL;

  // Store variables
  d_fftSize = fftsize;
  d_centerFreq = radioCenterFreq;
  d_minWidthHz = minWidthHz;
  d_maxWidthHz = maxWidthHz;
  d_sampleRate = sampleRate;
  d_holdUpSec = holdUpSec;
  d_enableDebug = enableDebug;
  d_framesToAvg = framesToAvg;

  // Init some attributes
  d_startInitialized = false;

  d_genSignalPDUs = genSignalPDUs;

  // Calc Duty Cycle
  float minDutyCycle = calcMinDutyCycle();

  // Create energy analyzer
  pEnergyAnalyzer = new EnergyAnalyzer(fftsize, squelchThreshold, minDutyCycle);
  d_detectionMethod = detectionMethod;

  // Make sure we have a multiple of fftsize coming in
  gr::block::set_output_multiple(fftsize * d_framesToAvg);

  // Set up PDUs
  message_port_register_in(pmt::mp("msgin"));
  set_msg_handler(pmt::mp("msgin"),
                  std::bind(&SignalDetector_impl::handleMsgIn, this, std::placeholders::_1));

  message_port_register_out(pmt::mp("signaldetect"));
  message_port_register_out(pmt::mp("signals"));
  message_port_register_out(pmt::mp("state"));
}

float SignalDetector_impl::calcMinDutyCycle() {
  double hzPerBucket = d_sampleRate / (double)d_fftSize;
  double binsForMinHz = d_minWidthHz / hzPerBucket;
  double minDutyCycle = binsForMinHz / d_fftSize;

  return minDutyCycle;
}

float SignalDetector_impl::getSquelch() const {
  return pEnergyAnalyzer->getThreshold();
}

void SignalDetector_impl::setSquelch(float newValue) {
  pEnergyAnalyzer->setThreshold(newValue);

  if (d_enableDebug)
    std::cout << "[Mesa Detector] Changing squelch to " << newValue
              << std::endl;
}

double SignalDetector_impl::getCenterFrequency() const { return d_centerFreq; }

void SignalDetector_impl::setCenterFrequency(double newValue) {
  gr::thread::scoped_lock guard(d_mutex);

  if (d_startInitialized) {
    pmt::pmt_t meta = pmt::make_dict();

    meta = pmt::dict_add(meta, pmt::mp("state"), pmt::mp(0));

    pmt::pmt_t pdu = pmt::cons(meta, pmt::PMT_NIL);
    message_port_pub(pmt::mp("signaldetect"), pdu);

    d_startInitialized = 0.0;
  }

  d_centerFreq = newValue;

  if (d_enableDebug)
    std::cout << "[Mesa Detector] Changing frequency to " << newValue
              << std::endl;
}

double SignalDetector_impl::getMinWidthHz() const { return d_minWidthHz; }

void SignalDetector_impl::setMinWidthHz(double newValue) {
  // Calc Duty Cycle
  d_minWidthHz = newValue;
  float minDutyCycle = calcMinDutyCycle();
  pEnergyAnalyzer->setDutyCycle(minDutyCycle);

  if (d_enableDebug)
    std::cout << "[Mesa Detector] Changing min width (Hz) to " << newValue
              << std::endl;
}

double SignalDetector_impl::getMaxWidthHz() const { return d_maxWidthHz; }

void SignalDetector_impl::setMaxWidthHz(double newValue) {
  d_maxWidthHz = newValue;
  if (d_enableDebug)
    std::cout << "[Mesa Detector] Changing max width (Hz) to " << newValue
              << std::endl;
}

bool SignalDetector_impl::stop() {
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
SignalDetector_impl::~SignalDetector_impl() { bool retVal = stop(); }

void SignalDetector_impl::handleMsgIn(pmt::pmt_t msg) {
  if (!d_genSignalPDUs)
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

void SignalDetector_impl::sendState(bool state) {
  int newState;
  if (state) {
    newState = 1;
  } else {
    newState = 0;
  }

  pmt::pmt_t pdu = pmt::cons(pmt::intern("state"), pmt::from_long(newState));

  message_port_pub(pmt::mp("state"), pdu);
}

int SignalDetector_impl::processData(int noutput_items, const gr_complex *in,
                                     gr_complex *out, pmt::pmt_t *pMetadata) {
  gr::thread::scoped_lock guard(d_mutex);
  // First get the max hold curve for this block
  FloatVector maxSpectrum;
  long samplesProcessed =
      pEnergyAnalyzer->maxHold(in, noutput_items, maxSpectrum, true);

  // Now look if we have signals
  int numSignals = 0;
  SignalOverviewVector signalVector;

  if (d_detectionMethod == SIGDETECTOR_METHOD_SEPARATESIGNALS) {
    // Last param says stop looking on the first detected signal.
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

  // If we have signals, deal with that (set PDU, etc.)
  // Tell runtime system how many output items we produced.
  if (numSignals > 0) {
    memcpy((void *)out, (void *)in, noutput_items * sizeof(gr_complex));

  } else {
    memset(out, 0, noutput_items * sizeof(gr_complex));
  }

  // start initialized tracks if we've picked up a signal and we're in a "high"
  // / have signal state
  bool justDetectedSignal = false; // first detection
  bool lostSignal = false;
  bool signalPresent = false;
  bool inHoldDown = false;

  if (numSignals > 0) { // We have a detection
    signalPresent = true;

    if (!d_startInitialized) {
      // Haven't seen a signal in a while, start the rising edge.
      startup = std::chrono::steady_clock::now();
      endup = startup;
      d_startInitialized = true;
      justDetectedSignal = true;
      if (d_enableDebug)
        std::cout << "[Mesa Detector] Just detected signal." << std::endl;
    } else {
      // We're continuing to see a signal.  Move the end indicator
      endup = std::chrono::steady_clock::now();
    }
  } else {                    // No Detection
    if (d_startInitialized) { // We had a signal so we can track losing it.
      // Before we say we've lost it, let's see if we're within our hold timer
      std::chrono::time_point<std::chrono::steady_clock> curTimestamp =
          std::chrono::steady_clock::now();
      std::chrono::duration<double> elapsed_seconds = curTimestamp - endup;
      if (elapsed_seconds.count() > (double)d_holdUpSec) {
        // No detection and we've exceeded our hold window.  Reset start
        // tracker.
        d_startInitialized = false;
        lostSignal = true;
        if (d_enableDebug)
          std::cout << "[Mesa Detector] Just lost signal." << std::endl;
      } else {
        inHoldDown = true;
      }
    }
  }

  // PDU Output:
  // signalState:
  // 1 - Signal just acquired
  // 0 - Signal just lost
  //
  // numSignals
  // Number of signals if signals were detected

  // If just detected signal, send new PDU
  if (justDetectedSignal) {
    // Find the max power signal:
    double maxCtrFreq = 0.0;
    double maxWidth = 0.0;
    float maxPower = -999.0;

    for (int i = 0; i < signalVector.size(); i++) {
      if (signalVector[i].maxPower > maxPower) {
        maxCtrFreq = signalVector[i].centerFreqHz;
        maxWidth = signalVector[i].widthHz;
        maxPower = signalVector[i].maxPower;
      }
    }

    pmt::pmt_t meta = pmt::make_dict();

    meta = pmt::dict_add(meta, pmt::mp("state"), pmt::mp(1));
    meta = pmt::dict_add(meta, pmt::mp("decisionvalue"),
                         pmt::mp((int)signalVector.size()));
    meta = pmt::dict_add(meta, pmt::mp("numsignals"),
                         pmt::mp((int)signalVector.size()));
    meta = pmt::dict_add(meta, pmt::mp("radioFreq"), pmt::mp(d_centerFreq));
    meta = pmt::dict_add(meta, pmt::mp("sampleRate"), pmt::mp(d_sampleRate));
    meta = pmt::dict_add(meta, pmt::mp("strongestCenterFreq"),
                         pmt::mp(maxCtrFreq));
    meta = pmt::dict_add(meta, pmt::mp("strongestWidthHz"), pmt::mp(maxWidth));
    meta = pmt::dict_add(meta, pmt::mp("strongestPower"), pmt::mp(maxPower));

    pmt::pmt_t pdu = pmt::cons(meta, pmt::PMT_NIL);
    message_port_pub(pmt::mp("signaldetect"), pdu);

    sendState(true);
  }
  // if Just lost signal, send PDU
  if (lostSignal) {
    pmt::pmt_t meta = pmt::make_dict();

    meta = pmt::dict_add(meta, pmt::mp("state"), pmt::mp(0));
    meta = pmt::dict_add(meta, pmt::mp("decisionvalue"), pmt::mp(0));
    meta = pmt::dict_add(meta, pmt::mp("radioFreq"), pmt::mp(d_centerFreq));
    meta = pmt::dict_add(meta, pmt::mp("sampleRate"), pmt::mp(d_sampleRate));

    pmt::pmt_t pdu = pmt::cons(meta, pmt::PMT_NIL);
    message_port_pub(pmt::mp("signaldetect"), pdu);

    sendState(false);
  }

  // This takes some processing, so we only do this if it's requested.
  if (d_genSignalPDUs) {
    pmt::pmt_t data_out(pmt::init_c32vector(noutput_items, in));

    for (int i = 0; i < signalVector.size(); i++) {
      pmt::pmt_t meta = pmt::make_dict();

      if (numSignals > 0) {
        if (!pMetadata) {
          meta =
              pmt::dict_add(meta, pmt::mp("radioFreq"), pmt::mp(d_centerFreq));
          meta =
              pmt::dict_add(meta, pmt::mp("sampleRate"), pmt::mp(d_sampleRate));
          meta = pmt::dict_add(meta, pmt::mp("signalCenterFreq"),
                               pmt::mp(signalVector[i].centerFreqHz));
          meta = pmt::dict_add(meta, pmt::mp("widthHz"),
                               pmt::mp(signalVector[i].widthHz));
          meta = pmt::dict_add(meta, pmt::mp("maxPower"),
                               pmt::mp(signalVector[i].maxPower));
        } else {
          if (!pmt::dict_has_key(*pMetadata, pmt::mp("radioFreq")))
            *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("radioFreq"),
                                       pmt::mp(d_centerFreq));
          if (!pmt::dict_has_key(*pMetadata, pmt::mp("sampleRate")))
            *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("sampleRate"),
                                       pmt::mp(d_sampleRate));
          if (!pmt::dict_has_key(*pMetadata, pmt::mp("signalCenterFreq")))
            *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("signalCenterFreq"),
                                       pmt::mp(signalVector[i].centerFreqHz));
          if (!pmt::dict_has_key(*pMetadata, pmt::mp("widthHz")))
            *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("widthHz"),
                                       pmt::mp(signalVector[i].widthHz));
          if (!pmt::dict_has_key(*pMetadata, pmt::mp("maxPower")))
            *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("maxPower"),
                                       pmt::mp(signalVector[i].maxPower));
        }
      } else {
        if (!pMetadata) {
          meta =
              pmt::dict_add(meta, pmt::mp("radioFreq"), pmt::mp(d_centerFreq));
          meta =
              pmt::dict_add(meta, pmt::mp("sampleRate"), pmt::mp(d_sampleRate));
          meta = pmt::dict_add(meta, pmt::mp("signalCenterFreq"), pmt::mp(0.0));
          meta = pmt::dict_add(meta, pmt::mp("widthHz"), pmt::mp(0.0));
          meta = pmt::dict_add(meta, pmt::mp("maxPower"), pmt::mp(0.0));
        } else {
          if (!pmt::dict_has_key(*pMetadata, pmt::mp("radioFreq")))
            *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("radioFreq"),
                                       pmt::mp(d_centerFreq));
          if (!pmt::dict_has_key(*pMetadata, pmt::mp("sampleRate")))
            *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("sampleRate"),
                                       pmt::mp(d_sampleRate));
          if (!pmt::dict_has_key(*pMetadata, pmt::mp("signalCenterFreq")))
            *pMetadata = pmt::dict_add(*pMetadata, pmt::mp("signalCenterFreq"),
                                       pmt::mp(0.0));
          if (!pmt::dict_has_key(*pMetadata, pmt::mp("widthHz")))
            *pMetadata =
                pmt::dict_add(*pMetadata, pmt::mp("widthHz"), pmt::mp(0.0));
          if (!pmt::dict_has_key(*pMetadata, pmt::mp("maxPower")))
            *pMetadata =
                pmt::dict_add(*pMetadata, pmt::mp("maxPower"), pmt::mp(0.0));
        }
      }

      if (!pMetadata) {
        pmt::pmt_t pdu = pmt::cons(meta, data_out);
        message_port_pub(pmt::mp("signals"), pdu);
      } else {
        pmt::pmt_t pdu = pmt::cons(*pMetadata, data_out);
        message_port_pub(pmt::mp("signals"), pdu);
      }
    }
  }

  return noutput_items;
}

int SignalDetector_impl::work(int noutput_items,
                              gr_vector_const_void_star &input_items,
                              gr_vector_void_star &output_items) {
  const gr_complex *in = (const gr_complex *)input_items[0];
  gr_complex *out = (gr_complex *)output_items[0];

  return processData(noutput_items, in, out, NULL);
} // end work

} /* namespace mesa */
} /* namespace gr */
