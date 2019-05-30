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

#ifndef INCLUDED_MESA_MESAENERGYDETECTOR_IMPL_H
#define INCLUDED_MESA_MESAENERGYDETECTOR_IMPL_H

#include <mesa/MesaEnergyDetector.h>
#include "signals_mesa.h"
#include <chrono>
#include <ctime>

using namespace MesaSignals;

namespace gr {
  namespace mesa {

    class MesaEnergyDetector_impl : public MesaEnergyDetector
    {
     protected:
        boost::mutex d_mutex;
    	EnergyAnalyzer *pEnergyAnalyzer;
        gr_complex *pMsgOutBuff;
        int msgBufferSize;

    	float d_sampleRate;
    	float d_centerFreq;
    	float d_minWidthHz;
    	float d_maxWidthHz;
    	int d_framesToAvg;

    	int d_fftSize;
    	bool d_enableDebug;

    	bool d_genSignalPDUs;

    	std::chrono::time_point<std::chrono::steady_clock> startup, endup;
    	bool d_startInitialized;
    	float d_holdUpSec;

    	// Methods
    	float calcMinDutyCycle();
		virtual int processData(int noutput_items,const gr_complex *in,gr_complex *out,pmt::pmt_t *pMetadata);
		void sendState(bool state);

     public:
      MesaEnergyDetector_impl(int fftsize, float squelchThreshold, float minWidthHz, float maxWidthHz, float radioCenterFreq, float sampleRate, float holdUpSec,
    		  	  	  	  	  int framesToAvg, bool genSignalPDUs, bool enableDebug);
      virtual ~MesaEnergyDetector_impl();

      virtual bool stop();

      void setup_rpc();
      void handleMsgIn(pmt::pmt_t msg);

      virtual float getSquelch() const;
      virtual void setSquelch(float newValue);

      virtual float getCenterFrequency() const;
      virtual void setCenterFrequency(float newValue);

      virtual float getMinWidthHz() const;
      virtual void setMinWidthHz(float newValue);

      virtual float getMaxWidthHz() const;
      virtual void setMaxWidthHz(float newValue);

      // Where all the action really happens
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);
    };

  } // namespace mesa
} // namespace gr

#endif /* INCLUDED_MESA_MESAENERGYDETECTOR_IMPL_H */

