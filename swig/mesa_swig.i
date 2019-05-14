/* -*- c++ -*- */

#define MESA_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "mesa_swig_doc.i"

%{
#include "mesa/MesaEnergyDetector.h"
#include "mesa/AutoDopplerCorrect.h"
#include "mesa/MaxPower.h"
#include "mesa/SourceSelector.h"
%}


%include "mesa/MesaEnergyDetector.h"
GR_SWIG_BLOCK_MAGIC2(mesa, MesaEnergyDetector);
%include "mesa/AutoDopplerCorrect.h"
GR_SWIG_BLOCK_MAGIC2(mesa, AutoDopplerCorrect);
%include "mesa/MaxPower.h"
GR_SWIG_BLOCK_MAGIC2(mesa, MaxPower);
%include "mesa/SourceSelector.h"
GR_SWIG_BLOCK_MAGIC2(mesa, SourceSelector);
