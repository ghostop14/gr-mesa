# Copyright 2020 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

########################################################################
# Check if there is C++ code at all
########################################################################
if(NOT mesa_sources)
    MESSAGE(STATUS "No C++ sources... skipping python bindings")
    return()
endif(NOT mesa_sources)

########################################################################
# Check for pygccxml
########################################################################
GR_PYTHON_CHECK_MODULE_RAW(
    "pygccxml"
    "import pygccxml"
    PYGCCXML_FOUND
    )

include(GrPybind)

########################################################################
# Python Bindings
########################################################################

list(APPEND mesa_python_files
    AutoDopplerCorrect_python.cc
    AvgToMsg_python.cc
    ioselector_python.cc
    LongTermIntegrator_python.cc
    MaxPower_python.cc
    phase_shift_python.cc
    SignalDetector_python.cc
    SourceSelector_python.cc python_bindings.cc)

GR_PYBIND_MAKE_OOT(mesa 
   ../..
   gr::mesa
   "${mesa_python_files}")

install(TARGETS mesa_python DESTINATION ${GR_PYTHON_DIR}/mesa COMPONENT pythonapi)
