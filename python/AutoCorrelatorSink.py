#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Copyright 2019 ghostop14.
# 
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# 
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
# 

from gnuradio import gr
from gnuradio import blocks
from mesa import AutoCorrelator
from mesa import Normalize
from gnuradio import qtgui
# SIP import lets us wrap our control as a pyQt widget
import sip
from PyQt5 import QtGui
from PyQt5.QtWidgets import QWidget

class AutoCorrelatorSink(gr.hier_block2):
    """
    docstring for block AutoCorrelatorSink
    """
    def __init__(self, sample_rate, fac_size, fac_decimation, title,  autoScale, grid, yMin, yMax,  useDB):
        gr.hier_block2.__init__(self,
            "AutoCorrelatorSink",
            gr.io_signature(1, 1, gr.sizeof_gr_complex),  # Input signature
            gr.io_signature(0, 0, 0)) # Output signature

        self.fac_size = fac_size
        self.fac_decimation = fac_decimation
        self.sample_rate = sample_rate
        
        autoCorr = AutoCorrelator(sample_rate, fac_size, fac_decimation,  useDB)
        vecToStream = blocks.vector_to_stream(gr.sizeof_float, self.fac_size)

        self.timeSink = qtgui.time_sink_f(self.fac_size/2, sample_rate, title, 1)
        self.timeSink.enable_grid(grid)
        self.timeSink.set_y_axis(yMin, yMax)
        self.timeSink.enable_autoscale(autoScale)
        self.timeSink.disable_legend()
        self.timeSink.set_update_time(0.1)

        if useDB:
            self.connect(self, autoCorr,  vecToStream,  self.timeSink)
        else:
            norm = Normalize.Normalize(self.fac_size)
            self.connect(self, autoCorr,  norm,  vecToStream,  self.timeSink)

        #pyQt  = self.timeSink.pyqwidget()
        #self.pyWin = sip.wrapinstance(pyQt, QtGui.QWidget)
        # self.pyWin.show()

    def getWidget(self):
        return sip.wrapinstance(self.timeSink.pyqwidget(), QWidget)
