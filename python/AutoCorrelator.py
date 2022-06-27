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
from gnuradio import blocks, fft, filter
import math

class AutoCorrelator(gr.hier_block2):
    """
    docstring for block AutoCorrelator
    """
    def __init__(self, sample_rate,  fac_size,fac_decimation,  useDB):
        gr.hier_block2.__init__(self,"AutoCorrelator",
            gr.io_signature(1, 1, gr.sizeof_gr_complex),  # Input signature
            gr.io_signature(1, 1, gr.sizeof_float*fac_size)) # Output signature

        self.fac_size = fac_size
        self.fac_decimation = fac_decimation
        self.sample_rate = sample_rate
            
        streamToVec = blocks.stream_to_vector(gr.sizeof_gr_complex, self.fac_size)
        # Make sure N is at least 1
        decimation =  int(self.sample_rate/self.fac_size/self.fac_decimation)
        self.one_in_n = blocks.keep_one_in_n(gr.sizeof_gr_complex * self.fac_size, max(1,decimation))

        # FFT Note: No windowing.
        fac = fft.fft_vcc(self.fac_size, True, ())

        complex2Mag = blocks.complex_to_mag(self.fac_size)
        self.avg = filter.single_pole_iir_filter_ff(1.0, self.fac_size)

        fac_fac   = fft.fft_vfc(self.fac_size, True, ())
        fac_c2mag = blocks.complex_to_mag(fac_size)

        # There's a note in Baz's block about needing to add 3 dB to each bin but the DC bin, however it was never implemented
        n = 20
        k =  -20*math.log10(self.fac_size)
        log = blocks.nlog10_ff(n, self.fac_size, k )

        if useDB:
            self.connect(self, streamToVec, self.one_in_n, fac, complex2Mag,  fac_fac, fac_c2mag, self.avg, log,  self)
        else:
            self.connect(self, streamToVec, self.one_in_n, fac, complex2Mag,  fac_fac, fac_c2mag, self.avg, self)
