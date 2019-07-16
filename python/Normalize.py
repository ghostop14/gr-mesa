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

class Normalize(gr.hier_block2):
	def __init__(self, vecsize=1024):
		gr.hier_block2.__init__(
			self, "Normalize",
			gr.io_signature(1, 1, gr.sizeof_float*vecsize),
			gr.io_signature(1, 1, gr.sizeof_float*vecsize),
		)

		##################################################
		# Parameters
		##################################################
		self.vecsize = vecsize

		##################################################
		# Blocks
		##################################################
		self.blocks_stream_to_vector_0 = blocks.stream_to_vector(gr.sizeof_float, vecsize)
		self.blocks_repeat_0 = blocks.repeat(gr.sizeof_float, vecsize)
		self.blocks_max_xx_0 = blocks.max_ff(vecsize)
		self.blocks_divide_xx_0 = blocks.divide_ff(vecsize)

		##################################################
		# Connections
		##################################################
		self.connect((self.blocks_divide_xx_0, 0), (self, 0))
		self.connect((self.blocks_stream_to_vector_0, 0), (self.blocks_divide_xx_0, 1))
		self.connect((self, 0), (self.blocks_max_xx_0, 0))
		self.connect((self.blocks_repeat_0, 0), (self.blocks_stream_to_vector_0, 0))
		self.connect((self.blocks_max_xx_0, 0), (self.blocks_repeat_0, 0))
		self.connect((self, 0), (self.blocks_divide_xx_0, 0))


	def get_vecsize(self):
		return self.vecsize

	def set_vecsize(self, vecsize):
		self.vecsize = vecsize

