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


import numpy
from gnuradio import gr
import pmt

import threading
import time
import sys

# Cheat to send messages across threads
from PyQt5.QtWidgets import QFrame
from PyQt5 import QtCore

class RotationThread(threading.Thread):
    def __init__(self, parentCtl, valueList,rotationTime):
        threading.Thread.__init__(self)
        
        self.stopThread = False
        self.threadRunning = False
        
        self.hold = False
        
        self.parentCtl = parentCtl
        self.valueList = valueList
        self.rotationTime = rotationTime
        self.curValue = valueList[0]
    
    def run(self):
        self.threadRunning = True
        
        curIndex = 0
        numEntries = len(self.valueList)
        
        time.sleep(0.5) # TODO: Find better way to know if init is all done
        
        while not self.stopThread:
            if not self.hold:
                self.curValue = self.valueList[curIndex]
                self.parentCtl.sigNewValue.emit(self.curValue,curIndex)
            
                # Need to shorten detection time for threadRunning
                time.sleep(self.rotationTime)
                
                curIndex = curIndex + 1
                if curIndex > (numEntries-1):
                    curIndex = 0
            else:
                # Just a short sleep waiting to see if the hold gets released
                time.sleep(0.002)
        
        self.threadRunning = False
    
class VariableRotator(gr.sync_block,QFrame):
    """
    docstring for block VariableRotator
    """
    sigNewValue = QtCore.pyqtSignal(float,int)
        
    def __init__(self, valueList,rotationTime,varCallback,parent):
        gr.sync_block.__init__(self,name="VariableRotator",in_sig=None,out_sig=None)
        QFrame.__init__(self,parent=parent)

        if (len(valueList) == 0):
            print("[ValueRotator] Error: Please provide a list of values")
            sys.exit(1)

        self.varCallback = varCallback
        self.sigNewValue.connect(self.newValue)
        
        self.lastValue = valueList[0]
        
        self.thread = RotationThread(self, valueList, rotationTime)
        self.thread.start()
        
        self.message_port_register_in(pmt.intern("hold"))
        self.set_msg_handler(pmt.intern("hold"), self.holdHandler)   
        
        self.message_port_register_out(pmt.intern("value"))
        self.message_port_register_out(pmt.intern("index"))

    def holdHandler(self, msg):
        try:    
            newVal = pmt.to_python(pmt.cdr(msg))

            if type(newVal) == int:
                if newVal == 0:
                    self.thread.hold = False
                else:
                    self.thread.hold = True
            else:
                print("[ValueRotator] Error: Value received was not an int: %s" % str(e))
                
        except Exception as e:
            print("[ValueRotator] Error with message conversion: %s" % str(e))
        
    def newValue(self,curValue,curIndex):
        p_val = pmt.from_float(curValue)
        p_index = pmt.from_long(curIndex)

        if curValue != self.lastValue:
            try:
                self.varCallback(curValue)
            except:
                pass
            
            self.message_port_pub(pmt.intern("value"),pmt.cons(pmt.intern("value"),p_val))
            
            self.lastValue = curValue
            
        self.message_port_pub(pmt.intern("index"),pmt.cons(pmt.intern("index"),p_index))
    
    def stop(self):
        self.thread.stopThread = True
        self.thread.join()
            
        return True
      

    def work(self, input_items, output_items):
        out = output_items[0]
        return len(output_items[0])

