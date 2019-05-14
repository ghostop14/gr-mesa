# gr-mesa - GNURadio Modules for Enhanced Signal Analysis

## Overview
gr-mesa is a project to incorporate some enhanced fundamental signal identification and processing blocks to assist with signal input select and downstream analysis.  Current blocks include:


WARNING: THESE MODULES ARE IN ACTIVE DEVELOPMENT, AND STILL NEED SOME MORE TESTING!  But folks had been asking for early access, so I've posted them.


1. Max Power Detection - This block will analyze the input signal and based on some parameters that control the length of time / averaging will output a max power message.  The block can also, given a threshold value, output a state change (max power above the threshold / below the threshold) when the threshold is crossed.  Holddown timers prevent bouncing.  This can then be used downstream if signal detection can be based solely on power levels.  (e.g. good signal filtering in a dedicated band where seeing a signal power above a noise floor is sufficient to activate downstream processing).  This block can also optionally be used to transition from stream inputs to message-based data outputs.
2. Source Selector - This block monitors the message data inputs for a meta tag "decisionvalue" associated with the input data.  Whichever input has the maximum decision value is the one whose data is sent to the output.  A hold-down timer is available to limit "bouncing".  This can be combined with the MaxPower block to select the input with the best signal strength to continue downstream for processing.  This block does provide some buffering to prevent jittery signals due to a lack of samples, however it will not by itself account for delays between signals due to the time variations between multiple receivers receiving the same signal.
3. Auto Doppler Correct - This block scans the input signal for a signal near the center frequency and attempts to keep the center frequency of the detected signal centered by automatically shifting the signal.  This can be useful if you have unkown or dynamic doppler shifting going on.  If you would like to switch to processing the output as a PDU at this block, enable PDU processing.  This will enable the msgout port.  
4. Signal Detector - This block scans the input signal looking for sub signals of the specified min/max width.  The block takes a max-hold average to inspect the spectrum, then determines any signals present in the spectrum.  Note that the FFT frames to average is configurable.  Too small and the detection is jittery, too high and too many samples will be held/processed so pick a number that works well (or stick with the default).  When a signal is detected, a PDU will be generated on the signaldetect connector with a 'state' metadata tag set to 1.  PDU's are only sent on state changes, so any downstream blocks should track their own state.  When no signals are present and the hold timer has expired, a PDU will be generated with 'state' set to 0.  If more downstream data processing is desired, 'Gen Signal PDUs' can be turned on.  In that case, for each detected signal, a PDU is generated along with some metadata (radio freq, sample rate, signal center freq, signal width, and signal max power) along with the full data block.  This can be used downstream to tune filters and/or shift the signal.  

## Building
gr-mesa has no core dependencies.  However if you will be using the state out ports, it is highly recommended to install gr-filerepeater as additional state blocks are included there.

``
cd <clone directory>

mkdir build

cd build

cmake ..

make

[sudo] make install

sudo ldconfig
``

If each step was successful (do not overlook the "sudo ldconfig" step if this is the first installation).

