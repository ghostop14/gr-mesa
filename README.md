# gr-mesa - GNURadio Modules for Enhanced Signal Analysis

## Overview
gr-mesa is a project to incorporate some enhanced fundamental signal identification and processing blocks to assist with signal input select and downstream analysis.  Current blocks include:


1. Max Power Detection - This block will analyze the input signal and based on some parameters that control the length of time / averaging will output a max power message.  The block can also, given a threshold value, output a state change (max power above the threshold / below the threshold) when the threshold is crossed.  Holddown timers prevent bouncing.  This can then be used downstream if signal detection can be based solely on power levels.  (e.g. good signal filtering in a dedicated band where seeing a signal power above a noise floor is sufficient to activate downstream processing).  This block can also optionally be used to transition from stream inputs to message-based data outputs.
2. Source Selector - This block monitors the message data inputs for a meta tag "decisionvalue" associated with the input data.  Whichever input has the maximum decision value is the one whose data is sent to the output.  A hold-down timer is available to limit "bouncing".  This can be combined with the MaxPower block to select the input with the best signal strength to continue downstream for processing.  This block does provide some buffering to prevent jittery signals due to a lack of samples, however it will not by itself account for delays between signals due to the time variations between multiple receivers receiving the same signal.
3. Auto Doppler Correct - This block scans the input signal for a signal near the center frequency and attempts to keep the center frequency of the detected signal centered by automatically shifting the signal.  This can be useful if you have unkown or dynamic doppler shifting going on.  If you would like to switch to processing the output as a PDU at this block, enable PDU processing.  This will enable the msgout port.
4. Signal Detector - This block scans the input signal looking for sub signals of the specified min/max width.  The block takes a max-hold average to inspect the spectrum, then determines any signals present in the spectrum.  Note that the FFT frames to average is configurable.  Too small and the detection is jittery, too high and too many samples will be held/processed so pick a number that works well (or stick with the default).  When a signal is detected, a PDU will be generated on the signaldetect connector with a 'state' metadata tag set to 1.  PDU's are only sent on state changes, so any downstream blocks should track their own state.  When no signals are present and the hold timer has expired, a PDU will be generated with 'state' set to 0.  If more downstream data processing is desired, 'Gen Signal PDUs' can be turned on.  In that case, for each detected signal, a PDU is generated along with some metadata (radio freq, sample rate, signal center freq, signal width, and signal max power) along with the full data block.  This can be used downstream to tune filters and/or shift the signal.
5. A QT GUI version of the Fast Auto-correlator (example in the examples directory).  This conversion makes this block GR 3.8/3.9-Ready.
6. A fast auto-correlator block that provides correlated vectors as output (example in the examples directory).
7. Normalize - Take an input vector and normalize all values to 1.0.
8. Phase Shift - Shift an incoming signal by shift_radians.  Shift can be controlled via variable or incoming float message
9. Average to Message - Take the average of an incoming float vector and output the scalar average as a message
10. Variable Rotator - While named and functioning more generically, the drive behind this block was a block that could rotate frequencies as part of a GNURadio-based scanner.  The block id can be used as a variable, and 2 messages get output: One is the current value, and one is the corresponding index from the list of provided values.  The index facilitates different downstream processing paths using the IO Selector for each value.  For instance, if the list is frequencies, f1 may be NBFM, f2 may be digital, etc.  The block also has a message input that can be used to lock/hold a frequency where activity is detected.  See the Scanner section below for more details.

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

## GNURadio-Based Scanner
One exciting solution that could be developed with gr-mesa (**note this also requires gr-lfast**) is a complete GNURadio-based scanner.  Two basic examples are included in the examples directory.  (Note the 3 frequencies is not a limitation, just a setting in the design of the flowgraph for this example)
1. The first flowgraph scans for voice NBFM signals on 3 different channels and allows for 3 different paths of decoding (examples/scanner_fm.grc). 
2. The second folowgraph uses the same track for all decodes. (examples/scanner_fm_single_decode.grc)
 
The key component behind enabling a scanner in GNURadio is the Variable Rotator block in this OOT module, which provides the fundamentals to iterate through a frequency list at set time intervals.  The rotator also outputs an index corresponding to the configured list so that different downstream processing paths can be taken for each frequency (if that's how you would like to use it).  The first Signal Detector block is then combined with this to detect when a signal is actually present.  This mimics the basic scanner function of "is there a signal present?  If so, stop here."  The state output from the Signal Detector goes high when a signal is detected, when matches up with the hold input of the rotator block creating the necessary feedback loop to hold on a channel when a signal is detected.

In the first example flowgraph, each path for 3 different frequencies is looking for a NBFM audio/analog signal.  Each path uses a separate signal detector such that when the processing holds on a channel, the individual channel signal detector goes high telling an Advanced File Sink from the **gr-filerepeater** OOT module to start recording the signal to a WAV file that can be played back with any WAV file player.  (Note: There were some recent updates to the Adv Sink block to support the scanning functionality, so if you already had it installed, please git pull and refresh it).  When the signal goes away on the active frequency and the variable rotator's hold is released and it goes to the next frequency, the individual channel detector will transition low after a hold period and close the file.  The net result of this whole process is individual recordings for each signal detection on each channel saved in files named and timestamped corresponding to their frequency.

The second example with a single track capitalizes on the Adv File Sink's capability to rotate files when the frequency changes.  If all scanned channels are the same type (in this example NBFM), this is a more efficient approach as it cuts out adding individual signal detectors.

Both examples use 2 blocks from the **gr-guiextra** OOT module for some better visualization to complete the scanner.  The first is a familiar frequency / digital number display.  The other is a push / toggle button.  Note gr-mesa's note about the digital display, make sure you 'sudo pip3 install pyqtchart' (or pip install if you're still on python2).  There's a version issue with the native apt version, even in Ubuntu 18.04 as it only includes version 5.9 and version 5.11 or better is required for a specific function.  If you have issues with this approach you can always remove the frequency display from the flowgraph and use a standard GNURadio control instead.  THe second control is a toggle button in gr-guiextra that will stay down when pressed and generate messages on state changes.  When combined with the State Message Or block from gr-filerepeater as demonstrated in the example flowgraphs, you can within the flowgraph dynamically HOLD or lock onto the current frequency.  This tells the rotator to not go to the next frequency until the hold is released and there is no signal. With all of that said, this is a very basic but functioning example of a scanner implemented in GNURadio, and the flowgraphs provide a framework for more advanced processing depending on your needs (this part's up to you).

From these examples, it's up to you how complex you make it.  A couple of suggestions to keep in mind:
1. Test a single processing track in an isolated flowgraph before thinking something isn't working.  And watch any decimations along the way if you integrate a number of stand-alone flowgraphs into one with frequency rotation.
2. The rotation most likely won't be timing-accurate enough to follow say FHSS, so if you try to put your own together to do that, it probably won't work.
3. Watch how different the frequencies are relative to the tuning of your antenna.  Antenna rules still apply.  An antenna tuned to 2m won't be optimal on 70cm, etc.
4. In the basic example flowgraph provided, you'll see different thresholds set in each track for the signal detector squelch thresholds.  This is specifically to adjust differences observed due to (3) with the example frequencies used.  You'll need to manually monitor and experiment with the best values here depending on the frequencies you select and overall design.
5. See the developer's notes below about non QtGUI flowgraphs with the rotator.

### Variable Rotator Technical / Developer Notes
There are a few "tricks" in the Variable Rotator block worth mentioning.  First, because a separate thread is used to control the scheduling of rotation and messages, sending pmt messages within the Qt GUI context runs into exceptions sending cross-thread.  As a result, the workaround was to create the block as a QFrame and leverage Qt's signaling mechanisms to queue it into the message queue of the primary thread with an emit() call.  So no GUI control is visually displayed, but one is used behind the scenes to allow for cross-thread behavior to work as expected.  While not tested, this COULD mean that using the variable rotator may not work in non-QtGUI flowgraphs.  Just something to keep in mind.

