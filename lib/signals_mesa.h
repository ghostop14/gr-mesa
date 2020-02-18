/*
 * signals_mesa.h
 *
 *      Copyright 2017, Michael Piscopo
 *
 *
 */

#ifndef LIB_SIGNALS_MESA_H_
#define LIB_SIGNALS_MESA_H_

#include "scomplex.h"
#include <boost/thread/mutex.hpp>
#include <fftw3.h>
#include <volk/volk.h>

typedef std::vector<float> FloatVector;
// These match the FFTW definitions
//#define FFTDIRECTION_FORWARD -1
// #define FFTDIRECTION_BACKWARD 1
#define FFTDIRECTION_FORWARD FFTW_FORWARD
#define FFTDIRECTION_BACKWARD FFTW_BACKWARD

#define WINDOWTYPE_NONE 0
#define WINDOWTYPE_HAMMING 1
#define WINDOWTYPE_BLACKMAN_HARRIS 2

#define SQUELCH_DISABLE -1000.0
#define NOISE_FLOOR -100.0

using namespace std;

namespace MesaSignals {
// global test function
void printArray(FloatVector &arr, string name);
void printArray(float *arr, int arrSize, string name);

/*
 * FFT Transforms
 */
class FFT {
public:
  // Actually it appears to get better performance, at least for 1024 sample
  // FFT's with 1 thread The fftw3f doc says multiple threads only really help
  // for large fft sizes, which may be why.
  FFT(int FFTDirection, int initFFTSize, int initNumThreads = 1);
  virtual ~FFT();

  inline int inputBufferLength() const { return fftSize; }
  inline int outputBufferLength() const { return fftSize; }

  inline SComplex *getInputBuffer() { return inputBuffer; };
  inline SComplex *getOutputBuffer() { return outputBuffer; };

  virtual void setWindow(int winType);
  // Note: For this setWindow, the length of newTaps should be fftSize
  virtual void setWindow(FloatVector &newTaps);

  virtual void clearWindow();

  // Execute computes the FFT and if a windowing function has been set for the
  // class, it is applied appropriately before executing the FFT
  virtual void execute(bool shift = false);

  // This execute applies the specified taps rather than the class taps.
  // NOTE: the length of pTaps must be fftSize.
  // Passing NULL will disable windowing for this run.
  inline void execute(float *pTaps, bool shift = false);

  // So PSD computes power in dBm which is basically the RSSI.
  // PowerSpectralDensity: Call Execute first, then call this function
  // to compute PSD and store it in the provided psdBuffer buffer
  // (note: psdBuffer must be at least fftSize in length and should be aligned
  // memory) Note: PowerSpectralDensity output is basically the same as RSSI
  // without squelch void PowerSpectralDensity(float *psdBuffer);
  void PowerSpectralDensity(float *psdBuffer,
                            float squelchThreshold = SQUELCH_DISABLE,
                            float onSquelchSetRSSI = NOISE_FLOOR);

  // RSSI is very similar to PSD in terms of the PSD calculations to dBm,
  // however you can set a squelch threshold Anything below that threshold is
  // set to onSquelchSetRSSI to to "zero" the bin.
  void rssi(float *psdBuffer, float squelchThreshold = SQUELCH_DISABLE,
            float onSquelchSetRSSI = NOISE_FLOOR);

protected:
  boost::mutex d_mutex;

  float log2To10Factor;
  float rssi_K_const;

  string wisdomFilename;
  int fftSize;
  int numThreads;
  void *fftPlan;
  int fftDirection;

  float *alignedWindowTaps;
  bool hasTaps = false;

  SComplex *inputBuffer;
  SComplex *outputBuffer;
  float *tmpBuff; // Used for swapping to center DC in spectrums
  int fftLenBytes;
  int halfFFTSizeBytes;
  int halfFFTSize;

  void init();
  void initThreads();
  void importWisdom();
  void exportWisdom();
};

/*
 * Waterfall and energy analyzers
 */
class SpectrumOverview {
public:
  SpectrumOverview(){};
  virtual ~SpectrumOverview(){};
  SpectrumOverview &operator=(const SpectrumOverview &other);

  float dutyCycle = 0.0;
  float maxPower = NOISE_FLOOR;
  float minPower = 1000.0;
  float centerAvgPower = NOISE_FLOOR;
  float avgPower = NOISE_FLOOR;
  float minPowerOverThreshold = NOISE_FLOOR;
  bool energyOverThreshold = false;
};

typedef std::vector<SpectrumOverview> SpectrumOverviewVector;

class SignalOverview {
public:
  SignalOverview(){};
  virtual ~SignalOverview(){};
  SignalOverview &operator=(const SignalOverview &other);

  double widthHz = 0.0;
  double centerFreqHz = 0.0;
  float maxPower = NOISE_FLOOR;

  void print();
};

typedef std::vector<SignalOverview> SignalOverviewVector;

/*
 * Waterfall Data Class
 */

class WaterfallData {
public:
  // Waterfall will be fftSize wide by numRows long
  double centerFrequency;
  int fftSize;
  long numRows;

  float *data;

  virtual void reserve(int newFFTSize, long newNumRows);
  virtual void clear();
  virtual bool isEmpty();

  WaterfallData();
  virtual ~WaterfallData();
};

typedef std::vector<SignalOverview> SignalOverviewVector;

/*
 * EnergyAnalyzer class
 */
class EnergyAnalyzer {
protected:
  int fftSize;
  int centerBucket;
  float squelchThreshold;
  float minDutyCycle;

  FFT *fftProc;
  float *psdSpectrum;

public:
  // squelch threshold should be a number like -75.0
  // min duty cycle should be a fractional percentage (e.g. cycle = 0.1 for 10%)
  EnergyAnalyzer(int initFFTSize, float initSquelchThreshold,
                 float initMinDutyCycle, bool useWindow = true);
  virtual ~EnergyAnalyzer();

  inline void setThreshold(float newThreshold) {
    squelchThreshold = newThreshold;
  };
  inline float getThreshold() { return squelchThreshold; };
  inline void setDutyCycle(float newDutyCycle) { minDutyCycle = newDutyCycle; };
  inline float getDutyCycle() { return minDutyCycle; };

  inline float getFFTSize() { return fftSize; };
  inline FFT *getFFTProcessor() { return fftProc; };

  // Analyze chunks through frame and for each FFTSize block returns a
  // spectrumOverview object which describes duty cycle, max power, avg power,
  // etc. Think of it as an analysis/summary of each row in a waterfall plot.
  // return value is the number of samples analyzed (will be an integer multiple
  // of fftSize and the results of each FFTSize block
  virtual long analyze(const SComplex *frame, long numSamples,
                       SpectrumOverviewVector &results);

  // maxHold computes the max spectrum curve for the given frame.  Return value
  // is the number of samples processed.
  virtual long maxHold(const SComplex *frame, long numSamples,
                       FloatVector &maxSpectrum, bool useSquelch = true);

  // maxPower will look through a spectrum (something from maxHold or psd/rssi
  // from FFT and find the max value
  float maxPower(FloatVector &maxSpectrum);

  // AnalyzeSpectrum is a bit more generic.  It does rely on the local fftSize
  // and squelch threshold however, it analyzes just a single spectrum line and
  // returns params about it. Think of it as a 1-line analysis of what analyze
  // does in bulk. This can be useful in analyzing say a max spectrum line after
  // maxHold is called.
  void analyzeSpectrum(const float *spectrum, float &dutyCycle, float &maxPower,
                       float &minPower, float &centerAvgPower, float &avgPower);

  // findSignals looks through spectrum for signals meeting min/max/edge
  // requirements assumes spectrum is fftSize long sampleRate is needed to
  // determine Hz/bin If you don't know the center frequency for the radio, just
  // use 0.0 then the signal center frequency will be relative to center (0)
  // representing the frequncy shift off center Note it does use the
  // squelchThreshold value set in the constructor as well. edgeDBDown defines
  // required drop-off on edges (e.g. 15 dB for valid signal) Note: It's best to
  // feed a squelched spectrum (say from FFT::PowerSpectralDensity or from
  // maxHold) to this.
  int findSignals(const float *spectrum, double sampleRate,
                  double centerFrequencyHz, double minWidthHz,
                  double maxWidthHz, SignalOverviewVector &signalVector,
                  bool stopOnFirst);
  // float edgeDBDown, SignalOverviewVector& signalVector, bool stopOnFirst);

  // findSingleSignal makes some assumptions, first that the potential input
  // signal has been filtered to the band that may contain the signal.  Second,
  // that there may only be a single signal.  Therefore a simple boxing method
  // can be applied to determine any signal that may be present.  The minWidthHz
  // is used to ensure whatever is found at least meets the minimum.
  int findSingleSignal(const float *spectrum, double sampleRate,
                       double centerFrequencyHz, double minWidthHz,
                       SignalOverview &signalOverview);

  long getWaterfall(const SComplex *frame, long numSamples,
                    WaterfallData &waterfallData);

  // power binary slicer treats the spectrum/waterfall like OOK given the
  // squelch threshold and duty cycle that have been set.  The resulting output
  // is the bits (on or off) based on duty cycle for a frame rssi will be the
  // max power observed in any of the blocks that meet the duty cycle
  // requirement
  long powerBinarySlicer(const SComplex *frame, long numSamples,
                         FloatVector &bits, float &rssi);

  // Energy present uses the set squelch threshold and duty cycle to look for an
  // FFT block with the specified amount of energy present.  As soon as it finds
  // one, it returns. rssi will be the max power observed in any of the blocks
  // that meet the duty cycle requirement
  bool energyPresent(const SComplex *frame, long numSamples, float &rssi);

  // countEnergyBlocks uses the set squelch threshold and duty cycle to look for
  // FFT blocks with the specified amount of energy present.  This function
  // returns the count of those blocks rssi will be the max power observed in
  // any of the blocks that meet the duty cycle requirement
  long countEnergyBlocks(const SComplex *frame, long numSamples, float &rssi);
};

} // namespace MesaSignals

#endif /* LIB_SIGNALS_MESA_H_ */
