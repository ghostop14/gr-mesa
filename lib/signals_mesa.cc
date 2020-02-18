/*
 * signals_mesa.cc
 *
 *      Copyright 2019, Michael Piscopo
 *
 */

#include "signals_mesa.h"

#include <boost/math/special_functions/round.hpp>
#include <cassert>
#include <fftw3.h>
#include <gnuradio/fft/window.h>
#include <iostream> // std::reverse

#define PRINTDEBUG

namespace MesaSignals {

static int fftInstanceCount = 0;

void printArray(FloatVector &arr, string name) {
  int arrSize = arr.size();

  if (name.length() > 0)
    std::cout << "[ Array '" << name << "' Size " << arrSize << "] ";
  else
    std::cout << "[Size " << arrSize << "] ";

  for (int i = 0; i < arrSize; i++) {
    if (i < (arrSize - 1))
      std::cout << arr[i] << ",";
    else
      std::cout << arr[i] << std::endl << std::endl;
  }
}

void printArray(float *arr, int arrSize, string name) {
  if (name.length() > 0)
    std::cout << "[ Array '" << name << "' Size " << arrSize << "] ";
  else
    std::cout << "[Size " << arrSize << "] ";

  for (int i = 0; i < arrSize; i++) {
    if (i < (arrSize - 1))
      std::cout << arr[i] << ",";
    else
      std::cout << arr[i] << std::endl << std::endl;
  }
}

// ------------------   FFT   ---------------------------------------
FFT::FFT(int FFTDirection, int initFFTSize, int initNumThreads)
    : fftDirection(FFTDirection), fftSize(initFFTSize),
      numThreads(initNumThreads) {
  init();
}

void FFT::init() {
  rssi_K_const =
      -10.0 * log10((float)fftSize) - 20.0 * log10(log2((float)fftSize));
  log2To10Factor = 10.0 / log2(10.0);

  fftPlan = NULL;
  fftLenBytes = fftSize * sizeof(float);
  halfFFTSizeBytes = fftLenBytes / 2;
  halfFFTSize = fftSize / 2;

  initThreads();
  importWisdom();

  // Requires FFT Size which is set in the constructor
  /*
  inputBuffer = new SComplex[inputBufferLength()];
  outputBuffer = new SComplex[outputBufferLength()];
  */
  // Aligned memory has better performance:
  size_t memAlignment = volk_get_alignment();
  inputBuffer = (SComplex *)volk_malloc(inputBufferLength() * sizeof(SComplex),
                                        memAlignment);
  outputBuffer = (SComplex *)volk_malloc(
      outputBufferLength() * sizeof(SComplex), memAlignment);
  tmpBuff =
      (float *)volk_malloc(outputBufferLength() * sizeof(float), memAlignment);

  alignedWindowTaps =
      (float *)volk_malloc(fftSize * sizeof(float), memAlignment);
  hasTaps = false;

  // fftDirection is set in the constructor

  // Some safety checks
  // Make sure the sizes are the same (they should be)
  assert(sizeof(fftwf_complex) == sizeof(SComplex));

  if (!inputBuffer)
    throw std::runtime_error("[FFT] input buffer allocation failed");

  if (!outputBuffer) {
    volk_free(inputBuffer);
    throw std::runtime_error("[FFT] output buffer allocation failed");
  }

  if (!tmpBuff) {
    volk_free(inputBuffer);
    volk_free(outputBuffer);
    throw std::runtime_error("[FFT] tmp buffer allocation failed");
  }

  fftPlan = fftwf_plan_dft_1d(
      fftSize, reinterpret_cast<fftwf_complex *>(inputBuffer),
      reinterpret_cast<fftwf_complex *>(outputBuffer),
      fftDirection, // fftDirection ? FFTW_FORWARD : FFTW_BACKWARD,
      FFTW_MEASURE);

  exportWisdom();
}

void FFT::initThreads() {
  if (fftInstanceCount == 0) {
    // Only need to initialize the thread system once
    fftwf_init_threads();
  }

  fftInstanceCount++;

  fftwf_plan_with_nthreads(numThreads);
}

void FFT::setWindow(int winType) {
  boost::mutex::scoped_lock scoped_lock(d_mutex);

  FloatVector windowTaps;
  size_t memAlignment = volk_get_alignment();

  switch (winType) {
  case WINDOWTYPE_NONE:
    hasTaps = false;
    break;
  case WINDOWTYPE_HAMMING:
    windowTaps = gr::fft::window::hamming(fftSize);
    memcpy(alignedWindowTaps, &windowTaps[0], fftSize * sizeof(float));
    hasTaps = true;
    break;

  case WINDOWTYPE_BLACKMAN_HARRIS:
    windowTaps = gr::fft::window::blackman_harris(fftSize);
    memcpy(alignedWindowTaps, &windowTaps[0], fftSize * sizeof(float));
    hasTaps = true;
    break;

  default:
    throw std::out_of_range("[FFT]: unknown window type.");
  }
}

void FFT::setWindow(FloatVector &newTaps) {
  boost::mutex::scoped_lock scoped_lock(d_mutex);

  if (newTaps.size() > 0) {
    hasTaps = true;

    if (newTaps.size() != fftSize)
      throw std::out_of_range("[FFT]: setWindow(newTaps) tap size " +
                              to_string(newTaps.size()) + " != fft size " +
                              to_string(fftSize));

    memcpy(alignedWindowTaps, &newTaps[0], fftSize * sizeof(float));

  } else
    hasTaps = false;
}

void FFT::clearWindow() { hasTaps = false; }

void FFT::importWisdom() {
  wisdomFilename = ".fftw_wisdom";

  FILE *pFile = fopen(wisdomFilename.c_str(), "r");
  if (pFile != NULL) {
    // File exists
    fclose(pFile);
    int retVal = fftwf_import_wisdom_from_filename(wisdomFilename.c_str());
  }
  /*
  else {
        std::cout << "[FFT] INFO: fftw wisdom file " << wisdomFilename << " does
  not exist.  We'll create it." << std::endl;
  }
  */
}

void FFT::exportWisdom() {
  int retVal = fftwf_export_wisdom_to_filename(wisdomFilename.c_str());
  /*
  if (retVal == 0) {
            std::cout << "[FFT] WARNING: Unable to export fftw wisdom file: " <<
  wisdomFilename << ".  Could be locked?" << std::endl;
  }
  */
}

inline void FFT::execute(bool shift) {
  boost::mutex::scoped_lock scoped_lock(d_mutex);

  if (hasTaps) {
    // Apply window function first
    volk_32fc_32f_multiply_32fc(inputBuffer, inputBuffer, alignedWindowTaps,
                                fftSize);
  }

  fftwf_execute((fftwf_plan)fftPlan);

  if (shift) {
    // DC center is at out[0] so need to swap halves first.
    // Move top half to tmp buffer
    memcpy(tmpBuff, &outputBuffer[halfFFTSize], halfFFTSizeBytes);
    // move lower half up
    memcpy(&outputBuffer[halfFFTSize], &outputBuffer[0], halfFFTSizeBytes);
    // put top half back in the lower half
    memcpy(&outputBuffer[0], &tmpBuff[0], halfFFTSizeBytes);
  }
}

inline void FFT::execute(float *pTaps, bool shift) {
  boost::mutex::scoped_lock scoped_lock(d_mutex);

  if (pTaps) {
    // Apply window function first
    volk_32fc_32f_multiply_32fc(inputBuffer, inputBuffer, pTaps, fftSize);
  }

  fftwf_execute((fftwf_plan)fftPlan);

  if (shift) {
    // DC center is at out[0] so need to swap halves first.
    // Move top half to tmp buffer
    memcpy(tmpBuff, &outputBuffer[halfFFTSize], halfFFTSizeBytes);
    // move lower half up
    memcpy(&outputBuffer[halfFFTSize], &outputBuffer[0], halfFFTSizeBytes);
    // put top half back in the lower half
    memcpy(&outputBuffer[0], &tmpBuff[0], halfFFTSizeBytes);
  }
}

inline void FFT::PowerSpectralDensity(float *psdBuffer, float squelchThreshold,
                                      float onSquelchSetRSSI) {
  return rssi(psdBuffer, squelchThreshold, onSquelchSetRSSI);
}

void FFT::rssi(float *psdBuffer, float squelchThreshold,
               float onSquelchSetRSSI) {
  // Note: using aligned memory can be notably faster than the unaligned
  // versions Calcs were slightly different using the 3-call approach.  Not sure
  // why. The math looks the same.
  /*
  // Complex to mag^2
  volk_32fc_magnitude_squared_32f(psdBuffer,outputBuffer,fftSize);

  // 10 Log10() - (10*math.log10(fftsize) - 20*math.log10(math.log(fftsize,2)))
  // See
  https://lists.gnu.org/archive/html/discuss-gnuradio/2015-09/msg00350.html for
  ref...
  // Calc 10 log10(x) as 10*log2(x)/log2(10) = (10/log2(10)) * log2(x)
  volk_32f_log2_32f(psdBuffer,psdBuffer,fftSize);

  // Finish the 10Log(x) calculation
  // Store results in place
  volk_32f_s32f_multiply_32f(psdBuffer,psdBuffer,log2To10Factor,fftSize);
  */

  volk_32fc_s32f_x2_power_spectral_density_32f(psdBuffer, outputBuffer, fftSize,
                                               1.0, fftSize);

  // DC center is at out[0] so need to swap halves first.
  // Move top half to tmp buffer
  memcpy(tmpBuff, &psdBuffer[halfFFTSize], halfFFTSizeBytes);
  // move lower half up
  memcpy(&psdBuffer[halfFFTSize], &psdBuffer[0], halfFFTSizeBytes);
  // put top half back in the lower half
  memcpy(&psdBuffer[0], &tmpBuff[0], halfFFTSizeBytes);

  if (squelchThreshold != SQUELCH_DISABLE) {
    for (int j = 0; j < fftSize; j++) {
      if (psdBuffer[j] <= squelchThreshold) // Squelch noise in the spectrum
        psdBuffer[j] = onSquelchSetRSSI;    // -100.0;
    }
  }
}

FFT::~FFT() {
  fftwf_destroy_plan((fftwf_plan)fftPlan);

  /*
      delete[] inputBuffer;
      delete[] outputBuffer;
      */
  volk_free(inputBuffer);
  volk_free(outputBuffer);
  volk_free(tmpBuff);

  if (alignedWindowTaps) {
    // In any case we need to clear what we had.
    volk_free(alignedWindowTaps);
    alignedWindowTaps = NULL;
  }

  fftInstanceCount--;

  if (fftInstanceCount == 0) {
    fftwf_cleanup_threads();
  }
}

// -----------------  End FFT ---------------------------------------

// -----------------  Start SpectrumOverview
// ---------------------------------------
SpectrumOverview &SpectrumOverview::operator=(const SpectrumOverview &other) {
  if (this != &other) { // self-assignment check expected
    dutyCycle = other.dutyCycle;
    maxPower = other.maxPower;
    energyOverThreshold = other.energyOverThreshold;
    avgPower = other.avgPower;
    minPowerOverThreshold = other.minPowerOverThreshold;
  }
  return *this;
}

// -----------------  End SpectrumOverview
// ---------------------------------------

// -----------------  Start Waterfall Data
// ---------------------------------------
WaterfallData::WaterfallData() {
  data = NULL;
  fftSize = 0;
  numRows = 0;
  centerFrequency = 0.0;
}

void WaterfallData::reserve(int newFFTSize, long newNumRows) {
  if ((fftSize * numRows) < (newFFTSize * newNumRows)) {
    if (data) {
      volk_free(data);
      data = NULL;
    }

    fftSize = newFFTSize;
    numRows = newNumRows;

    size_t memAlignment = volk_get_alignment();
    data =
        (float *)volk_malloc(fftSize * numRows * sizeof(float), memAlignment);
    memset(data, 0x00, numRows * fftSize * sizeof(float));
  }
}

WaterfallData::~WaterfallData() {
  if (data) {
    volk_free(data);
    data = NULL;
  }
}

void WaterfallData::clear() {
  if (data) {
    memset(data, 0x00, numRows * fftSize * sizeof(float));
  }
}

bool WaterfallData::isEmpty() {
  if (data) // if we have data, we're not empty
    return false;
  else
    return true;
}

// -----------------  End Waterfall Data ---------------------------------------

// -----------------  Start SignalOverview
// ---------------------------------------
SignalOverview &SignalOverview::operator=(const SignalOverview &other) {
  if (this != &other) { // self-assignment check expected
    widthHz = other.widthHz;
    maxPower = other.maxPower;
    centerFreqHz = other.centerFreqHz;
  }
  return *this;
}

void SignalOverview::print() {
  std::cout << "Center Frequency: " << centerFreqHz << std::endl;

  std::cout << "Width: " << widthHz << std::endl;
  std::cout << "Max Power: " << maxPower << std::endl;
}

// -----------------  End SignalOverview ---------------------------------------

// -----------------  Start Energy Analyzer
// ---------------------------------------
EnergyAnalyzer::EnergyAnalyzer(int initFFTSize, float initSquelchThreshold,
                               float initMinDutyCycle, bool useWindow) {
  fftSize = initFFTSize;
  centerBucket = fftSize / 2;

  squelchThreshold = initSquelchThreshold;
  minDutyCycle = initMinDutyCycle;

  fftProc = new FFT(FFTDIRECTION_FORWARD, fftSize);

  if (useWindow)
    fftProc->setWindow(WINDOWTYPE_BLACKMAN_HARRIS);

  size_t memAlignment = volk_get_alignment();
  psdSpectrum = (float *)volk_malloc(fftSize * sizeof(float), memAlignment);
}

EnergyAnalyzer::~EnergyAnalyzer() {
  delete fftProc;

  volk_free(psdSpectrum);
}

long EnergyAnalyzer::analyze(const SComplex *frame, long numSamples,
                             SpectrumOverviewVector &results) {
  long numBlocks = numSamples / fftSize;

  if (numBlocks <= 0 || (frame == NULL))
    return 0;

  SComplex *fftInput;
  long index;

  results.clear();
  // Allocate the memory all at one time.
  results.reserve(numBlocks);

  for (long i = 0; i < numBlocks; i++) {
    // Calculate the FFT for the current block
    index = i * fftSize;
    fftInput = fftProc->getInputBuffer();
    memcpy(fftInput, &frame[index], fftSize * sizeof(SComplex));
    fftProc->execute();

    // Get the PSD of the result with a squelch threshold
    fftProc->PowerSpectralDensity(psdSpectrum, squelchThreshold);

    // Now analyze the spectrum
    int bucketswithPower = 0;
    float maxPower = NOISE_FLOOR;
    float totalPower = 0.0;
    float minPower = 1000.0;

    for (int j = 0; j < fftSize; j++) {
      if (psdSpectrum[j] >= squelchThreshold) {
        bucketswithPower++;
      }

      if (psdSpectrum[j] < minPower) {
        minPower = psdSpectrum[j];
      }

      totalPower += psdSpectrum[j];

      if (psdSpectrum[j] > maxPower)
        maxPower = psdSpectrum[j];
    }

    if (minPower == 1000.0)
      minPower = NOISE_FLOOR;

    float curDutyCycle = (float)bucketswithPower / (float)fftSize;

    SpectrumOverview spectrumOverview;
    spectrumOverview.dutyCycle = curDutyCycle;
    spectrumOverview.maxPower = maxPower;
    spectrumOverview.minPower = minPower;
    spectrumOverview.centerAvgPower =
        (psdSpectrum[centerBucket] + psdSpectrum[centerBucket + 1]) / 2.0;
    spectrumOverview.avgPower = totalPower / (float)fftSize;
    spectrumOverview.minPowerOverThreshold = minPower;

    results.push_back(spectrumOverview);
  }

  return (numBlocks * fftSize);
}

void EnergyAnalyzer::analyzeSpectrum(const float *spectrum, float &dutyCycle,
                                     float &maxPower, float &minPower,
                                     float &centerAvgPower, float &avgPower) {
  int bucketswithPower = 0;
  maxPower = NOISE_FLOOR;
  float totalPower = 0.0;
  minPower = 1000.0;

  for (int j = 0; j < fftSize; j++) {
    if (spectrum[j] >= squelchThreshold) {
      bucketswithPower++;
    }

    if (spectrum[j] < minPower)
      minPower = spectrum[j];

    totalPower += spectrum[j];

    if (spectrum[j] > maxPower)
      maxPower = spectrum[j];
  }

  dutyCycle = (float)bucketswithPower / (float)fftSize;

  centerAvgPower = (spectrum[centerBucket] + spectrum[centerBucket + 1]) / 2.0;

  avgPower = totalPower / (float)fftSize;
}

long EnergyAnalyzer::getWaterfall(const SComplex *frame, long numSamples,
                                  WaterfallData &waterfallData) {
  long numBlocks = numSamples / fftSize;

  if (numBlocks <= 0 || (frame == NULL))
    return 0;

  SComplex *fftInput;
  long index;

  for (long i = 0; i < numBlocks; i++) {
    // Calculate the FFT for the current block
    index = i * fftSize;
    fftInput = fftProc->getInputBuffer();
    memcpy(fftInput, &frame[index], fftSize * sizeof(SComplex));
    fftProc->execute();

    // Get the PSD of the result with a squelch threshold
    fftProc->PowerSpectralDensity(psdSpectrum, squelchThreshold);

    memcpy(&waterfallData.data[index], psdSpectrum, fftSize * sizeof(float));
  }

  return (numBlocks * fftSize);
}

long EnergyAnalyzer::powerBinarySlicer(const SComplex *frame, long numSamples,
                                       FloatVector &bits, float &rssi) {
  long numBlocks = numSamples / fftSize;

  rssi = NOISE_FLOOR;

  if (numBlocks <= 0 || (frame == NULL))
    return 0;

  SComplex *fftInput;
  long index;

  if (bits.size() != numBlocks) {
    bits.resize(numBlocks);
  }

  float spectrumDutyCycle;
  float *pBit;
  float maxPower;
  int bucketswithPower;
  int j;

  pBit = &bits[0];

  for (long i = 0; i < numBlocks; i++) {
    // Calculate the FFT for the current block
    index = i * fftSize;
    fftInput = fftProc->getInputBuffer();
    memcpy(fftInput, &frame[index], fftSize * sizeof(SComplex));
    fftProc->execute();

    // Get the PSD of the result with a squelch threshold
    fftProc->PowerSpectralDensity(psdSpectrum, squelchThreshold);

    bucketswithPower = 0;
    maxPower = NOISE_FLOOR;

    for (j = 0; j < fftSize; j++) {
      if (psdSpectrum[j] >= squelchThreshold) {
        bucketswithPower++;

        if (psdSpectrum[j] > maxPower)
          maxPower = psdSpectrum[j];
      }
    }

    spectrumDutyCycle = (float)bucketswithPower / (float)fftSize;

    if (spectrumDutyCycle >= minDutyCycle) {
      *pBit++ = 1;

      if (maxPower > rssi)
        rssi = maxPower;
    } else {
      *pBit++ = 0;
    }
  }

  return (numBlocks * fftSize);
}

bool EnergyAnalyzer::energyPresent(const SComplex *frame, long numSamples,
                                   float &rssi) {
  long numBlocks = numSamples / fftSize;

  if (numBlocks <= 0 || (frame == NULL))
    return 0;

  SComplex *fftInput;
  long index;
  rssi = NOISE_FLOOR;

  for (long i = 0; i < numBlocks; i++) {
    // Calculate the FFT for the current block
    index = i * fftSize;
    fftInput = fftProc->getInputBuffer();
    memcpy(fftInput, &frame[index], fftSize * sizeof(SComplex));
    fftProc->execute();

    // Get the PSD of the result with a squelch threshold
    fftProc->PowerSpectralDensity(psdSpectrum, squelchThreshold);

    // Now analyze the spectrum
    int bucketswithPower = 0;

    for (int j = 0; j < fftSize; j++) {
      if (psdSpectrum[j] >= squelchThreshold) {
        bucketswithPower++;
      }

      if (psdSpectrum[j] > rssi)
        rssi = psdSpectrum[j];
    }

    float curDutyCycle = (float)bucketswithPower / (float)fftSize;

    if (curDutyCycle >= minDutyCycle)
      return true; // return immediately
  }

  return false;
}

long EnergyAnalyzer::countEnergyBlocks(const SComplex *frame, long numSamples,
                                       float &rssi) {
  long numBlocks = numSamples / fftSize;

  if (numBlocks <= 0 || (frame == NULL))
    return 0;

  SComplex *fftInput;
  long index;
  long numEnergyBlocks = 0;
  rssi = NOISE_FLOOR;
  float maxPower;

  for (long i = 0; i < numBlocks; i++) {
    // Calculate the FFT for the current block
    index = i * fftSize;
    fftInput = fftProc->getInputBuffer();
    memcpy(fftInput, &frame[index], fftSize * sizeof(SComplex));
    fftProc->execute();

    // Get the PSD of the result with a squelch threshold
    fftProc->PowerSpectralDensity(psdSpectrum, squelchThreshold);

    // Now analyze the spectrum
    int bucketswithPower = 0;
    maxPower = NOISE_FLOOR;

    for (int j = 0; j < fftSize; j++) {
      if (psdSpectrum[j] >= squelchThreshold) {
        bucketswithPower++;
      }
      if (psdSpectrum[j] > maxPower)
        maxPower = psdSpectrum[j];
    }

    float curDutyCycle = (float)bucketswithPower / (float)fftSize;

    if (curDutyCycle >= minDutyCycle) {
      if (maxPower > rssi)
        rssi = maxPower;

      numEnergyBlocks++;
    }
  }

  return numEnergyBlocks;
}

long EnergyAnalyzer::maxHold(const SComplex *frame, long numSamples,
                             FloatVector &maxSpectrum, bool useSquelch) {
  long numBlocks = numSamples / fftSize;

  if (numBlocks <= 0 || (frame == NULL))
    return 0;

  SComplex *fftInput;
  long index;

  if (maxSpectrum.size() != fftSize) {
    maxSpectrum.resize(fftSize);
  }

  // Set default to NOISE_FLOOR
  for (int i = 0; i < fftSize; i++) {
    maxSpectrum[i] = NOISE_FLOOR;
  }

  for (long i = 0; i < numBlocks; i++) {
    // Calculate the FFT for the current block
    index = i * fftSize;
    fftInput = fftProc->getInputBuffer();
    memcpy(fftInput, &frame[index], fftSize * sizeof(SComplex));
    fftProc->execute();

    // Get the PSD of the result with a squelch threshold
    if (useSquelch)
      fftProc->PowerSpectralDensity(psdSpectrum, squelchThreshold);
    else
      fftProc->PowerSpectralDensity(psdSpectrum, SQUELCH_DISABLE);

    for (int j = 0; j < fftSize; j++) {
      if (psdSpectrum[j] >= maxSpectrum[j]) {
        maxSpectrum[j] = psdSpectrum[j];
      }
    }
  }

  return (numBlocks * fftSize);
}

float EnergyAnalyzer::maxPower(FloatVector &maxSpectrum) {
  unsigned int index;

  volk_32f_index_max_32u(&index, &maxSpectrum[0], maxSpectrum.size());

  return maxSpectrum[index];
}

int EnergyAnalyzer::findSingleSignal(const float *spectrum, double sampleRate,
                                     double centerFrequencyHz,
                                     double minWidthHz,
                                     SignalOverview &signalOverview) {
  if (spectrum == NULL)
    return 0;

  double hzPerBucket = sampleRate / (float)fftSize;
  double minFrequency = centerFrequencyHz - (sampleRate / 2.0);
  float maxPower = NOISE_FLOOR;

  int fftStart = 0;
  int fftEnd = fftSize - 1;

  // Let's eliminate all of the squelched spectrum at the ends for the main
  // loop. NOTE: This shouldn't mean more FOR looping, it just breaks it into 3
  // chunks. Low side
  bool foundSignal = false;

  for (int i = 0; i < fftSize; i++) {
    if (spectrum[i] > squelchThreshold) {
      fftStart = i;
      foundSignal = true;
      break;
    }
  }

  if (!foundSignal)
    return 0;

  // High side
  for (int i = fftSize - 1; i >= 0; i--) {
    if (spectrum[i] > squelchThreshold) {
      fftEnd = i;

      break;
    }
  }

  for (int i = fftStart; i < (fftEnd + 1);
       i++) { // going -1 here since we're looking for the trailing edge.
    if (spectrum[i] > maxPower)
      maxPower = spectrum[i];
  }

  double widthHz = (float)(fftEnd - fftStart + 1) * hzPerBucket;

  if ((widthHz >= minWidthHz)) {
    int centerBin = fftStart + (fftEnd - fftStart) / 2;
    signalOverview.widthHz = widthHz;
    signalOverview.centerFreqHz = minFrequency + (float)centerBin * hzPerBucket;
    signalOverview.maxPower = maxPower;

    return 1;
  } else {
    return 0;
  }
}

int EnergyAnalyzer::findSignals(const float *spectrum, double sampleRate,
                                double centerFrequencyHz, double minWidthHz,
                                double maxWidthHz,
                                SignalOverviewVector &signalVector,
                                bool stopOnFirst) {
  //		float edgeDBDown,SignalOverviewVector& signalVector, bool
  //stopOnFirst) {

  if (spectrum == NULL)
    return 0;

  double hzPerBucket = sampleRate / (float)fftSize;
  double minFrequency = centerFrequencyHz - (sampleRate / 2.0);

  signalVector.clear();

  float bucketBefore =
      spectrum[0]; // one bucket back.  For first pass if we set it to noise
                   // floor, it'll trigger too frequently.
  int startBucket = -1;
  int endBucket = -1;
  float maxPower = NOISE_FLOOR;
  bool inSignal = false;

  int fftStart = 0;
  int fftEnd = fftSize - 1;

  // Let's eliminate all of the squelched spectrum at the ends for the main
  // loop. NOTE: This shouldn't mean more FOR looping, it just breaks it into 3
  // chunks. Low side
  bool foundSignal = false;

  for (int i = 0; i < fftSize; i++) {
    if (spectrum[i] > squelchThreshold) {
      fftStart = i;
      foundSignal = true;
      break;
    }
  }

  if (!foundSignal)
    return 0;

  // High side
  for (int i = fftSize - 1; i >= 0; i--) {
    if (spectrum[i] > squelchThreshold) {
      fftEnd = i;

      break;
    }
  }

  // backoff for the lookahead of i+1 below.
  if (fftEnd == (fftSize - 1))
    fftEnd--;

  // Main spectrum
  for (int i = fftStart; i < (fftEnd + 1);
       i++) { // going -1 here since we're looking for the trailing edge.
    if (inSignal) {
      if (spectrum[i] > maxPower)
        maxPower = spectrum[i];

      // We've already found a rising edge, now we're looking for the end
      if (spectrum[i] <= squelchThreshold &&
          (spectrum[i + 1] <= squelchThreshold)) {
        // looks like we found the far edge.
        endBucket = i;

        double widthHz = (double)(endBucket - startBucket + 1) * hzPerBucket;

        if ((widthHz >= minWidthHz) && (widthHz <= maxWidthHz)) {
          if (maxPower > squelchThreshold) {
            // Signal descriptor looks good.  Push back.
            SignalOverview signalOverview;

            int centerBin = startBucket + (endBucket - startBucket) / 2;
            signalOverview.widthHz = widthHz;
            signalOverview.centerFreqHz =
                minFrequency + (float)centerBin * hzPerBucket;
            signalOverview.maxPower = maxPower;

            signalVector.push_back(signalOverview);

            if (stopOnFirst)
              return signalVector.size();
          }
        }

        // Reset some params for next signal
        maxPower = NOISE_FLOOR;
        inSignal = false;
      }
    } else {
      // We're looking for the rising edge
      // NOTE: Taking squelch into account here can cause issues.  It can cause
      // an unintended sharp edge.
      if (spectrum[i] >
          squelchThreshold) { // (spectrum[i]-bucketBefore) >= edgeDBDown) {
        // found a leading edge.
        inSignal = true;
        startBucket = endBucket = i;
      }
    }

    bucketBefore = spectrum[i];
  }

  if (inSignal) {
    // Signal was still continuing at the high edge of the spectrum.  So test
    // last one.
    endBucket = fftSize - 1; // set to last bucket

    double widthHz = (double)(endBucket - startBucket) * hzPerBucket;

    // Just check if what we have so far is wide enough
    if (widthHz >= minWidthHz) {
      if (maxPower > squelchThreshold) {
        // Signal descriptor looks good.  Push back.
        SignalOverview signalOverview;

        int centerBin = (endBucket - startBucket) / 2;
        signalOverview.widthHz = widthHz;
        signalOverview.centerFreqHz =
            minFrequency + (float)centerBin * hzPerBucket;
        signalOverview.maxPower = maxPower;

        signalVector.push_back(signalOverview);
      }
    }
  }

  return signalVector.size();
}

} // namespace MesaSignals

// -----------------  End Energy Analyzer
// ---------------------------------------
