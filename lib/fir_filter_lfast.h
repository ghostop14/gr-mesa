/*
 * fir_filter_lfast.h
 *
 *      Author: ghostop14
 */

#ifndef INCLUDE_LFAST_FIR_FILTER_LFAST_H_
#define INCLUDE_LFAST_FIR_FILTER_LFAST_H_

#include <gnuradio/gr_complex.h>

#include <boost/thread/thread.hpp>
#include <volk/volk.h>
using namespace std;

namespace gr {
namespace lfast {
/*
 * Base Filter
 */
template <class io_type, class tap_type> class Filter {
protected:
  tap_type *alignedTaps;
  io_type *singlePointBuffer;
  std::vector<tap_type> d_taps;
  long numTaps;

public:
  Filter();
  Filter(const std::vector<tap_type> &newTaps);
  virtual ~Filter();

  virtual void setTaps(const std::vector<tap_type> &newTaps);
  // For compatibility
  inline virtual void set_taps(const std::vector<tap_type> &newTaps) {
    setTaps(newTaps);
  };
  virtual std::vector<tap_type> getTaps() const;
  inline virtual std::vector<tap_type> taps() const { return getTaps(); };
  inline virtual long ntaps() { return numTaps; };

  // Returns number of samples consumed / produced
  virtual long filter(io_type *outputBuffer, const io_type *inputBuffer,
                      long numSamples) {
    return 0;
  };
}; // end base filter template

/*
 * FIR Filter - Complex inputs, float taps
 */
class FIRFilterCCF : public Filter<gr_complex, float> {
public:
  FIRFilterCCF();
  FIRFilterCCF(const std::vector<float> &newTaps);
  virtual ~FIRFilterCCF();

  // NOTE: This routine is expecting numSamples to be an integer multiple of the
  // number of taps
  virtual long filterN(gr_complex *outputBuffer, const gr_complex *inputBuffer,
                       long numSamples);

  // OutputBuffer here should be at least numSamples/decimation in length
  virtual long filterNdec(gr_complex *outputBuffer,
                          const gr_complex *inputBuffer, long numSamples,
                          int decimation);

  virtual gr_complex filter(const gr_complex *inputBuffer);

  // This is for testing comparison.  This function does all calculations in 1
  // CPU-based thread.
  virtual long filterCPU(gr_complex *outputBuffer,
                         const gr_complex *inputBuffer, long numSamples);
};

/*
 * FIR Filter - float inputs, float taps
 */
class FIRFilterFFF : public Filter<float, float> {
public:
  FIRFilterFFF();
  FIRFilterFFF(const std::vector<float> &newTaps);
  virtual ~FIRFilterFFF();

  // NOTE: This routine is expecting numSamples to be an integer multiple of the
  // number of taps
  virtual long filterN(float *outputBuffer, const float *inputBuffer,
                       long numSamples);

  // OutputBuffer here should be at least numSamples/decimation in length
  virtual long filterNdec(float *outputBuffer, const float *inputBuffer,
                          long numSamples, int decimation);

  virtual gr_complex filter(const float *inputBuffer);

  // This is for testing comparison.  This function does all calculations in 1
  // CPU-based thread.
  virtual long filterCPU(float *outputBuffer, const float *inputBuffer,
                         long numSamples);
};

/*
 * FIR Filter - complex inputs, complex taps
 */
class FIRFilterCCC : public Filter<gr_complex, gr_complex> {
public:
  FIRFilterCCC();
  FIRFilterCCC(const std::vector<gr_complex> &newTaps);
  virtual ~FIRFilterCCC();

  // NOTE: This routine is expecting numSamples to be an integer multiple of the
  // number of taps
  virtual long filterN(gr_complex *outputBuffer, const gr_complex *inputBuffer,
                       long numSamples);

  // OutputBuffer here should be at least numSamples/decimation in length
  virtual long filterNdec(gr_complex *outputBuffer,
                          const gr_complex *inputBuffer, long numSamples,
                          int decimation);

  virtual gr_complex filter(const gr_complex *inputBuffer);

  // This is for testing comparison.  This function does all calculations in 1
  // CPU-based thread.
  virtual long filterCPU(gr_complex *outputBuffer,
                         const gr_complex *inputBuffer, long numSamples);
};

// -----------------------------------------------------------------
// ------  Multi-threaded filters ----------------------------------
// -----------------------------------------------------------------

/*
 * Multi-threaded base template.
 * Provides threading and buffers
 *
 */
template <class io_type> class MTBase {
protected:
  const io_type *pInputBuffer;
  io_type *pOutputBuffer;

  // Supports up to 16 threads
  boost::thread *threads[16];
  bool threadRunning[16];
  bool threadReady;
  int d_nthreads;
  int decimation;

public:
  MTBase(int nthreads = 4);
  virtual ~MTBase();

  int numThreads() { return d_nthreads; };

  inline virtual void setDecimation(int newDecimation) {
    decimation = newDecimation;
  };
  inline virtual int getDecimation() { return decimation; };
  inline virtual bool decimating() {
    if (decimation > 1)
      return true;
    else
      return false;
  };

  virtual long calcDecimationBlockSize(long numSamples);

  // NOTE: calcDecimationIndex assumes the index is an even multiple of a
  // decimation block size In other words, calcDecimationBlockSize was called to
  // get an appropriate block size, then blockStartIndex is an integer multiple
  // of that result.
  virtual long calcDecimationIndex(long blockStartIndex);

  // NOTE: This method IS NOT thread-safe.  Make sure to use a scoped_lock or be
  // sure no threads are running before changing the number of threads.
  virtual void setThreads(int nthreads);

  virtual bool anyThreadRunning();
};

// --------------------------------------------------
// Multi-threaded filter, complex data, float taps
// --------------------------------------------------
class FIRFilterCCF_MT : public MTBase<gr_complex>, public FIRFilterCCF {
protected:
  virtual void runThread1(long startIndex, long numSamples);
  virtual void runThread2(long startIndex, long numSamples);
  virtual void runThread3(long startIndex, long numSamples);
  virtual void runThread4(long startIndex, long numSamples);
  virtual void runThread5(long startIndex, long numSamples);
  virtual void runThread6(long startIndex, long numSamples);
  virtual void runThread7(long startIndex, long numSamples);
  virtual void runThread8(long startIndex, long numSamples);
  virtual void runThread9(long startIndex, long numSamples);
  virtual void runThread10(long startIndex, long numSamples);
  virtual void runThread11(long startIndex, long numSamples);
  virtual void runThread12(long startIndex, long numSamples);
  virtual void runThread13(long startIndex, long numSamples);
  virtual void runThread14(long startIndex, long numSamples);
  virtual void runThread15(long startIndex, long numSamples);
  virtual void runThread16(long startIndex, long numSamples);

public:
  FIRFilterCCF_MT(int nthreads);
  FIRFilterCCF_MT(const std::vector<float> &newTaps, int nthreads);
  virtual ~FIRFilterCCF_MT();

  // NOTE: This routine is expecting numSamples to be an integer multiple of the
  // number of taps
  virtual long filterN(gr_complex *outputBuffer, const gr_complex *inputBuffer,
                       long numSamples);

  // OutputBuffer here should be at least numSamples/decimation in length
  virtual long filterNdec(gr_complex *outputBuffer,
                          const gr_complex *inputBuffer, long numSamples,
                          int decimation);
};

// --------------------------------------------------
// Multi-threaded filter, complex data, float taps
// --------------------------------------------------
class FIRFilterFFF_MT : public MTBase<float>, public FIRFilterFFF {
protected:
  virtual void runThread1(long startIndex, long numSamples);
  virtual void runThread2(long startIndex, long numSamples);
  virtual void runThread3(long startIndex, long numSamples);
  virtual void runThread4(long startIndex, long numSamples);
  virtual void runThread5(long startIndex, long numSamples);
  virtual void runThread6(long startIndex, long numSamples);
  virtual void runThread7(long startIndex, long numSamples);
  virtual void runThread8(long startIndex, long numSamples);
  virtual void runThread9(long startIndex, long numSamples);
  virtual void runThread10(long startIndex, long numSamples);
  virtual void runThread11(long startIndex, long numSamples);
  virtual void runThread12(long startIndex, long numSamples);
  virtual void runThread13(long startIndex, long numSamples);
  virtual void runThread14(long startIndex, long numSamples);
  virtual void runThread15(long startIndex, long numSamples);
  virtual void runThread16(long startIndex, long numSamples);

public:
  FIRFilterFFF_MT(int nthreads);
  FIRFilterFFF_MT(const std::vector<float> &newTaps, int nthreads);
  virtual ~FIRFilterFFF_MT();

  // NOTE: This routine is expecting numSamples to be an integer multiple of the
  // number of taps
  virtual long filterN(float *outputBuffer, const float *inputBuffer,
                       long numSamples);

  // OutputBuffer here should be at least numSamples/decimation in length
  virtual long filterNdec(float *outputBuffer, const float *inputBuffer,
                          long numSamples, int decimation);
};

// --------------------------------------------------
// Multi-threaded filter, complex data, complex taps
// --------------------------------------------------
class FIRFilterCCC_MT : public MTBase<gr_complex>, public FIRFilterCCC {
protected:
  virtual void runThread1(long startIndex, long numSamples);
  virtual void runThread2(long startIndex, long numSamples);
  virtual void runThread3(long startIndex, long numSamples);
  virtual void runThread4(long startIndex, long numSamples);
  virtual void runThread5(long startIndex, long numSamples);
  virtual void runThread6(long startIndex, long numSamples);
  virtual void runThread7(long startIndex, long numSamples);
  virtual void runThread8(long startIndex, long numSamples);
  virtual void runThread9(long startIndex, long numSamples);
  virtual void runThread10(long startIndex, long numSamples);
  virtual void runThread11(long startIndex, long numSamples);
  virtual void runThread12(long startIndex, long numSamples);
  virtual void runThread13(long startIndex, long numSamples);
  virtual void runThread14(long startIndex, long numSamples);
  virtual void runThread15(long startIndex, long numSamples);
  virtual void runThread16(long startIndex, long numSamples);

public:
  FIRFilterCCC_MT(int nthreads);
  FIRFilterCCC_MT(const std::vector<gr_complex> &newTaps, int nthreads);
  virtual ~FIRFilterCCC_MT();

  // NOTE: This routine is expecting numSamples to be an integer multiple of the
  // number of taps
  virtual long filterN(gr_complex *outputBuffer, const gr_complex *inputBuffer,
                       long numSamples);

  // OutputBuffer here should be at least numSamples/decimation in length
  virtual long filterNdec(gr_complex *outputBuffer,
                          const gr_complex *inputBuffer, long numSamples,
                          int decimation);
};
} // namespace lfast
} // namespace gr

#endif /* INCLUDE_LFAST_FIR_FILTER_LFAST_H_ */
