/*
 * scomplex.h
 *
 *      Copyright 2017, Michael Piscopo
 *
 */

#ifndef LIB_SCOMPLEX_H_
#define LIB_SCOMPLEX_H_

#include <complex>
#include <vector>

// Common type definitions
typedef std::complex<float> SComplex;
typedef std::vector<std::complex<float>> ComplexVector;

// This structure version is used in some places for optimization to get direct
// access to the values without having to make function calls to get/set.
struct ComplexStruct {
  float real;
  float imag;
};

typedef struct ComplexStruct StructComplex;

#endif /* LIB_SCOMPLEX_H_ */
