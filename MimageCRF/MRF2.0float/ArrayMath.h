#ifndef __ARRAYMATH_H__
#define __ARRAYMATH_H__

#define FloatType float

namespace ArrayMath {

  /** Copy vector src to dst */
  void copyTo(FloatType* dst, FloatType* src, const int len);

  /** Copy sparse vector src to dst */
  void copyToSparse(FloatType* dst, FloatType* src, bool* sp, const int len);

  /** Adds vector src to accumulator dst */
  void addTo(FloatType* dst, FloatType* src, const int len);

  /** Adds sparse vector src to accumulator dst */
  void addToSparse(FloatType* dst, FloatType* src, bool* sp, const int len);

  /** Zeros the elements of the array. */
  void zero(FloatType* array, const int len);

  /** Convex update of an array. */
  void convexAdd(const FloatType alpha, 
                 FloatType* dst, FloatType* src, const int len);

  /** Convex update of an array in log-space */
  void convexLogAdd(const FloatType alpha, 
                 FloatType* dst, FloatType* src, const int len);


  /** Subtracts vector src from dst */
  void subtractFrom(FloatType* dst, FloatType* src, const int len);

  /** Subtracts vector src from dst */
  void subtractFromSparse(FloatType* dst, FloatType* src, bool* sp, const int len);
  
  /** Subtracts a constant from dst */
  void subtractFrom(FloatType* dst, FloatType sth, const int len);

  /** Subtracts a constant from sparse dst */
  void subtractFromSparse(FloatType* dst, FloatType sth, bool* sp, const int len);
  
  /** Sums an array stored in log-space */
  FloatType logSum(FloatType* array, const int len);  

  /** Sums a spasre array stored in log-space */
  FloatType logSumSparse(FloatType* array, bool* sp, const int len);  

  /** Normalizes an array stored in log-space */
  void logNormalize(FloatType* array, const int len);

  /** Normalizes a sparse array stored in log-space */
  void logNormalizeSparse(FloatType* array, bool* sp, const int len);

  /** Finds the index of the largest element */
  int argMax(FloatType* array, const int len);
}

  
  

#endif
