#include <iostream>
#include <math.h>
#include <string.h>
#include "ArrayMath.h"


/** Copy vector src to dst */
void ArrayMath::copyTo(FloatType* dst, FloatType* src, const int len)
{
  memcpy(dst,src,len*sizeof(FloatType));
}

/** Copy sparse vector src to dst */
void ArrayMath::copyToSparse(FloatType* dst, FloatType* src, bool* sp, 
                             const int len)
{
  for (int i=0 ; i<len ; i++, dst++, src++, sp++)
    if (*sp)
      *dst = *src;
}


/** Zeros the elements of the array. */
void ArrayMath::zero(FloatType* array, const int len)
{
  for (int i=0 ; i<len ; i++, array++)
    *array = 0;
}

/** Adds vector src to accumulator dst */
void ArrayMath::addTo(FloatType* dst, FloatType* src, const int len)
{
  for (int i=0 ; i<len ; i++, dst++,src++)
    *dst += *src;
}

/** Adds sparse vector src to accumulator dst */
void ArrayMath::addToSparse(FloatType* dst, FloatType* src, bool* sp, 
                            const int len)
{
  for (int i=0 ; i<len ; i++, dst++,src++,sp++)
    if (*sp)
      *dst += *src;
}

/** Convex update of an array.
 * 
 * Parameter alpha is the fraction of the new array to use in the update
 */
void ArrayMath::convexAdd(const FloatType alpha, 
               FloatType* dst, FloatType* src, const int len)
{
    
  FloatType calpha = 1-alpha; // Complement of alpha
    
  for (int i=0 ; i<len ; i++, dst++,src++)
    *dst = calpha * (*dst) + alpha * (*src);
    
}

/** Convex update of an array in log-space.
 * 
 * Parameter alpha is the fraction of the new array to use in the update
 */
void ArrayMath::convexLogAdd(const FloatType alpha, 
               FloatType* dst, FloatType* src, const int len)
{
    
  FloatType calpha = 1-alpha; // Complement of alpha
    
  for (int i=0 ; i<len ; i++, dst++,src++)
    if (*dst > *src)
      *dst = log( calpha + alpha * exp(*src-*dst) ) + *dst;
    else
      *dst = log( calpha * exp(*dst-*src) + alpha) + *src;
}

/** Subtracts vector src from dst */
void ArrayMath::subtractFrom(FloatType* dst, FloatType* src, const int len)
{
  for (int i=0 ; i<len ; i++, dst++,src++)
    *dst -= *src;
}

/** Subtracts sparse vector src from dst */
void ArrayMath::subtractFromSparse(FloatType* dst, FloatType* src, bool* sp, 
                                   const int len)
{
  for (int i=0 ; i<len ; i++, dst++,src++,sp++)
    if (*sp)
      *dst -= *src;
}

/** Subtracts a constant from dst */
void ArrayMath::subtractFrom(FloatType* dst, FloatType sth, const int len)
{
  for (int i=0 ; i<len ; i++, dst++)
    *dst -= sth; // Subtract the subtrahend!
}

/** Subtracts a constant from dst */
void ArrayMath::subtractFromSparse(FloatType* dst, FloatType sth, bool* sp, 
                                   const int len)
{
  for (int i=0 ; i<len ; i++, dst++, sp++)
    if (*sp)
      *dst -= sth; // Subtract the subtrahend!
}
  
/** Sums an array stored in log-space */
FloatType ArrayMath::logSum(FloatType* array, const int len) {
    
  register FloatType sum = 0;

  // Find maximum
  register FloatType mx = *array;
  FloatType* pa = array+1;

  for (int i=1; i<len ; i++, pa++)
    if (*pa>mx)
      mx = *pa;
    
  // Cumulatively add, removing the maximum from the exponent
  for (int i=0 ; i<len ; i++, array++)
    sum += exp(*array-mx);

  // std::cout <<" sum=" << sum <<" mx=" << mx << " len=" <<len <<std::endl;

  return log(sum)+mx;
}

/** Sums a sparse array stored in log-space */
FloatType ArrayMath::logSumSparse(FloatType* array, bool* sp, const int len) {

  register FloatType sum = 0;

  // Find maximum
  register FloatType mx = 0;
  FloatType* pa = array;
  bool* ps = sp;
  bool fmx = false;
  for (int i=0; i<len ; i++, pa++,ps++)
    if (*ps && ((fmx && *pa>mx) || !fmx))
      {
        mx = *pa;
        fmx = true;
      }
  
  // Cumulatively add, removing the maximum from the exponent
  for (int i=0 ; i<len ; i++, array++,sp++)
    if (*sp)
      sum += exp(*array-mx);
    
  return log(sum)+mx;
}
  
/** Normalizes an array stored in log-space */
void ArrayMath::logNormalize(FloatType* array, const int len)
{
  FloatType logZ = logSum(array, len);
  subtractFrom(array,logZ, len);
}

/** Normalizes a sparse array stored in log-space */
void ArrayMath::logNormalizeSparse(FloatType* array, bool* sp, const int len)
{
  FloatType logZ = logSumSparse(array, sp, len);
  subtractFromSparse(array,logZ, sp, len);
}

/** Finds the index of the largest element */
int ArrayMath::argMax(FloatType* array, const int len)
{

  int amx = 0;
  FloatType mx = *array;
  
  array++;

  for (int i=1 ; i<len ; i++, array++)
    if (mx < *array)
      {
        amx = i;
        mx = *array;
      }
  
  return amx;
}
