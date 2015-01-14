/**
 * Class outlining sparse sum-product belief
 * propagation on a grid.
 *
 * Calculates sparsity at a node using the current approximate
 * marginal. Summands in the message update are eliminated by
 * respecting this sparsity.
 *
 * Jerod Weinman
 * jerod@acm.org
 *
 * $Id: SparseSumProd.h,v 1.1.1.1 2007/12/01 01:01:15 cpal Exp $
 */

#ifndef __SPARSESUMPROD_H__
#define __SPARSESUMPROD_H__

#include "SumProd.h"


//////////////////////
// Optimizer Arguments
//////////////////////

// Epsilon threshold for KL divergence
#define PARAM_DIVTOL 3

// Minimum number of non-sparse components
#define PARAM_MINBEAMS 4

class SparseSumProd : public virtual SumProd {

 public:
  
  SparseSumProd(int width, int height, int nLabels, EnergyFunction *eng);
  
  SparseSumProd(int nPixels, int nLabels,EnergyFunction *eng);

  ~SparseSumProd();


  void setParameters(int , void *);

  int getNumNonSparseStates();	
  int getMaxNumNonSparseStates();	
  void printHistSparsity(int numBins = 10);

  //#ifndef JDEBUG
  void printAll();
  void printArray(bool* array);
  //#endif
  
 protected:


  // Compute the sparsely updated message between two nodes
  void computeMessage(const int pixel1, const int pixel2,
                      NodeMessages* message1, NodeMessages* message2,
                      bool* sparsity, int direction);
  
  
  // Array of (pointers to) sparsity indicators -- FALSE when the label is EXCLUDED
  bool** m_nodeSparsity;

  // Tolerance for compressing beliefs based on KL Divergence
  FloatType m_divTol;

  // Minimum number of beliefs to keep when compressing
  int m_minBeams;

  
  //
  // Array arithmetic
  //

  /** Adds sparse vector src to accumulator dst */
  inline void addToSparse (FloatType* dst, FloatType* src,bool* sp) 
    { ArrayMath::addToSparse(dst, src, sp, m_nLabels); };

  /** Copy sparse vetor src to dst */
  inline void copyToSparse (FloatType* dst, FloatType* src, bool* sp)
    { ArrayMath::copyToSparse(dst, src, sp, m_nLabels); }

  /** Subtracts sparse vector src to accumulator dst */
  inline void subtractFromSparse (FloatType* dst, FloatType* src, bool* sp)
    { ArrayMath::subtractFromSparse(dst, src, sp, m_nLabels); };

  /** Sums a sparse vector stored in log-spce */
  inline FloatType logSumSparse (FloatType *array, bool* sp)
    { return ArrayMath::logSumSparse(array, sp, m_nLabels); };



private:
  
  // Common initialization routine
  void initializeSparseSumProd();

  // Calculate sparsity at each node based on current approximate marginals
  void updateSparsity();
  void updateSparsity(int pixel);

  void minusLogReverseCumulLogSum (FloatType* dst,FloatType* src);
};

#endif
