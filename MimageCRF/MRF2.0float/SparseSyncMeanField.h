/**
 * Class implementing sparse synchronous mean-field on a grid.
 *
 * Calculates sparsity at a node using the current approximate
 * marginal. Summands in the belief update are eliminated by
 * respecting this sparsity.
 *
 * Jerod Weinman
 * jerod@acm.org
 */

#ifndef __SPARSESYNCMEANFIELD_H__
#define __SPARSESYNCMEANFIELD_H__

#include "SyncMeanField.h"

#define FloatType float
#define FLOATTYPE float

//////////////////////
// Optimizer Arguments
//////////////////////

// Epsilon threshold for KL divergence
#define PARAM_DIVTOL 3

// Minimum number of non-sparse components
#define PARAM_MINBEAMS 4

class SparseSyncMeanField : public SyncMeanField {

 public: 

  SparseSyncMeanField(int width, int height, int nLabels, EnergyFunction *eng);
  SparseSyncMeanField(int nPixels, int nLabels, EnergyFunction *eng);
  ~SparseSyncMeanField();

  void setParameters(int , void *);

  int getNumNonSparseStates();	
  void printHistSparsity (int numBins);

 protected:

  void initializeSparseMeanField();
  void initializeAlg();
  void optimizeAlg(int nIterations);

  // Add the mean smoothness cost to a vector using the current SPARSE beliefs
  void addMeanSmoothnessCost(const int fixedPixel, const int varPixel, 
                             FloatType* dst );
  
  // Updates sparsities based on current approximate marginals
  void updateSparsity();
  void updateSparsity(int pixel);

  inline void logNormalizeSparse (FloatType* array, bool* sp )
    { ArrayMath::logNormalizeSparse (array, sp, m_nLabels); }

  void minusLogReverseCumulLogSum (FloatType* dst,FloatType* src);
  FloatType logSum(const FloatType a, const FloatType b);

  // Array of (pointers to) sparsity indicators -- FALSE when the label is EXCLUDED
  bool** m_nodeSparsity;
  
  // Array of (pointers to) normalized, exponentiated sparse beliefs
  FloatType** m_spExpNodeBeliefs;
  
  // Tolerance for compressing beliefs based on KL Divergence
  FloatType m_divTol;

  // Minimum number of beliefs to keep when compressing
  int m_minBeams;

#ifndef NDEBUG
  void printAll();
  void printArray(bool* array);
#endif

};
  
#endif
