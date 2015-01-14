/**
 * Class implementing sparse mean-field on a grid.
 *
 * Calculates sparsity at a node using the current approximate
 * marginal. Summands in the belief update are eliminated by
 * respecting this sparsity.
 *
 * Jerod Weinman
 * jerod@acm.org
 */


#ifndef __SPARSEMEANFIELD_H__
#define __SPARSEMEANFIELD_H__

#include "MeanField.h"

#define FloatType float
#define FLOATTYPE float

//////////////////////
// Optimizer Arguments
//////////////////////

// Epsilon threshold for KL divergence
#define PARAM_DIVTOL 3

// Minimum number of non-sparse components
#define PARAM_MINBEAMS 4

class SparseMeanField : public MeanField {


 public: 

  SparseMeanField(int width, int height, int nLabels, EnergyFunction *eng);
  SparseMeanField(int nPixels, int nLabels, EnergyFunction *eng);
  ~SparseMeanField();

  void setParameters(int , void *);
  
  // Miscellaneous sparsity information probes
  int getNumNonSparseStates();	
  void printHistSparsity (int numBins);

  // Variational free energy calculated using the sparse distribution
  FloatType sparseFreeEnergy();

  // Variational entropy calculated using the sparse distribution
  FloatType getSparseTotalNodeEntropy();


  void initializeAlg(InitBelType init);

  void initializeAlg();

 protected:




  void optimizeAlg(int nIterations);
  FloatType computeUpdateBeliefs();

  // Add the mean smoothness cost to a vector using the current SPARSE beliefs
  void addMeanSmoothnessCost(const int fixedPixel, const int varPixel, 
                             FloatType* dst );
  
  // Updates sparsities based on current approximate marginals
  void updateSparsity();
  void updateSparsity(int pixel);

  // Reset all sparsities to full
  void resetSparsity();

  // Array of (pointers to) sparsity indicators -- FALSE when the label is EXCLUDED
  bool** m_nodeSparsity;
  
  // Array of (pointers to) normalized, exponentiated sparse beliefs
  FloatType** m_spExpNodeBeliefs;
  
  // Tolerance for compressing beliefs based on KL Divergence
  FloatType m_divTol;

  // Minimum number of beliefs to keep when compressing
  int m_minBeams;

  //
  // Array arithmetic: shorthanded and inlined
  //

  inline void logNormalizeSparse (FloatType* array, bool* sp )
    { ArrayMath::logNormalizeSparse (array, sp, m_nLabels); }

 private:

  void initializeSparseMeanField();
  
  // Adds two numbers stored in logspace in a numerically stable fashion
  FloatType logSum(const FloatType a, const FloatType b);
  
};


#endif
