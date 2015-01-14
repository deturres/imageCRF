/**
 * Class implementing sparse asynchronous mean-field on a grid.
 *
 * Jerod Weinman
 * jerod@acm.org
 */

#ifndef __SPARSEASYNCMEANFIELD_H__
#define __SPARSEASYNCMEANFIELD_H__

#include "SparseSyncMeanField.h"


class SparseAsyncMeanField : public SparseSyncMeanField {

 public: 

  SparseAsyncMeanField(int width, int height, int nLabels, EnergyFunction *eng);
  SparseAsyncMeanField(int nPixels, int nLabels, EnergyFunction *eng);
  ~SparseAsyncMeanField();

  FloatType sparseFreeEnergy();
  FloatType getSparseTotalNodeEntropy();

  void initializeAlg(InitBelType init);
 protected:

  void optimizeAlg(int nIterations);

  FloatType computeUpdateBeliefs();
};
  
#endif
