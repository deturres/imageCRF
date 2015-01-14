/**
 * Class implementing sparse asynchronous sum-product belief
 * propagation on a grid.
 *
 * Jerod Weinman
 * jerod@acm.org
 */

#ifndef __SPARSEASYNCSUMPROD_H__
#define __SPARSEASYNCSUMPROD_H__

#include "SparseSyncSumProd.h"

class SparseAsyncSumProd : public SparseSyncSumProd {

 public:
  
  SparseAsyncSumProd(int width, int height, int nLabels, EnergyFunction *eng);
  
  SparseAsyncSumProd(int nPixels, int nLabels, EnergyFunction *eng);

  ~SparseAsyncSumProd();


 protected:
  void optimizeAlg(int nIterations);
  



private:

  // Message update computation functions
  FloatType computeLeftMessages();
  FloatType computeRightMessages();
  FloatType computeUpMessages();
  FloatType computeDownMessages();


};

#endif
