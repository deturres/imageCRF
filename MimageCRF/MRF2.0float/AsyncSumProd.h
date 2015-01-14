/**
 * Class implementing asynchronous sum-product belief propagation on a grid.
 * Note that all operations are conducted and values stored in log space. While we 
 * use terms like "message product," it's actually the sum of log messages.
 *
 * Jerod Weinman
 * jerod@acm.org
 */

#ifndef __ASYNCSUMPROD_H__
#define __ASYNCSUMPROD_H__

#include "SyncSumProd.h"

class AsyncSumProd : public SyncSumProd {

 public:
  
  AsyncSumProd(int width, int height, int nLabels, EnergyFunction *eng);
  
  AsyncSumProd(int nPixels, int nLabels,EnergyFunction *eng);

  ~AsyncSumProd();
  
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
