/**
 * Class implementing asynchronous mean field inference on a grid.
 * Note that all operations are conducted and values stored in log space. While we 
 * use terms like "message product," it's actually the sum of log messages.
 *
 * Jerod Weinman
 * jerod@acm.org
 */
#ifndef __ASYNCMEANFIELD_H__
#define __ASYNCMEANFIELD_H__

#include "SyncMeanField.h"
#include "ArrayMath.h"

class AsyncMeanField : public SyncMeanField {

 public:
  
    
  AsyncMeanField(int width, int height, int nLabels, EnergyFunction *eng);
  
  AsyncMeanField(int nPixels, int nLabels,EnergyFunction *eng);

  ~AsyncMeanField();
  
 protected:

  FloatType computeUpdateBeliefs(); 


  void optimizeAlg(int nIterations);


 private:
  
  
};

#endif
