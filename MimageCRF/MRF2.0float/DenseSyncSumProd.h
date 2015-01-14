/**
 * Class implementing dense synchronous sum-product belief
 * propagation on a grid.
 *
 * Jerod Weinman
 * jerod@acm.org
 *
 * $Id: DenseSyncSumProd.h,v 1.1.1.1 2007/12/01 01:01:14 cpal Exp $
 *
 */

#ifndef __DENSESYNCSUMPROD_H__
#define __DENSESYNCSUMPROD_H__

#include "SyncSumProd.h"
#include "DenseSumProd.h"

class DenseSyncSumProd : public DenseSumProd, public SyncSumProd {

 public:
  
  DenseSyncSumProd(int width, int height, int nLabels, EnergyFunction *eng);
  
  DenseSyncSumProd(int nPixels, int nLabels,EnergyFunction *eng);

  ~DenseSyncSumProd();

  
 protected:

  



private:
  
};

#endif
