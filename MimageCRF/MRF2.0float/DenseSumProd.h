/**
 * Class outlining sum-product belief propagation on a grid.  Note
 * that all operations are conducted and values stored in log
 * space. While we use terms like "message product," it's actually the
 * sum of log messages.
 *
 * Jerod Weinman
 * jerod@acm.org
 *
 * $Id: DenseSumProd.h,v 1.1.1.1 2007/12/01 01:01:14 cpal Exp $
 */

#ifndef __DENSESUMPROD_H__
#define __DENSESUMPROD_H__

#include "SumProd.h"
#include "NodeMessages.h"
#include "ArrayMath.h"

#define FloatType float

class DenseSumProd : public virtual SumProd {

 public:
  
  DenseSumProd(int width, int height, int nLabels, EnergyFunction *eng);
  
  DenseSumProd(int nPixels, int nLabels,EnergyFunction *eng);

  ~DenseSumProd();
  
  // Required by MRF -- NOT SUPPORTED
  void setNeighbors(int pix1, int pix2, CostVal weight);

  void setParameters(int , void *);
 protected:
  
  // Compute the updated message between two nodes
  virtual void computeMessage(const int pixel1, const int pixel2,
                              NodeMessages* message1, NodeMessages* message2,
                              int direction);
    


private:


};

#endif
