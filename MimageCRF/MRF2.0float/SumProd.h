/**
 * Class outlining sum-product belief propagation on a grid.  Note
 * that all operations are conducted and values stored in log
 * space. While we use terms like "message product," it's actually the
 * sum of log messages.
 *
 * Derived classes will implement 
 *  - optimizeAlg (for synchronous/asynchronous scheduling)
 *  - computeMessage (for dense/sparse)
 *
 * Jerod Weinman
 * jerod@acm.org
 *
 * $Id: SumProd.h,v 1.1.1.1 2007/12/01 01:01:15 cpal Exp $
 */

#ifndef __SUMPROD_H__
#define __SUMPROD_H__

#include "MRFEnergy.h"
#include "NodeMessages.h"
#include "ArrayMath.h"

#define FloatType float
#define FLOATTYPE float

//////////////////////
// Optimizer Arguments
//////////////////////

// Message damping
#define PARAM_DAMPER 0 

// Message difference function (i.e., L1)
#define PARAM_MSGDFUN 1

// Message difference tolerance
#define PARAM_MSGDTOL 2

class SumProd : public virtual MRFEnergy {

 public:
  
  SumProd(int width, int height, int nLabels, EnergyFunction *eng);
  
  SumProd(int nPixels, int nLabels,EnergyFunction *eng);

  ~SumProd();
  
  // Required by MRF -- NOT SUPPORTED
  void setNeighbors(int pix1, int pix2, CostVal weight);

  void setParameters(int , void *);
	
  int getWidth();
  int getHeight();
  int getNumLabels();

  // Get the current approximate marginal for a pixel *in logspace*
  void getBelief(const int pixel, FloatType* dst);
  
  // Get the current approximate marginal for two neighboring pixels *in logspace*
  void getEdgeBelief(const int pixel1, const int pixel2, FloatType* dst);

  // Get the log probability for a configuration
  FloatType getLogProb(Label* labels);

  FloatType getTotalNodeEntropy();

  //#ifndef NDEBUG
  void printAll();
  void printArray(FloatType*);
  void printArray2(FloatType*);
  void printArray(int*);
  void printAllBeliefs();
  //#endif
  
 protected:
  
  SumProd();

  void initializeAlg();

  virtual void optimizeAlg(int nIterations)=0;
  
  // Compute the updated message between two nodes
  virtual void computeMessage(const int pixel1, const int pixel2,
                              NodeMessages* message1, NodeMessages* message2,
                              int direction)=0;
    

  // Array of (pointers to) messages into each node
  NodeMessages** m_nodeMessages;

  // Damping coefficient in a convex combination of new and old messages
  FloatType m_messageDamper;
  
  // Function to use for calculating message differences 
  NodeMessages::DiffType m_messageDiffFun;

  // Tolerance for convergence check on message differences
  FloatType m_messageDiffTol;

  // Sets answer to current Maximum Posterior Marginal (MPM) estimate
  void updateAnswer();

  // Array arithmetic
  inline int argMax(FloatType* array)
    { return ArrayMath::argMax(array, m_nLabels); };

  /** Copy vector src to dst */
  inline void copyTo(FloatType* dst, FloatType* src) 
    { ArrayMath::copyTo(dst, src, m_nLabels); };

  /** Adds vector src to accumulator dst */
  inline void addTo(FloatType* dst, FloatType* src)
    { ArrayMath::addTo(dst,src,m_nLabels); };

  /** Subtracts vector src from accumulator dst */
  inline void subtractFrom(FloatType* dst, FloatType* src) 
    { ArrayMath::subtractFrom(dst,src,m_nLabels); };

  /** Subtracts a constant from dst */
  inline void subtractFrom(FloatType* dst, FloatType sth)
    { ArrayMath::subtractFrom(dst,sth,m_nLabels);};

  /** Sums an array stored in log-space */
  inline FloatType logSum(FloatType* array)
    { return ArrayMath::logSum(array,m_nLabels); };
  
  /** Normalizes an array stored in log-space */
  inline void logNormalize(FloatType* array)
    { ArrayMath::logNormalize(array,m_nLabels); };

  /** Normalizes a matrix stored in log-space */
  inline void logNormalizeEdge(FloatType* array)
    { ArrayMath::logNormalize(array,m_nLabels*m_nLabels); };


private:

  void initializeSumProd();


};

#endif
