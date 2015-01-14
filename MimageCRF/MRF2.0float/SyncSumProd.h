/**
 * Class implementing synchronous sum-product belief propagation on a grid.
 * Note that all operations are conducted and values stored in log space. While we 
 * use terms like "message product," it's actually the sum of log messages.
 *
 * Jerod Weinman
 * jerod@acm.org
 */

#ifndef __SYNCSUMPROD_H__
#define __SYNCSUMPROD_H__

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

class SyncSumProd : public MRFEnergy {

 public:
  
  SyncSumProd(int width, int height, int nLabels, EnergyFunction *eng);
  
  SyncSumProd(int nPixels, int nLabels,EnergyFunction *eng);

  ~SyncSumProd();
  
  // Required by MRF -- NOT SUPPORTED
  void setNeighbors(int pix1, int pix2, CostVal weight);

  void setParameters(int , void *);
	
  int getWidth();
  int getHeight();
  int getNumLabels();

  // Get the current approximate marginal for a pixel *in logspace*
  void getBelief(const int pixel, FloatType* dst);
  
  FloatType getBelief(const int pixel, const Label label);

  // Get the current approximate marginal for two neighboring pixels *in logspace*
  void getEdgeBelief(const int pixel1, const int pixel2, FloatType* dst);

  // Get the log probability for a configuration
  FloatType getLogProb(Label* labels);

  FloatType getTotalNodeEntropy();

  
 protected:
  void initializeAlg();

  void optimizeAlg(int nIterations);
  
  // Compute the updated message between two nodes
  virtual void computeMessage(const int pixel1, const int pixel2,
                              NodeMessages* message1, NodeMessages* message2,
                              int direction);
    

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

  // Array arithmetic: shorthanded and inlined 

  /** Index of array max */
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

  // Message update computation functions
  void computeMessages();
  void computeLeftRightMessages();
  void computeUpDownMessages();


};

#endif
