/**
 * Class implementing synchronous mean field inference on a grid.
 * Note that all operations are conducted and values stored in log space. While we 
 * use terms like "message product," it's actually the sum of log messages.
 *
 * Jerod Weinman
 * jerod@acm.org
 */
#ifndef __SYNCMEANFIELD_H__
#define __SYNCMEANFIELD_H__

#include "MRFEnergy.h"
#include "ArrayMath.h"

#define FloatType float
#define FLOATTYPE float

//////////////////////
// Optimizer Arguments
//////////////////////

// Belief damping
#define PARAM_DAMPER 0 

// Belief difference function (i.e., L1)
#define PARAM_BELDFUN 1

// Belief difference tolerance
#define PARAM_BELDTOL 2

class SyncMeanField : public MRFEnergy {

 public:

  // Type for measuring the difference between current and updated beliefs
  typedef enum
    {
      NOP,
      L1, 
      L2,
      MAX,
      KLD,
      ENERGY // Variational free energy E_Q[log P] - H(Q) + C
    } DiffType;

  // Type for initializing beliefs
  typedef enum
    {
      UNIFORM,    // Uniform beliefs
      LOCAL,      // Uses only local potentials
      LABELPOINT  // Point beliefs from labeling
    } InitBelType;
  
    
  SyncMeanField(int width, int height, int nLabels, EnergyFunction *eng);
  
  SyncMeanField(int nPixels, int nLabels,EnergyFunction *eng);

  ~SyncMeanField();
  
  // Required by MRF -- NOT SUPPORTED
  void setNeighbors(int pix1, int pix2, CostVal weight);

  void setParameters(int , void *);

  int getWidth();
  int getHeight();
  int getNumLabels();

  // Entropy of node marginals (factorized distribution)
  FloatType getTotalNodeEntropy();
  FloatType getNodeEntropy(const int pixel);

  // Variational free energy E_Q[U+V] - H(Q)
  FloatType freeEnergy();

  // Get the current approximate marginal for a pixel *in logspace*
  inline void getBelief(const int pixel, FloatType* dst) 
    { ArrayMath::copyTo(dst,m_nodeBeliefs[pixel],m_nLabels); }
  
  inline FloatType getBelief(const int pixel, const Label label) 
    { return m_nodeBeliefs[pixel][label]; }
  
  // Get the current approximate marginal for two neighboring pixels *in logspace*
  void getEdgeBelief(const int pixel1, const int pixel2, FloatType* dst);

  // Get the log probability for a configuration
  FloatType getLogProb(Label* labels);

  void initializeAlg(InitBelType init);
 protected:

  void computeBeliefs(); 
  void computeBelief(int r, int c);

  void updateBeliefs();
  void updateBelief(int pixel);
  
  // Calculate the difference between the current and updated beliefs
  FloatType beliefDifference();
  FloatType beliefDifference(int pixel);

  FloatType combineDifference(FloatType diff1, FloatType diff2);

  void initializeAlg(); // specified by MRF



  void optimizeAlg(int nIterations);

  // Add the mean smoothness cost to a vector using the current beliefs
  virtual void addMeanSmoothnessCost(const int fixedPixel, const int varPixel, 
                                     FloatType* dst);

  
  // Array of (pointers to) beliefs at each node -- stored in log-space
  FloatType** m_nodeBeliefs;
  FloatType** m_newNodeBeliefs;
  
  // Array of (pointers to) beliefs at each node stored in regular space
  FloatType** m_expNodeBeliefs;

  FloatType m_beliefDamper;
  DiffType m_beliefDiffFun;
  FloatType m_beliefDiffTol;

  // Sets answer to current Maximum Posterior Marginal (MPM) estimate
  void updateAnswer();

#ifndef NDEBUG
  void printAll();
  void printArray(FloatType*);
  void printArray2(FloatType*);
  void printAllBeliefs();
#endif

  
  // Array arithmetic

  /** Copy vector src to dst */
  inline void copyTo(FloatType* dst, FloatType* src) 
    { ArrayMath::copyTo(dst, src, m_nLabels); };

  inline int argMax(FloatType* array)
    { return ArrayMath::argMax(array, m_nLabels); };

  /** Normalizes an array stored in log-space */
  inline void logNormalize(FloatType* array)
    { ArrayMath::logNormalize(array,m_nLabels); };

  /** Normalizes a matrix stored in log-space */
  inline void logNormalizeEdge(FloatType* array)
    { ArrayMath::logNormalize(array,m_nLabels*m_nLabels); };

  /** Convex update of an array. */
  inline void convexLogAdd(const FloatType alpha, FloatType* dst, FloatType* src)
    { ArrayMath::convexLogAdd(alpha,dst,src,m_nLabels); }

  // Convex log update, then difference  
  // a - log((1-alpha)*exp(a)+alpha*exp(b))
  FloatType convexLogDiff(const FloatType a, const FloatType b );


  void initAnswer();

 private:
  
  void initializeMeanField();

  
};

#endif
