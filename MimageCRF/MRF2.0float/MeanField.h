/**
 * Class implementing mean field inference on a grid. Note that all
 * operations are conducted and values stored in log space. While we
 * use terms like "message product," it's actually the sum of log
 * messages.
 *
 * Jerod Weinman
 * jerod@acm.org
 */
#ifndef __MEANFIELD_H__
#define __MEANFIELD_H__

#include "MRFEnergy.h"
#include "ArrayMath.h"
#include <ostream>

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





class MeanField : public MRFEnergy {

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
      CURRENT,        // Uses current beliefs, whatever they may be
      RESET,      // Resets any derived quantities, but keeps current messages
      UNIFORM,    // Uniform beliefs
      LOCAL,      // Uses only local potentials
      LABELPOINT  // Point beliefs from labeling
    } InitBelType;

  MeanField(int width, int height, int nLabels, EnergyFunction *eng);

  MeanField(int nPixels, int nLabels,EnergyFunction *eng);

  ~MeanField();

  // Required by MRF -- NOT SUPPORTED - now supported
//chet  void setNeighbors(int pix1, int pix2, CostVal weight);

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

  // chet
  void setLog2(std::ostream *stream) { m_log = stream; }

protected:

  void computeBelief(int r, int c);
  void computeBelief(int pixel);
  void updateBelief(int pixel);

  // Calculate the difference between the current and updated beliefs
  FloatType beliefDifference(int pixel);

  FloatType computeUpdateBeliefs();

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

  // Sets answer using only data cost
  void initAnswer();

  //
  // Array arithmetic: inlined and shorthanded
  //

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


 private:

  void initializeMeanField();

};

#endif
