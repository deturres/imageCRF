#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include <assert.h>
#include "ArrayMath.h"
#include "SyncMeanField.h"

/** 
 * You should be able to trace the general flow through:
 *  optimizeAlg() -> 
 *  computeMessages() -> 
 *  computeUpDownMessages(),computeLeftRightMessages() -> 
 *  computeMessage()
 */

#define PIXEL_IND(r,c) (r)*m_width+(c)

SyncMeanField::SyncMeanField(int width, int height, int nLabels, 
                             EnergyFunction *eng) :
  MRFEnergy(width,height,nLabels,eng)
{
  initializeMeanField();
}

SyncMeanField::SyncMeanField(int nPixels, int nLabels,EnergyFunction *eng) :
  MRFEnergy(nPixels,nLabels,eng)
{  
  std::cout << "Non grid graphs are not supported";
  assert(m_grid_graph);
}

SyncMeanField::~SyncMeanField()
{

  
  for (int i=0 ; i<m_nPixels ; i++)
    {
      delete[] m_nodeBeliefs[i];
      delete[] m_expNodeBeliefs[i];
      delete[] m_newNodeBeliefs[i];
    }
  
  delete[] m_nodeBeliefs;
  delete[] m_expNodeBeliefs;
  delete[] m_newNodeBeliefs;

}

/** Common initialization routine */
void SyncMeanField::initializeMeanField()
{

  m_beliefDamper = 1.0;
  m_beliefDiffFun = KLD;
  m_beliefDiffTol = .001;

  // Belief storage
  m_nodeBeliefs = new FloatType*[m_nPixels];
  m_newNodeBeliefs = new FloatType*[m_nPixels];
  m_expNodeBeliefs = new FloatType*[m_nPixels];

  //
  // Create beliefs
  //
  
  for (int i=0 ; i<m_nPixels ; i++)
    {
      m_nodeBeliefs[i] = new FloatType[m_nLabels];
      m_newNodeBeliefs[i] = new FloatType[m_nLabels];
      m_expNodeBeliefs[i] = new FloatType[m_nLabels];
    }

}


/** Initialization from parent classes.
 *
 * Start beliefs at uniform.
 */
void SyncMeanField::initializeAlg()
{
  initializeAlg(UNIFORM);
}

void SyncMeanField::initializeAlg(InitBelType init)
{

  switch (init)
    {
    case UNIFORM:
      {
        std::cout<<"SyncMeanField::initializeAlg(UNIFORM)\n";
        FloatType logZ = -log((FloatType)m_nLabels);
        FloatType Z = 1.0/((FloatType)m_nLabels);
        
        for (int i=0 ; i<m_nPixels ; i++)
          for (int b=0 ; b<m_nLabels ; b++)
            {
              m_nodeBeliefs[i][b] = logZ;
              m_expNodeBeliefs[i][b] = Z;
            }
        break;
      }
    case LOCAL:
      {
        std::cout<<"SyncMeanField::initializeAlg(LOCAL)\n";
        for (int i=0 ; i<m_nPixels ; i++)
          {
            copyDataCost (i,m_nodeBeliefs[i]);
            logNormalize(m_nodeBeliefs[i]);

            for (int b=0 ; b<m_nLabels ; b++)
              m_expNodeBeliefs[i][b] = exp(m_nodeBeliefs[i][b]);
      
          }
        break;
      }
    case LABELPOINT:
      {
        std::cout<<"SyncMeanField::initializeAlg(LABELPOINT)\n";
        //FloatType negInf = - std::numeric_limits<FloatType>::infinity();
        FloatType epsilon = 1e-8;
        FloatType logProb0 = log(epsilon) - log((FloatType)m_nLabels-1);
        FloatType logProb1 = log(1-epsilon);

        FloatType prob0 = epsilon/m_nLabels;
        FloatType prob1 = 1-epsilon;

        for (int i=0 ; i<m_nPixels ; i++)
          {
            for (int b=0 ; b<m_nLabels ; b++)
              {
                m_nodeBeliefs[i][b] = logProb0;
                m_expNodeBeliefs[i][b] = prob0;
              }

            m_nodeBeliefs[i][m_labels[i]] = logProb1;
            m_expNodeBeliefs[i][m_labels[i]] = prob1;
          }
        break;
      }
    }
  
  initAnswer();
}
/** Run the algorithm (synchronous plain-vanilla mean field)
 */     

void SyncMeanField::optimizeAlg (int nIterations)
{

  FloatType diff=0, energy;

#ifndef NDEBUG
  printAll();
#endif

  FloatType energy0 = freeEnergy();

  FloatType energyOld = energy0;

  for (int i=0 ; i<nIterations ; i++)
    {
      computeBeliefs();

#ifndef NDEBUG
      //      printAll();
#endif
      
      if (m_beliefDiffFun!=ENERGY)
        diff = beliefDifference();
      
      updateBeliefs();
      
      energy = freeEnergy();
      
      if (m_beliefDiffFun==ENERGY)
        {
          diff = 100*(energyOld-energy)/energyOld;
          energyOld = energy;
        }
      
      std::cout << "SyncMeanField: Iter = " << i << "\tDiff = " << diff << 
        "\tFree Energy = " << energy << "\n";

      if (diff < m_beliefDiffTol)
        {
          std::cout << "SyncMeanField converged in " << i << " iterations.\n";
          break;
        }
    }

  updateAnswer();

#ifndef NDEBUG
  printAll();
  printAllBeliefs();
#endif


}

/** Synchronously compute all new beliefs from the current beliefs.
 */
void SyncMeanField::computeBeliefs()
{
  for (int r=0 ; r<m_height ; r++)
    for (int c=0 ; c<m_width ; c++)
      computeBelief (r,c);
}


/** Compute a new belief at a pixel from the current beliefs.
 */
void SyncMeanField::computeBelief(int r, int c) 
{
  int pixC = PIXEL_IND(r,c);
  int pixL = PIXEL_IND(r,c-1);
  int pixR = PIXEL_IND(r,c+1);
  int pixU = PIXEL_IND(r-1,c);
  int pixD = PIXEL_IND(r+1,c);
        
  // Get the belief buffer
  FloatType* belief = m_newNodeBeliefs[pixC];
  
  // Start with -D(p)
  copyDataCost(pixC, belief);
  
  if (c>0)
    // Left neighbor
    addMeanSmoothnessCost(pixC, pixL, belief);
  
  if (c<m_width-1)
    // Right neighbor
    addMeanSmoothnessCost(pixC, pixR, belief);
  
  if (r>0)
    // Up neighbor
    addMeanSmoothnessCost(pixC, pixU, belief);
  
  if (r<m_height-1)
    // Down neighbor
    addMeanSmoothnessCost(pixC, pixD, belief);
  
  logNormalize(belief);
  
}

/** Add the mean smoothness cost to a vector using the current beliefs
 *
 * Parameters:
 *  - dst is a buffer over which fixedPixel indexes
 *
 * Output: M[i] = - sum_j b[j] * V(i,j)
 *
 */
void SyncMeanField::addMeanSmoothnessCost(const int fixedPixel, 
                                          const int varPixel, FloatType* dst)
{
  
  FloatType* belief0 = m_expNodeBeliefs[varPixel];
  FloatType* belief;

  register FloatType mean;

  for (int i=0 ; i<m_nLabels ; i++, dst++)
    {
      belief = belief0;
      
      mean = 0;

      for (int j=0 ; j<m_nLabels ; j++, belief++)
        {
          mean += (*belief) * getSmoothnessCost(fixedPixel, varPixel, i, j);

#ifndef NDEBUG
//           std::cout << "[" << fixedPixel << " : " << i << "]\t" <<
//             "[" << varPixel << " : " << j << "]\t" << 
//             "bel = " << *belief << "\t" <<
//             "V = " << getSmoothnessCost(fixedPixel, varPixel, i, j) << "\n";
#endif
          
        }

      *dst -= mean;
    }
}

/** Synchronous belief update once all have been computed */
void SyncMeanField::updateBeliefs()
{
  
  //
  // Update log-space beliefs
  //
  if (m_beliefDamper==1)
    {
      // Update by swapping pointers
      FloatType** tmp = m_nodeBeliefs;
      
      m_nodeBeliefs = m_newNodeBeliefs;
      m_newNodeBeliefs = tmp;
    }
  else
    {
      // Update by convex combination
      FloatType **belief = m_nodeBeliefs;
      FloatType **nbelief = m_newNodeBeliefs;

      for (int i=0 ; i<m_nPixels ; i++, belief++, nbelief++)
        convexLogAdd (m_beliefDamper, *belief, *nbelief);
    }
    
  //
  // Update/calculate exponentiated beliefs
  //
  FloatType *ebelief,*belief;

  for (int i=0 ; i<m_nPixels ; i++)
    {

      ebelief = m_expNodeBeliefs[i];
      belief = m_nodeBeliefs[i];

      for (int b=0 ; b<m_nLabels ; b++, ebelief++, belief++)
        *ebelief = exp(*belief);
    }
}


/** Asynchronous belief update */
void SyncMeanField::updateBelief(int pixel)
{
  
  //
  // Update log-space beliefs
  //
  if (m_beliefDamper==1)
    {
      // Update by swapping pointers
      FloatType* tmp = m_nodeBeliefs[pixel];
      
      m_nodeBeliefs[pixel] = m_newNodeBeliefs[pixel];
      m_newNodeBeliefs[pixel] = tmp;
    }
  else
    {
      // Update by convex combination
      FloatType *belief = m_nodeBeliefs[pixel];
      FloatType *nbelief = m_newNodeBeliefs[pixel];

      convexLogAdd (m_beliefDamper, belief, nbelief);
    }
    
  //
  // Update/calculate exponentiated beliefs
  //
  FloatType *ebelief = m_expNodeBeliefs[pixel];
  FloatType *belief = m_nodeBeliefs[pixel];
  
  for (int b=0 ; b<m_nLabels ; b++, ebelief++, belief++)
    *ebelief = exp(*belief);
  
}

/* Calculate the difference between the current and updated beliefs.
 *
 */

FloatType SyncMeanField::beliefDifference ()
{

  FloatType** ppb = m_nodeBeliefs;
  FloatType** ppe = m_expNodeBeliefs;
  FloatType** ppn = m_newNodeBeliefs;


  FloatType* pb;
  FloatType* pe;
  FloatType* pn;

  register FloatType diff0, diff = 0;

  switch (m_beliefDiffFun)
    {
    case NOP:
      break;
    case L1:

      for (int i=0 ; i<m_nPixels ; i++, ppb++, ppn++)
        {
          pb = *ppb;
          pn = *ppn;
          
          for (int b=0 ; b<m_nLabels ; b++, pb++, pn++)
            diff += fabs( convexLogDiff(*pb, *pn) );
        }
      
      
      break;
    case L2:

      for (int i=0 ; i<m_nPixels ; i++, ppb++, ppn++)
        {
            pb = *ppb;
            pn = *ppn;
            
            for (int b=0 ; b<m_nLabels ; b++, pb++, pn++)
              {
                diff0 = fabs( convexLogDiff(*pb, *pn) );
                diff += diff0*diff0;
              }
          }
        
        break;
        
    case MAX:

        diff = fabs ( **ppb - **ppn );

        for (int i=0 ; i<m_nPixels ; i++, ppb++, ppn++ )
          {
            pb = *ppb;
            pn = *ppn;
            
            for (int b=0 ; b<m_nLabels ; b++, pb++, pn++)
              { 
                diff0 = fabs( convexLogDiff(*pb, *pn) );
                if (diff0>diff)
                  diff = diff0;
              }
          }

        break;

    case KLD:

        for (int i=0 ; i<m_nPixels ; i++, ppb++, ppn++, ppe++)
          {
            pb = *ppb;
            pe = *ppe;
            pn = *ppn;
            
            for (int b=0 ; b<m_nLabels ; b++, pb++, pn++, pe++)
              diff += (*pe) * ( convexLogDiff(*pb, *pn) );
          }

        break;
    case ENERGY:

      diff = 0;
      break;

    default:
      std::cerr<<"SyncMeanField::beliefDifference " << m_beliefDiffFun << " not implemented yet. Ignoring.\n";
    }

  return diff;
}

/* Calculate the difference between the current and updated beliefs.
 *
 */

FloatType SyncMeanField::beliefDifference (int pixel)
{

  FloatType* pb = m_nodeBeliefs[pixel];
  FloatType* pe = m_expNodeBeliefs[pixel];
  FloatType* pn = m_newNodeBeliefs[pixel];


  register FloatType diff0, diff = 0;

  switch (m_beliefDiffFun)
    {
    case NOP:
      break;
    case L1:
      
      for (int b=0 ; b<m_nLabels ; b++, pb++, pn++)
        diff += fabs( *pb - *pn );
      
      break;
    case L2:
      
      for (int b=0 ; b<m_nLabels ; b++, pb++, pn++)
        {
          diff0 = fabs( *pb - *pn );
          diff += diff0*diff0;
        }
        
        break;
        
    case MAX:

      diff = fabs ( *pb - *pn );

      for (int b=0 ; b<m_nLabels ; b++, pb++, pn++)
        { 
          diff0 = fabs( *pb - *pn );
          if (diff0>diff)
            diff = diff0;
        }
      break;
      
    case KLD:
      
      for (int b=0 ; b<m_nLabels ; b++, pb++, pn++, pe++)
        diff += (*pe) * ( *pb - *pn );
      
      break;

    case ENERGY:

      diff = 0;
      break;
      
    default:
      std::cerr<<"SyncMeanField::beliefDifference " << m_beliefDiffFun << " not implemented yet. Ignoring.\n";
    }

  return diff;
}


/* Helper function for merging the message differences from two directions
 */
FloatType SyncMeanField::combineDifference(FloatType diff1, FloatType diff2)
{
  FloatType diff = 0;

switch (m_beliefDiffFun)
    {
    case NOP:
    case ENERGY:
      break;
    case L1:
    case L2:
    case KLD:
      diff = diff1+diff2;
      break;
    case MAX:
      diff = (diff2>diff1 ? diff2 : diff1);
      break;
      
    }

  return diff;
}


inline FloatType SyncMeanField::convexLogDiff(const FloatType a, 
                                              const FloatType b )

{
  if (m_beliefDamper==1)
    return a-b;
  else
    return a - log((1-m_beliefDamper)*exp(a)+m_beliefDamper*exp(b));
}
    

 /** Sets answer  to current Maximum Posterior Marginal (MPM) estimate
 */
void SyncMeanField::updateAnswer()
{
   Label* label = m_labels;

   FloatType** belief = m_nodeBeliefs;

   for (int i=0; i<m_nPixels ; i++, label++, belief++)
     *label = argMax (*belief);
            
}

void SyncMeanField::initAnswer()
{
  Label* label = m_labels;
  
  FloatType *cost = new FloatType[m_nLabels];

  for (int i=0 ; i<m_nPixels ; i++, label++)
    {
      copyDataCost(i, cost);
      *label = argMax (cost);
    }
  delete[] cost;
}

/** Calculate the current approximate marginal for a pair of neighboring pixels
 *  and write it to the destination array provided (must be of length m_nlabels^2.

 */
void SyncMeanField::getEdgeBelief(const int pixel1, const int pixel2, 
                                FloatType* dst)
{

  FloatType* bel1 = m_nodeBeliefs[pixel1];
  FloatType* bel2 = m_nodeBeliefs[pixel2];

  //
  // Pixel 2 
  //
  
  // Repeat the one dimensional belief over the ROWS of the destination

  FloatType *pd = dst;

  for (int b1=0 ; b1<m_nLabels ; b1++, pd+=m_nLabels)
    copyTo(pd, bel2);

  //
  // Pixel 1
  //
  
  // Repeat (add) the one dimensional belief over the COLS of the destination

  pd = dst;
  
  for (int b1=0 ; b1<m_nLabels ; b1++, bel1++)
    for (int b2=0 ; b2<m_nLabels ; b2++, pd++)
      *pd += *bel1;
  
}

/** Calculate the log probability of a labeling from the current beliefs
 */
FloatType SyncMeanField::getLogProb(Label* labels)
{
  
  register FloatType logprob = 0;

  FloatType** belief = m_nodeBeliefs;

  for (int i=0 ; i<m_nPixels ; i++, belief++, labels++)
    {
      logprob += *((*belief) + (*labels));
    }
  
  return logprob;
}

/** Parameter settings
 * 
 * See PARAM_* macros in SyncSumProd.h for options
 */
void SyncMeanField::setParameters (int param, void* value)
{
  
  switch (param)
    {
    case PARAM_DAMPER:
      {
        FloatType damp = *((FloatType*)value);
        
        if (damp>1 || damp<0)
          std::cerr << "SyncMeanField::setParameters: PARAM_DAMPER must be " <<
            "between zero and one. Ignoring.\n";
        else
          m_beliefDamper = *((FloatType*)value);
      
        break;
      }
    case PARAM_BELDFUN:
      m_beliefDiffFun = *((DiffType*)value);
      break;

    case PARAM_BELDTOL:
      m_beliefDiffTol = *((FloatType*)value);
      break;

    default:
      std::cerr << "SyncMeanField::setParameters: Unknown parameter index " <<
        param << '\n';
    }

}

 /** Add an edge to the graph.
 *
 * Note: Required by MRF, but NOT SUPPORTED
 */
void SyncMeanField::setNeighbors (int pix1, int pix2, CostVal weight)
{
  std::cerr << "SyncMeanField::setNeighbors not supported\n";
  assert(false);
}


FloatType SyncMeanField::freeEnergy()
{

  register FloatType energy = 0;

  int pixel,pixel1;
  
  FloatType *belief, *belief1;

  for (int r=0 ; r<m_height ; r++)
    for (int c=0 ; c<m_width ; c++)
      {
        pixel = PIXEL_IND(r,c);

        // Data terms
        belief = m_expNodeBeliefs[pixel];

        for (int b=0 ; b<m_nLabels ; b++ )
          energy += belief[b] * getDataCost(pixel,b);

        // Smoothness terms

        if (r<m_height-1)
          { // Down
            pixel1 = PIXEL_IND(r+1,c);

            belief = m_nodeBeliefs[pixel];
            belief1 = m_nodeBeliefs[pixel1];
            

            for (int b=0 ; b<m_nLabels ; b++ )
              for (int b1=0 ; b1<m_nLabels ; b1++)
                energy += exp(belief[b]+belief1[b1]) * 
                  getSmoothnessCost(pixel,pixel1,b,b1);
            
          }

        if (c<m_width-1)
          { // Right
            pixel1 = PIXEL_IND(r,c+1);

            belief = m_nodeBeliefs[pixel];
            belief1 = m_nodeBeliefs[pixel1];
            

            for (int b=0 ; b<m_nLabels ; b++ )
              for (int b1=0 ; b1<m_nLabels ; b1++)
                energy += exp(belief[b]+belief1[b1]) * 
                  getSmoothnessCost(pixel,pixel1,b,b1);
            
          }
            

      }

  return energy - getTotalNodeEntropy();
    

  
  
}


FloatType SyncMeanField::getTotalNodeEntropy()
{

  register FloatType entropy = 0;

  FloatType** beliefs = m_nodeBeliefs;
  FloatType** ebeliefs = m_expNodeBeliefs;

  FloatType *bel,*ebel;

  for (int i=0 ; i<m_nPixels ; i++, beliefs++, ebeliefs++)
    {
      bel = *beliefs;
      ebel = *ebeliefs;

      for (int b=0 ; b<m_nLabels ; b++, bel++, ebel++)
        entropy += (*ebel) * (*bel);
    }

  return -entropy;
}

FloatType SyncMeanField::getNodeEntropy(const int pixel)
{

  register FloatType entropy = 0;

  FloatType* belief = m_nodeBeliefs[pixel];
  FloatType* ebelief = m_expNodeBeliefs[pixel];

  for (int b=0 ; b<m_nLabels ; b++, belief++, ebelief++)
        entropy += (*ebelief) * (*belief);


  return -entropy;
}

#ifndef NDEBUG

void SyncMeanField::printAll()
{

//  MRFEnergy::printAll();

  std::cout << "SyncSumProd(" << m_beliefDiffFun << 
    ", " << m_beliefDiffTol << ")\n";

  for (int i=0 ; i<m_nPixels ; i++)
    {
      std::cout<<"Node " << i << "\n";
      printArray(m_nodeBeliefs[i]);
      printArray(m_expNodeBeliefs[i]);
      printArray(m_newNodeBeliefs[i]);
    }

  
  std::cout << "\n";
}

void SyncMeanField::printArray(FloatType* array)
{
  for (int i=0 ; i<m_nLabels ; i++, array++)
    std::cout << *array <<"\t";
  std::cout << "\n";

}

void SyncMeanField::printArray2(FloatType* array)
{
  for (int i=0 ; i<m_nLabels*m_nLabels ; i++, array++)
    std::cout << *array <<"\t";
  std::cout << "\n";

}

void SyncMeanField::printAllBeliefs()
{

  for (int i=0 ; i<m_nPixels ; i++)
    {
      std::cout << "Node " << i << ":\t";
      printArray(m_nodeBeliefs[i]);
      std::cout << "\n";
    }
  std::cout << "\n";

  FloatType *belief2 = new FloatType[m_nLabels*m_nLabels];

  // Horizontal Edge Beliefs
  for (int r=0 ; r<m_height ; r++)
    for (int c=1 ; c<m_width ; c++)
      {
        std::cout << "Edge ("<<r<<","<<c-1<<")->("<<r<<","<<c<<"):\t";
        getEdgeBelief(PIXEL_IND(r,c-1),PIXEL_IND(r,c),belief2);
        printArray2(belief2);
      }

  // Vertical Edge Beliefs
  for (int r=1 ; r<m_height ; r++)
    for (int c=0 ; c<m_width ; c++)
      {
        std::cout << "Edge ("<<r-1<<","<<c<<")->("<<r<<","<<c<<"):\t";
        getEdgeBelief(PIXEL_IND(r-1,c),PIXEL_IND(r,c),belief2);
        printArray2(belief2);
      }
	delete[] belief2;
}

#endif

