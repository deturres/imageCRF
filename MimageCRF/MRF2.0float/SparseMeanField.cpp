#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "ArrayMath.h"
#include "SparseMeanField.h"
#include <time.h>


using std::endl;

/** 
 * You should be able to trace the general flow through:
 *  optimizeAlg() -> 


 */

SparseMeanField::SparseMeanField(int width, int height, int nLabels, 
                             EnergyFunction *eng) :
  MeanField(width,height,nLabels,eng)
{
  initializeSparseMeanField();
}

SparseMeanField::SparseMeanField(int nPixels, int nLabels,
                                         EnergyFunction *eng) :
  MeanField(nPixels,nLabels,eng)
{  
  std::cerr << "Non grid graphs are not supported";
  assert(m_grid_graph);
}

SparseMeanField::~SparseMeanField()
{

  
  for (int i=0 ; i<m_nPixels ; i++)
    {
      delete[] m_nodeSparsity[i];
      delete[] m_spExpNodeBeliefs[i];
    }
  
  delete[] m_nodeSparsity;
  delete[] m_spExpNodeBeliefs;


}

/** Common initialization routine */
void SparseMeanField::initializeSparseMeanField()
{

  m_divTol = 0.001;
  m_minBeams = 1;


  // Sparsity storage
  m_nodeSparsity = new bool*[m_nPixels];
  m_spExpNodeBeliefs = new FloatType*[m_nPixels];

  // Create sparsity terms
  for (int i=0 ; i<m_nPixels ; i++)
    {
      m_nodeSparsity[i] = new bool[m_nLabels];
      m_spExpNodeBeliefs[i] = new FloatType[m_nLabels];
    }


}


/** Initialization from parent classes */
void SparseMeanField::initializeAlg()
{
  MeanField::initializeAlg();

  // Initialize all sparsities to true (since all are initially included)
  // and sparse beliefs to whatever the parent initializes them to
  SparseMeanField::resetSparsity();


}


void SparseMeanField::initializeAlg(InitBelType init)
{

  switch(init)
    {
    case CURRENT:
      break;
    case RESET:
      resetSparsity();
    case UNIFORM:
      SparseMeanField::initializeAlg ();
      break;
    case LOCAL:
      MeanField::initializeAlg(init);
      //updateSparsity();
      break;
    case LABELPOINT:
      {
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
                m_spExpNodeBeliefs[i][b] = 0;
                m_nodeSparsity[i][b] = false;
              }

            m_nodeBeliefs[i][m_labels[i]] = logProb1;
            m_expNodeBeliefs[i][m_labels[i]] = prob1;

            m_spExpNodeBeliefs[i][m_labels[i]] = 1.0;
            m_nodeSparsity[i][m_labels[i]] = true;
          }
        break;
      }
    }

  initAnswer();  

}

/** Run the algorithm (asynchronous sparse mean field)
 */     

void SparseMeanField::optimizeAlg (int nIterations)
{

  FloatType diff, energy = 0;

  FloatType energy0;

  if (m_beliefDiffFun==ENERGY)
    energy0 = sparseFreeEnergy();
  else
    energy0 = 0;

  FloatType energyOld = energy0;
  FloatType time = 0;
  
  // Log format is (tab separated): iteration energy time

  //// Print initial free energy and info ////
  *m_log << -1 << "\tNaN\t";

  m_log->precision(1);
  m_log->setf(std::ios::fixed,std::ios::floatfield);
  *m_log << energy0 <<"\t";

  m_log->precision(2);
  m_log->setf(std::ios::fixed,std::ios::floatfield);
  *m_log << time << endl;

  // Write intermediate image
  if (m_imWriter!=0)
    {
      updateAnswer();
      m_imWriter->write(m_labels);
    }


  clock_t start, end;

  for (int i=0 ; i<nIterations ; i++)
    {

      start = clock();

      diff = computeUpdateBeliefs();
  
      if (m_beliefDiffFun==ENERGY)
        {
          energy = sparseFreeEnergy();
          diff = 100*(energyOld-energy)/energyOld;
          energyOld = energy;
        }

      end = clock();

      time += (float) (((double)(end-start)) / CLOCKS_PER_SEC);
      
      //// Print iteration free energy and info ////
      *m_log << i << "\t";

      m_log->precision(5);
      m_log->setf(std::ios::fixed,std::ios::floatfield);
      *m_log << diff << "\t";

      m_log->precision(1);
      m_log->setf(std::ios::fixed,std::ios::floatfield);
      *m_log << energy << "\t";

      m_log->precision(2);
      m_log->setf(std::ios::fixed,std::ios::floatfield);
      *m_log << time << endl;
      ////

      if (m_imWriter!=0)
        {
          updateAnswer();
          m_imWriter->write(m_labels);
        }

      if (diff < m_beliefDiffTol)
        {
          break;
        }
    }

  updateAnswer();


}


/** Asynchronously compute all new beliefs from the current beliefs.
 */
FloatType SparseMeanField::computeUpdateBeliefs()
{
  
  FloatType diff,diff0;

  diff = 0;
  int pixel;

  for (int r=0 ; r<m_height ; r++)
    for (int c=0 ; c<m_width ; c++)
      {
        computeBelief (r,c);

        pixel = pixelIndex (r,c);

        diff0 = beliefDifference (pixel);

        diff = combineDifference (diff,diff0);

        updateBelief (pixel);
        updateSparsity (pixel);

      }

  return diff;
}

/** Variational free energy calculated using the sparse distribution
 */
FloatType SparseMeanField::sparseFreeEnergy()
{

  register FloatType energy = 0;

  int pixel,pixel1;
  
  FloatType *belief, *belief1;
  bool *sp, *sp1;

  for (int r=0 ; r<m_height ; r++)
    for (int c=0 ; c<m_width ; c++)
      {
        pixel = pixelIndex(r,c);

        // Data terms
        belief = m_spExpNodeBeliefs[pixel];
        sp = m_nodeSparsity[pixel];

        for (int b=0 ; b<m_nLabels ; b++ )
          if (sp[b])
            energy += belief[b] * getDataCost(pixel,b);

        // Smoothness terms

        if (r<m_height-1)
          { // Down
            pixel1 = pixelIndex(r+1,c);

            belief = m_spExpNodeBeliefs[pixel];
            belief1 = m_spExpNodeBeliefs[pixel1];
            
            sp = m_nodeSparsity[pixel];
            sp1 = m_nodeSparsity[pixel1];

            for (int b=0 ; b<m_nLabels ; b++ )
              for (int b1=0 ; b1<m_nLabels ; b1++)
                if (sp[b] && sp1[b1])
                  energy += belief[b] * belief1[b1] * 
                    getSmoothnessCost(pixel,pixel1,b,b1);
            
          }

        if (c<m_width-1)
          { // Right
            pixel1 = pixelIndex(r,c+1);

            belief = m_spExpNodeBeliefs[pixel];
            belief1 = m_spExpNodeBeliefs[pixel1];
            
            sp = m_nodeSparsity[pixel];
            sp1 = m_nodeSparsity[pixel1];
            
            for (int b=0 ; b<m_nLabels ; b++ )
              for (int b1=0 ; b1<m_nLabels ; b1++)
                if (sp[b] && sp1[b1])
                  energy += belief[b] * belief1[b1]  * 
                    getSmoothnessCost(pixel,pixel1,b,b1);
            
          }
            

      }

  return energy - getSparseTotalNodeEntropy();
    

  
  
}

/** Variational entropy calculated using the sparse distribution 
 */
FloatType SparseMeanField::getSparseTotalNodeEntropy()
{

  register FloatType entropy = 0;

  FloatType** beliefs = m_spExpNodeBeliefs;
  FloatType* bel;

  bool** sparsity = m_nodeSparsity;
  bool* sp;
  for (int i=0 ; i<m_nPixels ; i++, beliefs++, sparsity++)
    {
      bel = *beliefs;
      sp = *sparsity;

      for (int b=0 ; b<m_nLabels ; b++, bel++, sp++)
        if (*sp)
          entropy +=  (*bel) * log(*bel);
    }
  
  return -entropy;
}


/** Add the mean smoothness cost to a vector using the current SPARSE beliefs
 *
 * Parameters:
 *  - dst is a buffer over which fixedPixel indexes
 *
 * Output: M[i] = - sum_j s[j] * b[j] * V(i,j)
 *
 */
void SparseMeanField::addMeanSmoothnessCost(const int fixedPixel, 
                                               const int varPixel, FloatType* dst)
{

  bool* sparsity0 = m_nodeSparsity[varPixel];
  bool* sparsity;

  FloatType* belief0 = m_spExpNodeBeliefs[varPixel];
  FloatType* belief;

  register FloatType mean;

  for (int i=0 ; i<m_nLabels ; i++, dst++)
    {
      belief = belief0;
      sparsity = sparsity0;
      mean = 0;

      for (int j=0 ; j<m_nLabels ; j++, belief++, sparsity++)
        {
          if (*sparsity)
            {
              mean += (*belief) * getSmoothnessCost(fixedPixel, varPixel, i, j);
            }          
        }

      *dst -= mean;
    }
}

/** Resets the sparsity at all pixels
 */
void SparseMeanField::resetSparsity()
{

  bool *pSparse;

  for (int i=0 ; i<m_nPixels ; i++)
  {
    // Copy non-sparse exponentiated representation over to sparse storage
    copyTo (m_spExpNodeBeliefs[i], m_expNodeBeliefs[i]);
    
    pSparse = m_nodeSparsity[i];
    
    // Include all labels
    for (int j=0 ; j<m_nLabels ; j++, pSparse++)
        *pSparse = true;
  }
}

/** Updates the sparsity at a pixel by sorting the probabilities and 
 *  calculating the divergence from the original distribution as a result of 
 *  eliminating the states (zeroing their probabilities).
 */
void SparseMeanField::updateSparsity(int pixel)
{

  FloatType *belSort = new FloatType[m_nLabels];
  assert(belSort!=0);

  int kBeam = 1;   // Index of cut-off

  FloatType kThresh; // Belief threshold

  bool* pSparse;
  FloatType* pBelief;

  copyTo (belSort, m_nodeBeliefs[pixel]);
  
  // Sorts in ascending order (doh!)
  std::sort(belSort,belSort+m_nLabels);
      
  // Find index of cut-off: largest index (from end)
  // where divergence is less than the threshold
  FloatType *bel = belSort+m_nLabels-1;
  
  FloatType logcpf = *bel;

  bel--;
  
  for ( kBeam = 0 ; 
        kBeam < m_nLabels && kBeam <= m_nLabels-m_minBeams && -logcpf>=m_divTol; 
        kBeam++, bel--)
    {
      logcpf = logSum(logcpf,*bel);
    }


  // Update sparsity 
  pSparse = m_nodeSparsity[pixel];

  // We may have advanced and then failed the test
  // in this case, there can be no pruning
  if (kBeam==m_nLabels)
    {
      for (int b=0 ; b<m_nLabels ; b++, pSparse++ )
        *pSparse = true;
      
      // Update sparse exponentiated beliefs
      copyTo(m_spExpNodeBeliefs[pixel],m_expNodeBeliefs[pixel]);
      
    } 
  else
    {
      kThresh = belSort[m_nLabels-1-kBeam];
      
      pBelief = m_nodeBeliefs[pixel];
      
      for (int b=0 ; b<m_nLabels ; b++, pSparse++, pBelief++ )
        *pSparse = (*pBelief >= kThresh);
      
      // Update (normalize) sparse exponentiated beliefs
      pSparse = m_nodeSparsity[pixel];
      pBelief = m_spExpNodeBeliefs[pixel];
      
      copyTo (pBelief, m_nodeBeliefs[pixel]); 
      
      logNormalizeSparse (pBelief, pSparse); 
      
      // exponentiate 
      for (int b=0 ; b<m_nLabels ; b++, pBelief++, pSparse++)
        if (*pSparse)
          *pBelief = exp(*pBelief);

    }
  delete[] belSort;
}

/** Print a histogram of the number of states at the pixels */
void SparseMeanField::printHistSparsity(int numBins)
{

  register int sum=0;
  bool** ps = m_nodeSparsity;
  bool* psn;

  int *hist = new int[numBins];
  assert(hist!=0);

  for (int c=0 ; c<numBins; c++)
    hist[c] = 0;

  for (int i=0 ; i<m_nPixels ; i++, ps++)
    {
      psn = *ps;
      sum = 0;
      
      for (int b=0 ; b<m_nLabels ; b++, psn++)
        sum+= (*psn ? 1 : 0);

      hist[(int)floor( ((float)sum)/((float)(m_nLabels)) * numBins)]++;

      
    }

  *m_log<<"Sparsity Histogram\n";

  for (int c=0 ; c<numBins ; c++)
    *m_log<<hist[c]<<"\t";
  *m_log<<"\n";
  delete[] hist;
  
}

/** Parameter settings
 * 
 * See PARAM_* macros in SparseMeanField.h and MeanField.h for options
 */
void SparseMeanField::setParameters (int param, void* value)
{
  
  switch (param)
    {
    case PARAM_DIVTOL:
      {
        FloatType tol = *((FloatType*)value);
        
        if (tol<0)
          std::cerr<<"SparseMeanField::setParameters: PARAM_DIVTOL must " <<
            "be positive. Ignoring value " << tol << ".\n";
        else
          m_divTol = tol;
        
        break;
      }
    case PARAM_MINBEAMS:
      {
        int beams = *((int*)value);
        
        if (beams<1 || beams>m_nLabels)
          std::cerr<<"SparseMeanField::setParameters: PARAM_MINBEAMS must " <<
            "be between 1 and " << m_nLabels << ". Ignoring value " << beams << 
            ".\n";
        else
          m_minBeams = beams;
        
      break;
      }
    default:
      MeanField::setParameters(param,value);
    }

}


/** Total number of non-sparse states */
int SparseMeanField::getNumNonSparseStates()
{
  bool** ps = m_nodeSparsity;
  bool* psn;

  register int sum = 0;
  for (int i=0 ; i<m_nPixels ; i++, ps++)
    {
      psn = *ps;

      for (int b=0 ; b<m_nLabels ; b++,psn++)
        sum += (*psn ? 1 : 0);
    }
  
  return sum;
}

/** Adds two numbers stored in logspace in a numerically stable fashion */
inline FloatType SparseMeanField::logSum(const FloatType a, const FloatType b)
{
  if (a>b)
    return a+log(exp(b-a)+1);
  else
    return b+log(exp(a-b)+1);
}
