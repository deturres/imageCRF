#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "ArrayMath.h"
#include "SparseSyncMeanField.h"

#ifndef NDEBUG
#include <time.h>
#endif

using std::endl;

/** 
 * You should be able to trace the general flow through:
 *  optimizeAlg() -> 
 *  computeMessages() -> 
 *  computeUpDownMessages(),computeLeftRightMessages() -> 
 *  computeMessage()
 */


// $Id: SparseSyncMeanField.cpp,v 1.1.1.1 2007/12/01 01:01:15 cpal Exp $

SparseSyncMeanField::SparseSyncMeanField(int width, int height, int nLabels, 
                             EnergyFunction *eng) :
  SyncMeanField(width,height,nLabels,eng)
{
  initializeSparseMeanField();
}

SparseSyncMeanField::SparseSyncMeanField(int nPixels, int nLabels,
                                         EnergyFunction *eng) :
  SyncMeanField(nPixels,nLabels,eng)
{  
  *m_log << "Non grid graphs are not supported";
  assert(m_grid_graph);
}

SparseSyncMeanField::~SparseSyncMeanField()
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
void SparseSyncMeanField::initializeSparseMeanField()
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
void SparseSyncMeanField::initializeAlg()
{
  SyncMeanField::initializeAlg();

  // Initialize all sparsities to true (since all are initially included)
  // and sparse beliefs to whatever the parent initializes them to
  for (int i=0 ; i<m_nPixels ; i++)
    {
      copyTo (m_spExpNodeBeliefs[i], m_expNodeBeliefs[i]);

      for (int b=0 ; b<m_nLabels ; b++)
        m_nodeSparsity[i][b] = true;
    }


}

/** Run the algorithm (synchronous sparse mean field)
 */     

void SparseSyncMeanField::optimizeAlg (int nIterations)
{

  FloatType diff;

#ifndef NDEBUG
  printAll();
#endif

  for (int i=0 ; i<nIterations ; i++)
    {
      computeBeliefs();

#ifndef NDEBUG
      printAll();
#endif
      
      diff = beliefDifference();
      
      *m_log << "SyncMeanField: Iter = " << i << "\tDiff = " << diff <<"\n";
      
      updateBeliefs();

      updateSparsity();
      

      *m_log<< "Average sparsity = " << 
        ((float)getNumNonSparseStates())/((float)m_nPixels*m_nLabels)<<'\n';


      printHistSparsity(10);
      
      if (diff < m_beliefDiffTol)
        {
          *m_log << "SparseSyncMeanField converged in " << i << " iterations.\n";
          break;
        }
    }

  updateAnswer();

#ifndef NDEBUG
  printAll();
  printAllBeliefs();
#endif


}

/** Add the mean smoothness cost to a vector using the current SPARSE beliefs
 *
 * Parameters:
 *  - dst is a buffer over which fixedPixel indexes
 *
 * Output: M[i] = - sum_j s[j] * b[j] * V(i,j)
 *
 */
void SparseSyncMeanField::addMeanSmoothnessCost(const int fixedPixel, 
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

#ifndef NDEBUG
//           *m_log << "[" << fixedPixel << " : " << i << "]\t" <<
//             "[" << varPixel << " : " << j << "]\t" << 
//             "bel = " << *belief << "\t" <<
//             "V = " << getSmoothnessCost(fixedPixel, varPixel, i, j) << "\n";
#endif
            }          
        }

      *dst -= mean;
    }
}

/** Updates sparsities based on current approximate marginals */
void SparseSyncMeanField::updateSparsity()
{

  for (int i=0 ; i<m_nPixels ; i++)
    updateSparsity (i);

}

void SparseSyncMeanField::updateSparsity(int pixel)
{

  FloatType *belSort = new FloatType[m_nLabels];
  FloatType *divSort = new FloatType[m_nLabels];

  assert(belSort!=0);
  assert(divSort!=0);

  int kBeam = 1;   // Index of cut-off

  FloatType kThresh; // Belief threshold

  bool* pSparse;
  FloatType* pBelief;

 #ifndef NDEBUG
//   *m_log <<"Node "<<pixel<<":\t";
//   SyncMeanField::printArray(m_nodeBeliefs[pixel]);
#endif
  
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



#ifndef JDEBUG
  //*m_log<< "Node " << pixel << "\n";
  //*m_log << "Beliefs:\t";
  //printArray(belief);
  //*m_log << "Sorted:\t";
  //printArray(belSort);
  //*m_log << "Divergences\t";
  //printArray(divSort);
  //*m_log << "Beam=" << kBeam << "\tThresh=" << kThresh <<'\n';
  //*m_log << "Sparsity\t";
  //printArray(m_nodeSparsity[pixel]);
#endif
  
#ifndef NDEBUG
//   pSparse = m_nodeSparsity[pixel];
//   int ss = 0;
//   for (int b=0 ; b<m_nLabels ; b++,pSparse++)
//     if (*pSparse)
//       ss++;
  
//   *m_log<<"Total States:\t"<<ss<<"\n";
#endif
	delete[] belSort;
	delete[] divSort;
}

void SparseSyncMeanField::printHistSparsity(int numBins)
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
 * See PARAM_* macros in SparseSyncMeanField.h for options
 */
void SparseSyncMeanField::setParameters (int param, void* value)
{
  
  switch (param)
    {
    case PARAM_DIVTOL:
      {
        FloatType tol = *((FloatType*)value);
        
        if (tol<0)
          std::cerr<<"SparseSyncMeanField::setParameters: PARAM_DIVTOL must " <<
            "be positive. Ignoring value " << tol << ".\n";
        else
          m_divTol = tol;
        
        break;
      }
    case PARAM_MINBEAMS:
      {
        int beams = *((int*)value);
        
        if (beams<1 || beams>m_nLabels)
          std::cerr<<"SparseSyncMeanField::setParameters: PARAM_MINBEAMS must " <<
            "be between 1 and " << m_nLabels << ". Ignoring value " << beams << 
            ".\n";
        else
          m_minBeams = beams;
        
      break;
      }
    default:
      SyncMeanField::setParameters(param,value);
    }

}


/** Total number of non-sparse states */
int SparseSyncMeanField::getNumNonSparseStates()
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

/** Minus the log of the (reverse) cumulative sum of an array in log-space.
 *
 * In other words:
 *
 *  a[i] = -log(sum(exp(a[i:end])))
 *
 * Yes, I know noone is going to ever reuse this function, but it keeps the caller
 * tidy. :-)
 */
inline void SparseSyncMeanField::minusLogReverseCumulLogSum(FloatType* dst,
                                                            FloatType* src )
{
  register FloatType sum = 0;

  // Move to end so we can go backwards
  src = src + m_nLabels - 1;
  dst = dst + m_nLabels - 1;

  for (int i=0 ; i<m_nLabels ; i++, dst--,src--)
    {
      sum += exp(*src);
      *dst = -log(sum);
    }
}

inline FloatType SparseSyncMeanField::logSum(const FloatType a, const FloatType b)
{
  if (a>b)
    return a+log(exp(b-a)+1);
  else
    return b+log(exp(a-b)+1);
}

#ifndef NDEBUG
void SparseSyncMeanField::printAll()
{
  SyncMeanField::printAll();

  for (int i=0 ; i<m_nPixels ; i++)
    {
      *m_log << "Sparsity " << i << ":\t";
      printArray(m_nodeSparsity[i]);
      *m_log << "\n";

      *m_log << "SpExpBeliefs " << i << ":\t";
      SyncMeanField::printArray(m_spExpNodeBeliefs[i]);
      *m_log << "\n";

    }
}


void SparseSyncMeanField::printArray(bool* array)
{
  for (int i=0 ; i<m_nLabels ; i++, array++)
    *m_log << (*array?1:0) <<"\t";
  *m_log << "\n";
}

#endif
