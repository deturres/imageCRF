#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>
#include "ArrayMath.h"
#include "SparseSyncSumProd.h"

using std::endl;

/** 
 * You should be able to trace the general flow through:
 *  optimizeAlg() -> 
 *  computeMessages() -> 
 *  computeUpDownMessages(),computeLeftRightMessages() -> 
 *  computeMessage()
 *
 * The sparsity is indicated with a class-wide array for each node,
 * which is updated by this class.
 */


// $Id: SparseSyncSumProd.cpp,v 1.5 2007/12/07 15:02:03 weinman Exp $

SparseSyncSumProd::SparseSyncSumProd(int width, int height, int nLabels, 
                               EnergyFunction *eng) :
  SyncSumProd(width,height,nLabels,eng)
{
  initializeSparseSumProd();
}

SparseSyncSumProd::SparseSyncSumProd(int nPixels, int nLabels,EnergyFunction *eng) :
  SyncSumProd(nPixels,nLabels,eng)
{  
  *m_log << "Non grid graphs are not supported";
  assert(m_grid_graph);
}

SparseSyncSumProd::~SparseSyncSumProd()
{


  for (int i=0 ; i<m_nPixels ; i++) {
    delete[] m_nodeSparsity[i];
  }

  delete[] m_nodeSparsity;
}

/** Common initialization routine
 */
void SparseSyncSumProd::initializeSparseSumProd()
{

  m_divTol = 0.001;
  m_minBeams = 1;


  // Sparsity storage
  m_nodeSparsity = new bool*[m_nPixels];

  // Create sparsity terms
  for (int i=0 ; i<m_nPixels ; i++)
    m_nodeSparsity[i] = new bool[m_nLabels];


}

/** Initialization from parent classes */
void SparseSyncSumProd::initializeAlg()
{
  SyncSumProd::initializeAlg();

  // Initialize all sparsities according to initial beliefs
  //updateSparsity(); // doing this seems to worsen results

   // Initialize all sparsities to true (since all are initially included)
   for (int i=0 ; i<m_nPixels ; i++)
     for (int b=0 ; b<m_nLabels ; b++)
       m_nodeSparsity[i][b] = true;

}

/** Run the algorithm (sparse synchronous loopy sum product).
 */
void SparseSyncSumProd::optimizeAlg (int nIterations)
{

  FloatType diff;

  for (int i=0 ; i<nIterations ; i++)
    {
      computeMessages();

      diff = NodeMessages::messageDifference (m_messageDiffFun, 
                                              m_nodeMessages,
                                              m_nPixels);
      
      *m_log << "SparseSyncSumProd: Iter = " << i << "\tDiff = " << diff << '\n';

      NodeMessages::updateMessages (m_messageDamper, m_nodeMessages, m_nPixels);

      updateSparsity();

      if (diff < m_messageDiffTol )
        {
          *m_log << "SparseSyncSumProd converged in " << i << " iterations.\n";
          break;
        }
    }

  updateAnswer();

}

/** Compute all messages
 */
void SparseSyncSumProd::computeMessages ()
{

  computeLeftRightMessages();
  computeUpDownMessages();

}

/** Compute messages for all rows from left to right and vice-versa.
 */
void SparseSyncSumProd::computeLeftRightMessages() 
{

  int pixel;
  NodeMessages** msg;
  bool** sp;

  for (int r=0 ; r<m_height ; r++)
    {

      pixel = r*m_width;
      msg = m_nodeMessages + pixel;
      sp = m_nodeSparsity + pixel;
      
      for (int c=0 ; c<m_width-1; c++, pixel++, msg++,sp++)
        {
          // Left-to-Right
          computeMessage(pixel, pixel+1, *msg, *(msg+1), *sp, 
                         NodeMessages::MSG_LEFT);
          
          // Right-to-Left
          computeMessage(pixel+1, pixel, *(msg+1), *msg, *(sp+1),
                         NodeMessages::MSG_RIGHT);
        }
    }
}

/** Compute messages for all columns from up to down and vice-versa.
 */
void SparseSyncSumProd::computeUpDownMessages()
{

  int pixel,pixel1; // index for base, up, down 
  NodeMessages** msg; // message for up (top)
  NodeMessages** msg1; // message for down (bottom)

  bool** sp; // sparsity of top
  bool** sp1; // sparsity of bottom

  for (int r=0 ; r<m_height-1 ; r++)
    {

      pixel = r*m_width;
      pixel1 = pixel + m_width;

      msg = m_nodeMessages + pixel;
      msg1 = m_nodeMessages + pixel1;
      
      sp = m_nodeSparsity + pixel;
      sp1 = m_nodeSparsity + pixel1;

      for (int c=0 ; c<m_width; c++ , pixel++,pixel1++,msg++,msg1++,sp++,sp1++)
        {
          // Up-to-Down
          computeMessage(pixel, pixel1, *msg, *msg1, *sp, NodeMessages::MSG_UP);
          
          // Down-to-Up
          computeMessage(pixel1, pixel, *msg1, *msg, *sp1, NodeMessages::MSG_DOWN);
        }
    }
}

    
/** Compute the message between two nodes
 *
 * Parameter pixel1 is source node and pixel2 is the receiving node.
 *
 * Parameter direction is 
 *    - from pixel2 (dst) to pixel1 (src)
 *    - counter-intuitive?
 *    - and redundant
 *    but all of this should save some computation.
 * 
 */
void SparseSyncSumProd::computeMessage(const int pixel1, const int pixel2,
                                       NodeMessages* message1, 
                                       NodeMessages* message2,
                                       bool* sparsity, int direction)
{

  // NOTE -- could cache these two arrays between calls (i.e. pass them as 
  // arguments) rather than redeclaring each time. Maybe would save time?

  // Summand array for the source node
  CostVal *summand = new CostVal[m_nLabels];
  assert(summand!=0);

  // Data cost and message product at source node
  CostVal *dataCostMsgProd = new CostVal[m_nLabels];
  assert(dataCostMsgProd!=0);

  // Start with data cost of pixel1
  copySparseDataCost(pixel1,dataCostMsgProd, sparsity);  

  // Add message product of pixel1
  addToSparse(dataCostMsgProd,message1->getProduct(),sparsity);

  // Subtract the message from pixel2 to pixel1
  subtractFromSparse(dataCostMsgProd, 
                     message1->getMessage (NodeMessages::reverseDir(direction)),
                     sparsity);



  // Space for storing the updated message
  FloatType* dst = message2->getNewMessage(direction);
  FloatType* dst0 = dst;

  FloatType logsum;
  /* Here's the idea: 
   *   For each label at pixel2
   *     We will sum over labels for pixel1 a term that contains:
   *       D(pixel1,label1) + V(pixel1,pixel2,label1,label2) +
   *       sum_{i:pixel1 neighbors i} Message[from i to pixel1](label1) -
   *       Message[from pixel2 to pixel1](label1)
   */
  for (Label label2=0 ; label2<m_nLabels ; label2++, dst++)
    {
      // Copy the data cost and message product at pixel1
      copyToSparse(summand, dataCostMsgProd, sparsity);
      
      // Add smoothness terms (pixel2 is fixed at label2, 
      // pixel1 varies over summand array)
      addSparseSmoothnessCost(pixel2, pixel1,label2, summand, sparsity);

      // Finally, sum the summand vector (which is in log-space)
      logsum = logSumSparse(summand,sparsity);

      if (false)//isinf(logsum))
        {
          std::cerr<<"WARNING: ["<<pixel1<<"->"<<pixel2<<" yields -Inf. ";
          int ss=0;
          bool* ps=sparsity;
          for (int b=0  ; b<m_nLabels ; b++,ps++) 
            if (*ps) ss++;

          *m_log<<"|sparsity| = "<<ss<<"\n";

        }

      *dst = logsum;
    }

  // Normalize the message
  logNormalize(dst0);
  delete[] summand;
  delete[] dataCostMsgProd;

}

/** Updates sparsities based on current approximate marginals */
void SparseSyncSumProd::updateSparsity()
{
  for (int i=0 ; i<m_nPixels ; i++)
    updateSparsity(i);
}

void SparseSyncSumProd::updateSparsity(int pixel)
{
  FloatType *belief = new FloatType[m_nLabels];
  FloatType *belSort = new FloatType[m_nLabels];

  int kBeam = 1;   // Index of cut-off

  FloatType kThresh; // Belief threshold

  bool* pSparse;
  FloatType* pBelief;


  getBelief(pixel,belief);
  
  copyTo(belSort,belief);
  
  // Sorts in ascending order (doh!)
  std::sort(belSort,belSort+m_nLabels);
  
  // Find index of cut-off: largest index (from end)
  // where divergence is less than the threshold
  FloatType *bel = belSort+m_nLabels-1;
  
  FloatType logcpf = *bel;

  bel--;

  for ( kBeam = 0; 
        kBeam < m_nLabels && kBeam <= m_nLabels-m_minBeams && -logcpf >= m_divTol ;
        kBeam++, bel-- )
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
    } 
  else
    {
      kThresh = belSort[m_nLabels-1-kBeam];
      
      pBelief = belief;
      
      
      for (int b=0 ; b<m_nLabels ; b++, pSparse++, pBelief++ )
        *pSparse = (*pBelief >= kThresh);
    }
  delete[] belief;
  delete[] belSort;
}


/** Parameter settings
 * 
 * See PARAM_* macros in SyncSumProd.h for options
 */
void SparseSyncSumProd::setParameters (int param, void* value)
{
  
  switch (param)
    {
    case PARAM_DIVTOL:
      {
        FloatType tol = *((FloatType*)value);
        
        if (tol<0)
          std::cerr<<"SparseSyncSumProd::setParameters: PARAM_DIVTOL must " <<
            "be positive. Ignoring value " << tol << ".\n";
        else
          m_divTol = tol;
        
        break;
      }
    case PARAM_MINBEAMS:
      {
        int beams = *((int*)value);
        
        if (beams<1 || beams>m_nLabels)
          std::cerr<<"SparseSyncSumProd::setParameters: PARAM_MINBEAMS must" <<
            " be between 1 and " << m_nLabels << ". Ignoring value " << beams <<
            ".\n";
        else
          m_minBeams = beams;
        
      break;
      }
    default:
      SyncSumProd::setParameters(param,value);
    }

}


/** Total number of non-sparse states */
int SparseSyncSumProd::getNumNonSparseStates()
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

/** Total number of non-sparse states */
int SparseSyncSumProd::getNumNonSparseStates(const int pixel)
{
  
  register int sum = 0;
  bool* ps = m_nodeSparsity[pixel];

  for (int b=0 ; b<m_nLabels ; b++, ps++)
    sum += (*ps ? 1 : 0);
  
  return sum;
}

/** Maximum number of non-sparse states */
int SparseSyncSumProd::getMaxNumNonSparseStates()
{
  bool** ps = m_nodeSparsity;
  bool* psn;
  int mx = 0;
  register int sum = 0;
  for (int i=0 ; i<m_nPixels ; i++, ps++)
    {
      psn = *ps;
      sum = 0;

      for (int b=0 ; b<m_nLabels ; b++,psn++)
        sum += (*psn ? 1 : 0);

      if (sum>mx) mx = sum;
    }
  
  return mx;
}

/** Print a histogram of the number of states at the pixels */
void SparseSyncSumProd::printHistSparsity(int numBins)
{

  register int sum=0;
  bool** ps = m_nodeSparsity;
  bool* psn;

  int *hist = new int[numBins+1];

  for (int c=0 ; c<numBins+1; c++)
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

   for (int c=0 ; c<numBins+1 ; c++)
     *m_log<<hist[c]<<"\t";
   *m_log<<"\n";

   delete[] hist;
  
}

/** Adds two numbers stored in logspace in a numerically stable fashion */
inline FloatType SparseSyncSumProd::logSum(const FloatType a, const FloatType b)
{
  if (a>b)
    return a+log(exp(b-a)+1);
  else
    return b+log(exp(a-b)+1);
}
