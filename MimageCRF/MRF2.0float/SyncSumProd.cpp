#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "ArrayMath.h"
#include "SyncSumProd.h"

using std::endl;

/** 
 * You should be able to trace the general flow through:
 *  optimizeAlg() -> 
 *  computeMessages() -> 
 *  computeUpDownMessages(),computeLeftRightMessages() -> 
 *  computeMessage()
 */

// $Id: SyncSumProd.cpp,v 1.1.1.1 2007/12/01 01:01:15 cpal Exp $



SyncSumProd::SyncSumProd(int width, int height, int nLabels, EnergyFunction *eng) :
  MRFEnergy(width,height,nLabels,eng)
{
  initializeSumProd();
}

SyncSumProd::SyncSumProd(int nPixels, int nLabels,EnergyFunction *eng) :
  MRFEnergy(nPixels,nLabels,eng)
{  
  std::cerr << "Non grid graphs are not supported";
  assert(m_grid_graph);
}

SyncSumProd::~SyncSumProd()
{
  
  for (int i=0 ; i<m_nPixels ; i++)
    delete m_nodeMessages[i];
  
  delete[] m_nodeMessages;

}

void SyncSumProd::initializeSumProd()
{

  m_messageDamper = 1;//0.8;
  m_messageDiffFun = NodeMessages::L2;
  m_messageDiffTol = 1;//.001;

  // Message storage
  m_nodeMessages = new NodeMessages*[m_nPixels];

  //
  // Create messages
  //

  int nbr;

  for (int r=0; r<m_height ; r++)
    for (int c=0 ; c<m_width ; c++)
      {
        // Neighbor indicator
        nbr = 
          (r>0          ? NodeMessages::NBR_UP    : 0) | 
          (r<m_height-1 ? NodeMessages::NBR_DOWN  : 0) | 
          (c>0          ? NodeMessages::NBR_LEFT  : 0) | 
          (c<m_width-1  ? NodeMessages::NBR_RIGHT : 0);

        // Create node message struct 
        m_nodeMessages[pixelIndex(r,c)] = new NodeMessages ( m_nLabels, nbr );
      }

}

/** Initialization from parent classes */
void SyncSumProd::initializeAlg()
{

  // Reset all the messages and products
  NodeMessages** msg = m_nodeMessages;

  for (int i=0 ; i<m_nPixels ; i++,msg++)
    {
      (*msg)->resetMessages();
      (*msg)->resetProduct();
    }

  updateAnswer(); // Initialize to MPM

}

/** Run the algorithm (synchronous loopy sum product).
 */
void SyncSumProd::optimizeAlg (int nIterations)
{

  FloatType diff;

  for (int i=0 ; i<nIterations ; i++)
    {
      computeMessages();

      diff = NodeMessages::messageDifference (m_messageDiffFun, 
                                              m_nodeMessages,
                                              m_nPixels);
      
      *m_log << "SyncSumProd: Iter = " << i << "\tDiff = " << diff << '\n';

      NodeMessages::updateMessages (m_messageDamper, m_nodeMessages, m_nPixels);
      
      if (diff < m_messageDiffTol )
        {
          *m_log << "SyncSumProd converged in " << i << " iterations.\n";
          break;
        }
    }

  updateAnswer();

}

/** Compute all messages
 */
void SyncSumProd::computeMessages ()
{

  computeLeftRightMessages();
  computeUpDownMessages();

}

/** Compute messages for all rows from left to right and vice-versa.
 */
void SyncSumProd::computeLeftRightMessages() 
{

  int pixel;
  NodeMessages** msg;

  for (int r=0 ; r<m_height ; r++)
    {

      pixel = r*m_width;
      msg = m_nodeMessages + pixel;

      for (int c=0 ; c<m_width-1; c++, pixel++, msg++)
        {
          // Left-to-Right
          computeMessage(pixel, pixel+1, *msg, *(msg+1), NodeMessages::MSG_LEFT);
          
          // Right-to-Left
          computeMessage(pixel+1, pixel, *(msg+1), *msg, NodeMessages::MSG_RIGHT);
        }
    }
}

/** Compute messages for all columns from up to down and vice-versa.
 */
void SyncSumProd::computeUpDownMessages()
{

  int pixel,pixel1; // index for base, up, down 
  NodeMessages** msg; // message for up (top)
  NodeMessages** msg1; // message for down (bottom)

  for (int r=0 ; r<m_height-1 ; r++)
    {

      pixel = r*m_width;
      pixel1 = pixel + m_width;
      
      msg = m_nodeMessages + pixel;
      msg1 = m_nodeMessages + pixel1;
      
      for (int c=0 ; c<m_width; c++, pixel++, pixel1++,msg++,msg1++)
        {

          // Up-to-Down
          computeMessage(pixel, pixel1, *msg, *msg1, NodeMessages::MSG_UP);
          
          // Down-to-Up
          computeMessage(pixel1, pixel, *msg1, *msg, NodeMessages::MSG_DOWN);
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
void SyncSumProd::computeMessage(const int pixel1, const int pixel2,
                                 NodeMessages* message1, NodeMessages* message2,
                                 int direction)
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
  copyDataCost(pixel1,dataCostMsgProd);  

  // Add message product of pixel1
  addTo(dataCostMsgProd,message1->getProduct());

  // Subtract the message from pixel2 to pixel1
  subtractFrom(dataCostMsgProd, 
               message1->getMessage (NodeMessages::reverseDir(direction)));



  // Space for storing the updated message
  FloatType* dst = message2->getNewMessage(direction);
  FloatType* dst0 = dst;

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
      copyTo(summand,dataCostMsgProd);
      
      // Add smoothness terms (pixel2 is fixed at label2, 
      // pixel1 varies over summand array)
      addSmoothnessCost(pixel2, pixel1,label2,summand);

      // Finally, sum the summand vector (which is in log-space)
      *dst = logSum(summand);
    }

  // Normalize the message
  logNormalize(dst0);
  delete[] summand;
  delete[] dataCostMsgProd;

}

/** Calculate the current approximate marginal for a pixel and write
 * it to the destination array provided (must be of length m_nlabels.
 *
 * Note the answer is returned in logspace
 *
 */
void SyncSumProd::getBelief(const int pixel, FloatType* dst)
{
  // b_i(p) \propto exp(-D_i(p)) *  prod_j M_ij(p)

  copyDataCost(pixel,dst);
  addTo(dst, m_nodeMessages[pixel]->getProduct());
  logNormalize(dst);

}

FloatType SyncSumProd::getBelief(const int pixel, const Label label)
{
  std::cerr << "SyncSumProd::getBelief(pixel,label) not implemented." << endl;
  assert(false);
  return(0.0);
}

/** Calculate the current approximate marginal for a pair of neighboring pixels
 *  and write it to the destination array provided (must be of length m_nlabels^2.
 *
 *  Note: pixel1 should be above or left of pixel2. Then pixel1 is the row and
 *        pixel2 is the column in a row-major order of the destination array.
 */
void SyncSumProd::getEdgeBelief(const int pixel1, const int pixel2, FloatType* dst)
{

  // Calculate relative orientation of the nodes
  int direction1;
  int direction2;

  FloatType *dataCostMsgProd = new FloatType[m_nLabels];
  assert(dataCostMsgProd!=0);

  if (pixel1+1==pixel2 && m_width>1)
    {
      direction1 = NodeMessages::MSG_RIGHT;
      direction2 = NodeMessages::MSG_LEFT;
    }
  else
    {
      direction1 = NodeMessages::MSG_DOWN;
      direction2 = NodeMessages::MSG_UP;
    }

  //
  // Pixel 2 - Message + Data Cost
  //

  // Start with data cost of pixel2
  copyDataCost(pixel2, dataCostMsgProd);

  // Add product of messages into pixel2
  addTo (dataCostMsgProd, m_nodeMessages[pixel2]->getProduct());
  
  // Remove pixel1's message
  subtractFrom (dataCostMsgProd, m_nodeMessages[pixel2]->getMessage(direction2));

  // Repeat the one dimensional D(p2)+M(p2) over the ROWS of the destination
  FloatType* pd = dst;
  for (int b1=0 ; b1<m_nLabels; b1++, pd+=m_nLabels)
    copyTo(pd,dataCostMsgProd);
  

  // 
  // Pixel 1
  //

  // Start with data cost of pixel1
  copyDataCost(pixel1, dataCostMsgProd);

  // Product of messages into pixel1
  addTo (dataCostMsgProd, m_nodeMessages[pixel1]->getProduct());
  
  // Remove pixel2's message
  subtractFrom (dataCostMsgProd, m_nodeMessages[pixel1]->getMessage(direction1));
  
  // Repeat (add) the one dimensional D(p1)+M(p1) over the COLUMNS
  FloatType* ps = dataCostMsgProd;
  pd = dst;
  for (int b1=0 ; b1<m_nLabels ; b1++, ps++)
      for (int b2=0 ; b2<m_nLabels ; b2++, pd++)
        *pd += *ps;
  

  // Smoothness Cost
  addSmoothnessCost (pixel1,pixel2,dst);

  // Normalize
  logNormalizeEdge (dst);
  delete[] dataCostMsgProd;
  
}

/** Calculate the log probability of a labeling from the current beliefs
 *
 * Note: assumes labels is of length m_nPixels
 *  See Equation (39) in
 *
 * JS Yedidia, WT Freeman, Y Weiss. Constructing Free Energy
 * Approximations and Generalized Belief Propagation Algorithms. IEEE
 * Transactions on Information Theory, 2005
 */
FloatType SyncSumProd::getLogProb(Label* labels)
{
  
  assert(m_grid_graph);

  register FloatType logprob = 0;

  FloatType degree;

  FloatType *nodeBelief = new FloatType[m_nLabels];
  assert(nodeBelief!=0);

  FloatType *edgeBelief = new FloatType[m_nLabels*m_nLabels];
  assert(edgeBelief!=0);

  int pixC, pixR, pixD; 
  for (int r=0 ; r<m_height ; r++)
    for (int c=0 ; c<m_width ; c++)
      {

        pixC = pixelIndex(r,c);
        pixD = pixelIndex(r+1,c);
        pixR = pixelIndex(r,c+1);

        degree = 
          (r>0          ? 1 : 0 ) + // Up neighbor
          (r<m_height-1 ? 1 : 0 ) + // Down neighbor
          (c>0          ? 1 : 0 ) + // Left neighbor
          (c<m_width-1  ? 1 : 0 );  // Right neighbor


        getBelief(pixC, nodeBelief);

        // Subtract node term
        logprob -= (degree-1) * nodeBelief[labels[pixC]];

        // Add down edge
        if (r<m_height-1)
          {
            getEdgeBelief(pixC,pixD, edgeBelief);
            
            logprob += edgeBelief[labelIndex(labels[pixC],labels[pixD])];
          }

        // Add right edge
        if (c<m_width-1)
          {
            getEdgeBelief(pixC,pixR, edgeBelief);
            
            logprob += edgeBelief[labelIndex(labels[pixC],labels[pixR])];

          }

      }
  delete[] nodeBelief;
  delete[] edgeBelief;
  return logprob;
}

/** Total entropy of node marginals */
FloatType SyncSumProd::getTotalNodeEntropy()
{
  FloatType* belief;
  FloatType *belief0 = new FloatType[m_nLabels];
  assert(belief0 != 0);

  register FloatType entropy = 0;

  for (int i=0 ; i<m_nPixels ; i++)
    {
      getBelief(i,belief0);
      belief = belief0;
      for (int b=0 ; b<m_nLabels ; b++, belief++)
        entropy -= exp(*belief) * (*belief);
    }
  delete[] belief0;
  return entropy;
}

/** Sets answer  to current Maximum Posterior Marginal (MPM) estimate
 */
void SyncSumProd::updateAnswer()
{

  FloatType* belief = new FloatType[m_nLabels];

  Label* label = m_labels;

  for (int i=0; i<m_nPixels ; i++, label++)
    {
      getBelief (i, belief);
      *label = argMax (belief);
    }

}


/** Parameter settings
 * 
 * See PARAM_* macros in SyncSumProd.h for options
 */
void SyncSumProd::setParameters (int param, void* value)
{
  
  switch (param)
    {

    case PARAM_DAMPER:
      {
      FloatType damp = *((FloatType*)value);
      
      if (damp>1 || damp<0)
        std::cerr << "SyncSumProd::setParameters: PARAM_DAMPER must be " <<
          "between zero and one. Ignoring.\n";
      else
        m_messageDamper = damp;

      break;
      }

    case PARAM_MSGDFUN:
      m_messageDiffFun = *((NodeMessages::DiffType*)value);
      break;

    case PARAM_MSGDTOL:
      m_messageDiffTol = *((FloatType*)value);
      break;
    default:
      std::cerr << "SyncSumProd::setParameters: Unknown parameter index " <<
        param << '\n';
    }

}


/** Add an edge to the graph.
 *
 * Note: Required by MRF, but NOT SUPPORTED
 */
void SyncSumProd::setNeighbors (int pix1, int pix2, CostVal weight)
{
  std::cerr << "SyncSumProd::setNeighbors not supported\n";
  assert(false);
}
