#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "SparseAsyncSumProd.h"

/** 
 * You should be able to trace the general flow through:
 *  optimizeAlg() -> 
 *  computeMessages() -> 
 *  compute[Down,Left,Up,Right]Messages() ->
 *  computeMessage()
 *
 * The sparsity is indicated with a class-wide array for each node,
 * which is updated by this class.
 */


using std::endl;

// $Id: SparseAsyncSumProd.cpp,v 1.1.1.1 2007/12/01 01:01:15 cpal Exp $

SparseAsyncSumProd::SparseAsyncSumProd(int width, int height, int nLabels, 
                               EnergyFunction *eng) :
  SparseSyncSumProd(width,height,nLabels,eng)
{ }

SparseAsyncSumProd::SparseAsyncSumProd(int nPixels, int nLabels,EnergyFunction *eng) :
  SparseSyncSumProd(nPixels,nLabels,eng)
{  
  *m_log << "Non grid graphs are not supported";
  assert(m_grid_graph);
}

SparseAsyncSumProd::~SparseAsyncSumProd() { }

/** Run the algorithm (asynchronous loopy sum product).
 */
void SparseAsyncSumProd::optimizeAlg (int nIterations)
{

  FloatType diff, diff0;

  for (int i=0 ; i<nIterations ; i++)
    {
      
      diff = 0;

      diff0 = computeDownMessages();
      diff = NodeMessages::combineDifference( m_messageDiffFun, diff, diff0);


      diff0 = computeLeftMessages();
      diff = NodeMessages::combineDifference( m_messageDiffFun, diff, diff0);

      diff0 = computeUpMessages();
      diff = NodeMessages::combineDifference( m_messageDiffFun, diff, diff0);

      diff0 = computeRightMessages();
      diff = NodeMessages::combineDifference( m_messageDiffFun, diff, diff0);

      
      *m_log << "SparseAsyncSumProd: Iter = " << i << "\tDiff = " << diff << '\n';

      updateSparsity(); // Here or after each update?

      if (diff < m_messageDiffTol )
        {
          *m_log << "AsyncSumProd converged in " << i << " iterations.\n";
          break;
        }
    }

  updateAnswer();

}


/** Compute and update messages for all rows from left to right
 *
 * Returns total message difference
 */
FloatType SparseAsyncSumProd::computeLeftMessages() 
{

  int pixel;
  NodeMessages** msg;
  bool** sp;

  FloatType diff = 0;
  FloatType diff0;

  for (int r=0 ; r<m_height ; r++)
    {

      pixel = r*m_width;
      msg = m_nodeMessages + pixel;
      sp = m_nodeSparsity + pixel;

      for (int c=0 ; c<m_width-1; c++, pixel++, msg++, sp++)
        {
          // Calculate new message
          computeMessage(pixel, pixel+1, *msg, *(msg+1), *sp, 
                         NodeMessages::MSG_LEFT);
          
          // Calculate the message difference
          diff0 = (*(msg+1))->messageDifference(m_messageDiffFun,
                                                NodeMessages::MSG_LEFT);

          diff = NodeMessages::combineDifference(m_messageDiffFun, diff0, diff);
          
          // Update message
          (*(msg+1))->updateMessage(NodeMessages::MSG_LEFT,m_messageDamper);
          
        }
    }

  return diff;
}


/** Compute and update messages for all rows from right to left
 *
 * Returns total message difference
 */
FloatType SparseAsyncSumProd::computeRightMessages() 
{

  int pixel;
  NodeMessages** msg;
  bool** sp;

  FloatType diff = 0;
  FloatType diff0;
  
  for (int r=0 ; r<m_height ; r++)
    {

      pixel = (r+1)*m_width - 2;
      msg = m_nodeMessages + pixel;
      sp = m_nodeSparsity + pixel;

      for (int c=m_width-1 ; c>0; c--, pixel--, msg--, sp--)
        {
          // Calculate new message
          computeMessage(pixel+1, pixel, *(msg+1), *msg, *(sp+1),
                         NodeMessages::MSG_RIGHT);

          // Calculate the message difference
          diff0 = (*msg)->messageDifference(m_messageDiffFun,
                                            NodeMessages::MSG_RIGHT);
          
          diff = NodeMessages::combineDifference(m_messageDiffFun, diff0, diff);
          
          // Update message
          (*msg)->updateMessage(NodeMessages::MSG_RIGHT,m_messageDamper);

        }
    }
  return diff;
}

/** Compute and update messages for all columns from up to down
 *
 * Returns total message difference
 */
FloatType SparseAsyncSumProd::computeUpMessages()
{

  int pixel,pixel1; // index for base, up, down 

  NodeMessages** msg; // message for up (top)
  NodeMessages** msg1; // message for down (bottom)

  bool** sp; // sparsity of top

  FloatType diff = 0;
  FloatType diff0;

  for (int r=0 ; r<m_height-1 ; r++)
    {

      pixel = r*m_width;
      pixel1 = pixel + m_width;
      
      msg = m_nodeMessages + pixel;
      msg1 = m_nodeMessages + pixel1;
      
      sp = m_nodeSparsity + pixel;

      for (int c=0 ; c<m_width; c++, pixel++, pixel1++,msg++,msg1++,sp++)
        {
          // Calculate new message
          computeMessage(pixel, pixel1, *msg, *msg1, *sp, NodeMessages::MSG_UP);

          // Calculate the message difference
          diff0 = (*msg1)->messageDifference(m_messageDiffFun,
                                             NodeMessages::MSG_UP);
          
          diff = NodeMessages::combineDifference(m_messageDiffFun, diff0, diff);
          
          // Update message
          (*msg1)->updateMessage(NodeMessages::MSG_UP,m_messageDamper);

        }
    }
  return diff;

}

/** Compute and update messages for all columns from down to up
 *
 * Returns total message difference
 */
FloatType SparseAsyncSumProd::computeDownMessages()
{

  int pixel,pixel1; // index for base, up, down 
  NodeMessages** msg; // message for up (top)
  NodeMessages** msg1; // message for down (bottom)

  FloatType diff = 0;
  FloatType diff0;

  bool** sp1; // sparsity of bottom

  for (int r=m_height-2 ; r>=0 ; r--)
    {

      pixel = r*m_width;
      pixel1 = pixel + m_width;
      
      msg = m_nodeMessages + pixel;
      msg1 = m_nodeMessages + pixel1;
     
      sp1 = m_nodeSparsity + pixel1;

      for (int c=0 ; c<m_width; c++, pixel++, pixel1++,msg++,msg1++,sp1++)
        {
          // Calculate new message
          computeMessage(pixel1, pixel, *msg1, *msg, *sp1, 
                         NodeMessages::MSG_DOWN);

          // Calculate the message difference
          diff0 = (*msg)->messageDifference(m_messageDiffFun,
                                            NodeMessages::MSG_DOWN);
          
          diff = NodeMessages::combineDifference(m_messageDiffFun, diff0, diff);
          
          // Update message
          (*msg)->updateMessage(NodeMessages::MSG_DOWN,m_messageDamper);


        }
    }

  return diff;
}
