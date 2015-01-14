#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "AsyncSumProd.h"

using std::endl;

/** 
 * You should be able to trace the general flow through:
 *  optimizeAlg() -> 
 *  computeMessages() -> 
 *  computeUpDownMessages(),computeLeftRightMessages() -> 
 *  computeMessage()
 */

// $Id: AsyncSumProd.cpp,v 1.1.1.1 2007/12/01 01:01:14 cpal Exp $


AsyncSumProd::AsyncSumProd(int width, int height, int nLabels, EnergyFunction *eng) :
  SyncSumProd(width,height,nLabels,eng)
{ }

AsyncSumProd::AsyncSumProd(int nPixels, int nLabels,EnergyFunction *eng) :
  SyncSumProd(nPixels,nLabels,eng)
{  
  *m_log << "Non grid graphs are not supported";
  assert(m_grid_graph);
}

AsyncSumProd::~AsyncSumProd() {}


/** Run the algorithm (asynchronous loopy sum product).
 *
 * Although each direction is computed asynchronously, they are always computed
 * in the same order, so the schedule is not as random as it could be.
 */
void AsyncSumProd::optimizeAlg (int nIterations)
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

      
      *m_log << "AsyncSumProd: Iter = " << i << "\tDiff = " << diff << '\n';

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
FloatType AsyncSumProd::computeLeftMessages() 
{

  int pixel;
  NodeMessages** msg;
  
  FloatType diff = 0;
  FloatType diff0;

  for (int r=0 ; r<m_height ; r++)
    {

      pixel = r*m_width;
      msg = m_nodeMessages + pixel;

      for (int c=0 ; c<m_width-1; c++, pixel++, msg++)
        {
          // Calculate new message
          computeMessage(pixel, pixel+1, *msg, *(msg+1), NodeMessages::MSG_LEFT);
          
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
FloatType AsyncSumProd::computeRightMessages() 
{

  int pixel;
  NodeMessages** msg;

  FloatType diff = 0;
  FloatType diff0;
  
  for (int r=0 ; r<m_height ; r++)
    {

      pixel = (r+1)*m_width - 2;
      msg = m_nodeMessages + pixel;

      for (int c=m_width-1 ; c>0; c--, pixel--, msg--)
        {
          // Calculate new message
          computeMessage(pixel+1, pixel, *(msg+1), *msg, NodeMessages::MSG_RIGHT);

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
FloatType AsyncSumProd::computeUpMessages()
{

  int pixel,pixel1; // index for base, up, down 
  NodeMessages** msg; // message for up (top)
  NodeMessages** msg1; // message for down (bottom)

  FloatType diff = 0;
  FloatType diff0;

  for (int r=0 ; r<m_height-1 ; r++)
    {

      pixel = r*m_width;
      pixel1 = pixel + m_width;
      
      msg = m_nodeMessages + pixel;
      msg1 = m_nodeMessages + pixel1;
      
      for (int c=0 ; c<m_width; c++, pixel++, pixel1++,msg++,msg1++)
        {
          // Calculate new message
          computeMessage(pixel, pixel1, *msg, *msg1, NodeMessages::MSG_UP);

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
FloatType AsyncSumProd::computeDownMessages()
{

  int pixel,pixel1; // index for base, up, down 
  NodeMessages** msg; // message for up (top)
  NodeMessages** msg1; // message for down (bottom)

  FloatType diff = 0;
  FloatType diff0;

  for (int r=m_height-2 ; r>=0 ; r--)
    {

      pixel = r*m_width;
      pixel1 = pixel + m_width;
      
      msg = m_nodeMessages + pixel;
      msg1 = m_nodeMessages + pixel1;
      
      for (int c=0 ; c<m_width; c++, pixel++, pixel1++,msg++,msg1++)
        {
          // Calculate new message
          computeMessage(pixel1, pixel, *msg1, *msg, NodeMessages::MSG_DOWN);

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
