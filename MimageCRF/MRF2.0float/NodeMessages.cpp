#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <assert.h>
#include "NodeMessages.h"
#include "ArrayMath.h"

#define MAX(a,b)   (a>=b ? a : b)


int const NodeMessages::MSG_UP    = 0;
int const NodeMessages::MSG_DOWN  = 1;
int const NodeMessages::MSG_LEFT  = 2;
int const NodeMessages::MSG_RIGHT = 3;

int const NodeMessages::NBR_UP    = 1;
int const NodeMessages::NBR_DOWN  = 2;
int const NodeMessages::NBR_LEFT  = 4;
int const NodeMessages::NBR_RIGHT = 8;

NodeMessages::NodeMessages(int nLabels, 
                           int neighbors)
{


  if (nLabels<1)
    {
      std::cerr << "There must be a positive number of labels for nodes";
      assert(nLabels>0);
    }



  m_nLabels = nLabels;
  m_neighbors = neighbors;

  m_product = new FloatType[nLabels];
  
  if (!m_product)
    {
      std::cerr << "Insufficient memory for messages\n";
      exit(EXIT_FAILURE);
    }

  // Initialize messages and updated mesages
  if (NBR_UP & m_neighbors) {
    m_messages[MSG_UP] = new FloatType[nLabels];
    m_newMessages[MSG_UP] = new FloatType[nLabels];

    if (!m_messages[MSG_UP] || ! m_newMessages[MSG_UP])
      {
        std::cerr << "Insufficient memory for messages\n";
        exit(EXIT_FAILURE);
      }

  }

  if (NBR_DOWN & m_neighbors) {
    m_messages[MSG_DOWN] = new FloatType[nLabels];
    m_newMessages[MSG_DOWN] = new FloatType[nLabels];

    if (!m_messages[MSG_DOWN] || ! m_newMessages[MSG_DOWN])
      {
        std::cerr << "Insufficient memory for messages\n";
        exit(EXIT_FAILURE);
      }

  }

  if (NBR_LEFT & m_neighbors) {
    m_messages[MSG_LEFT] = new FloatType[nLabels];
    m_newMessages[MSG_LEFT] = new FloatType[nLabels];

    if (!m_messages[MSG_LEFT] || ! m_newMessages[MSG_LEFT])
      {
        std::cerr << "Insufficient memory for messages\n";
        exit(EXIT_FAILURE);
      }

  }
    
  if (NBR_RIGHT & m_neighbors) {
    m_messages[MSG_RIGHT] = new FloatType[nLabels];
    m_newMessages[MSG_RIGHT] = new FloatType[nLabels];

    if (!m_messages[MSG_RIGHT] || ! m_newMessages[MSG_RIGHT])
      {
        std::cerr << "Insufficient memory for messages\n";
        exit(EXIT_FAILURE);
      }

  }

  resetMessages();
  resetProduct();
  
}

NodeMessages::~NodeMessages()
{
  // Free messages and updated mesages
  if (NBR_UP & m_neighbors) {
    delete[] m_messages[MSG_UP];
    delete[] m_newMessages[MSG_UP];
  }

  if (NBR_DOWN & m_neighbors) {
    delete[] m_messages[MSG_DOWN];
    delete[] m_newMessages[MSG_DOWN];
  }

  if (NBR_LEFT & m_neighbors) {
    delete[] m_messages[MSG_LEFT];
    delete[] m_newMessages[MSG_LEFT];
  }
    
  if (NBR_RIGHT & m_neighbors) {
    delete[] m_messages[MSG_RIGHT];
    delete[] m_newMessages[MSG_RIGHT];
  }

  delete[] m_product;
}


/** Copy the new message for a direction into the current message.
 *
 * Note that cached products are updated.
 */
void NodeMessages::updateMessage(int direction)
{

  // Remove current message from product
  subtractFrom(m_product,m_messages[direction]);

  // Swap-out
  FloatType* tmp;

  tmp = m_messages[direction];
  m_messages[direction] = m_newMessages[direction];
  m_newMessages[direction] = tmp;

  // Add new message to product
  addTo(m_product,m_messages[direction]);
  
}



/** Convex combination of new messages with current messages. 
 *  Since this is an update the result is stored in the current message field.
 *
 *  Parameter damper is the fraction of new message to use.
 *  Note that cached products are updated.
 */
void NodeMessages::updateMessage(int direction, FloatType damper)
{

  if (damper==1)
    {
      updateMessage(direction);
      return;
    }
  
  // Remove current message from product
  subtractFrom(m_product,m_messages[direction]);

  // Update message
  convexLogAdd(damper, m_messages[direction], m_newMessages[direction]);

  // Add new message to product
  addTo(m_product,m_messages[direction]);
  

}
        
/** Copy the new messages into the current messages.
 *
 * This actually just swaps the pointers to re-use the space.
 * Note that cached products are updated.
 */
void NodeMessages::updateMessages()
{  

  FloatType* tmp;
  for (int i=0 ; i<4 ; i++)
    {
      tmp = m_messages[i];
      m_messages[i] = m_newMessages[i];
      m_newMessages[i] = tmp;
    }
  
  updateProduct();
}


/** Convex combination of new messages with current messages. 
 *  Since this is an update the result is stored in the current message field.
 *
 *  Parameter damper is the fraction of new message to use.
 *  Note that cached products are updated.
 */
void NodeMessages::updateMessages(FloatType damper)
{
  if (damper==1)
    {
      updateMessages();
      return;
    }

  if (NBR_UP & m_neighbors)
    convexLogAdd(damper, m_messages[MSG_UP], m_newMessages[MSG_UP]);

  if (NBR_DOWN & m_neighbors)
    convexLogAdd(damper, m_messages[MSG_DOWN], m_newMessages[MSG_DOWN]);

  if (NBR_LEFT & m_neighbors)
    convexLogAdd(damper, m_messages[MSG_LEFT], m_newMessages[MSG_LEFT]);

  if (NBR_RIGHT & m_neighbors)
    convexLogAdd(damper, m_messages[MSG_RIGHT], m_newMessages[MSG_RIGHT]);

  updateProduct();
}

/** Update messages on an array of nodes
 */
void NodeMessages::updateMessages(NodeMessages** nodeMsgs, int numNodes)
{

  for (int i=0 ; i<numNodes ; i++, nodeMsgs++)
    (*nodeMsgs)->updateMessages();

}

/** Convex update of messages on an array of nodes
 */
void NodeMessages::updateMessages(FloatType damper, 
                                         NodeMessages** nodeMsgs, 
                                         int numNodes)
{

  for (int i=0 ; i<numNodes ; i++, nodeMsgs++)
    (*nodeMsgs)->updateMessages(damper);

}
  
/* Erase all message (current and updated) by setting them to zeros.
 */
void NodeMessages::resetMessages()
{

  if (NBR_UP & m_neighbors)
    {
      zeroMessages(m_messages[MSG_UP]);
      zeroMessages(m_newMessages[MSG_UP]);
    }

  if (NBR_DOWN & m_neighbors)
    {
      zeroMessages(m_messages[MSG_DOWN]);
      zeroMessages(m_newMessages[MSG_DOWN]);
    }

  if (NBR_LEFT & m_neighbors) 
    {
      zeroMessages(m_messages[MSG_LEFT]);
      zeroMessages(m_newMessages[MSG_LEFT]);
    }
    
  if (NBR_RIGHT & m_neighbors)
    {
      zeroMessages(m_messages[MSG_RIGHT]);
      zeroMessages(m_newMessages[MSG_RIGHT]);
    }
}

/* Helper function for erasing message values.
 *
 * Assumes the array is of length m_nLabels.
 */
void NodeMessages::zeroMessages(FloatType* msg)
{

  for (int i=0 ; i<m_nLabels ; i++, msg++)
    *msg = 0;

}

/* Calculate the difference between the current and updated messages in 
 * a direction.
 *
 */
FloatType NodeMessages::messageDifference(const DiffType type, int direction)
{
  return messageDifference(type, m_messages[direction], 
                           m_newMessages[direction]);
}

/* Calculate the total difference between the current and updated messages 
 * for a node
 *
 */
FloatType NodeMessages::messageDifference(const DiffType type)
{
  
  bool isValid = false;
  FloatType diff = 0, diff1;

  if (NBR_UP & m_neighbors)
    {
      diff = messageDifference(type, m_messages[MSG_UP], m_newMessages[MSG_UP]);
      isValid = true;
    }

  if (NBR_DOWN & m_neighbors)
    {
      diff1 = messageDifference(type, m_messages[MSG_DOWN], 
                                m_newMessages[MSG_DOWN]);

      if (isValid)
        diff = combineDifference(type, diff,diff1);
      else
        {
          diff = diff1;
          isValid = true;
        }
    }

  if (NBR_LEFT & m_neighbors) 
    {
      diff1 = messageDifference(type, m_messages[MSG_LEFT],
                                m_newMessages[MSG_LEFT]);

      if (isValid)
        diff = combineDifference(type, diff,diff1);
      else
        {
          diff = diff1;
          isValid = true;
        }
    }
    
  if (NBR_RIGHT & m_neighbors)
    {
      diff1 = messageDifference(type, m_messages[MSG_RIGHT],
                                m_newMessages[MSG_RIGHT]);

      if (isValid)
        diff = combineDifference(type, diff,diff1);
      else
        {
          diff = diff1;
          isValid = true;
        }
    }

  if (type==L2)
    diff = sqrt(diff);

  return diff;

}

/** Actually calculates the difference between two vectors.
 * 
 * Note: Assumes the vectors have length m_nLabels
 */
FloatType NodeMessages::messageDifference(const DiffType type, 
                                          FloatType* msg, FloatType* nmsg)
{
  
  FloatType diff=0,diff1;

  switch(type)
    {
    case L1:
    case MAX:
        diff = fabs(*msg - *nmsg);
    case L2:
      diff = *msg - *nmsg;
      diff *= diff;
    case KL:
      diff = exp(*msg) * (*msg - *nmsg);
    }

  msg++;
  nmsg++;

  for (int i=1 ; i<m_nLabels ; i++, msg++, nmsg++)
    switch(type)
      {
      case L1:
        diff += fabs(*msg - *nmsg);
      case L2:
        diff1 = *msg - *nmsg;
        diff1 *= diff1;
        diff += diff1;
      case MAX:
        diff1 = fabs(*msg  - *nmsg);
        diff = MAX(diff, diff1);
      case KL:
        diff += exp(*msg) * (*msg - *nmsg);
      }

  return diff;
}

/* Helper function for merging the message differences from two directions
 */
FloatType NodeMessages::combineDifference(const DiffType type, 
                                          FloatType diff1, FloatType diff2)
{
  FloatType diff = 0;

  switch (type)
    {
    case L1:
    case L2:
    case KL:
      diff = diff1+diff2;
    case MAX:
      diff = MAX(diff1,diff2);
      
    }

  return diff;
}

/**
 * Calculate the total difference between current and updated messages on an
 * array of NodeMessages
 */
FloatType NodeMessages::messageDifference(const DiffType type,
                                          NodeMessages** nodeMsgs,
                                          int numNodes)
{
  FloatType diff = 0;

  if (numNodes<1)
    return diff;

  diff = (*nodeMsgs)->messageDifference(type);

  nodeMsgs++;

  for (int i=1 ; i<numNodes ; i++, nodeMsgs++)
    diff = combineDifference(type, diff, 
                             (*nodeMsgs)->messageDifference(type));
  
  return diff;
}
      
/** Calculate product of current messages
 *
 * Since messages are stored in log-space, this is actually the sum
 */
void NodeMessages::updateProduct()
{

  resetProduct();

  if (NBR_UP & m_neighbors)
    addTo(m_product,m_messages[MSG_UP]);

  if (NBR_DOWN & m_neighbors) 
    addTo(m_product,m_messages[MSG_DOWN]);
  
  if (NBR_LEFT & m_neighbors) 
    addTo(m_product,m_messages[MSG_LEFT]);
  
  if (NBR_RIGHT & m_neighbors) 
    addTo(m_product,m_messages[MSG_RIGHT]);
  
}


/** Reset the product (to zeros) 
 */
void NodeMessages::resetProduct()
{
  FloatType* prod = m_product;

  for (int i=0 ; i<m_nLabels ; i++, prod++)
    *prod=0;

}

/** Return the product array. 
 *
 * Dangerous. But play nice.
 */
FloatType* NodeMessages::getProduct()
{
  return m_product;
}

/** Return the message array from a particular direction.
 *
 * Dangerous. But play nice.
 */
FloatType* NodeMessages::getMessage(int direction)
{
  return m_messages[direction]; // brazen lack of error checking
}

/** Return the new message array from a particular direction.
 *
 * Dangerous. But play nice.
 */
FloatType* NodeMessages::getNewMessage(int direction)
{
  return m_newMessages[direction];
}

/** Gives the opposite of the message direction indicated
 */
int NodeMessages::reverseDir(int direction)
{

  assert(direction==MSG_UP || direction==MSG_DOWN || 
         direction==MSG_LEFT || direction==MSG_RIGHT);

  switch (direction)
    {
    case MSG_UP:
        return MSG_DOWN;
    case MSG_DOWN:
      return MSG_UP;
    case MSG_LEFT:
      return MSG_RIGHT;
    case MSG_RIGHT:
      return MSG_LEFT;
    }
}

