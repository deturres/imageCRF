/** Class for storing the (up-to) 4 messages coming into a node on a grid.
 *
 * Two versions of messages are stored: a "current" version from which
 * information may be drawn, and an "updated" or "new" version which
 * are used to build up updates from the current messages. This is
 * important for synchronous message updating schemes.
 *
 * In addition, the product of messages coming into a node is
 * calculated and cached.
 *
 * Jerod Weinman
 * jerod@acm.org
 */

#ifndef __NODEMESSAGES_H__
#define __NODEMESSAGES_H__

#include "ArrayMath.h"

#define FloatType float

class NodeMessages {

 public:
  
  // Indices for the messages from the stated direction
  static int const MSG_UP;
  static int const MSG_DOWN;
  static int const MSG_LEFT;
  static int const MSG_RIGHT;

  // Indicator bit for a neighbor in the stated direction
  static int const NBR_UP;
  static int const NBR_DOWN;
  static int const NBR_LEFT;
  static int const NBR_RIGHT;


  // Optional argument neighbors (see member definition below for details)
  NodeMessages(int nLabels, int neighbors = NBR_UP|NBR_DOWN|NBR_LEFT|NBR_RIGHT );
  //NodeMessages(int nLabels, int neighbors);
 
  ~NodeMessages();

 // Type for measuring the difference between current and updated messages
  typedef enum
    {
      L1, 
      L2,
      MAX,
      KL
    } DiffType;
  
  
  // "Copy" new messages to current messages
  void updateMessages();
  void updateMessage(int direction);

  // Take convex combination of new messages with current messages
  void updateMessages(FloatType damper);
  void updateMessage(int direction, FloatType damper);

  static void updateMessages(NodeMessages**, int numNodes);
  static void updateMessages(FloatType damper, NodeMessages**, int numNodes);

  // Calculate product of current messages
  void updateProduct();

  // Erases all messages (current and updated)
  void resetMessages();

  // Erase the product
  void resetProduct();

  // Get the product array. Play nice with it! It's not a copy!!
  FloatType* getProduct();

  // Get a message array.  Play nice with it! It's not a copy!!
  FloatType* getMessage(int direction);

  // Calculate the difference between the current and updated messages 
  // for the purposes of checking convergence 
  FloatType messageDifference(const DiffType type);
  FloatType messageDifference(const DiffType type, int direction);

  // Calculate the total difference between current and updated messages
  static FloatType messageDifference(const DiffType, NodeMessages**, 
                                     int numNodes);

  // Access for writing new messages
  FloatType* getNewMessage(int direction);

  // Opposite of the message direction indicated
  static int reverseDir(int direction);

  // Helper function for merging the message differences from two directions
  static FloatType combineDifference(const DiffType, FloatType, FloatType);

 protected:


  // Number of labels/states for the node
  int m_nLabels;

  // Neighbors of a node; indicated via bitwise OR with the constants 
  // NBR_UP, NBR_DOWN, NBR_LEFT, and NBR_RIGHT.
  int m_neighbors;

  // Current incoming messages
  FloatType* m_messages[4]; 

  // Storage for updated messages
  FloatType* m_newMessages[4];

  // Cached "product" of incoming messages
  // (Since everything is stored in log-space, this is actually a sum)
  FloatType* m_product;
  
 private:

  // Helper function for resetting message values
  void zeroMessages(FloatType*);
    
  // Calculates actual message differences
  FloatType messageDifference(const DiffType, FloatType*, FloatType*);

  // Adds src to accumulator dst
  inline void addTo(FloatType* dst, FloatType* src)
    { ArrayMath::addTo(dst,src,m_nLabels); };

  // Subtracts vector sth from dst
  inline void subtractFrom(FloatType* dst, FloatType* sth)
    { ArrayMath::subtractFrom(dst,sth,m_nLabels); };

  
  // Convex update of array: dst = log((1-alpha)*exp(dst) + alpha*exp(src))
  inline void convexLogAdd(const FloatType alpha, FloatType* dst, 
                           FloatType* src) 
    { ArrayMath::convexLogAdd(alpha,dst,src,m_nLabels); };


};

#endif
