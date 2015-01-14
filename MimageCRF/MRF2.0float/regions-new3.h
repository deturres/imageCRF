// (C2) 2006 Chris Pal, University of Massachusetts
// Modified from:
// (C) 2002 Marshall Tappen, MIT AI Lab
// Now used for Mean Field Message Passing

#ifndef _reg_h
#define _reg_h

#define FLOATTYPE float
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3
#include "MeanFieldMP.h"

class MeanFieldMP;
class OneNodeClusterMF
{
public:
  OneNodeClusterMF();
  
  //static const int UP = 0, DOWN = 1, LEFT = 2, RIGHT = 3;
  static int numStates;
  
  FLOATTYPE   *receivedMsgs[4],
              *nextRoundReceivedMsgs[4],
              *localEv;


  //   FLOATTYPE *psi[4];
  void ComputeMsgRight(FLOATTYPE *msgDest, int r, int c, MeanFieldMP *mrf);
  void ComputeMsgUp(FLOATTYPE *msgDest, int r, int c, MeanFieldMP *mrf);
  void ComputeMsgLeft(FLOATTYPE *msgDest, int r, int c, MeanFieldMP *mrf);
  void ComputeMsgDown(FLOATTYPE *msgDest, int r, int c, MeanFieldMP *mrf);

  // CPAL - added these 
  void ComputeMsgRightMF(FLOATTYPE *msgDest, int r, int c, MeanFieldMP *mrf);
  void ComputeMsgUpMF(FLOATTYPE *msgDest, int r, int c, MeanFieldMP *mrf);
  void ComputeMsgLeftMF(FLOATTYPE *msgDest, int r, int c, MeanFieldMP *mrf);
  void ComputeMsgDownMF(FLOATTYPE *msgDest, int r, int c, MeanFieldMP *mrf);


  void getBelief(FLOATTYPE *beliefVec);
  int getBeliefMaxInd();
  
  int msgChange(FLOATTYPE thresh);

  void deliverMsgs();
    
};

void initOneNodeMsgMem(OneNodeClusterMF *nodeArray, FLOATTYPE *memChunk, const int numNodes, 
                       const int msgChunkSize);

void computeMessagesUpDown(OneNodeClusterMF *nodeArray, const int numCols, const int numRows,
                           const int currCol, const FLOATTYPE alpha, MeanFieldMP *mrf);
void computeMessagesLeftRight(OneNodeClusterMF *nodeArray, const int numCols, const int numRows,
                              const int currRow, const FLOATTYPE alpha, MeanFieldMP *mrf);

void computeOneNodeMessagesPeriodic(OneNodeClusterMF *nodeTopArray, OneNodeClusterMF *nodeBotArray,
                                    const int numCols, const FLOATTYPE alpha);

// CPAL - added these for MF updates
void computeMessagesUpDownMF(OneNodeClusterMF *nodeArray, const int numCols, const int numRows,
                           const int currCol, const FLOATTYPE alpha, MeanFieldMP *mrf);
void computeMessagesLeftRightMF(OneNodeClusterMF *nodeArray, const int numCols, const int numRows,
                              const int currRow, const FLOATTYPE alpha, MeanFieldMP *mrf);




#endif
