// (C) 2002 Marshall Tappen, MIT AI Lab

#ifndef _reg_h
#define _reg_h

#define FLOATTYPE float
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3
#include "SumProdBP.h"

class SumProdBP;

class OneNodeClusterMarg
{
public:
  OneNodeClusterMarg();
  
  //static const int UP = 0, DOWN = 1, LEFT = 2, RIGHT = 3;
  static int numStates;
  
  FLOATTYPE   *receivedMsgs[4],
              *nextRoundReceivedMsgs[4],
              *localEv;


//   FLOATTYPE *psi[4];
  void ComputeMsgRight(FLOATTYPE *msgDest, int r, int c, SumProdBP *mrf);
  void ComputeMsgUp(FLOATTYPE *msgDest, int r, int c, SumProdBP *mrf);

  void ComputeMsgLeft(FLOATTYPE *msgDest, int r, int c, SumProdBP *mrf);

  void ComputeMsgDown(FLOATTYPE *msgDest, int r, int c, SumProdBP *mrf);

  void getBelief(FLOATTYPE *beliefVec);
  int getBeliefMaxInd();
  
  int msgChange(FLOATTYPE thresh);

  void deliverMsgs();

  // CPAL - added for 
  //FLOATTYPE sumNegLogProb(FLOATTYPE a, FLOATTYPE b);

    
};

void initOneNodeMsgMem(OneNodeClusterMarg *nodeArray, FLOATTYPE *memChunk, const int numNodes, 
                       const int msgChunkSize);

void computeMessagesUpDown(OneNodeClusterMarg *nodeArray, const int numCols, const int numRows,
                           const int currCol, const FLOATTYPE alpha, SumProdBP *mrf);
void computeMessagesLeftRight(OneNodeClusterMarg *nodeArray, const int numCols, const int numRows,
                              const int currRow, const FLOATTYPE alpha, SumProdBP *mrf);

void computeOneNodeMessagesPeriodic(OneNodeClusterMarg *nodeTopArray, OneNodeClusterMarg *nodeBotArray,
                                    const int numCols, const FLOATTYPE alpha);


#endif
