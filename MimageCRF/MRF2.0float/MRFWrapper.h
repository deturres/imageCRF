/*
 * MRFWrapper.h
 *
 *  Created on: May 6, 2012
 *      Author: bhole
 */

#ifndef MRFWRAPPER_H_
#define MRFWRAPPER_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include "mrf.h"
#include "LinkedBlockList.h"
#include "regions-new.h"  // for OneNodeCluster
#include "common_typedefs.h"

#define FloatType float
#define FLOATTYPE float

class OneNodeCluster;

class MRFWrapper : public MRF {
public:
  typedef CostVal REAL;
  MRFWrapper(int width, int height, int nLabels, EnergyFunction *eng);
  MRFWrapper(int nPixels, int nLabels,EnergyFunction *eng);
  ~MRFWrapper();
  void setNeighbors(int pix1, int pix2, CostVal weight);
  int traverseNeighbors(int pixel1);
  Label getLabel(int pixel){return(m_answer[pixel]);};
  void setLabel(int pixel,Label label){m_answer[pixel] = label;};
  Label* getAnswerPtr(){return(m_answer);};
  void clearAnswer();
  void setParameters(int , void *){printf("No optional parameters to set");}
  EnergyVal smoothnessEnergy();
  EnergyVal dataEnergy();
  EnergyFunction *getEnergyFunction();
  int getWidth();
  int getHeight();
  FLOATTYPE *getScratchMatrix();
  int getNLabels();
  bool varWeights();
  void setExpScale(int expScale);
  friend void getPsiMat(OneNodeCluster &cluster, FLOATTYPE *&destMatrix,
      int r, int c, MRF *mrf, int direction, FLOATTYPE &var_weight);

  InputType getSmoothType();
  FLOATTYPE getExpV(int i);
  FLOATTYPE *getExpV();

  CostVal getHorizWeight(int r, int c);
  CostVal getVertWeight(int r, int c);


  CostVal m_lambda;
  CostVal m_smoothMax;
  int m_smoothExp;

  enum //Borrowed from BP-S.h by vnk
    {
      NONE,
      L1,
      L2,
      FIXED_MATRIX,
      GENERAL,
      BINARY,
    } m_type;

protected:
  void setData(DataCostFn dcost);
  void setData(CostVal* data);
  void setSmoothness(SmoothCostGeneralFn cost);
  void setSmoothness(CostVal* V);
  void setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda);
  void setCues(CostVal* hCue, CostVal* vCue);
  void initializeAlg();
  void BPinitializeAlg();
  virtual void optimizeAlg(int nIterations);

  void extractEdgeList(std::vector<std::vector<std::vector<int> > > *node_edge_list,
                       std::vector<std::vector<int> > *edge_list);

  int findIdxOfNeighbor(int mainpix, int pixel1);
  int findNeighborWithIdx(int mainpix, int idx);

  int *m_nNeighbors;
  OneNodeCluster *nodeArray;

  DataCostFn m_dataFn;
  SmoothCostGeneralFn m_smoothFn;

  int checkSubmodularityEdge(const int &node1, const int &node2);

  // this function assumes working on current graph
  // graph_edge_list is indexed by edge number [edge no [node1 node2]]*
  void extractEdgeType(
      const std::vector<std::vector<int> > &graph_edge_list,
      std::vector<int> *edge_list_type);

private:
  Label *m_answer;
  CostVal *m_V;
  CostVal *m_D;
  CostVal *m_horizWeights;
  CostVal *m_vertWeights;
  FLOATTYPE m_exp_scale;

  bool m_needToFreeV;
  FLOATTYPE *m_scratchMatrix;
  FLOATTYPE *m_ExpData;
  FLOATTYPE *m_message_chunk;

  typedef struct NeighborStruct {
    int     to_node;
    CostVal weight;
  } Neighbor;

  LinkedBlockList *m_neighbors;

  int m_maxnNB;
  // int   m_messageArraySizeInBytes;

  // REAL** m_messagesNG; // non-grid message structure

};



#endif /* MRFWRAPPER_H_ */
