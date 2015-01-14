/*
 * MaxProdBPTree.h
 *
 *  Created on: Oct 3, 2011
 *      Author: bhole
 */

#ifndef MAXPRODBPTREE_H_
#define MAXPRODBPTREE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include "mrf.h"
#include "LinkedBlockList.h"
#include "regions-new.h"

#define FloatType float
#define FLOATTYPE float
class MaxProdBPTree;

class MaxProdBPTree : public MRF{
  public:
   typedef CostVal REAL;
   MaxProdBPTree(int nPixels, int nLabels,EnergyFunction *eng);
   ~MaxProdBPTree();
   void setNeighbors(int pix1, int pix2, CostVal weight);
   int traverseNeighbors(int pixel1);
   Label getLabel(int pixel){return(m_answer[pixel]);};
   void setLabel(int pixel,Label label){m_answer[pixel] = label;};
   Label* getAnswerPtr(){return(m_answer);};
   void clearAnswer();
   void setParameters(int , void *) {
     printf("No optional parameters to set");
   }

   virtual MRF::EnergyVal call_smoothFn(int pix_i, int pix_j, int label_i, int label_j);

   EnergyVal smoothnessEnergy();
   EnergyVal dataEnergy();
   EnergyFunction *getEnergyFunction();
   FLOATTYPE *getScratchMatrix();
   int getNLabels();
   bool varWeights();
   void setExpScale(int expScale);
   friend void getPsiMat(OneNodeCluster &cluster, FLOATTYPE *&destMatrix,
			int r, int c, MaxProdBP *mrf, int direction, FLOATTYPE &var_weight);

  InputType getSmoothType();
  FLOATTYPE getExpV(int i);
  FLOATTYPE *getExpV();

  void updateData(CostVal* data);

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
    SmoothCostGeneralFn getSmoothnessFn();
    void setSmoothness(CostVal* V);
    void setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda);
    void setCues(CostVal* hCue, CostVal* vCue){};
    void initializeAlg();
    void BPinitializeAlg();
    void optimizeAlg(int nIterations);

  private:
	void extractParentArrayUsingBFS(const int& root);
    Label *m_answer;
    CostVal *m_V;
    CostVal *m_D;
    CostVal *m_horizWeights;
    CostVal *m_vertWeights;
    FLOATTYPE m_exp_scale;
    DataCostFn m_dataFn;
    SmoothCostGeneralFn m_smoothFn;
    bool m_needToFreeV;
    FLOATTYPE *m_scratchMatrix;
    FLOATTYPE *m_ExpData;
    FLOATTYPE *m_message_chunk;
    OneNodeCluster *nodeArray;

    typedef struct NeighborStruct {
      int     to_node;
      CostVal weight;
    } Neighbor;

    LinkedBlockList *m_neighbors;

    int *m_nNeighbors;
    int m_maxnNB;
    int	m_messageArraySizeInBytes;

    REAL** m_messagesNG; // non-grid message structure

    int findIdxOfNeighbor(int mainpix, int pixel1);
    int findNeighborWithIdx(int mainpix, int idx);

    std::vector<std::vector<int> > parent_array_;
};



#endif /* MAXPRODBPTREE_H_ */
