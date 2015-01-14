#ifndef __BPS_H__
#define __BPS_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "mrf.h"
#include "LinkedBlockList.h"


class BPS : public MRF{
 public:
    typedef CostVal REAL;

    BPS(int width, int height, int nLabels, EnergyFunction *eng);
    BPS(int nPixels, int nLabels,EnergyFunction *eng);
    ~BPS();
    void setNeighbors(int /*pix1*/, int /*pix2*/, CostVal /*weight*/); //{printf("Not implemented"); exit(1);}
    int traverseNeighbors(int pixel1);
    Label getLabel(int pixel){return(m_answer[pixel]);};
    void setLabel(int pixel,Label label){m_answer[pixel] = label;};
    Label* getAnswerPtr(){return(m_answer);};
    void clearAnswer();
    void setParameters(int /*numParam*/, void * /*param*/){printf("No optional parameters to set"); exit(1);}
    EnergyVal smoothnessEnergy();
    EnergyVal dataEnergy();

 protected:
    void setData(DataCostFn dcost);
    void setData(CostVal* data);
    void setSmoothness(SmoothCostGeneralFn cost);
    void setSmoothness(CostVal* V);
    void setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda);
    void setCues(CostVal* hCue, CostVal* vCue);
    void Allocate();
    void initializeAlg();
    void optimizeAlg(int nIterations);

 private:

    enum
	{
	    NONE,
	    L1,
	    L2,
	    FIXED_MATRIX,
	    GENERAL,
	    BINARY,
	} m_type;

    CostVal m_smoothMax; // used only if
    CostVal m_lambda;    // m_type == L1 or m_type == L2

    Label *m_answer;
    CostVal *m_V; // points to array of size nLabels^2 (if type==FIXED_MATRIX) or of size nEdges*nLabels^2 (if type==GENERAL)
    CostVal *m_D;
    CostVal *m_DBinary; // valid if type == BINARY
    CostVal *m_horzWeights;
    CostVal *m_vertWeights;
    CostVal *m_horzWeightsBinary;
    CostVal *m_vertWeightsBinary;
    DataCostFn m_dataFn;
    SmoothCostGeneralFn m_smoothFn;
    bool m_needToFreeV;
    bool m_needToFreeD;

    REAL* m_messages; // size of one message: N = 1 if m_type == BINARY, N = K otherwise
    // message between edges (x,y)-(x+1,y): m_messages+(2*x+2*y*m_width)*N
    // message between edges (x,y)-(x,y+1): m_messages+(2*x+2*y*m_width+1)*N

    int	  m_messageArraySizeInBytes;

	REAL** m_messagesNG; // non-grid message structure

	//chet
	typedef struct NeighborStruct{
	 int  to_node;
	 CostVal weight;
	} Neighbor;
	LinkedBlockList *m_neighbors;

	int *m_nNeighbors;
	int m_maxnNB;

    void optimize_GRID_L1(int nIterations);
    void optimize_GRID_L2(int nIterations);
    void optimize_GRID_FIXED_MATRIX(int nIterations);
    void optimize_GRID_GENERAL(int nIterations);
    void optimize_GRID_BINARY(int nIterations);

//chet

    void optimize_NONGRID_FIXED_MATRIX_GENERAL(int nIterations);

    void printMat(REAL* M, int K);
    int findIdxOfNeighbor(int mainpix, int pixel1);
    int findNeighborWithIdx(int mainpix, int idx);
    BPS::REAL UpdateMessageFIXED_MATRIX_NG(BPS::REAL* M, BPS::REAL* Di_hat, int K, BPS::REAL gamma, MRF::CostVal lambda, MRF::CostVal* V, void* buf, BPS::REAL* N);
    BPS::REAL UpdateMessageGENERAL_NG(BPS::REAL* M, BPS::REAL* Di_hat, int K, BPS::REAL gamma, MRF::CostVal lambda, int m, int q, void* buf, BPS::REAL* N);

    void computeSolution_GRID_L1_L2_FIXED_MATRIX();
    void computeSolution_GRID_GENERAL();
	void computeSolution_GRID_BINARY();

	void computeSolution_NONGRID_FIXED_MATRIX_GENERAL();
	void fillAnswersK();

};

#endif /*  __BPS_H__ */
