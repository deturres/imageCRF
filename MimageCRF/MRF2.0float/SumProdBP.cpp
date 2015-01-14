#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "SumProdBP.h"
#include "regions-new2.h"

#define m_D(pix,l)  m_D[(pix)*m_nLabels+(l)]
#define m_V(l1,l2)  m_V[(l1)*m_nLabels+(l2)]


		

SumProdBP::SumProdBP(int width, int height, int nLabels,EnergyFunction *eng):MRF(width,height,nLabels,eng)
{
	m_needToFreeV = 0;
	BPinitializeAlg();
}
SumProdBP::SumProdBP(int nPixels, int nLabels,EnergyFunction *eng):MRF(nPixels,nLabels,eng)
{
	m_needToFreeV = 0;
	BPinitializeAlg();
}

SumProdBP::~SumProdBP()
{ 
	delete[] m_answer;
	if (m_message_chunk) delete[] m_message_chunk;
	if (!m_grid_graph) delete[] m_neighbors;
	if ( m_needToFreeV ) delete[] m_V;
}


void SumProdBP::initializeAlg()
{
}

void SumProdBP::BPinitializeAlg()
{
	m_answer = (Label *) new Label[m_nPixels];
	if ( !m_answer ){printf("\nNot enough memory, exiting");exit(0);}

	m_scratchMatrix = new FLOATTYPE[m_nLabels * m_nLabels];
	// MEMORY LEAK? where does this ever get deleted??

	nodeArray =     new OneNodeClusterMarg[m_nPixels];
	// MEMORY LEAK? where does this ever get deleted??

	OneNodeClusterMarg::numStates = m_nLabels;

	if (!m_grid_graph)
	{
	  assert(0);
	  // Only Grid Graphs are supported
		m_neighbors = (LinkedBlockList *) new LinkedBlockList[m_nPixels];
		if (!m_neighbors) {printf("Not enough memory,exiting");exit(0);};
	}
	else
	{
// 	  const int clen = 4*m_nPixels * m_nLabels + 4*m_nPixels * m_nLabels * m_nLabels;
	  const int clen = 4*m_nPixels * m_nLabels ;
	  //printf("clen:%d\n",clen/1024/1024);
	  m_message_chunk = (FloatType *) new FloatType[clen];
	  if ( !m_message_chunk ){printf("\nNot enough memory for messages, exiting");exit(0);}
	  for(int i = 0; i < clen; i++)
	    m_message_chunk[i]=0;
	  initOneNodeMsgMem(nodeArray, m_message_chunk, m_nPixels, m_nLabels);
	}
}

MRF::InputType SumProdBP::getSmoothType()
{
  return m_smoothType;
}

EnergyFunction *SumProdBP::getEnergyFunction()
{
  return m_e;
}

void SumProdBP::setExpScale(int expScale)
{
  m_exp_scale = (float)expScale;
}

int SumProdBP::getNLabels()
{
  return m_nLabels;
}

int SumProdBP::getWidth()
{
  return m_width;
}

int SumProdBP::getHeight()
{
  return m_height;
}
FLOATTYPE SumProdBP::getExpV(int i)
{
  return m_ExpData[i];
}

FLOATTYPE *SumProdBP::getExpV()
{
  return m_ExpData;
}


MRF::CostVal SumProdBP::getHorizWeight(int r, int c)
{
  int x = c;
  int y = r;
  int pix    = x+y*m_width;
  return m_varWeights ? m_horizWeights[pix] :  1;
}

MRF::CostVal SumProdBP::getVertWeight(int r, int c)
{
  int x = c;
  int y = r;
  int pix    = x+y*m_width;
  return  m_varWeights ? m_vertWeights[pix] :  1;
}

bool SumProdBP::varWeights()
{
  return m_varWeights;
}

FLOATTYPE *SumProdBP::getScratchMatrix()
{
  return m_scratchMatrix;
}

void SumProdBP::clearAnswer()
{
	memset(m_answer, 0, m_nPixels*sizeof(Label));
}


void SumProdBP::setNeighbors(int pixel1, int pixel2, CostVal weight)
{
  assert(0);
  //Only Grid Graphs are supported
	assert(!m_grid_graph);
	assert(pixel1 < m_nPixels && pixel1 >= 0 && pixel2 < m_nPixels && pixel2 >= 0);


	Neighbor *temp1 = (Neighbor *) new Neighbor;
	Neighbor *temp2 = (Neighbor *) new Neighbor;

	if ( !temp1 || ! temp2 ) {printf("\nNot enough memory, exiting");exit(0);}

	temp1->weight  = weight;
	temp1->to_node = pixel2;

	temp2->weight  = weight;
	temp2->to_node = pixel1;

	m_neighbors[pixel1].addFront(temp1);
	m_neighbors[pixel2].addFront(temp2);
}


MRF::EnergyVal SumProdBP::smoothnessEnergy()
{
    double eng =  0;
    double weight;
    int x,y,pix;
        if ( m_smoothType != FUNCTION  )
        {
            for ( y = 0; y < m_height; y++ )
                for ( x = 1; x < m_width; x++ )
                {
                    pix    = x+y*m_width;
                    weight = m_varWeights ? m_horizWeights[pix-1] :  1;
                    eng = eng + m_V(m_answer[pix],m_answer[pix-1])*weight;
                }

            for ( y = 1; y < m_height; y++ )
                for ( x = 0; x < m_width; x++ )
                {
                    pix = x+y*m_width;
                    weight = m_varWeights ? m_vertWeights[pix-m_width] :  1;
                    eng = eng + m_V(m_answer[pix],m_answer[pix-m_width])*weight;
                }
        }
        else
        {
            for ( y = 0; y < m_height; y++ )
                for ( x = 1; x < m_width; x++ )
                {
                    pix = x+y*m_width;
                    eng = eng + m_smoothFn(pix,pix-1,m_answer[pix],m_answer[pix-1]);
                }

            for ( y = 1; y < m_height; y++ )
                for ( x = 0; x < m_width; x++ )
                {
                    pix = x+y*m_width;
                    eng = eng + m_smoothFn(pix,pix-m_width,m_answer[pix],m_answer[pix-m_width]);
                }
        }
	return((EnergyVal)eng);

}



MRF::EnergyVal SumProdBP::dataEnergy()
{
    double eng =  0;

    
    if ( m_dataType == ARRAY) 
    {
        for ( int i = 0; i < m_nPixels; i++ )
            eng = eng + m_D(i,m_answer[i]);
    }
    else
    {
        for ( int i = 0; i < m_nPixels; i++ )
            eng = eng + m_dataFn(i,m_answer[i]);
    }
    return((EnergyVal) eng);

}


void SumProdBP::setData(DataCostFn dcost)
{
  m_dataFn = dcost;
  int i;
  int j;
  m_ExpData = new FloatType[m_nPixels * m_nLabels];
  // MEMORY LEAK? where does this ever get deleted??
  if(!m_ExpData)
  {
    exit(0);
  }
  
  
  m_exp_scale = 1;//FLOATTYPE(cmax)*4.0;
  FloatType *cData = m_ExpData;
  for ( i= 0; i < m_nPixels; i++)
  {
    nodeArray[i].localEv = cData;
    for( j = 0; j < m_nLabels; j++)
    {
      *cData = (float)m_dataFn(i,j);
      cData++;
    }
  }
}

void SumProdBP::setData(CostVal* data)
{
  int i;
  int j;
  m_D = data;
  m_ExpData = new FloatType[m_nPixels * m_nLabels];
  // MEMORY LEAK? where does this ever get deleted??
  if(!m_ExpData)
  {
    exit(0);
  }
  
  
  m_exp_scale = 1;//FLOATTYPE(cmax)*4.0;
  FloatType *cData = m_ExpData;
  for ( i= 0; i < m_nPixels; i++)
  {
    nodeArray[i].localEv = cData;
    for( j = 0; j < m_nLabels; j++)
    {
      *cData = (float)m_D(i,j);
      cData++;
    }
  }
}


void SumProdBP::setSmoothness(SmoothCostGeneralFn cost)
{
  m_smoothFn = cost;
}
void SumProdBP::setSmoothness(CostVal* V)
{
  m_type = FIXED_MATRIX;
  m_V = V;
}


void SumProdBP::setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda)
{
	int i, j;
	CostVal cost;
	m_type = (smoothExp == 1) ? L1 : L2; //borrowed from BP-S.cpp from vnk
	m_lambda = lambda;
	m_smoothMax = smoothMax;
	m_smoothExp = smoothExp;
	m_needToFreeV = 1;

	m_V = (CostVal *) new CostVal[m_nLabels*m_nLabels*sizeof(CostVal)];
	if (!m_V) { fprintf(stderr, "Not enough memory!\n"); exit(1); }


	for (i=0; i<m_nLabels; i++)
		for (j=i; j<m_nLabels; j++)
		{
			cost = (smoothExp == 1) ? j - i : (j - i)*(j - i);
   			if (cost > smoothMax) cost = smoothMax;
			m_V[i*m_nLabels + j] = m_V[j*m_nLabels + i] = cost*lambda;
		}

	
	
}


void SumProdBP::setCues(CostVal* hCue, CostVal* vCue)
{
	m_horizWeights = hCue;
	m_vertWeights  = vCue;
}


void SumProdBP::optimizeAlg(int nIterations)
{
    //int x, y, i, j, n;
    //Label* l;
    //CostVal* dataPtr;

	if ( !m_grid_graph) {printf("\nSumProdBP is not implemented for nongrids yet!");exit(1);}

	int numRows = getHeight();
	int numCols = getWidth();
	const FLOATTYPE alpha = 0.8f;
	for (int niter=0; niter < nIterations; niter++)
	{
	  for(int r = 0; r < numRows; r++)
	  {
	    computeMessagesLeftRight(nodeArray, numCols, numRows, r, alpha, this);
	  }
	  for(int c = 0; c < numCols; c++)
	  {
	    computeMessagesUpDown(nodeArray, numCols, numRows, c, alpha, this);
	  }
	}

      
      Label *currAssign = m_answer;
      for(int m = 0; m < numRows; m++)
      {
	for(int n = 0; n < numCols; n++)
	{
	  int maxInd = nodeArray[m*numCols+n].getBeliefMaxInd();
	  currAssign[m * numCols +n] = maxInd;
	}
      }

}

