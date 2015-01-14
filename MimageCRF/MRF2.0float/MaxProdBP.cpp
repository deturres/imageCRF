#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "MaxProdBP.h"
#include "regions-new.h"

#define m_D(pix,l)  m_D[(pix)*m_nLabels+(l)]
#define m_V(l1,l2)  m_V[(l1)*m_nLabels+(l2)]

#define MIN(a,b)  (((a) < (b)) ? (a) : (b))
#define MAX(a,b)  (((a) > (b)) ? (a) : (b))
#define TRUNCATE_MIN(a,b) { if ((a) > (b)) (a) = (b); }
#define TRUNCATE_MAX(a,b) { if ((a) < (b)) (a) = (b); }
#define TRUNCATE TRUNCATE_MIN


		

MaxProdBP::MaxProdBP(int width, int height, int nLabels,EnergyFunction *eng):MRF(width,height,nLabels,eng)
{
	m_needToFreeV = 0;
	BPinitializeAlg();
}
MaxProdBP::MaxProdBP(int nPixels, int nLabels,EnergyFunction *eng):MRF(nPixels,nLabels,eng)
{
	m_needToFreeV = 0;
	BPinitializeAlg();
}

MaxProdBP::~MaxProdBP()
{ 
	delete[] m_answer;
	//if (m_message_chunk) delete[] m_message_chunk;
	//if (!m_grid_graph) delete[] m_neighbors;
	if ( m_needToFreeV ) delete[] m_V;

	//chet
	// not sure if the following two should be here. Need to check
	delete[] m_scratchMatrix;
	delete[] nodeArray;

	//chet new addition
	delete[] m_ExpData;


	//chet - what about these??

/*
    if ( m_needToFreeD ) delete [] m_D;
    if ( m_needToFreeV ) delete [] m_V;
    if ( m_messages ) delete [] m_messages;
    if ( m_DBinary ) delete [] m_DBinary;
    if ( m_horzWeightsBinary ) delete [] m_horizWeightsBinary;
    if ( m_vertWeightsBinary ) delete [] m_vertWeightsBinary;

*/

	//chet
	if (m_grid_graph)
	{
		delete[] m_message_chunk;
	}
    if (!m_grid_graph) {
        // need to delete m_messagesNG
        // need to delete m_neighbors
    	for (int i=0; i<m_nPixels; i++) {
    		delete[] m_messagesNG[i];
    	}
    	delete[] m_messagesNG;

    	delete[] m_neighbors;
    	delete[] m_nNeighbors;
    }

}

int MaxProdBP::traverseNeighbors(int pixel1)
{

    assert(pixel1 < m_nPixels && pixel1 >= 0);
    assert(m_grid_graph == 0);

	LinkedBlockList *m_NBTR = &m_neighbors[pixel1];

	int cnter = 0;
	m_NBTR->setCursorFront();

	while(m_NBTR->hasNext()) {
		cnter++;
		m_NBTR->next();
	}

	return cnter;

//	printf("Pixel %d has neighbor count =  %d \n", pixel1, cnter);

}

void MaxProdBP::initializeAlg()
{

	if (!m_grid_graph)
	{
		//general graphs

		// allocate messages
		// this is message number per neighbor
		int messageNumPerNB = m_nPixels*m_nLabels;
		// true number of messages will be this multiplied with the number of neighbors per pixel-label
		// but in general graph this number of neighbors could be variable

		m_messagesNG = new REAL*[m_nPixels];
		m_nNeighbors = new int[m_nPixels];

		// need to traverse the neighbor list to get neighbors per pixel

		int totNB=0;

		m_maxnNB=0;

		// traverse Neighbors list
		for (int m=0; m < m_nPixels; m++) {
			int numNB=traverseNeighbors((int)m);

			m_nNeighbors[m] = numNB;

			if (m_maxnNB<numNB) {
				m_maxnNB=numNB;
			}

			m_messagesNG[m] = new REAL[numNB*m_nLabels];
			memset(m_messagesNG[m], 0, numNB*m_nLabels*sizeof(REAL));

			totNB +=numNB;
		}

		m_messageArraySizeInBytes = totNB*sizeof(REAL);
	}

}

void MaxProdBP::BPinitializeAlg()
{
	m_answer = (Label *) new Label[m_nPixels];
	if ( !m_answer ){printf("\nNot enough memory, exiting");exit(0);}

	m_scratchMatrix = new FLOATTYPE[m_nLabels * m_nLabels];
	// MEMORY LEAK? where does this ever get deleted??

	nodeArray =     new OneNodeCluster[m_nPixels];
	// MEMORY LEAK? where does this ever get deleted??

	OneNodeCluster::numStates = m_nLabels;

	if (!m_grid_graph)
	{
	  // assert(0);  //chetan
	  // Only Grid Graphs are supported
		m_neighbors = (LinkedBlockList *) new LinkedBlockList[m_nPixels];
		if (!m_neighbors) {printf("Not enough memory,exiting");exit(0);};
		m_varWeights = false;
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
	 // delete[] m_answer;  //chetan
	//delete[] m_scratchMatrix;
	//delete[] nodeArray;

}

MRF::InputType MaxProdBP::getSmoothType()
{
  return m_smoothType;
}

EnergyFunction *MaxProdBP::getEnergyFunction()
{
  return m_e;
}

void MaxProdBP::setExpScale(int expScale)
{
  m_exp_scale = (float)expScale;
}

int MaxProdBP::getNLabels()
{
  return m_nLabels;
}

int MaxProdBP::getWidth()
{
  return m_width;
}

int MaxProdBP::getHeight()
{
  return m_height;
}
FLOATTYPE MaxProdBP::getExpV(int i)
{
  return m_ExpData[i];
}

FLOATTYPE *MaxProdBP::getExpV()
{
  return m_ExpData;
}


MRF::CostVal MaxProdBP::getHorizWeight(int r, int c)
{
  int x = c;
  int y = r;
  int pix    = x+y*m_width;
  return m_varWeights ? m_horizWeights[pix] :  1;
}

MRF::CostVal MaxProdBP::getVertWeight(int r, int c)
{
  int x = c;
  int y = r;
  int pix    = x+y*m_width;
  return  m_varWeights ? m_vertWeights[pix] :  1;
}

bool MaxProdBP::varWeights()
{
  return m_varWeights;
}

FLOATTYPE *MaxProdBP::getScratchMatrix()
{
  return m_scratchMatrix;
}

void MaxProdBP::clearAnswer()
{
	memset(m_answer, 0, m_nPixels*sizeof(Label));
}


void MaxProdBP::setNeighbors(int pixel1, int pixel2, CostVal weight)
{
  //assert(0);  // chetan
  //Only Grid Graphs are supported
	// assert(!m_grid_graph);  //chetan

	assert(m_grid_graph == 0);
	assert(pixel1 < m_nPixels && pixel1 >= 0 && pixel2 < m_nPixels && pixel2 >= 0);


	Neighbor *temp1 = (Neighbor *) new Neighbor;
	Neighbor *temp2 = (Neighbor *) new Neighbor;

	if ( !temp1 || ! temp2 ) {printf("\nNot enough memory, exiting");exit(0);}

	temp1->weight  = weight;
	temp1->to_node = pixel2;

	temp2->weight  = weight;
	temp2->to_node = pixel1;

	m_neighbors[pixel1].addFront(temp1, 1);
	m_neighbors[pixel2].addFront(temp2, 1);
}


MRF::EnergyVal MaxProdBP::smoothnessEnergy()
{
    double eng =  0;
    double weight;
    int x,y,pix;

    if ( m_grid_graph )
	{
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
	}
    else // general graph
	{
	    if ( m_smoothType != FUNCTION  )
		{
		    for (int i=0; i<m_nPixels; i++)
			{
		    	LinkedBlockList *m_NBTR = &m_neighbors[i];
		    	m_NBTR->setCursorFront();
		    	while(m_NBTR->hasNext()) {

		    		Neighbor *temp1 = (Neighbor *) m_NBTR->next();
		    	    int pix2  = temp1->to_node;

		    		eng = eng + m_V(m_answer[i],m_answer[pix2]);
		    	}
			}
		    eng = eng/2;
		}
	    else // smooth is a function
		{
		    for (int i=0; i<m_nPixels; i++)
			{
		    	LinkedBlockList *m_NBTR = &m_neighbors[i];
		    	m_NBTR->setCursorFront();
		    	while(m_NBTR->hasNext()) {
		    		Neighbor *temp1 = (Neighbor *) m_NBTR->next();
		    	    int pix2  = temp1->to_node;
		    		eng = eng + m_smoothFn(i, pix2, m_answer[i], m_answer[pix2]);

		    	}
			}
		    eng = eng/2;
		}
	}

    return((EnergyVal) eng);

}



MRF::EnergyVal MaxProdBP::dataEnergy()
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


void MaxProdBP::setData(DataCostFn dcost)
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

void MaxProdBP::setData(CostVal* data)
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


void MaxProdBP::setSmoothness(SmoothCostGeneralFn cost)
{
  m_smoothFn = cost;
}
void MaxProdBP::setSmoothness(CostVal* V)
{
  m_type = FIXED_MATRIX;
  m_V = V;
}


void MaxProdBP::setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda)
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


void MaxProdBP::setCues(CostVal* hCue, CostVal* vCue)
{
	m_horizWeights = hCue;
	m_vertWeights  = vCue;
}


inline void SetVectorToZero(MaxProdBP::REAL* to, int K)
{
	MaxProdBP::REAL* to_finish = to + K;
    do
	{
	    *to ++ = 0;
	} while (to < to_finish);
}

inline void AddAllMessages(MaxProdBP::REAL* to, MaxProdBP::REAL* from, int K)
{
	MaxProdBP::REAL* to_finish = to + K;
    do
	{
	    *to ++ += *from ++;
	} while (to < to_finish);
}

inline void AddDataCost(MaxProdBP::REAL* to, MaxProdBP::REAL* from, int K)
{
	MaxProdBP::REAL* to_finish = to + K;
    do
	{
	    *to ++ -= *from ++;  // the negative is because we want to have -ve datacost added
	} while (to < to_finish);
}


int MaxProdBP::findIdxOfNeighbor(int mainpix, int pixel1)
{

    assert(pixel1 < m_nPixels && pixel1 >= 0);
    assert(mainpix < m_nPixels && mainpix >= 0);
    assert(m_grid_graph == 0);

	LinkedBlockList *m_NBTR = &m_neighbors[mainpix];

	m_NBTR->setCursorFront();

	int idx=-1;

	int cnter=0;
	while(m_NBTR->hasNext()) {
		Neighbor* nb = (Neighbor *)m_NBTR->next();
		if (nb->to_node == pixel1) {
			idx = cnter;
			break;
		}
		cnter++;
	}

	return idx;
}

int MaxProdBP::findNeighborWithIdx(int mainpix, int idx)
{

    assert(mainpix < m_nPixels && mainpix >= 0);
    assert(m_grid_graph == 0);

	LinkedBlockList *m_NBTR = &m_neighbors[mainpix];

	int pixel=-1;

	int cnter = 0;
	m_NBTR->setCursorFront();

	while(m_NBTR->hasNext()) {
		if (idx==cnter) {
			Neighbor* nb = (Neighbor *) m_NBTR->next();
			pixel = nb->to_node;
			break;
		} else {
			cnter++;
			m_NBTR->next();
		}
	}

	return pixel;
}




void MaxProdBP::optimizeAlg(int nIterations)
{
    //int x, y, i, j, n;
    //Label* l;
    //CostVal* dataPtr;

	// if ( !m_grid_graph) {printf("\nMaxProdBP is not implemented for nongrids yet!");exit(1);}

	if (m_grid_graph)
	{

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

	} else {
		// for non-grid graphs

		const FLOATTYPE alpha = 0.8f;

		REAL* Di;
		Di = new REAL[m_nLabels];

		REAL* Di2;
		Di2 = new REAL[m_nLabels];

		REAL* M_ptr, *N_ptr;
		REAL delta;

		CostVal lambda = 1;

		for (int niter=0; niter < nIterations; niter++)
		{

			// forward pass

			for (int m=0; m<m_nPixels; m++) {

				SetVectorToZero(Di, m_nLabels);
		    	M_ptr = m_messagesNG[m];
		    	for (int p=0; p<m_nNeighbors[m]; p++, M_ptr+=m_nLabels) {
		    		AddAllMessages(Di, M_ptr, m_nLabels);
		    	}
		    	AddDataCost(Di, nodeArray[m].localEv, m_nLabels);

          // for (int kix=0; kix<m_nLabels; kix++)
          //  printf(" %f ", Di[kix]);

		    	M_ptr = m_messagesNG[m];
		    	for (int p=0; p<m_nNeighbors[m]/2; p++, M_ptr+=m_nLabels) {

		    		// m <-- q is the message pointed by M_ptr

		    		int q=findNeighborWithIdx(m, p);
		    		if (q==-1) {
		    			printf("Error\n");
		    			exit(1);
		    		}
		    		int idx = findIdxOfNeighbor(q, m);
		    		if (idx==-1) {
		    			printf("Error\n");
		    			exit(1);
		    		}

	          printf(" : %d %d : ", m, q);

		    		N_ptr = m_messagesNG[q];
		    		N_ptr += idx*m_nLabels;

		    		// m --> q is the message pointed by N_ptr and we want to update this

		    		if ( m_smoothType != FUNCTION  ) { //assume ARRAY
		    			// not coded yet
		    			exit(0);
		    		} else { // this is smoothfunction is function()
		    			// UpdateMessageGENERAL_NG(M_ptr, Di, K, 1, lambda, m, q, buf, N_ptr);

		    		    for (int ki=0; ki<m_nLabels; ki++)
		    			{
		    			    Di2[ki] = Di[ki] - M_ptr[ki];
		    			}

		    		    for (int kj=0; kj<m_nLabels; kj++)
		    			{
		    			    N_ptr[kj] = Di2[0] - lambda*m_smoothFn(m, q, 0, kj);//[0];

		    			    for (int ki=1; ki<m_nLabels; ki++)
		    				{
		    			    	TRUNCATE_MAX(N_ptr[kj], Di2[ki] - lambda*m_smoothFn(m, q, ki, kj)); //V[0]);
		    				}
		    			}

              // for (int kix=0; kix<m_nLabels; kix++)
              //  printf(" %f ", N_ptr[kix]);

		    		  delta = N_ptr[0];
		    		  for (int kj=1; kj<m_nLabels; kj++) TRUNCATE(delta, N_ptr[kj]);
		    		  for (int kj=0; kj<m_nLabels; kj++) N_ptr[kj] -= delta;

	            // for (int kix=0; kix<m_nLabels; kix++)
	            //  printf(" %f ", N_ptr[kix]);
	            // printf("\n");

		    		}

		    	}


			}


			// backward pass


			for (int m=0; m<m_nPixels; m++) {

				SetVectorToZero(Di, m_nLabels);
		    	M_ptr = m_messagesNG[m];
		    	for (int p=0; p<m_nNeighbors[m]; p++, M_ptr+=m_nLabels) {
		    		AddAllMessages(Di, M_ptr, m_nLabels);
		    	}
		    	AddDataCost(Di, nodeArray[m].localEv, m_nLabels);

		    	M_ptr = m_messagesNG[m];

		    	for (int p=m_nNeighbors[m]/2; p<m_nNeighbors[m]; p++, M_ptr+=m_nLabels) {

		    		// m <-- q is the message pointed by M_ptr

		    		int q=findNeighborWithIdx(m, p);
		    		if (q==-1) {
		    			printf("Error\n");
		    			exit(1);
		    		}
		    		int idx = findIdxOfNeighbor(q, m);
		    		if (idx==-1) {
		    			printf("Error\n");
		    			exit(1);
		    		}

		    		// printf(" : %d %d : ", m, q);

		    		N_ptr = m_messagesNG[q];
		    		N_ptr += idx*m_nLabels;

		    		// m --> q is the message pointed by N_ptr and we want to update this

		    		if ( m_smoothType != FUNCTION  ) { //assume ARRAY
		    			// not coded yet
		    			exit(0);
		    		} else { // this is smoothfunction is function()
		    			// UpdateMessageGENERAL_NG(M_ptr, Di, K, 1, lambda, m, q, buf, N_ptr);

		    		    for (int ki=0; ki<m_nLabels; ki++)
		    			{
		    			    Di2[ki] = Di[ki] - M_ptr[ki];
		    			}

		    		    for (int kj=0; kj<m_nLabels; kj++)
		    			{
		    			    N_ptr[kj] = Di2[0] - lambda*m_smoothFn(m, q, 0, kj);//[0];

		    			    for (int ki=1; ki<m_nLabels; ki++)
		    				{
		    			    	TRUNCATE_MAX(N_ptr[kj], Di2[ki] - lambda*m_smoothFn(m, q, ki, kj)); //V[0]);
		    				}
		    			}

		    		    delta = N_ptr[0];
		    		    for (int kj=1; kj<m_nLabels; kj++) TRUNCATE(delta, N_ptr[kj]);
		    		    for (int kj=0; kj<m_nLabels; kj++) N_ptr[kj] -= delta;

	              //for (int kix=0; kix<m_nLabels; kix++)
	              //  printf(" %f ", N_ptr[kix]);
	              //printf("\n");


		    		}

		    	}


			}

		}


		for (int m=0; m<m_nPixels; m++)
		{

			SetVectorToZero(Di, m_nLabels);
	    	M_ptr = m_messagesNG[m];
	    	for (int p=0; p<m_nNeighbors[m]; p++, M_ptr+=m_nLabels) {
	    		AddAllMessages(Di, M_ptr, m_nLabels);
	    	}
	    	AddDataCost(Di, nodeArray[m].localEv, m_nLabels);

	    	REAL currBelief,bestBelief;
	    	int bestInd = 0;
	    	bestBelief = Di[0];

			for(int i = 1; i < m_nLabels; i++)
			{
				currBelief = Di[i];
        if(currBelief == bestBelief && rand() % 2 == 1)
        {
          bestInd=i;
          bestBelief = currBelief;
        }

				if(currBelief > bestBelief)
				{
				  bestInd=i;
				  bestBelief = currBelief;
				}
			}

			m_answer[m] = bestInd;

		}


		delete[] Di;
		delete[] Di2;

	}

}

