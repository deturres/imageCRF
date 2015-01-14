/*
 * MaxProdBPTree.cpp
 *
 *  Created on: Oct 3, 2011
 *      Author: bhole
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "MaxProdBPTree.h"
#include "regions-new.h"

#define m_D(pix,l)  m_D[(pix)*m_nLabels+(l)]
#define m_V(l1,l2)  m_V[(l1)*m_nLabels+(l2)]

#define MIN(a,b)  (((a) < (b)) ? (a) : (b))
#define MAX(a,b)  (((a) > (b)) ? (a) : (b))
#define TRUNCATE_MIN(a,b) { if ((a) > (b)) (a) = (b); }
#define TRUNCATE_MAX(a,b) { if ((a) < (b)) (a) = (b); }
#define TRUNCATE TRUNCATE_MIN




MaxProdBPTree::MaxProdBPTree(int nPixels, int nLabels,EnergyFunction *eng):MRF(nPixels,nLabels,eng)
{
	m_needToFreeV = 0;
	BPinitializeAlg();
}

MaxProdBPTree::~MaxProdBPTree()
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

int MaxProdBPTree::traverseNeighbors(int pixel1)
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

void MaxProdBPTree::initializeAlg()
{

  if (m_grid_graph) {
    printf(" MaxProdBPTree should not handle grid graphs \n");
    exit(1);
  }
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

void MaxProdBPTree::BPinitializeAlg()
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
	  printf(" MaxProdBPTree should not handle grid graphs \n");
	  exit(1);
	}

	// delete[] m_answer;  //chetan
	//delete[] m_scratchMatrix;
	//delete[] nodeArray;

}

MRF::InputType MaxProdBPTree::getSmoothType()
{
  return m_smoothType;
}


MRF::SmoothCostGeneralFn MaxProdBPTree::getSmoothnessFn()
{
  return m_smoothFn;
}

EnergyFunction *MaxProdBPTree::getEnergyFunction()
{
  return m_e;
}

void MaxProdBPTree::setExpScale(int expScale)
{
  m_exp_scale = (float)expScale;
}

int MaxProdBPTree::getNLabels()
{
  return m_nLabels;
}

FLOATTYPE MaxProdBPTree::getExpV(int i)
{
  return m_ExpData[i];
}

FLOATTYPE *MaxProdBPTree::getExpV()
{
  return m_ExpData;
}


bool MaxProdBPTree::varWeights()
{
  return m_varWeights;
}

FLOATTYPE *MaxProdBPTree::getScratchMatrix()
{
  return m_scratchMatrix;
}

void MaxProdBPTree::clearAnswer()
{
	memset(m_answer, 0, m_nPixels*sizeof(Label));
}


void MaxProdBPTree::setNeighbors(int pixel1, int pixel2, CostVal weight)
{
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



MRF::EnergyVal MaxProdBPTree::call_smoothFn(int pix_i, int pix_j, int label_i, int label_j)
{
  return m_smoothFn(pix_i, pix_j, label_i, label_j);
}


MRF::EnergyVal MaxProdBPTree::smoothnessEnergy()
{
    double eng =  0;
    double weight;
    int x,y,pix;

    if ( m_grid_graph )
    {
      printf(" MaxProdBPTree should not handle grid graphs \n");
      exit(1);
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
		    		eng = eng + call_smoothFn(i, pix2, m_answer[i], m_answer[pix2]);

		    	}
			}
		    eng = eng/2;
		}
	}

    return((EnergyVal) eng);

}



MRF::EnergyVal MaxProdBPTree::dataEnergy()
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
    return((EnergyVal)eng);

}


void MaxProdBPTree::setData(DataCostFn dcost)
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

void MaxProdBPTree::setData(CostVal* data)
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


void MaxProdBPTree::updateData(CostVal* data)
{

  // FloatType *cData = m_ExpData;
  for (int i= 0; i < m_nPixels; i++)
  {
    FloatType *cData = nodeArray[i].localEv; // = cData;
    for(int j = 0; j < m_nLabels; j++)
    {
      *cData = (float)m_D(i,j);
      cData++;
    }
  }
}




void MaxProdBPTree::setSmoothness(SmoothCostGeneralFn cost)
{
  m_smoothFn = cost;
}
void MaxProdBPTree::setSmoothness(CostVal* V)
{
  m_type = FIXED_MATRIX;
  m_V = V;
}


void MaxProdBPTree::setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda)
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


inline void SetVectorToZero(MaxProdBPTree::REAL* to, int K)
{
	MaxProdBPTree::REAL* to_finish = to + K;
    do
	{
	    *to ++ = 0;
	} while (to < to_finish);
}

inline void AddAllMessages(MaxProdBPTree::REAL* to, MaxProdBPTree::REAL* from, int K)
{
	MaxProdBPTree::REAL* to_finish = to + K;
    do
	{
	    *to ++ += *from ++;
	} while (to < to_finish);
}

inline void AddDataCost(MaxProdBPTree::REAL* to, MaxProdBPTree::REAL* from, int K)
{
	MaxProdBPTree::REAL* to_finish = to + K;
    do
	{
	    *to ++ -= *from ++;  // the negative is because we want to have -ve datacost added
	} while (to < to_finish);
}


int MaxProdBPTree::findIdxOfNeighbor(int mainpix, int pixel1)
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

int MaxProdBPTree::findNeighborWithIdx(int mainpix, int idx)
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


void MaxProdBPTree::extractParentArrayUsingBFS(const int& root)
{
  parent_array_.resize(m_nPixels, std::vector<int>(2, -1));

  std::vector<int> visited(m_nPixels, 0); // 0 means not visited, 1 means visited.

  // simple queue
  std::vector<int> queue(m_nPixels, -1); // -1 means that element is not filled yet.

  int cnt = 0;

  int m = root;

  // initialize root to have no parent.
  parent_array_[cnt][0] = -1;
  parent_array_[cnt][1] = m;
  cnt++;

  queue[0] = m;
  visited[0] = 1;
  int start_queue = 0;
  int end_queue = 0;

  while (end_queue >= start_queue)
  {
	m = queue[start_queue];
	start_queue++;
    for (int p=0; p<m_nNeighbors[m]; p++)
    {
      int q=findNeighborWithIdx(m, p);
      if (q==-1)
      {
	    printf("Error\n");
		exit(1);
	  }
      // just checking that the graph/tree is constructed symmetrically.
	  int idx = findIdxOfNeighbor(q, m);
	  if (idx==-1)
	  {
	    printf("Error\n");
	    exit(1);
	  }
	  if (visited[q] == 0)
	  {
	    end_queue++;
	    queue[end_queue] = q;
	    parent_array_[cnt][0] = m;
	    parent_array_[cnt][1] = q;
	    cnt++;
	    visited[q] = 1;
	  }
    }
  }
}


void MaxProdBPTree::optimizeAlg(int nIterations)
{
  if ( m_grid_graph )
  {
    printf(" MaxProdBPTree should not handle grid graphs \n");
    exit(1);
  } else {
    // for non-grid graphs

    // set root as first node
    // extract parent array
    // do BFS for extraction of parent array
    // (nodes are labeled from parent followed by all children of that parent and so on)
    // for i=leaves to root
    //   pass message to parents
    // collect score at root
    // get indices for optimal framework.

    const int root = 0;
    extractParentArrayUsingBFS(root);
    // printf("Extracted parent array \n");
    // for (int i=0; i<m_nPixels; i++)
    // {
    //   printf("%d %d\n", parent_array_[i][0], parent_array_[i][1]);
    // }

    /*
    // check that at least last element is child.
    // is [m_nPixels-1][0] is non -1
    // and [m_nPixels-1][1] is -1
    if (parent_array_[m_nPixels-1][0] == -1 || parent_array_[m_nPixels-1][1] != -1)
    {
      printf("Error : at least last element should be child \n ");
      exit(1);
    }
    */

    // check that root has parent -1
    if (parent_array_[0][0] != -1)
    {
      printf(" Error : root should not have a parent \n");
      exit(1);
    }

    REAL* Di;
    Di = new REAL[m_nLabels];

    REAL* Di2;
    Di2 = new REAL[m_nLabels];

    REAL* M_ptr, *N_ptr;
    REAL delta;

    std::vector<std::vector<int> > back_indexes(m_nPixels, std::vector<int>(m_nLabels, -1));

    const CostVal lambda = 1;

    // do for all nodes including root.
    for (int m = m_nPixels-1; m >= 0; m--)
    {
      // get child node and parent node.
      const int current = parent_array_[m][1];
      const int parent = parent_array_[m][0];
      SetVectorToZero(Di, m_nLabels);
      if (current != -1) // non-leaf node
      {
  	    M_ptr = m_messagesNG[current];
  	    for (int p = 0; p < m_nNeighbors[current]; p++, M_ptr+=m_nLabels)
  	    {
  		  // add messages from children only
  	      int q=findNeighborWithIdx(current, p);
  	      if (q != parent)
  	      {
  	        AddAllMessages(Di, M_ptr, m_nLabels);
  	      }
  	    }
  	  }
  	  AddDataCost(Di, nodeArray[current].localEv, m_nLabels);
  	  // printf("\nnode : %d\n", current);
  	  // for (int kk=0; kk<m_nLabels; kk++)
  	  // {
  		//   printf("Di[%d] = %f    ", kk, Di[kk]);
  	  // }

  	  // pass message to parent.
  	  if (parent != -1) // checking current is not the root.
  	  {
  	   	int idx = findIdxOfNeighbor(parent, current);
  	   	if (idx==-1)
  	   	{
  	   	  printf("Error\n");
  	   	  exit(1);
  	   	}
  	   	// current --> parent
  	   	N_ptr = m_messagesNG[parent];
  	   	N_ptr += idx*m_nLabels;

  	   	if ( m_smoothType != FUNCTION  )
  	   	{ //assume ARRAY
  	   	  printf(" Not coded yet !");
  	   	  exit(0);
  	   	} else { // this is smoothfunction is function()
  	   	  // UpdateMessageGENERAL_NG(M_ptr, Di, K, 1, lambda, m, q, buf, N_ptr);

  	   	  for (int ki=0; ki<m_nLabels; ki++)
  	   	  {
  	   	    Di2[ki] = Di[ki];
  	   	  }

  	   	  for (int kj=0; kj<m_nLabels; kj++)
  	   	  {
  	   	    N_ptr[kj] = Di2[0] - lambda*call_smoothFn(current, parent, 0, kj);
  	   	    int maxidx = 0;
  	   	    for (int ki=1; ki<m_nLabels; ki++)
  	   	    {
  	   	      // TRUNCATE_MAX(N_ptr[kj], Di2[ki] - lambda*m_smoothFn(current, parent, ki, kj)); //V[0]);
  	   	      if ((N_ptr[kj]) < (Di2[ki] - lambda*call_smoothFn(current, parent, ki, kj))) {
  	   	        (N_ptr[kj]) = (Di2[ki] - lambda*call_smoothFn(current, parent, ki, kj));
  	   	        maxidx = ki;
  	   	      }
  	   	    }
  	   	   back_indexes[current][kj] = maxidx;
  	   	  }
  	   	  //delta = N_ptr[0];
  	   	  //for (int kj=1; kj<m_nLabels; kj++) TRUNCATE(delta, N_ptr[kj]);
  	   	  //for (int kj=0; kj<m_nLabels; kj++) N_ptr[kj] -= delta;
  	   	}
  	  }
    }

    // find best score at root.
    // Di contains the root overall values because it will be the last one to be processed.
    int bestIdx = 0;
    double max_value = Di[0];
    for (int i = 1; i < m_nLabels; i++)
    {
      if (max_value < Di[i])
      {
        max_value = Di[i];
        bestIdx = i;
      }
    }
    m_answer[0] = bestIdx;

    delete[] Di;
    delete[] Di2;

    // back track and put solutions in mrf
    for (int m = 1; m < m_nPixels; m++)
    {
      // get child node and parent node.
      const int current = parent_array_[m][1];
      const int parent = parent_array_[m][0];
      // for (int mmm = 0; mmm < m_nLabels; mmm++)
      //   printf(" vVVv %d ", back_indexes[current][m_answer[parent]]);
      // printf("\n");
      m_answer[current] = back_indexes[current][m_answer[parent]];
    }
  }
}



