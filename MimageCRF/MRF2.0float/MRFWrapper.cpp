/*
 * MRFWrapper.cpp
 *
 *  Created on: May 6, 2012
 *      Author: bhole
 */

#include <time.h>
#include <math.h>

#include "MRFWrapper.h"

#define m_D(pix,l)  m_D[(pix)*m_nLabels+(l)]
#define m_V(l1,l2)  m_V[(l1)*m_nLabels+(l2)]

#define MIN(a,b)  (((a) < (b)) ? (a) : (b))
#define MAX(a,b)  (((a) > (b)) ? (a) : (b))
#define TRUNCATE_MIN(a,b) { if ((a) > (b)) (a) = (b); }
#define TRUNCATE_MAX(a,b) { if ((a) < (b)) (a) = (b); }
#define TRUNCATE TRUNCATE_MIN




MRFWrapper::MRFWrapper(int width, int height, int nLabels,EnergyFunction *eng):MRF(width,height,nLabels,eng)
{
  m_needToFreeV = 0;
  BPinitializeAlg();
}
MRFWrapper::MRFWrapper(int nPixels, int nLabels,EnergyFunction *eng):MRF(nPixels,nLabels,eng)
{
  m_needToFreeV = 0;
  BPinitializeAlg();
}

MRFWrapper::~MRFWrapper()
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
/*      for (int i=0; i<m_nPixels; i++) {
        delete[] m_messagesNG[i];
      }
      delete[] m_messagesNG;
*/
      delete[] m_neighbors;
      delete[] m_nNeighbors;
    }

}

int MRFWrapper::traverseNeighbors(int pixel1)
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
}

void MRFWrapper::initializeAlg()
{

  if (m_grid_graph)
  {
    printf(" MRFWrapper does not work with grid graphs interface \n");
    exit(1);
  }

  if (!m_grid_graph)
  {
    //general graphs

    // allocate messages
    // this is message number per neighbor
    //int messageNumPerNB = m_nPixels*m_nLabels;
    // true number of messages will be this multiplied with the number of neighbors per pixel-label
    // but in general graph this number of neighbors could be variable

    // m_messagesNG = new REAL*[m_nPixels];
    m_nNeighbors = new int[m_nPixels];

    // need to traverse the neighbor list to get neighbors per pixel

    // int totNB=0;

    m_maxnNB=0;

    // traverse Neighbors list
    for (int m=0; m < m_nPixels; m++) {
      int numNB=traverseNeighbors((int)m);

      m_nNeighbors[m] = numNB;

      if (m_maxnNB<numNB) {
        m_maxnNB=numNB;
      }

      //m_messagesNG[m] = new REAL[numNB*m_nLabels];
      //memset(m_messagesNG[m], 0, numNB*m_nLabels*sizeof(REAL));

      // totNB +=numNB;
    }

    //m_messageArraySizeInBytes = totNB*sizeof(REAL);
  }

}

void MRFWrapper::BPinitializeAlg()
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
    printf(" Not supported \n");
    exit(1);
  }
   // delete[] m_answer;  //chetan
  //delete[] m_scratchMatrix;
  //delete[] nodeArray;

}

MRF::InputType MRFWrapper::getSmoothType()
{
  return m_smoothType;
}

EnergyFunction *MRFWrapper::getEnergyFunction()
{
  return m_e;
}

void MRFWrapper::setExpScale(int expScale)
{
  m_exp_scale = (float)expScale;
}

int MRFWrapper::getNLabels()
{
  return m_nLabels;
}

int MRFWrapper::getWidth()
{
  return m_width;
}

int MRFWrapper::getHeight()
{
  return m_height;
}
FLOATTYPE MRFWrapper::getExpV(int i)
{
  return m_ExpData[i];
}

FLOATTYPE *MRFWrapper::getExpV()
{
  return m_ExpData;
}


MRF::CostVal MRFWrapper::getHorizWeight(int r, int c)
{
  int x = c;
  int y = r;
  int pix    = x+y*m_width;
  return m_varWeights ? m_horizWeights[pix] :  1;
}

MRF::CostVal MRFWrapper::getVertWeight(int r, int c)
{
  int x = c;
  int y = r;
  int pix    = x+y*m_width;
  return  m_varWeights ? m_vertWeights[pix] :  1;
}

bool MRFWrapper::varWeights()
{
  return m_varWeights;
}

FLOATTYPE *MRFWrapper::getScratchMatrix()
{
  return m_scratchMatrix;
}

void MRFWrapper::clearAnswer()
{
  memset(m_answer, 0, m_nPixels*sizeof(Label));
}


void MRFWrapper::setNeighbors(int pixel1, int pixel2, CostVal weight)
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


MRF::EnergyVal MRFWrapper::smoothnessEnergy()
{
  double eng = 0;
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
          while(m_NBTR->hasNext())
          {
            Neighbor *temp1 = (Neighbor *) m_NBTR->next();
            int pix2  = temp1->to_node;
            eng = eng + m_smoothFn(i, pix2, m_answer[i], m_answer[pix2]);
          }
      }
        eng = eng/2;
    }
  }

    return((EnergyVal)eng);

}



MRF::EnergyVal MRFWrapper::dataEnergy()
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


void MRFWrapper::setData(DataCostFn dcost)
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

void MRFWrapper::setData(CostVal* data)
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


void MRFWrapper::setSmoothness(SmoothCostGeneralFn cost)
{
  m_smoothFn = cost;
}

void MRFWrapper::setSmoothness(CostVal* V)
{
  m_type = FIXED_MATRIX;
  m_V = V;
}


void MRFWrapper::setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda)
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


void MRFWrapper::setCues(CostVal* hCue, CostVal* vCue)
{
  m_horizWeights = hCue;
  m_vertWeights  = vCue;
}


inline void SetVectorToZero(MRFWrapper::REAL* to, int K)
{
  MRFWrapper::REAL* to_finish = to + K;
    do
  {
      *to ++ = 0;
  } while (to < to_finish);
}

inline void AddAllMessages(MRFWrapper::REAL* to, MRFWrapper::REAL* from, int K)
{
  MRFWrapper::REAL* to_finish = to + K;
    do
  {
      *to ++ += *from ++;
  } while (to < to_finish);
}

inline void AddDataCost(MRFWrapper::REAL* to, MRFWrapper::REAL* from, int K)
{
  MRFWrapper::REAL* to_finish = to + K;
    do
  {
      *to ++ -= *from ++;  // the negative is because we want to have -ve datacost added
  } while (to < to_finish);
}


int MRFWrapper::findIdxOfNeighbor(int mainpix, int pixel1)
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

int MRFWrapper::findNeighborWithIdx(int mainpix, int idx)
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


void MRFWrapper::extractEdgeList(std::vector<std::vector<std::vector<int> > > *node_edge_list,
                                        std::vector<std::vector<int> > *edge_list)
{
  int edgenumber = 0;

  std::vector<int> temp_list(2,-1);

  for (int i = 0; i < m_nPixels; i++)
  {
    for (int p = 0; p < m_nNeighbors[i]; p++)
    {
      int q=findNeighborWithIdx(i, p);
      if (q >= i)
      {
        // store [other node, edge number]
        temp_list[0] = q;
        temp_list[1] = edgenumber;
        (*node_edge_list)[i].push_back(temp_list);

        // store [other node, edge number]
        temp_list[0] = i;
        (*node_edge_list)[q].push_back(temp_list);

        // store [first node, second node]
        temp_list[1] = q;
        edge_list->push_back(temp_list);

        edgenumber++;
      } // else it means it was put in already
    }
  }
}



// 1 for submodular, 0 for non submodular
int MRFWrapper::checkSubmodularityEdge(const int &node1, const int &node2)
{
  if ( m_grid_graph )
  {
    printf(" Grid graphs interface doesn't work \n");
    exit(1);
  }

  if (! m_grid_graph) // general graph
  {
    if ( m_smoothType != FUNCTION  )
    {
      printf(" Fixed array shouldn't work for now \n");
      exit(1);
    }
    else // smooth is a function
    {
      MRF::CostVal E00 = m_smoothFn(node1, node2, 0, 0);
      MRF::CostVal E01 = m_smoothFn(node1, node2, 0, 1);
      MRF::CostVal E10 = m_smoothFn(node1, node2, 1, 0);
      MRF::CostVal E11 = m_smoothFn(node1, node2, 1, 1);

      if (E00 + E11 <= E10 + E01)
        return 1;
      else
        return 0;
    }
  }
}

void MRFWrapper::extractEdgeType(
    const std::vector<std::vector<int> > &graph_edge_list,
    std::vector<int> *edge_list_type)
{
  // here we assume current graph is flipper

  for (unsigned int e=0; e < graph_edge_list.size(); e++)
  {
    int first_node = graph_edge_list[e][0];
    int second_node = graph_edge_list[e][1];

    if (checkSubmodularityEdge(first_node, second_node) == 1) // submodular
    {
      (*edge_list_type)[e] = 0;
    } else
    { // non submodular
      (*edge_list_type)[e] = 1;
    }
  }
}



void MRFWrapper::optimizeAlg(int nIterations)
{

}
