/*
 * example_dd_mix2.cpp
 *
 *  Created on: May 31, 2012
 *      Author: bhole
 */

// mix of submodular and  non submodular graph (both random and 50/50 graph)

#define FULLRANDOM 1  // set 1 for full random, 0 for 50/50 and 2 for object type, 4 for pairwise different for each edge

#include "mrf.h"
#include "DualDecompositionFlipper.h"
#include "DualDecompositionFlipperES.h"
#include "DualDecompositionSubGraph.h"
#include "DualDecompositionTree.h"
#include "DualDecompositionTreeES.h"
#include "TRW-S.h"
#include "MaxProdBP.h"
#include "SyncSumProd.h"
#include "BP-S.h"
#include "QPBO-v1.3.src/QPBO.cpp"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <new>
#include "common_typedefs.h"

#include "../randgen/stocc/randomc.h"                   // define classes for random number generators
#include "../randgen/stocc/stocc.h"                     // define random library classes

#ifndef MULTIFILE_PROJECT
// If compiled as a single file then include these cpp files,
// If compiled as a project then compile and link in these cpp files.
   #include "../randgen/randomc/mersenne.cpp"             // code for random number generator
   #include "../randgen/stocc/stoc1.cpp"                // random library source code
   #include "../randgen/randomc/userintf.cpp"             // define system specific user interface
#endif

StochasticLib1 *sto;          // make instance of random library



const int sizeX = 1000;
const int sizeY = 1000;
const int numPixels = 1000000;


const int outerIter = 100;
const int innerIter = 1;

const int numLabels = 2;

MRF::CostVal D[numPixels*numLabels];

std::vector<int> *flippedinfo;
std::vector<int> *subgraph_node_to_graph_node_g;

// For Flipper methods
// this will dynamically point to correct flipper graph when Flipper -> optimize() is called
std::vector<int> *flippedinfo_f;

// DualDecompositionFlipperES
myumap *edge_to_subgraphcount_p;

// used as [pixel1pixel2, value]
myumap *mappx;

float fncostarray[numLabels][numLabels];
float nonsubfncostarray[numLabels][numLabels];

myumapvf *mappxvf;

MRF::CostVal fnCost(int pix1, int pix2, int i, int j)
{

  MRF::CostVal answer = 0.0;
  int toggled = 0;

  int smallerpix = pix1;
  int largerpix = pix2;
  if (pix1 > pix2)
  {
    smallerpix = pix2;
    largerpix = pix1;
    toggled = 1;
  }

  if (FULLRANDOM == 0) // 50/50 and is unaffected by the setSubNonSub() function
  {
    if (smallerpix % sizeX <= (sizeX/2))
      answer = fncostarray[i][j];
    else
      answer = nonsubfncostarray[i][j];
  } else if (FULLRANDOM == 1 || FULLRANDOM == 2)
  {

    char buff[30];
    sprintf(buff, "%012d%012d", smallerpix, largerpix);
    std::string str(buff);

    myumap::const_iterator got = (*mappx).find (buff);

    int value = 0;
    if ( got == (*mappx).end() )
    {
      printf("\n Real problem here : smallerpix %d largerpix %d  buff %s\n", smallerpix, largerpix, buff);
      exit(1);
    }
    else
      value = got->second;

    if (value == 0)
      answer = fncostarray[i][j];
    else
      answer = nonsubfncostarray[i][j];

  }

  else if (FULLRANDOM == 4)
    {

      char buff[30];
      sprintf(buff, "%012d%012d", smallerpix, largerpix);
      std::string str(buff);

      myumapvf::const_iterator got = (*mappxvf).find (buff);

      std::vector<float> value;
      if ( got == (*mappxvf).end() )
      {
        printf("\n Real problem here : smallerpix %d largerpix %d  buff %s\n", smallerpix, largerpix, buff);
        exit(1);
      }
      else
        value = got->second;

      if (toggled == 0)
      {
        answer = value[i*2+j];
      } else
        answer = value[j*2+i];

    }

  return answer;
}

MRF::CostVal subfnCost(int pix1, int pix2, int i, int j)
{

  int pixg_i = (*subgraph_node_to_graph_node_g)[pix1];
  int pixg_j = (*subgraph_node_to_graph_node_g)[pix2];

  return fnCost(pixg_i, pixg_j, i, j);
}


// note that flippedinfo is indexed by subgraph node number. So directly use pix1 and pix2
MRF::CostVal nonsubfnCost(int pix1, int pix2, int i, int j)
{

  if ((*flippedinfo)[pix1] == (*flippedinfo)[pix2])
  {
    printf(" Error : Pairs must have one node as flipped !! \n");
    exit(1);
  }

  // we know it is binary for now so 1- solves the flipping needed
  if ((*flippedinfo)[pix1] == 1)
    i = 1-i;

  if ((*flippedinfo)[pix2] == 1)
    j = 1-j;

  int pixg_i = (*subgraph_node_to_graph_node_g)[pix1];
  int pixg_j = (*subgraph_node_to_graph_node_g)[pix2];

  return fnCost(pixg_i, pixg_j, i, j);
}


// this function will be called by DualDecompositionFlipper
// it needs to know the subgraph->graph node matching
MRF::CostVal fnCost_dd_flipper(int pix1, int pix2, int i, int j)
{
  // no flipping needs to be done by DDF

  int pixg_i = (*subgraph_node_to_graph_node_g)[pix1];
  int pixg_j = (*subgraph_node_to_graph_node_g)[pix2];

  // printf("subgraph %d %d, graph %d %d\n", pix1, pix2, pixg_i, pixg_j);

  return fnCost(pixg_i, pixg_j, i, j);
}

// this function will be called by Flipper
// it needs to know what node is flipped, and also needs to
// know the subgraph->graph node matching
MRF::CostVal fnCost_flipper(int pix1, int pix2, int i, int j)
{

  if ((*flippedinfo_f)[pix1] == (*flippedinfo_f)[pix2])
  {
    // should be a submodular edge
  } else if ((*flippedinfo_f)[pix1] != (*flippedinfo_f)[pix2])
  {
    // should be a nonsubmodular edge
  }

  // we know it is binary for now so 1- solves the flipping needed
  if ((*flippedinfo_f)[pix1] == 1)
    i = 1-i;

  if ((*flippedinfo_f)[pix2] == 1)
    j = 1-j;

  int pixg_i = (*subgraph_node_to_graph_node_g)[pix1];
  int pixg_j = (*subgraph_node_to_graph_node_g)[pix2];

  // printf("subgraph %d %d, graph %d %d\n", pix1, pix2, pixg_i, pixg_j);

  return fnCost(pixg_i, pixg_j, i, j);
}


int getSubgraphCount(int pixg_i, int pixg_j)
{
  int pix1, pix2; // pix1 smaller than pix2

  if (pixg_i <= pixg_j)
  {
    pix1 = pixg_i;
    pix2 = pixg_j;
  } else {
    pix2 = pixg_i;
    pix1 = pixg_j;
  }

  char buff[30];
  sprintf(buff, "%012d%012d", pix1, pix2);
  std::string str(buff);

  myumap::const_iterator got = (*edge_to_subgraphcount_p).find (buff);

  int subgraph_count = 0;
  if ( got == (*edge_to_subgraphcount_p).end() )
  {
    printf("\n Real problem here : pix1 %d pix2 %d  buff %s\n", pix1, pix2, buff);
    exit(1);
  }
  else
    subgraph_count = got->second;

  // printf("%d ", subgraph_count);
  return subgraph_count;

}





// this function will be called by DualDecompositionFlipperES
// it needs to know the subgraph->graph node matching
MRF::CostVal fnCost_dd_flipper_es(int pix1, int pix2, int i, int j)
{
  // no flipping needs to be done by DDF

  int pixg_i = (*subgraph_node_to_graph_node_g)[pix1];
  int pixg_j = (*subgraph_node_to_graph_node_g)[pix2];

  int subgraph_count = getSubgraphCount(pixg_i, pixg_j);

  return fnCost(pixg_i, pixg_j, i, j) / double(subgraph_count);
}

// this function will be called by Flipper
// it needs to know what node is flipped, and also needs to
// know the subgraph->graph node matching
MRF::CostVal fnCost_flipper_es(int pix1, int pix2, int i, int j)
{

  if ((*flippedinfo_f)[pix1] == (*flippedinfo_f)[pix2])
  {
    // should be a submodular edge
  } else if ((*flippedinfo_f)[pix1] != (*flippedinfo_f)[pix2])
  {
    // should be a nonsubmodular edge
  }

  // we know it is binary for now so 1- solves the flipping needed
  if ((*flippedinfo_f)[pix1] == 1)
    i = 1-i;

  if ((*flippedinfo_f)[pix2] == 1)
    j = 1-j;

  int pixg_i = (*subgraph_node_to_graph_node_g)[pix1];
  int pixg_j = (*subgraph_node_to_graph_node_g)[pix2];

  // printf("subgraph %d %d, graph %d %d\n", pix1, pix2, pixg_i, pixg_j);

  int subgraph_count = getSubgraphCount(pixg_i, pixg_j);

  return fnCost(pixg_i, pixg_j, i, j) / double(subgraph_count);
}





DataCost* computeDataCost()
{
  int dsiIndex = 0;

  for (int p = 0; p < numPixels; p++)
  {
    for (int d = 0; d < numLabels; d++)
    {
      //MRF::CostVal dsiValue = ((MRF::CostVal)(rand() % 100))/10 + 1;
      double r = 2.0*rand()/RAND_MAX - 1; // fabs(sto->Normal(0, 1));
      MRF::CostVal dsiValue = r;
      D[dsiIndex++] =  (float)dsiValue; //((MRF::CostVal)(rand() % 100))/10 + 1;
    }
  }
  return new DataCost (D);
}

// seems like unordered_map<pair, int> does not work
// and map would be log(n) and hence slow though it works

void setSubNonSubFullRandom()
{
  mappx = new(myumap);

  int pix1 = 0;
  int pix2 = 0;

  for (int j=0; j<sizeY; j++)
  {
    for (int i=0; i<sizeX-1; i++)
    {
      pix1 = j*sizeX+i;
      pix2 = j*sizeX+i+1;
      // note that pix1 is smaller than pix2
      char buff[30];
      sprintf(buff, "%012d%012d", pix1, pix2);
      //printf("%s\n", buff);
      std::string str(buff);

      int value = 0;
      if (rand() % 2 == 0)
        value = 1;
      std::pair<std::string, int> pairsi(buff, value);
      (*mappx).insert(pairsi);
    }
  }

  for (int j=0; j<sizeY-1; j++)
  {
    for (int i=0; i<sizeX; i++)
    {
      pix1 = j*sizeX+i;
      pix2 = (j+1)*sizeX+i;
      // note that pix1 is smaller than pix2
      char buff[30];
      sprintf(buff, "%012d%012d", pix1, pix2);
      //printf("%s\n", buff);
      std::string str(buff);

      int value = 0;
      if (rand() % 2 == 0)
        value = 1;
      std::pair<std::string, int> pairsi(buff, value);
      (*mappx).insert(pairsi);
    }
  }

}


// the function is assuming 1000x1000 matrix
void setSubNonSubObject()
{
  mappx = new(myumap);

  int pix1 = 0;
  int pix2 = 0;

  for (int j=0; j<sizeY; j++)
  {
    for (int i=0; i<sizeX-1; i++)
    {
      pix1 = j*sizeX+i;
      pix2 = j*sizeX+i+1;
      // note that pix1 is smaller than pix2
      char buff[30];
      sprintf(buff, "%012d%012d", pix1, pix2);
      std::string str(buff);

      int value = 0;
      if (i >= 251 && i<= 748 && j>=251 && j<=748) // inner part is submodular
        value = 0;
      else if (i >= 249 && i<= 751 && j>=249 && j<=751) // object boundary
        value = 1;
      else // outer background
        value = 0;

      std::pair<std::string, int> pairsi(buff, value);
      (*mappx).insert(pairsi);
    }
  }

  for (int j=0; j<sizeY-1; j++)
  {
    for (int i=0; i<sizeX; i++)
    {
      pix1 = j*sizeX+i;
      pix2 = (j+1)*sizeX+i;
      // note that pix1 is smaller than pix2
      char buff[30];
      sprintf(buff, "%012d%012d", pix1, pix2);
      std::string str(buff);

      int value = 0;
      if (i >= 251 && i<= 748 && j>=251 && j<=748) // inner part is submodular
        value = 0;
      else if (i >= 249 && i<= 751 && j>=249 && j<=751) // object boundary
        value = 1;
      else // outer background
        value = 0;

      std::pair<std::string, int> pairsi(buff, value);
      (*mappx).insert(pairsi);
    }
  }

}

void computePairwiseMatrix()
{

  for (int ii=0; ii<numLabels; ii++)
    for (int jj=ii; jj<numLabels; jj++)
    {
      // value was 3.7 fixed before
      double r = fabs(sto->Normal(0, 1)); // is always +ve and maintains all values as submodular
      fncostarray[ii][jj] = fncostarray[jj][ii] = (ii == jj) ? (MRF::CostVal)0 : (MRF::CostVal)r;
    }

  for (int ii=0; ii<numLabels; ii++)
    for (int jj=ii; jj<numLabels; jj++)
    {
      // value was 2.3 fixed before
      double r = fabs(sto->Normal(0, 1)); // is always +ve and maintains all values as non-submodular
      nonsubfncostarray[ii][jj] = nonsubfncostarray[jj][ii] = (ii == jj) ? (MRF::CostVal)r : (MRF::CostVal)0;
    }

}


void setFullRandomNotSame()
{
  mappx = new(myumap);
  mappxvf = new(myumapvf);

  int pix1 = 0;
  int pix2 = 0;

  for (int j=0; j<sizeY; j++)
  {
    for (int i=0; i<sizeX-1; i++)
    {
      pix1 = j*sizeX+i;
      pix2 = j*sizeX+i+1;
      // note that pix1 is smaller than pix2
      char buff[30];
      sprintf(buff, "%012d%012d", pix1, pix2);
      //printf("%s\n", buff);
      std::string str(buff);

      int value = 0;
      if (rand() % 2 == 0)
        value = 1;

      std::vector<float> EE(4, 0.0);
      float E00 = 2.0*rand()/RAND_MAX - 1;
      float E01 = 2.0*rand()/RAND_MAX - 1;
      float E10 = 2.0*rand()/RAND_MAX - 1;
      float E11 = 2.0*rand()/RAND_MAX - 1;

      if ((value == 0 && (E00 + E11 <= E01 + E10)) || (value==1 && (E00 + E11 > E01 + E10)))
      {
          EE[0] = E00;
          EE[1] = E01;
          EE[2] = E10;
          EE[3] = E11;
      } else {

          EE[0] = E01;
          EE[1] = E00;
          EE[2] = E11;
          EE[3] = E10;
      }

      std::pair<std::string, int> pairsi(buff, value);
      (*mappx).insert(pairsi);

      std::pair<std::string, std::vector<float> > pairsivf(buff, EE);
      (*mappxvf).insert(pairsivf);

    }
  }

  for (int j=0; j<sizeY-1; j++)
  {
    for (int i=0; i<sizeX; i++)
    {
      pix1 = j*sizeX+i;
      pix2 = (j+1)*sizeX+i;
      // note that pix1 is smaller than pix2
      char buff[30];
      sprintf(buff, "%012d%012d", pix1, pix2);
      //printf("%s\n", buff);
      std::string str(buff);

      int value = 0;
      if (rand() % 2 == 0)
        value = 1;

      std::vector<float> EE(4, 0.0);
      float E00 = 2.0*rand()/RAND_MAX - 1;
      float E01 = 2.0*rand()/RAND_MAX - 1;
      float E10 = 2.0*rand()/RAND_MAX - 1;
      float E11 = 2.0*rand()/RAND_MAX - 1;

      if ((value == 0 && (E00 + E11 <= E01 + E10)) || (value==1 && (E00 + E11 > E01 + E10)))
      {
          EE[0] = E00;
          EE[1] = E01;
          EE[2] = E10;
          EE[3] = E11;
      } else {

          EE[0] = E01;
          EE[1] = E00;
          EE[2] = E11;
          EE[3] = E10;
      }


      std::pair<std::string, int> pairsi(buff, value);
      (*mappx).insert(pairsi);

      std::pair<std::string, std::vector<float> > pairsivf(buff, EE);
      (*mappxvf).insert(pairsivf);


    }
  }

}



typedef double REAL;

void QPBOAddTerms(QPBO<REAL>* q, EnergyFunction *energy) {
  // Assuming the data is ARRAY type i.e. ( m_dataType == ARRAY)
  for (int i = 0; i < numPixels; ++i) {
    double E0 = energy->m_dataCost->getDataCostRaw(2*i);
    double E1 = energy->m_dataCost->getDataCostRaw(2*i+1);
    // printf(" Node : %d , Energy : %f %f \n", i, E0, E1);
    q->AddUnaryTerm(i, E0, E1);
  }

  for (int j=0; j<sizeY; j++)
  {
    for (int i=0; i<sizeX-1; i++)
    {
      double E00 = energy->m_smoothCost->getSmoothnessCostRaw(j*sizeX+i, j*sizeX+i+1, 0, 0);
      double E01 = energy->m_smoothCost->getSmoothnessCostRaw(j*sizeX+i, j*sizeX+i+1, 0, 1);
      double E10 = energy->m_smoothCost->getSmoothnessCostRaw(j*sizeX+i, j*sizeX+i+1, 1, 0);
      double E11 = energy->m_smoothCost->getSmoothnessCostRaw(j*sizeX+i, j*sizeX+i+1, 1, 1);
      // printf(" Node : %d %d , Energy : %f %f %f %f \n", j*sizeX+i, j*sizeX+i+1, E00, E01, E10, E11);
      q->AddPairwiseTerm(j*sizeX+i, j*sizeX+i+1, E00, E01, E10, E11);
    }
  }

  for (int j=0; j<sizeY-1; j++)
  {
    for (int i=0; i<sizeX; i++)
    {
      double E00 = energy->m_smoothCost->getSmoothnessCostRaw(j*sizeX+i, (j+1)*sizeX+i, 0, 0);
      double E01 = energy->m_smoothCost->getSmoothnessCostRaw(j*sizeX+i, (j+1)*sizeX+i, 0, 1);
      double E10 = energy->m_smoothCost->getSmoothnessCostRaw(j*sizeX+i, (j+1)*sizeX+i, 1, 0);
      double E11 = energy->m_smoothCost->getSmoothnessCostRaw(j*sizeX+i, (j+1)*sizeX+i, 1, 1);
      // printf(" Node : %d %d , Energy : %f %f %f %f \n", j*sizeX+i, (j+1)*sizeX+i, E00, E01, E10, E11);
      q->AddPairwiseTerm(j*sizeX+i, (j+1)*sizeX+i, E00, E01, E10, E11);
    }
  }
}



int main()
{
  //srand(time(NULL));
  srand(100);
  sto = new StochasticLib1((int) 0) ; //(int) time(NULL));

  computePairwiseMatrix();

  DataCost *dcost = computeDataCost();
  SmoothnessCost *scost = new SmoothnessCost((MRF::SmoothCostGeneralFn) fnCost);
  // SmoothnessCost *nsscost = new SmoothnessCost((MRF::SmoothCostGeneralFn) nonsubfnCost); // nonsubmodular
  EnergyFunction *energy = new EnergyFunction(dcost, scost);

  // generating a graph where the submod and nonsubmod functions will be fixed so the same
  // graph is used for all algorithms
  if (FULLRANDOM == 1)
  {
    setSubNonSubFullRandom();
  }
  else if (FULLRANDOM == 0)
  {} // 50/50
  else if (FULLRANDOM == 2) // object
  {
    setSubNonSubObject();
  }
  else if (FULLRANDOM == 4)
  {
    setFullRandomNotSame();
  }


  bool runDDF = false;
  bool runDDG = false;
  bool runDDT = false;
  bool runDDTes = false;
  bool runDDFes = false;
  bool runTRWS = true;
  bool runMaxProdBP = false;
  bool runBPS       = false;
  bool runSSP       = false;
  bool runQPBO  = true;



  if (runDDF) {

    MRF::EnergyVal E;
    float t,tot_t;

    MRF* mrf = new DualDecompositionFlipper(numPixels, numLabels, energy);

    ((DualDecompositionFlipper*)mrf)->setSplitType(TREESADD);

    // you need to set DDFlipper function
    ((DualDecompositionFlipper*)mrf)->setDDFlipperFunction((MRF::SmoothCostGeneralFn) fnCost_dd_flipper);
    // you need to set DDFlipper function
    ((DualDecompositionFlipper*)mrf)->setFlipperFunction((MRF::SmoothCostGeneralFn) fnCost_flipper);
    // you need to set *flippedinfo before calling that function
    ((DualDecompositionFlipper*)mrf)->setflippedinfo(&flippedinfo_f);
    // you need to set the subgraph node to graph node mapping
    ((DualDecompositionFlipper*)mrf)->setNodeMappingStructure(&subgraph_node_to_graph_node_g);



    for (int j=0; j<sizeY; j++)
    {
      for (int i=0; i<sizeX-1; i++)
      {
        mrf->setNeighbors(j*sizeX+i, j*sizeX+i+1, 1);
      }
    }

    for (int j=0; j<sizeY-1; j++)
    {
      for (int i=0; i<sizeX; i++)
      {
        mrf->setNeighbors(j*sizeX+i, (j+1)*sizeX+i, 1);
      }
    }

    printf("\n*******  Started Dual Decomposition Flipper *****\n");
    mrf->initialize();
    mrf->clearAnswer();

    E = mrf->totalEnergy();
    printf("Energy at the Start= %g (%g,%g)\n", (float)E,
      (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

    tot_t = 0;
    mrf->optimize(1, t);
    E = mrf->totalEnergy();


    tot_t = tot_t + t ;
    printf("energy = %g (%f secs)\n", (float)E, tot_t);



//  for (int i=0; i<numPixels; i++)
//  {
//    unsigned char rowvalue = mrf->getLabel(i);
//    printf(" Node : %d  state : %d \n", i, (int)rowvalue);
//  }



    delete mrf;


  }





  if (runDDFes) {

    MRF::EnergyVal E;
    float t,tot_t;

    MRF* mrf = new DualDecompositionFlipperES(numPixels, numLabels, energy);

    ((DualDecompositionFlipperES*)mrf)->setSplitType(DDFES::TREESADD_EDGESHARE);

    // you need to set DDFlipper function
    ((DualDecompositionFlipperES*)mrf)->setDDFlipperFunction((MRF::SmoothCostGeneralFn) fnCost_dd_flipper_es);
    // you need to set DDFlipper function
    ((DualDecompositionFlipperES*)mrf)->setFlipperFunction((MRF::SmoothCostGeneralFn) fnCost_flipper_es);
    // you need to set *flippedinfo before calling that function
    ((DualDecompositionFlipperES*)mrf)->setflippedinfo(&flippedinfo_f);
    // you need to set the subgraph node to graph node mapping
    ((DualDecompositionFlipperES*)mrf)->setNodeMappingStructure(&subgraph_node_to_graph_node_g);
    // you need to set the edge sharing map datastructure
    ((DualDecompositionFlipperES*)mrf)->set_edge_to_subgraphcount(&edge_to_subgraphcount_p);



    for (int j=0; j<sizeY; j++)
    {
      for (int i=0; i<sizeX-1; i++)
      {
        //if (j*sizeX+i == 2 || j*sizeX+i == 12){} // || j*sizeX+i == 22) {}
        //else
        mrf->setNeighbors(j*sizeX+i, j*sizeX+i+1, 1);
      }
    }

    for (int j=0; j<sizeY-1; j++)
    {
      for (int i=0; i<sizeX; i++)
      {
        mrf->setNeighbors(j*sizeX+i, (j+1)*sizeX+i, 1);
      }
    }

    printf("\n*******  Started Dual Decomposition subgraph *****\n");
    mrf->initialize();
    mrf->clearAnswer();

    E = mrf->totalEnergy();
    printf("Energy at the Start= %g (%g,%g)\n", (float)E,
      (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

    tot_t = 0;
    mrf->optimize(1, t);
    E = mrf->totalEnergy();


    tot_t = tot_t + t ;
    printf("energy = %g (%f secs)\n", (float)E, tot_t);


/*
  for (int i=0; i<numPixels; i++)
  {
    unsigned char rowvalue = mrf->getLabel(i);
    printf(" Node : %d  state : %d \n", i, (int)rowvalue);
  }
*/
    delete mrf;
  }






  if (runDDG) {

    MRF::EnergyVal E;
    float t,tot_t;

    MRF* mrf = new DualDecompositionSubGraph(numPixels, numLabels, energy);

    // the functions below are so that the functions work on the subgraphs
    // you need to set submodular function
    ((DualDecompositionSubGraph*)mrf)->setSubmodularFunction((MRF::SmoothCostGeneralFn) subfnCost);
    // you need to set non submodular function
    ((DualDecompositionSubGraph*)mrf)->setNonSubmodularFunction((MRF::SmoothCostGeneralFn) nonsubfnCost);
    // you need to set the subgraph node to graph node mapping
    ((DualDecompositionSubGraph*)mrf)->setNodeMappingStructure(&subgraph_node_to_graph_node_g);
    // you also need to make sure energy calls it correctly ie during computation of total cost in the inside subgraph
    // you need to set *flippedinfo before calling that function
    ((DualDecompositionSubGraph*)mrf)->setflippedinfo(&flippedinfo);

    for (int j=0; j<sizeY; j++)
    {
      for (int i=0; i<sizeX-1; i++)
      {
        mrf->setNeighbors(j*sizeX+i, j*sizeX+i+1, 1);
      }
    }

    for (int j=0; j<sizeY-1; j++)
    {
      for (int i=0; i<sizeX; i++)
      {
        mrf->setNeighbors(j*sizeX+i, (j+1)*sizeX+i, 1);
      }
    }

    printf("\n*******  Started Dual Decomposition subgraph *****\n");
    mrf->initialize();
    mrf->clearAnswer();

    E = mrf->totalEnergy();
    printf("Energy at the Start= %g (%g,%g)\n", (float)E,
      (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

    tot_t = 0;
    mrf->optimize(1, t);
    E = mrf->totalEnergy();
    tot_t = tot_t + t ;
    printf("energy = %g (%f secs)\n", (float)E, tot_t);

/*
  for (int i=0; i<numPixels; i++)
  {
    unsigned char rowvalue = mrf->getLabel(i);
    printf(" Node : %d  state : %d \n", i, (int)rowvalue);
  }
*/


    delete mrf;

  }


  ///////////////////////////////////////////////
  // DD Tree
  //////////////////////////////////////////////

  if (runDDT)
  {
    MRF::EnergyVal E;
    float t,tot_t;

    printf("\n*******Started Dual Decomposition *****\n");
    MRF* mrf = new DualDecompositionTree(sizeX*sizeY, numLabels, energy);

    for (int j=0; j<sizeY; j++)
    {
      for (int i=0; i<sizeX-1; i++)
      {
        mrf->setNeighbors(j*sizeX+i, j*sizeX+i+1, 1);
      }
    }

    for (int j=0; j<sizeY-1; j++)
    {
      for (int i=0; i<sizeX; i++)
      {
        mrf->setNeighbors(j*sizeX+i, (j+1)*sizeX+i, 1);
      }
    }

    mrf->initialize();
    mrf->clearAnswer();

    E = mrf->totalEnergy();
    printf("Energy at the Start= %g (%g,%g)\n", (float)E,
        (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

    tot_t = 0;

    mrf->optimize(innerIter, t);

    E = mrf->totalEnergy();
    tot_t = tot_t + t ;
    printf("energy = %g (%f secs)\n", (float)E, tot_t);


    delete mrf;

  }


  ///////////////////////////////////////////////
  // DD Tree Edge sharing
  //////////////////////////////////////////////

  if (runDDTes)
  {
    MRF::EnergyVal E;
    float t,tot_t;

    printf("\n*******Started Dual Decomposition *****\n");
    MRF* mrf = new DualDecompositionTreeES(sizeX*sizeY, numLabels, energy);

    for (int j=0; j<sizeY; j++)
    {
      for (int i=0; i<sizeX-1; i++)
      {
        mrf->setNeighbors(j*sizeX+i, j*sizeX+i+1, 1);
      }
    }

    for (int j=0; j<sizeY-1; j++)
    {
      for (int i=0; i<sizeX; i++)
      {
        mrf->setNeighbors(j*sizeX+i, (j+1)*sizeX+i, 1);
      }
    }

    mrf->initialize();
    mrf->clearAnswer();

    E = mrf->totalEnergy();
    printf("Energy at the Start= %g (%g,%g)\n", (float)E,
        (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

    tot_t = 0;

    mrf->optimize(innerIter, t);

    E = mrf->totalEnergy();
    tot_t = tot_t + t ;
    printf("energy = %g (%f secs)\n", (float)E, tot_t);


    delete mrf;

  }






  ////////////////////////////////////////////////
  //                  TRW-S                     //
  ////////////////////////////////////////////////
  if (runTRWS)
  {
    printf("\n*******Started TRW-S *****\n");
    double lowerBound;
    MRF* mrf = new TRWS(sizeX, sizeY, numLabels, energy);

    mrf->initialize();
    mrf->clearAnswer();

    MRF::EnergyVal E;
    float t,tot_t;

    E = mrf->totalEnergy();
    printf("Energy at the Start= %g (%g,%g)\n", (float)E,
        (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

    tot_t = 0;
    // 200
    for (int iter=0; iter<10; iter++)
    {
      mrf->optimize(10, t);

      E = mrf->totalEnergy();
      lowerBound = mrf->lowerBound();
      tot_t = tot_t + t ;
      printf("energy = %g, lower bound = %f (%f secs)\n", (float)E, lowerBound, tot_t);
    }

/*
      for (int i=0; i<numPixels; i++)
      {
        unsigned char rowvalue = mrf->getLabel(i);
        printf(" Node : %d  state : %d \n", i, (int)rowvalue);
      }
*/

    delete mrf;
  }



  ////////////////////////////////////////////////
  //          Belief Propagation                //
  ////////////////////////////////////////////////
  if (runMaxProdBP)
  {
    MRF::EnergyVal E;
    float t,tot_t;

    printf("\n*******  Started MaxProd Belief Propagation *****\n");
    MRF*  mrf = new MaxProdBP(sizeX,sizeY,numLabels,energy);
    mrf->initialize();
    mrf->clearAnswer();

    E = mrf->totalEnergy();
    printf("Energy at the Start= %g (%g,%g)\n", (float)E,
       (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

    tot_t = 0;
    for (int iter=0; iter < 400; iter++)
    {
      mrf->optimize(1, t);

      E = mrf->totalEnergy();
      tot_t = tot_t + t ;
      printf("energy = %g (%f secs)\n", (float)E, tot_t);
    }

    delete mrf;
  }


  ////////////////////////////////////////////////
  //                  BP-S                     //
  ////////////////////////////////////////////////
  if (runBPS)
  {
    MRF::EnergyVal E;
    float t,tot_t;

    printf("\n*******Started BP-S *****\n");
    MRF*  mrf = new BPS(sizeX,sizeY,numLabels,energy);

    // can disable caching of values of general smoothness function:
    //mrf->dontCacheSmoothnessCosts();

    mrf->initialize();
    mrf->clearAnswer();

    E = mrf->totalEnergy();
    printf("Energy at the Start= %g (%g,%g)\n", (float)E,
       (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

    tot_t = 0;
    for (int iter=0; iter<200; iter++)
    {
      mrf->optimize(10, t);

      E = mrf->totalEnergy();
      tot_t = tot_t + t ;
      printf("energy = %g (%f secs)\n", (float)E, tot_t);
    }

    delete mrf;
  }

  ////////////////////////////////////////////////
  //                  SyncSumProd                     //
  ////////////////////////////////////////////////
  if (runSSP)
  {
    MRF::EnergyVal E;
    float t,tot_t;

    printf("\n*******Started SyncSumProd *****\n");
    MRF*  mrf = new SyncSumProd(sizeX,sizeY,numLabels,energy);

    // can disable caching of values of general smoothness function:
    //mrf->dontCacheSmoothnessCosts();

    mrf->initialize();
    mrf->clearAnswer();

    E = mrf->totalEnergy();
    printf("Energy at the Start= %g (%g,%g)\n", (float)E,
       (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

    tot_t = 0;
    for (int iter=0; iter<outerIter; iter++)
    {
      mrf->optimize(innerIter, t);

      E = mrf->totalEnergy();
      tot_t = tot_t + t ;
      printf("energy = %g (%f secs)\n", (float)E, tot_t);
    }

    delete mrf;
  }


  if (runQPBO)
  {

    printf("\n*******Started QPBO *****\n");
    QPBO<REAL>* q;
    // start timer
    clock_t start = clock();
    // grid is sizeX by size Y
    // so number of edges are sizeX*(sizeY-1) + sizeY*(sizeX-1)
    int numEdges = sizeX*(sizeY-1) + sizeY*(sizeX-1);
    q = new QPBO<REAL>(numPixels, numEdges); // max number of nodes & edges
    q->AddNode(numPixels); // add two nodes

    QPBOAddTerms(q, energy);

    q->Solve();

    int count1 = 0;
    for (int i=0; i<numPixels; i++)
    {
      int rowvalue = q->GetLabel(i);
      if (rowvalue < 0)
        count1++;
      // printf(" Node : %d  state : %d \n", i, (int)rowvalue);
    }
    printf(" count = %d \n", count1);

    q->ComputeWeakPersistencies();

    // stop timer
    clock_t finish = clock();
    float timed = (float) (((double)(finish - start)) / CLOCKS_PER_SEC);

    int count = 0;
    for (int i=0; i<numPixels; i++)
    {
      int rowvalue = q->GetLabel(i);
      if (rowvalue < 0)
        count++;
      // printf(" Node : %d  state : %d \n", i, (int)rowvalue);
    }
    printf(" count = %d \n", count);

    printf("energy = %g, lower bound = %f (%f secs)\n",
        (float)q->ComputeTwiceEnergy()/2.0, q->ComputeTwiceLowerBound()/2.0, timed);
  }



  delete energy;

  return 0;
}


