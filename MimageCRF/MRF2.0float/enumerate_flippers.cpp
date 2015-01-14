/*
 * enumerate_flippers.cpp
 *
 *  Created on: Apr 13, 2013
 *      Author: bhole
 */

#include "mrf.h"
#include "FlipperGeneral.h"
#include "DualDecompositionSubGraph.h"
#include "DualDecompositionTree.h"
#include "DualDecompositionTreeES.h"
#include "TRW-S.h"
#include "MaxProdBP.h"
#include "SyncSumProd.h"
#include "BP-S.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <new>

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

const double primaldualtol = 0.1;

//const int sizeX = 5;
//const int sizeY = 4;
//const int numPixels = 20; //5x4

const int sizeX = 4;   // width
const int sizeY = 4; // height
const int numPixels = 16; // conversion is pix = h*sizeX + w  (h is along height, w along width

//const int sizeX = 2;
//const int sizeY = 1;
//const int numPixels = 2;

const int outerIter = 100;
const int innerIter = 1;

const int numLabels = 2;

MRF::CostVal D[numPixels*numLabels];

// For DD methods
std::vector<int> *flippedinfo;
std::vector<int> *subgraph_node_to_graph_node_g;

float nonsubfncostarray[numLabels][numLabels];

// For Flipper methods
std::vector<int> *flippedinfo_f;

MRF::CostVal fnCost(int pix1, int pix2, int i, int j)
{
/*
  float fncostarray[5][numLabels][numLabels] = {{{15, 20}, {20, 15}},
                                          {{15, 20}, {20, 15}},
                                          {{15, 20}, {20, 15}},
                                          {{15, 20}, {20, 15}},
                                          {{15, 20}, {20, 15}}
                                              };
*/



/*
   nonsubfncostarray[0][0] = (MRF::CostVal)2.3;
   nonsubfncostarray[1][1] = (MRF::CostVal)2.3;
   nonsubfncostarray[0][1] = (MRF::CostVal) 0;
   nonsubfncostarray[1][0] = (MRF::CostVal) 0;
*/


  MRF::CostVal answer = nonsubfncostarray[i][j];

  if (pix1 == 0 && pix2 == 1)
    answer = -answer;

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

  return fnCost(pix1, pix2, i, j);
}



DataCost* computeDataCost()
{
  int dsiIndex = 0;

  for (int p = 0; p < numPixels; p++)
  {
    for (int d = 0; d < numLabels; d++)
    {
      //MRF::CostVal dsiValue = ((MRF::CostVal)(rand() % 100))/10 + 1;
      double r = fabs(sto->Normal(0, 1));
      MRF::CostVal dsiValue = r;
      D[dsiIndex++] =  (float)dsiValue; //((MRF::CostVal)(rand() % 100))/10 + 1;
    }
  }

  return new DataCost (D);
}


void computePairwiseMatrix()
{

  for (int ii=0; ii<numLabels; ii++)
    for (int jj=ii; jj<numLabels; jj++)
    {
      // value was 2.3 fixed before
      double r = fabs(sto->Normal(0, 3)); // is always +ve and maintains all values as non-submodular
      //double r = 100;
      nonsubfncostarray[ii][jj] = nonsubfncostarray[jj][ii] = (ii == jj) ?  (MRF::CostVal)0 : (MRF::CostVal)r;
    }

}




int main()
{
  srand(0);
  sto = new StochasticLib1((int) 0); // time(NULL) ==> returns secs since 1970

  computePairwiseMatrix();

  DataCost *dcost = computeDataCost();
  SmoothnessCost *scost = new SmoothnessCost((MRF::SmoothCostGeneralFn) fnCost);
  EnergyFunction *energy = new EnergyFunction(dcost, scost);

  MRF* mrf = new FlipperGeneral(numPixels, numLabels, energy);

  // you need to set submodular/nonsubmodular function
  ((FlipperGeneral*)mrf)->setFlipperFunction((MRF::SmoothCostGeneralFn) fnCost_flipper);
  // you need to set *flippedinfo before calling that function
  ((FlipperGeneral*)mrf)->setflippedinfo(&flippedinfo_f);


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

  MRF::EnergyVal E;
  float t,tot_t;

  printf("\n*******  Started Flipper *****\n");
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


//     for (int i=0; i<numPixels; i++)
//     {
//     unsigned char rowvalue = mrf->getLabel(i);
//     printf(" Node : %d  state : %d \n", i, (int)rowvalue);
//     }


  delete mrf;



  delete energy;
  delete scost;
  delete dcost;

  return 0;
}


