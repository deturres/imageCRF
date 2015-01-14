/*
 * example_dd_subgraph.cpp
 *
 *  Created on: Dec 28, 2011
 *      Author: bhole
 */

#include "mrf.h"
#include "DualDecompositionSubGraph.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <new>


const int numPixels = 20; //5x4
const int numLabels = 2;

MRF::CostVal D[numPixels*numLabels];

std::vector<int> *flippedinfo;
std::vector<int> *subgraph_node_to_graph_node_g;

MRF::CostVal fnCost(int pix1, int pix2, int i, int j)
{

  double fncostarray[5][numLabels][numLabels] = {{{15, 20}, {20, 15}},
		  	  	  	  	  	  	  	  	  	  {{15, 20}, {20, 15}},
		  	  	  	  	  	  	  	  	  	  {{15, 20}, {20, 15}},
		  	  	  	  	  	  	  	  	  	  {{15, 20}, {20, 15}},
		  	  	  	  	  	  	  	  	  	  {{15, 20}, {20, 15}}
  	  	  	  	  	  	  	  	  	  	  	  };

  double nonsubfncostarray[5][numLabels][numLabels] = {{{20, 15}, {15, 20}},
                                          {{20, 15}, {15, 20}},
                                          {{20, 15}, {15, 20}},
                                          {{20, 15}, {15, 20}},
                                          {{20, 15}, {15, 20}}
                                              };


  int arrayidx = (pix1+pix2) % 5;

  int smallerpix = pix1;
  if (pix1 > pix2)
  {
    smallerpix = pix2;
  }

  MRF::CostVal answer = 0.0;

  if (smallerpix % 5 <= 2)
    answer = fncostarray[arrayidx][i][j];
  else
    answer = nonsubfncostarray[arrayidx][i][j];

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





DataCost* computeDataCost()
{
  double datac[numPixels][numLabels] = {{10,15}, {20,15}, {15,10}, {20,15}, {10,15},
      {10,15}, {20,15}, {15,10}, {20,15}, {10,15},
      {10,15}, {20,15}, {15,10}, {20,15}, {10,15},
      {10,15}, {20,15}, {15,10}, {20,15}, {10,15}};
  int dsiIndex = 0;

  for (int p = 0; p < numPixels; p++)
  {
    for (int d = 0; d < numLabels; d++)
    {
      // MRF::CostVal dsiValue = (((p-d)*(p-d)) % 4)*(float(p+d) / 3) * 0.5 + 1.5;
      MRF::CostVal dsiValue = datac[p][d];
      D[dsiIndex++] =  (float)dsiValue; //((MRF::CostVal)(rand() % 100))/10 + 1;
      // printf(" p = %d , d = %d, dsiValue = %f \n", p, d, D[dsiIndex-1]);
    }
  }
  return new DataCost (D);
}



int main()
{
  DataCost *dcost = computeDataCost();
  SmoothnessCost *scost = new SmoothnessCost((MRF::SmoothCostGeneralFn) fnCost);
  // SmoothnessCost *nsscost = new SmoothnessCost((MRF::SmoothCostGeneralFn) nonsubfnCost); // nonsubmodular
  EnergyFunction *energy = new EnergyFunction(dcost, scost);
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


  mrf->setNeighbors(0, 1, 1);
  mrf->setNeighbors(1, 2, 1);
  mrf->setNeighbors(2, 3, 1);
  mrf->setNeighbors(3, 4, 1);
  mrf->setNeighbors(5, 6, 1);
  mrf->setNeighbors(6, 7, 1);
  mrf->setNeighbors(7, 8, 1);
  mrf->setNeighbors(8, 9, 1);
  mrf->setNeighbors(10, 11, 1);
  mrf->setNeighbors(11, 12, 1);
  mrf->setNeighbors(12, 13, 1);
  mrf->setNeighbors(13, 14, 1);
  mrf->setNeighbors(15, 16, 1);
  mrf->setNeighbors(16, 17, 1);
  mrf->setNeighbors(17, 18, 1);
  mrf->setNeighbors(18, 19, 1);

  mrf->setNeighbors(0, 5, 1);
  mrf->setNeighbors(1, 6, 1);
  mrf->setNeighbors(2, 7, 1);
  mrf->setNeighbors(3, 8, 1);
  mrf->setNeighbors(4, 9, 1);
  mrf->setNeighbors(5, 10, 1);
  mrf->setNeighbors(6, 11, 1);
  mrf->setNeighbors(7, 12, 1);
  mrf->setNeighbors(8, 13, 1);
  mrf->setNeighbors(9, 14, 1);
  mrf->setNeighbors(10, 15, 1);
  mrf->setNeighbors(11, 16, 1);
  mrf->setNeighbors(12, 17, 1);
  mrf->setNeighbors(13, 18, 1);
  mrf->setNeighbors(14, 19, 1);


  MRF::EnergyVal E;
  float t,tot_t;

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
  for (int i=0; i<6; i++)
  {
    unsigned char rowvalue = mrf->getLabel(i);
    printf(" Node : %d  state : %d \n", i, (int)rowvalue);
  }
  */


  delete energy;
//  delete mrf;

  return 0;
}


