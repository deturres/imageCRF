/*
 * example_dd_tree_es.cpp
 *
 *  Created on: May 27, 2012
 *      Author: bhole
 */

// edge sharing enabled.

#include "mrf.h"
#include "DualDecompositionTreeES.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <new>


const int numPixels = 20; //5x4
const int numLabels = 3;

MRF::CostVal D[numPixels*numLabels];

MRF::CostVal fnCost(int pix1, int pix2, int i, int j)
{

  double fncostarray[5][numLabels][numLabels] = {{{15, 20, 10}, {20, 15, 20}, {10, 20, 15}},
                                          {{15, 20, 20}, {20, 15, 10}, {20, 10, 15}},
                                          {{15, 20, 20}, {20, 15, 20}, {20, 20, 10}},
                                          {{15, 20, 20}, {20, 15, 10}, {20, 10, 15}},
                                          {{15, 20, 10}, {20, 15, 20}, {10, 20, 15}}
                                              };

  int arrayidx = (pix1+pix2) % 5;

  MRF::CostVal answer = fncostarray[arrayidx][i][j];
  return answer;
}


DataCost* computeDataCost()
{
  double datac[numPixels][numLabels] = {{10,15,20}, {20,15,10}, {15,10,20}, {20,15,10}, {10,15,20},
      {10,15,20}, {20,15,10}, {15,10,20}, {20,15,10}, {10,15,20},
      {10,15,20}, {20,15,10}, {15,10,20}, {20,15,10}, {10,15,20},
      {10,15,20}, {20,15,10}, {15,10,20}, {20,15,10}, {10,15,20}};
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
  EnergyFunction *energy = new EnergyFunction(dcost, scost);
  MRF* mrf = new DualDecompositionTreeES(numPixels, numLabels, energy);

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

  printf("\n*******  Started Dual Decomposition Tree *****\n");
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
  delete mrf;

  return 0;
}


