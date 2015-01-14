/*
 * example_tree.cpp
 *
 *  Created on: Oct 3, 2011
 *      Author: bhole
 */

// sample to run MaxProdBPTree

#include "mrf.h"
#include "MaxProdBPTree.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <new>


#define numPixels 6
#define numEdges 5
#define numLabels 3

MRF::CostVal D[numPixels*numLabels];

MRF::CostVal fnCost(int pix1, int pix2, int i, int j)
{

  /*
  double fncostarray[numEdges][numLabels][numLabels] = {{{1, 1.3, 1.4}, {1.3, 1, 1.7}, {1.4, 1.7, 1}},
		  	  	  	  	  	  	  	  	  	  {{1, 2.3, 2.4}, {2.3, 1, 2.7}, {2.4, 2.7, 1}},
		  	  	  	  	  	  	  	  	  	  {{1, 3.3, 3.4}, {3.3, 1, 3.7}, {3.4, 3.7, 1}},
		  	  	  	  	  	  	  	  	  	  {{1, 4.3, 4.4}, {4.3, 1, 4.7}, {4.4, 4.7, 1}},
		  	  	  	  	  	  	  	  	  	  {{1, 5.3, 5.4}, {5.3, 1, 5.7}, {5.4, 5.7, 1}}
  	  	  	  	  	  	  	  	  	  	  	  };
  */
  double fncostarray[numEdges][numLabels][numLabels] = {{{15, 20, 10}, {20, 15, 20}, {10, 20, 15}},
                                          {{15, 20, 20}, {20, 15, 10}, {20, 10, 15}},
                                          {{15, 20, 20}, {20, 15, 20}, {20, 20, 10}},
                                          {{15, 20, 20}, {20, 15, 10}, {20, 10, 15}},
                                          {{15, 20, 10}, {20, 15, 20}, {10, 20, 15}}
                                              };


  int arrayidx = 99;
  if ((pix1==0 && pix2==1) || (pix1==1 && pix2==0))
  {
    arrayidx = 0;
  }
  if ((pix1==2 && pix2==1) || (pix1==1 && pix2==2))
  {
    arrayidx = 1;
  }
  if ((pix1==4 && pix2==1) || (pix1==1 && pix2==4))
  {
    arrayidx = 2;
  }
  if ((pix1==3 && pix2==4) || (pix1==4 && pix2==3))
  {
    arrayidx = 3;
  }
  if ((pix1==4 && pix2==5) || (pix1==5 && pix2==4))
  {
    arrayidx = 4;
  }

  MRF::CostVal answer = fncostarray[arrayidx][i][j];
  return answer;
}


DataCost* computeDataCost()
{
//  double datac[numPixels][numLabels] = {{1,2,3}, {5,4,1}, {4,1,1}, {1,5,5}, {1,4,1}, {2,3,2}};
  double datac[numPixels][numLabels] = {{10,15,20}, {20,15,10}, {15,10,20}, {20,15,10}, {10,15,20}, {10, 15, 20}};


  int dsiIndex = 0;

  for (int p = 0; p < numPixels; p++)
  {
    for (int d = 0; d < numLabels; d++)
    {
      // MRF::CostVal dsiValue = (((p-d)*(p-d)) % 4)*(float(p+d) / 3) * 0.5 + 1.5;
      MRF::CostVal dsiValue = datac[p][d];
      D[dsiIndex++] =  (float)dsiValue; //((MRF::CostVal)(rand() % 100))/10 + 1;
      printf(" p = %d , d = %d, dsiValue = %f \n", p, d, D[dsiIndex-1]);
    }
  }
  return new DataCost (D);
}



int main()
{
  DataCost *dcost = computeDataCost();
  SmoothnessCost *scost = new SmoothnessCost((MRF::SmoothCostGeneralFn) fnCost);
  EnergyFunction *energy = new EnergyFunction(dcost, scost);
  MRF* mrf = new MaxProdBPTree(numPixels, numLabels, energy);

  mrf->setNeighbors(0, 1, 1);
  mrf->setNeighbors(1, 2, 1);
  mrf->setNeighbors(3, 4, 1);
  mrf->setNeighbors(4, 5, 1);
  mrf->setNeighbors(1, 4, 1);

  MRF::EnergyVal E;
  float t,tot_t;

  printf("\n*******  Started MaxProdTree Belief Propagation *****\n");
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

  for (int i=0; i<numPixels; i++)
  {
    unsigned char rowvalue = mrf->getLabel(i);
    printf(" Node : %d  state : %d \n", i, (int)rowvalue);
  }


  delete energy;
  delete mrf;

  return 0;
}





