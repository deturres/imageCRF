/*
 * example_dd_flipper.cpp
 *
 *  Created on: May 6, 2012
 *      Author: bhole
 */

#include "mrf.h"
#include "DualDecompositionFlipper.h"
#include "TRW-S.h"
#include "Flipper.h"

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

const int sizeX = 100;
const int sizeY = 100;
const int numPixels = 10000;

const int outerIter = 100;
const int innerIter = 1;

const int numLabels = 2;

MRF::CostVal D[numPixels*numLabels];

float fncostarray[numLabels][numLabels];
float nonsubfncostarray[numLabels][numLabels];

// For Flipper methods
// this will dynamically point to correct flipper graph when Flipper -> optimize() is called
std::vector<int> *flippedinfo_f;
std::vector<int> *subgraph_node_to_graph_node_g;

MRF::CostVal fnCost(int pix1, int pix2, int i, int j)
{

  MRF::CostVal answer = 0.0;


  int smallerpix = pix1;
  int largerpix = pix2;
  if (pix1 > pix2)
  {
    smallerpix = pix2;
    largerpix = pix1;
  }


  if (smallerpix % sizeX <= (sizeX/2))
    answer = fncostarray[i][j];
  else
    answer = nonsubfncostarray[i][j];
/*

  if (smallerpix == 0 && largerpix == 1)
    answer = nonsubfncostarray[i][j];
  else if (smallerpix == 0 && largerpix == 2)
    answer = nonsubfncostarray[i][j];
  else if (smallerpix == 1 && largerpix == 3)
    answer = nonsubfncostarray[i][j];
  else if (smallerpix == 2 && largerpix == 3)
    answer = fncostarray[i][j];
*/

  return answer;
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


MRF::CostVal fnCost_only_flipper(int pix1, int pix2, int i, int j)
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




void computePairwiseMatrix()
{

  for (int ii=0; ii<numLabels; ii++)
    for (int jj=ii; jj<numLabels; jj++)
    {
      // value was 3.7 fixed before
      double r = fabs(sto->Normal(0, 3)); // is always +ve and maintains all values as submodular
      fncostarray[ii][jj] = fncostarray[jj][ii] = (ii == jj) ? (MRF::CostVal)0 : (MRF::CostVal)r;
    }

  for (int ii=0; ii<numLabels; ii++)
    for (int jj=ii; jj<numLabels; jj++)
    {
      // value was 2.3 fixed before
      double r = fabs(sto->Normal(0, 3)); // is always +ve and maintains all values as non-submodular
      nonsubfncostarray[ii][jj] = nonsubfncostarray[jj][ii] = (ii == jj) ? (MRF::CostVal)r : (MRF::CostVal)0;
    }

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
      // if (p == 2 || p == 3 || p == 12 || p == 13 || p == 22 || p == 23)
      //  printf(" %d %f \n", p, r);
      D[dsiIndex++] =  (float)dsiValue; //((MRF::CostVal)(rand() % 100))/10 + 1;
    }
  }
  return new DataCost (D);
}




int main()
{

  srand(0);
  sto = new StochasticLib1((int) 0) ;

  computePairwiseMatrix();

  DataCost *dcost = computeDataCost();
  SmoothnessCost *scost = new SmoothnessCost((MRF::SmoothCostGeneralFn) fnCost);
  EnergyFunction *energy = new EnergyFunction(dcost, scost);

  bool runF = false;
  bool runDDF = true;
  bool runTRWS = true;


  if (runF)
  {

    MRF* mrf = new Flipper(numPixels, numLabels, energy);

    // you need to set submodular/nonsubmodular function
    ((Flipper*)mrf)->setFlipperFunction((MRF::SmoothCostGeneralFn) fnCost_only_flipper);
    // you need to set *flippedinfo before calling that function
    ((Flipper*)mrf)->setflippedinfo(&flippedinfo_f);


    for (int j=0; j<sizeY; j++)
    {
      for (int i=0; i<sizeX-1; i++)
      {
        if (j*sizeX+i == 2 || j*sizeX+i == 12 || j*sizeX+i == 22) {}
        else
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

/*
     for (int i=0; i<numPixels; i++)
     {
     unsigned char rowvalue = mrf->getLabel(i);
     printf(" Node : %d  state : %d \n", i, (int)rowvalue);
     }
*/

    delete mrf;

  }










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
    for (int iter=0; iter<200; iter++)
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



  return 0;
}
