/*
 * example_flipper2.cpp
 *
 *  Created on: May 7, 2012
 *      Author: bhole
 */

#include "mrf.h"
#include "Flipper.h"
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

const int sizeX = 1000;   // width
const int sizeY = 1000; // height
const int numPixels = 1000000; // conversion is pix = h*sizeX + w  (h is along height, w along width

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

myumapvf *mappx;



// For Flipper methods
std::vector<int> *flippedinfo_f;

MRF::CostVal fnCost(int pix1, int pix2, int i, int j)
{
  MRF::CostVal answer = 0.0;

  int smallerpix = pix1;
  int largerpix = pix2;
  int labeli = i;
  int labelj = j;

  if (pix1 > pix2)
  {
    smallerpix = pix2;
    largerpix = pix1;
    labeli = j;
    labelj = i;
  }
  char buff[30];
  sprintf(buff, "%012d%012d", smallerpix, largerpix);
  std::string str(buff);

  myumapvf::const_iterator got = (*mappx).find (buff);
  Epair ep;
  if ( got == (*mappx).end() )
  {
    printf("\n Real problem here : smallerpix %d largerpix %d  buff %s\n", smallerpix, largerpix, buff);
    exit(1);
  }
  else
    ep = got->second;

  if (labeli == 0 && labelj == 0)
    answer = ep[0];
  else if (labeli == 0 && labelj == 1)
    answer = ep[1];
  else if (labeli == 1 && labelj == 0)
    answer = ep[2];
  else if (labeli == 1 && labelj == 1)
    answer = ep[3];

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




// assuming the size is 1000x1000
void setSubNonSubPairwiseEnergies()
{
  mappx = new(myumapvf);

  int pix1 = 0;
  int pix2 = 0;

  // Y = height, X = width
  // (1) have values for 0*0 to Y*299 as submodular
  // (2) have values for Y*650 to Y*999 as submodular
  // (3) have values for Y*300 to Y*649 as nonsubmodular
  // (4) have alternate edges Y*299 to Y*300 as submod and nonsubmod
  // (5) have alternate edges Y*649 to Y*650 as nonsubmod and submod

  double r = 0.0;

  // Epair is flattened as e(0,0), e(0,1), e(1,0), e(1,1)

  // (1) //////////////////////////////////////////////////////////////

  for (int j=0; j<sizeY; j++)
  {
    for (int i=0; i<299; i++)
    {
      pix1 = j*sizeX+i;
      pix2 = j*sizeX+i+1;
      // note that pix1 is smaller than pix2
      char buff[30];
      sprintf(buff, "%012d%012d", pix1, pix2);
      std::string str(buff);

      // submodular
      Epair ep(4, 0.0); // default ep(0) == ep(3) == 0
      r = fabs(sto->Normal(0, 3));
      // ep(0,1) = r;
      ep[1] = r;
      r = fabs(sto->Normal(0, 3));
      //ep(1,0) = r;
      ep[2] = r;

      std::pair<std::string, Epair> pairsi(buff, ep);
      (*mappx).insert(pairsi);
    }
  }

  for (int j=0; j<sizeY-1; j++)
  {
    for (int i=0; i<300; i++)
    {
      pix1 = j*sizeX+i;
      pix2 = (j+1)*sizeX+i;
      // note that pix1 is smaller than pix2
      char buff[30];
      sprintf(buff, "%012d%012d", pix1, pix2);
      std::string str(buff);

      // submodular
      Epair ep(4, 0.0); // default ep(0) == ep(3) == 0
      r = fabs(sto->Normal(0, 3));
      // ep(0,1) = r;
      ep[1] = r;
      r = fabs(sto->Normal(0, 3));
      //ep(1,0) = r;
      ep[2] = r;

      std::pair<std::string, Epair> pairsi(buff, ep);
      (*mappx).insert(pairsi);
    }
  }


  // (2) //////////////////////////////////////////////////////////////

  for (int j=0; j<sizeY; j++)
  {
    for (int i=650; i<999; i++)
    {
      pix1 = j*sizeX+i;
      pix2 = j*sizeX+i+1;
      // note that pix1 is smaller than pix2
      char buff[30];
      sprintf(buff, "%012d%012d", pix1, pix2);
      std::string str(buff);

      // submodular
      Epair ep(4, 0.0); // default ep(0) == ep(3) == 0
      r = fabs(sto->Normal(0, 3));
      // ep(0,1) = r;
      ep[1] = r;
      r = fabs(sto->Normal(0, 3));
      //ep(1,0) = r;
      ep[2] = r;

      std::pair<std::string, Epair> pairsi(buff, ep);
      (*mappx).insert(pairsi);
    }
  }

  for (int j=0; j<sizeY-1; j++)
  {
    for (int i=650; i<1000; i++)
    {
      pix1 = j*sizeX+i;
      pix2 = (j+1)*sizeX+i;
      // note that pix1 is smaller than pix2
      char buff[30];
      sprintf(buff, "%012d%012d", pix1, pix2);
      std::string str(buff);

      // submodular
      Epair ep(4, 0.0); // default ep(0) == ep(3) == 0
      r = fabs(sto->Normal(0, 3));
      // ep(0,1) = r;
      ep[1] = r;
      r = fabs(sto->Normal(0, 3));
      //ep(1,0) = r;
      ep[2] = r;

      std::pair<std::string, Epair> pairsi(buff, ep);
      (*mappx).insert(pairsi);
    }
  }


  // (3) //////////////////////////////////////////////////////////////

  for (int j=0; j<sizeY; j++)
  {
    for (int i=300; i<649; i++)
    {
      pix1 = j*sizeX+i;
      pix2 = j*sizeX+i+1;
      // note that pix1 is smaller than pix2
      char buff[30];
      sprintf(buff, "%012d%012d", pix1, pix2);
      std::string str(buff);

      // nonsubmodular
      Epair ep(4, 0.0); // default ep(1) == ep(2) == 0
      r = fabs(sto->Normal(0, 3));
      // ep(0,0) = r;
      ep[0] = r;
      r = fabs(sto->Normal(0, 3));
      //ep(1,1) = r;
      ep[3] = r;

      std::pair<std::string, Epair> pairsi(buff, ep);
      (*mappx).insert(pairsi);
    }
  }

  for (int j=0; j<sizeY-1; j++)
  {
    for (int i=300; i<650; i++)
    {
      pix1 = j*sizeX+i;
      pix2 = (j+1)*sizeX+i;
      // note that pix1 is smaller than pix2
      char buff[30];
      sprintf(buff, "%012d%012d", pix1, pix2);
      std::string str(buff);

      // nonsubmodular
      Epair ep(4, 0.0); // default ep(1) == ep(2) == 0
      r = fabs(sto->Normal(0, 3));
      // ep(0,0) = r;
      ep[0] = r;
      r = fabs(sto->Normal(0, 3));
      //ep(1,1) = r;
      ep[3] = r;

      std::pair<std::string, Epair> pairsi(buff, ep);
      (*mappx).insert(pairsi);
    }
  }


  // (4) //////////////////////////////////////////////////////////////
  for (int j=0; j<sizeY; j++)
  {
    int i = 299;
    pix1 = j*sizeX+i;
    pix2 = j*sizeX+i+1;
    // note that pix1 is smaller than pix2
    char buff[30];
    sprintf(buff, "%012d%012d", pix1, pix2);
    std::string str(buff);

    if (j%2 == 0)
    {
      // submodular
      Epair ep(4, 0.0); // default ep(0) == ep(3) == 0
      r = fabs(sto->Normal(0, 3));
      // ep(0,1) = r;
      ep[1] = r;
      r = fabs(sto->Normal(0, 3));
      //ep(1,0) = r;
      ep[2] = r;

      std::pair<std::string, Epair> pairsi(buff, ep);
      (*mappx).insert(pairsi);
    } else
    {
      // nonsubmodular
      Epair ep(4, 0.0); // default ep(1) == ep(2) == 0
      r = fabs(sto->Normal(0, 3));
      // ep(0,0) = r;
      ep[0] = r;
      r = fabs(sto->Normal(0, 3));
      //ep(1,1) = r;
      ep[3] = r;

      std::pair<std::string, Epair> pairsi(buff, ep);
      (*mappx).insert(pairsi);
    }
  }



  // (5) //////////////////////////////////////////////////////////////
  for (int j=0; j<sizeY; j++)
  {
    int i = 649;
    pix1 = j*sizeX+i;
    pix2 = j*sizeX+i+1;
    // note that pix1 is smaller than pix2
    char buff[30];
    sprintf(buff, "%012d%012d", pix1, pix2);
    std::string str(buff);

    if (j%2 == 1)
    {
      // submodular
      Epair ep(4, 0.0); // default ep(0) == ep(3) == 0
      r = fabs(sto->Normal(0, 3));
      // ep(0,1) = r;
      ep[1] = r;
      r = fabs(sto->Normal(0, 3));
      //ep(1,0) = r;
      ep[2] = r;

      std::pair<std::string, Epair> pairsi(buff, ep);
      (*mappx).insert(pairsi);
    } else
    {
      // nonsubmodular
      Epair ep(4, 0.0); // default ep(1) == ep(2) == 0
      r = fabs(sto->Normal(0, 3));
      // ep(0,0) = r;
      ep[0] = r;
      r = fabs(sto->Normal(0, 3));
      //ep(1,1) = r;
      ep[3] = r;

      std::pair<std::string, Epair> pairsi(buff, ep);
      (*mappx).insert(pairsi);
    }
  }
}




int main()
{
  srand(0);
  sto = new StochasticLib1((int) 0); // time(NULL) ==> returns secs since 1970

  setSubNonSubPairwiseEnergies();

  DataCost *dcost = computeDataCost();
  SmoothnessCost *scost = new SmoothnessCost((MRF::SmoothCostGeneralFn) fnCost);
  EnergyFunction *energy = new EnergyFunction(dcost, scost);

  bool runF = true;
  bool runDDG = false;
  bool runDD = false;
  bool runTRWS = true;
  bool runMaxProdBP = false;
  bool runBPS       = false;
  bool runSSP       = false;

  bool runDDTes = true;

  if (runF)
  {

    MRF* mrf = new Flipper(numPixels, numLabels, energy);

    // you need to set submodular/nonsubmodular function
    ((Flipper*)mrf)->setFlipperFunction((MRF::SmoothCostGeneralFn) fnCost_flipper);
    // you need to set *flippedinfo before calling that function
    ((Flipper*)mrf)->setflippedinfo(&flippedinfo_f);


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

    /*
     for (int i=0; i<numPixels; i++)
     {
     unsigned char rowvalue = mrf->getLabel(i);
     printf(" Node : %d  state : %d \n", i, (int)rowvalue);
     }
     */

    delete mrf;

  }


  if (runDDG)
  {

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


  if (runDD)
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
  // DD Tree ES
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
    for (int iter=0; iter<400; iter++)
    {
      mrf->optimize(10, t);

      E = mrf->totalEnergy();
      lowerBound = mrf->lowerBound();
      tot_t = tot_t + t ;
      printf("energy = %g, lower bound = %f (%f secs)\n", (float)E, lowerBound, tot_t);

      if (fabs((double)E - lowerBound) < primaldualtol)
        break;

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
    for (int iter=0; iter<100; iter++)
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

  delete energy;
  delete scost;
  delete dcost;

  return 0;
}
