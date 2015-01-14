/*
 * example_dd_realworld2.cpp
 *
 *  Created on: May 9, 2012
 *      Author: bhole
 */



// this program is used to denoise examples

// mix of submodular and  non submodular edges on real world images

#include "mrf.h"
#include "DualDecompositionSubGraph.h"
#include "DualDecompositionTree.h"
#include "TRW-S.h"
#include "MaxProdBP.h"
#include "SyncSumProd.h"
#include "BP-S.h"

#include "imageLib.h"

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


int sizeX ;  // width
int sizeY ;  // height
int numPixels ;


const int outerIter = 100;
const int innerIter = 1;

const int numLabels = 2;

MRF::CostVal *D; //[numPixels*numLabels];

std::vector<int> *flippedinfo;
std::vector<int> *subgraph_node_to_graph_node_g;


double a = 1;
double b = -1.1;
double c = 0.5;
double d = -1;

CByteImage im;

const double primaldualtol = 1e-2;

// images are grayscale
int getIntensity(int pix)
{
  // extract x,y coordinates from pix
  int yrow = pix/sizeX;
  int xcol = pix % sizeX;

  unsigned char *pix1t = 0;
  pix1t = &im.Pixel(xcol, yrow, 0);

  return ((int)pix1t[0]);
}

MRF::CostVal fnCost(int pix1, int pix2, int labeli, int labelj)
{

  MRF::CostVal answer = 0.0;

  if (labeli == labelj)
  {}
  else
  {
    double p1val = getIntensity(pix1)/255.0;
    double p2val = getIntensity(pix2)/255.0;

    answer = a + b*fabs(p1val - p2val);
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






DataCost* computeDataCost()
{
  int dsiIndex = 0;

  for (int p = 0; p < numPixels; p++)
  {
    for (int m = 0;  m < numLabels; m++)
    {
      double r = 0;
      if (m == 1)
      {
        int intensity = getIntensity(p);
        r = c + d*(intensity/255.0);
      }

      MRF::CostVal dsiValue = r;
      D[dsiIndex++] =  (float)dsiValue;
    }
  }

  // handling interactive points by maxing local costs extreme so they would favor
  // one label over the other
/*
  // let (0,0) label be 1
  // let (99,99) label be 0
  D[0] = -100.0;
  D[1] = +100.0;

//  D[numPixels*numLabels - 2] = -100.0;
//  D[numPixels*numLabels - 1] = +100.0;
*/

  // let (0,99) == pix 99 ie [198==-ve 199==+ve]   label be 0
  // let (99,0) == pix 9900 ie [19800==+ve 19801==-ve] label be 1
//  D[198] = +100.0;
//  D[199] = -100.0;

//  D[19800] = -100.0;
//  D[19801] = +100.0;

// no interaction needed for these examples

  return new DataCost (D);
}


int writeToImage(MRF *mrf, CByteImage *imo)
{
  int n = 0;

  for (int y = 0; y < sizeY; y++)
  {
    unsigned char *row = &imo->Pixel(0, y, 0);

    for (int x = 0; x < sizeX; x++)
    {
      row[x] = mrf->getLabel(n++);
    }
  }

  return 0;
}



int main(int argc, char* argv[])
{
  //srand(time(NULL));
  srand(0);
  sto = new StochasticLib1((int) 0) ; //(int) time(NULL));

  if (argc != 3 && argc != 7)
  {
    printf(" Please supply input image and output location (with image prefix) on commandline ! \n");
    printf(" OR Please supply input image and output location (with image prefix), a, b, c and d on commandline ! \n");
    printf(" Note that trws, ddg, ddt, bps, maxp will be attached to the output name \n");
    exit(1);
  }

  string fname(argv[1]);
  string ofname(argv[2]);

  if (argc == 7)
  {
    string astr(argv[3]);
    string bstr(argv[4]);
    string cstr(argv[5]);
    string dstr(argv[6]);
    // a = std::stod(astr);
    // b = std::stod(bstr);
    // c = std::stod(cstr);
    // d = std::stod(dstr);
    a = atof(astr.c_str());
    b = atof(bstr.c_str());
    c = atof(cstr.c_str());
    d = atof(dstr.c_str());
  }

  char outnameW[500];

  // read input image
  ReadImage(im, fname.c_str());

  sizeX = im.Shape().width;
  sizeY = im.Shape().height;

  numPixels = sizeX * sizeY;

  D = new MRF::CostVal[numPixels*numLabels];

  CByteImage imo;
  CShape sh(sizeX, sizeY, 1);
  imo.ReAllocate(sh);

  printf("%s _ w = %d _ h = %d _ w*h = %d \n", fname.c_str(), sizeX, sizeY, numPixels);

  DataCost *dcost = computeDataCost();
  SmoothnessCost *scost = new SmoothnessCost((MRF::SmoothCostGeneralFn) fnCost);
  EnergyFunction *energy = new EnergyFunction(dcost, scost);

  bool runDDG = true;
  bool runDD = true;
  bool runTRWS = true;
  bool runMaxProdBP = false;
  bool runBPS       = true;
  bool runSSP       = false;

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
    ((DualDecompositionSubGraph*)mrf)->setToleranceParam(primaldualtol);

    // initialize labels to previous disps
/*    int n = 0;

    for (int y = 0; y < sizeY; y++)
    {
      for (int x = 0; x < sizeX; x++)
      {
        int labelv = 0;
        if (rand() % 2 == 1)
          labelv = 1;
        printf("%d ", labelv);
        mrf->setLabel(n++, labelv);
      }
    }
*/




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

    writeToImage(mrf, &imo);
    sprintf(outnameW, "%s_ddg.png", ofname.c_str());
    CByteImage disp2;
    ScaleAndOffset(imo, disp2, 255.0, 0);
    WriteImage(disp2, outnameW);


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
    ((DualDecompositionTree*)mrf)->setToleranceParam(primaldualtol);

    E = mrf->totalEnergy();
    printf("Energy at the Start= %g (%g,%g)\n", (float)E,
        (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

    tot_t = 0;

    mrf->optimize(innerIter, t);

    E = mrf->totalEnergy();
    tot_t = tot_t + t ;
    printf("energy = %g (%f secs)\n", (float)E, tot_t);

    writeToImage(mrf, &imo);
    sprintf(outnameW, "%s_ddt.png", ofname.c_str());
    CByteImage disp2;
    ScaleAndOffset(imo, disp2, 255.0, 0);
    WriteImage(disp2, outnameW);

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
    for (int iter=0; iter<200; iter++)
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
    writeToImage(mrf, &imo);
    sprintf(outnameW, "%s_trws.png", ofname.c_str());
    CByteImage disp2;
    ScaleAndOffset(imo, disp2, 255.0, 0);
    WriteImage(disp2, outnameW);


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
    for (int iter=0; iter < 200; iter++)
    {
      mrf->optimize(1, t);

      E = mrf->totalEnergy();
      tot_t = tot_t + t ;
      printf("energy = %g (%f secs)\n", (float)E, tot_t);
    }

    writeToImage(mrf, &imo);
    sprintf(outnameW, "%s_map.png", ofname.c_str());
    CByteImage disp2;
    ScaleAndOffset(imo, disp2, 255.0, 0);
    WriteImage(disp2, outnameW);

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

    writeToImage(mrf, &imo);
    sprintf(outnameW, "%s_bps.png", ofname.c_str());
    CByteImage disp2;
    ScaleAndOffset(imo, disp2, 255.0, 0);
    WriteImage(disp2, outnameW);


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
  delete D;

  return 0;
}



