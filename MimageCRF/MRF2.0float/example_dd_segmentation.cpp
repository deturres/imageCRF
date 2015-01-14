/*
 * example_dd_segmentation.cpp
 *
 *  Created on: May 22, 2012
 *      Author: bhole
 */

// this program assumes you pass the folder name as a parameter
// and Euni.txt as the local cost (2 colums), Ehor.txt as hotizontal and Ever as
// vertical energies (4 columns)

#include "mrf.h"
#include "Flipper.h"
#include "DualDecompositionFlipper.h"
#include "DualDecompositionSubGraph.h"
#include "DualDecompositionTree.h"
#include "TRW-S.h"
#include "MaxProdBP.h"
#include "SyncSumProd.h"
#include "BP-S.h"
#include "QPBO-v1.3.src/QPBO.cpp"

#include "imageLib.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <new>
#include <fstream>

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

VVVF Euni;
VVVVF Ehor;
VVVVF Ever;

MRF::CostVal *D; //[numPixels*numLabels];



std::vector<int> *flippedinfo;
std::vector<int> *flippedinfo_f;
std::vector<int> *subgraph_node_to_graph_node_g;


const double primaldualtol = 1e-2;




// we assume the pixels will be from a grid
MRF::CostVal fnCost(int pix1, int pix2, int labeli, int labelj)
{

  MRF::CostVal answer = 0.0;

  int px1 = pix1/sizeX;
  int py1 = pix1%sizeX;

  int px2 = pix2/sizeX;
  int py2 = pix2%sizeX;

  if (abs(px1-px2)==1 && py1==py2)
  {
    // use vertical energies as x is along height
    if (px1 < px2)
      answer = Ever[px1][py1][labeli][labelj];  // labeli is for pix1
    else if (px2 < px1)
      answer = Ever[px2][py2][labelj][labeli];  // labelj is for pix2

  } else if (abs(py1-py2)==1 && px1==px2)
  {
    // use horizontal energies as y is along width
    if (py1 < py2)
      answer = Ehor[px1][py1][labeli][labelj];  // labeli is for pix1
    else if (py2 < py1)
      answer = Ehor[px2][py2][labelj][labeli];  // labelj is for pix2
  } else
  {
    printf("%d(%d, %d) %d(%d,%d)\n",pix1, px1, py1, pix2, px2, py2);
    printf("Error not possible \n");
    exit(1);
  }


  return -answer;
}


//////////////////////////////////////////////////////////////
/////////////// METHODS NEEDED BY DDS ////////////////////////


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

//////////////////////END OF METHODS NEEDED BY DDS ///////////////////////










//////////////////////////////////////////////////////////////
////////// METHODS NEEDED BY DDF /////////////////////////////
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

//////////////////////END OF METHODS NEEDED BY DDF ///////////////////////



////////////// FOR PURE FLIPPER GRAPH //////////////////////////////

MRF::CostVal fnCost_pure_flipper(int pix1, int pix2, int i, int j)
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


//////////////// END OF PURE FLIPPER GRAPH //////////////////////////









DataCost* computeDataCost()
{
  int dsiIndex = 0;

  int row=0, col=0;

  for (int p = 0; p < numPixels; p++)
  {
    for (int m = 0;  m < numLabels; m++)
    {
      D[dsiIndex++] =  -Euni[row][col][m];
    }
    col++;
    if (col>=sizeX)
    {
      col=0;
      row++;
    }
  }

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


typedef double REAL;


int writeToImage(QPBO<REAL> *q, CByteImage *imo)
{
  int n = 0;

  for (int y = 0; y < sizeY; y++)
  {
    unsigned char *row = &imo->Pixel(0, y, 0);

    for (int x = 0; x < sizeX; x++)
    {
      // do i need to type cast to char explicitly ?
      row[x] = q->GetLabel(n++);
    }
  }

  return 0;
}




void readEunifile(const char* foldername)
{
  char Eunif[500];
  sprintf(Eunif, "%s/Euni.txt", foldername);

  int nolines=0;
  ifstream myfile (Eunif);

  float E0, E1;
  int row=0, col=0;

  if (myfile.is_open())
  {
    while ( myfile.good() )
    {
      myfile >> E0;
      if (myfile.eof())
        break;
      myfile >> E1;
      // printf("%f %f\n", E0, E1);
      Euni[row][col][0] = E0;
      Euni[row][col][1] = E1;
      nolines++;
      col++;
      if (col>=sizeX)
      {
        col=0;
        row++;
      }

    }

    if (nolines != sizeX*sizeY)
      printf("All lines not read Euni %d %d\n", nolines, sizeX*sizeY);

    myfile.close();
  }
  else printf("Unable to open Euni file\n");
}


void readEhorfile(const char* foldername)
{
  char Ehorf[500];
  sprintf(Ehorf, "%s/Ehor.txt", foldername);

  int nolines=0;
  ifstream myfile (Ehorf);

  float E00, E01, E10, E11;
  int row=0, col=0;

  if (myfile.is_open())
  {
    while ( myfile.good() )
    {
      myfile >> E00;
      if (myfile.eof())
        break;
      myfile >> E01;
      myfile >> E10;
      myfile >> E11;

      // printf("%f %f\n", E0, E1);
      Ehor[row][col][0][0] = E00;
      Ehor[row][col][0][1] = E01;
      Ehor[row][col][1][0] = E10;
      Ehor[row][col][1][1] = E11;
      nolines++;
      col++;
      if (col>=sizeX-1)
      {
        col=0;
        row++;
      }

    }

    if (nolines != (sizeX-1)*sizeY)
      printf("All lines not read Ehor %d %d\n", nolines, sizeX*sizeY);

    myfile.close();
  }
  else printf("Unable to open Ehor file\n");
}




void readEverfile(const char* foldername)
{
  char Everf[500];
  sprintf(Everf, "%s/Ever.txt", foldername);

  int nolines=0;
  ifstream myfile (Everf);

  float E00, E01, E10, E11;
  int row=0, col=0;

  if (myfile.is_open())
  {
    while ( myfile.good() )
    {
      myfile >> E00;
      if (myfile.eof())
        break;
      myfile >> E01;
      myfile >> E10;
      myfile >> E11;

      // printf("%f %f\n", E0, E1);
      Ever[row][col][0][0] = E00;
      Ever[row][col][0][1] = E01;
      Ever[row][col][1][0] = E10;
      Ever[row][col][1][1] = E11;
      nolines++;
      col++;
      if (col>=sizeX)
      {
        col=0;
        row++;
      }

    }

    if (nolines != sizeX*(sizeY-1))
      printf("All lines not read Ever %d %d\n", nolines, sizeX*sizeY);

    myfile.close();
  }
  else printf("Unable to open Ever file\n");
}




void readEnergyTextFiles(const char* foldername)
{

  readEunifile(foldername);

  readEhorfile(foldername);

  readEverfile(foldername);

}

// Y is for height, X is for width
void allocateEnergyVectors()
{
  Euni.resize(sizeY);
  Ehor.resize(sizeY);
  Ever.resize(sizeY-1); // 1 row less

  for (int i=0; i<sizeY; i++)
  {
    Euni[i].resize(sizeX);
    Ehor[i].resize(sizeX-1); // 1 column less
  }

  for (int i=0; i<sizeY-1; i++)
  {
    Ever[i].resize(sizeX);
  }


  for (int i=0; i<sizeY; i++)
    for (int j=0; j<sizeX; j++)
      Euni[i][j].resize(2);

  for (int i=0; i<sizeY; i++)
    for (int j=0; j<sizeX-1; j++)
    {
      Ehor[i][j].resize(2);
      for (int k=0; k<2; k++)
      {
        Ehor[i][j][k].resize(2);
      }
    }

  for (int i=0; i<sizeY-1; i++)
    for (int j=0; j<sizeX; j++)
    {
      Ever[i][j].resize(2);
      for (int k=0; k<2; k++)
      {
        Ever[i][j][k].resize(2);
      }
    }
}





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



int main(int argc, char* argv[])
{
  //srand(time(NULL));
  srand(0);
  sto = new StochasticLib1((int) 0) ; //(int) time(NULL));

  if (argc != 3)
  {
    printf(" Please supply input folder and output location (with image prefix) on commandline ! \n");
    exit(1);
  }

  char imfilename[500];
  sprintf(imfilename, "%s/im.png", argv[1]);

  string ofname(argv[2]);

  char outnameW[500];

  CByteImage im;
  // read input image
  ReadImage(im, imfilename);

  sizeX = im.Shape().width;
  sizeY = im.Shape().height;

  numPixels = sizeX * sizeY;

  D = new MRF::CostVal[numPixels*numLabels];

  CByteImage imo;
  CShape sh(sizeX, sizeY, 1);
  imo.ReAllocate(sh);

  printf("%s _ w = %d _ h = %d _ w*h = %d \n", imfilename, sizeX, sizeY, numPixels);


  // allocate energy vectors
  // printf("Allocating memory for vectors\n");
  allocateEnergyVectors();
  // read energy text files
  // energies get stored in Euni, Ehor, Ever (which are global vectors)
  // printf("Reading energy into vectors\n");
  readEnergyTextFiles(argv[1]);

  // printf("Computing datacost\n");
  DataCost *dcost = computeDataCost();


  SmoothnessCost *scost = new SmoothnessCost((MRF::SmoothCostGeneralFn) fnCost);
  EnergyFunction *energy = new EnergyFunction(dcost, scost);

  bool runF = false;
  bool runDDF = false;
  bool runDDG = false;
  bool runDD = false;
  bool runTRWS = true;
  bool runMaxProdBP = false;
  bool runBPS       = false;
  bool runSSP       = false;
  bool runQPBO  = true;


  if (runF)
  {

    MRF* mrf = new Flipper(numPixels, numLabels, energy);

    // you need to set submodular/nonsubmodular function
    ((Flipper*)mrf)->setFlipperFunction((MRF::SmoothCostGeneralFn) fnCost_pure_flipper);
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


    writeToImage(mrf, &imo);
    sprintf(outnameW, "%s_ddf.png", ofname.c_str());
    CByteImage disp2;
    ScaleAndOffset(imo, disp2, 255.0, 0);
    WriteImage(disp2, outnameW);



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
    ((DualDecompositionSubGraph*)mrf)->setToleranceParam(primaldualtol);


    E = mrf->totalEnergy();

    printf("Energy at the Start= %g (%g,%g)\n", (float)E,
      (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());



    tot_t = 0;
    mrf->optimize(1, t);
    E = mrf->totalEnergy();
    tot_t = tot_t + t ;
    printf("energy = %g (%f secs)\n", (float)E, tot_t);


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




  ////////////////////////////////////////////////
  //                  QPBO                      //
  ////////////////////////////////////////////////

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
    q->ComputeWeakPersistencies();

    // stop timer
    clock_t finish = clock();
    float timed = (float) (((double)(finish - start)) / CLOCKS_PER_SEC);


    printf("energy = %g, lower bound = %f (%f secs)\n",
        (float)q->ComputeTwiceEnergy()/2.0, q->ComputeTwiceLowerBound()/2.0, timed);

/*
    for (int i=0; i<numPixels; i++)
    {
      int rowvalue = q->GetLabel(i);
      printf(" Node : %d  state : %d \n", i, (int)rowvalue);
    }
*/


    writeToImage(q, &imo);
    sprintf(outnameW, "%s_qpbo.png", ofname.c_str());
    CByteImage disp2;
    ScaleAndOffset(imo, disp2, 255.0, 0);
    WriteImage(disp2, outnameW);

  }





  delete energy;
  delete D;

  return 0;
}




