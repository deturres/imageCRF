/* crfmodel.cpp
 *  Heavily modified for image segmentation by Chetan Bhole
 *
 * $Id: crfstereo.cpp,v 1.10 2007/12/08 14:11:48 weinman Exp $
 * 11/13/2006 initial version by Daniel Scharstein
 *
 *
 */



// maximum number of iterations for graph cut / BP
// (can be big, since termination is controlled via closeEnoughPercent parameter)
#define MAXITER 1

#include "crfmodel.h"
#include "features.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include "GCoptimization.h"
#include "BP-S.h"
#include "MaxProdBP.h"
#include "SyncSumProd.h"
#include "SparseSyncSumProd.h"
#include "MeanField.h"
#include "SparseMeanField.h"
#include "TRW-S.h"
#include <vector>
#include <math.h>
#include <algorithm>
#include "Global.h"

extern GlobalPairwiseParam globalP;

using std::vector;

int verbose;
FILE *debugfile;

// these parameters are also defined in expectations.cpp
// make sure to change them there too
int loccubex = 4, loccubey = 4, loccubez = 4;

int skipXYZlocklr = 1;

// global variables to keep track of what needs to be freed later
MRF::CostVal *dsiArray = NULL;
std::vector<MRF::CostVal *> dsiArray2;
std::vector<int> globalNpixels;
int globalNdisps = 0;
bool randomv = false;

unsigned int cueIndex = 0;
std::vector<MRF::CostVal *> hCueArrayV;
std::vector<MRF::CostVal *> vCueArrayV;
std::vector<MRF::CostVal *> dCueArrayV;

std::vector<unsigned int *> hCueArrayVG;
std::vector<unsigned int *> vCueArrayVG;
std::vector<unsigned int *> dCueArrayVG;

std::vector<MRF::CostVal *> hBOcclArrayV;
std::vector<MRF::CostVal *> vBOcclArrayV;
std::vector<MRF::CostVal *> hSOcclArrayV;
std::vector<MRF::CostVal *> vSOcclArrayV;

std::vector<MRF::CostVal *> dsiArrayV;
std::vector<MRF::CostVal *> dsiArrayVGen;

std::vector<MRF::CostVal *> hpairwiseProbV;
std::vector<MRF::CostVal *> wpairwiseProbV;
std::vector<MRF::CostVal *> dpairwiseProbV;

std::vector<MRF::CostVal *> wpairwiseProbZ;




bool pairwise = false;
bool localOccl = false;
bool pairwiseInteraction = false;

std::vector<int> globalwidth; // for debugging
//float globallambda1 = 1;
//float globallambda2 = 1;

int inferencer = USE_BP;

MRF::CostVal occluded1 = 0.0;
MRF::CostVal occluded2 = 0.0;

void setGlobalNpixels(std::vector<int>& w, std::vector<int>& h){
  for(unsigned int i = 0; i < w.size(); ++i)
    globalNpixels.push_back(w[i]*h[i]);
}

void setGlobalNpixels(std::vector<int>& d, std::vector<int>& w, std::vector<int>& h){
  for(unsigned int i = 0; i < w.size(); ++i)
    globalNpixels.push_back(w[i]*h[i]*d[i]);
}

void setGlobalWidth(std::vector<int>& v){
  globalwidth = v;
}
// make sure MRF library is compiled with float costs
void checkForFloatCosts() {
  if (typeid(MRF::EnergyVal) != typeid(float) ||
      typeid(MRF::CostVal) != typeid(float)) {
	fprintf(stderr, "Need to define MRF::EnergyVal and MRF::CostVal as float in mrf.h\n");
	exit(1);
  }
}

double getLogSumExp(std::vector<double> v)
{

  std::sort(v.begin(), v.end());

  if (v.size()==1)
    return v[0];

  double ktemp = 1;

  for (int ii=1; ii<v.size(); ii++)
    {
      ktemp = 1 + ktemp*exp(v[ii-1] - v[ii]);
    }  

  ktemp = v[v.size()-1] + log(ktemp);

  return ktemp;

}

/*


void computeAppearanceVectors(std::vector <std::vector <CImage2> > im1,
                              int nD,
                              fvec thetaA,
                              Appearance** appclass,
                              std::vector <std::vector <std::vector <CByteImage> > > &appdirImage,
                              int app)
{
  checkForFloatCosts();

  int nA = thetaA.size();

  int nT;
  if (app==1)
    nT = nA/nD; //this gives centers per class
  else
    nT = appclass[0]->getNTextons();  // since nD = 1 when we pass app==2

  if (nA<=0)
    return;



  for (unsigned int j = 0; j < im1.size(); ++j) {

    int depth = im1[j].size();

    CShape sh = im1[j][0].Shape();
    int width = sh.width, height = sh.height, nB = sh.nBands;
    int tempglobalNpixels = depth * width * height;
    int nColors = __min(3, nB);

    for(unsigned int i = 0; i < im1[j].size(); ++i){

      DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);

      sh.nBands = 1;

      for (int ss=0; ss<nD; ss++)
        appdirImage[j][i][ss].ReAllocate(sh);

      int dsiIndex = 0;
      for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
          for (int d = 0; d < nD; d++) {

            uchar *appIm = &appdirImage[j][i][d].Pixel(x, y, 0);

            // assume the patch size is the same for all classes
            int patchSize = appclass[d]->getPatchSize();
            int nTextons = appclass[d]->getNTextons();

            if (nTextons!=nT)
              {
                printf(" Mismatch in number of textons per class \n");
                exit(1);
              }

            float** textons = appclass[d]->getTextons();

            if (x>(int)patchSize/2 && x<width-patchSize/2 && y>(int)patchSize/2 && y<height-patchSize/2)
              {
                // select best texton center
                int bestT=-1;
                double bestDist=9e99;

                for (int mm=0; mm<nTextons; mm++)
                  {
                    double distT=0;
                    int cnter=0;
                    for (int qq=-(int)patchSize/2 + y; qq<=(int)patchSize/2 + y; qq++)
                      for (int pp=-(int)patchSize/2 + x; pp<=(int)patchSize/2 + x; pp++)
                        {
                          distT += (textons[mm][cnter]-(float)im1[j][i].Pixel(pp, qq, 0)) * (textons[mm][cnter]-(float)im1[j][i].Pixel(pp, qq, 0));
                          cnter++;
                        }

                    if (bestDist>distT)
                      {
                        bestDist = distT;
                        bestT = mm;
                      }

                  }

                if (bestT==-1)
                  {
                    printf(" some error in appearance model \n");
                    exit(1);
                  }

                *appIm = (unsigned char)bestT;
              }
          }
        }
      }
    }
  }
}



// this uses Yuv (for now only uv) and shrinks image to half for computation reasons
// remember k-means and/or EM are trained using this half images to find the centers too
void computeAppearanceVectors(std::vector <std::vector <CByteImage> > im1,
                              int nD,
                              fvec thetaA,
                              Appearance** appclass,
                              std::vector <std::vector <std::vector <CByteImage> > > &appdirImage,
                              int app)
{
  checkForFloatCosts();

  int nA = thetaA.size();

  int nT;
  if (app==1)
    nT = nA/nD; //this gives centers per class
  else
    nT = appclass[0]->getNTextons();  // since nD = 1 when we pass app==2

  if (nA<=0)
    return;



  for (unsigned int j = 0; j < im1.size(); ++j) {

    int depth = im1[j].size();

    CShape sh = im1[j][0].Shape();
    int width = sh.width, height = sh.height, nB = sh.nBands;
    int tempglobalNpixels = depth * width * height;
    int nColors = __min(3, nB);

    for(unsigned int i = 0; i < im1[j].size(); ++i){

      DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);

      sh.nBands = 1;

      for (int ss=0; ss<nD; ss++)
        appdirImage[j][i][ss].ReAllocate(sh);

      // create image buffer of half size
      std::vector < std::vector <std::vector <unsigned char> > > halfImage;
      halfImage.resize(height/2);
      for (int kk=0; kk< height/2; kk++)
      {
    	  halfImage[kk].resize(width/2);
    	  for (int mm=0; mm<width/2; mm++)
    	  {
    		  halfImage[kk][mm].resize(3);
    	  }
      }

      // copy image to make half and convert to Yuv
      for (int kk=0; kk< height/2; kk++)
      {
    	  for (int mm=0; mm<width/2; mm++)
    	  {
    		  double tempR=0, tempG=0, tempB=0;
    		  int rgb[3];
    		  for (int nn=kk*2; nn<=kk*2+1; nn++)
    			  for (int ss=mm*2; ss<=mm*2+1; ss++)
    			  {
    				  if (nn>=height || ss>=width)
    					  continue;

					  tempB += (float)im1[j][i].Pixel(ss, nn, 0);
					  tempG += (float)im1[j][i].Pixel(ss, nn, 1);
					  tempR += (float)im1[j][i].Pixel(ss, nn, 2);
    			  }

			  tempR = tempR/4;
			  tempG = tempG/4;
			  tempB = tempB/4;



			  rgb[0]=round(tempB);
			  rgb[1]=round(tempG);
			  rgb[2]=round(tempR);
			  int Yuv[3];
			  rgb2Yuv(rgb, Yuv);

			  int r=rgb[2], g=rgb[1], b=rgb[0];

			  halfImage[kk][mm][0] = Yuv[0];  //v
			  halfImage[kk][mm][1] = Yuv[1];  //u
			  halfImage[kk][mm][2] = Yuv[2];  //Y
    	  }
      }


      int dsiIndex = 0;
      for (int y = 0; y < height/2; y++) {
        for (int x = 0; x < width/2; x++) {
          for (int d = 0; d < nD; d++) {

            // assume the patch size is the same for all classes
            int patchSize = appclass[d]->getPatchSize();
            int nTextons = appclass[d]->getNTextons();

            if (nTextons!=nT)
              {
                printf(" Mismatch in number of textons per class \n");
                exit(1);
              }

            float** textons = appclass[d]->getTextons();
            // textons are arranged all u then v for uv for a given patch
            // or all Y and then all u and then v for a given patch for Yuv

            // here using only uv for now

            if (x>(int)patchSize/2 && x<width/2-patchSize/2 && y>(int)patchSize/2 && y<height/2-patchSize/2)
              {
                // select best texton center
                int bestT=-1;
                double bestDist=9e99;

                for (int mm=0; mm<nTextons; mm++)
                  {
                    double distT=0;
                    int cnter=0;

                    // u => note texton is stored as [u v]
                    for (int qq=-(int)patchSize/2 + y; qq<=(int)patchSize/2 + y; qq++)
                      for (int pp=-(int)patchSize/2 + x; pp<=(int)patchSize/2 + x; pp++)
                        {
                          distT += (textons[mm][cnter]-halfImage[qq][pp][1]) * (textons[mm][cnter]-halfImage[qq][pp][1]);
                          cnter++;
                        }

                    // v
                    for (int qq=-(int)patchSize/2 + y; qq<=(int)patchSize/2 + y; qq++)
                      for (int pp=-(int)patchSize/2 + x; pp<=(int)patchSize/2 + x; pp++)
                        {
                          distT += (textons[mm][cnter]-halfImage[qq][pp][0]) * (textons[mm][cnter]-halfImage[qq][pp][0]);
                          cnter++;
                        }
                    if (bestDist>distT)
                      {
                        bestDist = distT;
                        bestT = mm;
                      }

                  }

                if (bestT==-1)
                  {
                    printf(" some error in appearance model \n");
                    exit(1);
                  }

        		for (int yy=2*y; yy<=2*y+1; yy++)
        		{
					for(int xx=2*x;xx<=2*x+1; xx++)
					{
            			if (xx>=width || yy>=height)
            				continue;
            			uchar *appIm = &appdirImage[j][i][d].Pixel(xx, yy, 0);
            			*appIm = (unsigned char)bestT;
            		}
            	}

              }
          }
        }
      }
    }
  }
}






void computeAppearanceProb(std::vector <std::vector <CImage2> > im1,
                           int nD,
                           std::vector <std::vector <std::vector <matrixB<double> > > > &appdirMVProb,
                           AppearanceMV** AppMVclass
                           )
{



  for (unsigned int j = 0; j < im1.size(); ++j) {

    int depth = im1[j].size();

    CShape sh = im1[j][0].Shape();
    int width = sh.width, height = sh.height, nB = sh.nBands;

    matrixB<double> improb(height, width);

    for(unsigned int i = 0; i < im1[j].size(); ++i){

      DEBUG_OUT4(verbose, debugfile, "Appearance: Image [%d][%d] \n", j, i);

      for (int d = 0; d < nD; d++) {
        for (int y = 0; y < height; y++) {
          for (int x = 0; x < width; x++) {

            float appcost;
            float *pii;
            float **mu;
            float ***SigmaInv;
            float *normConst;
            int dim;

            dim = AppMVclass[d]->getDimension();
            pii = AppMVclass[d]->getPi();
            mu = AppMVclass[d]->getMu();
            SigmaInv = AppMVclass[d]->getSigmaInv();
            normConst = AppMVclass[d]->getNormConst();


            //					  AppMVclass[d]->printParameters();

            double appval = 0.0;

            int patchSize = sqrt(dim);

            if (x>(int)patchSize/2 && x<width-patchSize/2 && y>(int)patchSize/2 && y<height-patchSize/2)
              {
                for (int c=0; c<AppMVclass[d]->getNClusters(); c++)
                  {
                    int nnn = 0;
                    double tempval=0;
                    for (int xn=x-(int)patchSize/2; xn<=x+(int)patchSize/2; xn++)
                      {
                        for (int yn=y-(int)patchSize/2; yn<=y+(int)patchSize/2; yn++)
                          {
                            unsigned short pixt = im1[j][i].Pixel(xn, yn, 0);

                            double xs = (double)pixt-mu[c][nnn];

                            tempval += xs*xs*SigmaInv[c][nnn][nnn];

                          }
                      }
                    appval += exp(log(pii[c])+(-tempval/2)-log(normConst[c]));
                  }

                if (appval==0.0)
                  appval = pow(10.0, -300.0);  //WORK-THIS

                appcost = -log(appval);

                (improb)(y, x) = appcost;

              } else {
              (improb)(y,x) = -log(1/nD); //uniform distribution
            }
          }
        } //finished filling matrix(y,x)

        appdirMVProb[j][i].push_back(improb);

      }
    }
  }

}



void readHoGProb(std::vector <std::vector <CImage2> > im1, int nD, std::vector <std::vector <std::vector <matrixB<double> > > > &hogdirMVProb, std::vector<int> testDirIndexV)
{


  for (unsigned int j = 0; j < im1.size(); ++j) {

    int depth = im1[j].size();

    CShape sh = im1[j][0].Shape();
    int width = sh.width, height = sh.height, nB = sh.nBands;

    //WORK-THIS - how to take into account multiple patients etc in an elegant way

    //WORK-THIS assume the patient is pa13, will need to modify this code
    // to account for multiple patients

    char cname[50];
    // need to work on testDirIndex
    if (testDirIndexV[j]==1)
      strcpy(cname, "test/");
    else {

      if (j==0)
        strcpy(cname, "pa13/");
      else if (j==1)
        strcpy(cname, "test/"); //WORK-THIS will need to change depend on the data

    }

    for (int d = 0; d < nD; d++) {

      char fname[200];
      strcpy(fname, "../art/HoGgen30/" );
      strcat(fname, cname);

      if (d==0)
        strcat(fname,"hogprobbgnd.txt");
      else if (d==1)
        strcat(fname,"hogprobliv.txt");
      else if (d==2)
        strcat(fname,"hogprobrk.txt");
      else if (d==3)
        strcat(fname,"hogproblk.txt");
      else if (d==4)
        strcat(fname,"hogprobgb.txt");
      else if (d==5)
        strcat(fname,"hogprobspl.txt");

      std::cout << "filename : "<<fname<<endl;

      ifstream myfile (fname);

      if (myfile.is_open())
        {
          string line;

          int z=0, y=0, x=0;  //y is for height or row, x is for width or column
          matrixB<double> improb(height, width);

          while( getline( myfile, line ) )
            {
              tokenizer<escaped_list_separator<char> > tok(line);
              for(tokenizer<escaped_list_separator<char> >::iterator beg=tok.begin(); beg!=tok.end(); ++beg){
                improb(y,x) = boost::lexical_cast<double>(*beg);
                x++;
                if (x>width) {
                  std::cout <<" Error in the data in hog param file \n";
                  exit(1);
                }
              }
              y++; x=0;
              if (y>=height) {
                hogdirMVProb[j][z].push_back(improb);
                z++; y=0; x=0;
              }
            }
          myfile.close();
        } else {
        std::cout << "Unable to open hogprob parameter file";
        exit(1);
      }

    }

  }
}




*/


//this is the initialization using location
void initializeDataCostVector(int genparam, std::vector <std::vector <CImage2> > im1,
                              std::vector <std::vector <CByteImage> > hogIm,
                              int nD,               // number of disparities
                              int numTrainingPats,
                              std::vector <std::vector <CByteImage> > &wta,      // winner-take-all disparities
                              fvec thetaU, fvec thetaA, fvec thetaH, fvec thetaL,
                              int featureCode, int featureCodeKlr,
                              std::vector <std::vector <std::vector <CByteImage> > > appdirImage,
                              Appearance** appclass, LocationMV** LocationMVclass, AppearanceMV** AppMVclass,
                              std::vector <std::vector <std::vector  <matrixB<double> > > > appdirMVProb,
                              intClass **intObj,
                              std::vector <std::vector <CByteImage> > gtIm,
                              std::vector <std::vector <std::vector <CImage2> > > locdirImage,
                              std::vector <std::vector <std::vector  <matrixB<double> > > > hogdirMVProb,
                              intensityNB* intNBclass,
                              matrix<DBL_TYPE> wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM, DBL_TYPE lambda, DBL_TYPE p1,
                              std::vector<int> startSliceNo, std::vector<int> endSliceNo,
                              int nbinsNB, int hogdimklr, int appdimklr, int locdimklr, int biasklr
                              )
{
  checkForFloatCosts();

  globalNdisps = nD;

  int intensity = (int)featureCode/10000;
  int hog = (int) (featureCode/1000)%10;
  int app = (int) (featureCode/100)%10;
  int loc = (int) (featureCode/10)%10;

  int klr = 0;

  int intensityKlr = (int)featureCodeKlr/10000;
  int hogKlr = (int) (featureCodeKlr/1000)%10;
  int appKlr = (int) (featureCodeKlr/100)%10;
  int locKlr = (int) (featureCodeKlr/10)%10;

  if (intensityKlr!=0 || appKlr!=0 || hogKlr!=0 || locKlr!=0)
    klr = 1;


  if (klr==1)
    {
      if (intensityKlr>0)
  		intensity=0;
      if (appKlr>0)
  		app=0;
      if (locKlr>0)
  		loc=0;
      if (hogKlr>0)
  		hog=0;
    }


  int dim = 0;
  IVM* ivm;


  if (intensityKlr==2 || hogKlr==2 || appKlr==2 || locKlr==2)
    {
	  klr = 1;
	  ivm = new IVM(&Xtrain, &wparam, RBF, lambda, nD, p1);
	  dim = Xtrain.size2();
    }



  int nU = thetaU.size();
  int nbins = nU/nD;

  // assuming that we have png data that goes upto 65535.
  // since i have the data in the range of 850 to 1250 which is shifted to 0 to 400
  double binner = 256.0/nbins;
  // so for any point. divide it by binner to get bin number (floor int value) it is in
  // ranging from bin0 to bin{nbins-1}
  // since >=400 will be special case, put it in the last bin if there.

  int nA = thetaA.size();

  int nT;

  nT = nA/nD; //this gives centers per class
  // this will work for app==2 and app==1 because in app==2, we have nD*the number of total clusters as the number of parameters
  // while in app==1 it is nD*number of clusters of that class (only that clusters are also = number of params)



  int nH = thetaH.size(); // number of HoG parameters
  int nHpC = nH/nD; // HoG vocabulary size

  int nL = thetaL.size(); // number of location parameters
  int nLpC = nL/nD;  // this variable will be useful only for cube setting (loc==1) or loc==9

  MRF::CostVal badcost = 1e20;

  for(int j=0; j<nU; j++)
    badcost += thetaU[j];

  if (nA > 0)
    for(int j=0; j<nA; j++)
      badcost += thetaA[j];

  if (nH > 0)
    for(int j=0; j<nH; j++)
      badcost += thetaH[j];

  if (nL > 0)
    for(int j=0; j<nL; j++)
      badcost += thetaL[j];




  if (genparam==1) {
    int numGens=0;
    if (intensity==3 || intensity==6)
      numGens++;
    if (loc==3 || loc==8 || loc==9)
      numGens++;
    if (app==3)
      numGens++;
    if (hog==3)
      numGens++;


    //WORK-THIS
    // assuming worst cost after taking log is -log(10^-99)

    badcost += numGens*(99*3);  //that 3 accounts for 2.7 value of e (natural log)

  }




  for (unsigned int j = 0; j < im1.size(); ++j) {

    int depth = im1[j].size();
    CShape sh = im1[j][0].Shape();
    int width = sh.width, height = sh.height, nB = sh.nBands;
    int tempglobalNpixels = depth * width * height;
    int nColors = __min(3, nB);

    if (genparam==0)
      dsiArrayV.push_back(new MRF::CostVal[tempglobalNpixels * globalNdisps]);
    else if (genparam==1)
      dsiArrayVGen.push_back(new MRF::CostVal[tempglobalNpixels * globalNdisps]);

    int zslice = 0;
    zslice = startSliceNo[j];


    int dsiIndex = 0;
    for(unsigned int i = 0; i < im1[j].size(); ++i){

      DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);
      sh.nBands = 1;
      wta[j][i].ReAllocate(sh);

      matrix<DBL_TYPE> hogdescs;

      if (hogKlr==2)
        {
          // read file
          char hogfile[100];

          sprintf(hogfile, "/p/imageseg/medicalSeg/parameters/hog/hog2dcsv_90n/hog2ddata_%d.csv", i+zslice );

          // read this file in matrix
          if(!load_one_data(hogfile, hogdescs)){
            std::cout << "csv read error: Bad data file filename : "<< hogfile << endl;
            exit(1);
          }

        }

      // need to load appfile if appklr==2
      matrix<DBL_TYPE> appdescs;

      if (appKlr==2)
        {
          // read file
          char appfile[100];

          sprintf(appfile, "/p/imageseg/medicalSeg/parameters/app/appcsv_81n/app2ddatan_%d.csv", i+zslice );

          // read this file in matrix
          if(!load_one_data(appfile, appdescs)){
            std::cout << "csv read error: Bad data file filename : "<< appfile << endl;
            exit(1);
          }

        }

      for (int y = 0; y < height; y++) {
        uchar *WTArow = &wta[j][i].Pixel(0, y, 0);
        for (int x = 0; x < width; x++) {
          unsigned short pix1 = im1[j][i].Pixel(x, y, 0);

          uchar *gtpix = &gtIm[j][i].Pixel(x, y, 0);

          uchar *hogval;
          if (nH > 0)
            hogval = &hogIm[j][i].Pixel(x, y, 0);

          ublas::vector<DBL_TYPE> f;
          MRF::CostVal bestval = badcost;
          int* bestd = new int [nD];
          int numbest = 0;

          int thetaid;

          if (klr==1)
            {
              ublas::vector<DBL_TYPE> Xtest(dim);

              int idxx = 0;

              if (intensityKlr==2) // for now only intensity case
                {
                  getCostIntensityKlr(Xtest, (double)pix1, minmaxM, idxx);
                }

              if (locKlr==2)
                {
                  // 25 values. // this should be a parameter  // location features
                  // x y z <prob of belonging to each cluster>*16 <prob of belong to each class>*6
                  getCostLocationKlr(Xtest, i, zslice, x, y, nD, minmaxM, LocationMVclass, idxx, skipXYZlocklr);

                }

              if (appKlr==2)
                {
                  getCostAppKlr(Xtest, x, y, width, minmaxM, idxx, appdimklr, &appdescs);
                }


              if (hogKlr==2)
                {
                  getCostHoGKlr(Xtest, x, y, width, minmaxM, idxx, hogdimklr, &hogdescs);
                }

              std::cout<<gtpix[0]+1;
              for (int kop=0; kop<dim; kop++)
              {
            	  std::cout<<" "<<kop+1<<":"<<Xtest(kop);
              }
              std::cout<<"\n";


              f = ivm->classify_klrexp(&Xtest);

            }




          for (int d = 0; d < nD; d++) {

            MRF::CostVal dsiValue = 0;

            if (klr==1) {
              dsiValue += -f(d);
              dsiValue += biasklr;
            }


            if (nU>0 && intensity==1)  // intensity as histogram bins
              {
                dsiValue += getCostIntensityDiscBins(pix1, nbins, d, thetaU, thetaid);

              } else if (intensity==3 && genparam==1) // generative setting with Gaussians
              {
                dsiValue += getCostIntensityGenGaussian(intObj[d], pix1);

              } else if (intensity == 6 && genparam==1) // using NBayes and bins
              {
                dsiValue += getCostIntensityGenNB(pix1, nbins, d, intNBclass);
              }


            if (nA>0 && app==1) // we have appearance features
              {
                // different features for different classes
					   
                dsiValue += getCostAppDiscPatch((int)appdirImage[j][i][d].Pixel(x, y, 0), appclass[d]->getPatchSize(), x, y, width, height, thetaA, d, nT, thetaid);

              } else if (nA>0 && app==2) // we have appearance features
              {
                // common features for all classes
					    
                dsiValue += getCostAppDiscPatch((int)appdirImage[j][i][0].Pixel(x, y, 0), appclass[0]->getPatchSize(), x, y, width, height, thetaA, d, nT, thetaid);

              }
            else if (app==3 && genparam==1) // generative setting
              {

                dsiValue += getCostAppGenPatch((appdirMVProb[j][i][d])(y,x), AppMVclass[d]->getDimension(), x, y, width, height);

              }



            if (nH > 0 && hog==1) {
					    
              dsiValue += getCostHoGDiscBins((int) *hogval, d, nHpC, thetaH, thetaid);

            } else if (hog==3 && genparam==1) // generative setting
              {

                int startPatch = 8; // index starting from 0
                // this variable indicates where the Hog descriptors are defined
                // since i have generated them using the entire patient body, don't need to
                // worry about the z axis and can use all slices here.
                //WORK-THIS need to change above stuff for flawless working

                dsiValue += getCostHoGGenBins((hogdirMVProb[j][i][d])(y,x), x, y, width, height, startPatch);

              }

            if (nL>0 && loc==1) { // CRF cube code

              // this function is completely out-dated
              // and for older series of data, will need to modify significantly for use
					    
              dsiValue += getCostLocationDiscCube((unsigned short) locdirImage[j][d][i].Pixel(x,y,0), x, y, i, loccubex, loccubey, loccubez, numTrainingPats, nLpC, thetaL, width, height, depth, d, thetaid);

            } else if (nL>0 && (loc==6 || loc==8) )  //CRF hard or soft assignment
              {
					    
                dsiValue += getCostLocationDiscGaussian(LocationMVclass, i, zslice, x, y, thetaL, d, loc, genparam, thetaid);

              } else if (nL>0 && (loc==7 || loc==9) )  //CRF soft or hard assignment redefined
              {
					    
                dsiValue += getCostLocationDiscGaussian2(LocationMVclass, i, zslice, x, y, thetaL, d, loc, genparam, nD, thetaid);

              } else if (loc==3 && genparam==1) // generative setting
              {

                dsiValue += getCostLocationGenGaussian(LocationMVclass[d], i, zslice, x, y);

              }

            // The cost of pixel p and label l is stored at dsiArray[p*nLabels+l]
            if (genparam==0)
              (dsiArrayV[j])[dsiIndex++] = (float)dsiValue;
            else if (genparam==1)
              (dsiArrayVGen[j])[dsiIndex++] = (float)dsiValue;

            if (dsiValue < bestval) {
              bestval = dsiValue;
              bestd[0] = d;
              numbest = 1;
            } else if (dsiValue == bestval) {
              bestd[numbest] = d;
              numbest++;
            }

            if (pix1 >= 256) {
              bestd[0] = 0; // assume 0 is background
              numbest = 1;
            }
          }

          int curr = rand() % numbest;
          WTArow[x] = bestd[curr];

          delete [] bestd;
        }
      }
    }
  }



  if (klr==1)
    delete ivm;


}


void initializeDataCostVector(int genparam, std::vector <std::vector <CImage2> > im1,
                              std::vector <std::vector <CByteImage> > hogIm,
                              std::vector <std::vector <CByteImage> > motIm,
                              int numTrainingPats,
                              std::vector <std::vector <CByteImage> > &wta,      // winner-take-all disparities
                              std::vector <std::vector <std::vector <CByteImage> > > appdirImage,
                              Appearance** appclass, LocationMV** LocationMVclass, AppearanceMV** AppMVclass,
                              std::vector <std::vector <std::vector  <matrixB<double> > > > appdirMVProb,
                              intClass **intObj,
                              std::vector <std::vector <CByteImage> > gtIm,
                              std::vector <std::vector <std::vector <CImage2> > > locdirImage,
                              std::vector <std::vector <std::vector  <matrixB<double> > > > hogdirMVProb,
                              intensityNB* intNBclass,
                              matrix<DBL_TYPE> wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM,
                              std::vector<int> startSliceNo, std::vector<int> endSliceNo, Parameters prm
                              )
{

  int nD = prm.nD;
  std::vector<float> thetaU, thetaA, thetaH, thetaL, thetaM;
    // copy thetas
  vecCopy(prm.featurepv.thetaU, thetaU);
  vecCopy(prm.featurepv.thetaA, thetaA);
  vecCopy(prm.featurepv.thetaH, thetaH);
  vecCopy(prm.featurepv.thetaM, thetaM);
  vecCopy(prm.featurepv.thetaL, thetaL);

  int featureCode = prm.featurepv.featureCode;
  int featureCodeKlr = prm.featurepv.featureCodeKlr;

  DBL_TYPE lambda = prm.featurepv.klrpv.lambda;
  DBL_TYPE p1 = prm.featurepv.klrpv.p1;

  int nbinsNB = prm.featurepv.nbinsintNB;
  int hogdimklr = prm.featurepv.klrpv.hogdimklr;
  int appdimklr = prm.featurepv.klrpv.appdimklr;
  int locdimklr = prm.featurepv.klrpv.locdimklr;
  int biasklr = prm.featurepv.klrpv.biasklr;


  checkForFloatCosts();

  globalNdisps = nD;

  int intensity = (int)featureCode/1000000;
  int hog = (int) (featureCode/100000)%10;
  int app = (int) (featureCode/10000)%10;
  int loc = (int) (featureCode/1000)%10;
  int mot = (int) (featureCode/100)%10;
  int rbm = (int) (featureCode/10)%10;



  int klr = 0;

  int intensityKlr = (int)featureCodeKlr/10000;
  int hogKlr = (int) (featureCodeKlr/1000)%10;
  int appKlr = (int) (featureCodeKlr/100)%10;
  int locKlr = (int) (featureCodeKlr/10)%10;

  if (intensityKlr!=0 || appKlr!=0 || hogKlr!=0 || locKlr!=0)
    klr = 1;


  if (klr==1)
    {
      if (intensityKlr>0)
  		intensity=0;
      if (appKlr>0)
  		app=0;
      if (locKlr>0)
  		loc=0;
      if (hogKlr>0)
  		hog=0;
    }


  int dim = 0;
  IVM* ivm;


  if (intensityKlr==2 || hogKlr==2 || appKlr==2 || locKlr==2)
    {
	  klr = 1;
	  ivm = new IVM(&Xtrain, &wparam, RBF, lambda, nD, p1);
	  dim = Xtrain.size2();
    }



  int nU = thetaU.size();
  int nbins = nU/nD;

  // assuming that we have png data that goes upto 65535.
  // since i have the data in the range of 850 to 1250 which is shifted to 0 to 400
  double binner = 256.0/nbins;
  // so for any point. divide it by binner to get bin number (floor int value) it is in
  // ranging from bin0 to bin{nbins-1}
  // since >=400 will be special case, put it in the last bin if there.

  int nA = thetaA.size();

  int nT;

  nT = nA/nD; //this gives centers per class
  // this will work for app==2 and app==1 because in app==2, we have nD*the number of total clusters as the number of parameters
  // while in app==1 it is nD*number of clusters of that class (only that clusters are also = number of params)



  int nH = thetaH.size(); // number of HoG parameters
  int nHpC = nH/nD; // HoG vocabulary size


  int nM = thetaM.size(); // number of HoG parameters
  int nMpC = nM/nD; // HoG vocabulary size


  int nL = thetaL.size(); // number of location parameters
  int nLpC = nL/nD;  // this variable will be useful only for cube setting (loc==1) or loc==9

  MRF::CostVal badcost = 1e20;

  for(int j=0; j<nU; j++)
    badcost += thetaU[j];

  if (nA > 0)
    for(int j=0; j<nA; j++)
      badcost += thetaA[j];

  if (nH > 0)
    for(int j=0; j<nH; j++)
      badcost += thetaH[j];

  if (nM > 0)
    for(int j=0; j<nM; j++)
      badcost += thetaM[j];

  if (nL > 0)
    for(int j=0; j<nL; j++)
      badcost += thetaL[j];




  if (genparam==1) {
    int numGens=0;
    if (intensity==3 || intensity==6)
      numGens++;
    if (loc==3 || loc==8 || loc==9)
      numGens++;
    if (app==3)
      numGens++;
    if (hog==3)
      numGens++;

    // mot not implemented for generative model


    //WORK-THIS
    // assuming worst cost after taking log is -log(10^-99)

    badcost += numGens*(99*3);  //that 3 accounts for 2.7 value of e (natural log)

  }




  for (unsigned int j = 0; j < im1.size(); ++j) {

    int depth = im1[j].size();
    CShape sh = im1[j][0].Shape();
    int width = sh.width, height = sh.height, nB = sh.nBands;
    int tempglobalNpixels = depth * width * height;
    int nColors = __min(3, nB);

    if (genparam==0)
      dsiArrayV.push_back(new MRF::CostVal[tempglobalNpixels * globalNdisps]);
    else if (genparam==1)
      dsiArrayVGen.push_back(new MRF::CostVal[tempglobalNpixels * globalNdisps]);

    int zslice = 0;
    zslice = startSliceNo[j];


    int dsiIndex = 0;
    for(unsigned int i = 0; i < im1[j].size(); ++i){

      DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);
      sh.nBands = 1;
      wta[j][i].ReAllocate(sh);

      matrix<DBL_TYPE> hogdescs;

      if (hogKlr==2)
        {
          // read file
          char hogfile[100];

          sprintf(hogfile, "/p/imageseg/medicalSeg/parameters/hog/hog2dcsv_90n/hog2ddata_%d.csv", i+zslice );

          // read this file in matrix
          if(!load_one_data(hogfile, hogdescs)){
            std::cout << "csv read error: Bad data file filename : "<< hogfile << endl;
            exit(1);
          }

        }

      // need to load appfile if appklr==2
      matrix<DBL_TYPE> appdescs;

      if (appKlr==2)
        {
          // read file
          char appfile[100];

          sprintf(appfile, "/p/imageseg/medicalSeg/parameters/app/appcsv_81n/app2ddatan_%d.csv", i+zslice );

          // read this file in matrix
          if(!load_one_data(appfile, appdescs)){
            std::cout << "csv read error: Bad data file filename : "<< appfile << endl;
            exit(1);
          }

        }

      for (int y = 0; y < height; y++) {
        uchar *WTArow = &wta[j][i].Pixel(0, y, 0);
        for (int x = 0; x < width; x++) {
          unsigned short pix1 = im1[j][i].Pixel(x, y, 0);

          uchar *gtpix = &gtIm[j][i].Pixel(x, y, 0);

          uchar *hogval;
          if (nH > 0)
            hogval = &hogIm[j][i].Pixel(x, y, 0);

          uchar *motval;
          if (nM > 0)
            motval = &motIm[j][i].Pixel(x, y, 0);


          ublas::vector<DBL_TYPE> f;
          MRF::CostVal bestval = badcost;
          int* bestd = new int [nD];
          int numbest = 0;

          int thetaid;

          if (klr==1)
            {
              ublas::vector<DBL_TYPE> Xtest(dim);

              int idxx = 0;

              if (intensityKlr==2) // for now only intensity case
                {
                  getCostIntensityKlr(Xtest, (double)pix1, minmaxM, idxx);
                }

              if (locKlr==2)
                {
                  // 25 values. // this should be a parameter  // location features
                  // x y z <prob of belonging to each cluster>*16 <prob of belong to each class>*6
                  getCostLocationKlr(Xtest, i, zslice, x, y, nD, minmaxM, LocationMVclass, idxx, skipXYZlocklr);

                }

              if (appKlr==2)
                {
                  getCostAppKlr(Xtest, x, y, width, minmaxM, idxx, appdimklr, &appdescs);
                }


              if (hogKlr==2)
                {
                  getCostHoGKlr(Xtest, x, y, width, minmaxM, idxx, hogdimklr, &hogdescs);
                }

              std::cout<<gtpix[0]+1;
              for (int kop=0; kop<dim; kop++)
              {
            	  std::cout<<" "<<kop+1<<":"<<Xtest(kop);
              }
              std::cout<<"\n";


              f = ivm->classify_klrexp(&Xtest);

            }




          for (int d = 0; d < nD; d++) {

            MRF::CostVal dsiValue = 0;

            if (klr==1) {
              dsiValue += -f(d);
              dsiValue += biasklr;
            }


            if (nU>0 && intensity==1)  // intensity as histogram bins
              {
                dsiValue += getCostIntensityDiscBins(pix1, nbins, d, thetaU, thetaid);

              } else if (intensity==3 && genparam==1) // generative setting with Gaussians
              {
                dsiValue += getCostIntensityGenGaussian(intObj[d], pix1);

              } else if (intensity == 6 && genparam==1) // using NBayes and bins
              {
                dsiValue += getCostIntensityGenNB(pix1, nbinsNB, d, intNBclass);
              }


            if (nA>0 && app==1) // we have appearance features
              {
                // different features for different classes

                dsiValue += getCostAppDiscPatch((int)appdirImage[j][i][d].Pixel(x, y, 0), appclass[d]->getPatchSize(), x, y, width, height, thetaA, d, nT, thetaid);

              } else if (nA>0 && app==2) // we have appearance features
              {
                // common features for all classes

                dsiValue += getCostAppDiscPatch((int)appdirImage[j][i][0].Pixel(x, y, 0), appclass[0]->getPatchSize(), x, y, width, height, thetaA, d, nT, thetaid);

              }
            else if (app==3 && genparam==1) // generative setting
              {

                dsiValue += getCostAppGenPatch((appdirMVProb[j][i][d])(y,x), AppMVclass[d]->getDimension(), x, y, width, height);

              }



            if (nH > 0 && hog==1) {

              dsiValue += getCostHoGDiscBins((int) *hogval, d, nHpC, thetaH, thetaid);

            } else if (hog==3 && genparam==1) // generative setting
              {

                int startPatch = 8; // index starting from 0
                // this variable indicates where the Hog descriptors are defined
                // since i have generated them using the entire patient body, don't need to
                // worry about the z axis and can use all slices here.
                //WORK-THIS need to change above stuff for flawless working

                dsiValue += getCostHoGGenBins((hogdirMVProb[j][i][d])(y,x), x, y, width, height, startPatch);

              }



            if (nM > 0 && mot==1) {

              dsiValue += getCostmotDiscBins((int) *motval, d, nMpC, thetaM, thetaid);

            } else if (mot==3 && genparam==1) // generative setting
            {
            	std::cout << " mot generative not implemented! \n";
            	exit(-1);
            }





            if (nL>0 && loc==1) { // CRF cube code

              // this function is completely out-dated
              // and for older series of data, will need to modify significantly for use

              dsiValue += getCostLocationDiscCube((unsigned short) locdirImage[j][d][i].Pixel(x,y,0), x, y, i, loccubex, loccubey, loccubez, numTrainingPats, nLpC, thetaL, width, height, depth, d, thetaid);

            } else if (nL>0 && (loc==6 || loc==8) )  //CRF hard or soft assignment
              {

                dsiValue += getCostLocationDiscGaussian(LocationMVclass, i, zslice, x, y, thetaL, d, loc, genparam, thetaid);

              } else if (nL>0 && (loc==7 || loc==9) )  //CRF soft or hard assignment redefined
              {

                dsiValue += getCostLocationDiscGaussian2(LocationMVclass, i, zslice, x, y, thetaL, d, loc, genparam, nD, thetaid);

              } else if (loc==3 && genparam==1) // generative setting
              {

                dsiValue += getCostLocationGenGaussian(LocationMVclass[d], i, zslice, x, y);

              }

            //std::cout << " j, i, y, x, d, dsi : " << j << "  " <<  i << "  " << y << "  " << x << "  " << d << "  " << dsiValue << std::endl;

            // The cost of pixel p and label l is stored at dsiArray[p*nLabels+l]
            if (genparam==0)
              (dsiArrayV[j])[dsiIndex++] = (float)dsiValue;
            else if (genparam==1)
              (dsiArrayVGen[j])[dsiIndex++] = (float)dsiValue;

            if (dsiValue < bestval) {
              bestval = dsiValue;
              bestd[0] = d;
              numbest = 1;
            } else if (dsiValue == bestval) {
              bestd[numbest] = d;
              numbest++;
            }

            if (pix1 >= 256) {
              bestd[0] = 0; // assume 0 is background
              numbest = 1;
            }
          }

          if (numbest == 0)
        	  numbest = nD;

          int curr = rand() % numbest;
          WTArow[x] = bestd[curr];

          delete [] bestd;
        }
      }
    }
  }



  if (klr==1)
    delete ivm;


}



void initializeDataCostVector(int genparam, std::vector <std::vector <CByteImage> > im1,
                              std::vector <std::vector <CByteImage> > hogIm,
                              std::vector <std::vector <CByteImage> > motIm,
                              int numTrainingPats,
                              std::vector <std::vector <CByteImage> > &wta,      // winner-take-all disparities
                              std::vector <std::vector <std::vector <CByteImage> > > appdirImage,
                              Appearance** appclass, LocationMV** LocationMVclass, AppearanceMV** AppMVclass,
                              std::vector <std::vector <std::vector  <matrixB<double> > > > appdirMVProb,
                              intClass **intObj,
                              std::vector <std::vector <CByteImage> > gtIm,
                              std::vector <std::vector <std::vector <CImage2> > > locdirImage,
                              std::vector <std::vector <std::vector  <matrixB<double> > > > hogdirMVProb,
                              intensityNB* intNBclass,
                              matrix<DBL_TYPE> wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM,
                              std::vector<int> startSliceNo, std::vector<int> endSliceNo, Parameters prm
                              )
{

  int nD = prm.nD;
  std::vector<float> thetaU, thetaA, thetaH, thetaM, thetaL;
    // copy thetas
  vecCopy(prm.featurepv.thetaU, thetaU);
  vecCopy(prm.featurepv.thetaA, thetaA);
  vecCopy(prm.featurepv.thetaH, thetaH);
  vecCopy(prm.featurepv.thetaM, thetaM);
  vecCopy(prm.featurepv.thetaL, thetaL);

  int featureCode = prm.featurepv.featureCode;
  int featureCodeKlr = prm.featurepv.featureCodeKlr;

  DBL_TYPE lambda = prm.featurepv.klrpv.lambda;
  DBL_TYPE p1 = prm.featurepv.klrpv.p1;

  int nbinsNB = prm.featurepv.nbinsintNB;
  int hogdimklr = prm.featurepv.klrpv.hogdimklr;
  int appdimklr = prm.featurepv.klrpv.appdimklr;
  int locdimklr = prm.featurepv.klrpv.locdimklr;
  int biasklr = prm.featurepv.klrpv.biasklr;


  checkForFloatCosts();

  globalNdisps = nD;

  int intensity = (int)featureCode/1000000;
  int hog = (int) (featureCode/100000)%10;
  int app = (int) (featureCode/10000)%10;
  int loc = (int) (featureCode/1000)%10;
  int mot = (int) (featureCode/100)%10;
  int rbm = (int) (featureCode/10)%10;




  int klr = 0;

  int intensityKlr = (int)featureCodeKlr/10000;
  int hogKlr = (int) (featureCodeKlr/1000)%10;
  int appKlr = (int) (featureCodeKlr/100)%10;
  int locKlr = (int) (featureCodeKlr/10)%10;

  if (intensityKlr!=0 || appKlr!=0 || hogKlr!=0 || locKlr!=0)
    klr = 1;


  if (klr==1)
    {
      if (intensityKlr>0)
  		intensity=0;
      if (appKlr>0)
  		app=0;
      if (locKlr>0)
  		loc=0;
      if (hogKlr>0)
  		hog=0;
    }


  int dim = 0;
  IVM* ivm;


  if (intensityKlr==2 || hogKlr==2 || appKlr==2 || locKlr==2)
    {
	  klr = 1;
	  ivm = new IVM(&Xtrain, &wparam, RBF, lambda, nD, p1);
	  dim = Xtrain.size2();
    }



  int nU = thetaU.size();
  int nbins = nU/nD;

  // assuming that we have png data that goes upto 65535.
  // since i have the data in the range of 850 to 1250 which is shifted to 0 to 400
  double binner = 256.0/nbins;
  // so for any point. divide it by binner to get bin number (floor int value) it is in
  // ranging from bin0 to bin{nbins-1}
  // since >=400 will be special case, put it in the last bin if there.

  int nA = thetaA.size();

  int nT;

  nT = nA/nD; //this gives centers per class
  // this will work for app==2 and app==1 because in app==2, we have nD*the number of total clusters as the number of parameters
  // while in app==1 it is nD*number of clusters of that class (only that clusters are also = number of params)



  int nH = thetaH.size(); // number of HoG parameters
  int nHpC = nH/nD; // HoG vocabulary size


  int nM = thetaM.size(); // number of mot parameters
  int nMpC = nM/nD; // HoG vocabulary size


  int nL = thetaL.size(); // number of location parameters
  int nLpC = nL/nD;  // this variable will be useful only for cube setting (loc==1) or loc==9

  MRF::CostVal badcost = 1e20;

  for(int j=0; j<nU; j++)
    badcost += thetaU[j];

  if (nA > 0)
    for(int j=0; j<nA; j++)
      badcost += thetaA[j];

  if (nH > 0)
    for(int j=0; j<nH; j++)
      badcost += thetaH[j];

  if (nM > 0)
    for(int j=0; j<nM; j++)
      badcost += thetaM[j];


  if (nL > 0)
    for(int j=0; j<nL; j++)
      badcost += thetaL[j];




  if (genparam==1) {
    int numGens=0;
    if (intensity==3 || intensity==6)
      numGens++;
    if (loc==3 || loc==8 || loc==9)
      numGens++;
    if (app==3)
      numGens++;
    if (hog==3)
      numGens++;
    // not coded for mot - generative model


    //WORK-THIS
    // assuming worst cost after taking log is -log(10^-99)

    badcost += numGens*(99*3);  //that 3 accounts for 2.7 value of e (natural log)

  }




  for (unsigned int j = 0; j < im1.size(); ++j) {

    int depth = im1[j].size();
    CShape sh = im1[j][0].Shape();
    int width = sh.width, height = sh.height, nB = sh.nBands;
    int tempglobalNpixels = depth * width * height;
    int nColors = __min(3, nB);

    if (genparam==0)
      dsiArrayV.push_back(new MRF::CostVal[tempglobalNpixels * globalNdisps]);
    else if (genparam==1)
      dsiArrayVGen.push_back(new MRF::CostVal[tempglobalNpixels * globalNdisps]);

    int zslice = 0;
    zslice = startSliceNo[j];

//    if (j==1)
//    	exit(0);

    int dsiIndex = 0;
    for(unsigned int i = 0; i < im1[j].size(); ++i){

      DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);
      sh.nBands = 1;
      wta[j][i].ReAllocate(sh);

      matrix<DBL_TYPE> hogdescs;

      if (hogKlr==2)
        {
          // read file
          char hogfile[100];

          sprintf(hogfile, "/p/imageseg/medicalSeg/parameters/hog/hog2dcsv_90n/hog2ddata_%d.csv", i+zslice );

          // read this file in matrix
          if(!load_one_data(hogfile, hogdescs)){
            std::cout << "csv read error: Bad data file filename : "<< hogfile << endl;
            exit(1);
          }

        }

      // need to load appfile if appklr==2
      matrix<DBL_TYPE> appdescs;

      if (appKlr==2)
        {
          // read file
          char appfile[100];

          sprintf(appfile, "/p/imageseg/medicalSeg/parameters/app/appcsv_81n/app2ddatan_%d.csv", i+zslice );

          // read this file in matrix
          if(!load_one_data(appfile, appdescs)){
            std::cout << "csv read error: Bad data file filename : "<< appfile << endl;
            exit(1);
          }

        }

      for (int y = 0; y < height; y++) {
        uchar *WTArow = &wta[j][i].Pixel(0, y, 0);
        for (int x = 0; x < width; x++) {
          uchar *pix1 = &im1[j][i].Pixel(x, y, 0);

 //         std::cout<<j<<" "<<i<<" "<<(int)pix1[0]<<" "<<(int)pix1[1]<<" "<<(int)pix1[2]<<"\n";

          uchar *gtpix = &gtIm[j][i].Pixel(x, y, 0);

          uchar *hogval;
          if (nH > 0)
            hogval = &hogIm[j][i].Pixel(x, y, 0);


          uchar *motval;
          if (nM > 0)
            motval = &motIm[j][i].Pixel(x, y, 0);


          ublas::vector<DBL_TYPE> f;
          MRF::CostVal bestval = badcost;
          int* bestd = new int [nD];
          int numbest = 0;

          int thetaid;
          int thetaidrgb[3]={0,0,0};
          int thetaidYuv[2]={0,0};

          if (klr==1)
            {
              ublas::vector<DBL_TYPE> Xtest(dim);

              int idxx = 0;

              if (intensityKlr==2) // for now only intensity case
                {
                  getCostIntensityKlrYuv(Xtest, pix1, minmaxM, idxx);
                }

              if (locKlr==2)
                {
                  // 25 values. // this should be a parameter  // location features
                  // x y z <prob of belonging to each cluster>*16 <prob of belong to each class>*6
                  getCostLocationKlr(Xtest, i, zslice, x, y, nD, minmaxM, LocationMVclass, idxx, skipXYZlocklr);

                }

              if (appKlr==2)
                {
                  getCostAppKlr(Xtest, x, y, width, minmaxM, idxx, appdimklr, &appdescs);
                }


              if (hogKlr==2)
                {
                  getCostHoGKlr(Xtest, x, y, width, minmaxM, idxx, hogdimklr, &hogdescs);
                }

/*              std::cout<<gtpix[0]+1;
              for (int kop=0; kop<dim; kop++)
              {
            	  std::cout<<" "<<kop+1<<":"<<Xtest(kop);
              }
              std::cout<<"\n";
*/

              f = ivm->classify_klrexp(&Xtest);

            }


//          std::cout << "  \n";


          for (int d = 0; d < nD; d++) {

            MRF::CostVal dsiValue = 0;

            if (klr==1) {
              dsiValue += -f(d);
              dsiValue += biasklr;
            }

//            std::cout << dsiValue << " ";


            if (nU>0 && intensity==1)  // intensity as histogram bins
              {
                //dsiValue += getCostIntensityDiscBins(pix1, nbins, d, thetaU, thetaid);
                // dsiValue += getCostIntensityDiscBins(pix1, nbins, d, thetaU, thetaidrgb);
                dsiValue += getCostIntensityDiscBinsYuv(pix1, nbins, d, thetaU, thetaidYuv);

              } else if (intensity==3 && genparam==1) // generative setting with Gaussians
              {
                //dsiValue += getCostIntensityGenGaussian(intObj[d], pix1);

              } else if (intensity == 6 && genparam==1) // using NBayes and bins
              {
                //dsiValue += getCostIntensityGenNB(pix1, nbins, d, intNBclass);
              }


            if (nA>0 && app==1) // we have appearance features
              {
                // different features for different classes
            	  if ((int)appdirImage[j][i][d].Pixel(x, y, 0) >= nT)
            	  {
            		  // this is to handle the possibility that some of the pixels in the image borders might have invalid values
            	  } else
            	  {
            		  dsiValue += getCostAppDiscPatch((int)appdirImage[j][i][d].Pixel(x, y, 0), appclass[d]->getPatchSize(), x, y, width, height, thetaA, d, nT, thetaid);
            	  }
              } else if (nA>0 && app==2) // we have appearance features
              {
                // common features for all classes
               	  if ((int)appdirImage[j][i][0].Pixel(x, y, 0) >= nA)
               	  {
                		  // this is to handle the possibility that some of the pixels in the image borders might have invalid values
               	  } else
               	  {
               		  dsiValue += getCostAppDiscPatch((int)appdirImage[j][i][0].Pixel(x, y, 0), appclass[0]->getPatchSize(), x, y, width, height, thetaA, d, nT, thetaid);
               	  }
              }
            else if (app==3 && genparam==1) // generative setting
              {

                dsiValue += getCostAppGenPatch((appdirMVProb[j][i][d])(y,x), AppMVclass[d]->getDimension(), x, y, width, height);

              }



            if (nH > 0 && hog==1) {

              dsiValue += getCostHoGDiscBins((int) *hogval, d, nHpC, thetaH, thetaid);

            } else if (hog==3 && genparam==1) // generative setting
              {

                int startPatch = 8; // index starting from 0
                // this variable indicates where the Hog descriptors are defined
                // since i have generated them using the entire patient body, don't need to
                // worry about the z axis and can use all slices here.
                //WORK-THIS need to change above stuff for flawless working

                dsiValue += getCostHoGGenBins((hogdirMVProb[j][i][d])(y,x), x, y, width, height, startPatch);

              }


            if (nM > 0 && mot==1) {

              dsiValue += getCostmotDiscBins((int) *motval, d, nMpC, thetaM, thetaid);

            } else if (mot==3 && genparam==1) // generative setting
            {
            	std::cout << " mog generative not implemented \n";
				exit(-1);
            }



            if (nL>0 && loc==1) { // CRF cube code

              // this function is completely out-dated
              // and for older series of data, will need to modify significantly for use

              dsiValue += getCostLocationDiscCube((unsigned short) locdirImage[j][d][i].Pixel(x,y,0), x, y, i, loccubex, loccubey, loccubez, numTrainingPats, nLpC, thetaL, width, height, depth, d, thetaid);

            } else if (nL>0 && (loc==6 || loc==8) )  //CRF hard or soft assignment
              {

                dsiValue += getCostLocationDiscGaussian(LocationMVclass, i, zslice, x, y, thetaL, d, loc, genparam, thetaid);

              } else if (nL>0 && (loc==7 || loc==9) )  //CRF soft or hard assignment redefined
              {

                dsiValue += getCostLocationDiscGaussian2(LocationMVclass, i, zslice, x, y, thetaL, d, loc, genparam, nD, thetaid);

              } else if (loc==3 && genparam==1) // generative setting
              {

                dsiValue += getCostLocationGenGaussian(LocationMVclass[d], i, zslice, x, y);

              }

            // The cost of pixel p and label l is stored at dsiArray[p*nLabels+l]
            if (genparam==0)
              (dsiArrayV[j])[dsiIndex++] = (float)dsiValue;
            else if (genparam==1)
              (dsiArrayVGen[j])[dsiIndex++] = (float)dsiValue;

            if (dsiValue < bestval) {
              bestval = dsiValue;
              bestd[0] = d;
              numbest = 1;
            } else if (dsiValue == bestval) {
              bestd[numbest] = d;
              numbest++;
            }

          }

          if (numbest == 0)
        	  numbest = nD;

          int curr = rand() % numbest;
          WTArow[x] = bestd[curr];

          delete [] bestd;
        }
      }
    }
  }



  if (klr==1)
    delete ivm;


}









// create data cost using parameters thetaU
DataCost *computeDataCost(int generative, std::vector<CImage2> set1Im, std::vector <CByteImage> hogIm, int numTrainingPats, std::vector<float> thetaU, std::vector<float> thetaA, std::vector<float> thetaH, std::vector<float> thetaL, int featureCode, int featureCodeKlr, std::vector <std::vector <CByteImage> >  appImage, Appearance** appclass, int index, LocationMV** LocationMVclass, AppearanceMV** AppMVclass, std::vector <std::vector  <matrixB<double> > > appdirMVProb, intClass** intObj, std::vector <std::vector <CImage2> > locdirImage, std::vector <std::vector  <matrixB<double> > > hogdirMVProb, intensityNB* intNBclass, matrix<DBL_TYPE> wparam,
                          matrix<DBL_TYPE> Xtrain,
                          matrix<DBL_TYPE> minmaxM, DBL_TYPE lambda, DBL_TYPE p1,
                          std::vector<int> startSliceNo, std::vector<int> endSliceNo,
                          int nbinsNB, int hogdimklr, int appdimklr, int locdimklr, int biasklr)
{
  return new DataCost ( computeDataCostArray (generative, set1Im, hogIm, numTrainingPats, thetaU, thetaA, thetaH, thetaL, featureCode, featureCodeKlr, appImage, appclass, index, LocationMVclass, AppMVclass, appdirMVProb, intObj, locdirImage, hogdirMVProb, intNBclass, wparam, Xtrain, minmaxM, lambda, p1, startSliceNo, endSliceNo, nbinsNB, hogdimklr, appdimklr, locdimklr, biasklr));
}

MRF::CostVal* computeDataCostArray(int generative, std::vector<CImage2> im1, std::vector <CByteImage> hogIm, int numTrainingPats, std::vector<float> thetaU, std::vector<float> thetaA, std::vector<float> thetaH, std::vector<float> thetaL, int featureCode, int featureCodeKlr, std::vector <std::vector <CByteImage> >  appImage, Appearance** appclass, int index, LocationMV** LocationMVclass, AppearanceMV** AppMVclass, std::vector <std::vector  <matrixB<double> > >  appdirMVProb, intClass** intObj, std::vector <std::vector <CImage2> > locdirImage, std::vector <std::vector  <matrixB<double> > > hogdirMVProb, intensityNB* intNBclass, 		matrix<DBL_TYPE> wparam,
                                   matrix<DBL_TYPE> Xtrain,
                                   matrix<DBL_TYPE> minmaxM, DBL_TYPE lambda, DBL_TYPE p1,
                                   std::vector<int> startSliceNo, std::vector<int> endSliceNo,
                                   int nbinsNB, int hogdimklr, int appdimklr, int locdimklr, int biasklr)
{

  int intensity = (int)featureCode/10000;
  int hog = (int) (featureCode/1000)%10;
  int app = (int) (featureCode/100)%10;
  int loc = (int) (featureCode/10)%10;

  int nD = globalNdisps;


  int klr = 0;

  int intensityKlr = (int)featureCodeKlr/10000;
  int hogKlr = (int) (featureCodeKlr/1000)%10;
  int appKlr = (int) (featureCodeKlr/100)%10;
  int locKlr = (int) (featureCodeKlr/10)%10;

  if (intensityKlr!=0 || appKlr!=0 || hogKlr!=0 || locKlr!=0)
    klr = 1;


  if (klr==1)
    {
      if (intensityKlr>0)
        intensity=0;
      if (appKlr>0)
        app=0;
      if (locKlr>0)
        loc=0;
      if (hogKlr>0)
        hog=0;
    }

  int dim = 0;
  IVM* ivm;



  if (intensityKlr==2 || hogKlr==2 || appKlr==2 || locKlr==2)
    {
      klr = 1;

      ivm = new IVM(&Xtrain, &wparam, RBF, lambda, nD, p1);

      dim = Xtrain.size2();
    }




  int zslice = 0;

  if (generative==0) {
    if (dsiArrayV.empty() == true)
      throw CError("call initializeDataCost first");
  } else if (generative==1) {
    if (dsiArrayVGen.empty() == true)
      throw CError("call initializeDataCostLoc first");
  }

  // currently used for slice location in a training subpart of volume

  zslice = startSliceNo[index];

  int nU = (int)thetaU.size();
  int nbins = nU/globalNdisps;

  //	double binner = 400.0/nbins;
  int depth = im1.size();
  CShape sh = im1[0].Shape();
  int width = sh.width, height = sh.height;
  //int nB = sh.nBands;

  int nA = thetaA.size();
  int nT = nA/globalNdisps; //this gives centers per class

  int nH = thetaH.size(); // number of HoG parameters
  int nHpC = nH/globalNdisps; // HoG vocabulary size

  int nL = thetaL.size();
  int nLpC = nL/globalNdisps;

  int dsiIndex = 0;



  for (int z=0; z < depth; z++) {


    matrix<DBL_TYPE> hogdescs;

    if (hogKlr==2)
      {
        // read file
        char hogfile[100];
        sprintf(hogfile, "/p/imageseg/medicalSeg/parameters/hog/hog2dcsv_90n/hog2ddata_%d.csv", z+zslice );  //change this

        // read this file in matrix

        if(!load_one_data(hogfile, hogdescs)){
          std::cout << "csv read error: Bad data file filename : "<< hogfile << endl;
          exit(1);
        }

      }


    matrix<DBL_TYPE> appdescs;

    if (appKlr==2)
      {
        // read file
        char appfile[100];

        sprintf(appfile, "/p/imageseg/medicalSeg/parameters/app/appcsv_81n/app2ddatan_%d.csv", z+zslice );

        // read this file in matrix

        if(!load_one_data(appfile, appdescs)){
          std::cout << "csv read error: Bad data file filename : "<< appfile << endl;
          exit(1);
        }
      }




    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {

        unsigned short pix1 = im1[z].Pixel(x, y, 0);

        ublas::vector<DBL_TYPE> f;
        int thetaid;

        uchar *hogval;
        if (nH > 0)
          hogval = &hogIm[z].Pixel(x, y, 0);


        if (klr==1)
          {
            ublas::vector<DBL_TYPE> Xtest(dim);

            int idxx = 0;

            if (intensityKlr==2) // for now only intensity case
              {
                getCostIntensityKlr(Xtest, (double)pix1, minmaxM, idxx);
              }

            if (locKlr==2)
              {
                getCostLocationKlr(Xtest, z, zslice, x, y, nD, minmaxM, LocationMVclass, idxx, skipXYZlocklr);
              }

            if (appKlr==2)
              {
                getCostAppKlr(Xtest, x, y, width, minmaxM, idxx, appdimklr, &appdescs);
              }


            if (hogKlr==2)
              {
                getCostHoGKlr(Xtest, x, y, width, minmaxM, idxx, hogdimklr, &hogdescs);
              }


            f = ivm->classify_klrexp(&Xtest);
          }


        for (int d = 0; d < globalNdisps; d++) {

          MRF::CostVal dsiValue = 0;

          if (klr==1) {
            dsiValue += -f(d);
            dsiValue += biasklr;
          }


          if (nU>0 && intensity==1)
            {
              dsiValue += getCostIntensityDiscBins(pix1, nbins, d, thetaU, thetaid);

            } else if (intensity==3 && generative==1) // generative setting
            {
              dsiValue += getCostIntensityGenGaussian(intObj[d], pix1);

            } else if (intensity == 6 && generative==1) // using NBayes and bins
            {
              dsiValue += getCostIntensityGenNB(pix1, nbins, d, intNBclass);
            }



          if (nA>0 && app==1) // we have appearance features
            {
              // different features for different classes
              dsiValue += getCostAppDiscPatch((int)appImage[z][d].Pixel(x, y, 0), appclass[d]->getPatchSize(), x, y, width, height, thetaA, d, nT, thetaid);

            } else if (nA>0 && app==2) // we have appearance features
            {
              // common features for all classes
              dsiValue += getCostAppDiscPatch((int)appImage[z][0].Pixel(x, y, 0), appclass[0]->getPatchSize(), x, y, width, height, thetaA, d, nT, thetaid);
            }
          else if (app==3  && generative==1) // generative setting
            {
              dsiValue += getCostAppGenPatch((appdirMVProb[z][d])(y,x), AppMVclass[d]->getDimension(), x, y, width, height);
            }



          if (nH > 0 && hog==1)
            {
              dsiValue += getCostHoGDiscBins((int) *hogval, d, nHpC, thetaH, thetaid);

            } else if (hog==3  && generative==1) // generative setting
            {

              int startPatch = 8; // index starting from 0
              // this variable indicates where the Hog descriptors are defined
              // since i have generated them using the entire patient body, don't need to
              // worry about the z axis and can use all slices here.
              //WORK-THIS need to change above stuff for flawless working
              dsiValue += getCostHoGGenBins((hogdirMVProb[z][d])(y,x), x, y, width, height, startPatch);

            }


          if (nL>0 && loc==1) { //WORK-THIS check boundary conditions

            // this function is completely outdated
            // and for older series of data, will need to modify significantly for use
            dsiValue += getCostLocationDiscCube((unsigned short) locdirImage[d][z].Pixel(x,y,0), x, y, z, loccubex, loccubey, loccubez, numTrainingPats, nLpC, thetaL, width, height, depth, d, thetaid);

          } else if (nL>0 && (loc==6 || loc==8)) // CRF hard or soft
            {

              dsiValue += getCostLocationDiscGaussian(LocationMVclass, z, zslice, x, y, thetaL, d, loc, generative, thetaid);

            } else if (nL>0 && (loc==7 || loc==9) )  //CRF soft or hard assignment redefined
            {

              dsiValue += getCostLocationDiscGaussian2(LocationMVclass, z, zslice, x, y, thetaL, d, loc, generative, nD, thetaid);

            }
          else if (loc==3 && generative==1) // generative setting
            {

              dsiValue += getCostLocationGenGaussian(LocationMVclass[d], z, zslice, x, y);

            }


          if (generative==0) {
            (dsiArrayV[index])[dsiIndex++] = (float)dsiValue;
          } else if (generative==1) {
            (dsiArrayVGen[index])[dsiIndex++] = (float)dsiValue;
          }

        }
      }
    }
  }
  if (generative==0)
    return dsiArrayV[index];
  else if (generative==1)
    return dsiArrayVGen[index];


  if (klr==1)
    delete ivm;

}



DataCost *computeDataCost(int generative, std::vector<CImage2> set1Im, std::vector <CByteImage> hogIm, int numTrainingPats, std::vector <std::vector <CByteImage> >  appImage, Appearance** appclass, int index, LocationMV** LocationMVclass, AppearanceMV** AppMVclass, std::vector <std::vector  <matrixB<double> > > appdirMVProb, intClass** intObj, std::vector <std::vector <CImage2> > locdirImage, std::vector <std::vector  <matrixB<double> > > hogdirMVProb, intensityNB* intNBclass, matrix<DBL_TYPE> wparam,
                          matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM,
                          std::vector<int> startSliceNo, std::vector<int> endSliceNo, Parameters prm)
{
  return new DataCost ( computeDataCostArray (generative, set1Im, hogIm, numTrainingPats, appImage, appclass, index, LocationMVclass, AppMVclass, appdirMVProb, intObj, locdirImage, hogdirMVProb, intNBclass, wparam, Xtrain, minmaxM, startSliceNo, endSliceNo, prm));
}

MRF::CostVal* computeDataCostArray(int generative, std::vector<CImage2> im1, std::vector <CByteImage> hogIm, int numTrainingPats, std::vector <std::vector <CByteImage> >  appImage, Appearance** appclass, int index, LocationMV** LocationMVclass, AppearanceMV** AppMVclass, std::vector <std::vector  <matrixB<double> > >  appdirMVProb, intClass** intObj, std::vector <std::vector <CImage2> > locdirImage, std::vector <std::vector  <matrixB<double> > > hogdirMVProb, intensityNB* intNBclass,
								 matrix<DBL_TYPE> wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM,
                                   std::vector<int> startSliceNo, std::vector<int> endSliceNo, Parameters prm)
{

  std::vector<float> thetaU, thetaA, thetaH, thetaL;
	    // copy thetas
  vecCopy(prm.featurepv.thetaU, thetaU);
  vecCopy(prm.featurepv.thetaA, thetaA);
  vecCopy(prm.featurepv.thetaH, thetaH);
  vecCopy(prm.featurepv.thetaL, thetaL);

  int featureCode = prm.featurepv.featureCode;
  int featureCodeKlr = prm.featurepv.featureCodeKlr;

  DBL_TYPE lambda = prm.featurepv.klrpv.lambda;
  DBL_TYPE p1 = prm.featurepv.klrpv.p1;

  int nbinsNB = prm.featurepv.nbinsintNB;
  int hogdimklr = prm.featurepv.klrpv.hogdimklr;
  int appdimklr = prm.featurepv.klrpv.appdimklr;
  int locdimklr = prm.featurepv.klrpv.locdimklr;
  int biasklr = prm.featurepv.klrpv.biasklr;



  int intensity = (int)featureCode/10000;
  int hog = (int) (featureCode/1000)%10;
  int app = (int) (featureCode/100)%10;
  int loc = (int) (featureCode/10)%10;

  int nD = globalNdisps;


  int klr = 0;

  int intensityKlr = (int)featureCodeKlr/10000;
  int hogKlr = (int) (featureCodeKlr/1000)%10;
  int appKlr = (int) (featureCodeKlr/100)%10;
  int locKlr = (int) (featureCodeKlr/10)%10;

  if (intensityKlr!=0 || appKlr!=0 || hogKlr!=0 || locKlr!=0)
    klr = 1;


  if (klr==1)
    {
      if (intensityKlr>0)
        intensity=0;
      if (appKlr>0)
        app=0;
      if (locKlr>0)
        loc=0;
      if (hogKlr>0)
        hog=0;
    }

  int dim = 0;
  IVM* ivm;



  if (intensityKlr==2 || hogKlr==2 || appKlr==2 || locKlr==2)
    {
      klr = 1;

      ivm = new IVM(&Xtrain, &wparam, RBF, lambda, nD, p1);

      dim = Xtrain.size2();
    }




  int zslice = 0;

  if (generative==0) {
    if (dsiArrayV.empty() == true)
      throw CError("call initializeDataCost first");
  } else if (generative==1) {
    if (dsiArrayVGen.empty() == true)
      throw CError("call initializeDataCostLoc first");
  }

  // currently used for slice location in a training subpart of volume

  zslice = startSliceNo[index];

  int nU = (int)thetaU.size();
  int nbins = nU/globalNdisps;

  //	double binner = 400.0/nbins;
  int depth = im1.size();
  CShape sh = im1[0].Shape();
  int width = sh.width, height = sh.height;
  //int nB = sh.nBands;

  int nA = thetaA.size();
  int nT = nA/globalNdisps; //this gives centers per class

  int nH = thetaH.size(); // number of HoG parameters
  int nHpC = nH/globalNdisps; // HoG vocabulary size

  int nL = thetaL.size();
  int nLpC = nL/globalNdisps;

  int dsiIndex = 0;



  for (int z=0; z < depth; z++) {


    matrix<DBL_TYPE> hogdescs;

    if (hogKlr==2)
      {
        // read file
        char hogfile[100];
        sprintf(hogfile, "/p/imageseg/medicalSeg/parameters/hog/hog2dcsv_90n/hog2ddata_%d.csv", z+zslice );  //change this

        // read this file in matrix

        if(!load_one_data(hogfile, hogdescs)){
          std::cout << "csv read error: Bad data file filename : "<< hogfile << endl;
          exit(1);
        }

      }


    matrix<DBL_TYPE> appdescs;

    if (appKlr==2)
      {
        // read file
        char appfile[100];

        sprintf(appfile, "/p/imageseg/medicalSeg/parameters/app/appcsv_81n/app2ddatan_%d.csv", z+zslice );

        // read this file in matrix

        if(!load_one_data(appfile, appdescs)){
          std::cout << "csv read error: Bad data file filename : "<< appfile << endl;
          exit(1);
        }
      }




    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {

        unsigned short pix1 = im1[z].Pixel(x, y, 0);

        ublas::vector<DBL_TYPE> f;
        int thetaid;

        uchar *hogval;
        if (nH > 0)
          hogval = &hogIm[z].Pixel(x, y, 0);


        if (klr==1)
          {
            ublas::vector<DBL_TYPE> Xtest(dim);

            int idxx = 0;

            if (intensityKlr==2) // for now only intensity case
              {
                getCostIntensityKlr(Xtest, (double)pix1, minmaxM, idxx);
              }

            if (locKlr==2)
              {
                getCostLocationKlr(Xtest, z, zslice, x, y, nD, minmaxM, LocationMVclass, idxx, skipXYZlocklr);
              }

            if (appKlr==2)
              {
                getCostAppKlr(Xtest, x, y, width, minmaxM, idxx, appdimklr, &appdescs);
              }


            if (hogKlr==2)
              {
                getCostHoGKlr(Xtest, x, y, width, minmaxM, idxx, hogdimklr, &hogdescs);
              }


            f = ivm->classify_klrexp(&Xtest);
          }


        for (int d = 0; d < globalNdisps; d++) {

          MRF::CostVal dsiValue = 0;

          if (klr==1) {
            dsiValue += -f(d);
            dsiValue += biasklr;
          }


          if (nU>0 && intensity==1)
            {
              dsiValue += getCostIntensityDiscBins(pix1, nbins, d, thetaU, thetaid);

            } else if (intensity==3 && generative==1) // generative setting
            {
              dsiValue += getCostIntensityGenGaussian(intObj[d], pix1);

            } else if (intensity == 6 && generative==1) // using NBayes and bins
            {
              dsiValue += getCostIntensityGenNB(pix1, nbinsNB, d, intNBclass);
            }



          if (nA>0 && app==1) // we have appearance features
            {
              // different features for different classes
              dsiValue += getCostAppDiscPatch((int)appImage[z][d].Pixel(x, y, 0), appclass[d]->getPatchSize(), x, y, width, height, thetaA, d, nT, thetaid);

            } else if (nA>0 && app==2) // we have appearance features
            {
              // common features for all classes
              dsiValue += getCostAppDiscPatch((int)appImage[z][0].Pixel(x, y, 0), appclass[0]->getPatchSize(), x, y, width, height, thetaA, d, nT, thetaid);
            }
          else if (app==3  && generative==1) // generative setting
            {
              dsiValue += getCostAppGenPatch((appdirMVProb[z][d])(y,x), AppMVclass[d]->getDimension(), x, y, width, height);
            }



          if (nH > 0 && hog==1)
            {
              dsiValue += getCostHoGDiscBins((int) *hogval, d, nHpC, thetaH, thetaid);

            } else if (hog==3  && generative==1) // generative setting
            {

              int startPatch = 8; // index starting from 0
              // this variable indicates where the Hog descriptors are defined
              // since i have generated them using the entire patient body, don't need to
              // worry about the z axis and can use all slices here.
              //WORK-THIS need to change above stuff for flawless working
              dsiValue += getCostHoGGenBins((hogdirMVProb[z][d])(y,x), x, y, width, height, startPatch);

            }


          if (nL>0 && loc==1) { //WORK-THIS check boundary conditions

            // this function is completely outdated
            // and for older series of data, will need to modify significantly for use
            dsiValue += getCostLocationDiscCube((unsigned short) locdirImage[d][z].Pixel(x,y,0), x, y, z, loccubex, loccubey, loccubez, numTrainingPats, nLpC, thetaL, width, height, depth, d, thetaid);

          } else if (nL>0 && (loc==6 || loc==8)) // CRF hard or soft
            {

              dsiValue += getCostLocationDiscGaussian(LocationMVclass, z, zslice, x, y, thetaL, d, loc, generative, thetaid);

            } else if (nL>0 && (loc==7 || loc==9) )  //CRF soft or hard assignment redefined
            {

              dsiValue += getCostLocationDiscGaussian2(LocationMVclass, z, zslice, x, y, thetaL, d, loc, generative, nD, thetaid);

            }
          else if (loc==3 && generative==1) // generative setting
            {

              dsiValue += getCostLocationGenGaussian(LocationMVclass[d], z, zslice, x, y);

            }


          if (generative==0) {
            (dsiArrayV[index])[dsiIndex++] = (float)dsiValue;
          } else if (generative==1) {
            (dsiArrayVGen[index])[dsiIndex++] = (float)dsiValue;
          }

        }
      }
    }
  }
  if (generative==0)
    return dsiArrayV[index];
  else if (generative==1)
    return dsiArrayVGen[index];


  if (klr==1)
    delete ivm;

}



DataCost *computeDataCost(int generative, std::vector<CByteImage> set1Im, std::vector <CByteImage> hogIm, std::vector <CByteImage> motIm, int numTrainingPats, std::vector <std::vector <CByteImage> >  appImage, Appearance** appclass, int index, LocationMV** LocationMVclass, AppearanceMV** AppMVclass, std::vector <std::vector  <matrixB<double> > > appdirMVProb, intClass** intObj, std::vector <std::vector <CImage2> > locdirImage, std::vector <std::vector  <matrixB<double> > > hogdirMVProb, intensityNB* intNBclass, matrix<DBL_TYPE> wparam,
                          matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM,
                          std::vector<int> startSliceNo, std::vector<int> endSliceNo, Parameters prm)
{
  return new DataCost ( computeDataCostArray (generative, set1Im, hogIm, motIm, numTrainingPats, appImage, appclass, index, LocationMVclass, AppMVclass, appdirMVProb, intObj, locdirImage, hogdirMVProb, intNBclass, wparam, Xtrain, minmaxM, startSliceNo, endSliceNo, prm));
}

MRF::CostVal* computeDataCostArray(int generative, std::vector<CByteImage> im1, std::vector <CByteImage> hogIm, std::vector <CByteImage> motIm, int numTrainingPats, std::vector <std::vector <CByteImage> >  appImage, Appearance** appclass, int index, LocationMV** LocationMVclass, AppearanceMV** AppMVclass, std::vector <std::vector  <matrixB<double> > >  appdirMVProb, intClass** intObj, std::vector <std::vector <CImage2> > locdirImage, std::vector <std::vector  <matrixB<double> > > hogdirMVProb, intensityNB* intNBclass,
								 matrix<DBL_TYPE> wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM,
                                   std::vector<int> startSliceNo, std::vector<int> endSliceNo, Parameters prm)
{

  std::vector<float> thetaU, thetaA, thetaH, thetaL, thetaM;
	    // copy thetas
  vecCopy(prm.featurepv.thetaU, thetaU);
  vecCopy(prm.featurepv.thetaA, thetaA);
  vecCopy(prm.featurepv.thetaH, thetaH);
  vecCopy(prm.featurepv.thetaM, thetaM);
  vecCopy(prm.featurepv.thetaL, thetaL);

  int featureCode = prm.featurepv.featureCode;
  int featureCodeKlr = prm.featurepv.featureCodeKlr;

  DBL_TYPE lambda = prm.featurepv.klrpv.lambda;
  DBL_TYPE p1 = prm.featurepv.klrpv.p1;

  int nbinsNB = prm.featurepv.nbinsintNB;
  int hogdimklr = prm.featurepv.klrpv.hogdimklr;
  int appdimklr = prm.featurepv.klrpv.appdimklr;
  int locdimklr = prm.featurepv.klrpv.locdimklr;
  int biasklr = prm.featurepv.klrpv.biasklr;



  int intensity = (int)featureCode/1000000;
  int hog = (int) (featureCode/100000)%10;
  int app = (int) (featureCode/10000)%10;
  int loc = (int) (featureCode/1000)%10;
  int mot = (int) (featureCode/100)%10;
  int rbm = (int) (featureCode/10)%10;

  int nD = globalNdisps;


  int klr = 0;

  int intensityKlr = (int)featureCodeKlr/10000;
  int hogKlr = (int) (featureCodeKlr/1000)%10;
  int appKlr = (int) (featureCodeKlr/100)%10;
  int locKlr = (int) (featureCodeKlr/10)%10;

  if (intensityKlr!=0 || appKlr!=0 || hogKlr!=0 || locKlr!=0)
    klr = 1;


  if (klr==1)
    {
      if (intensityKlr>0)
        intensity=0;
      if (appKlr>0)
        app=0;
      if (locKlr>0)
        loc=0;
      if (hogKlr>0)
        hog=0;
    }

  int dim = 0;
  IVM* ivm;



  if (intensityKlr==2 || hogKlr==2 || appKlr==2 || locKlr==2)
    {
      klr = 1;

      ivm = new IVM(&Xtrain, &wparam, RBF, lambda, nD, p1);

      dim = Xtrain.size2();
    }




  int zslice = 0;

  if (generative==0) {
    if (dsiArrayV.empty() == true)
      throw CError("call initializeDataCost first");
  } else if (generative==1) {
    if (dsiArrayVGen.empty() == true)
      throw CError("call initializeDataCostLoc first");
  }

  // currently used for slice location in a training subpart of volume

  zslice = startSliceNo[index];

  int nU = (int)thetaU.size();
  int nbins = nU/globalNdisps;

  //	double binner = 400.0/nbins;
  int depth = im1.size();
  CShape sh = im1[0].Shape();
  int width = sh.width, height = sh.height;
  //int nB = sh.nBands;

  int nA = thetaA.size();
  int nT = nA/globalNdisps; //this gives centers per class

  int nH = thetaH.size(); // number of HoG parameters
  int nHpC = nH/globalNdisps; // HoG vocabulary size

  int nM = thetaM.size(); // number of mot parameters
  int nMpC = nM/globalNdisps; // mot vocabulary size

  int nL = thetaL.size();
  int nLpC = nL/globalNdisps;

  int dsiIndex = 0;



  for (int z=0; z < depth; z++) {


    matrix<DBL_TYPE> hogdescs;

    if (hogKlr==2)
      {
        // read file
        char hogfile[100];
        sprintf(hogfile, "/p/imageseg/medicalSeg/parameters/hog/hog2dcsv_90n/hog2ddata_%d.csv", z+zslice );  //change this

        // read this file in matrix

        if(!load_one_data(hogfile, hogdescs)){
          std::cout << "csv read error: Bad data file filename : "<< hogfile << endl;
          exit(1);
        }

      }


    matrix<DBL_TYPE> appdescs;

    if (appKlr==2)
      {
        // read file
        char appfile[100];

        sprintf(appfile, "/p/imageseg/medicalSeg/parameters/app/appcsv_81n/app2ddatan_%d.csv", z+zslice );

        // read this file in matrix

        if(!load_one_data(appfile, appdescs)){
          std::cout << "csv read error: Bad data file filename : "<< appfile << endl;
          exit(1);
        }
      }




    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {

        unsigned char *pix1 = &im1[z].Pixel(x, y, 0);

        ublas::vector<DBL_TYPE> f;
        int thetaid;
        int thetaidrgb[3]={0,0,0};
        int thetaidYuv[2]={0, 0};

        uchar *hogval;
        if (nH > 0)
          hogval = &hogIm[z].Pixel(x, y, 0);

        uchar *motval;
        if (nM > 0)
          motval = &motIm[z].Pixel(x, y, 0);


        if (klr==1)
          {
            ublas::vector<DBL_TYPE> Xtest(dim);

            int idxx = 0;

            if (intensityKlr==2) // for now only intensity case
              {
                //getCostIntensityKlr(Xtest, (double)pix1, minmaxM, idxx);
              }

            if (locKlr==2)
              {
                getCostLocationKlr(Xtest, z, zslice, x, y, nD, minmaxM, LocationMVclass, idxx, skipXYZlocklr);
              }

            if (appKlr==2)
              {
                getCostAppKlr(Xtest, x, y, width, minmaxM, idxx, appdimklr, &appdescs);
              }


            if (hogKlr==2)
              {
                getCostHoGKlr(Xtest, x, y, width, minmaxM, idxx, hogdimklr, &hogdescs);
              }


            f = ivm->classify_klrexp(&Xtest);
          }


        for (int d = 0; d < globalNdisps; d++) {

          MRF::CostVal dsiValue = 0;

          if (klr==1) {
            dsiValue += -f(d);
            dsiValue += biasklr;
          }


          if (nU>0 && intensity==1)
            {
              // dsiValue += getCostIntensityDiscBins(pix1, nbins, d, thetaU, thetaid);
              // dsiValue += getCostIntensityDiscBins(pix1, nbins, d, thetaU, thetaidrgb);
              dsiValue += getCostIntensityDiscBinsYuv(pix1, nbins, d, thetaU, thetaidYuv);

            } else if (intensity==3 && generative==1) // generative setting
            {
              //dsiValue += getCostIntensityGenGaussian(intObj[d], pix1);

            } else if (intensity == 6 && generative==1) // using NBayes and bins
            {
              //dsiValue += getCostIntensityGenNB(pix1, nbins, d, intNBclass);
            }


          if (nA>0 && app==1) // we have appearance features
            {
              // different features for different classes
        	  if ((int)appImage[z][d].Pixel(x, y, 0) >= nT)
        	  {
        		  // this is to handle the possibility that some of the pixels in the image borders might have invalid values
        	  } else
        	  {
        		  dsiValue += getCostAppDiscPatch((int)appImage[z][d].Pixel(x, y, 0), appclass[d]->getPatchSize(), x, y, width, height, thetaA, d, nT, thetaid);
        	  }
            } else if (nA>0 && app==2) // we have appearance features
            {
              // common features for all classes
            	if ((int)appImage[z][0].Pixel(x, y, 0) >= nA)
            	{
            		  // this is to handle the possibility that some of the pixels in the image borders might have invalid values
            	} else
            	{
            		dsiValue += getCostAppDiscPatch((int)appImage[z][0].Pixel(x, y, 0), appclass[0]->getPatchSize(), x, y, width, height, thetaA, d, nT, thetaid);
            	}
            }
          else if (app==3  && generative==1) // generative setting
            {
              dsiValue += getCostAppGenPatch((appdirMVProb[z][d])(y,x), AppMVclass[d]->getDimension(), x, y, width, height);
            }



          if (nH > 0 && hog==1)
            {
              dsiValue += getCostHoGDiscBins((int) *hogval, d, nHpC, thetaH, thetaid);

            } else if (hog==3  && generative==1) // generative setting
            {

              int startPatch = 8; // index starting from 0
              // this variable indicates where the Hog descriptors are defined
              // since i have generated them using the entire patient body, don't need to
              // worry about the z axis and can use all slices here.
              //WORK-THIS need to change above stuff for flawless working
              dsiValue += getCostHoGGenBins((hogdirMVProb[z][d])(y,x), x, y, width, height, startPatch);

            }



			if (nM > 0 && mot==1)
			{
			  dsiValue += getCostmotDiscBins((int) *motval, d, nMpC, thetaM, thetaid);

			} else if (mot==3  && generative==1) // generative setting
			{
				std::cout << " mot generative not implemented! \n";
				exit(-1);
			}


          if (nL>0 && loc==1) { //WORK-THIS check boundary conditions

            // this function is completely outdated
            // and for older series of data, will need to modify significantly for use
            dsiValue += getCostLocationDiscCube((unsigned short) locdirImage[d][z].Pixel(x,y,0), x, y, z, loccubex, loccubey, loccubez, numTrainingPats, nLpC, thetaL, width, height, depth, d, thetaid);

          } else if (nL>0 && (loc==6 || loc==8)) // CRF hard or soft
            {

              dsiValue += getCostLocationDiscGaussian(LocationMVclass, z, zslice, x, y, thetaL, d, loc, generative, thetaid);

            } else if (nL>0 && (loc==7 || loc==9) )  //CRF soft or hard assignment redefined
            {

              dsiValue += getCostLocationDiscGaussian2(LocationMVclass, z, zslice, x, y, thetaL, d, loc, generative, nD, thetaid);

            }
          else if (loc==3 && generative==1) // generative setting
            {

              dsiValue += getCostLocationGenGaussian(LocationMVclass[d], z, zslice, x, y);

            }


          if (generative==0) {
            (dsiArrayV[index])[dsiIndex++] = (float)dsiValue;
          } else if (generative==1) {
            (dsiArrayVGen[index])[dsiIndex++] = (float)dsiValue;
          }

        }
      }
    }
  }
  if (generative==0)
    return dsiArrayV[index];
  else if (generative==1)
    return dsiArrayVGen[index];


  if (klr==1)
    delete ivm;

}









//
void getCostArray(float* &costs, int &nPixels, int &nDisps)
{
  costs = (dsiArray2.empty() == true)? dsiArrayV[cueIndex] : dsiArray2[cueIndex];
  nPixels = globalNpixels[cueIndex];
  nDisps = globalNdisps;
}

// filter definitions
enum FTLR { esimple2, esimple4, esimple2sq, esimple4sq, eprewitt3, eprewitt4, esobel3, esobel5, ewilson3, ewilson4 };
// simple2
double simple2[2] = {-1, +1};
// simple4
double simple4[4] = {-1, -1, +1, +1};
// simple2sq
double simple2sq[2] = {-1, +1};
// simple4sq
double simple4sq[4] = {-1, -1, +1, +1};
// prewitt3x3
double prewitt3x3[3][3] = {{-1, 0, +1}, {-1, 0, +1}, {-1, 0, +1}};
double prewitt3y3[3][3] = {{-1, -1, -1}, {0, 0, 0}, {+1, +1, +1}};
// prewitt5x5
double prewitt4x4[4][4] = {{-3, -1, +1, +3}, {-3, -1, +1, +3}, {-3, -1, +1, +3}, {-3, -1, +1, +3}};
double prewitt4y4[4][4] = {{-3, -3, -3, -3}, {-1, -1, -1, -1}, {+1, +1, +1, +1}, {+3, +3, +3, +3}};
// sobel3x3
double sobel3x3[3][3] = {{-1, 0, +1}, {-2, 0, +2}, {-1, 0, +1}};
double sobel3y3[3][3] = {{-1, -2, -1}, {0, 0, 0}, {+1, +2, +1}};
// sobel5x5
double sobel5x5[5][5] = {{-1, -2, 0, +2, +1}, {-4, -8, 0, +8, +4}, {-6, -12, 0 , +12, +6}, {-4, -8, 0, +8, +4}, {-1, -2, 0, +2, +1}};
double sobel5y5[5][5] = {{-1, -4, -6, -4, -1}, {-2, -8, -12, -8, -2}, {0, 0, 0 , 0, 0}, {+2, +8, +12, +8, +2}, {+1, +4, +6, +4, +1}};
// canny9x9
//?
// wilson3x3
double wilson3x3[3][3] = {{-0.0741, -0.0955, 0}, {-0.0955, 0, 0.0955}, {0, 0.0955, 0.0741}};
double wilson3y3[3][3] = {{0, -0.0955, -0.0741}, {+0.0955, 0, -0.0955}, {0.0741, 0.0955, 0}};
// wilson4x4
double wilson4x4[4][4] = {{-0.0107, -0.0496, -0.0277, 0}, {-0.0496, -0.1292, 0, 0.0277}, {-0.0277, 0, 0.1292, 0.0496}, {0, 0.0277, 0.0496, 0.0107}};
double wilson4y4[4][4] = {{0, -0.0277,-0.0496, -0.0107}, {0.0277, 0, -0.1292, -0.0496}, {0.0496, 0.1292, 0, -0.0277}, {0.0107, 0.0496, 0.0277, 0}};



void computeQuantizedGradients(std::vector<CImage2> im1, std::vector<CByteImage> &im1grad, std::vector<int> gradThresh)
  /*
    use nK = gradThresh.size() levels of quantization.
    quantization is done as follows:
    compute absolute gradient (RMS difference of color bands) in x and y
    given gradient g, assign largest k such that gradThresh[k] <= g
    the last entry (gradThresh[nK-1]) is ignored and is assumed to be infinity
    in im1grad, use band 0 for x gradient, and band 1 for y gradients
  */
{
  int depth = im1.size();
  CShape sh = im1[0].Shape();
  int width = sh.width, height = sh.height, nB = sh.nBands;
  int nColors = __min(3, nB);

  enum FTLR filter = esimple2;

  int nK = (int)gradThresh.size();
  nK = nK + 1; // adding a 1 because will push maxlimit (inf) in gradThresh
  int k;

  DEBUG_OUT1(verbose, debugfile, "%d Gradient bins:\t0-", nK);
  for (k = 0; k < nK-1; k++)
	DEBUG_OUT2(verbose, debugfile, "%d\t%d-", gradThresh[k]-1, gradThresh[k]);
  DEBUG_OUT0(verbose, debugfile, "inf\n");

  sh.nBands = 3;  // band 0: x gradient, band 1: y gradient, band 2: z gradient

  for (int z = 0; z<depth; z++) {
	  im1grad[z].ReAllocate(sh);
  }


  // maxlimit depends on the filter.
  // since i am using the range - 0 to 400, using 450 as the maximal difference between two pixels for the simple filters

  // gradThresh.push_back(450 * 450 * nColors); // this actually depends on the filter type and size

  // the constant factor is the sum of the +ve numbers in the filter matrix
  double upbound = 0, fctor=0;

  if (filter == esimple2) {
	  gradThresh.push_back(450 * nColors);
	  fctor = 1;
  } else if (filter == esimple4) {
	  gradThresh.push_back((450*2) * nColors);  // the 2 is because there are two +1 and two -1
	  fctor = 2;
  } else if (filter == esimple2sq) {
	  gradThresh.push_back(450 * 450 * nColors);
	  fctor = 1;
  } else if (filter == esimple4sq) {
	  gradThresh.push_back(450 * 450 * 4 * nColors);
	  fctor = 4;
  } else if (filter == eprewitt3) {
	  gradThresh.push_back((450*3) * nColors);
	  fctor = 3;
  } else if (filter == eprewitt4) {
	  gradThresh.push_back((450*(3+1)*4) * nColors);
	  fctor = 16;
  } else if (filter == esobel3) {
	  gradThresh.push_back((450*4) * nColors);
	  fctor = 4;
  } else if (filter == esobel5) {
	  gradThresh.push_back((450*48) * nColors);
	  fctor = 48;
  } else if (filter == ewilson3) {
	  // gradThresh.push_back((450*0.2651) * nColors);
	  gradThresh.push_back((450) * nColors);
	  fctor = 1;
  } else if (filter == ewilson4) {
	  // gradThresh.push_back((450*0.2945) * nColors);
	  gradThresh.push_back((450) * nColors);
	  fctor = 1;
  }
  // need to fill in for the other filters

  upbound = 450 * fctor;


  if (filter==esimple2sq || filter==esimple4sq) {
	  for (k = 0; k < nK-1; k++) // compute sum of squared color differences below,
		  gradThresh[k] *= gradThresh[k] * fctor * nColors;  //  so need to adjust thresholds to get RMS distance
  } else if (filter==esimple2 || filter==esimple4 || filter==eprewitt3 || filter==esobel3 || filter==ewilson3 || filter==eprewitt4 || filter==ewilson4 || filter==esobel5) {
	  for (k = 0; k < nK-1; k++) // compute sum of squared color differences below,
		  gradThresh[k] *= fctor * nColors;  //  so need to adjust thresholds to get RMS distance
  }


  double *filx3, *fily3, *filx4, *fily4, *filx5, *fily5;

  if (filter==eprewitt3) {
	  filx3=&prewitt3x3[0][0];
	  fily3=&prewitt3y3[0][0];
  } else if (filter==esobel3) {
	  filx3=&sobel3x3[0][0];
	  fily3=&sobel3y3[0][0];
  } else if (filter==ewilson3) {
	  filx3=&wilson3x3[0][0];
	  fily3=&wilson3y3[0][0];
  } else if (filter==eprewitt4) {
	  filx4=&prewitt4x4[0][0];
	  fily4=&prewitt4y4[0][0];
  } else if (filter==ewilson4) {
	  filx4=&wilson4x4[0][0];
	  fily4=&wilson4y4[0][0];
  } else if (filter==esobel5) {
	  filx5=&sobel5x5[0][0];
	  fily5=&sobel5y5[0][0];
  }



// need to decide whether we use the 2D filter or 3D filter
// for now i only have 2D filters
  // and so no restriction on the depth for loop
  for (int z = 0; z < depth; z++) {
	  for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {

		  uchar *grad   = &im1grad[z].Pixel(x, y, 0);
		  int sx = 0, sy = 0, sz = 0; // for now sz will not really be used
		  int flagc = 0;

		  if (filter==esimple2 || filter==esimple2sq || filter==eprewitt3 || filter==esobel3 || filter==ewilson3) {
			  if (x==0 || y==0 || x==width-1 || y==height-1) {
				  sx=0; sy=0; sz=0;
				  flagc = 1;
			  }
		  } else if (filter==esimple4 || filter==esimple4sq || filter==eprewitt4 || filter==ewilson4 || filter==esobel5) {
			  if (x==0 || x==1 || y==0 || y==1 || x==width-1 || x==width-2 || y==height-1 || y==height-2) {
				  sx = 0; sy = 0; sz = 0;
				  flagc = 1;
			  }
		  }

		  if (flagc==0) { // which means i still need to compute sx, sy, sz

			  double dx, dy, dz;

			  for (int b = 0; b < nColors; b++) {

				  if (filter==esimple2 || filter==esimple2sq) {
					  dx = simple2[0]*im1[z].Pixel(x-1, y, b) + simple2[1]*im1[z].Pixel(x, y, b);
					  dy = simple2[0]*im1[z].Pixel(x, y-1, b) + simple2[1]*im1[z].Pixel(x, y, b);
					  // need to take care of dz too
					  dz = 0;
				  }

				  if (filter==esimple4 || filter==esimple4sq) {
					  dx = simple4[0]*im1[z].Pixel(x-2, y, b) + simple4[1]*im1[z].Pixel(x-1, y, b) + simple4[2]*im1[z].Pixel(x, y, b) + simple4[3]*im1[z].Pixel(x+1, y, b);
					  dy = simple4[0]*im1[z].Pixel(x, y-2, b) + simple4[1]*im1[z].Pixel(x, y-1, b) + simple4[2]*im1[z].Pixel(x, y, b) + simple4[3]*im1[z].Pixel(x, y+1, b);
					  // need to take care of dz too
					  dz = 0;
				  }

				  if (filter==eprewitt3 || filter==esobel3 || filter==ewilson3) {
					  dx = 0; dy = 0; dz = 0;
					  for (int m=0; m<3; m++) {
						  for (int n=0; n<3; n++) {
							  dx += filx3[3*m+n]*(double)im1[z].Pixel(x+m-1, y+n-1, b);
							  dy += fily3[3*m+n]*(double)im1[z].Pixel(x+m-1, y+n-1, b);
							  dz = 0; // for now
						  }
					  }
				  }

				  if (filter==eprewitt4 || filter==ewilson4) {
					  dx = 0; dy = 0; dz = 0;
					  for (int m=0; m<4; m++) {
						  for (int n=0; n<4; n++) {
							  dx += filx4[4*m+n]*(double)im1[z].Pixel(x+m-2, y+n-2, b);
							  dy += fily4[4*m+n]*(double)im1[z].Pixel(x+m-2, y+n-2, b);
							  dz = 0; // for now
						  }
					  }
				  }

				  if (filter==esobel5) {
					  dx = 0; dy = 0; dz = 0;
					  for (int m=0; m<5; m++) {
						  for (int n=0; n<5; n++) {
							  dx += filx5[5*m+n]*(double)im1[z].Pixel(x+m-2, y+n-2, b);
							  dy += fily5[5*m+n]*(double)im1[z].Pixel(x+m-2, y+n-2, b);
							  dz = 0; // for now
						  }
					  }
				  }

				  if (dx<0)
					  dx=-dx;
				  if (dy<0)
					  dy=-dy;

				  if (dx>upbound)
					  dx = upbound;
				  if (dy>upbound)
					  dy = upbound;


				  if (filter==esimple2 || filter==esimple4 || filter==eprewitt3 || filter==esobel3 || filter==ewilson3 || filter==eprewitt4 || filter==ewilson4 || filter==esobel5) {
					  sx = (int) dx;
					  sy = (int) dy;
					  sz = (int) dz;
				  } else if (filter==esimple2sq || filter==esimple4sq) {
					  sx += (int)(dx * dx);
					  sy += (int)(dy * dy);
					  sz += (int)(dz * dz);
				  }

			  }
		  }
		  for (k = 0; gradThresh[k] <= sx; k++)
			; // find lowest bin for x gradient
		  assert(k < nK);
		  grad[0] = k;

		  for (k = 0; gradThresh[k] <= sy; k++)
			; // find lowest bin for y gradient
		  assert(k < nK);
		  grad[1] = k;

		  for (k = 0; gradThresh[k] <= sz; k++)
			; // find lowest bin for z gradient
		  assert(k < nK);
		  grad[2] = k;

		}
	  }
  }
}



// grayscale (x-y)
void computeQuantizedGradients00(std::vector<CImage2> im1, std::vector<CByteImage> &im1grad, std::vector<int> gradThresh)
  /*
    use nK = gradThresh.size() levels of quantization.
    quantization is done as follows:
    compute absolute gradient (RMS difference of color bands) in x and y
    given gradient g, assign largest k such that gradThresh[k] <= g
    the last entry (gradThresh[nK-1]) is ignored and is assumed to be infinity
    in im1grad, use band 0 for x gradient, and band 1 for y gradients
  */
{
  int depth = im1.size();
  CShape sh = im1[0].Shape();
  int width = sh.width, height = sh.height, nB = sh.nBands;
  int nColors = __min(3, nB);
  int nColorsg = 1;

  enum FTLR filter = esimple2;

  int nK = (int)gradThresh.size();
  nK = nK + 1; // adding a 1 because will push maxlimit (inf) in gradThresh
  int k;

  DEBUG_OUT1(verbose, debugfile, "%d Gradient bins:\t0-", nK);
  for (k = 0; k < nK-1; k++)
	DEBUG_OUT2(verbose, debugfile, "%d\t%d-", gradThresh[k]-1, gradThresh[k]);
  DEBUG_OUT0(verbose, debugfile, "inf\n");

  sh.nBands = 3;  // band 0: x gradient, band 1: y gradient, band 2: z gradient

  for (int z = 0; z<depth; z++) {
	  im1grad[z].ReAllocate(sh);
  }

  // the constant factor is the sum of the +ve numbers in the filter matrix
  double upbound = 0, fctor=0;

  if (filter == esimple2) {
	  gradThresh.push_back(260 * nColorsg);
	  fctor = 1;
  }
  // need to fill in for the other filters

  upbound = 260 * fctor;

// need to decide whether we use the 2D filter or 3D filter
// for now i only have 2D filters
  // and so no restriction on the depth for loop
  for (int z = 0; z < depth; z++) {
	  for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {

		  uchar *grad   = &im1grad[z].Pixel(x, y, 0);
		  int sx = 0, sy = 0, sz = 0;

		  if (filter==esimple2) {

			  int pix1 = im1[z].Pixel(x, y, 0);
			  int pix2 = 0;
			  int pix3 = 0;
			  int pix4 = 0;

			  double dx=0.0, dy=0.0, dz=0.0;

			  if (x==width-1)
				  sx = 0;
			  else {
				  pix2 = im1[z].Pixel(x+1, y, 0);
				  dx = simple2[0]*pix2 + simple2[1]*pix1;
			  }

			  if (y==height-1)
				  sy = 0;
			  else {
				  pix3 = im1[z].Pixel(x, y+1, 0);
				  dy = simple2[0]*pix3 + simple2[1]*pix1;
			  }


			  if (dx<0)
				  dx=-dx;
			  if (dy<0)
				  dy=-dy;

			  if (dx>upbound)
				  dx = upbound;
			  if (dy>upbound)
				  dy = upbound;

			  sx = (int) dx;
			  sy = (int) dy;
			  sz = (int) 0;

		  }

		  for (k = 0; gradThresh[k] <= sx; k++)
			; // find lowest bin for x gradient
		  assert(k < nK);
		  grad[0] = k;

		  for (k = 0; gradThresh[k] <= sy; k++)
			; // find lowest bin for y gradient
		  assert(k < nK);
		  grad[1] = k;

		  grad[2] = 0;

		}
	  }
  }
}





//(x-y-z)
// grayscale (for medical images)
void computeQuantizedGradients0(std::vector<CImage2> im1, std::vector<CByteImage> &im1grad, std::vector<int> gradThresh)
  /*
    use nK = gradThresh.size() levels of quantization.
    quantization is done as follows:
    compute absolute gradient (RMS difference of color bands) in x and y
    given gradient g, assign largest k such that gradThresh[k] <= g
    the last entry (gradThresh[nK-1]) is ignored and is assumed to be infinity
    in im1grad, use band 0 for x gradient, and band 1 for y gradients
  */
{
  int depth = im1.size();
  CShape sh = im1[0].Shape();
  int width = sh.width, height = sh.height, nB = sh.nBands;
  int nColors = __min(3, nB);
  int nColorsg = 1;

  enum FTLR filter = esimple2;

  int nK = (int)gradThresh.size();
  nK = nK + 1; // adding a 1 because will push maxlimit (inf) in gradThresh
  int k;

  DEBUG_OUT1(verbose, debugfile, "%d Gradient bins:\t0-", nK);
  for (k = 0; k < nK-1; k++)
	DEBUG_OUT2(verbose, debugfile, "%d\t%d-", gradThresh[k]-1, gradThresh[k]);
  DEBUG_OUT0(verbose, debugfile, "inf\n");

  sh.nBands = 3;  // band 0: x gradient, band 1: y gradient, band 2: z gradient

  for (int z = 0; z<depth; z++) {
	  im1grad[z].ReAllocate(sh);
  }

  // the constant factor is the sum of the +ve numbers in the filter matrix
  double upbound = 0, fctor=0;

  if (filter == esimple2) {
	  gradThresh.push_back(260 * nColorsg);
	  fctor = 1;
  }
  // need to fill in for the other filters

  upbound = 260 * fctor;

// need to decide whether we use the 2D filter or 3D filter
// for now i only have 2D filters
  // and so no restriction on the depth for loop
  for (int z = 0; z < depth; z++) {
	  for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {

		  uchar *grad   = &im1grad[z].Pixel(x, y, 0);
		  int sx = 0, sy = 0, sz = 0;

		  if (filter==esimple2) {

			  int pix1 = im1[z].Pixel(x, y, 0);
			  int pix2 = 0;
			  int pix3 = 0;
			  int pix4 = 0;

			  double dx=0.0, dy=0.0, dz=0.0;

			  if (x==width-1)
				  sx = 0;
			  else {
				  pix2 = im1[z].Pixel(x+1, y, 0);
				  dx = simple2[0]*pix2 + simple2[1]*pix1;
			  }

			  if (y==height-1)
				  sy = 0;
			  else {
				  pix3 = im1[z].Pixel(x, y+1, 0);
				  dy = simple2[0]*pix3 + simple2[1]*pix1;
			  }

			  if (z==depth-1)
				  sz = 0;
			  else {
				  pix4 = im1[z+1].Pixel(x, y, 0);
				  dz = simple2[0]*pix4 + simple2[1]*pix1;
			  }

			  if (dx<0)
				  dx=-dx;
			  if (dy<0)
				  dy=-dy;
			  if (dz<0)
				  dz=-dz;

			  if (dx>upbound)
				  dx = upbound;
			  if (dy>upbound)
				  dy = upbound;
			  if (dz>upbound)
				  dz = upbound;

			  sx = (int) dx;
			  sy = (int) dy;
			  sz = (int) dz;

		  }

		  for (k = 0; gradThresh[k] <= sx; k++)
			; // find lowest bin for x gradient
		  assert(k < nK);
		  grad[0] = k;

		  for (k = 0; gradThresh[k] <= sy; k++)
			; // find lowest bin for y gradient
		  assert(k < nK);
		  grad[1] = k;

		  for (k = 0; gradThresh[k] <= sz; k++)
			; // find lowest bin for z gradient
		  assert(k < nK);
		  grad[2] = k;

		}
	  }
  }
}



// RGB converted to grayscale
void computeQuantizedGradients(std::vector<CByteImage> im1, std::vector<CByteImage> &im1grad, std::vector<int> gradThresh)
  /*
    use nK = gradThresh.size() levels of quantization.
    quantization is done as follows:
    compute absolute gradient (RMS difference of color bands) in x and y
    given gradient g, assign largest k such that gradThresh[k] <= g
    the last entry (gradThresh[nK-1]) is ignored and is assumed to be infinity
    in im1grad, use band 0 for x gradient, and band 1 for y gradients
  */
{
  int depth = im1.size();
  CShape sh = im1[0].Shape();
  int width = sh.width, height = sh.height, nB = sh.nBands;
  int nColors = __min(3, nB);
  int nColorsg = 1;

  enum FTLR filter = esimple2;

  int nK = (int)gradThresh.size();
  nK = nK + 1; // adding a 1 because will push maxlimit (inf) in gradThresh
  int k;

  DEBUG_OUT1(verbose, debugfile, "%d Gradient bins:\t0-", nK);
  for (k = 0; k < nK-1; k++)
	DEBUG_OUT2(verbose, debugfile, "%d\t%d-", gradThresh[k]-1, gradThresh[k]);
  DEBUG_OUT0(verbose, debugfile, "inf\n");

  sh.nBands = 3;  // band 0: x gradient, band 1: y gradient, band 2: z gradient

  for (int z = 0; z<depth; z++) {
	  im1grad[z].ReAllocate(sh);
  }

  // the constant factor is the sum of the +ve numbers in the filter matrix
  double upbound = 0, fctor=0;

  if (filter == esimple2) {
	  gradThresh.push_back(256 * nColorsg);
	  fctor = 1;
  }
  // need to fill in for the other filters

  upbound = 256 * fctor;

// need to decide whether we use the 2D filter or 3D filter
// for now i only have 2D filters
  // and so no restriction on the depth for loop
  for (int z = 0; z < depth; z++) {
	  for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {

		  uchar *grad   = &im1grad[z].Pixel(x, y, 0);
		  int sx = 0, sy = 0, sz = 0;

		  if (filter==esimple2) {

			  unsigned char *pix1 = &im1[z].Pixel(x, y, 0);
			  unsigned char *pix2 = 0;
			  unsigned char *pix3 = 0;
			  unsigned char *pix4 = 0;

			  double dx=0.0, dy=0.0, dz=0.0;

			  unsigned int pix1g, pix2g, pix3g, pix4g;
			  pix1g = (int)(0.11*pix1[0] + 0.59*pix1[1] + 0.3*pix1[2]);

			  if (x==width-1)
				  sx = 0;
			  else {
				  pix2 = &im1[z].Pixel(x+1, y, 0);
				  pix2g = (int)(0.11*pix2[0] + 0.59*pix2[1] + 0.3*pix2[2]);
				  dx = simple2[0]*pix2g + simple2[1]*pix1g;
			  }

			  if (y==height-1)
				  sy = 0;
			  else {
				  pix3 = &im1[z].Pixel(x, y+1, 0);
				  pix3g = (int)(0.11*pix3[0] + 0.59*pix3[1] + 0.3*pix3[2]);
				  dy = simple2[0]*pix3g + simple2[1]*pix1g;
			  }

			  if (z==depth-1)
				  sz = 0;
			  else {
				  pix4 = &im1[z+1].Pixel(x, y, 0);
				  pix4g = (int)(0.11*pix4[0] + 0.59*pix4[1] + 0.3*pix4[2]);
				  dz = simple2[0]*pix4g + simple2[1]*pix1g;
			  }

			  if (dx<0)
				  dx=-dx;
			  if (dy<0)
				  dy=-dy;
			  if (dz<0)
				  dz=-dz;

			  if (dx>upbound)
				  dx = upbound;
			  if (dy>upbound)
				  dy = upbound;
			  if (dz>upbound)
				  dz = upbound;

			  sx = (int) dx;
			  sy = (int) dy;
			  sz = (int) dz;

		  }

		  for (k = 0; gradThresh[k] <= sx; k++)
			; // find lowest bin for x gradient
		  assert(k < nK);
		  grad[0] = k;

		  for (k = 0; gradThresh[k] <= sy; k++)
			; // find lowest bin for y gradient
		  assert(k < nK);
		  grad[1] = k;

		  for (k = 0; gradThresh[k] <= sz; k++)
			; // find lowest bin for z gradient
		  assert(k < nK);
		  grad[2] = k;

		}
	  }
  }
}




// RGB converted to grayscale (x-y and z)
void computeQuantizedGradients(std::vector<CByteImage> im1, std::vector<CByteImage> &im1grad, std::vector<int> gradThresh, std::vector<int> gradThreshZ)
  /*
    use nK = gradThresh.size() levels of quantization.
    quantization is done as follows:
    compute absolute gradient (RMS difference of color bands) in x and y
    given gradient g, assign largest k such that gradThresh[k] <= g
    the last entry (gradThresh[nK-1]) is ignored and is assumed to be infinity
    in im1grad, use band 0 for x gradient, and band 1 for y gradients
  */
{
  int depth = im1.size();
  CShape sh = im1[0].Shape();
  int width = sh.width, height = sh.height, nB = sh.nBands;
  int nColors = __min(3, nB);
  int nColorsg = 1;

  enum FTLR filter = esimple2;

  int nK = (int)gradThresh.size();
  nK = nK + 1; // adding a 1 because will push maxlimit (inf) in gradThresh
  int k;

  int nKz = (int)gradThreshZ.size();
  nKz = nKz + 1; // adding a 1 because will push maxlimit (inf) in gradThresh
  int kz;


  DEBUG_OUT1(verbose, debugfile, "%d Gradient bins:\t0-", nK);
  for (k = 0; k < nK-1; k++)
	DEBUG_OUT2(verbose, debugfile, "%d\t%d-", gradThresh[k]-1, gradThresh[k]);
  DEBUG_OUT0(verbose, debugfile, "inf\n");

  DEBUG_OUT1(verbose, debugfile, "%d GradientZ bins:\t0-", nKz);
  for (kz = 0; kz < nKz-1; kz++)
	DEBUG_OUT2(verbose, debugfile, "%d\t%d-", gradThreshZ[kz]-1, gradThreshZ[kz]);
  DEBUG_OUT0(verbose, debugfile, "inf\n");



  sh.nBands = 3;  // band 0: x gradient, band 1: y gradient, band 2: z gradient

  for (int z = 0; z<depth; z++) {
	  im1grad[z].ReAllocate(sh);
  }

  // the constant factor is the sum of the +ve numbers in the filter matrix
  double upbound = 0, fctor=0;

  if (filter == esimple2) {
	  gradThresh.push_back(256 * nColorsg);
	  fctor = 1;
  }
  // need to fill in for the other filters

  upbound = 256 * fctor;

// need to decide whether we use the 2D filter or 3D filter
// for now i only have 2D filters
  // and so no restriction on the depth for loop
  for (int z = 0; z < depth; z++) {
	  for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {

		  uchar *grad   = &im1grad[z].Pixel(x, y, 0);
		  int sx = 0, sy = 0, sz = 0;

		  if (filter==esimple2) {

			  unsigned char *pix1 = &im1[z].Pixel(x, y, 0);
			  unsigned char *pix2 = 0;
			  unsigned char *pix3 = 0;
			  unsigned char *pix4 = 0;

			  double dx=0.0, dy=0.0, dz=0.0;

			  unsigned int pix1g, pix2g, pix3g, pix4g;
			  pix1g = (int)(0.11*pix1[0] + 0.59*pix1[1] + 0.3*pix1[2]);

			  if (x==width-1)
				  sx = 0;
			  else {
				  pix2 = &im1[z].Pixel(x+1, y, 0);
				  pix2g = (int)(0.11*pix2[0] + 0.59*pix2[1] + 0.3*pix2[2]);
				  dx = simple2[0]*pix2g + simple2[1]*pix1g;
			  }

			  if (y==height-1)
				  sy = 0;
			  else {
				  pix3 = &im1[z].Pixel(x, y+1, 0);
				  pix3g = (int)(0.11*pix3[0] + 0.59*pix3[1] + 0.3*pix3[2]);
				  dy = simple2[0]*pix3g + simple2[1]*pix1g;
			  }

			  if (z==depth-1)
				  sz = 0;
			  else {
				  pix4 = &im1[z+1].Pixel(x, y, 0);
				  pix4g = (int)(0.11*pix4[0] + 0.59*pix4[1] + 0.3*pix4[2]);
				  dz = simple2[0]*pix4g + simple2[1]*pix1g;
			  }

			  if (dx<0)
				  dx=-dx;
			  if (dy<0)
				  dy=-dy;
			  if (dz<0)
				  dz=-dz;

			  if (dx>upbound)
				  dx = upbound;
			  if (dy>upbound)
				  dy = upbound;
			  if (dz>upbound)
				  dz = upbound;

			  sx = (int) dx;
			  sy = (int) dy;
			  sz = (int) dz;

		  }

		  for (k = 0; gradThresh[k] <= sx; k++)
			; // find lowest bin for x gradient
		  assert(k < nK);
		  grad[0] = k;

		  for (k = 0; gradThresh[k] <= sy; k++)
			; // find lowest bin for y gradient
		  assert(k < nK);
		  grad[1] = k;

		  for (kz = 0; gradThreshZ[k] <= sz; kz++)
			; // find lowest bin for z gradient
		  assert(kz < nKz);
		  grad[2] = kz;

		}
	  }
  }
}



// compute quantized absolute color gradients of input image
void computeQuantizedGradients(CByteImage im1, CByteImage &im1grad, std::vector<int> gradThresh)
  /*
    use nK = gradThresh.size() levels of quantization.
    quantization is done as follows:
    compute absolute gradient (RMS difference of color bands) in x and y
    given gradient g, assign largest k such that gradThresh[k] <= g
    the last entry (gradThresh[nK-1]) is ignored and is assumed to be infinity
    in im1grad, use band 0 for x gradient, and band 1 for y gradients
  */
{
  CShape sh = im1.Shape();
  int width = sh.width, height = sh.height, nB = sh.nBands;
  int nColors = __min(3, nB);


  if ((int)gradThresh.size() == 0)
	  return; // this is the case for logistic regression with no edge boundaries.


  gradThresh.push_back(256 * 256 * nColors); // make sure last threshold is big enough

  int nK = (int)gradThresh.size();
  int k;

  DEBUG_OUT1(verbose, debugfile, "%d Gradient bins:\t0-", nK);
  for (k = 0; k < nK-1; k++)
	DEBUG_OUT2(verbose, debugfile, "%d\t%d-", gradThresh[k]-1, gradThresh[k]);
  DEBUG_OUT0(verbose, debugfile, "inf\n");

  for (k = 0; k < nK-1; k++) // compute sum of squared color differences below,
	gradThresh[k] *= gradThresh[k] * nColors;  //  so need to adjust thresholds to get RMS distance


  sh.nBands = 2;  // band 0: x gradient, band 1: y gradient
  im1grad.ReAllocate(sh);
  for (int y = 0; y < height; y++) {
	//printf("y=%d\n",y);
	for (int x = 0; x < width; x++) {
      uchar *pix   = &im1.Pixel(x, y, 0);
      uchar *pix1x = &im1.Pixel(x + (x < width-1), y, 0);
      uchar *pix1y = &im1.Pixel(x, y + (y < height-1), 0);
      uchar *grad   = &im1grad.Pixel(x, y, 0);
      int sx = 0, sy = 0;
      for (int b = 0; b < nColors; b++) {
		int dx = pix[b] - pix1x[b];
		int dy = pix[b] - pix1y[b];
		sx += dx * dx;
		sy += dy * dy;
      }
      for (k = 0; gradThresh[k] <= sx; k++)
		; // find lowest bin for x gradient
      assert(k < nK);
      grad[0] = k;
      for (k = 0; gradThresh[k] <= sy; k++)
		; // find lowest bin for y gradient
      assert(k < nK);
      grad[1] = k;
	}
  }
  if (0) {
	CByteImage g;
	BandSelect(im1grad, g, 0, 0);
	ScaleAndOffset(g, g, (float)256.0/nK, 0);
	WriteImageVerb(g, "grad0.png", true);
	BandSelect(im1grad, g, 1, 0);
	ScaleAndOffset(g, g, (float)256.0/nK, 0);
	WriteImageVerb(g, "grad1.png", true);
	//exit(1);
  }
}


// function-based smoothness cost computation
MRF::CostVal fnCost(int pix1, int pix2, int di, int dj)
{

  if (di == dj)
		return 0;

  MRF::CostVal gradcost = 0;
  if (pix1 < pix2) {
	if (pix2 - pix1 == 1)
      gradcost = hCueArrayV[cueIndex][pix1];
	else if (pix2 - pix1 == globalwidth[cueIndex]) {
      gradcost = vCueArrayV[cueIndex][pix1];
	} else {
		gradcost = dCueArrayV[cueIndex][pix1];
	}

  } else { // (pix2 < pix1)
	if (pix1 - pix2 == 1)
      gradcost = hCueArrayV[cueIndex][pix2];
	else if (pix1 - pix2 == globalwidth[cueIndex]) {
      gradcost = vCueArrayV[cueIndex][pix2];
	} else {
		gradcost = dCueArrayV[cueIndex][pix2];
	}
  }
  return gradcost;
}

fvec thetaVG, thetaZG, thetaCG;

int nDG;
int nGG;
int nGzG;


void copyThetaVZC(fvec thetaV, fvec thetaZ, fvec thetaC, int nD)
{

	thetaVG.erase(thetaVG.begin(), thetaVG.end());
	thetaZG.erase(thetaZG.begin(), thetaZG.end());
	thetaCG.erase(thetaCG.begin(), thetaCG.end());

	int nK = (int)thetaV.size();
	int nZ = (int)thetaZ.size();
	int nC = (int)thetaC.size();

	for (int k = 0; k < nK; k++) {
		thetaVG.push_back(thetaV[k]);
	}

	for (int k = 0; k < nZ; k++) {
		thetaZG.push_back(thetaZ[k]);
	}

	for (int k = 0; k < nC; k++) {
		thetaCG.push_back(thetaC[k]);
	}

	nDG = nD;
}





void setGlobalValues(int nD, int nG, int nGz)
{
  nDG = nD;
  nGG = nG;
  nGzG = nGz;
}


// // function-based smoothness cost computation
MRF::CostVal fnCostContext(int pix1, int pix2, int di, int dj)
{

  int mind, maxd;

  if (di == dj)
		return 0;

  MRF::CostVal smoothcost = 0;


  // gradient cost
  if (pix1 < pix2) {
	if (pix2 - pix1 == 1)
      smoothcost = thetaVG[hCueArrayVG[cueIndex][pix1]];
	else if (pix2 - pix1 == globalwidth[cueIndex]) {
      smoothcost = thetaVG[vCueArrayVG[cueIndex][pix1]];
	} else {
	  smoothcost = thetaVG[dCueArrayVG[cueIndex][pix1]];
	}

  } else { // (pix2 < pix1)
	if (pix1 - pix2 == 1)
      smoothcost = thetaVG[hCueArrayVG[cueIndex][pix2]];
	else if (pix1 - pix2 == globalwidth[cueIndex]) {
      smoothcost = thetaVG[vCueArrayVG[cueIndex][pix2]];
	} else {
	  smoothcost = thetaVG[dCueArrayVG[cueIndex][pix2]];
	}
  }

  //context cost
  if (di < dj) {
	  mind = di;
	  maxd = dj;
  } else {
	  mind = dj;
	  maxd = di;
  }

  int cnt = 0;
  for (int k=mind-1; k>=0; k--) {
	  cnt += (nDG-(k+1));
  }
  cnt += (maxd-mind-1);


  smoothcost += thetaCG[cnt];

  return smoothcost;
}

// function-based smoothness cost computation
MRF::CostVal fnCostGradContext(int pix1, int pix2, int di, int dj)
{

  // int mind, maxd;

//  if (di == dj)
//		return 0;

//  int nb = 0;
//  float tri = nDG*(nDG-1)/2;
//  nb = (int) thetaVG.size() / tri;
  // nb captures the number of bins per pair of distinct organs (background inclusive)

  //context
  // if (di < dj) {
//	  mind = di;
//	  maxd = dj;
  // } else {
//	  mind = dj;
//	  maxd = di;
//  }


  int gradval = 0;


  // gradient cost
  if (pix1 < pix2) {
	if (pix2 - pix1 == 1)
      gradval = hCueArrayVG[cueIndex][pix1];
	else if (pix2 - pix1 == globalwidth[cueIndex]) {
      gradval = vCueArrayVG[cueIndex][pix1];
	} else {
	  gradval = dCueArrayVG[cueIndex][pix1];
	}

  } else { // (pix2 < pix1)
	if (pix1 - pix2 == 1)
      gradval = hCueArrayVG[cueIndex][pix2];
	else if (pix1 - pix2 == globalwidth[cueIndex]) {
      gradval = vCueArrayVG[cueIndex][pix2];
	} else {
	  gradval = dCueArrayVG[cueIndex][pix2];
	}
  }

  int d1V=di, d2V=dj;

  if (d1V>d2V) {
    d1V = dj;
    d2V = di;
  }

  return thetaVG[(d1V*nDG + d2V - (d1V*(d1V+1))/2)*nGG + gradval];

/*
              distV[(d1*nD + d2 - d1*(d1+1)/2)*nG + g] ++;

 */
//  int cnt = 0;
//  for (int k=mind-1; k>=0; k--) {
//	  cnt += (nDG-(k+1))*nb;
//  }
//  cnt += (maxd-mind-1)*nb + gradval;

//  return thetaVG[cnt];
}


// for (x-y)
// function-based smoothness cost computation
MRF::CostVal fnCostGradContext00(int pix1, int pix2, int di, int dj)
{

  int gradval = 0;


  // gradient cost
  if (pix1 < pix2) {
	if (pix2 - pix1 == 1)
      gradval = hCueArrayVG[cueIndex][pix1];
	else if (pix2 - pix1 == globalwidth[cueIndex]) {
      gradval = vCueArrayVG[cueIndex][pix1];
	} else {

	}

  } else { // (pix2 < pix1)
	if (pix1 - pix2 == 1)
      gradval = hCueArrayVG[cueIndex][pix2];
	else if (pix1 - pix2 == globalwidth[cueIndex]) {
      gradval = vCueArrayVG[cueIndex][pix2];
	} else {

	}
  }

  int d1V=di, d2V=dj;

  if (d1V>d2V) {
    d1V = dj;
    d2V = di;
  }

  return thetaVG[(d1V*nDG + d2V - (d1V*(d1V+1))/2)*nGG + gradval];

}



// for x-y and z mode
MRF::CostVal fnCostGradContextZ(int pix1, int pix2, int di, int dj)
{


  int gradval = 0;
  int zflag = 0;

  // gradient cost
  if (pix1 < pix2) {
	if (pix2 - pix1 == 1)
      gradval = hCueArrayVG[cueIndex][pix1];
	else if (pix2 - pix1 == globalwidth[cueIndex]) {
      gradval = vCueArrayVG[cueIndex][pix1];
	} else {
		zflag = 1;
		gradval = dCueArrayVG[cueIndex][pix1];
	}

  } else { // (pix2 < pix1)
	if (pix1 - pix2 == 1)
      gradval = hCueArrayVG[cueIndex][pix2];
	else if (pix1 - pix2 == globalwidth[cueIndex]) {
      gradval = vCueArrayVG[cueIndex][pix2];
	} else {
		zflag = 1;
		gradval = dCueArrayVG[cueIndex][pix2];
	}
  }

  int d1V=di, d2V=dj;

  if (d1V>d2V) {
    d1V = dj;
    d2V = di;
  }

  if (zflag == 1)
	  return thetaZG[(d1V*nDG + d2V - (d1V*(d1V+1))/2)*nGzG + gradval];
  else
	  return thetaVG[(d1V*nDG + d2V - (d1V*(d1V+1))/2)*nGG + gradval];

}







SmoothnessCost *computeSmoothnessCost(CByteImage im1grad, std::vector<float> theta, int i)
{
  std::vector<float> thetaP;

  return computeSmoothnessCost(im1grad,theta,thetaP, i);
}


SmoothnessCost *computeSmoothnessCost(std::vector<CByteImage> im1grad, std::vector<float> theta, int i)
{
  return new SmoothnessCost( computeSmoothnessCostFunction (im1grad, theta, i) );
}

SmoothnessCost *computeSmoothnessCost(std::vector<CByteImage> im1grad, std::vector<float> thetaV, std::vector<float> thetaC, int i, int gradContext)
{
  return new SmoothnessCost( computeSmoothnessCostFunction (im1grad, thetaV, thetaC, i, gradContext) );
}

SmoothnessCost *computeSmoothnessCost00(std::vector<CByteImage> im1grad, std::vector<float> thetaV, std::vector<float> thetaC, int i, int gradContext)
{
  return new SmoothnessCost( computeSmoothnessCostFunction00 (im1grad, thetaV, thetaC, i, gradContext) );
}


SmoothnessCost *computeSmoothnessCost(std::vector<CByteImage> im1grad, std::vector<float> thetaV, std::vector<float> thetaZ, std::vector<float> thetaC, int i, int gradContext)
{
  return new SmoothnessCost( computeSmoothnessCostFunction (im1grad, thetaV, thetaZ, thetaC, i, gradContext) );
}


MRF::SmoothCostGeneralFn computeSmoothnessCostFunction(std::vector<CByteImage> im1grad,
													   std::vector<float> theta,
                                                       int i)
{
  checkForFloatCosts();
  int nK = (int)theta.size();

  int depth = im1grad.size();
  CShape sh = im1grad[0].Shape();
  int width = sh.width, height = sh.height;

  if (nK == 0) // run without smoothness params
	theta.push_back(1.0f);

  DEBUG_OUT1(verbose, debugfile, "%d Gradient weights:\t", nK);
  for (int k = 0; k < nK; k++)
	DEBUG_OUT1(verbose, debugfile, "%g\t", theta[k]);
  DEBUG_OUT0(verbose, debugfile, "\n");


  //DEBUG_OUT(verbose, debugfile, "Smoothness term: L%d norm, truncated at %d\n", smoothexp, smoothmax);

  // allocate the arrays when called the first time
  if (hCueArrayV.size() <= i)
	  hCueArrayV.push_back(new MRF::CostVal[depth * width * height]);
  if (vCueArrayV.size() <= i)
	  vCueArrayV.push_back(new MRF::CostVal[depth * width * height]);
  if (dCueArrayV.size() <= i)
	  dCueArrayV.push_back(new MRF::CostVal[depth * width * height]);

  int n = 0;
  for (int z = 0; z < depth; z++) {
	  for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
		  uchar* grad = &im1grad[z].Pixel(x, y, 0);
		  hCueArrayV[i][n] = theta[grad[0]];
		  vCueArrayV[i][n] = theta[grad[1]];
		  dCueArrayV[i][n] = theta[grad[2]];
		  n++;
		}
	  }
  }


  return fnCost;
}




MRF::SmoothCostGeneralFn computeSmoothnessCostFunction(std::vector<CByteImage> im1grad,
													   std::vector<float> thetaV,
													   std::vector<float> thetaC,
                                                       int i, int gradContext)
{
  checkForFloatCosts();

  int depth = im1grad.size();
  CShape sh = im1grad[0].Shape();
  int width = sh.width, height = sh.height;

  // allocate the arrays when called the first time
  if (hCueArrayVG.size() <= i)
	  hCueArrayVG.push_back(new unsigned int[depth * width * height]);
  if (vCueArrayVG.size() <= i)
	  vCueArrayVG.push_back(new unsigned int[depth * width * height]);
  if (dCueArrayVG.size() <= i)
	  dCueArrayVG.push_back(new unsigned int[depth * width * height]);

  int n = 0;
  for (int z = 0; z < depth; z++) {
	  for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
		  uchar* grad = &im1grad[z].Pixel(x, y, 0);
//		  DEBUG_OUT3(verbose, debugfile, " %d  %d  %d \t\t", grad[0], grad[1], grad[2]);
		  hCueArrayVG[i][n] = grad[0];
		  vCueArrayVG[i][n] = grad[1];
		  dCueArrayVG[i][n] = grad[2];

		  n++;
		}
	  }
  }

  if (gradContext==0)
	  return fnCostContext;
  else if(gradContext==1)
	  return fnCostGradContext;

}



MRF::SmoothCostGeneralFn computeSmoothnessCostFunction00(std::vector<CByteImage> im1grad,
													   std::vector<float> thetaV,
													   std::vector<float> thetaC,
                                                       int i, int gradContext)
{
  checkForFloatCosts();

  int depth = im1grad.size();
  CShape sh = im1grad[0].Shape();
  int width = sh.width, height = sh.height;

  // allocate the arrays when called the first time
  if (hCueArrayVG.size() <= i)
	  hCueArrayVG.push_back(new unsigned int[depth * width * height]);
  if (vCueArrayVG.size() <= i)
	  vCueArrayVG.push_back(new unsigned int[depth * width * height]);

  int n = 0;
  for (int z = 0; z < depth; z++) {
	  for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
		  uchar* grad = &im1grad[z].Pixel(x, y, 0);
//		  DEBUG_OUT3(verbose, debugfile, " %d  %d  %d \t\t", grad[0], grad[1], grad[2]);
		  hCueArrayVG[i][n] = grad[0];
		  vCueArrayVG[i][n] = grad[1];

		  n++;
		}
	  }
  }

  if (gradContext==0) // not coded for this case
	  return fnCostContext;
  else if(gradContext==1)
	  return fnCostGradContext00;

}



MRF::SmoothCostGeneralFn computeSmoothnessCostFunction(std::vector<CByteImage> im1grad,
													   std::vector<float> thetaV, std::vector<float> thetaZ,
													   std::vector<float> thetaC,
                                                       int i, int gradContext)
{
  checkForFloatCosts();

  int depth = im1grad.size();
  CShape sh = im1grad[0].Shape();
  int width = sh.width, height = sh.height;

  // allocate the arrays when called the first time
  if (hCueArrayVG.size() <= i)
	  hCueArrayVG.push_back(new unsigned int[depth * width * height]);
  if (vCueArrayVG.size() <= i)
	  vCueArrayVG.push_back(new unsigned int[depth * width * height]);
  if (dCueArrayVG.size() <= i)
	  dCueArrayVG.push_back(new unsigned int[depth * width * height]);

  int n = 0;
  for (int z = 0; z < depth; z++) {
	  for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
		  uchar* grad = &im1grad[z].Pixel(x, y, 0);
//		  DEBUG_OUT3(verbose, debugfile, " %d  %d  %d \t\t", grad[0], grad[1], grad[2]);
		  hCueArrayVG[i][n] = grad[0];
		  vCueArrayVG[i][n] = grad[1];
		  dCueArrayVG[i][n] = grad[2];

		  n++;
		}
	  }
  }

  if (gradContext==0) // need to work on or remove this condition
	  return fnCostContext;
  else if(gradContext==1)
	  return fnCostGradContextZ;

}







SmoothnessCost *computeSmoothnessCost(CByteImage im1grad, std::vector<float> theta, std::vector<float> thetaP, int i)
{

  return new SmoothnessCost( computeSmoothnessCostFunction (im1grad, theta, thetaP, i) );
}

MRF::SmoothCostGeneralFn computeSmoothnessCostFunction(CByteImage im1grad,
													   std::vector<float> theta, int i)
{
  std::vector<float> thetaP;
  return computeSmoothnessCostFunction(im1grad,theta,thetaP, i);
}

MRF::SmoothCostGeneralFn computeSmoothnessCostFunction(CByteImage im1grad,
													   std::vector<float> theta,
                                                       std::vector<float> thetaP, int i)
{
  checkForFloatCosts();
  int nK = (int)theta.size();
  int nP = (int)thetaP.size();

  CShape sh = im1grad.Shape();
  int width = sh.width, height = sh.height;

  if (nK == 0) // run without smoothness params
	theta.push_back(1.0f);

  DEBUG_OUT1(verbose, debugfile, "%d Gradient weights:\t", nK);
  for (int k = 0; k < nK; k++)
	DEBUG_OUT1(verbose, debugfile, "%g\t", theta[k]);
  DEBUG_OUT0(verbose, debugfile, "\n");


  //DEBUG_OUT(verbose, debugfile, "Smoothness term: L%d norm, truncated at %d\n", smoothexp, smoothmax);

  // allocate the arrays when called the first time
  if (hCueArrayV.size() <= i)
	  hCueArrayV.push_back(new MRF::CostVal[width * height]);
  if (vCueArrayV.size() <= i)
	  vCueArrayV.push_back(new MRF::CostVal[width * height]);

  //globallambda1 = theta[0];
  //globallambda2 = theta[1];
  //DEBUG_OUT(verbose, debugfile, "Smoothness term: lambda1 = %g, lambda2 = %g\n", globallambda1, globallambda2);

  int n = 0;
  for (int y = 0; y < height; y++) {
	for (int x = 0; x < width; x++) {
      uchar *grad = &im1grad.Pixel(x, y, 0);
      hCueArrayV[i][n] = theta[grad[0]];
      vCueArrayV[i][n] = theta[grad[1]];
      n++;
	}
  }


  if( nP == 4*nK){
	n = 0;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			uchar *grad = &im1grad.Pixel(x, y, 0);

			hBOcclArrayV[i][n] = thetaP[grad[0]];
			vBOcclArrayV[i][n] = thetaP[nK + grad[1]];

			hSOcclArrayV[i][n] = thetaP[2*nK + grad[0]];
			vSOcclArrayV[i][n] = thetaP[3*nK + grad[1]];
			n++;
		}
	}
  }

  // FUNCTION WRAPPING BREAKS THE ABILITY FOR THIS RETURN TYPE - JJW
//   if (0) { // define smoothnest cost using truncated linear
// 	MRF::CostVal lambda = 1;
// 	int smoothexp = 1;
// 	MRF::CostVal smoothmax = 1;
// 	// Potts model
// 	return new SmoothnessCost(smoothexp, smoothmax, lambda, hCueArray, vCueArray);
//   }
//   else { // define smoothness cost using general function
	return fnCost;
    //  }
}


// delete any arrays allocated for data and smoothness costs
void deleteGlobalArrays() {
  delete [] dsiArray;


  for(std::vector<MRF::CostVal *>::reverse_iterator b = dsiArray2.rbegin(); b != dsiArray2.rend(); ++b){
	  free(*b);
  }

  for(std::vector<MRF::CostVal *>::reverse_iterator b = hCueArrayV.rbegin(); b != hCueArrayV.rend(); ++b){
	  free(*b);
  }

  for(std::vector<MRF::CostVal *>::reverse_iterator b = vCueArrayV.rbegin(); b != vCueArrayV.rend(); ++b){
	  free(*b);
  }

  for(std::vector<MRF::CostVal *>::reverse_iterator b = hBOcclArrayV.rbegin(); b != hBOcclArrayV.rend(); ++b){
	  free(*b);
  }

  for(std::vector<MRF::CostVal *>::reverse_iterator b = vBOcclArrayV.rbegin(); b != vBOcclArrayV.rend(); ++b){
	  free(*b);
  }

  for(std::vector<MRF::CostVal *>::reverse_iterator b = hSOcclArrayV.rbegin(); b != hSOcclArrayV.rend(); ++b){
	  free(*b);
  }

  for(std::vector<MRF::CostVal *>::reverse_iterator b = vSOcclArrayV.rbegin(); b != vSOcclArrayV.rend(); ++b){
	  free(*b);
  }

  for(std::vector<MRF::CostVal *>::reverse_iterator b = dsiArrayV.rbegin(); b != dsiArrayV.rend(); ++b){
	  // free(*b);
	  delete [] *b;
  }

  for(std::vector<unsigned int *>::reverse_iterator b = hCueArrayVG.rbegin(); b != hCueArrayVG.rend(); ++b){
	  delete [] *b;
  }

  for(std::vector<unsigned int *>::reverse_iterator b = vCueArrayVG.rbegin(); b != vCueArrayVG.rend(); ++b){
	  delete [] *b;
  }

  for(std::vector<unsigned int *>::reverse_iterator b = dCueArrayVG.rbegin(); b != dCueArrayVG.rend(); ++b){
	  delete [] *b;
  }

  /*delete [][] hCueArray;
  delete [][] vCueArray;
  delete [][] hBOcclArray;
  delete [][] vBOcclArray;
  delete [][] hSOcclArray;
  delete [][] vSOcclArray;*/
}

// if WTAdisp is a valid image, it is used to initialize the labels

void crfmodel(int width, int height, int nD, DataCost *dcost, SmoothnessCost *scost,
               CByteImage WTAdisp, CByteImage &disp, float closeEnoughPercent)
{
  crfmodel (width, height, nD, dcost, scost, WTAdisp, disp, closeEnoughPercent,
             USE_GC,1);
}



MRF* crfmodel(int depth, int width, int height, int nD, DataCost *dcost, SmoothnessCost *scost,
               std::vector<CByteImage> &disp, float closeEnoughPercent,
               int inferencer, int innerIter, std::vector<CByteImage> interImage, std::vector<CByteImage> truedisp, int interactive, std::ofstream *logStream, int testv)
{

  if (innerIter == -1)
    innerIter = 2*width*height*depth;

  int outerIter = 1; //50; // also defined as the default while calling crfmodel()

  EnergyFunction *energy = new EnergyFunction(dcost, scost);

  MRF* mrf;

  unsigned int size3D = width*height*depth;
  unsigned int size2D = width*height;

  FloatType damper = 1;
  FloatType msgTol = 0.001;
//  MeanField::DiffType msgDiffFun = MeanField::ENERGY;
  MeanField::DiffType msgDiffFun = MeanField::KLD;

  switch (inferencer)
  {
  case USE_TRWS:
    //mrf = new TRWS(width, height, nD, energy);
	DEBUG_OUT0(verbose, debugfile, "NOT Creating TRWS\n");
	exit(1);
	break;
  case USE_MF:
  {
	mrf = new MeanField(size3D, nD, energy);

	mrf->setParameters(PARAM_DAMPER, &damper);
	mrf->setParameters(PARAM_MSGDTOL, &msgTol);
	mrf->setParameters(PARAM_MSGDFUN, &msgDiffFun);


	mrf->setLog2(logStream);


	innerIter = 20;
	outerIter = 5;

	DEBUG_OUT2(verbose, debugfile, "damper = %0.1f  tol = %g\n",
			 damper,msgTol);

    DEBUG_OUT0(verbose, debugfile, "Creating MeanField\n");
    break;
  }
  case USE_SMF:
//    mrf = new SparseMeanField(width, height, nD, energy);
    DEBUG_OUT0(verbose, debugfile, "Not Creating SparseMeanField\n");
    exit(1);
    break;
  case USE_BP:
//    mrf = new SyncSumProd(width, height, nD, energy);
    DEBUG_OUT0(verbose, debugfile, "NOT Creating SyncSumProd\n");
    exit(1);
    break;
  case USE_SBP:
    {
      DEBUG_OUT0(verbose, debugfile, "NOT Creating SparseSyncSumProd\n");
      exit(1);
      break;
    }
  case USE_GC:
    {
	  mrf = new Expansion(size3D, nD, energy);

	  DEBUG_OUT0(verbose, debugfile, "Creating Graph Cuts / Expansion\n");
	  std::cout << " in Graph cuts \n";

	  break;
    }
  case USE_BPS: {

	  mrf = new BPS(size3D, nD, energy);

	  DEBUG_OUT0(verbose, debugfile, "Creating BP-scaline \n");

	  innerIter = 100;


	  break;
  }
  case USE_MP: {

	  mrf = new MaxProdBP(size3D, nD, energy);

	  DEBUG_OUT0(verbose, debugfile, "Creating Max Product \n");

	  innerIter = 10;

	  break;
  }
  default:
	  fprintf(stderr, "Unknown inferencer type");
	  exit(1);
  }

  // set the neighborhood of the general graph

  if (globalP.getOpflowConnection() == 1)
  {

    int volidx = globalP.getOpflowVolumeIndex();
    std::map <std::pair <int, int>, std::pair <int, MRF::CostVal> >::iterator iter = globalP.edgeGlobalF[volidx].begin();

    for (; iter != globalP.edgeGlobalF[volidx].end(); iter++)
    {
      int node1 = iter->first.first;
      int node2 = iter->first.second;

      if (interactive == 0 || (interactive == 1 && testv == 0))
      {
        mrf->setNeighbors(node1, node2, 1);
      } else
      {
        int z1 = node1 / (width*height);
        int z2 = node2 / (width*height);

        int rem1 = node1 - z1*(width*height);
        int rem2 = node2 - z2*(width*height);

        int y1 = rem1 / width;
        int y2 = rem2 / width;

        int x1 = rem1 - y1*width;
        int x2 = rem2 - y2*width;

        uchar *inimp = &interImage[z1].Pixel(x1, y1, 0);
        uchar *inimq = &interImage[z2].Pixel(x2, y2, 0);

        if (*inimp==0 && *inimq==0)
          mrf->setNeighbors(node1, node2, 1);

      }
    }
  } else { // opflowConnection = 0

    if (interactive == 0 || (interactive == 1 && testv == 0))
    {

      for (int m=0; m < depth-1; m++) {
        for (unsigned int k=0; k < size2D; k++) {
          mrf->setNeighbors(m*size2D+k, (m+1)*size2D+k, 1);
        }
      }

      for (unsigned int k=0; k < size3D; k++) {
        unsigned int kd = k%size2D;

        if (kd%width == width-1 && kd/width == height-1) {
          // last pix in 2D slice, do nothing
        } else if (kd/width == height-1) {
          // last row
          mrf->setNeighbors(k, k+1, 1);
        } else if (kd%width == width-1) {
          // last column
          mrf->setNeighbors(k, k+width, 1);
        } else {
          mrf->setNeighbors(k, k+1, 1);
          mrf->setNeighbors(k, k+width, 1);
        }
      }

    } else
    {

      for (int m=0; m < depth-1; m++) {
        for (unsigned int k=0; k < size2D; k++) {
          int xx = k%width;
          int yy = k/width;

          uchar *inimp = &interImage[m].Pixel(xx, yy, 0);
          uchar *inimq = &interImage[m+1].Pixel(xx, yy, 0);

          if (*inimp==0 && *inimq==0)
            mrf->setNeighbors(m*size2D+k, (m+1)*size2D+k, 1);
        }
      }

      for (unsigned int k=0; k < size3D; k++) {
        unsigned int kd = k%size2D;
        unsigned int m = k/size2D;
        int xx = kd%width;
        int yy = kd/width;

        if (kd%width == width-1 && kd/width == height-1) {
          // last pix in 2D slice, do nothing
        } else if (kd/width == height-1) {
          // last row

          uchar *inimp = &interImage[m].Pixel(xx, yy, 0);
          uchar *inimq = &interImage[m].Pixel(xx+1, yy, 0);

          if (*inimp==0 && *inimq==0)
            mrf->setNeighbors(k, k+1, 1);

        } else if (kd%width == width-1) {
          // last column

          uchar *inimp = &interImage[m].Pixel(xx, yy, 0);
          uchar *inimq = &interImage[m].Pixel(xx, yy+1, 0);

          if (*inimp==0 && *inimq==0)
            mrf->setNeighbors(k, k+width, 1);

        } else {

          uchar *inimp = &interImage[m].Pixel(xx, yy, 0);
          uchar *inimq = &interImage[m].Pixel(xx+1, yy, 0);

          if (*inimp==0 && *inimq==0)
            mrf->setNeighbors(k, k+1, 1);

          inimp = &interImage[m].Pixel(xx, yy, 0);
          inimq = &interImage[m].Pixel(xx, yy+1, 0);

          if (*inimp==0 && *inimq==0)
            mrf->setNeighbors(k, k+width, 1);
        }
      }
    }
  }



  if(inferencer != USE_TRWS && inferencer != USE_BPS && inferencer != USE_MF && inferencer != USE_MP)
    mrf->setParameters(1,&randomv);

  mrf->initialize();
  mrf->clearAnswer();

  // initialize labels to previous disps
  int n = 0;
  for (int z = 0; z < depth; z++)
  {
    for (int y = 0; y < height; y++)
    {
      uchar *row = &disp[z].Pixel(0, y, 0);
      for (int x = 0; x < width; x++)
      {
        mrf->setLabel(n++, (int)row[x]);
      }
    }
  }

  delete energy;

  // std::cout << " just before inner crfmodel call \n";
  crfmodel(depth, width, height, disp, mrf, closeEnoughPercent, innerIter, outerIter);
  // std::cout << " just after inner crfmodel call \n";

  if (interactive == 1 && testv == 1)
  {
	  // copy the gtruth results to the mrf for the pixels that are manually labeled in interacrtive mode.

    int n = 0;

    for (int z = 0; z < depth; z++)
    {
		  for (int y = 0; y < height; y++)
			{
			  uchar *row = &truedisp[z].Pixel(0, y, 0);
			  uchar *inim = &interImage[z].Pixel(0, y, 0);

			  for (int x = 0; x < width; x++)
				{
				  if (inim[x]!=0) // it is interactive, then update result with gtruth
				  {
					  mrf->setLabel(n, row[x]); // set mrf with row[x] which has gtruth

					  uchar *rowv = &disp[z].Pixel(x, y, 0); // also set the disp image
					  *rowv = row[x];
				  }
				  n++;
				}
			}
    }
  }

  return mrf;

}




void crfmodel(int depth, int width, int height, std::vector<CByteImage> &disp, MRF* mrf,
               float closeEnoughPercent, int innerIter, int maxIter )
{
  MRF::EnergyVal E = 0, Ed = 0, Es = 0, Eold = 0;
  Ed = mrf->dataEnergy();
  Es = mrf->smoothnessEnergy();
  E = Ed + Es; // mrf->totalEnergy();

  //std::cout << "*&* E = " << E << "    Ed = " << Ed << "  Es = " << Es << std::endl;

  //DEBUG_OUT3(verbose, debugfile, "Energy = %f (Ed=%g, Es=%g) at start\n", E, Ed, Es);
  fprintf(debugfile, "Energy = %f (Ed=%g, Es=%g) at start\n", E, Ed, Es);

  Eold = E;

  float t, tot_t = 0;

  for (int iter = 0; iter < maxIter; iter++)
  {
    // std::cout << " Iter : " << iter << "\n";
    mrf->optimize(innerIter, t);
    tot_t += t ;
    // std::cout << " After optimize \n";
    Ed = mrf->dataEnergy();
    Es = mrf->smoothnessEnergy();
    E = Ed + Es; // mrf->totalEnergy();

    // std::cout << " just after adding energies \n";

    if (fabs(Eold) < 1e-30)
      Eold = 1e-30;

    float percentDecrease;

    if (fabs(E) < 1e-30)
      percentDecrease = 0.0;
    else
      percentDecrease = (float)(100.0*((Eold - E)/fabs(Eold)));

    std::cout << "Energy = " << E << " (Ed= " << Ed << ", Es= " << Es <<" ), " << percentDecrease  << " decrease \n";
    DEBUG_OUT5(verbose, debugfile, "Energy = %f (Ed=%g, Es=%g), %.3f secs, %.3f%% decrease\n", E, Ed, Es, tot_t, percentDecrease);

    if (percentDecrease < closeEnoughPercent)
      break;

    if (E > Eold)
      fprintf(stderr, "Warning: energy is increasing!\n");

    Eold = E;
  }

  // get disparities

  CShape sh(width, height, 1);
  for (int z=0; z<depth; z++)
	  disp[z].ReAllocate(sh);

  int n = 0;
  for (int z=0; z<depth; z++) {
	  for (int y = 0; y < height; y++)
		{
		  uchar *row = &disp[z].Pixel(0, y, 0);
		  for (int x = 0; x < width; x++)
			{
			  row[x] = mrf->getLabel(n++);
			}
		}
  }
}





void crfmodel(int width, int height, int nD, DataCost *dcost, SmoothnessCost *scost,
               CByteImage WTAdisp, CByteImage &disp, float closeEnoughPercent,
               int inferencer, int innerIter)
{

  if (innerIter == -1)
    innerIter = 5*width*height;

  EnergyFunction *energy = new EnergyFunction(dcost, scost);

  MRF* mrf;

  switch (inferencer){
  case USE_TRWS:
    mrf = new TRWS(width, height, nD, energy);
	DEBUG_OUT0(verbose, debugfile, "Creating TRWS\n");
	break;
  case USE_MF:
    mrf = new MeanField(width, height, nD, energy);
    DEBUG_OUT0(verbose, debugfile, "Creating MeanField\n");
    break;
  case USE_SMF:
    mrf = new SparseMeanField(width, height, nD, energy);
    DEBUG_OUT0(verbose, debugfile, "Creating SparseMeanField\n");
    break;
  case USE_BP:
    mrf = new SyncSumProd(width, height, nD, energy);
    DEBUG_OUT0(verbose, debugfile, "Creating SyncSumProd\n");
    break;
  case USE_SBP:
    {
      mrf = new SparseSyncSumProd(width, height, nD, energy);
      FloatType damper = 0.9;
      mrf->setParameters(PARAM_DAMPER, &damper);

      FloatType msgTol = 0.0001 * 2 * (width*(height-1) + (width-1)*height);
      mrf->setParameters(PARAM_MSGDTOL, &msgTol);

      FloatType divTol = 0.001;
      mrf->setParameters(PARAM_DIVTOL, &divTol);

      DEBUG_OUT0(verbose, debugfile, "Creating SparseSyncSumProd\n");
      break;
    }
  case USE_GC:
    mrf = new Expansion(width, height, nD, energy);
    // chetos - use other expansion constructor and need to also call setNeighbours
    DEBUG_OUT0(verbose, debugfile, "Creating Graph Cuts / Expansion\n");
    std::cout << " Inside inner crfmodel : Creating Graph Cuts \n";
    break;
  default:
    fprintf(stderr, "Unknown inferencer type");
    exit(1);
  }

  if(inferencer != USE_TRWS)
	mrf->setParameters(1,&randomv);
  mrf->initialize();

  if (WTAdisp.Shape().width == width)
    {
      // initialize labels to WTA disps
      int n = 0;

      for (int y = 0; y < height; y++)
        {
          uchar *row = &WTAdisp.Pixel(0, y, 0);

          for (int x = 0; x < width; x++)
            {
              mrf->setLabel(n++, row[x]);
            }
        }
    }
  else
    {
      mrf->clearAnswer();
    }

  delete energy;

  crfmodel(width, height, disp, mrf, closeEnoughPercent, innerIter );

  delete mrf;
}

void crfmodel(int width, int height, CByteImage &disp, MRF* mrf,
               float closeEnoughPercent, int innerIter, int maxIter )
{
  MRF::EnergyVal E, Ed, Es, Eold;
  Ed = mrf->dataEnergy();
  Es = mrf->smoothnessEnergy();
  E = Ed + Es; // mrf->totalEnergy();

  DEBUG_OUT3(verbose, debugfile, "Energy = %f (Ed=%g, Es=%g) at start\n", E, Ed, Es);
  std::cout << " Inner inner crf : Energy computed before optimize called \n";


  Eold = E;

  float t, tot_t = 0;

  for (int iter = 0; iter < maxIter; iter++)
    {
      // std::cout << " Iter : " << iter << "\n";
      mrf->optimize(innerIter, t);
      // std::cout << " After optimize \n";
      tot_t += t ;

      Ed = mrf->dataEnergy();
      Es = mrf->smoothnessEnergy();
      E = Ed + Es; // mrf->totalEnergy();

      std::cout << " Energy = " << Ed <<"  " << Es << "  " << E << "\n";

      float percentDecrease = (float)(100.0*(1.0 - E/Eold));

      DEBUG_OUT5(verbose, debugfile, "Energy = %f (Ed=%g, Es=%g), %.3f secs, %.3f%% decrease\n", E, Ed, Es, tot_t, percentDecrease);

      if (percentDecrease < closeEnoughPercent)
	    break;

      if (E > Eold)
	    fprintf(stderr, "Warning: energy is increasing!\n");

      Eold = E;
    }

  // get disparities
  CShape sh(width, height, 1);
  disp.ReAllocate(sh);

  int n = 0;
  for (int y = 0; y < height; y++)
    {
      uchar *row = &disp.Pixel(0, y, 0);
      for (int x = 0; x < width; x++)
        {
          row[x] = mrf->getLabel(n++);
        }
    }

  std::cout << " End of inner inner crfmodel \n";
}

void writeDisparities(CByteImage disp, int outscale, char *outstem)
{
  CByteImage disp2;
//  DEBUG_OUT1(verbose, debugfile, "scaling disparities by %d\n", outscale);
  ScaleAndOffset(disp, disp2, (float)outscale, 0);

  char dispname[1000];
  sprintf(dispname, "%s.png", outstem);
//c  DEBUG_OUT1(verbose, debugfile, "Writing image %s\n", dispname);
  WriteImage(disp2, dispname);
}

void writeDisparities(std::vector<CByteImage> disp, int outscale, char *outstem)
{
  CByteImage disp2;
//  DEBUG_OUT1(verbose, debugfile, "scaling disparities by %d\n", outscale);
  for (int z=0; z<disp.size(); z++) {
	  ScaleAndOffset(disp[z], disp2, (float)outscale, 0);

	  char dispname[1000];
	  sprintf(dispname, "%s_%d.png", outstem,z);
//c	  DEBUG_OUT1(verbose, debugfile, "Writing image %s\n", dispname);
	  WriteImage(disp2, dispname);
  }
}

void setPairwise(bool v){
	pairwise = v;
}

void setPairwiseInteraction(bool v){
	pairwiseInteraction = v;
}

void setLocal(bool v){
	localOccl = v;
}

bool getLocal(){
	return localOccl;
}


void updateThetaP(std::vector<float> theta){
	int nP = theta.size();
	if(nP == 2){
		occluded1 = theta[1]; // Only one of the pair is occluded.
		occluded2 = theta[0]; // Both Pair are occluded.
	}
}

int getNumStates(){
	return globalNdisps;
}

float rangeNormailize(float rawMin, float rawMax, float min, float max, float real){
	if(real > rawMax){
		return max;
	}
	else if(real < rawMin){
		return min;
	}
	else{
		return (max - min) * (real - rawMin )/ ( rawMax - rawMin) + min;
	}
}

void updateCueIndex(unsigned int i){
	cueIndex = i;
}

void setRandom(bool rv){
	randomv = rv;
}

// this function is no longer in use

double estimateLogLikelihood(std::vector<CByteImage> gtImage, std::vector<float> thetaU, std::vector<float> thetaV, int index, int nD ) {

	double logLikelihood = 0.0;

	int depth = gtImage.size();
	CShape sh = gtImage[0].Shape();
	int width = sh.width, height = sh.height, nB = sh.nBands;

	double dataCostVal = 0.0, smoothCostVal = 0.0;

	unsigned int n=0;

	for (int z=0; z < depth; z++) {
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				uchar *pix   = &gtImage[z].Pixel(x, y, 0);

				  // The cost of pixel p and label l is stored at dsiArray[p*nLabels+l]

				dataCostVal += (dsiArrayV[index])[n*nD + pix[0]];

				smoothCostVal += hCueArrayV[index][n] + vCueArrayV[index][n] + dCueArrayV[index][n];

				n++;

			}
		}
	}

	logLikelihood = -dataCostVal-smoothCostVal;

	return logLikelihood;
}







void logisticRegressionPairwise(std::vector <CImage2> im1, int nD, int numTrainingPats, std::vector <CByteImage> &wta, fvec thetaV, std::vector<int> gradThresh, std::vector<CByteImage> im1grad, std::vector<int> startSliceNo, std::vector<int> endSliceNo, int index, double &loglikelihood, std::vector <CByteImage> gtIm)
{

  checkForFloatCosts();

  int nV = thetaV.size();
  int nG = gradThresh.size()+1;


  // for now using only the widthwise pixel pairwise terms
  if (wpairwiseProbV.size() == 0)
    for (int ii=0; ii<nG; ii++)
      wpairwiseProbV.push_back(new MRF::CostVal[nD*nD]);


  MRF::CostVal badcost = 1e20;

  for(int j=0; j<nV; j++)
	  badcost += thetaV[j];

  int thetaid;

  for (int tid = 0; tid < nG; tid++) {
    std::vector<double> sumP(nD*nD);
    int ih = 0;
    //    double sumP = 0.0;
    for (int d1 = 0; d1 < nD; d1++) {
      for (int d2 = d1; d2 < nD; d2++) {
        
        // d1 will always be less than or equal to d2

        wpairwiseProbV[tid][d1*nD+d2] = - thetaV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + tid];  // the -ve sign because we are modeling the exponent as E_s or cost
        wpairwiseProbV[tid][d2*nD+d1] = wpairwiseProbV[tid][d1*nD+d2];

        if (d1==d2)
          sumP[ih] = wpairwiseProbV[tid][d2*nD+d1];
        else {
          sumP[ih] = wpairwiseProbV[tid][d2*nD+d1];
          ih++;
          sumP[ih] = wpairwiseProbV[tid][d2*nD+d1];
        }
        ih++;

      }
    }

    double logSumP  = getLogSumExp(sumP);


    for (int d1 = 0; d1 < nD; d1++) {
      for (int d2 = d1; d2 < nD; d2++) {
        wpairwiseProbV[tid][d1*nD+d2] = wpairwiseProbV[tid][d1*nD+d2] - logSumP;
        wpairwiseProbV[tid][d2*nD+d1] = wpairwiseProbV[tid][d1*nD+d2];
      }
    }



  }


  int depth = im1.size();
  CShape sh = im1[0].Shape();
  int width = sh.width, height = sh.height, nB = sh.nBands;
  int tempglobalNpixels = depth * width * height;
  int nColors = __min(3, nB);

  int zslice = 0;
  zslice = startSliceNo[index];

  int dsiIndex = 0;

  int* bestd1 = new int[nD*nD];
  int* bestd2 = new int[nD*nD];

  for(unsigned int i = 0; i < im1.size(); ++i){

    DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);
    sh.nBands = 1;
    wta[i].ReAllocate(sh);

    for (int y = 0; y < height; y++) {
      uchar *WTArow = &wta[i].Pixel(0, y, 0);
      for (int x = 0; x < width-1; x=x+2) {
			
        //unsigned short pix1 = im1[i].Pixel(x, y, 0);
        // unsigned short pix2 = im1[i].Pixel(x+1, y, 0);

        uchar *grad = &im1grad[i].Pixel(x, y, 0);

        int tidcurr = grad[0];
                  
        // std::cout<< tidcurr << "\n";

        int numbest = 0;
        MRF::CostVal bestval = badcost;
                 
        for (int d1 = 0; d1 < nD; d1++) {
          for (int d2 = d1; d2 < nD; d2++) {

            MRF::CostVal dsiValue = wpairwiseProbV[tidcurr][d1*nD+d2];
            // dsiValue is now the probability,
            // and we want the one with max prob.
            if (dsiValue > bestval || numbest == 0) {
              bestval = dsiValue;
              bestd1[0] = d1;
              bestd2[0] = d2;
              numbest = 1;
            } else if (dsiValue == bestval) {
              bestd1[numbest] = d1;
              bestd2[numbest] = d2;
              numbest++;
            } 
          }
        }
                  
        // selecting the best candidates from the list uniformly

        if (numbest == 0)
        	numbest = nD;

        int curr = rand() % numbest;
                  
        WTArow[x] = bestd1[curr];
        WTArow[x+1] = bestd2[curr];
                 
        uchar *gtpix1 = &gtIm[i].Pixel(x, y, 0);
        uchar *gtpix2 = &gtIm[i].Pixel(x+1, y, 0);

        loglikelihood += wpairwiseProbV[tidcurr][gtpix1[0]*nD+gtpix2[0]];

      }
              
    }
  }

  delete [] bestd1;
  delete [] bestd2;

}


void logisticRegressionPairwise(int nD, int numTrainingPats, std::vector <CByteImage> &wta, fvec thetaV, std::vector<int> gradThresh, std::vector<CByteImage> im1grad, std::vector<int> startSliceNo, std::vector<int> endSliceNo, int index, double &loglikelihood, std::vector <CByteImage> gtIm)
{

  checkForFloatCosts();

  int nV = thetaV.size();
  int nG = gradThresh.size()+1;


  // for now using only the widthwise pixel pairwise terms
  if (wpairwiseProbV.size() == 0)
    for (int ii=0; ii<nG; ii++)
      wpairwiseProbV.push_back(new MRF::CostVal[nD*nD]);


  MRF::CostVal badcost = 1e20;

  for(int j=0; j<nV; j++)
	  badcost += thetaV[j];

  int thetaid;

  for (int tid = 0; tid < nG; tid++) {
    std::vector<double> sumP(nD*nD);
    int ih = 0;
    //    double sumP = 0.0;
    for (int d1 = 0; d1 < nD; d1++) {
      for (int d2 = d1; d2 < nD; d2++) {

        // d1 will always be less than or equal to d2

        wpairwiseProbV[tid][d1*nD+d2] = - thetaV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + tid];  // the -ve sign because we are modeling the exponent as E_s or cost
        wpairwiseProbV[tid][d2*nD+d1] = wpairwiseProbV[tid][d1*nD+d2];

        if (d1==d2)
          sumP[ih] = wpairwiseProbV[tid][d2*nD+d1];
        else {
          sumP[ih] = wpairwiseProbV[tid][d2*nD+d1];
          ih++;
          sumP[ih] = wpairwiseProbV[tid][d2*nD+d1];
        }
        ih++;

      }
    }

    double logSumP  = getLogSumExp(sumP);


    for (int d1 = 0; d1 < nD; d1++) {
      for (int d2 = d1; d2 < nD; d2++) {
        wpairwiseProbV[tid][d1*nD+d2] = wpairwiseProbV[tid][d1*nD+d2] - logSumP;
        wpairwiseProbV[tid][d2*nD+d1] = wpairwiseProbV[tid][d1*nD+d2];
      }
    }



  }


  int depth = gtIm.size();
  CShape sh = gtIm[0].Shape();
  int width = sh.width, height = sh.height, nB = sh.nBands;
  int tempglobalNpixels = depth * width * height;
  int nColors = __min(3, nB);

  int zslice = 0;
  zslice = startSliceNo[index];

  int dsiIndex = 0;

  int* bestd1 = new int[nD*nD];
  int* bestd2 = new int[nD*nD];

  for(unsigned int i = 0; i < gtIm.size(); ++i){

    DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);
    sh.nBands = 1;
    wta[i].ReAllocate(sh);

    for (int y = 0; y < height; y++) {
      uchar *WTArow = &wta[i].Pixel(0, y, 0);
      for (int x = 0; x < width-1; x=x+2) {

        //unsigned short pix1 = im1[i].Pixel(x, y, 0);
        // unsigned short pix2 = im1[i].Pixel(x+1, y, 0);

        uchar *grad = &im1grad[i].Pixel(x, y, 0);

        int tidcurr = grad[0];

        // std::cout<< tidcurr << "\n";

        int numbest = 0;
        MRF::CostVal bestval = badcost;

        for (int d1 = 0; d1 < nD; d1++) {
          for (int d2 = d1; d2 < nD; d2++) {

            MRF::CostVal dsiValue = wpairwiseProbV[tidcurr][d1*nD+d2];
            // dsiValue is now the probability,
            // and we want the one with max prob.
            if (dsiValue > bestval || numbest == 0) {
              bestval = dsiValue;
              bestd1[0] = d1;
              bestd2[0] = d2;
              numbest = 1;
            } else if (dsiValue == bestval) {
              bestd1[numbest] = d1;
              bestd2[numbest] = d2;
              numbest++;
            }
          }
        }

        // selecting the best candidates from the list uniformly
        if (numbest == 0)
        	numbest = nD;
        int curr = rand() % numbest;

        WTArow[x] = bestd1[curr];
        WTArow[x+1] = bestd2[curr];

        uchar *gtpix1 = &gtIm[i].Pixel(x, y, 0);
        uchar *gtpix2 = &gtIm[i].Pixel(x+1, y, 0);

        loglikelihood += wpairwiseProbV[tidcurr][gtpix1[0]*nD+gtpix2[0]];

      }

    }
  }

  delete [] bestd1;
  delete [] bestd2;

}



// for x-y and z
void logisticRegressionPairwise(int nD, int numTrainingPats, std::vector <CByteImage> &wta, fvec thetaV, std::vector<int> gradThresh, fvec thetaZ, std::vector<int> gradThreshZ, std::vector<CByteImage> im1grad, std::vector<int> startSliceNo, std::vector<int> endSliceNo, int index, double &loglikelihood, std::vector <CByteImage> gtIm)
{

  checkForFloatCosts();

  int nV = thetaV.size();
  int nG = gradThresh.size()+1;

  int nZ = thetaZ.size();
  int nGz = gradThreshZ.size()+1;


  // for now using only the widthwise pixel pairwise terms
  if (wpairwiseProbV.size() == 0)
    for (int ii=0; ii<nG; ii++)
      wpairwiseProbV.push_back(new MRF::CostVal[nD*nD]);

  if (wpairwiseProbZ.size() == 0)
    for (int ii=0; ii<nGz; ii++)
      wpairwiseProbZ.push_back(new MRF::CostVal[nD*nD]);


  MRF::CostVal badcost = 1e20;

  for(int j=0; j<nV; j++)
	  badcost += thetaV[j];

//  int thetaid;

  for (int tid = 0; tid < nG; tid++) {
    std::vector<double> sumP(nD*nD);
    int ih = 0;
    //    double sumP = 0.0;
    for (int d1 = 0; d1 < nD; d1++) {
      for (int d2 = d1; d2 < nD; d2++) {

        // d1 will always be less than or equal to d2

        wpairwiseProbV[tid][d1*nD+d2] = - thetaV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + tid];  // the -ve sign because we are modeling the exponent as E_s or cost
        wpairwiseProbV[tid][d2*nD+d1] = wpairwiseProbV[tid][d1*nD+d2];

        if (d1==d2)
          sumP[ih] = wpairwiseProbV[tid][d2*nD+d1];
        else {
          sumP[ih] = wpairwiseProbV[tid][d2*nD+d1];
          ih++;
          sumP[ih] = wpairwiseProbV[tid][d2*nD+d1];
        }
        ih++;

      }
    }

    double logSumP  = getLogSumExp(sumP);


    for (int d1 = 0; d1 < nD; d1++) {
      for (int d2 = d1; d2 < nD; d2++) {
        wpairwiseProbV[tid][d1*nD+d2] = wpairwiseProbV[tid][d1*nD+d2] - logSumP;
        wpairwiseProbV[tid][d2*nD+d1] = wpairwiseProbV[tid][d1*nD+d2];
      }
    }
  }


  badcost = 1e20;

  for(int j=0; j<nZ; j++)
	  badcost += thetaZ[j];

  for (int tid = 0; tid < nGz; tid++) {
    std::vector<double> sumP(nD*nD);
    int ih = 0;
    //    double sumP = 0.0;
    for (int d1 = 0; d1 < nD; d1++) {
      for (int d2 = d1; d2 < nD; d2++) {

        // d1 will always be less than or equal to d2

        wpairwiseProbZ[tid][d1*nD+d2] = - thetaZ[(d1*nD + d2 - (d1*(d1+1))/2)*nGz + tid];  // the -ve sign because we are modeling the exponent as E_s or cost
        wpairwiseProbZ[tid][d2*nD+d1] = wpairwiseProbZ[tid][d1*nD+d2];

        if (d1==d2)
          sumP[ih] = wpairwiseProbZ[tid][d2*nD+d1];
        else {
          sumP[ih] = wpairwiseProbZ[tid][d2*nD+d1];
          ih++;
          sumP[ih] = wpairwiseProbZ[tid][d2*nD+d1];
        }
        ih++;

      }
    }

    double logSumP  = getLogSumExp(sumP);


    for (int d1 = 0; d1 < nD; d1++) {
      for (int d2 = d1; d2 < nD; d2++) {
        wpairwiseProbZ[tid][d1*nD+d2] = wpairwiseProbZ[tid][d1*nD+d2] - logSumP;
        wpairwiseProbZ[tid][d2*nD+d1] = wpairwiseProbZ[tid][d1*nD+d2];
      }
    }
  }



  int depth = gtIm.size();
  CShape sh = gtIm[0].Shape();
  int width = sh.width, height = sh.height, nB = sh.nBands;
  int tempglobalNpixels = depth * width * height;
  int nColors = __min(3, nB);

//  int zslice = 0;
//  zslice = startSliceNo[index];

  int dsiIndex = 0;

  int* bestd1 = new int[nD*nD];
  int* bestd2 = new int[nD*nD];

  for(unsigned int i = 0; i < gtIm.size(); ++i){

    DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);
    sh.nBands = 1;
    wta[i].ReAllocate(sh);

    for (int y = 0; y < height; y++) {
      uchar *WTArow = &wta[i].Pixel(0, y, 0);
      for (int x = 0; x < width-1; x=x+2) {

        //unsigned short pix1 = im1[i].Pixel(x, y, 0);
        // unsigned short pix2 = im1[i].Pixel(x+1, y, 0);

        uchar *grad = &im1grad[i].Pixel(x, y, 0);

        int tidcurr = grad[0];

        // std::cout<< tidcurr << "\n";

        int numbest = 0;
        MRF::CostVal bestval = badcost;

        for (int d1 = 0; d1 < nD; d1++) {
          for (int d2 = d1; d2 < nD; d2++) {

            MRF::CostVal dsiValue = wpairwiseProbV[tidcurr][d1*nD+d2];
            // dsiValue is now the probability,
            // and we want the one with max prob.
            if (dsiValue > bestval || numbest == 0) {
              bestval = dsiValue;
              bestd1[0] = d1;
              bestd2[0] = d2;
              numbest = 1;
            } else if (dsiValue == bestval) {
              bestd1[numbest] = d1;
              bestd2[numbest] = d2;
              numbest++;
            }
          }
        }

        // selecting the best candidates from the list uniformly
        if (numbest == 0)
        	numbest = nD;

        int curr = rand() % numbest;

        WTArow[x] = bestd1[curr];
        WTArow[x+1] = bestd2[curr];

        uchar *gtpix1 = &gtIm[i].Pixel(x, y, 0);
        uchar *gtpix2 = &gtIm[i].Pixel(x+1, y, 0);

        loglikelihood += wpairwiseProbV[tidcurr][gtpix1[0]*nD+gtpix2[0]];

      }

    }
  }

  delete [] bestd1;
  delete [] bestd2;

}






void logisticRegressionPairwise(int nD, int numTrainingPats, std::vector <CByteImage> &wta, fvec thetaV, std::vector<int> gradThresh, fvec thetaZ, std::vector<int> gradThreshZ, std::vector<CByteImage> im1grad, std::vector<int> startSliceNo, std::vector<int> endSliceNo, int index, double &loglikelihood, std::vector <CByteImage> gtIm, std::vector <std::vector <CByteImage> > &dirDispXYZ)
{

  checkForFloatCosts();

  int nV = thetaV.size();
  int nG = gradThresh.size()+1;

  int nZ = thetaZ.size();
  int nGz = gradThreshZ.size()+1;


  // for now using only the widthwise pixel pairwise terms
  if (wpairwiseProbV.size() == 0)
    for (int ii=0; ii<nG; ii++)
      wpairwiseProbV.push_back(new MRF::CostVal[nD*nD]);

  if (wpairwiseProbZ.size() == 0)
    for (int ii=0; ii<nGz; ii++)
      wpairwiseProbZ.push_back(new MRF::CostVal[nD*nD]);


  MRF::CostVal badcost = 1e20;

  for(int j=0; j<nV; j++)
	  badcost += thetaV[j];

//  int thetaid;

  for (int tid = 0; tid < nG; tid++) {
    std::vector<double> sumP(nD*nD);
    int ih = 0;
    //    double sumP = 0.0;
    for (int d1 = 0; d1 < nD; d1++) {
      for (int d2 = d1; d2 < nD; d2++) {

        // d1 will always be less than or equal to d2

        wpairwiseProbV[tid][d1*nD+d2] = - thetaV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + tid];  // the -ve sign because we are modeling the exponent as E_s or cost
        wpairwiseProbV[tid][d2*nD+d1] = wpairwiseProbV[tid][d1*nD+d2];

        if (d1==d2)
          sumP[ih] = wpairwiseProbV[tid][d2*nD+d1];
        else {
          sumP[ih] = wpairwiseProbV[tid][d2*nD+d1];
          ih++;
          sumP[ih] = wpairwiseProbV[tid][d2*nD+d1];
        }
        ih++;

      }
    }

    double logSumP  = getLogSumExp(sumP);


    for (int d1 = 0; d1 < nD; d1++) {
      for (int d2 = d1; d2 < nD; d2++) {
        wpairwiseProbV[tid][d1*nD+d2] = wpairwiseProbV[tid][d1*nD+d2] - logSumP;
        wpairwiseProbV[tid][d2*nD+d1] = wpairwiseProbV[tid][d1*nD+d2];
      }
    }
  }


  badcost = 1e20;

  for(int j=0; j<nZ; j++)
	  badcost += thetaZ[j];

  for (int tid = 0; tid < nGz; tid++) {
    std::vector<double> sumP(nD*nD);
    int ih = 0;
    //    double sumP = 0.0;
    for (int d1 = 0; d1 < nD; d1++) {
      for (int d2 = d1; d2 < nD; d2++) {

        // d1 will always be less than or equal to d2

        wpairwiseProbZ[tid][d1*nD+d2] = - thetaZ[(d1*nD + d2 - (d1*(d1+1))/2)*nGz + tid];  // the -ve sign because we are modeling the exponent as E_s or cost
        wpairwiseProbZ[tid][d2*nD+d1] = wpairwiseProbZ[tid][d1*nD+d2];

        if (d1==d2)
          sumP[ih] = wpairwiseProbZ[tid][d2*nD+d1];
        else {
          sumP[ih] = wpairwiseProbZ[tid][d2*nD+d1];
          ih++;
          sumP[ih] = wpairwiseProbZ[tid][d2*nD+d1];
        }
        ih++;

      }
    }

    double logSumP  = getLogSumExp(sumP);


    for (int d1 = 0; d1 < nD; d1++) {
      for (int d2 = d1; d2 < nD; d2++) {
        wpairwiseProbZ[tid][d1*nD+d2] = wpairwiseProbZ[tid][d1*nD+d2] - logSumP;
        wpairwiseProbZ[tid][d2*nD+d1] = wpairwiseProbZ[tid][d1*nD+d2];
      }
    }
  }



  int depth = gtIm.size();
  CShape sh = gtIm[0].Shape();
  int width = sh.width, height = sh.height, nB = sh.nBands;
  int tempglobalNpixels = depth * width * height;
  int nColors = __min(3, nB);

//  int zslice = 0;
//  zslice = startSliceNo[index];

  int dsiIndex = 0;

  int* bestd1 = new int[nD*nD];
  int* bestd2 = new int[nD*nD];

  // x0-x1, x2-x3, x4-x5 etc

  for(unsigned int i = 0; i < gtIm.size(); ++i){

    DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);
    sh.nBands = 1;
    dirDispXYZ[0][i].ReAllocate(sh);

    for (int y = 0; y < height; y++) {
      uchar *WTArow = &dirDispXYZ[0][i].Pixel(0, y, 0);
      for (int x = 0; x < width-1; x=x+2) {

        //unsigned short pix1 = im1[i].Pixel(x, y, 0);
        // unsigned short pix2 = im1[i].Pixel(x+1, y, 0);

        uchar *grad = &im1grad[i].Pixel(x, y, 0);

        int tidcurr = grad[0];

        // std::cout<< tidcurr << "\n";

        int numbest = 0;
        MRF::CostVal bestval = badcost;

        for (int d1 = 0; d1 < nD; d1++) {
          for (int d2 = d1; d2 < nD; d2++) {

            MRF::CostVal dsiValue = wpairwiseProbV[tidcurr][d1*nD+d2];
            // dsiValue is now the probability,
            // and we want the one with max prob.
            if (dsiValue > bestval || numbest == 0) {
              bestval = dsiValue;
              bestd1[0] = d1;
              bestd2[0] = d2;
              numbest = 1;
            } else if (dsiValue == bestval) {
              bestd1[numbest] = d1;
              bestd2[numbest] = d2;
              numbest++;
            }
          }
        }

        // selecting the best candidates from the list uniformly
        if (numbest == 0)
        	numbest = nD;

        int curr = rand() % numbest;

        WTArow[x] = bestd1[curr];
        WTArow[x+1] = bestd2[curr];

        uchar *gtpix1 = &gtIm[i].Pixel(x, y, 0);
        uchar *gtpix2 = &gtIm[i].Pixel(x+1, y, 0);

        loglikelihood += wpairwiseProbV[tidcurr][gtpix1[0]*nD+gtpix2[0]];

      }

    }
  }

  // x1-x2, x3-x4, x5-x6 etc
  for(unsigned int i = 0; i < gtIm.size(); ++i){

    DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);
    sh.nBands = 1;
    dirDispXYZ[1][i].ReAllocate(sh);

    for (int y = 0; y < height; y++) {
      uchar *WTArow = &dirDispXYZ[1][i].Pixel(0, y, 0);
      for (int x = 1; x < width-1; x=x+2) {

        //unsigned short pix1 = im1[i].Pixel(x, y, 0);
        // unsigned short pix2 = im1[i].Pixel(x+1, y, 0);

        uchar *grad = &im1grad[i].Pixel(x, y, 0);

        int tidcurr = grad[0];

        // std::cout<< tidcurr << "\n";

        int numbest = 0;
        MRF::CostVal bestval = badcost;

        for (int d1 = 0; d1 < nD; d1++) {
          for (int d2 = d1; d2 < nD; d2++) {

            MRF::CostVal dsiValue = wpairwiseProbV[tidcurr][d1*nD+d2];
            // dsiValue is now the probability,
            // and we want the one with max prob.
            if (dsiValue > bestval || numbest == 0) {
              bestval = dsiValue;
              bestd1[0] = d1;
              bestd2[0] = d2;
              numbest = 1;
            } else if (dsiValue == bestval) {
              bestd1[numbest] = d1;
              bestd2[numbest] = d2;
              numbest++;
            }
          }
        }

        // selecting the best candidates from the list uniformly
        if (numbest == 0)
        	numbest = nD;

        int curr = rand() % numbest;

        WTArow[x] = bestd1[curr];
        WTArow[x+1] = bestd2[curr];

        uchar *gtpix1 = &gtIm[i].Pixel(x, y, 0);
        uchar *gtpix2 = &gtIm[i].Pixel(x+1, y, 0);

        loglikelihood += wpairwiseProbV[tidcurr][gtpix1[0]*nD+gtpix2[0]];

      }

    }
  }

  // y0-y1, y2-y3, etc
  for(unsigned int i = 0; i < gtIm.size(); ++i){

    DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);
    sh.nBands = 1;
    dirDispXYZ[2][i].ReAllocate(sh);

    for (int y = 0; y < height-1; y=y+2) {
      uchar *WTArow = &dirDispXYZ[2][i].Pixel(0, y, 0);
      uchar *WTArowy = &dirDispXYZ[2][i].Pixel(0, y+1, 0);
      for (int x = 0; x < width; x++) {

        //unsigned short pix1 = im1[i].Pixel(x, y, 0);
        // unsigned short pix2 = im1[i].Pixel(x+1, y, 0);

        uchar *grad = &im1grad[i].Pixel(x, y, 0);

        int tidcurr = grad[1];

        // std::cout<< tidcurr << "\n";

        int numbest = 0;
        MRF::CostVal bestval = badcost;

        for (int d1 = 0; d1 < nD; d1++) {
          for (int d2 = d1; d2 < nD; d2++) {

            MRF::CostVal dsiValue = wpairwiseProbV[tidcurr][d1*nD+d2];
            // dsiValue is now the probability,
            // and we want the one with max prob.
            if (dsiValue > bestval || numbest == 0) {
              bestval = dsiValue;
              bestd1[0] = d1;
              bestd2[0] = d2;
              numbest = 1;
            } else if (dsiValue == bestval) {
              bestd1[numbest] = d1;
              bestd2[numbest] = d2;
              numbest++;
            }
          }
        }

        // selecting the best candidates from the list uniformly

        if (numbest == 0)
        	numbest = nD;
        int curr = rand() % numbest;

        WTArow[x] = bestd1[curr];
        WTArowy[x] = bestd2[curr];

        uchar *gtpix1 = &gtIm[i].Pixel(x, y, 0);
        uchar *gtpix2 = &gtIm[i].Pixel(x, y+1, 0);

        loglikelihood += wpairwiseProbV[tidcurr][gtpix1[0]*nD+gtpix2[0]];

      }

    }
  }

  // y1-y2, y3-y4, etc
  for(unsigned int i = 0; i < gtIm.size(); ++i){

    DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);
    sh.nBands = 1;
    dirDispXYZ[3][i].ReAllocate(sh);

    for (int y = 1; y < height-1; y+=2) {
      uchar *WTArow = &dirDispXYZ[3][i].Pixel(0, y, 0);
      uchar *WTArowy = &dirDispXYZ[3][i].Pixel(0, y+1, 0);
      for (int x = 0; x < width; x++) {

        //unsigned short pix1 = im1[i].Pixel(x, y, 0);
        // unsigned short pix2 = im1[i].Pixel(x+1, y, 0);

        uchar *grad = &im1grad[i].Pixel(x, y, 0);

        int tidcurr = grad[1];

        // std::cout<< tidcurr << "\n";

        int numbest = 0;
        MRF::CostVal bestval = badcost;

        for (int d1 = 0; d1 < nD; d1++) {
          for (int d2 = d1; d2 < nD; d2++) {

            MRF::CostVal dsiValue = wpairwiseProbV[tidcurr][d1*nD+d2];
            // dsiValue is now the probability,
            // and we want the one with max prob.
            if (dsiValue > bestval || numbest == 0) {
              bestval = dsiValue;
              bestd1[0] = d1;
              bestd2[0] = d2;
              numbest = 1;
            } else if (dsiValue == bestval) {
              bestd1[numbest] = d1;
              bestd2[numbest] = d2;
              numbest++;
            }
          }
        }

        // selecting the best candidates from the list uniformly

        if (numbest == 0)
        	numbest = nD;

        int curr = rand() % numbest;

        WTArow[x] = bestd1[curr];
        WTArowy[x] = bestd2[curr];

        uchar *gtpix1 = &gtIm[i].Pixel(x, y, 0);
        uchar *gtpix2 = &gtIm[i].Pixel(x, y+1, 0);

        loglikelihood += wpairwiseProbV[tidcurr][gtpix1[0]*nD+gtpix2[0]];

      }

    }
  }

// z0-z1, z2-z3 etc
  for(unsigned int i = 0; i < gtIm.size()-1; i=i+2){

    DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);
    sh.nBands = 1;
    dirDispXYZ[4][i].ReAllocate(sh);

    for (int y = 0; y < height; y++) {
      uchar *WTArow = &dirDispXYZ[4][i].Pixel(0, y, 0);
      uchar *WTArowz = &dirDispXYZ[4][i+1].Pixel(0, y, 0);
      for (int x = 0; x < width; x++) {

        //unsigned short pix1 = im1[i].Pixel(x, y, 0);
        // unsigned short pix2 = im1[i].Pixel(x+1, y, 0);

        uchar *grad = &im1grad[i].Pixel(x, y, 0);

        int tidcurr = grad[2];

        // std::cout<< tidcurr << "\n";

        int numbest = 0;
        MRF::CostVal bestval = badcost;

        for (int d1 = 0; d1 < nD; d1++) {
          for (int d2 = d1; d2 < nD; d2++) {

            MRF::CostVal dsiValue = wpairwiseProbZ[tidcurr][d1*nD+d2];
            // dsiValue is now the probability,
            // and we want the one with max prob.
            if (dsiValue > bestval || numbest == 0) {
              bestval = dsiValue;
              bestd1[0] = d1;
              bestd2[0] = d2;
              numbest = 1;
            } else if (dsiValue == bestval) {
              bestd1[numbest] = d1;
              bestd2[numbest] = d2;
              numbest++;
            }
          }
        }

        // selecting the best candidates from the list uniformly
        if (numbest == 0)
        	numbest = nD;
        int curr = rand() % numbest;

        WTArow[x] = bestd1[curr];
        WTArowz[x] = bestd2[curr];

        uchar *gtpix1 = &gtIm[i].Pixel(x, y, 0);
        uchar *gtpix2 = &gtIm[i+1].Pixel(x, y, 0);

        loglikelihood += wpairwiseProbZ[tidcurr][gtpix1[0]*nD+gtpix2[0]];

      }

    }
  }

// z1-z2, z3-z4 etc
  for(unsigned int i = 1; i < gtIm.size() - 1; i=i+2){

    DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);
    sh.nBands = 1;
    dirDispXYZ[5][i].ReAllocate(sh);

    for (int y = 0; y < height; y++) {
      uchar *WTArow = &dirDispXYZ[5][i].Pixel(0, y, 0);
      uchar *WTArowz = &dirDispXYZ[5][i+1].Pixel(0, y, 0);
      for (int x = 0; x < width; x++) {

        //unsigned short pix1 = im1[i].Pixel(x, y, 0);
        // unsigned short pix2 = im1[i].Pixel(x+1, y, 0);

        uchar *grad = &im1grad[i].Pixel(x, y, 0);

        int tidcurr = grad[2];

        // std::cout<< tidcurr << "\n";

        int numbest = 0;
        MRF::CostVal bestval = badcost;

        for (int d1 = 0; d1 < nD; d1++) {
          for (int d2 = d1; d2 < nD; d2++) {

            MRF::CostVal dsiValue = wpairwiseProbZ[tidcurr][d1*nD+d2];
            // dsiValue is now the probability,
            // and we want the one with max prob.
            if (dsiValue > bestval || numbest == 0) {
              bestval = dsiValue;
              bestd1[0] = d1;
              bestd2[0] = d2;
              numbest = 1;
            } else if (dsiValue == bestval) {
              bestd1[numbest] = d1;
              bestd2[numbest] = d2;
              numbest++;
            }
          }
        }

        // selecting the best candidates from the list uniformly
        if (numbest == 0)
        	numbest = nD;
        int curr = rand() % numbest;

        WTArow[x] = bestd1[curr];
        WTArowz[x] = bestd2[curr];

        uchar *gtpix1 = &gtIm[i].Pixel(x, y, 0);
        uchar *gtpix2 = &gtIm[i+1].Pixel(x, y, 0);

        loglikelihood += wpairwiseProbZ[tidcurr][gtpix1[0]*nD+gtpix2[0]];

      }

    }
  }


  delete [] bestd1;
  delete [] bestd2;

  std::vector <int> vote;
  vote.resize(nD);

// voting from the 6 output values to make a decision
  for(unsigned int i = 0; i < gtIm.size(); ++i){

    sh.nBands = 1;
    wta[i].ReAllocate(sh);

    for (int y = 1; y < height-1; y++) {
      for (int x = 1; x < width-1; x++) {

    	  for (int kk=0; kk<nD; kk++)
    		  vote[kk] = 0;

    	  for (int jj=0; jj<4; jj++)
    	  {
    		  uchar *WTArow = &dirDispXYZ[jj][i].Pixel(x, y, 0);

    		  vote[WTArow[0]] ++;
    	  }

    	  // 5th
    	  if (gtIm.size() % 2 == 0) // even
    	  {
    		  uchar *WTArow = &dirDispXYZ[4][i].Pixel(x, y, 0);
    		  vote[WTArow[0]] ++;

    	  } else { // odd

    		  if (i==gtIm.size() - 1)
    		  {
    			  // skip
    		  } else {
        		  uchar *WTArow = &dirDispXYZ[4][i].Pixel(x, y, 0);
        		  vote[WTArow[0]] ++;
    		  }
    	  }

    	  // 6th
    	  if (i!=0)
    	  {
			  if (gtIm.size() % 2 == 0) // even
			  {
				  if (i==gtIm.size() - 1)
				  {
				      			  // skip
				  } else
				  {
	        		  uchar *WTArow = &dirDispXYZ[5][i].Pixel(x, y, 0);
	        		  vote[WTArow[0]] ++;
				  }

			  } else { // odd

				  uchar *WTArow = &dirDispXYZ[5][i].Pixel(x, y, 0);
				  vote[WTArow[0]] ++;

			  }
    	  }

    	  // this code is biased to take the first in case there is a tie of best points
    	  int maxvote = 0, maxvoteval = vote[0];

    	  for (int kk=1; kk<nD; kk++)
    	  {
    		  if (vote[kk]>maxvoteval)
    		  {
    			  maxvoteval = vote[kk];
    			  maxvote = kk;
    		  }
    	  }

    	  uchar *WTArow = &wta[i].Pixel(x, y, 0);
    	  WTArow[0] = maxvote;


      }
    }
  }



}



















// pairwise with datacost

void logisticRegressionPairwiseD(int genparam, std::vector <CImage2> im1, std::vector <CByteImage>  hogIm,
      int nD, int numTrainingPats, std::vector <CByteImage> &wta,      // winner-take-all disparities
      fvec thetaU, fvec thetaA, fvec thetaH, fvec thetaL, int featureCode, std::vector <std::vector <CByteImage> > appdirImage,
      Appearance** appclass, LocationMV** LocationMVclass, AppearanceMV** AppMVclass,
      std::vector <std::vector  <matrixB<double> > >  appdirMVProb, intClass **intObj, std::vector <std::vector <CImage2> > locdirImage,
      std::vector <std::vector  <matrixB<double> > > hogdirMVProb, intensityNB* intNBclass, 
      std::vector<int> startSliceNo, std::vector<int> endSliceNo, int nbinsNB, int index, double &loglikelihood, 
      std::vector <CByteImage> gtIm, fvec thetaV, std::vector<int> gradThresh, std::vector<CByteImage> im1grad)
{

  checkForFloatCosts();

  // dsiArrayV and dsiArrayVGen will be used to store the probability values in logistic regression instead 
  // will need to see how to process dsiArrayVGen (derive equations etc)
  if (genparam==0) {
    if (dsiArrayV.empty() == true)
      throw CError("call initializeDataCost first");
  } else if (genparam==1) {
    if (dsiArrayVGen.empty() == true)
      throw CError("call initializeDataCostLoc first");
  }

  // need to make sure reInitializeDataCostVectorPairwise() has been called.


  int intensity = (int)featureCode/10000;
  int hog = (int) (featureCode/1000)%10;
  int app = (int) (featureCode/100)%10;
  int loc = (int) (featureCode/10)%10;


  int nU = thetaU.size();
  int nbins = nU/nD;

  // assuming that we have png data that goes upto 65535.
  // since i have the data in the range of 850 to 1250 which is shifted to 0 to 256
  double binner = 256.0/nbins;
  // so for any point. divide it by binner to get bin number (floor int value) it is in
  // ranging from bin0 to bin{nbins-1}
  // since >=256 will be special case, put it in the last bin if there.

  int nA = thetaA.size();

  int nT;

  nT = nA/nD; //this gives centers per class
  // this will work for app==2 and app==1 because in app==2, we have nD*the number of total clusters as the number of parameters
  // while in app==1 it is nD*number of clusters of that class (only that clusters are also = number of params)

  int nH = thetaH.size(); // number of HoG parameters
  int nHpC = nH/nD; // HoG vocabulary size

  int nL = thetaL.size(); // number of location parameters
  int nLpC = nL/nD;  // this variable will be useful only for cube setting (loc==1) or loc==9


  int nV = thetaV.size();
  int nG = gradThresh.size()+1;


  MRF::CostVal badcost = 1e20;

  for(int j=0; j<nU; j++)
    badcost += thetaU[j];

  if (nA > 0)
    for(int j=0; j<nA; j++)
      badcost += thetaA[j];

  if (nH > 0)
    for(int j=0; j<nH; j++)
      badcost += thetaH[j];

  if (nL > 0)
    for(int j=0; j<nL; j++)
      badcost += thetaL[j];

  for (int k=0; k<nV; k++) 
    for(int j=0; j<nV; j++)
      badcost += thetaV[j];

  int thetaid;

  if (genparam==1) {
    int numGens=0;
    if (intensity==3 || intensity==6)
      numGens++;
    if (loc==3 || loc==8 || loc==9)
      numGens++;
    if (app==3)
      numGens++;
    if (hog==3)
      numGens++;


    //WORK-THIS
    // assuming worst cost after taking log is -log(10^-99)

    badcost += numGens*(99*3);  //that 3 accounts for 2.7 value of e (natural log)

  }

  //  int* bestd = new int[nD];
  int* bestd1 = new int[nD*nD];
  int* bestd2 = new int[nD*nD];



  // for (unsigned int j = 0; j < im1.size(); ++j) {

  int depth = im1.size();
  CShape sh = im1[0].Shape();
  int width = sh.width, height = sh.height, nB = sh.nBands;
  int tempglobalNpixels = depth * width * height;
  int nColors = __min(3, nB);

  int zslice = 0;
  zslice = startSliceNo[index];

  int dsiIndex = 0;
  for(unsigned int i = 0; i < im1.size(); ++i){

    DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);
    sh.nBands = 1;
    wta[i].ReAllocate(sh);

    for (int y = 0; y < height; y++) {

      uchar *WTArow = &wta[i].Pixel(0, y, 0);

      for (int x = 0; x < width-1; x=x+2) {

        unsigned short pix1p = im1[i].Pixel(x, y, 0);
        unsigned short pix1q = im1[i].Pixel(x+1, y, 0);

        uchar *hogvalp, *hogvalq;
        if (nH > 0) {
          hogvalp = &hogIm[i].Pixel(x, y, 0);
          hogvalq = &hogIm[i].Pixel(x+1, y, 0);
        }

        //        uchar *gtpixp = &gtIm[i].Pixel(x, y, 0);
        //uchar *gtpixq = &gtIm[i].Pixel(x+1, y, 0);

        uchar *grad = &im1grad[i].Pixel(x, y, 0);

        int tidcurr = grad[0];



        //int bestd = 0;
        int numbest = 0;
        MRF::CostVal bestval = badcost;

        std::vector<double> dsiValueSum(nD*nD);

        for (int d1 = 0; d1 < nD; d1++) {
          for (int d2 = 0; d2 < nD; d2++) {

            MRF::CostVal dsiValuep = 0, dsiValueq = 0, dsiValuepq = 0;

            if (nU>0 && intensity==1)  // intensity as histogram bins
              {
                dsiValuep += getCostIntensityDiscBins(pix1p, nbins, d1, thetaU, thetaid);
                dsiValueq += getCostIntensityDiscBins(pix1q, nbins, d2, thetaU, thetaid);

              } else if (intensity==3 && genparam==1) // generative setting with Gaussians
              {
                dsiValuep += getCostIntensityGenGaussian(intObj[d1], pix1p);
                dsiValueq += getCostIntensityGenGaussian(intObj[d2], pix1q);
              } else if (intensity == 6 && genparam==1) // using NBayes and bins
              {
                dsiValuep += getCostIntensityGenNB(pix1p, nbins, d1, intNBclass);
                dsiValueq += getCostIntensityGenNB(pix1q, nbins, d2, intNBclass);
              }


            if (nA>0 && app==1) // we have appearance features
              {
                // different features for different classes
					   
                dsiValuep += getCostAppDiscPatch((int)appdirImage[i][d1].Pixel(x, y, 0), appclass[d1]->getPatchSize(), x, y, width, height, thetaA, d1, nT, thetaid);
                dsiValueq += getCostAppDiscPatch((int)appdirImage[i][d2].Pixel(x+1, y, 0), appclass[d2]->getPatchSize(), x+1, y, width, height, thetaA, d2, nT, thetaid);

              } else if (nA>0 && app==2) // we have appearance features
              {
                // common features for all classes
					    
                dsiValuep += getCostAppDiscPatch((int)appdirImage[i][0].Pixel(x, y, 0), appclass[0]->getPatchSize(), x, y, width, height, thetaA, d1, nT, thetaid);
                dsiValueq += getCostAppDiscPatch((int)appdirImage[i][0].Pixel(x+1, y, 0), appclass[0]->getPatchSize(), x+1, y, width, height, thetaA, d2, nT, thetaid);

              }
            else if (app==3 && genparam==1) // generative setting
              {
                dsiValuep += getCostAppGenPatch((appdirMVProb[i][d1])(y,x), AppMVclass[d1]->getDimension(), x, y, width, height);
                dsiValueq += getCostAppGenPatch((appdirMVProb[i][d2])(y,x+1), AppMVclass[d2]->getDimension(), x+1, y, width, height);
              }


            if (nH > 0 && hog==1) {
					    
              dsiValuep += getCostHoGDiscBins((int) *hogvalp, d1, nHpC, thetaH, thetaid);
              dsiValueq += getCostHoGDiscBins((int) *hogvalq, d2, nHpC, thetaH, thetaid);

            } else if (hog==3 && genparam==1) // generative setting
              {

                int startPatch = 8; // index starting from 0
                // this variable indicates where the Hog descriptors are defined
                // since i have generated them using the entire patient body, don't need to
                // worry about the z axis and can use all slices here.
                //WORK-THIS need to change above stuff for flawless working

                dsiValuep += getCostHoGGenBins((hogdirMVProb[i][d1])(y,x), x, y, width, height, startPatch);
                dsiValueq += getCostHoGGenBins((hogdirMVProb[i][d2])(y,x+1), x+1, y, width, height, startPatch);

              }

            if (nL>0 && loc==1) { // CRF cube code

              // this function is completely out-dated
              // and for older series of data, will need to modify significantly for use
					    
              dsiValuep += getCostLocationDiscCube((unsigned short) locdirImage[d1][i].Pixel(x,y,0), x, y, i, loccubex, loccubey, loccubez, numTrainingPats, nLpC, thetaL, width, height, depth, d1, thetaid);
              dsiValueq += getCostLocationDiscCube((unsigned short) locdirImage[d2][i].Pixel(x+1,y,0), x+1, y, i, loccubex, loccubey, loccubez, numTrainingPats, nLpC, thetaL, width, height, depth, d2, thetaid);


            } else if (nL>0 && (loc==6 || loc==8) )  //CRF hard or soft assignment
              {
					    
                dsiValuep += getCostLocationDiscGaussian(LocationMVclass, i, zslice, x, y, thetaL, d1, loc, genparam, thetaid);
                dsiValueq += getCostLocationDiscGaussian(LocationMVclass, i, zslice, x+1, y, thetaL, d2, loc, genparam, thetaid);


              } else if (nL>0 && (loc==7 || loc==9) )  //CRF soft or hard assignment redefined
              {
					    
                dsiValuep += getCostLocationDiscGaussian2(LocationMVclass, i, zslice, x, y, thetaL, d1, loc, genparam, nD, thetaid);
                dsiValueq += getCostLocationDiscGaussian2(LocationMVclass, i, zslice, x+1, y, thetaL, d2, loc, genparam, nD, thetaid);


              } else if (loc==3 && genparam==1) // generative setting
              {

                dsiValuep += getCostLocationGenGaussian(LocationMVclass[d1], i, zslice, x, y);
                dsiValueq += getCostLocationGenGaussian(LocationMVclass[d2], i, zslice, x+1, y);

              }

            if (d2>d1)
              dsiValuepq =  thetaV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + tidcurr];
            else
              dsiValuepq =  thetaV[(d2*nD + d1 - (d2*(d2+1))/2)*nG + tidcurr];

            // The cost of pixel p and label l is stored at dsiArray[p*nLabels+l]
            // Since in logistic regression, I am not using E_d ie e(-E_d), I negate the values of dsiValue so i can work on prob()
            // directly.

            double dsiValue =  (double)dsiValuep + (double)dsiValueq + (double)dsiValuepq;

            if (genparam==0)
              (dsiArrayV[index])[dsiIndex++] = - dsiValue;
            else if (genparam==1)
              (dsiArrayVGen[index])[dsiIndex++] = - dsiValue;

            dsiValueSum[d1*nD + d2] = - dsiValue;


            // dsiValue here is the cost, so we want the minimum
            if (dsiValue < bestval) {
              bestval = dsiValue;
              bestd1[0] = d1;
              bestd2[0] = d2;
              numbest = 1;
            } else if (dsiValue == bestval) {
              bestd1[numbest] = d1;
              bestd2[numbest] = d2;
              numbest++;
            }

            if (pix1p >= 256 && pix1q >= 256) {
              bestd1[0] = 0; // assume 0 is background
              bestd2[0] = 0;
              numbest = 1;
            }

            // these 2 ifs are not the best way to force only one node to have some disp value
            if (pix1p >= 256) {
              bestd1[0] = 0; // assume 0 is background
              numbest = 1;
            }
            if (pix1q >= 256) {
              bestd2[0] = 0; // assume 0 is background
              numbest = 1;
            }

          }
        }
        if (numbest == 0)
        	numbest = nD;
        int curr = rand() % numbest;
        WTArow[x] = bestd1[curr];
        WTArow[x+1] = bestd2[curr];
                 
        uchar *gtpix1 = &gtIm[i].Pixel(x, y, 0);
        uchar *gtpix2 = &gtIm[i].Pixel(x+1, y, 0);

        double dsiValueLogSum  = getLogSumExp(dsiValueSum);

        for (int d = nD*nD; d > 0; d--) {
          if (genparam==0)
            (dsiArrayV[index])[dsiIndex-d] = (double)(dsiArrayV[index])[dsiIndex-d] - dsiValueLogSum;
          else if (genparam==1)
            (dsiArrayVGen[index])[dsiIndex-d] = (double)(dsiArrayVGen[index])[dsiIndex-d] - dsiValueLogSum;                    
        }

        
        // dsiIndex is already pointing to the new pixel and output level index and so need to rewind
        int tempid = dsiIndex-nD*nD+gtpix1[0]*nD+gtpix2[0]; //bestd[curr];
        if (genparam == 0) {
          loglikelihood += (dsiArrayV[index])[tempid];
          // std::cout << " logSum = "<< log(dsiValueSum) << "  LL =  " << (dsiArrayV[index])[tempid] << "\n";
        }
        else if (genparam==1) {
          loglikelihood += (dsiArrayVGen[index])[tempid];
        }

      }
    }
  }

  //  std::cout << "  LL =  " << loglikelihood  << "\n";
         
  delete [] bestd1;
  delete [] bestd2;

}


void reInitializeDataCostVectorPairwise(int numPatients, int genparam, int nD, int power)
{

  dsiArrayV.resize(0);
  dsiArrayVGen.resize(0);

  int nDp = 1;
  for (int i = 0; i < power; i++)
    nDp *= nDp*nD;

  for (unsigned int j = 0; j < numPatients; ++j) {

    if (genparam==0)
      dsiArrayV.push_back(new MRF::CostVal[int(globalNpixels[j]/2 + 1) * nDp]);  // since we have pairs now, we have one number for 2 pixels
    else if (genparam==1)
      dsiArrayVGen.push_back(new MRF::CostVal[int(globalNpixels[j]/2 + 1) * nDp]);

  }

}





void logisticRegressionPairwiseD4(int genparam, std::vector <CImage2> im1, std::vector <CByteImage>  hogIm,
      int nD, int numTrainingPats, std::vector <CByteImage> &wta,      // winner-take-all disparities
      fvec thetaU, fvec thetaA, fvec thetaH, fvec thetaL, int featureCode, std::vector <std::vector <CByteImage> > appdirImage,
      Appearance** appclass, LocationMV** LocationMVclass, AppearanceMV** AppMVclass,
      std::vector <std::vector  <matrixB<double> > >  appdirMVProb, intClass **intObj, std::vector <std::vector <CImage2> > locdirImage,
      std::vector <std::vector  <matrixB<double> > > hogdirMVProb, intensityNB* intNBclass, 
      std::vector<int> startSliceNo, std::vector<int> endSliceNo, int nbinsNB, int index, double &loglikelihood, 
      std::vector <CByteImage> gtIm, fvec thetaV, std::vector<int> gradThresh, std::vector<CByteImage> im1grad)
{

  checkForFloatCosts();

  // dsiArrayV and dsiArrayVGen will be used to store the probability values in logistic regression instead 
  // will need to see how to process dsiArrayVGen (derive equations etc)
  if (genparam==0) {
    if (dsiArrayV.empty() == true)
      throw CError("call initializeDataCost first");
  } else if (genparam==1) {
    if (dsiArrayVGen.empty() == true)
      throw CError("call initializeDataCostLoc first");
  }

  // need to make sure reInitializeDataCostVectorPairwise() has been called.


  int intensity = (int)featureCode/10000;
  int hog = (int) (featureCode/1000)%10;
  int app = (int) (featureCode/100)%10;
  int loc = (int) (featureCode/10)%10;


  int nU = thetaU.size();
  int nbins = nU/nD;

  // assuming that we have png data that goes upto 65535.
  // since i have the data in the range of 850 to 1250 which is shifted to 0 to 256
  double binner = 256.0/nbins;
  // so for any point. divide it by binner to get bin number (floor int value) it is in
  // ranging from bin0 to bin{nbins-1}
  // since >=256 will be special case, put it in the last bin if there.

  int nA = thetaA.size();

  int nT;

  nT = nA/nD; //this gives centers per class
  // this will work for app==2 and app==1 because in app==2, we have nD*the number of total clusters as the number of parameters
  // while in app==1 it is nD*number of clusters of that class (only that clusters are also = number of params)

  int nH = thetaH.size(); // number of HoG parameters
  int nHpC = nH/nD; // HoG vocabulary size

  int nL = thetaL.size(); // number of location parameters
  int nLpC = nL/nD;  // this variable will be useful only for cube setting (loc==1) or loc==9


  int nV = thetaV.size();
  int nG = gradThresh.size()+1;


  MRF::CostVal badcost = 1e20;

  for(int j=0; j<nU; j++)
    badcost += thetaU[j];

  if (nA > 0)
    for(int j=0; j<nA; j++)
      badcost += thetaA[j];

  if (nH > 0)
    for(int j=0; j<nH; j++)
      badcost += thetaH[j];

  if (nL > 0)
    for(int j=0; j<nL; j++)
      badcost += thetaL[j];

  for (int k=0; k<nV; k++) 
    for(int j=0; j<nV; j++)
      badcost += thetaV[j];

  int thetaid;

  if (genparam==1) {
    int numGens=0;
    if (intensity==3 || intensity==6)
      numGens++;
    if (loc==3 || loc==8 || loc==9)
      numGens++;
    if (app==3)
      numGens++;
    if (hog==3)
      numGens++;

    //WORK-THIS
    // assuming worst cost after taking log is -log(10^-99)
    badcost += numGens*(99*3);  //that 3 accounts for 2.7 value of e (natural log)

  }

  //  int* bestd = new int[nD];
  int* bestd1p = new int[nD*nD];
  int* bestd2p = new int[nD*nD];
  int* bestd1q = new int[nD*nD];
  int* bestd2q = new int[nD*nD];


  // for (unsigned int j = 0; j < im1.size(); ++j) {

  int depth = im1.size();
  CShape sh = im1[0].Shape();
  int width = sh.width, height = sh.height, nB = sh.nBands;
  int tempglobalNpixels = depth * width * height;
  int nColors = __min(3, nB);

  int zslice = 0;
  zslice = startSliceNo[index];

  int dsiIndex = 0;
  for(unsigned int i = 0; i < im1.size(); ++i) {

    DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);
    sh.nBands = 1;
    wta[i].ReAllocate(sh);
    int hv = 1;
    if (height%2 == 1) // if it is odd
      hv = 2;

    int wv = 1;
    if (width%2 == 1)
      wv = 2;

    for (int y = 1; y < height-hv; y+=2) {

      uchar *WTArow1 = &wta[i].Pixel(0, y, 0);
      uchar *WTArow2 = &wta[i].Pixel(0, y+1, 0);

      for (int x = 1; x < width-wv; x=x+2) {

        // 1p -- 1q
        // |      |
        // 2p -- 2q

        unsigned short pix1p = im1[i].Pixel(x, y, 0);
        unsigned short pix1q = im1[i].Pixel(x+1, y, 0);
        unsigned short pix2p = im1[i].Pixel(x, y+1, 0);
        unsigned short pix2q = im1[i].Pixel(x+1, y+1, 0);


        uchar *hogval1p, *hogval1q, *hogval2p, *hogval2q;
        if (nH > 0) {
          hogval1p = &hogIm[i].Pixel(x, y, 0);
          hogval1q = &hogIm[i].Pixel(x+1, y, 0);
          hogval2p = &hogIm[i].Pixel(x, y+1, 0);
          hogval2q = &hogIm[i].Pixel(x+1, y+1, 0);      
        }

        //        uchar *gtpixp = &gtIm[i].Pixel(x, y, 0);
        //uchar *gtpixq = &gtIm[i].Pixel(x+1, y, 0);

        uchar *grad1p = &im1grad[i].Pixel(x, y, 0);
        uchar *grad1q = &im1grad[i].Pixel(x+1, y, 0);
        uchar *grad2p = &im1grad[i].Pixel(x, y+1, 0);
        

        // up, left to right
        int tuplr = grad1p[0];
        // left, up to down
        int tlud = grad1p[1];
        // right, up to down
        int trud = grad1q[1];
        // down, left to right
        int tdlr = grad2p[0];



        int numbest = 0;
        MRF::CostVal bestval = badcost;

        std::vector<double> dsiValueSum(nD*nD*nD*nD);

        for (int d1p = 0; d1p < nD; d1p++) {
          for (int d1q = 0; d1q < nD; d1q++) {
            for (int d2p = 0; d2p < nD; d2p++) {
              for (int d2q = 0; d2q < nD; d2q++) {

                MRF::CostVal dsiValue1p = 0, dsiValue1q = 0, dsiValue2p = 0, dsiValue2q = 0;
                MRF::CostVal dsiValue1p1q = 0, dsiValue1p2p = 0, dsiValue1q2q = 0, dsiValue2p2q = 0;

                if (nU>0 && intensity==1)  // intensity as histogram bins
                  {
                    dsiValue1p += getCostIntensityDiscBins(pix1p, nbins, d1p, thetaU, thetaid);
                    dsiValue1q += getCostIntensityDiscBins(pix1q, nbins, d1q, thetaU, thetaid);
                    dsiValue2p += getCostIntensityDiscBins(pix2p, nbins, d2p, thetaU, thetaid);
                    dsiValue2q += getCostIntensityDiscBins(pix2q, nbins, d2q, thetaU, thetaid);

                  } else if (intensity==3 && genparam==1) // generative setting with Gaussians
                  {
                    dsiValue1p += getCostIntensityGenGaussian(intObj[d1p], pix1p);
                    dsiValue1q += getCostIntensityGenGaussian(intObj[d1q], pix1q);
                    dsiValue2p += getCostIntensityGenGaussian(intObj[d2p], pix2p);
                    dsiValue2q += getCostIntensityGenGaussian(intObj[d2q], pix2q);
                  
                  } else if (intensity == 6 && genparam==1) // using NBayes and bins
                  {
                    dsiValue1p += getCostIntensityGenNB(pix1p, nbins, d1p, intNBclass);
                    dsiValue1q += getCostIntensityGenNB(pix1q, nbins, d1q, intNBclass);
                    dsiValue2p += getCostIntensityGenNB(pix2p, nbins, d2p, intNBclass);
                    dsiValue2q += getCostIntensityGenNB(pix2q, nbins, d2q, intNBclass);
                  }


                if (nA>0 && app==1) // we have appearance features
                  {
                    // different features for different classes
					   
                    dsiValue1p += getCostAppDiscPatch((int)appdirImage[i][d1p].Pixel(x, y, 0), appclass[d1p]->getPatchSize(), x, y, width, height, thetaA, d1p, nT, thetaid);
                    dsiValue1q += getCostAppDiscPatch((int)appdirImage[i][d1q].Pixel(x+1, y, 0), appclass[d1q]->getPatchSize(), x+1, y, width, height, thetaA, d1q, nT, thetaid);
                    dsiValue2p += getCostAppDiscPatch((int)appdirImage[i][d2p].Pixel(x, y+1, 0), appclass[d2p]->getPatchSize(), x, y+1, width, height, thetaA, d2p, nT, thetaid);
                    dsiValue2q += getCostAppDiscPatch((int)appdirImage[i][d2q].Pixel(x+1, y+1, 0), appclass[d2q]->getPatchSize(), x+1, y+1, width, height, thetaA, d2q, nT, thetaid);

                  } else if (nA>0 && app==2) // we have appearance features
                  {
                    // common features for all classes
					    
                    dsiValue1p += getCostAppDiscPatch((int)appdirImage[i][0].Pixel(x, y, 0), appclass[0]->getPatchSize(), x, y, width, height, thetaA, d1p, nT, thetaid);
                    dsiValue1q += getCostAppDiscPatch((int)appdirImage[i][0].Pixel(x+1, y, 0), appclass[0]->getPatchSize(), x+1, y, width, height, thetaA, d1q, nT, thetaid);
                    dsiValue2p += getCostAppDiscPatch((int)appdirImage[i][0].Pixel(x, y+1, 0), appclass[0]->getPatchSize(), x, y+1, width, height, thetaA, d2p, nT, thetaid);
                    dsiValue2q += getCostAppDiscPatch((int)appdirImage[i][0].Pixel(x+1, y+1, 0), appclass[0]->getPatchSize(), x+1, y+1, width, height, thetaA, d2q, nT, thetaid);
                
                  }
                else if (app==3 && genparam==1) // generative setting
                  {
                    dsiValue1p += getCostAppGenPatch((appdirMVProb[i][d1p])(y,x), AppMVclass[d1p]->getDimension(), x, y, width, height);
                    dsiValue1q += getCostAppGenPatch((appdirMVProb[i][d1q])(y,x+1), AppMVclass[d1q]->getDimension(), x+1, y, width, height);
                    dsiValue2p += getCostAppGenPatch((appdirMVProb[i][d2p])(y+1,x), AppMVclass[d2p]->getDimension(), x, y+1, width, height);
                    dsiValue2q += getCostAppGenPatch((appdirMVProb[i][d2q])(y+1,x+1), AppMVclass[d2q]->getDimension(), x+1, y+1, width, height);
                  }

                if (nH > 0 && hog==1) {
					    
                  dsiValue1p += getCostHoGDiscBins((int) *hogval1p, d1p, nHpC, thetaH, thetaid);
                  dsiValue1q += getCostHoGDiscBins((int) *hogval1q, d1q, nHpC, thetaH, thetaid);
                  dsiValue2p += getCostHoGDiscBins((int) *hogval2p, d2p, nHpC, thetaH, thetaid);
                  dsiValue2q += getCostHoGDiscBins((int) *hogval2q, d2q, nHpC, thetaH, thetaid);

                } else if (hog==3 && genparam==1) // generative setting
                  {

                    int startPatch = 8; // index starting from 0
                    // this variable indicates where the Hog descriptors are defined
                    // since i have generated them using the entire patient body, don't need to
                    // worry about the z axis and can use all slices here.
                    //WORK-THIS need to change above stuff for flawless working

                    dsiValue1p += getCostHoGGenBins((hogdirMVProb[i][d1p])(y,x), x, y, width, height, startPatch);
                    dsiValue1q += getCostHoGGenBins((hogdirMVProb[i][d1q])(y,x+1), x+1, y, width, height, startPatch);
                    dsiValue2p += getCostHoGGenBins((hogdirMVProb[i][d2p])(y+1,x), x, y+1, width, height, startPatch);
                    dsiValue2q += getCostHoGGenBins((hogdirMVProb[i][d2q])(y+1,x+1), x+1, y+1, width, height, startPatch);

                  }

                if (nL>0 && loc==1) { // CRF cube code

                  // this function is completely out-dated
                  // and for older series of data, will need to modify significantly for use
					    
                  dsiValue1p += getCostLocationDiscCube((unsigned short) locdirImage[d1p][i].Pixel(x,y,0), x, y, i, loccubex, loccubey, loccubez, numTrainingPats, nLpC, thetaL, width, height, depth, d1p, thetaid);
                  dsiValue1q += getCostLocationDiscCube((unsigned short) locdirImage[d1q][i].Pixel(x+1,y,0), x+1, y, i, loccubex, loccubey, loccubez, numTrainingPats, nLpC, thetaL, width, height, depth, d1q, thetaid);
                  dsiValue2p += getCostLocationDiscCube((unsigned short) locdirImage[d2p][i].Pixel(x,y+1,0), x, y+1, i, loccubex, loccubey, loccubez, numTrainingPats, nLpC, thetaL, width, height, depth, d2p, thetaid);
                  dsiValue2q += getCostLocationDiscCube((unsigned short) locdirImage[d2q][i].Pixel(x+1,y+1,0), x+1, y+1, i, loccubex, loccubey, loccubez, numTrainingPats, nLpC, thetaL, width, height, depth, d2q, thetaid);


                } else if (nL>0 && (loc==6 || loc==8) )  //CRF hard or soft assignment
                  {
					    
                    dsiValue1p += getCostLocationDiscGaussian(LocationMVclass, i, zslice, x, y, thetaL, d1p, loc, genparam, thetaid);
                    dsiValue1q += getCostLocationDiscGaussian(LocationMVclass, i, zslice, x+1, y, thetaL, d1q, loc, genparam, thetaid);
                    dsiValue2p += getCostLocationDiscGaussian(LocationMVclass, i, zslice, x, y+1, thetaL, d2p, loc, genparam, thetaid);
                    dsiValue2q += getCostLocationDiscGaussian(LocationMVclass, i, zslice, x+1, y+1, thetaL, d2q, loc, genparam, thetaid);

                  } else if (nL>0 && (loc==7 || loc==9) )  //CRF soft or hard assignment redefined
                  {
					    
                    dsiValue1p += getCostLocationDiscGaussian2(LocationMVclass, i, zslice, x, y, thetaL, d1p, loc, genparam, nD, thetaid);
                    dsiValue1q += getCostLocationDiscGaussian2(LocationMVclass, i, zslice, x+1, y, thetaL, d1q, loc, genparam, nD, thetaid);
                    dsiValue2p += getCostLocationDiscGaussian2(LocationMVclass, i, zslice, x, y+1, thetaL, d2p, loc, genparam, nD, thetaid);
                    dsiValue2q += getCostLocationDiscGaussian2(LocationMVclass, i, zslice, x+1, y+1, thetaL, d2q, loc, genparam, nD, thetaid);

                  } else if (loc==3 && genparam==1) // generative setting
                  {

                    dsiValue1p += getCostLocationGenGaussian(LocationMVclass[d1p], i, zslice, x, y);
                    dsiValue1q += getCostLocationGenGaussian(LocationMVclass[d1q], i, zslice, x+1, y);
                    dsiValue2p += getCostLocationGenGaussian(LocationMVclass[d2p], i, zslice, x, y+1);
                    dsiValue2q += getCostLocationGenGaussian(LocationMVclass[d2q], i, zslice, x+1, y+1);
                  
                  }

                if (d1q>d1p)
                  dsiValue1p1q =  thetaV[(d1p*nD + d1q - (d1p*(d1p+1))/2)*nG + tuplr];
                else
                  dsiValue1p1q =  thetaV[(d1q*nD + d1p - (d1q*(d1q+1))/2)*nG + tuplr];

                if (d2q>d2p)
                  dsiValue2p2q =  thetaV[(d2p*nD + d2q - (d2p*(d2p+1))/2)*nG + tdlr];
                else
                  dsiValue2p2q =  thetaV[(d2q*nD + d2p - (d2q*(d2q+1))/2)*nG + tdlr];

                if (d2p>d1p)
                  dsiValue1p2p =  thetaV[(d1p*nD + d2p - (d1p*(d1p+1))/2)*nG + tlud];
                else
                  dsiValue1p2p =  thetaV[(d2p*nD + d1p - (d2p*(d2p+1))/2)*nG + tlud];

                if (d1q>d2q)
                  dsiValue1q2q =  thetaV[(d2q*nD + d1q - (d2q*(d2q+1))/2)*nG + trud];
                else
                  dsiValue1q2q =  thetaV[(d1q*nD + d2q - (d1q*(d1q+1))/2)*nG + trud];


                // The cost of pixel p and label l is stored at dsiArray[p*nLabels+l]
                // Since in logistic regression, I am not using E_d ie e(-E_d), I negate the values of dsiValue so i can work on prob()
                // directly.

                double dsiValue =  (double)dsiValue1p + (double)dsiValue1q + (double)dsiValue2p + (double)dsiValue2q +  (double)dsiValue1p1q + (double)dsiValue1p2p + (double)dsiValue1q2q + (double)dsiValue2p2q;

                if (genparam==0)
                  (dsiArrayV[index])[dsiIndex++] = - dsiValue;
                else if (genparam==1)
                  (dsiArrayVGen[index])[dsiIndex++] = - dsiValue;

                dsiValueSum[d1p*nD*nD*nD + d1q*nD*nD + d2p*nD + d2q] = - dsiValue;

                // dsiValue here is the cost, so we want the minimum
                if (dsiValue < bestval) {
                  bestval = dsiValue;
                  bestd1p[0] = d1p;
                  bestd2p[0] = d2p;
                  bestd1q[0] = d1q;
                  bestd2q[0] = d2q;
                  numbest = 1;
                } else if (dsiValue == bestval) {
                  bestd1p[numbest] = d1p;
                  bestd2p[numbest] = d2p;
                  bestd1q[numbest] = d1q;
                  bestd2q[numbest] = d2q;
                  numbest++;
                }

                if (pix1p >= 256 && pix1q >= 256 && pix2p >= 256 && pix2q >= 256) {
                  bestd1p[0] = 0; // assume 0 is background
                  bestd2p[0] = 0;
                  bestd1q[0] = 0; // assume 0 is background
                  bestd2q[0] = 0;
                  numbest = 1;
                }

                // these 2 ifs are not the best way to force only one node to have some disp value
                if (pix1p >= 256) {
                  bestd1p[0] = 0; // assume 0 is background
                  numbest = 1;
                }
                if (pix1q >= 256) {
                  bestd1q[0] = 0; // assume 0 is background
                  numbest = 1;
                }
                if (pix2p >= 256) {
                  bestd2p[0] = 0; // assume 0 is background
                  numbest = 1;
                }
                if (pix2q >= 256) {
                  bestd2q[0] = 0; // assume 0 is background
                  numbest = 1;
                }
              }
            }
          }
        }
        if (numbest == 0)
        	numbest = nD;
        int curr = rand() % numbest;
        WTArow1[x] = bestd1p[curr];
        WTArow1[x+1] = bestd1q[curr];
        WTArow2[x] = bestd2p[curr];
        WTArow2[x+1] = bestd2q[curr];

                 
        uchar *gtpix1p = &gtIm[i].Pixel(x, y, 0);
        uchar *gtpix1q = &gtIm[i].Pixel(x+1, y, 0);
        uchar *gtpix2p = &gtIm[i].Pixel(x, y+1, 0);
        uchar *gtpix2q = &gtIm[i].Pixel(x+1, y+1, 0);


        double dsiValueLogSum  = getLogSumExp(dsiValueSum);

        for (int d = nD*nD*nD*nD; d > 0; d--) {
          if (genparam==0)
            (dsiArrayV[index])[dsiIndex-d] = (double)(dsiArrayV[index])[dsiIndex-d] - dsiValueLogSum;
          else if (genparam==1)
            (dsiArrayVGen[index])[dsiIndex-d] = (double)(dsiArrayVGen[index])[dsiIndex-d] - dsiValueLogSum;                    
        }

        
        // dsiIndex is already pointing to the new pixel and output level index and so need to rewind
        int tempid = dsiIndex - nD*nD*nD*nD + gtpix1p[0]*nD*nD*nD + gtpix1q[0]*nD*nD + gtpix2p[0]*nD + gtpix2q[0]; //bestd[curr];
        if (genparam == 0) {
          loglikelihood += (dsiArrayV[index])[tempid];
          // std::cout << " logSum = "<< log(dsiValueSum) << "  LL =  " << (dsiArrayV[index])[tempid] << "\n";
        }
        else if (genparam==1) {
          loglikelihood += (dsiArrayVGen[index])[tempid];
        }

      }
    }
  }

  //  std::cout << "  LL =  " << loglikelihood  << "\n";
         
  delete [] bestd1p;
  delete [] bestd2p;
  delete [] bestd1q;
  delete [] bestd2q;


}




// computeEmpirDistV 


void computeEmpirDistVlpd4(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, std::vector<CByteImage> im1grad, fvec &distV, int nD, int ignoreVal, int gradContext, int logregpaird4, std::vector<int> gradThreshVec, int empir, int index, int paramlreg) 
{

  // have neglected anything about gradContext

  if (logregpaird4 != 1) {
    std::cout << " Error in compute lpd4 \n ";
    exit(1);
  }

  int depth = validdisp.size();
  CShape sh = disp[0].Shape();

  int width = sh.width, height = sh.height;

  int nV = (int)distV.size();
  int nG = gradThreshVec.size() + 1;

  if (nV == 0)
	return; // nothing to do

  //  int nb = 0;
  //float tri = nD*(nD-1)/2;
  //nb = (int) nV / tri;
  
  // nb captures the number of bins per pair of distinct organs (background inclusive)

  // nb useful only of gradContext = 1.

  for (int k=0; k < nV; k++)
	distV[k] = 0;

  int hv = 1;
  if (height%2 == 1) // if it is odd
    hv = 2;

  int wv = 1;
  if (width%2 == 1)
    wv = 2;


  for (int z=0; z < depth - 1; z++) {

    for (int y = 0; y < height-hv; y+=2)
      {

        for (int x = 0; x < width-wv; x+=2) {

          uchar *gtpix1p = &disp[z].Pixel(x, y, 0);
          uchar *gtpix1q = &disp[z].Pixel(x+1, y, 0);
          uchar *gtpix2p = &disp[z].Pixel(x, y+1, 0);
          uchar *gtpix2q = &disp[z].Pixel(x+1, y+1, 0);

          int d1p = (int) gtpix1p[0];
          int d1q = (int) gtpix1q[0];
          int d2p = (int) gtpix2p[0];
          int d2q = (int) gtpix2q[0];


          // 1p -- 1q
          // |      |
          // 2p -- 2q

          //              if (paramlreg == 5  && (y==0 || x==0 || (y==height-1 && height%2==0) || (y>=height-2 && height%2==1) ||  (x==width-1 && width%2==0) || (x>=width-2 && width%2==1)  ))
          //  continue;

          uchar *grad1p = &im1grad[z].Pixel(x, y, 0);
          uchar *grad1q = &im1grad[z].Pixel(x+1, y, 0);
          uchar *grad2p = &im1grad[z].Pixel(x, y+1, 0);
        

          // up, left to right
          int tuplr = grad1p[0];
          // left, up to down
          int tlud = grad1p[1];
          // right, up to down
          int trud = grad1q[1];
          // down, left to right
          int tdlr = grad2p[0];

          int idx1p1q, idx2p2q, idx1p2p, idx1q2q;

          if (d1q>d1p)
            idx1p1q =  (d1p*nD + d1q - (d1p*(d1p+1))/2)*nG + tuplr;
          else
            idx1p1q =  (d1q*nD + d1p - (d1q*(d1q+1))/2)*nG + tuplr;

          if (d2q>d2p)
            idx2p2q =  (d2p*nD + d2q - (d2p*(d2p+1))/2)*nG + tdlr;
          else
            idx2p2q =  (d2q*nD + d2p - (d2q*(d2q+1))/2)*nG + tdlr;

          if (d2p>d1p)
            idx1p2p =  (d1p*nD + d2p - (d1p*(d1p+1))/2)*nG + tlud;
          else
            idx1p2p =  (d2p*nD + d1p - (d2p*(d2p+1))/2)*nG + tlud;

          if (d1q>d2q)
            idx1q2q =  (d2q*nD + d1q - (d2q*(d2q+1))/2)*nG + trud;
          else
            idx1q2q =  (d1q*nD + d2q - (d1q*(d1q+1))/2)*nG + trud;


          if (empir == 0) { // this is model
            // up, left to right
            distV[idx1p1q] += exp(dsiArrayV[index][d1p*nD*nD*nD + d1q*nD*nD + d2p*nD + d2q]); // This value is already calculated at g 
            // left, up to down
            distV[idx1p2p] += exp(dsiArrayV[index][d1p*nD*nD*nD + d1q*nD*nD + d2p*nD + d2q]);
            // right, up to down
            distV[idx1q2q] += exp(dsiArrayV[index][d1p*nD*nD*nD + d1q*nD*nD + d2p*nD + d2q]);
            // down, left to right
            distV[idx2p2q] += exp(dsiArrayV[index][d1p*nD*nD*nD + d1q*nD*nD + d2p*nD + d2q]);

          } else {// this is empirical 
            distV[idx1p1q] ++;
            distV[idx2p2q] ++;
            distV[idx1p2p] ++;
            distV[idx1q2q] ++;
          }

        }
      }
  }
}



void copyImSet(std::vector <std::vector <CByteImage> > &dirDisp, std::vector <std::vector <CByteImage> > &dirWTA)
{

  for (unsigned int j=0; j < dirDisp.size(); j++) 
  {
    for(unsigned int i = 0; i < dirDisp[j].size(); i++) 
    {
      int ddpt1 = dirDisp[j].size();
      int ddpt2 = dirWTA[j].size();
      CShape sh1 = dirDisp[j][i].Shape();
      CShape sh2 = dirWTA[j][i].Shape();
      if (ddpt1 != ddpt2 || sh1.width != sh2.width || sh1.height != sh2.height || sh1.nBands != sh2.nBands) 
      {
        throw CError(" Depth and Shape dimensions of patient images should match ");
      } else 
      {
        // CopyPixels(dirDisp[j][i], dirWTA[j][i]);
        for (int m=0; m<sh1.height; m++) 
        {
          for (int n=0; n<sh1.width; n++) 
          {
            for (int p=0; p<sh1.nBands; p++) 
            {
              uchar *ptrp = &dirWTA[j][i].Pixel(n,m,p);
              *ptrp = dirDisp[j][i].Pixel(n,m,p);
            }
          }
        }

      }
    }
  }

}

/*

int checkNeighborPixel(int xp, int yp, int zp, int xq, int yq, int zq, int width, int height, int depth)
{
	if (xp==0 && xq<0)
		return -1;

	if (yp==0 && yq<0)
		return -1;

	if (zp==0 && zq<0)
		return -1;

	if (xp==width-1 && xq>width-1)
		return -1;

	if (yp==height-1 && yq>height-1)
		return -1;

	if (zp==depth-1 && zq>depth-1)
		return -1;

	return 0;

}


// using the 3D graph structure
// traverse thru all points.
// if the point is interactively labeled and all its neighbors are interactively labeled, don't update dcost
// if point is interactively labeled and its neighbor is not interactively labeled, update dcost of neighbor
// as dcost(d) + vcost(const, d) where const is the label of the interactively labeled point
// and d is the variable indexing possible labels for a point


void computeInterDataCost(DataCost *dcost, SmoothnessCost *scost, std::vector<CByteImage> interImage, std::vector<CByteImage> truedisp)
{

  int nD = globalNdisps;

  CShape sh = truedisp[0].Shape();
  int width = sh.width, height = sh.height;
  int depth = truedisp.size();

  int pIndex = 0;
  for (int z = 0; z < depth; z++)
  {
	  for (int y = 0; y < height; y++)
	  {
		  uchar *tdi = &truedisp[z].Pixel(0, y, 0);
		  uchar *inim = &interImage[z].Pixel(0, y, 0);

		  for (int x = 0; x < width; x++)
		  {

			  if (inim[x]==0)
			  {
				  // skip
			  } else
			  {
				  // list all neighbors of p(x,y,z)

				  // q(x+1, y, z)
				  if (checkNeighborPixel(x, y, z, x+1, y, z, width, height, depth)== 0)
				  { // fine
					  int qIndex = pIndex + 1;
					  if (inim[x+1] == 256)
					  {
						  // skip because it is labeled too
					  } else
					  {
						  for (int d=0; d<nD; d++)
						  {
							  dcost->updateDataCostRaw(pIndex*nD + d, dcost->getDataCostRaw(pIndex*nD + d) + scost->getSmoothnessCostRaw(pIndex, qIndex, tdi[x], d));
						  }
					  }
				  }

				  // q(x-1, y, z)
				  if (checkNeighborPixel(x, y, z, x-1, y, z, width, height, depth)== 0)
				  { // fine
					  int qIndex = pIndex - 1;
					  if (inim[x-1] == 256)
					  {
						  // skip because it is labeled too
					  } else
					  {
						  for (int d=0; d<nD; d++)
						  {
							  dcost->updateDataCostRaw(pIndex*nD + d, dcost->getDataCostRaw(pIndex*nD + d) + scost->getSmoothnessCostRaw(pIndex, qIndex, tdi[x], d));
						  }
					  }
				  }


				  // q(x, y+1, z)
				  if (checkNeighborPixel(x, y, z, x, y+1, z, width, height, depth)== 0)
				  { // fine

					  int qIndex = pIndex + width;
					  uchar *inimy = &interImage[z].Pixel(0, y+1, 0);
					  if (inimy[x] == 256)
					  {
						  // skip because it is labeled too
					  } else
					  {
						  for (int d=0; d<nD; d++)
						  {
							  dcost->updateDataCostRaw(pIndex*nD + d, dcost->getDataCostRaw(pIndex*nD + d) + scost->getSmoothnessCostRaw(pIndex, qIndex, tdi[x], d));
						  }
					  }
				  }

				  // q(x, y-1, z)
				  if (checkNeighborPixel(x, y, z, x, y-1, z, width, height, depth)== 0)
				  { // fine

					  int qIndex = pIndex - width;
					  uchar *inimy = &interImage[z].Pixel(0, y-1, 0);
					  if (inimy[x] == 256)
					  {
						  // skip because it is labeled too
					  } else
					  {
						  for (int d=0; d<nD; d++)
						  {
							  dcost->updateDataCostRaw(pIndex*nD + d, dcost->getDataCostRaw(pIndex*nD + d) + scost->getSmoothnessCostRaw(pIndex, qIndex, tdi[x], d));
						  }
					  }
				  }

				  // q(x, y, z+1)
				  if (checkNeighborPixel(x, y, z, x, y, z+1, width, height, depth)== 0)
				  { // fine

					  int qIndex = pIndex + width*height;
					  uchar *inimz = &interImage[z+1].Pixel(0, y, 0);
					  if (inimz[x] == 256)
					  {
						  // skip because it is labeled too
					  } else
					  {
						  for (int d=0; d<nD; d++)
						  {
							  dcost->updateDataCostRaw(pIndex*nD + d, dcost->getDataCostRaw(pIndex*nD + d) + scost->getSmoothnessCostRaw(pIndex, qIndex, tdi[x], d));
						  }
					  }
				  }

				  // q(x, y, z-1)
				  if (checkNeighborPixel(x, y, z, x, y, z-1, width, height, depth)== 0)
				  { // fine

					  int qIndex = pIndex - width*height;
					  uchar *inimz = &interImage[z-1].Pixel(0, y, 0);
					  if (inimz[x] == 256)
					  {
						  // skip because it is labeled too
					  } else
					  {
						  for (int d=0; d<nD; d++)
						  {
							  dcost->updateDataCostRaw(pIndex*nD + d, dcost->getDataCostRaw(pIndex*nD + d) + scost->getSmoothnessCostRaw(pIndex, qIndex, tdi[x], d));
						  }
					  }
				  }


			  }

			  pIndex++;
		  }
	  }
  }

}

*/

void checkReadImage(std::vector <std::vector <CByteImage> > interdirImage)
{

	int nump = 3;

	for (int i=2; i<3; i++)
	{

		  CShape sh = interdirImage[i][0].Shape();
		  int width = sh.width, height = sh.height;
		  int depth = interdirImage[i].size();

	      int n = 0;

	      for (int z = 0; z < depth; z++)
	      {
			  for (int y = 0; y < height; y++)
				{
				  for (int x = 0; x < width; x++)
					{
					  uchar *inim = &interdirImage[i][z].Pixel(x, y, 0);

					  std::cout << " " << (int)*inim << std::endl;

					  n++;
					}
				}
	      }
	  }


}















//this is the initialization using location
void initializeDataCostVectorPar(int genparam, std::vector <std::vector <CImage2> > im1,
                              std::vector <std::vector <CByteImage> > hogIm,
                              int nD,               // number of disparities
                              int numTrainingPats,
                              std::vector <std::vector <CByteImage> > &wta,      // winner-take-all disparities
                              fvec thetaU, fvec thetaA, fvec thetaH, fvec thetaL,
                              int featureCode, int featureCodeKlr,
                              std::vector <std::vector <std::vector <CByteImage> > > appdirImage,
                              Appearance** appclass, LocationMV** LocationMVclass, AppearanceMV** AppMVclass,
                              std::vector <std::vector <std::vector  <matrixB<double> > > > appdirMVProb,
                              intClass **intObj,
                              std::vector <std::vector <CByteImage> > gtIm,
                              std::vector <std::vector <std::vector <CImage2> > > locdirImage,
                              std::vector <std::vector <std::vector  <matrixB<double> > > > hogdirMVProb,
                              intensityNB* intNBclass,
                              matrix<DBL_TYPE> wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM, DBL_TYPE lambda, DBL_TYPE p1,
                              std::vector<int> startSliceNo, std::vector<int> endSliceNo,
                              int nbinsNB, int hogdimklr, int appdimklr, int locdimklr, int biasklr
                              )
{
  checkForFloatCosts();

  globalNdisps = nD;

  int intensity = (int)featureCode/10000;
  int hog = (int) (featureCode/1000)%10;
  int app = (int) (featureCode/100)%10;
  int loc = (int) (featureCode/10)%10;

  int klr = 0;

  int intensityKlr = (int)featureCodeKlr/10000;
  int hogKlr = (int) (featureCodeKlr/1000)%10;
  int appKlr = (int) (featureCodeKlr/100)%10;
  int locKlr = (int) (featureCodeKlr/10)%10;

  if (intensityKlr!=0 || appKlr!=0 || hogKlr!=0 || locKlr!=0)
    klr = 1;


  if (klr==1)
    {
      if (intensityKlr>0)
  		intensity=0;
      if (appKlr>0)
  		app=0;
      if (locKlr>0)
  		loc=0;
      if (hogKlr>0)
  		hog=0;
    }


  int dim = 0;
  IVM* ivm;


  if (intensityKlr==2 || hogKlr==2 || appKlr==2 || locKlr==2)
    {
	  klr = 1;
	  ivm = new IVM(&Xtrain, &wparam, RBF, lambda, nD, p1);
	  dim = Xtrain.size2();
    }



  int nU = thetaU.size();
  int nbins = nU/nD;

  // assuming that we have png data that goes upto 65535.
  // since i have the data in the range of 850 to 1250 which is shifted to 0 to 256
  double binner = 256.0/nbins;
  // so for any point. divide it by binner to get bin number (floor int value) it is in
  // ranging from bin0 to bin{nbins-1}
  // since >=256 will be special case, put it in the last bin if there.

  int nA = thetaA.size();

  int nT;

  nT = nA/nD; //this gives centers per class
  // this will work for app==2 and app==1 because in app==2, we have nD*the number of total clusters as the number of parameters
  // while in app==1 it is nD*number of clusters of that class (only that clusters are also = number of params)



  int nH = thetaH.size(); // number of HoG parameters
  int nHpC = nH/nD; // HoG vocabulary size

  int nL = thetaL.size(); // number of location parameters
  int nLpC = nL/nD;  // this variable will be useful only for cube setting (loc==1) or loc==9

  MRF::CostVal badcost = 1e20;

  for(int j=0; j<nU; j++)
    badcost += thetaU[j];

  if (nA > 0)
    for(int j=0; j<nA; j++)
      badcost += thetaA[j];

  if (nH > 0)
    for(int j=0; j<nH; j++)
      badcost += thetaH[j];

  if (nL > 0)
    for(int j=0; j<nL; j++)
      badcost += thetaL[j];




  if (genparam==1) {
    int numGens=0;
    if (intensity==3 || intensity==6)
      numGens++;
    if (loc==3 || loc==8 || loc==9)
      numGens++;
    if (app==3)
      numGens++;
    if (hog==3)
      numGens++;


    //WORK-THIS
    // assuming worst cost after taking log is -log(10^-99)

    badcost += numGens*(99*3);  //that 3 accounts for 2.7 value of e (natural log)

  }




  for (unsigned int j = 0; j < im1.size(); ++j) {

    int depth = im1[j].size();
    CShape sh = im1[j][0].Shape();
    int width = sh.width, height = sh.height, nB = sh.nBands;
    int tempglobalNpixels = depth * width * height;
    int nColors = __min(3, nB);

    if (genparam==0)
      dsiArrayV.push_back(new MRF::CostVal[tempglobalNpixels * globalNdisps]);
    else if (genparam==1)
      dsiArrayVGen.push_back(new MRF::CostVal[tempglobalNpixels * globalNdisps]);

    int zslice = 0;
    zslice = startSliceNo[j];


    int dsiIndex = 0;
    for(unsigned int i = 0; i < im1[j].size(); ++i){

      DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);
      sh.nBands = 1;
      wta[j][i].ReAllocate(sh);

      matrix<DBL_TYPE> hogdescs;

      if (hogKlr==2)
        {
          // read file
          char hogfile[100];

          sprintf(hogfile, "/p/imageseg/medicalSeg/parameters/hog/hog2dcsv_90n/hog2ddata_%d.csv", i+zslice );

          // read this file in matrix
          if(!load_one_data(hogfile, hogdescs)){
            std::cout << "csv read error: Bad data file filename : "<< hogfile << endl;
            exit(1);
          }

        }

      // need to load appfile if appklr==2
      matrix<DBL_TYPE> appdescs;

      if (appKlr==2)
        {
          // read file
          char appfile[100];

          sprintf(appfile, "/p/imageseg/medicalSeg/parameters/app/appcsv_81n/app2ddatan_%d.csv", i+zslice );

          // read this file in matrix
          if(!load_one_data(appfile, appdescs)){
            std::cout << "csv read error: Bad data file filename : "<< appfile << endl;
            exit(1);
          }

        }
/*
#pragma omp parallel for
      for (int y = 0; y < height; y++) {
        uchar *WTArow = &wta[j][i].Pixel(0, y, 0);
        for (int x = 0; x < width; x++) {
          unsigned short pix1 = im1[j][i].Pixel(x, y, 0);

          uchar *gtpix = &gtIm[j][i].Pixel(x, y, 0);

          uchar *hogval;
          if (nH > 0)
            hogval = &hogIm[j][i].Pixel(x, y, 0);

          ublas::vector<DBL_TYPE> f;
          MRF::CostVal bestval = badcost;
          int* bestd = new int [nD];
          int numbest = 0;

          int thetaid;

          if (klr==1)
            {
              ublas::vector<DBL_TYPE> Xtest(dim);

              int idxx = 0;

              if (intensityKlr==2) // for now only intensity case
                {
                  getCostIntensityKlr(Xtest, (double)pix1, minmaxM, idxx);
                }

              if (locKlr==2)
                {
                  // 25 values. // this should be a parameter  // location features
                  // x y z <prob of belonging to each cluster>*16 <prob of belong to each class>*6
                  getCostLocationKlr(Xtest, i, zslice, x, y, nD, minmaxM, LocationMVclass, idxx, skipXYZlocklr);

                }

              if (appKlr==2)
                {
                  getCostAppKlr(Xtest, x, y, width, minmaxM, idxx, appdimklr, &appdescs);
                }


              if (hogKlr==2)
                {
                  getCostHoGKlr(Xtest, x, y, width, minmaxM, idxx, hogdimklr, &hogdescs);
                }


              f = ivm->classify_klrexp(&Xtest);

            }




          for (int d = 0; d < nD; d++) {

            MRF::CostVal dsiValue = 0;

            if (klr==1) {
              dsiValue += -f(d);
              dsiValue += biasklr;
            }


            if (nU>0 && intensity==1)  // intensity as histogram bins
              {
                dsiValue += getCostIntensityDiscBins(pix1, nbins, d, thetaU, thetaid);

              } else if (intensity==3 && genparam==1) // generative setting with Gaussians
              {
                dsiValue += getCostIntensityGenGaussian(intObj[d], pix1);

              } else if (intensity == 6 && genparam==1) // using NBayes and bins
              {
                dsiValue += getCostIntensityGenNB(pix1, nbins, d, intNBclass);
              }


            if (nA>0 && app==1) // we have appearance features
              {
                // different features for different classes

                dsiValue += getCostAppDiscPatch((int)appdirImage[j][i][d].Pixel(x, y, 0), appclass[d]->getPatchSize(), x, y, width, height, thetaA, d, nT, thetaid);

              } else if (nA>0 && app==2) // we have appearance features
              {
                // common features for all classes

                dsiValue += getCostAppDiscPatch((int)appdirImage[j][i][0].Pixel(x, y, 0), appclass[0]->getPatchSize(), x, y, width, height, thetaA, d, nT, thetaid);

              }
            else if (app==3 && genparam==1) // generative setting
              {

                dsiValue += getCostAppGenPatch((appdirMVProb[j][i][d])(y,x), AppMVclass[d]->getDimension(), x, y, width, height);

              }



            if (nH > 0 && hog==1) {

              dsiValue += getCostHoGDiscBins((int) *hogval, d, nHpC, thetaH, thetaid);

            } else if (hog==3 && genparam==1) // generative setting
              {

                int startPatch = 8; // index starting from 0
                // this variable indicates where the Hog descriptors are defined
                // since i have generated them using the entire patient body, don't need to
                // worry about the z axis and can use all slices here.
                //WORK-THIS need to change above stuff for flawless working

                dsiValue += getCostHoGGenBins((hogdirMVProb[j][i][d])(y,x), x, y, width, height, startPatch);

              }

            if (nL>0 && loc==1) { // CRF cube code

              // this function is completely out-dated
              // and for older series of data, will need to modify significantly for use

              dsiValue += getCostLocationDiscCube((unsigned short) locdirImage[j][d][i].Pixel(x,y,0), x, y, i, loccubex, loccubey, loccubez, numTrainingPats, nLpC, thetaL, width, height, depth, d, thetaid);

            } else if (nL>0 && (loc==6 || loc==8) )  //CRF hard or soft assignment
              {

                dsiValue += getCostLocationDiscGaussian(LocationMVclass, i, zslice, x, y, thetaL, d, loc, genparam, thetaid);

              } else if (nL>0 && (loc==7 || loc==9) )  //CRF soft or hard assignment redefined
              {

                dsiValue += getCostLocationDiscGaussian2(LocationMVclass, i, zslice, x, y, thetaL, d, loc, genparam, nD, thetaid);

              } else if (loc==3 && genparam==1) // generative setting
              {

                dsiValue += getCostLocationGenGaussian(LocationMVclass[d], i, zslice, x, y);

              }

            // The cost of pixel p and label l is stored at dsiArray[p*nLabels+l]
            int dsiPar = i*height*width*nD + y*width*nD + x*nD + d;
            if (genparam==0)
              (dsiArrayV[j])[dsiPar] = (float)dsiValue;
            else if (genparam==1)
              (dsiArrayVGen[j])[dsiPar] = (float)dsiValue;

            if (dsiValue < bestval) {
              bestval = dsiValue;
              bestd[0] = d;
              numbest = 1;
            } else if (dsiValue == bestval) {
              bestd[numbest] = d;
              numbest++;
            }

            if (pix1 >= 256) {
              bestd[0] = 0; // assume 0 is background
              numbest = 1;
            }
          }
          if (numbest == 0)
          	numbest = nD;
          int curr = rand() % numbest;
          WTArow[x] = bestd[curr];

          delete [] bestd;
        }
      }
*/

    }
  }



  if (klr==1)
    delete ivm;


}



DataCost *computeDataCost(int generative, int index)
{
  return new DataCost ( computeDataCostArray(generative, index));
}


MRF::CostVal* computeDataCostArray(int generative, int index)
{
	if (generative==0)
		return dsiArrayV[index];
	else if (generative==1)
		return dsiArrayVGen[index];
}

/*

// create data cost using parameters thetaU
DataCost *computeDataCostPar(int generative, std::vector<CImage2> set1Im, std::vector <CByteImage> hogIm, int numTrainingPats, std::vector<float> thetaU, std::vector<float> thetaA, std::vector<float> thetaH, std::vector<float> thetaL, int featureCode, int featureCodeKlr, std::vector <std::vector <CByteImage> >  appImage, Appearance** appclass, int index, LocationMV** LocationMVclass, AppearanceMV** AppMVclass, std::vector <std::vector  <matrixB<double> > > appdirMVProb, intClass** intObj, std::vector <std::vector <CImage2> > locdirImage, std::vector <std::vector  <matrixB<double> > > hogdirMVProb, intensityNB* intNBclass, matrix<DBL_TYPE> wparam,
                          matrix<DBL_TYPE> Xtrain,
                          matrix<DBL_TYPE> minmaxM, DBL_TYPE lambda, DBL_TYPE p1,
                          std::vector<int> startSliceNo, std::vector<int> endSliceNo,
                          int nbinsNB, int hogdimklr, int appdimklr, int locdimklr, int biasklr)
{
  return new DataCost ( computeDataCostArrayPar (generative, set1Im, hogIm, numTrainingPats, thetaU, thetaA, thetaH, thetaL, featureCode, featureCodeKlr, appImage, appclass, index, LocationMVclass, AppMVclass, appdirMVProb, intObj, locdirImage, hogdirMVProb, intNBclass, wparam, Xtrain, minmaxM, lambda, p1, startSliceNo, endSliceNo, nbinsNB, hogdimklr, appdimklr, locdimklr, biasklr));
}

MRF::CostVal* computeDataCostArrayPar(int generative, std::vector<CImage2> im1, std::vector <CByteImage> hogIm, int numTrainingPats, std::vector<float> thetaU, std::vector<float> thetaA, std::vector<float> thetaH, std::vector<float> thetaL, int featureCode, int featureCodeKlr, std::vector <std::vector <CByteImage> >  appImage, Appearance** appclass, int index, LocationMV** LocationMVclass, AppearanceMV** AppMVclass, std::vector <std::vector  <matrixB<double> > >  appdirMVProb, intClass** intObj, std::vector <std::vector <CImage2> > locdirImage, std::vector <std::vector  <matrixB<double> > > hogdirMVProb, intensityNB* intNBclass, 		matrix<DBL_TYPE> wparam,
                                   matrix<DBL_TYPE> Xtrain,
                                   matrix<DBL_TYPE> minmaxM, DBL_TYPE lambda, DBL_TYPE p1,
                                   std::vector<int> startSliceNo, std::vector<int> endSliceNo,
                                   int nbinsNB, int hogdimklr, int appdimklr, int locdimklr, int biasklr)
{

  int intensity = (int)featureCode/10000;
  int hog = (int) (featureCode/1000)%10;
  int app = (int) (featureCode/100)%10;
  int loc = (int) (featureCode/10)%10;

  int nD = globalNdisps;


  int klr = 0;

  int intensityKlr = (int)featureCodeKlr/10000;
  int hogKlr = (int) (featureCodeKlr/1000)%10;
  int appKlr = (int) (featureCodeKlr/100)%10;
  int locKlr = (int) (featureCodeKlr/10)%10;

  if (intensityKlr!=0 || appKlr!=0 || hogKlr!=0 || locKlr!=0)
    klr = 1;


  if (klr==1)
    {
      if (intensityKlr>0)
        intensity=0;
      if (appKlr>0)
        app=0;
      if (locKlr>0)
        loc=0;
      if (hogKlr>0)
        hog=0;
    }

  int dim = 0;
  IVM* ivm;



  if (intensityKlr==2 || hogKlr==2 || appKlr==2 || locKlr==2)
    {
      klr = 1;

      ivm = new IVM(&Xtrain, &wparam, RBF, lambda, nD, p1);

      dim = Xtrain.size2();
    }




  int zslice = 0;

  if (generative==0) {
    if (dsiArrayV.empty() == true)
      throw CError("call initializeDataCost first");
  } else if (generative==1) {
    if (dsiArrayVGen.empty() == true)
      throw CError("call initializeDataCostLoc first");
  }

  // currently used for slice location in a training subpart of volume

  zslice = startSliceNo[index];

  int nU = (int)thetaU.size();
  int nbins = nU/globalNdisps;

  //	double binner = 256.0/nbins;
  int depth = im1.size();
  CShape sh = im1[0].Shape();
  int width = sh.width, height = sh.height;
  //int nB = sh.nBands;

  int nA = thetaA.size();
  int nT = nA/globalNdisps; //this gives centers per class

  int nH = thetaH.size(); // number of HoG parameters
  int nHpC = nH/globalNdisps; // HoG vocabulary size

  int nL = thetaL.size();
  int nLpC = nL/globalNdisps;

  int dsiIndex = 0;



  for (int z=0; z < depth; z++) {


    matrix<DBL_TYPE> hogdescs;

    if (hogKlr==2)
      {
        // read file
        char hogfile[100];
        sprintf(hogfile, "/p/imageseg/medicalSeg/parameters/hog/hog2dcsv_90n/hog2ddata_%d.csv", z+zslice );  //change this

        // read this file in matrix

        if(!load_one_data(hogfile, hogdescs)){
          std::cout << "csv read error: Bad data file filename : "<< hogfile << endl;
          exit(1);
        }

      }


    matrix<DBL_TYPE> appdescs;

    if (appKlr==2)
      {
        // read file
        char appfile[100];

        sprintf(appfile, "/p/imageseg/medicalSeg/parameters/app/appcsv_81n/app2ddatan_%d.csv", z+zslice );

        // read this file in matrix

        if(!load_one_data(appfile, appdescs)){
          std::cout << "csv read error: Bad data file filename : "<< appfile << endl;
          exit(1);
        }
      }



#pragma omp parallel for
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {

        unsigned short pix1 = im1[z].Pixel(x, y, 0);

        ublas::vector<DBL_TYPE> f;
        int thetaid;

        uchar *hogval;
        if (nH > 0)
          hogval = &hogIm[z].Pixel(x, y, 0);


        if (klr==1)
          {
            ublas::vector<DBL_TYPE> Xtest(dim);

            int idxx = 0;

            if (intensityKlr==2) // for now only intensity case
              {
                getCostIntensityKlr(Xtest, (double)pix1, minmaxM, idxx);
              }

            if (locKlr==2)
              {
                getCostLocationKlr(Xtest, z, zslice, x, y, nD, minmaxM, LocationMVclass, idxx, skipXYZlocklr);
              }

            if (appKlr==2)
              {
                getCostAppKlr(Xtest, x, y, width, minmaxM, idxx, appdimklr, &appdescs);
              }


            if (hogKlr==2)
              {
                getCostHoGKlr(Xtest, x, y, width, minmaxM, idxx, hogdimklr, &hogdescs);
              }


            f = ivm->classify_klrexp(&Xtest);
          }


        for (int d = 0; d < globalNdisps; d++) {

          MRF::CostVal dsiValue = 0;

          if (klr==1) {
            dsiValue += -f(d);
            dsiValue += biasklr;
          }


          if (nU>0 && intensity==1)
            {
              dsiValue += getCostIntensityDiscBins(pix1, nbins, d, thetaU, thetaid);

            } else if (intensity==3 && generative==1) // generative setting
            {
              dsiValue += getCostIntensityGenGaussian(intObj[d], pix1);

            } else if (intensity == 6 && generative==1) // using NBayes and bins
            {
              dsiValue += getCostIntensityGenNB(pix1, nbins, d, intNBclass);
            }



          if (nA>0 && app==1) // we have appearance features
            {
              // different features for different classes
              dsiValue += getCostAppDiscPatch((int)appImage[z][d].Pixel(x, y, 0), appclass[d]->getPatchSize(), x, y, width, height, thetaA, d, nT, thetaid);

            } else if (nA>0 && app==2) // we have appearance features
            {
              // common features for all classes
              dsiValue += getCostAppDiscPatch((int)appImage[z][0].Pixel(x, y, 0), appclass[0]->getPatchSize(), x, y, width, height, thetaA, d, nT, thetaid);
            }
          else if (app==3  && generative==1) // generative setting
            {
              dsiValue += getCostAppGenPatch((appdirMVProb[z][d])(y,x), AppMVclass[d]->getDimension(), x, y, width, height);
            }



          if (nH > 0 && hog==1)
            {
              dsiValue += getCostHoGDiscBins((int) *hogval, d, nHpC, thetaH, thetaid);

            } else if (hog==3  && generative==1) // generative setting
            {

              int startPatch = 8; // index starting from 0
              // this variable indicates where the Hog descriptors are defined
              // since i have generated them using the entire patient body, don't need to
              // worry about the z axis and can use all slices here.
              //WORK-THIS need to change above stuff for flawless working
              dsiValue += getCostHoGGenBins((hogdirMVProb[z][d])(y,x), x, y, width, height, startPatch);

            }


          if (nL>0 && loc==1) { //WORK-THIS check boundary conditions

            // this function is completely outdated
            // and for older series of data, will need to modify significantly for use
            dsiValue += getCostLocationDiscCube((unsigned short) locdirImage[d][z].Pixel(x,y,0), x, y, z, loccubex, loccubey, loccubez, numTrainingPats, nLpC, thetaL, width, height, depth, d, thetaid);

          } else if (nL>0 && (loc==6 || loc==8)) // CRF hard or soft
            {

              dsiValue += getCostLocationDiscGaussian(LocationMVclass, z, zslice, x, y, thetaL, d, loc, generative, thetaid);

            } else if (nL>0 && (loc==7 || loc==9) )  //CRF soft or hard assignment redefined
            {

              dsiValue += getCostLocationDiscGaussian2(LocationMVclass, z, zslice, x, y, thetaL, d, loc, generative, nD, thetaid);

            }
          else if (loc==3 && generative==1) // generative setting
            {

              dsiValue += getCostLocationGenGaussian(LocationMVclass[d], z, zslice, x, y);

            }

          int dsiPar = z*height*width*nD + y*width*nD + x*nD + d;
          if (generative==0) {
            (dsiArrayV[index])[dsiPar] = (float)dsiValue;
          } else if (generative==1) {
            (dsiArrayVGen[index])[dsiPar] = (float)dsiValue;
          }

        }
      }
    }
  }
  if (generative==0)
    return dsiArrayV[index];
  else if (generative==1)
    return dsiArrayVGen[index];


  if (klr==1)
    delete ivm;

}
*/




