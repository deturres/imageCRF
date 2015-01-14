// expectations.cpp
// modified by Chetan Bhole for image segmentation
//
// $Id: expectations.cpp,v 1.11 2007/12/11 13:28:48 weinman Exp $
//
// code to compute the empirical and model expections for the gradient update
// Daniel Scharstein

#include "crfmodel.h"
#include "expectations.h"
#include <assert.h>
#include <math.h>
#include <iostream>
#include <iomanip>
using std::cout;


#define ABS(x) ((x)>0? (x) : (-(x)))

// these parameters are also defined in crfmodel.cpp
// make sure to change them there too
extern int loccubex, loccubey, loccubez;
extern int skipXYZlocklr;


void computeMaskedImage(CByteImage disp, CByteImage truedisp, int ignoreVal)
{

  if (ignoreVal<0)
    return;

  CShape sh = disp.Shape();
  int width = sh.width, height = sh.height;

  for (int y = 0; y < height; y++) {
	uchar *di = &disp.Pixel(0, y, 0);
	uchar *ti = &truedisp.Pixel(0, y, 0);

	for (int x = 0; x < width; x++) {
      if (ti[x]==ignoreVal) { // only evaluate non-ignored pixels
		di[x]=ignoreVal;
      }
	}
  }
}

float computeMaskedLikelihood(CByteImage disp, SyncMeanField *mrf)
{

  CShape sh = disp.Shape();
  int width = sh.width, height = sh.height;

  register float loglik = 0;

  for (int y = 0; y < height; y++) {

	uchar *ti = &disp.Pixel(0, y, 0);

	for (int x = 0; x < width; x++) {

      if (ti[x]!=0)
        loglik += mrf->getBelief(mrf->pixelIndex(y,x),ti[x]);
    }
  }

  return loglik;
}

// Occlusion model
// return data cost at occluded pixels
void computeEmpirDistU(CByteImage disp, fvec &distU, uchar occlVal, int ignoreVal)
{
  CShape sh = disp.Shape();

  int width = sh.width, height = sh.height;

  int nU = (int)distU.size();

  if (nU == 0)
	return; // nothing to do


  float *costs;
  int nPixels, nStates;

  getCostArray(costs, nPixels, nStates);

  assert(nPixels == width * height);

  float totCost = 0;
  float totOccl = 0;

  for (int y = 0; y < height; y++) {

	uchar *di = &disp.Pixel(0, y, 0);

	for (int x = 0; x < width; x++) {

      if ((int)(di[x])!=ignoreVal) {

        if (di[x]!=occlVal){
          totCost += costs[di[x]];
        }
        else{
          // cost of occluded (delta function)
          ++totOccl;

        }
      }

      // move costs pointer to next pixel
      costs += nStates;
	}
  }

  distU[0] = totOccl;

  if(nU == 2)
	distU[1] = totCost;

}


// No occlusion model
// return data cost at non-occluded pixels

// grayscale***
void computeEmpirDistU(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distU, int ignoreVal, std::vector<CImage2> im1, int nD, std::vector<float> thetaU)
{
  int depth = validdisp.size();
  CShape sh = validdisp[0].Shape();

  int width = sh.width, height = sh.height;

  int nU = (int)distU.size();

  for (int kk=0; kk<nU; kk++) {
	  distU[kk] = 0;
  }

  if (nU == 0)
	return; // nothing to do

  int nbins = nU/nD;
  double binner = 256.0/nbins;

  for (int z = 0; z < depth; z++) {

	  for (int y = 0; y < height; y++) {

		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *tdi = &validdisp[z].Pixel(0, y, 0);

		for (int x = 0; x < width; x++) {

		  unsigned short pix1 = im1[z].Pixel(x, y, 0);

		  int thetaid;
		  getCostIntensityDiscBins(pix1, nbins, di[x], thetaU, thetaid);

	// for us the costs are the theta param, but in differentiation, they vanish,
	// so its just whether the feature is enabled or not.
	// move costs pointer to next pixel
		  distU[thetaid] ++;
		}
	  }
  }
}

//RGB***
void computeEmpirDistU(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distU, int ignoreVal, std::vector<CByteImage> im1, int nD, std::vector<float> thetaU, std::vector<CByteImage> interImage, int interactive)
{
  int depth = validdisp.size();
  CShape sh = validdisp[0].Shape();

  int width = sh.width, height = sh.height;

  int nU = (int)distU.size();

  for (int kk=0; kk<nU; kk++) {
	  distU[kk] = 0;
  }

  if (nU == 0)
	return; // nothing to do

  int nbins = nU/nD;
  double binner = 256.0/nbins;

  for (int z = 0; z < depth; z++) {

	  for (int y = 0; y < height; y++) {

		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *tdi = &validdisp[z].Pixel(0, y, 0);

		for (int x = 0; x < width; x++) {

		  unsigned char *pix1 = &im1[z].Pixel(x, y, 0);

			if (interactive==1)
			{
				uchar *inimp = &interImage[z].Pixel(x, y, 0);
				if (*inimp!=0)
					continue;
			}


		  int thetaidrgb[3]={0,0,0};
		  int thetaidYuv[2]={0,0};
		  // getCostIntensityDiscBins(pix1, nbins, di[x], thetaU, thetaidrgb);
		  getCostIntensityDiscBins(pix1, nbins, di[x], thetaU, thetaidYuv);

	// for us the costs are the theta param, but in differentiation, they vanish,
	// so its just whether the feature is enabled or not.
	// move costs pointer to next pixel
		  // distU[thetaidrgb[0]] ++;
		  // distU[thetaidrgb[1]] ++;
		  // distU[thetaidrgb[2]] ++;

		  distU[thetaidYuv[0]] ++;
		  distU[thetaidYuv[1]] ++;

		}
	  }
  }
}




void computeEmpirDistH(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distH, int ignoreVal, std::vector<CByteImage> hogIm, int nD, std::vector<float> thetaH)
{
  int depth = validdisp.size();
  CShape sh = validdisp[0].Shape();

  int width = sh.width, height = sh.height;

  int nH = (int)distH.size();

  if (nH == 0)
	return; // nothing to do

  for (int kk=0; kk<nH; kk++) {
	  distH[kk] = 0;
  }

  int nHpC = nH/nD;

  for (int z = 0; z < depth; z++) {

	  for (int y = 0; y < height; y++) {

		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *tdi = &validdisp[z].Pixel(0, y, 0);

		for (int x = 0; x < width; x++) {

          uchar *hogval = &hogIm[z].Pixel(x, y, 0);

          int thetaid;
          getCostHoGDiscBins((int) *hogval, di[x], nHpC, thetaH, thetaid);

          /*
			int hogbin = (int) *hogval;
			int thetaid = di[x]*nHpC + hogbin;
          */
          // for us the costs are the theta param, but in differentiation, they vanish,
          // so its just whether the feature is enabled or not.
          // move costs pointer to next pixel
          distH[thetaid] ++;

		}
	  }
  }
}


void computeEmpirDistM(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distM, int ignoreVal, std::vector<CByteImage> motIm, int nD, std::vector<float> thetaM)
{
  int depth = validdisp.size();
  CShape sh = validdisp[0].Shape();

  int width = sh.width, height = sh.height;

  int nM = (int)distM.size();

  if (nM == 0)
	return; // nothing to do

  for (int kk=0; kk<nM; kk++) {
	  distM[kk] = 0;
  }

  int nMpC = nM/nD;

  for (int z = 0; z < depth; z++) {

	  for (int y = 0; y < height; y++) {

		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *tdi = &validdisp[z].Pixel(0, y, 0);

		for (int x = 0; x < width; x++) {

          uchar *motval = &motIm[z].Pixel(x, y, 0);

          int thetaid;
          getCostmotDiscBins((int) *motval, di[x], nMpC, thetaM, thetaid);

          /*
			int hogbin = (int) *hogval;
			int thetaid = di[x]*nHpC + hogbin;
          */
          // for us the costs are the theta param, but in differentiation, they vanish,
          // so its just whether the feature is enabled or not.
          // move costs pointer to next pixel
          distM[thetaid] ++;

		}
	  }
  }
}







void computeEmpirDistL(int loc, int index, std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distL, int ignoreVal, int nD, std::vector <std::vector <CImage2> > locdirImage, int numTrainingPats, LocationMV** LocationMVclass,  std::vector<float> thetaL, int genparam, std::vector<int> startSliceNo)
{
  int depth = validdisp.size();
  CShape sh = validdisp[0].Shape();

  int width = sh.width, height = sh.height;

  int nL = (int)distL.size();

  for (int kk=0; kk<nL; kk++) {
	  distL[kk] = 0;
  }

  if (nL == 0)
	return; //nothing to do

  int nLpC = nL/nD;

  int zslice=0;
  zslice = startSliceNo[index];

    /*  if (index==0)
	zslice = 195;
  else if (index==1)
	zslice = 225;
    */

  for (int z = 0; z < depth; z++) {

	  for (int y = 0; y < height; y++) {

		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *tdi = &validdisp[z].Pixel(0, y, 0);

		for (int x = 0; x < width; x++) {

          if (loc==1) {

            unsigned short locpix = locdirImage[di[x]][z].Pixel(x,y,0);

            int thetaid;
            getCostLocationDiscCube(locpix, x, y, z, loccubex, loccubey, loccubez, numTrainingPats, nLpC, thetaL, width, height, depth, di[x], thetaid);
			  
            double locpixd = (double)locpix/ (loccubex * loccubey * loccubez * numTrainingPats);
            // for us the costs are the theta param, but in differentiation, they vanish,
            // so its just the value of the feature.
            // move costs pointer to next pixel
            distL[thetaid] += (1-locpixd);

          } else if (loc==6 || loc==8) { // CRF hard assignment and soft

            int thetaid;
            getCostLocationDiscGaussian(LocationMVclass, z, zslice, x, y, thetaL, di[x], loc, genparam, thetaid);

            distL[thetaid]++;

          } else if (loc==9 || loc==7) {  //CRF soft or hard assignment redefined
            
              int thetaid;
              getCostLocationDiscGaussian2(LocationMVclass, z, zslice, x, y, thetaL, di[x], loc, genparam, nD, thetaid);

              distL[thetaid]++;
          }

		}
	  }
  }
}




void computeEmpirDistW(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distW,  matrix<DBL_TYPE> &wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM, int ignoreVal, std::vector<CImage2> im1, int nD, DBL_TYPE lambda, DBL_TYPE p1, int featureCodeKlr, int index, LocationMV** LocationMVclass, std::vector<int> startSliceNo, std::vector<int> endSliceNo, int appdimklr, int hogdimklr)
{
  int depth = validdisp.size();
  CShape sh = validdisp[0].Shape();

  int width = sh.width, height = sh.height;


  int klr = 1; // if this function is called this is true

  int intensityKlr = (int)featureCodeKlr/10000;
  int hogKlr = (int) (featureCodeKlr/1000)%10;
  int appKlr = (int) (featureCodeKlr/100)%10;
  int locKlr = (int) (featureCodeKlr/10)%10;

  int dimW = wparam.size2();

  int dim = Xtrain.size2();
  int Ntr = wparam.size1();


  int nW = (int)distW.size();

  for (int kk=0; kk<nW; kk++) {
	  distW[kk] = 0;
  }

  int zslice = 0;
  zslice = startSliceNo[index];

  if (nW == 0)
	return; // nothing to do


  //std::cout<<" Xtrain  "<<Xtrain<<std::endl;

  IVM* ivm = new IVM(&Xtrain, &wparam, RBF, lambda, nD, p1);;


  for (int z = 0; z < depth; z++) {


	  matrix<DBL_TYPE> hogdescs;

	  if (hogKlr==2)
	  {
		  // read file
		  char hogfile[100];

		  sprintf(hogfile, "/p/imageseg/medicalSeg/parameters/hog/hog2dcsv_90n/hog2ddata_%d.csv", z+zslice );

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

		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *tdi = &validdisp[z].Pixel(0, y, 0);

		for (int x = 0; x < width; x++) {

			unsigned short pix1 = im1[z].Pixel(x, y, 0);

			  ublas::vector<DBL_TYPE> Xtest(dim);

			  int idxx = 0;

			  if (intensityKlr==2) // for now only intensity case
			  {
                  getCostIntensityKlr(Xtest, (double)pix1, minmaxM, idxx);
			  }

			  if (locKlr==2)
			  {
				  // 25 values. // this should be a parameter
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

			  for (int j=0; j<Ntr; j++)
			  {
				int thetaid = j*dimW + di[x];

/*				std::cout<<" Xtest out  "<<Xtest<<std::endl;
				if (di[x]>0)
					std::cout<<" Xtest "<<Xtest<<std::endl;
*/
				double Kval = ivm->classify_klrexp(&Xtest, j);
				distW[thetaid] += Kval;
			  }
			}
		  }
	  }

  delete ivm;
}







void computeEmpirDistWPar(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distW,  matrix<DBL_TYPE> &wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM, int ignoreVal, std::vector<CImage2> im1, int nD, DBL_TYPE lambda, DBL_TYPE p1, int featureCodeKlr, int index, LocationMV** LocationMVclass, std::vector<int> startSliceNo, std::vector<int> endSliceNo, int appdimklr, int hogdimklr)
{
  int depth = validdisp.size();
  CShape sh = validdisp[0].Shape();

  int width = sh.width, height = sh.height;


  int klr = 1; // if this function is called this is true

  int intensityKlr = (int)featureCodeKlr/10000;
  int hogKlr = (int) (featureCodeKlr/1000)%10;
  int appKlr = (int) (featureCodeKlr/100)%10;
  int locKlr = (int) (featureCodeKlr/10)%10;

  int dimW = wparam.size2();

  int dim = Xtrain.size2();
  int Ntr = wparam.size1();


  int nW = (int)distW.size();

  for (int kk=0; kk<nW; kk++) {
	  distW[kk] = 0;
  }

  int zslice = 0;
  zslice = startSliceNo[index];

  if (nW == 0)
	return; // nothing to do

  IVM* ivm = new IVM(&Xtrain, &wparam, RBF, lambda, nD, p1);;


  for (int z = 0; z < depth; z++) {


	  matrix<DBL_TYPE> hogdescs;

	  if (hogKlr==2)
	  {
		  // read file
		  char hogfile[100];

		  sprintf(hogfile, "/p/imageseg/medicalSeg/parameters/hog/hog2dcsv_90n/hog2ddata_%d.csv", z+zslice );

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

		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *tdi = &validdisp[z].Pixel(0, y, 0);

		for (int x = 0; x < width; x++) {

			unsigned short pix1 = im1[z].Pixel(x, y, 0);

			  ublas::vector<DBL_TYPE> Xtest(dim);

			  int idxx = 0;

			  if (intensityKlr==2) // for now only intensity case
			  {
                  getCostIntensityKlr(Xtest, (double)pix1, minmaxM, idxx);
			  }

			  if (locKlr==2)
			  {
				  // 25 values. // this should be a parameter
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

			  for (int j=0; j<Ntr; j++)
			  {
				int thetaid = j*dimW + di[x];

				double Kval = ivm->classify_klrexp(&Xtest, j);

				#pragma omp critical
				distW[thetaid] += Kval;
			  }
			}
		  }
	  }

  delete ivm;
}



void computeEmpirDistWPar(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distW,  matrix<DBL_TYPE> &wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM, int ignoreVal, std::vector<CImage2> im1, int nD, DBL_TYPE lambda, DBL_TYPE p1, int featureCodeKlr, int index, LocationMV** LocationMVclass, std::vector<int> startSliceNo, std::vector<int> endSliceNo, int appdimklr, int hogdimklr, MRFEnergy* mrf)
{
  int depth = validdisp.size();
  CShape sh = validdisp[0].Shape();

  int width = sh.width, height = sh.height;


  int klr = 1; // if this function is called this is true

  int intensityKlr = (int)featureCodeKlr/10000;
  int hogKlr = (int) (featureCodeKlr/1000)%10;
  int appKlr = (int) (featureCodeKlr/100)%10;
  int locKlr = (int) (featureCodeKlr/10)%10;

  int dimW = wparam.size2();

  int dim = Xtrain.size2();
  int Ntr = wparam.size1();


  int nW = (int)distW.size();

  for (int kk=0; kk<nW; kk++) {
	  distW[kk] = 0;
  }

  int zslice = 0;
  zslice = startSliceNo[index];

  if (nW == 0)
	return; // nothing to do

  IVM* ivm = new IVM(&Xtrain, &wparam, RBF, lambda, nD, p1);;


  for (int z = 0; z < depth; z++) {


	  matrix<DBL_TYPE> hogdescs;

	  if (hogKlr==2)
	  {
		  // read file
		  char hogfile[100];

		  sprintf(hogfile, "/p/imageseg/medicalSeg/parameters/hog/hog2dcsv_90n/hog2ddata_%d.csv", z+zslice );

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

		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *tdi = &validdisp[z].Pixel(0, y, 0);

		for (int x = 0; x < width; x++) {

			unsigned short pix1 = im1[z].Pixel(x, y, 0);

			float *belief = new float[nD];

			int pixel = z*width*height + y*width + x;
			mrf->getBelief(pixel,belief);


			  ublas::vector<DBL_TYPE> Xtest(dim);

			  int idxx = 0;

			  if (intensityKlr==2) // for now only intensity case
			  {
                  getCostIntensityKlr(Xtest, (double)pix1, minmaxM, idxx);
			  }

			  if (locKlr==2)
			  {
				  // 25 values. // this should be a parameter
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


			  for (int j=0; j<Ntr; j++)
			  {
				  double Kval = ivm->classify_klrexp(&Xtest, j);

				  for (int d=0; d<nD; d++)
				  {
					int thetaid = j*dimW + d;

					#pragma omp critical
					distW[thetaid] += exp(belief[d])*Kval;
				  }
			  }
			}
		  }
	  }

  delete ivm;
}


void computeEmpirDistW(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distW,  matrix<DBL_TYPE> &wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM, int ignoreVal, std::vector<CImage2> im1, int nD, DBL_TYPE lambda, DBL_TYPE p1, int featureCodeKlr, int index, LocationMV** LocationMVclass, std::vector<int> startSliceNo, std::vector<int> endSliceNo, int appdimklr, int hogdimklr, MRFEnergy* mrf)
{
  int depth = validdisp.size();
  CShape sh = validdisp[0].Shape();

  int width = sh.width, height = sh.height;


  int klr = 1; // if this function is called this is true

  int intensityKlr = (int)featureCodeKlr/10000;
  int hogKlr = (int) (featureCodeKlr/1000)%10;
  int appKlr = (int) (featureCodeKlr/100)%10;
  int locKlr = (int) (featureCodeKlr/10)%10;

  int dimW = wparam.size2();

  int dim = Xtrain.size2();
  int Ntr = wparam.size1();


  int nW = (int)distW.size();

  for (int kk=0; kk<nW; kk++) {
	  distW[kk] = 0;
  }

  int zslice = 0;
  zslice = startSliceNo[index];

  if (nW == 0)
	return; // nothing to do

  IVM* ivm = new IVM(&Xtrain, &wparam, RBF, lambda, nD, p1);;


  for (int z = 0; z < depth; z++) {


	  matrix<DBL_TYPE> hogdescs;

	  if (hogKlr==2)
	  {
		  // read file
		  char hogfile[100];

		  sprintf(hogfile, "/p/imageseg/medicalSeg/parameters/hog/hog2dcsv_90n/hog2ddata_%d.csv", z+zslice );

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

		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *tdi = &validdisp[z].Pixel(0, y, 0);

		for (int x = 0; x < width; x++) {

			unsigned short pix1 = im1[z].Pixel(x, y, 0);

			float *belief = new float[nD];

			int pixel = z*width*height + y*width + x;
			mrf->getBelief(pixel,belief);


			  ublas::vector<DBL_TYPE> Xtest(dim);

			  int idxx = 0;

			  if (intensityKlr==2) // for now only intensity case
			  {
                  getCostIntensityKlr(Xtest, (double)pix1, minmaxM, idxx);
			  }

			  if (locKlr==2)
			  {
				  // 25 values. // this should be a parameter
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


			  for (int j=0; j<Ntr; j++)
			  {
				  double Kval = ivm->classify_klrexp(&Xtest, j);

				  for (int d=0; d<nD; d++)
				  {
					int thetaid = j*dimW + d;

					distW[thetaid] += exp(belief[d])*Kval;
				  }
			  }
			}
		  }
	  }

  delete ivm;
}






// No occlusion model
// return data cost at non-occluded pixels

// grayscale
void computeEmpirDistU(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distU, int ignoreVal, std::vector<CImage2> im1, int nD, MRFEnergy* mrf, std::vector<float> thetaU)
{
  int depth = validdisp.size();
  CShape sh = validdisp[0].Shape();

  int width = sh.width, height = sh.height;

  int nU = (int)distU.size();

  for (int kk=0; kk<nU; kk++) {
	  distU[kk] = 0;
  }

  if (nU == 0)
	return; // nothing to do

  float *belief = new float[nD];

  int nbins = nU/nD;
  //  double binner = 256.0/nbins;

  int pixel=0;

  for (int z = 0; z < depth; z++) {

	  for (int y = 0; y < height; y++) {

		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *tdi = &validdisp[z].Pixel(0, y, 0);

		for (int x = 0; x < width; x++) {

		  unsigned short pix1 = im1[z].Pixel(x, y, 0);

		  mrf->getBelief(pixel,belief);

		  for (int d=0; d<nD; d++)
		  {

            int thetaid;
            getCostIntensityDiscBins(pix1, nbins, d, thetaU, thetaid);
            /*
//			  std::cout << " %f  " << exp(belief[d]);
			  int thetaid = d*nbins + binid;
            */
			  // for us the costs are the theta param, but in differentiation, they vanish,
			  // so its just whether the feature is enabled or not.
			  // move costs pointer to next pixel

			  distU[thetaid] += exp(belief[d]); //*costs[d];
		  }

		  pixel++;

		}
	  }
  }

  delete [] belief;

}


// RGB
void computeEmpirDistU(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distU, int ignoreVal, std::vector<CByteImage> im1, int nD, MRFEnergy* mrf, std::vector<float> thetaU, std::vector<CByteImage> interImage, int interactive)
{
  int depth = validdisp.size();
  CShape sh = validdisp[0].Shape();

  int width = sh.width, height = sh.height;

  int nU = (int)distU.size();

  for (int kk=0; kk<nU; kk++) {
	  distU[kk] = 0;
  }

  if (nU == 0)
	return; // nothing to do

  float *belief = new float[nD];

  int nbins = nU/nD;
  //  double binner = 256.0/nbins;

  int pixel=0;

  for (int z = 0; z < depth; z++) {

	  for (int y = 0; y < height; y++) {

		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *tdi = &validdisp[z].Pixel(0, y, 0);

		for (int x = 0; x < width; x++) {

		  unsigned char *pix1 = &im1[z].Pixel(x, y, 0);

			if (interactive==1)
			{
				uchar *inimp = &interImage[z].Pixel(x, y, 0);
				if (*inimp!=0)
					continue;
			}


		  mrf->getBelief(pixel,belief);

		  for (int d=0; d<nD; d++)
		  {

            int thetaidrgb[3]={0,0,0};
            int thetaidYuv[2]={0,0};
            // getCostIntensityDiscBins(pix1, nbins, d, thetaU, thetaidrgb);
            getCostIntensityDiscBins(pix1, nbins, d, thetaU, thetaidYuv);
            /*
//			  std::cout << " %f  " << exp(belief[d]);
			  int thetaid = d*nbins + binid;
            */
			  // for us the costs are the theta param, but in differentiation, they vanish,
			  // so its just whether the feature is enabled or not.
			  // move costs pointer to next pixel

			  // distU[thetaidrgb[0]] += exp(belief[d]); //*costs[d];
			  // distU[thetaidrgb[1]] += exp(belief[d]);
			  // distU[thetaidrgb[2]] += exp(belief[d]);

			  distU[thetaidYuv[0]] += exp(belief[d]);
			  distU[thetaidYuv[1]] += exp(belief[d]);

		  }

		  pixel++;

		}
	  }
  }


  delete [] belief;

}



void computeEmpirDistH(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distH, int ignoreVal, std::vector<CByteImage> hogIm, int nD, MRFEnergy* mrf, std::vector<float> thetaH)
{
  int depth = validdisp.size();
  CShape sh = validdisp[0].Shape();

  int width = sh.width, height = sh.height;

  int nH = (int)distH.size();

  for (int kk=0; kk<nH; kk++) {
	  distH[kk] = 0;
  }

  if (nH == 0)
	return; // nothing to do

  float *belief = new float[nD];

  int nHpC = nH/nD;

  int pixel=0;

  for (int z = 0; z < depth; z++) {

	  for (int y = 0; y < height; y++) {

		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *tdi = &validdisp[z].Pixel(0, y, 0);

		for (int x = 0; x < width; x++) {

   			uchar *hogval = &hogIm[z].Pixel(x, y, 0);
          //			int hogbin = (int) *hogval;
			mrf->getBelief(pixel,belief);

			for (int d=0; d<nD; d++)
			{
              int thetaid;
              getCostHoGDiscBins((int) *hogval, d, nHpC, thetaH, thetaid);


              //				int thetaid = d*nHpC + hogbin;

				  // for us the costs are the theta param, but in differentiation, they vanish,
				  // so its just whether the feature is enabled or not.
				  // move costs pointer to next pixel
				distH[thetaid] += exp(belief[d]); //*costs[d];
			}
			pixel++;
		}
	  }
  }

  delete[] belief;

}



void computeEmpirDistM(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distM, int ignoreVal, std::vector<CByteImage> motIm, int nD, MRFEnergy* mrf, std::vector<float> thetaM)
{
  int depth = validdisp.size();
  CShape sh = validdisp[0].Shape();

  int width = sh.width, height = sh.height;

  int nM = (int)distM.size();

  for (int kk=0; kk<nM; kk++) {
	  distM[kk] = 0;
  }

  if (nM == 0)
	return; // nothing to do

  float *belief = new float[nD];

  int nMpC = nM/nD;

  int pixel=0;

  for (int z = 0; z < depth; z++) {

	  for (int y = 0; y < height; y++) {

		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *tdi = &validdisp[z].Pixel(0, y, 0);

		for (int x = 0; x < width; x++) {

   			uchar *motval = &motIm[z].Pixel(x, y, 0);
          //			int hogbin = (int) *hogval;
			mrf->getBelief(pixel,belief);

			for (int d=0; d<nD; d++)
			{
              int thetaid;
              getCostmotDiscBins((int) *motval, d, nMpC, thetaM, thetaid);


              //				int thetaid = d*nHpC + hogbin;

				  // for us the costs are the theta param, but in differentiation, they vanish,
				  // so its just whether the feature is enabled or not.
				  // move costs pointer to next pixel
				distM[thetaid] += exp(belief[d]); //*costs[d];
			}
			pixel++;
		}
	  }
  }

  delete[] belief;

}








void computeEmpirDistL(int loc, int index, std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distL, int ignoreVal, int nD, std::vector <std::vector <CImage2> > locdirImage, int numTrainingPats, MRFEnergy* mrf, LocationMV** LocationMVclass, std::vector<float> thetaL, int genparam, std::vector<int> startSliceNo)
{
  int depth = validdisp.size();
  CShape sh = validdisp[0].Shape();

  int width = sh.width, height = sh.height;

  int nL = (int)distL.size();

  for (int kk=0; kk<nL; kk++) {
	  distL[kk] = 0;
  }

  if (nL == 0)
	return; // nothing to do

  int nLpC = nL/nD;

  int zslice=0;
  zslice = startSliceNo[index];

  float *belief = new float[nD];

  int pixel=0;
  int thetaid;

  for (int z = 0; z < depth; z++) {

	  for (int y = 0; y < height; y++) {

		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *tdi = &validdisp[z].Pixel(0, y, 0);

		for (int x = 0; x < width; x++) {

			if (loc==1) {

				  mrf->getBelief(pixel,belief);

                  for (int d=0; d<nD; d++)
				  {

						unsigned short locpix = locdirImage[d][z].Pixel(x,y,0);

                        getCostLocationDiscCube((unsigned short) locdirImage[d][z].Pixel(x,y,0), x, y, z, loccubex, loccubey, loccubez, numTrainingPats, nLpC, thetaL, width, height, depth, d, thetaid);

						double locpixd = (double)locpix/ (loccubex * loccubey * loccubez * numTrainingPats);

						// need to check this for boundary conditions
						// especially when loccube? not multiple of size
						// or when it is actually multiple of size of image volume
						
                        // int thetaid = nLpC*d + hno*maxwidth*maxheight + rowno*maxwidth + colno;

						// the exp(belief[d]) are normalized and so are a prob distribution
						distL[thetaid] += exp(belief[d]) * (1-locpixd);

				  }
			} else if (loc==6 || loc==8) {

				mrf->getBelief(pixel,belief);

				  for (int d=0; d<nD; d++)
				  {

                    getCostLocationDiscGaussian(LocationMVclass, z, zslice, x, y, thetaL, d, loc, genparam, thetaid);
             
						// the exp(belief[d]) are normalized and so are a prob distribution
					  distL[thetaid] += exp(belief[d]);
				  }

			} else if (loc==7 || loc==9)  //CRF soft or hard assignment redefined
			{
				  mrf->getBelief(pixel,belief);

             	  for (int d=0; d<nD; d++)
				  {
					  // find the theta corresponding to c
                    //					  int thetaid = totClusters*d + maxclass;
                    getCostLocationDiscGaussian2(LocationMVclass, z, zslice, x, y, thetaL, d, loc, genparam, nD, thetaid);
  
                    // WORK-THIS
                    //distL[thetaid]++;
                    distL[thetaid] += exp(belief[d]);
				  }

			} else if (loc==7 || loc==9)  //CRF soft or hard assignment redefined
			{
				  mrf->getBelief(pixel,belief);

             	  for (int d=0; d<nD; d++)
				  {
					  // find the theta corresponding to c
                    //					  int thetaid = totClusters*d + maxclass;
                    getCostLocationDiscGaussian2(LocationMVclass, z, zslice, x, y, thetaL, d, loc, genparam, nD, thetaid);
  
                    // WORK-THIS
                    //distL[thetaid]++;
                    
                    distL[thetaid] += exp(belief[d]);

				  }

			  }

		  pixel++;

		}
	  }
  }

  delete[] belief;

}



// Occlusion model
// return expected data cost
void computeModelDistU(CByteImage disp, fvec &distU, MRFEnergy* mrf, int ignoreVal)
{

  CShape sh = disp.Shape();


  int width = sh.width, height = sh.height;

  int nU = (int)distU.size();

  if (nU == 0)
	return; // nothing to do


  float *costs;
  int nPixels, nStates;

  getCostArray(costs, nPixels, nStates);

  assert(nPixels == width * height);

  float totCost = 0;
  float occCost = 0;

  float *belief = new float[nStates];

  for (int y = 0; y < height; y++) {

	uchar *di = &disp.Pixel(0, y, 0);

	for (int x = 0; x < width; x++) {

      if ((int)(di[x])!=ignoreVal) {

        mrf->getBelief(mrf->pixelIndex(y,x),belief);

        // cost for non-occluded pixel
        for (int b=0 ; b<nStates -1 ; b++)
          totCost += exp(belief[b])*costs[b];

        // cost for occluded pixel
        occCost += exp(belief[nStates-1]);

      }

      // move costs pointer to next pixel
      costs += nStates;
	}
  }

  delete[] belief;

  distU[0] = occCost;

  if(nU == 2)
    distU[1] = totCost;

}

// No occlusion model
void computeModelDistU(CByteImage validdisp, CByteImage disp, fvec &distU, MRFEnergy* mrf, int ignoreVal)
{

  CShape sh = validdisp.Shape();


  int width = sh.width, height = sh.height;

  int nU = (int)distU.size();

  if (nU == 0)
	return; // nothing to do


  float *costs;
  int nPixels, nStates;

  getCostArray(costs, nPixels, nStates);

  assert(nPixels == width * height);

  float totCost = 0;
  float occCost = 0;

  float *belief = new float[nStates];

  for (int y = 0; y < height; y++) {

    uchar *tdi = &validdisp.Pixel(0, y, 0);

	for (int x = 0; x < width; x++) {

      if ((int)(tdi[x])!=ignoreVal && tdi[x]!=0) { // skip occluded pixels


        mrf->getBelief(mrf->pixelIndex(y,x),belief);

        // cost for non-occluded pixel
        for (int b=0 ; b<nStates -1 ; b++)
          totCost += exp(belief[b])*costs[b];

        // cost for occluded pixel
        occCost += exp(belief[nStates-1]);
      }

      // move costs pointer to next pixel
      costs += nStates;
	}
  }

  delete[] belief;

  distU[0] = occCost;

  if(nU == 2)
    distU[1] = totCost;

}


// return smoothness cost component for each thetaV parameter at nonoccluded pixels
void computeEmpirDistVPI(CByteImage disp, CByteImage im1grad, fvec &distV, fvec &distP, uchar occlVal, int ignoreVal)
{
	CShape sh = disp.Shape();
	int width = sh.width, height = sh.height;
	int nV = (int)distV.size();
	int nP = (int)distP.size();

	if (nV == 0 && nP == 0)
		return; // nothing to do

	for (int k=0; k < nV; k++)
		distV[k] = 0;

	for (int k=0; k < nP; k++)
		distP[k] = 0;

	for (int y = 0; y < height-1; y++) {

		uchar *di = &disp.Pixel(0, y, 0);
		uchar *di2 = &disp.Pixel(0, y+1, 0);
		uchar *gr = &im1grad.Pixel(0, y, 0);

		for (int x = 0; x < width-1; x++) {

			int g = gr[x]; // quantized gradient at the pixel

            if ((int)(di[x])==ignoreVal)
              continue;

			// Horizontal Edge

			if (di[x] == occlVal && di[x+1] == occlVal){ // both are occluded
              distP[g]++;
			}
			else if (di[x] == occlVal || di[x+1] == occlVal){ // one is occluded
              distP[2*nV + g]++;
			}
			else if ((int)(di[x+1])!=ignoreVal) { // neither is occluded

					int dh = di[x] - di[x+1]; // horizontal disparity
					dh = ABS(dh); // disparity magnitude
					if (dh > 0){   // count if STRICTLY positive
						distV[g]++;

					}
			}

			// Vertical Edge

			if (di[x] == occlVal && di2[x] == occlVal){ // both are occluded
              distP[nV + g]++;
			}
			else if(di[x] == occlVal || di2[x] == occlVal){ // one is occluded
              distP[3*nV + g]++;
			}
			else if ((int)(di2[x])!=ignoreVal) { // neither is occluded

					int dv = di[x] - di2[x]; // vertical disparity
					dv = ABS(dv); // disparity magnitude
					if (dv > 0){ // count if STRICTLY positive
						distV[g]++;
					}
			}
		} // end for x
    } // for y
}


// Simple occlusion model
// return smoothness cost component for each thetaV parameter at nonoccluded pixels
void computeEmpirDistV(CByteImage disp, CByteImage im1grad, fvec &distV, fvec &distP, uchar occlVal, int ignoreVal)
{
	CShape sh = disp.Shape();
	int width = sh.width, height = sh.height;
	int nV = (int)distV.size();
	int nP = (int)distP.size();

	if (nV == 0 && nP==0)
		return; // nothing to do

	for (int k=0; k < nV; k++)
		distV[k] = 0;

	for (int k=0; k < nP; k++)
		distP[k] = 0;

    int cnt = 0;
	for (int y = 0; y < height-1; y++) {
		uchar *di = &disp.Pixel(0, y, 0);
		uchar *di2 = &disp.Pixel(0, y+1, 0);
		uchar *gr = &im1grad.Pixel(0, y, 0);

		for (int x = 0; x < width-1; x++) {

          if ((int)(di[x])==ignoreVal)
            continue;

			int g = gr[x]; // quantized gradient at the pixel

			// Horizontal Edge
			if (di[x] == occlVal && di[x+1] == occlVal){ // both are occluded
				distP[0] ++;
			}
			else if (di[x] == occlVal || di[x+1] == occlVal){ // one is occluded
				distP[1] ++;
			}
			else if ((int)(di[x+1])!=ignoreVal) { // neither is occluded
					int dh = di[x] - di[x+1]; // horizontal disparity
					dh = ABS(dh); // disparity magnitude
					if (dh > 0){   // count if STRICTLY positive
						distV[g]++;
					}
                    else
                      cnt++;
			}

			// Vertical Edge

			if (di[x] == occlVal && di2[x] == occlVal){ // both are occluded
		 		distP[0] ++;
			}
			else if(di[x] == occlVal || di2[x] == occlVal){ // one is occluded
				distP[1] ++;
			}
			else if ((int)(di2[x])!=ignoreVal){
					int dv = di[x] - di2[x]; // vertical disparity
					dv = ABS(dv); // disparity magnitude
					if (dv > 0){ // count if STRICTLY positive
						distV[g]++;
					}
                    else
                      cnt++;

			}
		} // end for x
    } // for y
    std::cout<<"count: " << cnt << "\n";
}


void computeEmpirDistV(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, std::vector<CByteImage> im1grad, fvec &distV, int nD, int ignoreVal, int gradContext, std::vector<int> gradThreshVec, std::vector<CByteImage> interImage, int interactive)
{
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

  for (int z=0; z < depth - 1; z++) {
    for (int y = 0; y < height-1; y++)
      {
        uchar *di = &disp[z].Pixel(0, y, 0);
        uchar *di2 = &disp[z].Pixel(0, y+1, 0);
        uchar *di3 = &disp[z+1].Pixel(0, y, 0);

        for (int x = 0; x < width-1; x++)
          {
            int gh, gv, gd;
            gh	= im1grad[z].Pixel(x, y, 0);
            gv	= im1grad[z].Pixel(x, y, 1);
            gd	= im1grad[z].Pixel(x, y, 2);

    		uchar *inimp = 0;
    		uchar *inimph = 0;
    		uchar *inimp2 = 0;
    		uchar *inimp3 = 0;

    		if (interactive==1)
    		{
    			inimp = &interImage[z].Pixel(x, y, 0);
    			inimph = &interImage[z].Pixel(x+1, y, 0);
    			inimp2 = &interImage[z].Pixel(x, y+1, 0);
    			inimp3 = &interImage[z+1].Pixel(x, y, 0);
    		}

				      
            //CHECK need to make sure it has the gradContext option too
            // i am not fixing this gradOnly option
            if (gradContext==0) 
			  {
                // Horizontal Edge
                if (di[x] != di[x+1])   // this is the border of the parts
                  distV[gh]++;

                // Vertical Edge
                if (di[x] != di2[x]) // this is the border of the parts
                  distV[gv]++;

                // depth edge
                if (di[x] != di3[x]) // this is the border of the parts
                  distV[gd]++;

			  } else { //gradContext option

              // Horizontal Edge
              int d1 = (int) di[x];
              int d2 = (int) di[x+1];

              if (d1>d2) {
                int swap = d2;
                d2 = d1;
                d1 = swap;
              }
            
              if (interactive==0)
			  {
				distV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + gh] ++;
			  } else if (interactive==1)
			  {
				  if (*inimp==0 && *inimph==0)
					  distV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + gh] ++;
			  }



              // Vertical Edge
              d1 = (int) di[x];
              d2 = (int) di2[x];

              if (d1>d2) {
                int swap = d2;
                d2 = d1;
                d1 = swap;
              }

              if (interactive==0)
			  {
            	  distV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + gv] ++;
			  } else if (interactive==1)
			  {
				  if (*inimp==0 && *inimph==0)
					  distV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + gv] ++;
			  }


              // Depth Edge
              d1 = (int) di[x];
              d2 = (int) di3[x];

              if (d1>d2) {
                int swap = d2;
                d2 = d1;
                d1 = swap;
              }


              if (interactive==0)
			  {
            	  distV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + gd] ++;
			  } else if (interactive==1)
			  {
				  if (*inimp==0 && *inimph==0)
					  distV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + gd] ++;
			  }


            }

          }
      }
  }
}



void computeEmpirDistV00(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, std::vector<CByteImage> im1grad, fvec &distV, int nD, int ignoreVal, int gradContext, std::vector<int> gradThreshVec, std::vector<CByteImage> interImage, int interactive)
{
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

  for (int z=0; z < depth ; z++) {
    for (int y = 0; y < height-1; y++)
      {
        uchar *di = &disp[z].Pixel(0, y, 0);
        uchar *di2 = &disp[z].Pixel(0, y+1, 0);

        for (int x = 0; x < width-1; x++)
          {
            int gh, gv, gd;
            gh	= im1grad[z].Pixel(x, y, 0);
            gv	= im1grad[z].Pixel(x, y, 1);

    		uchar *inimp = 0;
    		uchar *inimph = 0;
    		uchar *inimp2 = 0;

    		if (interactive==1)
    		{
    			inimp = &interImage[z].Pixel(x, y, 0);
    			inimph = &interImage[z].Pixel(x+1, y, 0);
    			inimp2 = &interImage[z].Pixel(x, y+1, 0);
    		}


            //CHECK need to make sure it has the gradContext option too
            // i am not fixing this gradOnly option
            if (gradContext==0)
			  {
                // Horizontal Edge
                if (di[x] != di[x+1])   // this is the border of the parts
                  distV[gh]++;

                // Vertical Edge
                if (di[x] != di2[x]) // this is the border of the parts
                  distV[gv]++;

			  } else { //gradContext option

              // Horizontal Edge
              int d1 = (int) di[x];
              int d2 = (int) di[x+1];

              if (d1>d2) {
                int swap = d2;
                d2 = d1;
                d1 = swap;
              }

              if (interactive==0)
			  {
				distV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + gh] ++;
			  } else if (interactive==1)
			  {
				  if (*inimp==0 && *inimph==0)
					  distV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + gh] ++;
			  }



              // Vertical Edge
              d1 = (int) di[x];
              d2 = (int) di2[x];

              if (d1>d2) {
                int swap = d2;
                d2 = d1;
                d1 = swap;
              }

              if (interactive==0)
			  {
            	  distV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + gv] ++;
			  } else if (interactive==1)
			  {
				  if (*inimp==0 && *inimph==0)
					  distV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + gv] ++;
			  }

            }

          }
      }
  }
}














void computeEmpirDistVZ(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, std::vector<CByteImage> im1grad, fvec &distV, fvec &distZ, int nD, int ignoreVal, int gradContext, std::vector<int> gradThreshVec, std::vector<int> gradThreshVecZ, std::vector<CByteImage> interImage, int interactive)
{
  int depth = validdisp.size();
  CShape sh = disp[0].Shape();

  int width = sh.width, height = sh.height;

  int nV = (int)distV.size();
  int nG = gradThreshVec.size() + 1;

  int nZ = (int)distZ.size();
  int nGz = gradThreshVecZ.size() + 1;


  if (nV + nZ == 0)
	return; // nothing to do

  // nb captures the number of bins per pair of distinct organs (background inclusive)
  // nb useful only of gradContext = 1.

  for (int k=0; k < nV; k++)
	distV[k] = 0;

  for (int k=0; k < nZ; k++)
	distZ[k] = 0;


  if (depth==1)
  {

	int z = 0;
	for (int y = 0; y < height-1; y++)
	  {
		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *di2 = &disp[z].Pixel(0, y+1, 0);

		for (int x = 0; x < width-1; x++)
		  {
			int gh, gv, gd;
			gh	= im1grad[z].Pixel(x, y, 0);
			gv	= im1grad[z].Pixel(x, y, 1);

			uchar *inimp = 0;
			uchar *inimph = 0;
			uchar *inimp2 = 0;

			if (interactive==1)
			{
				inimp = &interImage[z].Pixel(x, y, 0);
				inimph = &interImage[z].Pixel(x+1, y, 0);
				inimp2 = &interImage[z].Pixel(x, y+1, 0);
			}


			//CHECK need to make sure it has the gradContext option too
			// i am not fixing this gradOnly option
			if (gradContext==0) // left this out for now
			{

			} else { //gradContext option

			  // Horizontal Edge
			  int d1 = (int) di[x];
			  int d2 = (int) di[x+1];

			  if (d1>d2) {
				int swap = d2;
				d2 = d1;
				d1 = swap;
			  }

			  if (interactive==0)
			  {
				distV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + gh] ++;
			  } else if (interactive==1)
			  {
				  if (*inimp==0 && *inimph==0)
					  distV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + gh] ++;
			  }


			  // Vertical Edge
			  d1 = (int) di[x];
			  d2 = (int) di2[x];

			  if (d1>d2) {
				int swap = d2;
				d2 = d1;
				d1 = swap;
			  }

			  if (interactive==0)
			  {
				  distV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + gv] ++;
			  } else if (interactive==1)
			  {
				  if (*inimp==0 && *inimph==0)
					  distV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + gv] ++;
			  }

			}

		  }
	  }


  } else {

	  for (int z=0; z < depth - 1; z++) {
		for (int y = 0; y < height-1; y++)
		  {
			uchar *di = &disp[z].Pixel(0, y, 0);
			uchar *di2 = &disp[z].Pixel(0, y+1, 0);
			uchar *di3 = &disp[z+1].Pixel(0, y, 0);

			for (int x = 0; x < width-1; x++)
			  {
				int gh, gv, gd;
				gh	= im1grad[z].Pixel(x, y, 0);
				gv	= im1grad[z].Pixel(x, y, 1);
				gd	= im1grad[z].Pixel(x, y, 2);

				uchar *inimp = 0;
				uchar *inimph = 0;
				uchar *inimp2 = 0;
				uchar *inimp3 = 0;

				if (interactive==1)
				{
					inimp = &interImage[z].Pixel(x, y, 0);
					inimph = &interImage[z].Pixel(x+1, y, 0);
					inimp2 = &interImage[z].Pixel(x, y+1, 0);
					inimp3 = &interImage[z+1].Pixel(x, y, 0);
				}


				//CHECK need to make sure it has the gradContext option too
				// i am not fixing this gradOnly option
				if (gradContext==0)
				  {
					// Horizontal Edge
					if (di[x] != di[x+1])   // this is the border of the parts
					  distV[gh]++;

					// Vertical Edge
					if (di[x] != di2[x]) // this is the border of the parts
					  distV[gv]++;

					// depth edge
					if (di[x] != di3[x]) // this is the border of the parts
					  distZ[gd]++;

				  } else { //gradContext option

				  // Horizontal Edge
				  int d1 = (int) di[x];
				  int d2 = (int) di[x+1];

				  if (d1>d2) {
					int swap = d2;
					d2 = d1;
					d1 = swap;
				  }

				  if (interactive==0)
				  {
					distV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + gh] ++;
				  } else if (interactive==1)
				  {
					  if (*inimp==0 && *inimph==0)
						  distV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + gh] ++;
				  }



				  // Vertical Edge
				  d1 = (int) di[x];
				  d2 = (int) di2[x];

				  if (d1>d2) {
					int swap = d2;
					d2 = d1;
					d1 = swap;
				  }

				  if (interactive==0)
				  {
					  distV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + gv] ++;
				  } else if (interactive==1)
				  {
					  if (*inimp==0 && *inimph==0)
						  distV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + gv] ++;
				  }


				  // Depth Edge
				  d1 = (int) di[x];
				  d2 = (int) di3[x];

				  if (d1>d2) {
					int swap = d2;
					d2 = d1;
					d1 = swap;
				  }


				  if (interactive==0)
				  {
					  distZ[(d1*nD + d2 - (d1*(d1+1))/2)*nGz + gd] ++;
				  } else if (interactive==1)
				  {
					  if (*inimp==0 && *inimph==0)
						  distZ[(d1*nD + d2 - (d1*(d1+1))/2)*nGz + gd] ++;
				  }

				}

			  }
		  }
	  }

  }

}







void computeEmpirDistA(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distA, std::vector <std::vector <CByteImage> >  appImage, Appearance** appclass, int ignoreVal, int nD, int app, std::vector<float> thetaA)
{

	  int depth = validdisp.size();
	  CShape sh = disp[0].Shape();

	  int width = sh.width, height = sh.height;

	  int nA = (int)distA.size();

	  if (nA == 0)
		return; // nothing to do

	  for (int k=0; k < nA; k++)
		distA[k] = 0;

	int nT = nA/nD; //this gives centers per class

	for (int z=0; z < depth; z++)
	{
		for (int y = 0; y < height; y++)
		{
			uchar *di = &disp[z].Pixel(0, y, 0);

			for (int x = 0; x < width; x++) {
              
              int d = di[x];

              int dm = d;
              if (app==2)
                dm=0;

              int thetaid;
              if (((int)appImage[z][dm].Pixel(x, y, 0) >= nT && app==1) || ((int)appImage[z][dm].Pixel(x, y, 0) >= nA && app==2))
              {
          		continue;
              }
              getCostAppDiscPatch((int)appImage[z][dm].Pixel(x, y, 0), appclass[dm]->getPatchSize(), x, y, width, height, thetaA, d, nT, thetaid);

              // assume the patch size is the same for all classes
              int patchSize = appclass[dm]->getPatchSize();

              if (x>(int)patchSize/2 && x<width-patchSize/2 && y>(int)patchSize/2 && y<height-patchSize/2)
                {
                  // get thetaA id
                  distA[thetaid] ++;
                }
			}
		}
	}
}


void computeEmpirDistA(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distA, std::vector <std::vector <CByteImage> >  appImage, Appearance** appclass, MRFEnergy* mrf, int ignoreVal, int nD, int app, std::vector<float> thetaA)
{

	int depth = validdisp.size();
	CShape sh = disp[0].Shape();

	int width = sh.width, height = sh.height;

	int nA = (int)distA.size();

	if (nA == 0)
		return; // nothing to do

	for (int k=0; k < nA; k++)
		distA[k] = 0;

	float *belief = new float[nD];

	int nT = nA/nD; //this gives centers per class

	int pixel = 0;

	for (int z=0; z < depth; z++)
	{
		for (int y = 0; y < height; y++)
		{
			uchar *di = &disp[z].Pixel(0, y, 0);

			for (int x = 0; x < width; x++, pixel++)
			{
				  mrf->getBelief(pixel,belief);

				  for (int d=0; d<nD; d++)
				  {
					  int dm = d;
					  if (app==2)
						  dm=0;

                        int thetaid;

                        if (((int)appImage[z][dm].Pixel(x, y, 0) >= nT && app==1) || ((int)appImage[z][dm].Pixel(x, y, 0) >= nA && app==2))
                        {
                    		continue;
                        }

                        getCostAppDiscPatch((int)appImage[z][dm].Pixel(x, y, 0), appclass[dm]->getPatchSize(), x, y, width, height, thetaA, d, nT, thetaid);

                        //                        int clustNo = (int)appImage[z][dm].Pixel(x, y, 0);

					  // assume the patch size is the same for all classes
					  int patchSize = appclass[dm]->getPatchSize();

					  if (x>(int)patchSize/2 && x<width-patchSize/2 && y>(int)patchSize/2 && y<height-patchSize/2)
					  {
						  // get thetaA id and add to cost
						  distA[thetaid] += exp(belief[d]);
					  }
				  }
                  //pixel++;
			}
		}
	}

}


void computeEmpirDistC(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distC, int ignoreVal, int nD)
{
  int depth = validdisp.size();
  CShape sh = disp[0].Shape();

  int width = sh.width, height = sh.height;

  int nC = (int)distC.size();

  if (nC == 0)
	return; // nothing to do

  for (int k=0; k < nC; k++)
	distC[k] = 0;

  for (int z=0; z < depth - 1; z++) {
	  for (int y = 0; y < height-1; y++)
		{
		  uchar *di = &disp[z].Pixel(0, y, 0);
		  uchar *di2 = &disp[z].Pixel(0, y+1, 0);
		  uchar *di3 = &disp[z+1].Pixel(0, y, 0);

		  for (int x = 0; x < width-1; x++)
			{

			  int mind, maxd;
			  int g = 0;
			  int dix, djx;

			  // di[x] and di[x+1]
			  dix = di[x];
			  djx = di[x+1];

			  //context cost
			  if (dix < djx) {
				  mind = dix;
				  maxd = djx;
			  } else {
				  mind = djx;
				  maxd = dix;
			  }

			  for (int k=mind-1; k>=0; k--) {
				  g += (nD-(k+1));
			  }
			  g += (maxd-mind-1);

			  distC[g]++;



			  g = 0;
			  // di[x] and di2[x]
			  dix = di[x];
			  djx = di2[x];

			  //context cost
			  if (dix < djx) {
				  mind = dix;
				  maxd = djx;
			  } else {
				  mind = djx;
				  maxd = dix;
			  }

			  for (int k=mind-1; k>=0; k--) {
				  g += (nD-(k+1));
			  }
			  g += (maxd-mind-1);

			  distC[g]++;




			  g = 0;
			  // di[x] and di3[x]
			  dix = di[x];
			  djx = di3[x];

			  //context cost
			  if (dix < djx) {
				  mind = dix;
				  maxd = djx;
			  } else {
				  mind = djx;
				  maxd = dix;
			  }

			  for (int k=mind-1; k>=0; k--) {
				  g += (nD-(k+1));
			  }
			  g += (maxd-mind-1);

			  distC[g]++;



			}
		}
  }
}

void computeEmpirDistC(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distC, MRFEnergy* mrf, int ignoreVal, int nD)
{
  int depth = validdisp.size();
  CShape sh = disp[0].Shape();

  int width = sh.width, height = sh.height;

  int nC = (int)distC.size();

  if (nC == 0)
	return; // nothing to do

  for (int k=0; k < nC; k++)
	distC[k] = 0;


  int nStates = mrf->getNumLabels();

  float totBelEqDisp, belNotEq;

  float *belEqDisp;

  float *belief = new float[nStates*nStates];

  for (int z=0; z < depth - 1; z++)
  {
	  for (int y = 0; y < height-1; y++)
		{
		  for (int x = 0; x < width-1; x++)
			{
			  int p1, p2;

			  // Horizontal Edge
			  p1 = z*(width*height) + y*width + x;
			  p2 = z*(width*height) + y*width + x+1;
			  // don't i need to check the border conditions? as well as the original place
			  // i copied from: CHECK
			  mrf->getEdgeBelief( p1, p2, belief);

	          for (int ii=1; ii<nStates; ii++)
	          {
	        	  for (int jj=ii+1; jj<nStates; jj++)
	        	  {
			          int cnt = 0;
			          for (int k=ii-1; k>=0; k--) {
			        	  cnt += (nStates-(k+1));
			          }
			          cnt += (jj-ii-1);

			          distC[cnt] += exp(belief[ii*nStates + jj]) + exp(belief[jj*nStates + ii]);
	        	  }
	          }


			  // Vertical Edge
			  p1 = z*(width*height) + y*width + x;
			  p2 = z*(width*height) + (y+1)*width + x;
			  // don't i need to check the border conditions? as well as the original place
			  // i copied from: CHECK
			  mrf->getEdgeBelief( p1, p2, belief);

	          for (int ii=1; ii<nStates; ii++)
	          {
	        	  for (int jj=ii+1; jj<nStates; jj++)
	        	  {

			          int cnt = 0;
			          for (int k=ii-1; k>=0; k--) {
			        	  cnt += (nStates-(k+1));
			          }
			          cnt += (jj-ii-1);

			          distC[cnt] += exp(belief[ii*nStates + jj]) + exp(belief[jj*nStates + ii]);
	        	  }
	          }

			  // Depth Edge
			  p1 = z*(width*height) + y*width + x;
			  p2 = (z+1)*(width*height) + y*width + x;
			  // don't i need to check the border conditions? as well as the original place
			  // i copied from: CHECK
			  mrf->getEdgeBelief( p1, p2, belief);

	          for (int ii=1; ii<nStates; ii++)
	          {
	        	  for (int jj=ii+1; jj<nStates; jj++)
	        	  {

			          int cnt = 0;
			          for (int k=ii-1; k>=0; k--) {
			        	  cnt += (nStates-(k+1));
			          }
			          cnt += (jj-ii-1);

			          distC[cnt] += exp(belief[ii*nStates + jj]) + exp(belief[jj*nStates + ii]);
	        	  }
	          }
			}
		}
  }

  delete[] belief;

}



// No occlusion model
// return smoothness cost component for each thetaV parameter at nonoccluded pixels
void computeEmpirDistV(CByteImage validdisp, CByteImage disp, CByteImage im1grad, fvec &distV, int ignoreVal)
{
  CShape sh = disp.Shape();

  int width = sh.width, height = sh.height;

  int nV = (int)distV.size();

  if (nV == 0)
	return; // nothing to do

  for (int k=0; k < nV; k++)
	distV[k] = 0;


  for (int y = 0; y < height-1; y++)
    {
      uchar *di = &disp.Pixel(0, y, 0);
      uchar *di2 = &disp.Pixel(0, y+1, 0);
//chet      uchar *tdi = &validdisp.Pixel(0, y, 0);
//chet      uchar *tdi2 = &validdisp.Pixel(0, y+1, 0);

      uchar *gr = &im1grad.Pixel(0, y, 0);

      for (int x = 0; x < width-1; x++)
        {

//chet          if ((int)(tdi[x])==ignoreVal || tdi[x] == 0) // skip occluded pixels
//chet            continue;

          int g = gr[x]; // quantized gradient at the pixel

          //
          // Horizontal Edge
          //

          // neighbor right is not occluded
//chet          if ((int)(tdi[x+1])!=ignoreVal && tdi[x+1] != 0)
            {

              int dh = di[x] - di[x+1]; // horizontal disparity

              dh = ABS(dh); // disparity magnitude

              if (dh > 0)   // count if STRICTLY positive
                {
                  distV[g]++;
                }
            }

          //
          // Vertical Edge
          //

          // neighbor down is not occluded
//chet          if ((int)(tdi2[x])!=ignoreVal && tdi2[x] != 0)
            {
              int dv = di[x] - di2[x]; // vertical disparity
              dv = ABS(dv); // disparity magnitude

              if (dv > 0) // count if STRICTLY positive
                {
                  distV[g]++;
                }
            }
        }
    }
}

// Gradient-modulated occlusion model
// return expected smoothness cost component for each thetaV parameter
void computeModelDistVPI(CByteImage disp, CByteImage im1grad, fvec &distV, fvec &distP, MRFEnergy* mrf, int ignoreVal)
{
  CShape sh = disp.Shape();

  int width = sh.width, height = sh.height;

  int nV = (int)distV.size();
  int nP = (int)distP.size();

  if (nV == 0 & nP == 0)
	return; // nothing to do

  for (int k=0; k < nV; k++)
	distV[k] = 0;

  for (int k=0; k < nP; k++)
    distP[k] = 0;

  int nStates = mrf->getNumLabels();

  float totBelEqDisp,totBel1Occ,belBothOcc, belNotEq;

  float *belEqDisp;

  float *belief = new float[nStates*nStates];

  for (int y = 0; y < height-1; y++)
    {
      uchar *gr = &im1grad.Pixel(0, y, 0);
      uchar *di = &disp.Pixel(0, y, 0);
      uchar *di2 = &disp.Pixel(0,y+1,0);

      for (int x = 0; x < width-1; x++)
        {

          if ((int)(di[x])==ignoreVal)
            continue;

          int g = gr[x]; // quantized gradient at the pixel

          //
          // Horizontal Edge
          //

          if ((int)(di[x+1])!=ignoreVal) {

            mrf->getEdgeBelief( mrf->pixelIndex(y,x), mrf->pixelIndex(y,x+1),
                                belief);

            totBelEqDisp = 0;
            totBel1Occ = 0;
            belBothOcc = 0;

            belEqDisp = belief;

            // Sum Probability along diagonal to get Pr[d_p==d_q]
            // for non-occluded states
            for (int b=0 ; b<nStates-1 ; b++, belEqDisp+=(nStates+1))
            {
              totBelEqDisp += exp(*belEqDisp);
            }

            // Probability of both occluded
            belBothOcc = exp(*belEqDisp);

            // Sum probability of either pixel being occluded
            // JW: could speed with pointer arithmetic. i'm lazy and error prone.
            for (int b=0 ; b<nStates-1 ; b++)
            {
              totBel1Occ += exp(belief[b*(nStates)+nStates-1]); // [b][nDisps]
              totBel1Occ += exp(belief[(nStates-1)*(nStates)+b]); //[nDisps][b]
            }

            belNotEq = 1-totBelEqDisp-belBothOcc-totBel1Occ;


            if (belBothOcc>0)
              distP[g] += (belBothOcc>1 ? 1 : belBothOcc);


            if (totBel1Occ>0)
              distP[2*nV + g] += (totBel1Occ>1 ? 1 : totBel1Occ);


            if (belNotEq>0)
              distV[g] += (belNotEq>1 ? 1 : belNotEq);
          }


          //
          // Vertical Edge
          //

          if ((int)(di2[x])!=ignoreVal) {

            mrf->getEdgeBelief( mrf->pixelIndex(y,x), mrf->pixelIndex(y+1,x),
                                belief);

            totBelEqDisp = 0;
            totBel1Occ = 0;
            belBothOcc = 0;

            belEqDisp = belief;

            // Sum Probability along diagonal to get Pr[d_p==d_q]
            // for non-occluded states
            for (int b=0 ; b<nStates-1 ; b++, belEqDisp+=(nStates+1))
            {
              totBelEqDisp += exp(*belEqDisp);
            }

            // Probability of both occluded
            belBothOcc = exp(*belEqDisp);

            // Sum probability of either pixel being occluded
            // JW: could speed with pointer arithmetic. i'm lazy and error prone.
            for (int b=0 ; b<nStates-1 ; b++)
            {
              totBel1Occ += exp(belief[b*(nStates)+nStates-1]); // [b][nDisps]
              totBel1Occ += exp(belief[(nStates-1)*(nStates)+b]); //[nDisps][b]
            }

            belNotEq = 1-totBelEqDisp-belBothOcc-totBel1Occ;


          if (belBothOcc>0)
            distP[nV + g] += (belBothOcc>1 ? 1 : belBothOcc);

          if (totBel1Occ>0)
            distP[3*nV + g] += (totBel1Occ>1 ? 1 : totBel1Occ);

          if (belNotEq>0)
            distV[g] += (belNotEq>1 ? 1 : belNotEq);
          }
        }
    }
  delete[] belief;
}

// Simple occlusion model
// return expected smoothness cost component for each thetaV parameter

void computeModelDistV(CByteImage disp, CByteImage im1grad, fvec &distV, fvec &distP, MRFEnergy* mrf, int ignoreVal)
{
  CShape sh = disp.Shape();

  int width = sh.width, height = sh.height;

  int nV = (int)distV.size();
  int nP = (int)distP.size();

  if (nV == 0 & nP == 0)
	return; // nothing to do

  for (int k=0; k < nV; k++)
	distV[k] = 0;

  for (int k=0; k < nP; k++)
    distP[k] = 0;

  int nStates = mrf->getNumLabels();

  float totBelEqDisp,totBel1Occ,belBothOcc, belNotEq;

  float *belEqDisp;

  float *belief = new float[nStates*nStates];

  for (int y = 0; y < height-1; y++)
    {
      uchar *gr = &im1grad.Pixel(0, y, 0);
      uchar *di = &disp.Pixel(0,y,0);
      uchar *di2 = &disp.Pixel(0,y+1,0);

      for (int x = 0; x < width-1; x++)
        {

          if ((int)(di[x])==ignoreVal)
            continue;

          int g = gr[x]; // quantized gradient at the pixel

          //
          // Horizontal Edge
          //

          if ((int)(di[x+1])!=ignoreVal) {

            mrf->getEdgeBelief( mrf->pixelIndex(y,x), mrf->pixelIndex(y,x+1),
                                belief);

            totBelEqDisp = 0;
            totBel1Occ = 0;
            belBothOcc = 0;

            belEqDisp = belief;

            // Sum Probability along diagonal to get Pr[d_p==d_q]
            // for non-occluded states
            for (int b=0 ; b<nStates-1 ; b++, belEqDisp+=(nStates+1))
            {
              totBelEqDisp += exp(*belEqDisp);
            }

            // Probability of both occluded
            belBothOcc = exp(*belEqDisp);

            // Sum probability of either pixel being occluded
            // JW: could speed with pointer arithmetic. i'm lazy and error prone.
            for (int b=0 ; b<nStates-1 ; b++)
            {
              totBel1Occ += exp(belief[b*(nStates)+nStates-1]); // [b][nDisps]
              totBel1Occ += exp(belief[(nStates-1)*(nStates)+b]); //[nDisps][b]
            }

            belNotEq = 1-totBelEqDisp-belBothOcc-totBel1Occ;


            if (belBothOcc>0)
              distP[0] += (belBothOcc>1 ? 1 : belBothOcc);

            if (totBel1Occ>0)
            distP[1] += (totBel1Occ>1 ? 1 : totBel1Occ);

            if (belNotEq>0)
              distV[g] += (belNotEq>1 ? 1 : belNotEq);
          }

          //
          // Vertical Edge
          //

          if ((int)(di2[x])!=ignoreVal) {

            mrf->getEdgeBelief( mrf->pixelIndex(y,x), mrf->pixelIndex(y+1,x),
                                belief);

            totBelEqDisp = 0;
            totBel1Occ = 0;
            belBothOcc = 0;

            belEqDisp = belief;

            // Sum Probability along diagonal to get Pr[d_p==d_q]
          // for non-occluded states
            for (int b=0 ; b<nStates-1 ; b++, belEqDisp+=(nStates+1))
            {
              totBelEqDisp += exp(*belEqDisp);
            }

            // Probability of both occluded
            belBothOcc = exp(*belEqDisp);

            // Sum probability of either pixel being occluded
            // JW: could speed with pointer arithmetic. i'm lazy and error prone.
            for (int b=0 ; b<nStates-1 ; b++)
            {
              totBel1Occ += exp(belief[b*(nStates)+nStates-1]); // [b][nDisps]
              totBel1Occ += exp(belief[(nStates-1)*(nStates)+b]); //[nDisps][b]
            }

            belNotEq = 1-totBelEqDisp-belBothOcc-totBel1Occ;

            if (belBothOcc>0)
              distP[0] += (belBothOcc>1 ? 1 : belBothOcc);

            if (totBel1Occ>0)
              distP[1] += (totBel1Occ>1 ? 1 : totBel1Occ);

            if (belNotEq>0)
              distV[g] += (belNotEq>1 ? 1 : belNotEq);
          }
          }
        }

  delete[] belief;
}

// No occlusion model
// return expected smoothness cost component for each thetaV parameter
// at nonoccluded pixels
void computeModelDistV(CByteImage validdisp, CByteImage disp, CByteImage im1grad, fvec &distV, MRFEnergy* mrf, int ignoreVal)
{
  CShape sh = validdisp.Shape();

  int width = sh.width, height = sh.height;

  int nV = (int)distV.size();

  if (nV == 0)
	return; // nothing to do

  for (int k=0; k < nV; k++)
	distV[k] = 0;

  int nStates = mrf->getNumLabels();

  float totBelEqDisp,belNotEq;

  float *belEqDisp;

  float *belief = new float[nStates*nStates];

  for (int y = 0; y < height-1; y++)
    {
      uchar *gr = &im1grad.Pixel(0, y, 0);

      uchar *tdi = &validdisp.Pixel(0, y, 0);
      uchar *tdi2 = &validdisp.Pixel(0, y+1, 0);

      for (int x = 0; x < width-1; x++)
        {

          int g = gr[x]; // quantized gradient at the pixel

          if ((int)(tdi[x])==ignoreVal || tdi[x]==0) // skip occluded pixels
            continue;

          //
          // Horizontal Edge
          //

          // neighbor right is not occluded
          if ((int)(tdi[x+1])!=ignoreVal && tdi[x+1]!=0) {

            mrf->getEdgeBelief( mrf->pixelIndex(y,x), mrf->pixelIndex(y,x+1), belief);

            totBelEqDisp = 0;

            belEqDisp = belief;

            // Sum Probability along diagonal to get Pr[d_p==d_q]
            for (int b=0 ; b<nStates ; b++, belEqDisp+=(nStates+1))
              {
                totBelEqDisp += exp(*belEqDisp);
              }

            belNotEq = 1-totBelEqDisp;

            if (belNotEq>0)
              distV[g] += (belNotEq>1 ? 1 : belNotEq);
          }

          //
          // Vertical Edge
          //

          // neighbor down is not occluded
          if ((int)(tdi2[x])!=ignoreVal && tdi2[x]!=0) {

            mrf->getEdgeBelief( mrf->pixelIndex(y,x), mrf->pixelIndex(y+1,x), belief);

            totBelEqDisp = 0;

            belEqDisp = belief;

            // Sum Probability along diagonal to get Pr[d_p==d_q]
            for (int b=0 ; b<nStates-1 ; b++, belEqDisp+=(nStates+1))
            {
              totBelEqDisp += exp(*belEqDisp);
            }


            belNotEq = 1-totBelEqDisp;

            if (belNotEq>0)
              distV[g] += (belNotEq>1 ? 1 : belNotEq);

          }
        }
    }

  delete[] belief;
}



// No occlusion model
// return expected smoothness cost component for each thetaV parameter
// at nonoccluded pixels

void computeEmpirDistV(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, std::vector<CByteImage> im1grad, fvec &distV, MRFEnergy* mrf, int ignoreVal, int gradContext, std::vector<int> gradThreshVec, std::vector<CByteImage> interImage, int interactive)
{

  int depth = validdisp.size();
  CShape sh = disp[0].Shape();

  int width = sh.width, height = sh.height;

  int nV = (int)distV.size();
  int nG = gradThreshVec.size() + 1;


  if (nV == 0)
	return; // nothing to do

  for (int k=0; k < nV; k++)
	distV[k] = 0;

  int nD = mrf->getNumLabels();

  float totBelEqDisp, belNotEq;

  float *belEqDisp;

  float *belief = new float[nD*nD];


  for (int z=0; z < depth - 1; z++)
  {
	  for (int y = 0; y < height-1; y++)
	  {
		  for (int x = 0; x < width-1; x++)
		  {

            int gh, gv, gd; // quantized gradient at the pixel

            gh	= im1grad[z].Pixel(x, y, 0);
			gv	= im1grad[z].Pixel(x, y, 1);
			gd	= im1grad[z].Pixel(x, y, 2);

    		uchar *inimp = 0;
    		uchar *inimph = 0;
    		uchar *inimp2 = 0;
    		uchar *inimp3 = 0;

    		if (interactive==1)
    		{
    			inimp = &interImage[z].Pixel(x, y, 0);
    			inimph = &interImage[z].Pixel(x+1, y, 0);
    			inimp2 = &interImage[z].Pixel(x, y+1, 0);
    			inimp3 = &interImage[z+1].Pixel(x, y, 0);
    		}
	  
            // not working on this, will need to be fixed
            if (gradContext==0)
			  {
                // common grad bins for all pair of labels

                int p1, p2;

                // Horizontal Edge
                p1 = z*(width*height) + y*width + x;
                p2 = z*(width*height) + y*width + x+1;
                // don't i need to check the border conditions? as well as the original place
                // i copied from: CHECK
                mrf->getEdgeBelief( p1, p2, belief);

                totBelEqDisp = 0;
                belEqDisp = belief;

                // Sum Probability along diagonal to get Pr[d_p==d_q]
                for (int b=0 ; b<nD ; b++, belEqDisp+=(nD+1))
                  totBelEqDisp += exp(*belEqDisp);

                belNotEq = 1-totBelEqDisp;

                if (belNotEq>0) // meaning that the two pixels with different labels
                  distV[gh] += (belNotEq>1 ? 1 : belNotEq);

                // Vertical Edge
                p1 = z*(width*height) + y*width + x;
                p2 = z*(width*height) + (y+1)*width + x;
                // don't i need to check the border conditions? as well as the original place
                // i copied from: CHECK
                mrf->getEdgeBelief( p1, p2, belief);

                totBelEqDisp = 0;
                belEqDisp = belief;

                // Sum Probability along diagonal to get Pr[d_p==d_q]
                for (int b=0 ; b<nD ; b++, belEqDisp+=(nD+1))
                  totBelEqDisp += exp(*belEqDisp);

                belNotEq = 1-totBelEqDisp;

                if (belNotEq>0) // meaning that the two pixels with different labels
                  distV[gv] += (belNotEq>1 ? 1 : belNotEq);


                // depth edge
				  
                p1 = z*(width*height) + y*width + x;
                p2 = (z+1)*(width*height) + y*width + x;
                // don't i need to check the border conditions? as well as the original place
                // i copied from: CHECK
                mrf->getEdgeBelief( p1, p2, belief);

                totBelEqDisp = 0;
                belEqDisp = belief;

                // Sum Probability along diagonal to get Pr[d_p==d_q]
                for (int b=0 ; b<nD ; b++, belEqDisp+=(nD+1))
                  totBelEqDisp += exp(*belEqDisp);

                belNotEq = 1-totBelEqDisp;

                if (belNotEq>0) // meaning that the two pixels with different labels
                  distV[gd] += (belNotEq>1 ? 1 : belNotEq);

			  } else {

              //different grad bins for each pair of labels - CHECK

              int p1, p2;

              // Horizontal Edge
			
              p1 = z*(width*height) + y*width + x;
              p2 = z*(width*height) + y*width + x+1;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrf->getEdgeBelief( p1, p2, belief);

              /*
              for (int ii=1; ii<nStates; ii++)
                {
                  for (int jj=ii+1; jj<nStates; jj++)
                    {

                      int cnt = 0;
                      for (int k=ii-1; k>=0; k--) {
                        cnt += (nStates-(k+1))*nb;
                      }
                      cnt += (jj-ii-1)*nb + g;

                      distV[cnt] += exp(belief[ii*nStates + jj]) + exp(belief[jj*nStates + ii]);
                    }
                }
              */

              if (interactive==0)
			  {
				  for (int d1=0; d1<nD; d1++) {
					  for (int d2=0; d2<nD; d2++) {

						int d1V = d1;
						int d2V = d2;

						if (d1>d2) {
						  d2V = d1;
						  d1V = d2;
						}

						distV[(d1V*nD + d2V - (d1V*(d1V+1))/2)*nG + gh] +=  exp(belief[d1*nD + d2]);
					  }
				  }
			  } else if (interactive == 1)
			  {
				  if (*inimp==0 && *inimph==0)
				  {
					  for (int d1=0; d1<nD; d1++) {
						  for (int d2=0; d2<nD; d2++) {

							int d1V = d1;
							int d2V = d2;

							if (d1>d2) {
							  d2V = d1;
							  d1V = d2;
							}

							distV[(d1V*nD + d2V - (d1V*(d1V+1))/2)*nG + gh] +=  exp(belief[d1*nD + d2]);
						  }
					  }
				  }
			  }


              // Vertical Edge
			
              p1 = z*(width*height) + y*width + x;
              p2 = z*(width*height) + (y+1)*width + x;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrf->getEdgeBelief( p1, p2, belief);

              /*              for (int ii=1; ii<nStates; ii++)
                {
                  for (int jj=ii+1; jj<nStates; jj++)
                    {

                      int cnt = 0;
                      for (int k=ii-1; k>=0; k--) {
                        cnt += (nStates-(k+1))*nb;
                      }

                      cnt += (jj-ii-1)*nb + g;

                      distV[cnt] += exp(belief[ii*nStates + jj]) + exp(belief[jj*nStates + ii]);
                    }
                }
              */

              if (interactive==0)
			  {
				  for (int d1=0; d1<nD; d1++) {
					  for (int d2=0; d2<nD; d2++) {

						int d1V = d1;
						int d2V = d2;

						if (d1>d2) {
						  d2V = d1;
						  d1V = d2;
						}

						distV[(d1V*nD + d2V - (d1V*(d1V+1))/2)*nG + gv] +=  exp(belief[d1*nD + d2]);
					  }
				  }
			  } else if (interactive==1)
			  {
				  if (*inimp==0 && *inimp2==0)
				  {
					  for (int d1=0; d1<nD; d1++) {
						  for (int d2=0; d2<nD; d2++) {

							int d1V = d1;
							int d2V = d2;

							if (d1>d2) {
							  d2V = d1;
							  d1V = d2;
							}

							distV[(d1V*nD + d2V - (d1V*(d1V+1))/2)*nG + gv] +=  exp(belief[d1*nD + d2]);
						  }
					  }
				  }
			  }



              // CHECK if depth info causes problems in the distV
              // Depth Edge
			
              p1 = z*(width*height) + y*width + x;
              p2 = (z+1)*(width*height) + y*width + x;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrf->getEdgeBelief( p1, p2, belief);

              /*
              for (int ii=1; ii<nStates; ii++)
                {
                  for (int jj=ii+1; jj<nStates; jj++)
                    {

                      int cnt = 0;
                      for (int k=ii-1; k>=0; k--) {
                        cnt += (nStates-(k+1))*nb;
                      }
                      cnt += (jj-ii-1)*nb + g;

                      distV[cnt] += exp(belief[ii*nStates + jj]) + exp(belief[jj*nStates + ii]);
                    }
                }
              */
              if (interactive==0)
			  {
				  for (int d1=0; d1<nD; d1++) {
					for (int d2=0; d2<nD; d2++) {

					  int d1V = d1;
					  int d2V = d2;

					  if (d1>d2) {
						d2V = d1;
						d1V = d2;
					  }

					  distV[(d1V*nD + d2V - (d1V*(d1V+1))/2)*nG + gd] +=  exp(belief[d1*nD + d2]);
					}
				  }
			  } else if (interactive==1)
			  {
				  if (*inimp==0 && *inimp3==0)
				  {
					  for (int d1=0; d1<nD; d1++) {
						for (int d2=0; d2<nD; d2++) {

						  int d1V = d1;
						  int d2V = d2;

						  if (d1>d2) {
							d2V = d1;
							d1V = d2;
						  }

						  distV[(d1V*nD + d2V - (d1V*(d1V+1))/2)*nG + gd] +=  exp(belief[d1*nD + d2]);
						}
					  }
				  }
			  }

            }
          }
      }
  }

  delete[] belief;
}




void computeEmpirDistV00(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, std::vector<CByteImage> im1grad, fvec &distV, MRFEnergy* mrf, int ignoreVal, int gradContext, std::vector<int> gradThreshVec, std::vector<CByteImage> interImage, int interactive)
{

  int depth = validdisp.size();
  CShape sh = disp[0].Shape();

  int width = sh.width, height = sh.height;

  int nV = (int)distV.size();
  int nG = gradThreshVec.size() + 1;


  if (nV == 0)
	return; // nothing to do

  for (int k=0; k < nV; k++)
	distV[k] = 0;

  int nD = mrf->getNumLabels();

  float totBelEqDisp, belNotEq;

  float *belEqDisp;

  float *belief = new float[nD*nD];


  for (int z=0; z < depth; z++)
  {
	  for (int y = 0; y < height-1; y++)
	  {
		  for (int x = 0; x < width-1; x++)
		  {

            int gh, gv; // quantized gradient at the pixel

            gh	= im1grad[z].Pixel(x, y, 0);
			gv	= im1grad[z].Pixel(x, y, 1);

    		uchar *inimp = 0;
    		uchar *inimph = 0;
    		uchar *inimp2 = 0;


    		if (interactive==1)
    		{
    			inimp = &interImage[z].Pixel(x, y, 0);
    			inimph = &interImage[z].Pixel(x+1, y, 0);
    			inimp2 = &interImage[z].Pixel(x, y+1, 0);

    		}

            // not working on this, will need to be fixed
            if (gradContext==0)
			  {
                // common grad bins for all pair of labels

                int p1, p2;

                // Horizontal Edge
                p1 = z*(width*height) + y*width + x;
                p2 = z*(width*height) + y*width + x+1;
                // don't i need to check the border conditions? as well as the original place
                // i copied from: CHECK
                mrf->getEdgeBelief( p1, p2, belief);

                totBelEqDisp = 0;
                belEqDisp = belief;

                // Sum Probability along diagonal to get Pr[d_p==d_q]
                for (int b=0 ; b<nD ; b++, belEqDisp+=(nD+1))
                  totBelEqDisp += exp(*belEqDisp);

                belNotEq = 1-totBelEqDisp;

                if (belNotEq>0) // meaning that the two pixels with different labels
                  distV[gh] += (belNotEq>1 ? 1 : belNotEq);

                // Vertical Edge
                p1 = z*(width*height) + y*width + x;
                p2 = z*(width*height) + (y+1)*width + x;
                // don't i need to check the border conditions? as well as the original place
                // i copied from: CHECK
                mrf->getEdgeBelief( p1, p2, belief);

                totBelEqDisp = 0;
                belEqDisp = belief;

                // Sum Probability along diagonal to get Pr[d_p==d_q]
                for (int b=0 ; b<nD ; b++, belEqDisp+=(nD+1))
                  totBelEqDisp += exp(*belEqDisp);

                belNotEq = 1-totBelEqDisp;

                if (belNotEq>0) // meaning that the two pixels with different labels
                  distV[gv] += (belNotEq>1 ? 1 : belNotEq);


			  } else {

              //different grad bins for each pair of labels - CHECK

              int p1, p2;

              // Horizontal Edge

              p1 = z*(width*height) + y*width + x;
              p2 = z*(width*height) + y*width + x+1;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrf->getEdgeBelief( p1, p2, belief);


              if (interactive==0)
			  {
				  for (int d1=0; d1<nD; d1++) {
					  for (int d2=0; d2<nD; d2++) {

						int d1V = d1;
						int d2V = d2;

						if (d1>d2) {
						  d2V = d1;
						  d1V = d2;
						}

						distV[(d1V*nD + d2V - (d1V*(d1V+1))/2)*nG + gh] +=  exp(belief[d1*nD + d2]);
					  }
				  }
			  } else if (interactive == 1)
			  {
				  if (*inimp==0 && *inimph==0)
				  {
					  for (int d1=0; d1<nD; d1++) {
						  for (int d2=0; d2<nD; d2++) {

							int d1V = d1;
							int d2V = d2;

							if (d1>d2) {
							  d2V = d1;
							  d1V = d2;
							}

							distV[(d1V*nD + d2V - (d1V*(d1V+1))/2)*nG + gh] +=  exp(belief[d1*nD + d2]);
						  }
					  }
				  }
			  }


              // Vertical Edge

              p1 = z*(width*height) + y*width + x;
              p2 = z*(width*height) + (y+1)*width + x;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrf->getEdgeBelief( p1, p2, belief);



              if (interactive==0)
			  {
				  for (int d1=0; d1<nD; d1++) {
					  for (int d2=0; d2<nD; d2++) {

						int d1V = d1;
						int d2V = d2;

						if (d1>d2) {
						  d2V = d1;
						  d1V = d2;
						}

						distV[(d1V*nD + d2V - (d1V*(d1V+1))/2)*nG + gv] +=  exp(belief[d1*nD + d2]);
					  }
				  }
			  } else if (interactive==1)
			  {
				  if (*inimp==0 && *inimp2==0)
				  {
					  for (int d1=0; d1<nD; d1++) {
						  for (int d2=0; d2<nD; d2++) {

							int d1V = d1;
							int d2V = d2;

							if (d1>d2) {
							  d2V = d1;
							  d1V = d2;
							}

							distV[(d1V*nD + d2V - (d1V*(d1V+1))/2)*nG + gv] +=  exp(belief[d1*nD + d2]);
						  }
					  }
				  }
			  }

            }
          }
      }
  }

  delete[] belief;
}














void computeEmpirDistVZ(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, std::vector<CByteImage> im1grad, fvec &distV, fvec &distZ, MRFEnergy* mrf, int ignoreVal, int gradContext, std::vector<int> gradThreshVec,  std::vector<int> gradThreshVecZ, std::vector<CByteImage> interImage, int interactive)
{

  int depth = validdisp.size();
  CShape sh = disp[0].Shape();

  int width = sh.width, height = sh.height;

  int nV = (int)distV.size();
  int nG = gradThreshVec.size() + 1;

  int nZ = (int)distZ.size();
  int nGz = gradThreshVecZ.size() + 1;

  if (nV + nZ == 0)
	return; // nothing to do

  for (int k=0; k < nV; k++)
	distV[k] = 0;

  for (int k=0; k < nZ; k++)
	distZ[k] = 0;


  int nD = mrf->getNumLabels();

  float totBelEqDisp, belNotEq;

  float *belEqDisp;

  float *belief = new float[nD*nD];


  for (int z=0; z < depth - 1; z++)
  {
	  for (int y = 0; y < height-1; y++)
	  {
		  for (int x = 0; x < width-1; x++)
		  {

            int gh, gv, gd; // quantized gradient at the pixel

            gh	= im1grad[z].Pixel(x, y, 0);
			gv	= im1grad[z].Pixel(x, y, 1);
			gd	= im1grad[z].Pixel(x, y, 2);

    		uchar *inimp = 0;
    		uchar *inimph = 0;
    		uchar *inimp2 = 0;
    		uchar *inimp3 = 0;

    		if (interactive==1)
    		{
    			inimp = &interImage[z].Pixel(x, y, 0);
    			inimph = &interImage[z].Pixel(x+1, y, 0);
    			inimp2 = &interImage[z].Pixel(x, y+1, 0);
    			inimp3 = &interImage[z+1].Pixel(x, y, 0);
    		}

            // not working on this, will need to be fixed
            if (gradContext==0)
			  {
                // common grad bins for all pair of labels

                int p1, p2;

                // Horizontal Edge
                p1 = z*(width*height) + y*width + x;
                p2 = z*(width*height) + y*width + x+1;
                // don't i need to check the border conditions? as well as the original place
                // i copied from: CHECK
                mrf->getEdgeBelief( p1, p2, belief);

                totBelEqDisp = 0;
                belEqDisp = belief;

                // Sum Probability along diagonal to get Pr[d_p==d_q]
                for (int b=0 ; b<nD ; b++, belEqDisp+=(nD+1))
                  totBelEqDisp += exp(*belEqDisp);

                belNotEq = 1-totBelEqDisp;

                if (belNotEq>0) // meaning that the two pixels with different labels
                  distV[gh] += (belNotEq>1 ? 1 : belNotEq);

                // Vertical Edge
                p1 = z*(width*height) + y*width + x;
                p2 = z*(width*height) + (y+1)*width + x;
                // don't i need to check the border conditions? as well as the original place
                // i copied from: CHECK
                mrf->getEdgeBelief( p1, p2, belief);

                totBelEqDisp = 0;
                belEqDisp = belief;

                // Sum Probability along diagonal to get Pr[d_p==d_q]
                for (int b=0 ; b<nD ; b++, belEqDisp+=(nD+1))
                  totBelEqDisp += exp(*belEqDisp);

                belNotEq = 1-totBelEqDisp;

                if (belNotEq>0) // meaning that the two pixels with different labels
                  distV[gv] += (belNotEq>1 ? 1 : belNotEq);


                // depth edge

                p1 = z*(width*height) + y*width + x;
                p2 = (z+1)*(width*height) + y*width + x;
                // don't i need to check the border conditions? as well as the original place
                // i copied from: CHECK
                mrf->getEdgeBelief( p1, p2, belief);

                totBelEqDisp = 0;
                belEqDisp = belief;

                // Sum Probability along diagonal to get Pr[d_p==d_q]
                for (int b=0 ; b<nD ; b++, belEqDisp+=(nD+1))
                  totBelEqDisp += exp(*belEqDisp);

                belNotEq = 1-totBelEqDisp;

                if (belNotEq>0) // meaning that the two pixels with different labels
                  distZ[gd] += (belNotEq>1 ? 1 : belNotEq);

			  } else {

              //different grad bins for each pair of labels - CHECK

              int p1, p2;

              // Horizontal Edge

              p1 = z*(width*height) + y*width + x;
              p2 = z*(width*height) + y*width + x+1;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrf->getEdgeBelief( p1, p2, belief);

              /*
              for (int ii=1; ii<nStates; ii++)
                {
                  for (int jj=ii+1; jj<nStates; jj++)
                    {

                      int cnt = 0;
                      for (int k=ii-1; k>=0; k--) {
                        cnt += (nStates-(k+1))*nb;
                      }
                      cnt += (jj-ii-1)*nb + g;

                      distV[cnt] += exp(belief[ii*nStates + jj]) + exp(belief[jj*nStates + ii]);
                    }
                }
              */

              if (interactive==0)
			  {
				  for (int d1=0; d1<nD; d1++) {
					  for (int d2=0; d2<nD; d2++) {

						int d1V = d1;
						int d2V = d2;

						if (d1>d2) {
						  d2V = d1;
						  d1V = d2;
						}

						distV[(d1V*nD + d2V - (d1V*(d1V+1))/2)*nG + gh] +=  exp(belief[d1*nD + d2]);
					  }
				  }
			  } else if (interactive == 1)
			  {
				  if (*inimp==0 && *inimph==0)
				  {
					  for (int d1=0; d1<nD; d1++) {
						  for (int d2=0; d2<nD; d2++) {

							int d1V = d1;
							int d2V = d2;

							if (d1>d2) {
							  d2V = d1;
							  d1V = d2;
							}

							distV[(d1V*nD + d2V - (d1V*(d1V+1))/2)*nG + gh] +=  exp(belief[d1*nD + d2]);
						  }
					  }
				  }
			  }


              // Vertical Edge

              p1 = z*(width*height) + y*width + x;
              p2 = z*(width*height) + (y+1)*width + x;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrf->getEdgeBelief( p1, p2, belief);

              /*              for (int ii=1; ii<nStates; ii++)
                {
                  for (int jj=ii+1; jj<nStates; jj++)
                    {

                      int cnt = 0;
                      for (int k=ii-1; k>=0; k--) {
                        cnt += (nStates-(k+1))*nb;
                      }

                      cnt += (jj-ii-1)*nb + g;

                      distV[cnt] += exp(belief[ii*nStates + jj]) + exp(belief[jj*nStates + ii]);
                    }
                }
              */

              if (interactive==0)
			  {
				  for (int d1=0; d1<nD; d1++) {
					  for (int d2=0; d2<nD; d2++) {

						int d1V = d1;
						int d2V = d2;

						if (d1>d2) {
						  d2V = d1;
						  d1V = d2;
						}

						distV[(d1V*nD + d2V - (d1V*(d1V+1))/2)*nG + gv] +=  exp(belief[d1*nD + d2]);
					  }
				  }
			  } else if (interactive==1)
			  {
				  if (*inimp==0 && *inimp2==0)
				  {
					  for (int d1=0; d1<nD; d1++) {
						  for (int d2=0; d2<nD; d2++) {

							int d1V = d1;
							int d2V = d2;

							if (d1>d2) {
							  d2V = d1;
							  d1V = d2;
							}

							distV[(d1V*nD + d2V - (d1V*(d1V+1))/2)*nG + gv] +=  exp(belief[d1*nD + d2]);
						  }
					  }
				  }
			  }



              // CHECK if depth info causes problems in the distV
              // Depth Edge

              p1 = z*(width*height) + y*width + x;
              p2 = (z+1)*(width*height) + y*width + x;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrf->getEdgeBelief( p1, p2, belief);

              /*
              for (int ii=1; ii<nStates; ii++)
                {
                  for (int jj=ii+1; jj<nStates; jj++)
                    {

                      int cnt = 0;
                      for (int k=ii-1; k>=0; k--) {
                        cnt += (nStates-(k+1))*nb;
                      }
                      cnt += (jj-ii-1)*nb + g;

                      distV[cnt] += exp(belief[ii*nStates + jj]) + exp(belief[jj*nStates + ii]);
                    }
                }
              */
              if (interactive==0)
			  {
				  for (int d1=0; d1<nD; d1++) {
					for (int d2=0; d2<nD; d2++) {

					  int d1V = d1;
					  int d2V = d2;

					  if (d1>d2) {
						d2V = d1;
						d1V = d2;
					  }

					  distZ[(d1V*nD + d2V - (d1V*(d1V+1))/2)*nGz + gd] +=  exp(belief[d1*nD + d2]);
					}
				  }
			  } else if (interactive==1)
			  {
				  if (*inimp==0 && *inimp3==0)
				  {
					  for (int d1=0; d1<nD; d1++) {
						for (int d2=0; d2<nD; d2++) {

						  int d1V = d1;
						  int d2V = d2;

						  if (d1>d2) {
							d2V = d1;
							d1V = d2;
						  }

						  distZ[(d1V*nD + d2V - (d1V*(d1V+1))/2)*nGz + gd] +=  exp(belief[d1*nD + d2]);
						}
					  }
				  }
			  }

            }
          }
      }
  }

  delete[] belief;
}







// std::vector element-wise product: prod = a .* b
void vecEltWiseProd(fvec a, fvec b, fvec &prod)
{
  int nK = (int)a.size();
  assert((int)b.size() == nK);
  prod.resize(nK);
  for (int k = 0; k < nK; k++)
	prod[k] = a[k] * b[k];
}

// std::vector element-wise quotient: quot = a ./ b
void vecEltWiseQuot(fvec a, fvec b, fvec &quot)
{
  int nK = (int)a.size();
  assert((int)b.size() == nK);
  quot.resize(nK);
  for (int k = 0; k < nK; k++) {
	  if (a[k]>0 && b[k]==0) // pretend b[k] is a count of 1
		  quot[k] = a[k];
	  else if (a[k]==0 && b[k]==0)
		  quot[k] = 0;
	  else
		  quot[k] = a[k] / b[k];
//		  quot[k] = (a[k] == 0) ? 0 // don't divide by b if a is already 0 (avoids 0/0 case)
//      : a[k] / b[k];
  }
}

// std::vector difference: diff = a - b
void vecDiff(fvec a, fvec b, fvec &diff)
{
  int nK = (int)a.size();
  assert((int)b.size() == nK);
  diff.resize(nK);
  for (int k = 0; k < nK; k++)
	diff[k] = a[k] - b[k];
}

// std::vector sum: sum = a + b
void vecSum(fvec a, fvec b, fvec &sum)
{
  int nK = (int)a.size();
  assert((int)b.size() == nK);
  sum.resize(nK);
  for (int k = 0; k < nK; k++)
	sum[k] = a[k] + b[k];
}


// std::vector sum: sum = a + b
void vecSum(matrix<DBL_TYPE> &wparam, fvec b)
{
	int Ntr = wparam.size1();
	int dim = wparam.size2();

	int cnt = 0;

	for (int i=0; i<Ntr; i++)
		for (int j=0; j<dim; j++) {
			wparam(i,j) += b[cnt];
			cnt++;
		}
}


// std::vector sum: sum = a + b
void matrixSum(matrix<DBL_TYPE> &wparam, matrix<DBL_TYPE> oldwparam)
{
	int Ntr = wparam.size1();
	int dim = wparam.size2();

	for (int i=0; i<Ntr; i++)
		for (int j=0; j<dim; j++) {
			wparam(i,j) += oldwparam(i,j);
		}
}



// std::vector copy: dst = src
void vecCopy(fvec src, fvec &dst)
{
  int nK = (int)src.size();
  dst.resize(nK);
  for (int k = 0; k < nK; k++)
	dst[k] = (float)src[k];
}


void matrixCopy(matrix<DBL_TYPE> src, matrix<DBL_TYPE> &dst)
{

	int dim1 = src.size1();
	int dim2 = src.size2();

	dst.resize(dim1, dim2, 1);

	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
			dst(i,j) = src(i,j);

}


void matrixToVectorCopy(matrix<DBL_TYPE> src, fvec &dst)
{

	int dim1 = src.size1();
	int dim2 = src.size2();

	dst.resize(dim1*dim2);

	int vi = 0;
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
			dst[vi++] = src(i,j);

}


void overwriteVectorToMatrix(fvec src, matrix<DBL_TYPE> &dst)
{

	int dim1 = dst.size1();
	int dim2 = dst.size2();

	int vi = 0;
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
			dst(i,j) = src[vi++];

}



// scale std::vector: t = s* t
void vecScale(fvec &t, float s)
{
  int nK = (int)t.size();
  for (int k = 0; k < nK; k++)
	t[k] *= s;
}


// scale std::vector: t = s* t
void matrixScale(matrix<DBL_TYPE> &t, float s)
{
	int dim1 = t.size1();
	int dim2 = t.size2();

	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
			t(i,j) *= s;
}


// std::vector L2 norm
float vecNorm(fvec a)
{
  float ssum = 0;
  int nK = (int)a.size();
  for (int k = 0; k < nK; k++)
	ssum += a[k] * a[k];

  return sqrt(ssum);
}

float vecNorm(fvec a, int start, int end)
{
  float ssum = 0;
  int nK = (int)a.size();

  assert(end < nK);
  assert(start >= 0);

  for (int k = start; k <= end; k++)
	ssum += a[k] * a[k];

  return sqrt(ssum);
}


// limit std::vector elements to specified range
void vecBound(fvec &a, float minval, float maxval)
{
    int nK = (int)a.size();
    bool flagmin = false;
    bool flagmax = false;
    for (int k = 0; k < nK; k++) {
	if (a[k] < minval) {
	    a[k] = minval;
	    flagmin = true;
	}
	if (a[k] > maxval) {
	    a[k] = maxval;
	    flagmax = true;
	}
    }
    if (flagmin)
	fprintf(stderr, "vecBound: set to minval %g\n", minval);
    if (flagmax)
	fprintf(stderr, "vecBound: set to maxval %g\n", maxval);
}


// print std::vector to stderr
void printfVec(std::vector<int> a, char *name)
{
  int nK = (int)a.size();
  DEBUG_OUT2(verbose, debugfile, "%s: (%s", name, nK==0? ")\n" : "");
  for (int k = 0; k < nK; k++)
	DEBUG_OUT2(verbose, debugfile, "%d%s", a[k], k==nK-1? ")\n" : ",");
}

void printfVec(fvec a, char *name)
{
  int nK = (int)a.size();
  DEBUG_OUT2(verbose, debugfile, "%s: (%s", name, nK==0? ")\n" : "");
  for (int k = 0; k < nK; k++)
	DEBUG_OUT2(verbose, debugfile, "%g%s", a[k], k==nK-1? ")\n" : ",");
}

//negate vectors
void negateVec(fvec &abc)
{
  int nA = (int)abc.size();
  int k;
  for (k = 0; k < nA; k++)
	abc[k] = -abc[k];
}

// concatenate std::vectors
void concatVec(fvec a, fvec b, fvec c, fvec &abc)
{
  int nA = (int)a.size();
  int nB = (int)b.size();
  int nC = (int)c.size();
  abc.resize(nA + nB + nC);
  int k, i = 0;
  for (k = 0; k < nA; k++)
	abc[i++] = a[k];
  for (k = 0; k < nB; k++)
	abc[i++] = b[k];
  for (k = 0; k < nC; k++)
	abc[i++] = c[k];
}

// concatenate std::vectors
void concatVec(fvec a, fvec b, fvec c, fvec d, fvec &abcd)
{
  int nA = (int)a.size();
  int nB = (int)b.size();
  int nC = (int)c.size();
  int nD = (int)d.size();

  abcd.resize(nA + nB + nC + nD);
  int k, i = 0;
  for (k = 0; k < nA; k++)
	abcd[i++] = a[k];
  for (k = 0; k < nB; k++)
	abcd[i++] = b[k];
  for (k = 0; k < nC; k++)
	abcd[i++] = c[k];
  for (k = 0; k < nD; k++)
	abcd[i++] = d[k];
}

void concatVec(fvec a, fvec b, fvec c, fvec d, fvec e, fvec &abcde)
{
  int nA = (int)a.size();
  int nB = (int)b.size();
  int nC = (int)c.size();
  int nD = (int)d.size();
  int nE = (int)e.size();

  abcde.resize(nA + nB + nC + nD + nE);
  int k, i = 0;
  for (k = 0; k < nA; k++)
	abcde[i++] = a[k];
  for (k = 0; k < nB; k++)
	abcde[i++] = b[k];
  for (k = 0; k < nC; k++)
	abcde[i++] = c[k];
  for (k = 0; k < nD; k++)
	abcde[i++] = d[k];
  for (k = 0; k < nE; k++)
	abcde[i++] = e[k];
}



void concatVec(fvec a, fvec b, fvec c, fvec d, fvec e, fvec f, fvec &abcdef)
{
  int nA = (int)a.size();
  int nB = (int)b.size();
  int nC = (int)c.size();
  int nD = (int)d.size();
  int nE = (int)e.size();
  int nF = (int)f.size();

  abcdef.resize(nA + nB + nC + nD + nE + nF);
  int k, i = 0;
  for (k = 0; k < nA; k++)
	abcdef[i++] = a[k];
  for (k = 0; k < nB; k++)
	abcdef[i++] = b[k];
  for (k = 0; k < nC; k++)
	abcdef[i++] = c[k];
  for (k = 0; k < nD; k++)
	abcdef[i++] = d[k];
  for (k = 0; k < nE; k++)
	abcdef[i++] = e[k];
  for (k = 0; k < nF; k++)
	abcdef[i++] = f[k];
}


void concatVec(fvec a, fvec b, fvec c, fvec d, fvec e, fvec f, fvec g, fvec &abcdefg)
{
  int nA = (int)a.size();
  int nB = (int)b.size();
  int nC = (int)c.size();
  int nD = (int)d.size();
  int nE = (int)e.size();
  int nF = (int)f.size();
  int nG = (int)g.size();

  abcdefg.resize(nA + nB + nC + nD + nE + nF + nG);
  int k, i = 0;
  for (k = 0; k < nA; k++)
	abcdefg[i++] = a[k];
  for (k = 0; k < nB; k++)
	abcdefg[i++] = b[k];
  for (k = 0; k < nC; k++)
	abcdefg[i++] = c[k];
  for (k = 0; k < nD; k++)
	abcdefg[i++] = d[k];
  for (k = 0; k < nE; k++)
	abcdefg[i++] = e[k];
  for (k = 0; k < nF; k++)
	abcdefg[i++] = f[k];
  for (k = 0; k < nG; k++)
	abcdefg[i++] = g[k];
}

void concatVec(fvec a, fvec b, fvec c, fvec d, fvec e, fvec f, fvec g, fvec h, fvec &abcdefgh)
{
  int nA = (int)a.size();
  int nB = (int)b.size();
  int nC = (int)c.size();
  int nD = (int)d.size();
  int nE = (int)e.size();
  int nF = (int)f.size();
  int nG = (int)g.size();
  int nH = (int)h.size();

  abcdefgh.resize(nA + nB + nC + nD + nE + nF + nG + nH);
  int k, i = 0;
  for (k = 0; k < nA; k++)
	abcdefgh[i++] = a[k];
  for (k = 0; k < nB; k++)
	abcdefgh[i++] = b[k];
  for (k = 0; k < nC; k++)
	abcdefgh[i++] = c[k];
  for (k = 0; k < nD; k++)
	abcdefgh[i++] = d[k];
  for (k = 0; k < nE; k++)
	abcdefgh[i++] = e[k];
  for (k = 0; k < nF; k++)
	abcdefgh[i++] = f[k];
  for (k = 0; k < nG; k++)
	abcdefgh[i++] = g[k];
  for (k = 0; k < nH; k++)
	abcdefgh[i++] = h[k];
}


void concatVec(fvec a, fvec b, fvec c, fvec d, fvec e, fvec f, fvec g, fvec h, fvec j, fvec &abcdefghj)
{
  int nA = (int)a.size();
  int nB = (int)b.size();
  int nC = (int)c.size();
  int nD = (int)d.size();
  int nE = (int)e.size();
  int nF = (int)f.size();
  int nG = (int)g.size();
  int nH = (int)h.size();
  int nJ = (int)j.size();

  abcdefghj.resize(nA + nB + nC + nD + nE + nF + nG + nH + nJ);
  int k, i = 0;
  for (k = 0; k < nA; k++)
	abcdefghj[i++] = a[k];
  for (k = 0; k < nB; k++)
	abcdefghj[i++] = b[k];
  for (k = 0; k < nC; k++)
	abcdefghj[i++] = c[k];
  for (k = 0; k < nD; k++)
	abcdefghj[i++] = d[k];
  for (k = 0; k < nE; k++)
	abcdefghj[i++] = e[k];
  for (k = 0; k < nF; k++)
	abcdefghj[i++] = f[k];
  for (k = 0; k < nG; k++)
	abcdefghj[i++] = g[k];
  for (k = 0; k < nH; k++)
	abcdefghj[i++] = h[k];
  for (k = 0; k < nJ; k++)
	abcdefghj[i++] = j[k];


}










// split std::vectors
void splitVec(fvec abc, fvec &a, fvec &b, fvec &c)
{
  int nA = (int)a.size();
  int nB = (int)b.size();
  int nC = (int)c.size();

  assert((int)abc.size() == nA + nB + nC);
  int k, i = 0;
  for (k = 0; k < nA; k++)
	a[k] = abc[i++];
  for (k = 0; k < nB; k++)
	b[k] = abc[i++];
  for( k = 0; k < nC; k++)
    c[k] = abc[i++];
}

// split std::vectors
void splitVec(fvec abcd, fvec &a, fvec &b, fvec &c, fvec &d)
{
  int nA = (int)a.size();
  int nB = (int)b.size();
  int nC = (int)c.size();
  int nD = (int)d.size();

  assert((int)abcd.size() == nA + nB + nC + nD);
  int k, i = 0;
  for (k = 0; k < nA; k++)
	a[k] = abcd[i++];
  for (k = 0; k < nB; k++)
	b[k] = abcd[i++];
  for( k = 0; k < nC; k++)
    c[k] = abcd[i++];
  for( k = 0; k < nD; k++)
    d[k] = abcd[i++];
}



void splitVec(fvec abcde, fvec &a, fvec &b, fvec &c, fvec &d, fvec &e)
{
  int nA = (int)a.size();
  int nB = (int)b.size();
  int nC = (int)c.size();
  int nD = (int)d.size();
  int nE = (int)e.size();

  assert((int)abcde.size() == nA + nB + nC + nD + nE);
  int k, i = 0;
  for (k = 0; k < nA; k++)
	a[k] = abcde[i++];
  for (k = 0; k < nB; k++)
	b[k] = abcde[i++];
  for( k = 0; k < nC; k++)
    c[k] = abcde[i++];
  for( k = 0; k < nD; k++)
    d[k] = abcde[i++];
  for( k = 0; k < nE; k++)
    e[k] = abcde[i++];
}


void splitVec(fvec abcdef, fvec &a, fvec &b, fvec &c, fvec &d, fvec &e, fvec &f)
{
  int nA = (int)a.size();
  int nB = (int)b.size();
  int nC = (int)c.size();
  int nD = (int)d.size();
  int nE = (int)e.size();
  int nF = (int)f.size();

  assert((int)abcdef.size() == nA + nB + nC + nD + nE + nF);
  int k, i = 0;
  for (k = 0; k < nA; k++)
	a[k] = abcdef[i++];
  for (k = 0; k < nB; k++)
	b[k] = abcdef[i++];
  for( k = 0; k < nC; k++)
    c[k] = abcdef[i++];
  for( k = 0; k < nD; k++)
    d[k] = abcdef[i++];
  for( k = 0; k < nE; k++)
    e[k] = abcdef[i++];
  for( k = 0; k < nF; k++)
    f[k] = abcdef[i++];
}


void splitVec(fvec abcdefg, fvec &a, fvec &b, fvec &c, fvec &d, fvec &e, fvec &f, fvec &g)
{
  int nA = (int)a.size();
  int nB = (int)b.size();
  int nC = (int)c.size();
  int nD = (int)d.size();
  int nE = (int)e.size();
  int nF = (int)f.size();
  int nG = (int)g.size();

  assert((int)abcdefg.size() == nA + nB + nC + nD + nE + nF + nG);
  int k, i = 0;
  for (k = 0; k < nA; k++)
	a[k] = abcdefg[i++];
  for (k = 0; k < nB; k++)
	b[k] = abcdefg[i++];
  for( k = 0; k < nC; k++)
    c[k] = abcdefg[i++];
  for( k = 0; k < nD; k++)
    d[k] = abcdefg[i++];
  for( k = 0; k < nE; k++)
    e[k] = abcdefg[i++];
  for( k = 0; k < nF; k++)
    f[k] = abcdefg[i++];
  for( k = 0; k < nG; k++)
    g[k] = abcdefg[i++];


}


void splitVec(fvec abcdefgh, fvec &a, fvec &b, fvec &c, fvec &d, fvec &e, fvec &f, fvec &g, fvec &h)
{
  int nA = (int)a.size();
  int nB = (int)b.size();
  int nC = (int)c.size();
  int nD = (int)d.size();
  int nE = (int)e.size();
  int nF = (int)f.size();
  int nG = (int)g.size();
  int nH = (int)h.size();

  assert((int)abcdefgh.size() == nA + nB + nC + nD + nE + nF + nG + nH);
  int k, i = 0;
  for (k = 0; k < nA; k++)
	a[k] = abcdefgh[i++];
  for (k = 0; k < nB; k++)
	b[k] = abcdefgh[i++];
  for( k = 0; k < nC; k++)
    c[k] = abcdefgh[i++];
  for( k = 0; k < nD; k++)
    d[k] = abcdefgh[i++];
  for( k = 0; k < nE; k++)
    e[k] = abcdefgh[i++];
  for( k = 0; k < nF; k++)
    f[k] = abcdefgh[i++];
  for( k = 0; k < nG; k++)
    g[k] = abcdefgh[i++];
  for( k = 0; k < nH; k++)
    h[k] = abcdefgh[i++];

}




void splitVec(fvec abcdefghj, fvec &a, fvec &b, fvec &c, fvec &d, fvec &e, fvec &f, fvec &g, fvec &h, fvec &j)
{
  int nA = (int)a.size();
  int nB = (int)b.size();
  int nC = (int)c.size();
  int nD = (int)d.size();
  int nE = (int)e.size();
  int nF = (int)f.size();
  int nG = (int)g.size();
  int nH = (int)h.size();
  int nJ = (int)j.size();

  assert((int)abcdefghj.size() == nA + nB + nC + nD + nE + nF + nG + nH + nJ);
  int k, i = 0;
  for (k = 0; k < nA; k++)
	a[k] = abcdefghj[i++];
  for (k = 0; k < nB; k++)
	b[k] = abcdefghj[i++];
  for( k = 0; k < nC; k++)
    c[k] = abcdefghj[i++];
  for( k = 0; k < nD; k++)
    d[k] = abcdefghj[i++];
  for( k = 0; k < nE; k++)
    e[k] = abcdefghj[i++];
  for( k = 0; k < nF; k++)
    f[k] = abcdefghj[i++];
  for( k = 0; k < nG; k++)
    g[k] = abcdefghj[i++];
  for( k = 0; k < nH; k++)
    h[k] = abcdefghj[i++];
  for( k = 0; k < nJ; k++)
    j[k] = abcdefghj[i++];


}





void initializeVecToZero(fvec &a)
{
  int nA = (int)a.size();

  for (int k = 0; k < nA; k++)
	a[k] = 0;

}






// concatenate std::vectors
void concatVec(fvec a, fvec b, fvec &ab)
{
  int nA = (int)a.size();
  int nB = (int)b.size();

  ab.resize(nA + nB );
  int k, i = 0;
  for (k = 0; k < nA; k++)
	ab[i++] = a[k];
  for (k = 0; k < nB; k++)
	ab[i++] = b[k];
}

// split std::vectors
void splitVec(fvec ab, fvec &a, fvec &b)
{
  int nA = (int)a.size();
  int nB = (int)b.size();

  assert((int)ab.size() == nA + nB);
  int k, i = 0;
  for (k = 0; k < nA; k++)
	a[k] = ab[i++];
  for (k = 0; k < nB; k++)
	b[k] = ab[i++];
}


// std::vector difference: diff = a - b




// dump parameters

void dumpParameters(fvec thetaU, int nU, int intensity, fvec thetaH, int nH, int hog, fvec thetaA, int nA, int app, fvec thetaL, int nL, int loc, fvec thetaC, int nC, int context, fvec thetaV, int nV, fvec thetaZ, int nZ)
{

  if (nU>0 && intensity==1) 
  {
    for (int ii=0; ii<nU; ii++)
      DEBUG_OUT1(verbose, debugfile, "-u %g ", thetaU[ii]);
    DEBUG_OUT0(verbose, debugfile, "\n");
  }


  if (nH>0 && hog==1) 
  {
    for (int ii=0; ii<nH; ii++)
      DEBUG_OUT1(verbose, debugfile, "-h %g ", thetaH[ii]);
    DEBUG_OUT0(verbose, debugfile, "\n");
  }


  if (nA>0 && app==1) 
  {
    for (int ii=0; ii<nA; ii++)
      DEBUG_OUT1(verbose, debugfile, "-s %g ", thetaA[ii]);
    DEBUG_OUT0(verbose, debugfile, "\n");
  }

  if (nL>0 && (loc==1 || loc==6 || loc==7))
  {
    for (int ii=0; ii<nL; ii++)
      DEBUG_OUT1(verbose, debugfile, "-l %g ", thetaL[ii]);
    DEBUG_OUT0(verbose, debugfile, "\n");
  }

  if (nC>0 && context==1) 
  {
    for (int ii=0; ii<nC; ii++)
      DEBUG_OUT1(verbose, debugfile, "-c %g ", thetaC[ii]);
    DEBUG_OUT0(verbose, debugfile, "\n");
  }

  if (nV>0) 
  {
    for (int ii=0; ii<nV; ii++)
      DEBUG_OUT1(verbose, debugfile, "-v %g ", thetaV[ii]);
    DEBUG_OUT0(verbose, debugfile, "\n");
  }

}


void dumpParameters(Parameters prm)
{



	 int intensity = prm.featurepv.intensity;
	 int hog = prm.featurepv.hog;
	 int mot = prm.featurepv.mot;
	 int app = prm.featurepv.app;
	 int loc =  prm.featurepv.loc;
	 int context = prm.featurepv.context;


	  std::vector<float> thetaU, thetaA, thetaH, thetaM, thetaL, thetaC, thetaV, thetaZ;
		    // copy thetas
	  vecCopy(prm.featurepv.thetaU, thetaU);
	  int nU = prm.featurepv.nU;
	  vecCopy(prm.featurepv.thetaA, thetaA);
	  int nA = prm.featurepv.nA;
	  vecCopy(prm.featurepv.thetaH, thetaH);
	  int nH = prm.featurepv.nH;
	  vecCopy(prm.featurepv.thetaM, thetaM);
	  int nM = prm.featurepv.nM;
	  vecCopy(prm.featurepv.thetaL, thetaL);
	  int nL = prm.featurepv.nL;
	  vecCopy(prm.featurepv.thetaC, thetaC);
	  int nC = prm.featurepv.nC;
	  vecCopy(prm.featurepv.thetaV, thetaV);
	  int nV = prm.featurepv.nV;
	  vecCopy(prm.featurepv.thetaZ, thetaZ);
	  int nZ = prm.featurepv.nZ;




  if (nU>0 && intensity==1)
  {
    for (int ii=0; ii<nU; ii++)
      DEBUG_OUT1(verbose, debugfile, "-u %g ", thetaU[ii]);
    DEBUG_OUT0(verbose, debugfile, "\n");
  }


  if (nH>0 && hog==1)
  {
    for (int ii=0; ii<nH; ii++)
      DEBUG_OUT1(verbose, debugfile, "-h %g ", thetaH[ii]);
    DEBUG_OUT0(verbose, debugfile, "\n");
  }

  if (nM>0 && mot==1)
  {
    for (int ii=0; ii<nM; ii++)
      DEBUG_OUT1(verbose, debugfile, "-q %g ", thetaM[ii]);
    DEBUG_OUT0(verbose, debugfile, "\n");
  }


  if (nA>0 && app==1)
  {
    for (int ii=0; ii<nA; ii++)
      DEBUG_OUT1(verbose, debugfile, "-s %g ", thetaA[ii]);
    DEBUG_OUT0(verbose, debugfile, "\n");
  }

  if (nL>0 && (loc==1 || loc==6 || loc==7))
  {
    for (int ii=0; ii<nL; ii++)
      DEBUG_OUT1(verbose, debugfile, "-l %g ", thetaL[ii]);
    DEBUG_OUT0(verbose, debugfile, "\n");
  }


  if (nC>0 && context==1)
  {
    for (int ii=0; ii<nC; ii++)
      DEBUG_OUT1(verbose, debugfile, "-c %g ", thetaC[ii]);
    DEBUG_OUT0(verbose, debugfile, "\n");
  }


  if (nV>0)
  {
    for (int ii=0; ii<nV; ii++)
      DEBUG_OUT1(verbose, debugfile, "-v %g ", thetaV[ii]);
    DEBUG_OUT0(verbose, debugfile, "\n");
  }

  if (nZ>0)
  {
    for (int ii=0; ii<nZ; ii++)
      DEBUG_OUT1(verbose, debugfile, "-d %g ", thetaZ[ii]);
    DEBUG_OUT0(verbose, debugfile, "\n");
  }

}



int checkRegularityCondition(fvec thetaV, std::vector<int> gradThreshVec, int nD)
{

  int cRC = 1;

  int nV = (int)thetaV.size();
  int nG = gradThreshVec.size() + 1;

  for (int tid = 0; tid < nG; tid++) 
  {
    for (int d1 = 0; d1 < nD-1; d1++) 
    {
      for (int d2 = d1+1; d2 < nD; d2++) 
      {

        double gam00 = thetaV[(d1*nD + d1 - (d1*(d1+1))/2)*nG + tid];
        double gam11 = thetaV[(d2*nD + d2 - (d2*(d2+1))/2)*nG + tid];

        // assuming symmetry
        double gam01 = thetaV[(d1*nD + d2 - (d1*(d1+1))/2)*nG + tid];
        double gam10 = gam01;

        if (gam00 + gam11 > gam10 + gam01) // + 1e-5) kind of reduces cRC failures
        {
          std::cout << "(g, d1, d2)(g00, g11, g01, g10) = (" <<tid <<", " << d1 <<", " << d2 << ")(" << std::scientific << gam00 << ", " << gam11 << ", " << gam01 <<", " << gam10 << ")\n";
          cRC = 0;
        }
      }
    }
  }

  return cRC;
}
