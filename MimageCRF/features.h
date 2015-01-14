/*
 * features.h
 *
 *  Created on: Jul 16, 2010
 *      Author: bhole
 */

#ifndef FEATURES_H_
#define FEATURES_H_

#include <iostream>
#include <vector>
#include "locationMV.h"
#include "appearance.h"
#include "appearanceMV.h"
#include "intensityClass.h"
#include "intensityNB.h"
#include <math.h>
#include "readIVMparam.h"
#include "ivm.hpp"
#include "crfmodel.h"
#include "Parameters.h"









// color models and spaces
int rgb2Lab(int rgb[], double Lab[]);
int Lab2rgb(int rgb[], double Lab[]);
int rgb2Yuv(int rgb[], int Yuv[]);
int Yuv2rgb(int rgb[], int Yuv[]);



// intensity features
//grayscale
double getCostIntensityDiscBins(unsigned int pix1, int nbins, double binner, int range, int d, std::vector<float> thetaU, int &thetaid);
double getCostIntensityDiscBins(unsigned int pix1, int nbins, int d, std::vector<float> thetaU, int &thetaid);
double getCostIntensityGenGaussian(intClass *intObj, unsigned short pix1);
double getCostIntensityGenNB(unsigned int pix1, int nbinsNB, int d, intensityNB* intNBclass);

//RGB
double getCostIntensityDiscBins(unsigned char* pix1, int nbins, int d, std::vector<float> thetaU, int &thetaid);
double getCostIntensityDiscBins(unsigned char* pix1, int nbins, int d, std::vector<float> thetaU, int thetaid[]);

// Yuv
double getCostIntensityDiscBinsYuv(unsigned char* pix1, int nbins, int d, std::vector<float> thetaU, int thetaid[]);
double getCostIntensityDiscBinsuv(unsigned char* pix1, int nbins, int d, std::vector<float> thetaU, int thetaid[]);

// Lab
double getCostIntensityDiscBinsLab(unsigned char* pix1, int nbins, int d, std::vector<float> thetaU, int thetaid[]);



// appearance features

double getCostAppDiscPatch(int clustNo , int patchSize, int x, int y, int width, int height, std::vector<float> thetaA, int d, int nT, int &thetaid);
double getCostAppGenPatch(double appcostp, int dim, int x, int y, int width, int height);

// hog features
double getCostHoGDiscBins(int hogbin, int d, int nHpC, std::vector<float> thetaH, int &thetaid);
double getCostHoGGenBins(double hogcostP, int x, int y, int width, int height, int startPatch);

double getCostHoGDiscBinsTemp(int hogbin, int d, std::vector<float> thetaH, int &thetaid);

// bias features
double getCostBias(int d, std::vector<float> thetaB, int &thetaid);

// volume features
double getCostInverseClassSize(int d, std::vector<float> inverseClassSize, std::vector<float> thetaE, int &thetaid);


// motion features

double getCostmotDiscBins(int motbin, int d, int nMpC, std::vector<float> thetaM, int &thetaid);
double getCostmotGenBins(double motcostP, int x, int y, int width, int height, int startPatch);

// location features

double getCostLocationDiscGaussian(LocationMV** LocationMVclass, unsigned int i, int zslice, int x, int y, std::vector<float> thetaL, int d, int loc, int genparam, int &thetaid);
double getCostLocationDiscGaussian2(LocationMV** LocationMVclass, unsigned int z, int zslice, int x, int y, std::vector<float> thetaL, int d, int loc, int genparam, int nD, int &thetaid);
double getCostLocationGenGaussian(LocationMV* LocationMVclass, unsigned int z, int zslice, int x, int y);
double getCostLocationDiscCube(unsigned short locpix, int x, int y, int i, int loccubex, int loccubey, int loccubez, int numTrainingPats, int nLpC, std::vector<float> thetaL, int width, int height, int depth, int d, int &thetaid);
double getCostLocationDT(unsigned short locpixval, std::vector<float> thetaL, int d, int &thetaid);
double getCostLocationDTgen(unsigned short locpixval, std::vector<float> thetaL, int d, int &thetaid);


// rbm features


// daisy features

// optical flow features
void getIndicatorCostOpticalFlow(int depth, int z, int d, std::vector<float> thetaO, int thetaid[], float indcost[], int numframes, int frameidx);
double getCostOpflowBinsuv(double u, double v, int nbins, int d, std::vector<float> thetaO, int &thetaid, int numframes, int frameidx);


//klr features

void readklrParameters(std::string &klrfName, int &hogdimklr, int &appdimklr, int &locdimklr, int &biasklr, double &lambda, std::string &xfile, std::string &wfile, std::string &mmfile, std::string &kernel, double &p1, kernelType &ker);
void printklrParameters(std::string &klrfName, int hogdimklr, int appdimklr, int locdimklr, int biasklr, double lambda, std::string xfile, std::string wfile, std::string mmfile, std::string kernel, double p1, kernelType ker);
void readklrParameters(klrparams &klrpv);
void printklrParameters(klrparams klrpv);

void getCostLocationKlr(ublas::vector<DBL_TYPE> &Xloctest, unsigned int i, int zslice, int x, int y, int nD, matrix<DBL_TYPE> minmaxM, LocationMV** LocationMVclass, int &idxx, int skipXYZlocklr);
void getCostIntensityKlr(ublas::vector<DBL_TYPE> &Xtest, double pix1, matrix<DBL_TYPE> minmaxM, int &idxx);
void getCostAppKlr(ublas::vector<DBL_TYPE> &Xtest, int x, int y, int width, matrix<DBL_TYPE> minmaxM, int &idxx, int appdim, matrix<DBL_TYPE> *appdescs);
void getCostHoGKlr(ublas::vector<DBL_TYPE> &Xtest, int x, int y, int width, matrix<DBL_TYPE> minmaxM, int &idxx, int hogdim, matrix<DBL_TYPE> *hogdescs);

void getCostIntensityKlrYuv(ublas::vector<DBL_TYPE> &Xtest, unsigned char* pix1, matrix<DBL_TYPE> minmaxM, int &idxx);


// read parameter files
void readLocationMVParameters(LocationMV** &LocationMVclass, int nD);
void readApprearanceParametersCPC(Appearance** &appclass, int nD, std::string &appfName, int colordims);
void readApprearanceParameters(Appearance** &appclass, int nD, std::string &appfName);
void readAppearanceMVParameter(AppearanceMV** &AppMVclass, int nD, std::string &appfName);
void readintClassParameter(intClass** &intObj, int nD); 


// delete
void deleteLocationMVclass(int loc, int nD, LocationMV** &LocationMVclass);
void deleteAppMVclass(int app, int nD, AppearanceMV** &AppMVclass);
void deleteappclass(int app, int nD, Appearance** &appclass);
void deleteintObj(int intensity, int nD, intClass** &intObj);



#endif /* FEATURES_H_ */
