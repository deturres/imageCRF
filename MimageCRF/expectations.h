#ifndef EXPECTATIONS_H
#define EXPECTATIONS_H

//#include <omp.h>
#include "features.h"
#include "appearance.h"
#include "locationMV.h"
#include "MRFEnergy.h"
#include "SyncMeanField.h"
#include <vector>
 
#include "readIVMparam.h"
#include "ivm.hpp"

#include "Parameters.h"


#include "boostHeader.h"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

using boost::tokenizer;
using boost::lexical_cast;
using boost::escaped_list_separator;
using namespace boost::numeric::ublas;



typedef std::vector<float> fvec;


void computeMaskedImage(CByteImage disp, CByteImage trueimage, int ignoreVal=0);

float computeMaskedLikelihood(CByteImage disp, SyncMeanField *mrf);

// Using occlusion model
void computeEmpirDistU(CByteImage disp, fvec &distU, uchar occlVal=0, int ignoreVal=-1);
void computeModelDistU(CByteImage disp, fvec &distU, MRFEnergy* mrf, int ignoreVal=-1);
void computeModelDistU(CByteImage truedisp, CByteImage disp, fvec &distU, MRFEnergy* mrf, int ignoreVal=-1);

// No occlusion model
//void computeEmpirDistU(CByteImage truedisp, CByteImage disp, fvec &distU,
//                       int ignoreVal, CByteImage im1, int nD);

void computeEmpirDistU(std::vector<CByteImage> truedisp, std::vector<CByteImage> disp, fvec &distU,
                       int ignoreVal, std::vector<CImage2> im1, int nD, std::vector<float> thetaU);

void computeEmpirDistU(std::vector<CByteImage> truedisp, std::vector<CByteImage> disp, fvec &distU,
                       int ignoreVal, std::vector<CByteImage> im1, int nD, std::vector<float> thetaU, std::vector<CByteImage> interImage, int interactive);


void computeEmpirDistU(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distU, int ignoreVal, std::vector<CImage2> im1, int nD, MRFEnergy* mrf, std::vector<float> thetaU);

void computeEmpirDistU(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distU, int ignoreVal, std::vector<CByteImage> im1, int nD, MRFEnergy* mrf, std::vector<float> thetaU, std::vector<CByteImage> interImage, int interactive);


void computeEmpirDistW(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distW,  matrix<DBL_TYPE> &wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM, int ignoreVal, std::vector<CImage2> im1, int nD, DBL_TYPE lambda, DBL_TYPE p1, int featureCodeKlr, int index, LocationMV** LocationMVclass, std::vector<int> startSliceNo, std::vector<int> endSliceNo, int appdimklr, int hogdimklr);

void computeEmpirDistWPar(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distW,  matrix<DBL_TYPE> &wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM, int ignoreVal, std::vector<CImage2> im1, int nD, DBL_TYPE lambda, DBL_TYPE p1, int featureCodeKlr, int index, LocationMV** LocationMVclass, std::vector<int> startSliceNo, std::vector<int> endSliceNo, int appdimklr, int hogdimklr);

void computeEmpirDistWPar(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distW,  matrix<DBL_TYPE> &wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM, int ignoreVal, std::vector<CImage2> im1, int nD, DBL_TYPE lambda, DBL_TYPE p1, int featureCodeKlr, int index, LocationMV** LocationMVclass, std::vector<int> startSliceNo, std::vector<int> endSliceNo, int appdimklr, int hogdimklr, MRFEnergy* mrf);

void computeEmpirDistW(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distW,  matrix<DBL_TYPE> &wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM, int ignoreVal, std::vector<CImage2> im1, int nD, DBL_TYPE lambda, DBL_TYPE p1, int featureCodeKlr, int index, LocationMV** LocationMVclass, std::vector<int> startSliceNo, std::vector<int> endSliceNo, int appdimklr, int hogdimklr, MRFEnergy* mrf);


// Gradient-modulated Occlusion model
void computeEmpirDistVPI(CByteImage validdisp, CByteImage im1grad,
                         fvec &distV, fvec &distP,
                         uchar occlVal=0, int ignoreVal=-1);
void computeModelDistVPI(CByteImage disp, CByteImage im1grad,
                         fvec &distV, fvec &distP, MRFEnergy* mrf,
                         int ignoreVal=-1);

// Simple occlusion model
void computeEmpirDistV(CByteImage disp, CByteImage im1grad,
                       fvec &distV, fvec &distP,
                       uchar occlVal=0, int ignoreVal=-1);
void computeModelDistV(CByteImage disp, CByteImage im1grad,
                       fvec &distV, fvec &distP,
                       MRFEnergy* mrf, int ignoreVal=-1);

// No occlusion model
void computeEmpirDistV(CByteImage truedisp, CByteImage disp, CByteImage im1grad,
                       fvec &distV, int ignoreVal=-1);
void computeEmpirDistV(std::vector<CByteImage> truedisp, std::vector<CByteImage> disp, std::vector<CByteImage> im1grad, fvec &distV, int nD, int ignoreVal, int gradContext, std::vector<int> gradThreshVec, std::vector<CByteImage> interImage, int interactive);


void computeEmpirDistV00(std::vector<CByteImage> truedisp, std::vector<CByteImage> disp, std::vector<CByteImage> im1grad, fvec &distV, int nD, int ignoreVal, int gradContext, std::vector<int> gradThreshVec, std::vector<CByteImage> interImage, int interactive);
void computeEmpirDistV00(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, std::vector<CByteImage> im1grad,
						fvec &distV, MRFEnergy* mrf, int ignoreVal, int gradContext, std::vector<int> gradThreshVec, std::vector<CByteImage> interImage, int interactive);


void computeModelDistV(CByteImage truedisp, CByteImage disp, CByteImage im1grad,
                       fvec &distV, MRFEnergy* mrf, int ignoreVal=-1);
void computeEmpirDistV(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, std::vector<CByteImage> im1grad,
						fvec &distV, MRFEnergy* mrf, int ignoreVal, int gradContext, std::vector<int> gradThreshVec, std::vector<CByteImage> interImage, int interactive);
// x-y and z
void computeEmpirDistVZ(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, std::vector<CByteImage> im1grad, fvec &distV, fvec &distZ, int nD, int ignoreVal, int gradContext, std::vector<int> gradThreshVec, std::vector<int> gradThreshVecZ, std::vector<CByteImage> interImage, int interactive);

void computeEmpirDistVZ(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, std::vector<CByteImage> im1grad,
						fvec &distV, fvec &distZ, MRFEnergy* mrf, int ignoreVal, int gradContext, std::vector<int> gradThreshVec, std::vector<int> gradThreshVecZ, std::vector<CByteImage> interImage, int interactive);


void computeEmpirDistC(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp,
		fvec &distC, int ignoreVal, int nD);
void computeEmpirDistC(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp,
		fvec &distC, MRFEnergy* mrf, int ignoreVal, int nD);

void computeEmpirDistA(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp,
		fvec &distA, std::vector <std::vector <CByteImage> >  appImage, Appearance** appclass,
                       int ignoreVal, int nD, int app, std::vector<float> thetaA);
void computeEmpirDistA(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp,
		fvec &distA, std::vector <std::vector <CByteImage> >  appImage, Appearance** appclass, 
                       MRFEnergy* mrf, int ignoreVal, int nD, int app, std::vector<float> thetaA);


void computeEmpirDistH(std::vector<CByteImage> truedisp, std::vector<CByteImage> disp, fvec &distH,
                       int ignoreVal, std::vector<CByteImage> hogIm, int nD, std::vector<float> thetaH);

void computeEmpirDistH(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distH,
                       int ignoreVal, std::vector<CByteImage> hogIm, int nD, MRFEnergy* mrf, std::vector<float> thetaH);

void computeEmpirDistM(std::vector<CByteImage> truedisp, std::vector<CByteImage> disp, fvec &distM,
                       int ignoreVal, std::vector<CByteImage> motIm, int nD, std::vector<float> thetaM);

void computeEmpirDistM(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distM,
                       int ignoreVal, std::vector<CByteImage> motIm, int nD, MRFEnergy* mrf, std::vector<float> thetaM);




void computeEmpirDistL(int loc, int index, std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distL, 
                       int ignoreVal, int nD, std::vector <std::vector <CImage2> > locdirImage, int numTrainingPats, LocationMV** LocationMVclass, std::vector<float> thetaL, int genparam, std::vector<int> startSliceNo);

void computeEmpirDistL(int loc, int index, std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, fvec &distL,
		int ignoreVal, int nD, std::vector <std::vector <CImage2> > locdirImage, int numTrainingPats,
                       MRFEnergy* mrf, LocationMV** LocationMVclass, std::vector<float> thetaL, int genparam, std::vector<int> startSliceNo);



// std::vector element-wise product: prod = a .* b
void vecEltWiseProd(fvec a, fvec b, fvec &prod);

// std::vector element-wise quotient: quot = a ./ b
void vecEltWiseQuot(fvec a, fvec b, fvec &quot);

// std::vector difference: diff = a - b
void vecDiff(fvec a, fvec b, fvec &diff);

// std::vector copy: dst = src
void vecCopy(fvec src, fvec &dst);
void matrixCopy(matrix<DBL_TYPE> src, matrix<DBL_TYPE> &dst);

void matrixToVectorCopy(matrix<DBL_TYPE> src, fvec &dst);
void overwriteVectorToMatrix(fvec src, matrix<DBL_TYPE> &dst);

// std::vector sum: sum = a + b
void vecSum(fvec a, fvec b, fvec &sum);
void vecSum(matrix<DBL_TYPE> &wparam, fvec b);
void matrixSum(matrix<DBL_TYPE> &wparam, matrix<DBL_TYPE> oldwparam);

// scale std::vector: t = s * t
void vecScale(fvec &t, float s);
void matrixScale(matrix<DBL_TYPE> &t, float s);

// std::vector L2 norm
float vecNorm(fvec a);
float vecNorm(fvec a, int start, int end);

// limit std::vector elements to specified range
void vecBound(fvec &a, float minval, float maxval);

// print std::vector to stderr
void printfVec(std::vector<int> a, char *name);
void printfVec(fvec a, char *name);

// negate vectors
void negateVec(fvec &abc);

// concatenate std::vectors
void concatVec(fvec a, fvec b, fvec &ab);
void concatVec(fvec a, fvec b, fvec c, fvec &abc);
void concatVec(fvec a, fvec b, fvec c, fvec d, fvec &abcd);
void concatVec(fvec a, fvec b, fvec c, fvec d, fvec e, fvec &abcde);
void concatVec(fvec a, fvec b, fvec c, fvec d, fvec e, fvec f, fvec &abcdef);
void concatVec(fvec a, fvec b, fvec c, fvec d, fvec e, fvec f, fvec g, fvec &abcdefg);
void concatVec(fvec a, fvec b, fvec c, fvec d, fvec e, fvec f, fvec g, fvec h, fvec &abcdefgh);
void concatVec(fvec a, fvec b, fvec c, fvec d, fvec e, fvec f, fvec g, fvec h, fvec j, fvec &abcdefghj);

// split std::vectors
void splitVec(fvec ab, fvec &a, fvec &b);
void splitVec(fvec abc, fvec &a, fvec &b, fvec &c);
void splitVec(fvec abcd, fvec &a, fvec &b, fvec &c, fvec &d);
void splitVec(fvec abcde, fvec &a, fvec &b, fvec &c, fvec &d, fvec &e);
void splitVec(fvec abcdef, fvec &a, fvec &b, fvec &c, fvec &d, fvec &e, fvec &f);
void splitVec(fvec abcdefg, fvec &a, fvec &b, fvec &c, fvec &d, fvec &e, fvec &f, fvec &g);
void splitVec(fvec abcdefgh, fvec &a, fvec &b, fvec &c, fvec &d, fvec &e, fvec &f, fvec &g, fvec &h);
void splitVec(fvec abcdefghj, fvec &a, fvec &b, fvec &c, fvec &d, fvec &e, fvec &f, fvec &g, fvec &h, fvec &j);

void initializeVecToZero(fvec &a);


void dumpParameters(fvec thetaU, int nU, int intensity, fvec thetaH, int nH, int hog, fvec thetaA, int nA, int app, fvec thetaL, int nL, int loc, fvec thetaC, int nC, int context, fvec thetaV, int nV);
void dumpParameters(Parameters prm);


int checkRegularityCondition(fvec thetaV, std::vector<int> gradThreshVec, int nD);

#endif

