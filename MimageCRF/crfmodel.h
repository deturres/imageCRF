#ifndef CRFMODEL_H
#define CRFMODEL_H

//#include <omp.h>
#include <libconfig.h++>

#include "mrf.h"
#include "imageLib.h"
#include <vector>
#include "helperLib.h"
#include "locationMV.h"
#include "appearance.h"
#include "appearanceMV.h"
#include "intensityClass.h"
#include "intensityNB.h"

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


using namespace libconfig;



typedef std::vector<float> fvec;

class ImageWriter;

// global variables for debugging output
extern int verbose;
extern FILE *debugfile;

#ifdef _MSC_VER        // Visual C++
// unfortunately, variable argument macros not supported until VCC 2005...
#define DEBUG_OUT0(verbose, debugfile, fmt) { \
  if (verbose) fprintf(stderr, fmt); \
  if (debugfile) fprintf(debugfile, fmt); }
#define DEBUG_OUT1(verbose, debugfile, fmt, arg1) { \
  if (verbose) fprintf(stderr, fmt, arg1); \
  if (debugfile) fprintf(debugfile, fmt, arg1); }
#define DEBUG_OUT2(verbose, debugfile, fmt, arg1, arg2) { \
  if (verbose) fprintf(stderr, fmt, arg1, arg2); \
  if (debugfile) fprintf(debugfile, fmt, arg1, arg2); }
#define DEBUG_OUT3(verbose, debugfile, fmt, arg1, arg2, arg3) { \
  if (verbose) fprintf(stderr, fmt, arg1, arg2, arg3); \
  if (debugfile) fprintf(debugfile, fmt, arg1, arg2, arg3); }
#define DEBUG_OUT4(verbose, debugfile, fmt, arg1, arg2, arg3, arg4) { \
  if (verbose) fprintf(stderr, fmt, arg1, arg2, arg3, arg4); \
  if (debugfile) fprintf(debugfile, fmt, arg1, arg2, arg3, arg4); }
#define DEBUG_OUT5(verbose, debugfile, fmt, arg1, arg2, arg3, arg4, arg5) { \
  if (verbose) fprintf(stderr, fmt, arg1, arg2, arg3, arg4, arg5); \
  if (debugfile) fprintf(debugfile, fmt, arg1, arg2, arg3, arg4, arg5); }
#else
// easier in g++:
#define DEBUG_OUT(verbose, debugfile, args...) { \
  if (verbose) fprintf(stderr, args); \
  if (debugfile) fprintf(debugfile, args); }
#define DEBUG_OUT0 DEBUG_OUT
#define DEBUG_OUT1 DEBUG_OUT
#define DEBUG_OUT2 DEBUG_OUT
#define DEBUG_OUT3 DEBUG_OUT
#define DEBUG_OUT4 DEBUG_OUT
#define DEBUG_OUT5 DEBUG_OUT
#endif


double getLogSumExp(std::vector<double> v);
// prototypes


void setGlobalNpixels(std::vector<int>& w, std::vector<int>& h);
void setGlobalNpixels(std::vector<int>& d, std::vector<int>& w, std::vector<int>& h);
void setGlobalWidth(std::vector<int>& v);

// make sure MRF library is compiled with float costs
void checkForFloatCosts();

//klr related functions
//void copyToMatrix(matrixB<DBL_TYPE> &Yt, matrixB<DBL_TYPE> &Xt, vector <vector <CImage2> > indirImage, vector <vector <CByteImage> > gtdirImage, int testdirindex, int nD, int nbins);


// compute global data cost array - matching cost volume (absolute color diffs)
// also compute best data cost index per pixel (i.e. WTA disparities)
//void initializeDataCostVector(vector <vector<CImage2> >, vector <vector<CByteImage> >, int, int, vector <vector<CByteImage> > &, fvec, fvec, fvec, fvec, vector <vector <vector <CByteImage> > >, Appearance**, vector <vector <vector <CImage2> > >);
void initializeDataCostVector(int genparam, std::vector <std::vector<CImage2> >, std::vector <std::vector<CByteImage> >, int, int, std::vector <std::vector<CByteImage> > &, fvec, fvec, fvec, fvec, int featureCode, int featureCodeKlr,  std::vector <std::vector <std::vector <CByteImage> > >, Appearance**, LocationMV**, AppearanceMV**, std::vector <std::vector <std::vector  <matrixB<double> > > > appdirMVProb, intClass**,  std::vector <std::vector<CByteImage> >, std::vector <std::vector <std::vector <CImage2> > >, std::vector <std::vector <std::vector  <matrixB<double> > > >, intensityNB* intNBclass,
		matrix<DBL_TYPE> wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM, DBL_TYPE lambda, DBL_TYPE p1,
        std::vector<int> startSliceNo, std::vector<int> endSliceNo, int nbinsNB, int hogdimklr, int appdimklr, int locdimklr, int biasklr);

void initializeDataCostVector(int genparam, std::vector <std::vector<CImage2> >, std::vector <std::vector<CByteImage> >, std::vector <std::vector<CByteImage> >, int, std::vector <std::vector<CByteImage> > &, std::vector <std::vector <std::vector <CByteImage> > >, Appearance**, LocationMV**, AppearanceMV**, std::vector <std::vector <std::vector  <matrixB<double> > > > appdirMVProb, intClass**,  std::vector <std::vector<CByteImage> >, std::vector <std::vector <std::vector <CImage2> > >, std::vector <std::vector <std::vector  <matrixB<double> > > >, intensityNB* intNBclass,
		matrix<DBL_TYPE> wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM, std::vector<int> startSliceNo, std::vector<int> endSliceNo, Parameters prm);

void initializeDataCostVector(int genparam, std::vector <std::vector<CByteImage> >, std::vector <std::vector<CByteImage> >, std::vector <std::vector<CByteImage> >, int, std::vector <std::vector<CByteImage> > &, std::vector <std::vector <std::vector <CByteImage> > >, Appearance**, LocationMV**, AppearanceMV**, std::vector <std::vector <std::vector  <matrixB<double> > > > appdirMVProb, intClass**,  std::vector <std::vector<CByteImage> >, std::vector <std::vector <std::vector <CImage2> > >, std::vector <std::vector <std::vector  <matrixB<double> > > >, intensityNB* intNBclass,
		matrix<DBL_TYPE> wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM, std::vector<int> startSliceNo, std::vector<int> endSliceNo, Parameters prm);



void initializeDataCostVectorPar(int genparam, std::vector <std::vector<CImage2> >, std::vector <std::vector<CByteImage> >, int, int, std::vector <std::vector<CByteImage> > &, fvec, fvec, fvec, fvec, int featureCode, int featureCodeKlr,  std::vector <std::vector <std::vector <CByteImage> > >, Appearance**, LocationMV**, AppearanceMV**, std::vector <std::vector <std::vector  <matrixB<double> > > > appdirMVProb, intClass**,  std::vector <std::vector<CByteImage> >, std::vector <std::vector <std::vector <CImage2> > >, std::vector <std::vector <std::vector  <matrixB<double> > > >, intensityNB* intNBclass,
		matrix<DBL_TYPE> wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM, DBL_TYPE lambda, DBL_TYPE p1,
        std::vector<int> startSliceNo, std::vector<int> endSliceNo, int nbinsNB, int hogdimklr, int appdimklr, int locdimklr, int biasklr);


/*
void computeAppearanceVectors(std::vector <std::vector <CImage2> >, int, fvec, Appearance**, std::vector <std::vector <std::vector <CByteImage> > > &, int);
void computeAppearanceVectors(std::vector <std::vector <CByteImage> > im1, int nD, fvec thetaA, Appearance** appclass,
                              std::vector <std::vector <std::vector <CByteImage> > > &appdirImage, int app);

void computeAppearanceProb(std::vector <std::vector <CImage2> > im1,
						int nD,
						std::vector <std::vector <std::vector <matrixB<double> > > > &appdirMVProb,
						AppearanceMV** AppMVclass
						);

void readHoGProb(std::vector <std::vector <CImage2> > im1, int nD, std::vector <std::vector <std::vector <matrixB<double> > > > &hogdirMVProb, std::vector<int> testDirIndexV);
*/

void copyThetaVZC(fvec thetaV, fvec thetaZ, fvec thetaC, int nD);



/*
// create data cost from global cost array
DataCost *computeDataCost(int index);
MRF::CostVal *computeDataCostArray(int index);
*/

/*
// create data cost using parameters thetaU
DataCost *computeDataCost(CByteImage set1Im, std::vector<float> thetaU, int index);
MRF::CostVal *computeDataCostArray(CByteImage set1Im, std::vector<float> thetaU, int index);
*/

/*
DataCost *computeDataCost(std::vector<CImage2> set1Im, std::vector <CByteImage> hogIm, std::vector<float> thetaU, std::vector<float> thetaA, std::vector<float> thetaH, std::vector<float> thetaL, std::vector <std::vector <CByteImage> > appImage, Appearance** appclass, int index, int loc);
MRF::CostVal *computeDataCostArray(std::vector<CImage2> set1Im, std::vector <CByteImage> hogIm, std::vector<float> thetaU, std::vector<float> thetaA, std::vector<float> thetaH, std::vector<float> thetaL, std::vector <std::vector <CByteImage> > appImage, Appearance** appclass, int index, int loc);
*/

DataCost *computeDataCost(int generative, int index);
MRF::CostVal* computeDataCostArray(int generative, int index);


DataCost *computeDataCost(int generative, std::vector<CImage2> set1Im, std::vector <CByteImage> hogIm, int numTrainingPats, std::vector<float> thetaU, std::vector<float> thetaA, std::vector<float> thetaH, std::vector<float> thetaL, int featureCode, int featureCodeKlr, std::vector <std::vector <CByteImage> > appImage, Appearance** appclass, int index, LocationMV**, AppearanceMV**, std::vector <std::vector  <matrixB<double> > > , intClass**, std::vector <std::vector <CImage2> > locdirImage, std::vector <std::vector  <matrixB<double> > >, intensityNB* intNBclass,
		matrix<DBL_TYPE> wparam,
        matrix<DBL_TYPE> Xtrain,
        matrix<DBL_TYPE> minmaxM, DBL_TYPE lambda, DBL_TYPE p1,
        std::vector<int> startSliceNo, std::vector<int> endSliceNo,
        int nbinsNB, int hogdimklr, int appdimklr, int locdimklr, int biasklr);
MRF::CostVal *computeDataCostArray(int generative, std::vector<CImage2> set1Im, std::vector <CByteImage> hogIm, int numTrainingPats, std::vector<float> thetaU, std::vector<float> thetaA, std::vector<float> thetaH, std::vector<float> thetaL, int featureCode, int featureCodeKlr, std::vector <std::vector <CByteImage> > appImage, Appearance** appclass, int index, LocationMV**, AppearanceMV**, std::vector <std::vector  <matrixB<double> > >, intClass**, std::vector <std::vector <CImage2> > locdirImage, std::vector <std::vector  <matrixB<double> > >, intensityNB* intNBclass,
		matrix<DBL_TYPE> wparam,
        matrix<DBL_TYPE> Xtrain,
        matrix<DBL_TYPE> minmaxM, DBL_TYPE lambda, DBL_TYPE p1,
        std::vector<int> startSliceNo, std::vector<int> endSliceNo,
        int nbinsNB, int hogdimklr, int appdimklr, int locdimklr, int biasklr);

DataCost *computeDataCost(int generative, std::vector<CImage2> set1Im, std::vector <CByteImage> hogIm, int numTrainingPats, std::vector <std::vector <CByteImage> > appImage, Appearance** appclass, int index, LocationMV**, AppearanceMV**, std::vector <std::vector  <matrixB<double> > > , intClass**, std::vector <std::vector <CImage2> > locdirImage, std::vector <std::vector  <matrixB<double> > >, intensityNB* intNBclass,
		matrix<DBL_TYPE> wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM,
        std::vector<int> startSliceNo, std::vector<int> endSliceNo, Parameters prm);

MRF::CostVal *computeDataCostArray(int generative, std::vector<CImage2> set1Im, std::vector <CByteImage> hogIm, int numTrainingPats, std::vector <std::vector <CByteImage> > appImage, Appearance** appclass, int index, LocationMV**, AppearanceMV**, std::vector <std::vector  <matrixB<double> > >, intClass**, std::vector <std::vector <CImage2> > locdirImage, std::vector <std::vector  <matrixB<double> > >, intensityNB* intNBclass,
		matrix<DBL_TYPE> wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM,
        std::vector<int> startSliceNo, std::vector<int> endSliceNo, Parameters prm);

DataCost *computeDataCost(int generative, std::vector<CByteImage> set1Im, std::vector <CByteImage> hogIm, std::vector <CByteImage> motIm, int numTrainingPats, std::vector <std::vector <CByteImage> > appImage, Appearance** appclass, int index, LocationMV**, AppearanceMV**, std::vector <std::vector  <matrixB<double> > > , intClass**, std::vector <std::vector <CImage2> > locdirImage, std::vector <std::vector  <matrixB<double> > >, intensityNB* intNBclass,
		matrix<DBL_TYPE> wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM,
        std::vector<int> startSliceNo, std::vector<int> endSliceNo, Parameters prm);

MRF::CostVal *computeDataCostArray(int generative, std::vector<CByteImage> set1Im, std::vector <CByteImage> hogIm, std::vector <CByteImage> motIm, int numTrainingPats, std::vector <std::vector <CByteImage> > appImage, Appearance** appclass, int index, LocationMV**, AppearanceMV**, std::vector <std::vector  <matrixB<double> > >, intClass**, std::vector <std::vector <CImage2> > locdirImage, std::vector <std::vector  <matrixB<double> > >, intensityNB* intNBclass,
		matrix<DBL_TYPE> wparam, matrix<DBL_TYPE> Xtrain, matrix<DBL_TYPE> minmaxM,
        std::vector<int> startSliceNo, std::vector<int> endSliceNo, Parameters prm);






DataCost *computeDataCostPar(int generative, std::vector<CImage2> set1Im, std::vector <CByteImage> hogIm, int numTrainingPats, std::vector<float> thetaU, std::vector<float> thetaA, std::vector<float> thetaH, std::vector<float> thetaL, int featureCode, int featureCodeKlr, std::vector <std::vector <CByteImage> > appImage, Appearance** appclass, int index, LocationMV**, AppearanceMV**, std::vector <std::vector  <matrixB<double> > > , intClass**, std::vector <std::vector <CImage2> > locdirImage, std::vector <std::vector  <matrixB<double> > >, intensityNB* intNBclass,
		matrix<DBL_TYPE> wparam,
        matrix<DBL_TYPE> Xtrain,
        matrix<DBL_TYPE> minmaxM, DBL_TYPE lambda, DBL_TYPE p1,
        std::vector<int> startSliceNo, std::vector<int> endSliceNo,
        int nbinsNB, int hogdimklr, int appdimklr, int locdimklr, int biasklr);
MRF::CostVal *computeDataCostArrayPar(int generative, std::vector<CImage2> set1Im, std::vector <CByteImage> hogIm, int numTrainingPats, std::vector<float> thetaU, std::vector<float> thetaA, std::vector<float> thetaH, std::vector<float> thetaL, int featureCode, int featureCodeKlr, std::vector <std::vector <CByteImage> > appImage, Appearance** appclass, int index, LocationMV**, AppearanceMV**, std::vector <std::vector  <matrixB<double> > >, intClass**, std::vector <std::vector <CImage2> > locdirImage, std::vector <std::vector  <matrixB<double> > >, intensityNB* intNBclass,
		matrix<DBL_TYPE> wparam,
        matrix<DBL_TYPE> Xtrain,
        matrix<DBL_TYPE> minmaxM, DBL_TYPE lambda, DBL_TYPE p1,
        std::vector<int> startSliceNo, std::vector<int> endSliceNo,
        int nbinsNB, int hogdimklr, int appdimklr, int locdimklr, int biasklr);





// get pointer to cost array and its dimensions
void getCostArray(float* &costs, int &nPixels, int &nDisps);


//logistic regression dist functions


// compute quantized absolute color gradients of input image
void computeQuantizedGradients(CByteImage im1, CByteImage &im1Grad, std::vector<int> gradThresh);
/*
  use nK = gradThresh.size()+1 levels of quantization.
  quantization is done as follows:
  compute absolute gradient (Euclidean distance in color space) in x and y
  given gradient g, assign largest k such that g <= gradThresh[k]
  entry gradThresh[nK] is assumed to be infinity
  in im1Grad, use band 0 for x gradient, and band 1 for y gradients
*/

void computeQuantizedGradients(std::vector<CImage2> im1, std::vector<CByteImage> &im1Grad, std::vector<int> gradThresh);
void computeQuantizedGradients(std::vector<CByteImage> im1, std::vector<CByteImage> &im1Grad, std::vector<int> gradThresh);
void computeQuantizedGradients(std::vector<CByteImage> im1, std::vector<CByteImage> &im1Grad, std::vector<int> gradThresh, std::vector<int> gradThreshZ);

void computeQuantizedGradients00(std::vector<CImage2> im1, std::vector<CByteImage> &im1grad, std::vector<int> gradThresh);
void computeQuantizedGradients0(std::vector<CImage2> im1, std::vector<CByteImage> &im1grad, std::vector<int> gradThresh);

void computeQuantizedGradients016(std::vector<CImage2> im1, std::vector<CImage2> &im1grad);

SmoothnessCost *computeSmoothnessCost(CByteImage im1grad, std::vector<float> theta, int i);
MRF::SmoothCostGeneralFn computeSmoothnessCostFunction(CByteImage im1grad, std::vector<float> theta, int i);
// I don't know how to use a default object argument, so I leave the above two functions in  - JJW
SmoothnessCost *computeSmoothnessCost(CByteImage im1grad, std::vector<float> theta, std::vector<float> thetaP, int i);
MRF::SmoothCostGeneralFn computeSmoothnessCostFunction(CByteImage im1grad, std::vector<float> theta, std::vector<float> thetaP, int i);


SmoothnessCost *computeSmoothnessCost(std::vector<CByteImage> im1grad, std::vector<float> theta, int i);
MRF::SmoothCostGeneralFn computeSmoothnessCostFunction(std::vector<CByteImage> im1grad, std::vector<float> theta, int i);

SmoothnessCost *computeSmoothnessCost(std::vector<CByteImage> im1grad, std::vector<float> thetaV, std::vector<float> thetaC, int i, int gradContext);
MRF::SmoothCostGeneralFn computeSmoothnessCostFunction(std::vector<CByteImage> im1grad, std::vector<float> thetaV, std::vector<float> thetaC, int i, int gradContext);

SmoothnessCost *computeSmoothnessCost00(std::vector<CByteImage> im1grad, std::vector<float> thetaV, std::vector<float> thetaC, int i, int gradContext);
MRF::SmoothCostGeneralFn computeSmoothnessCostFunction00(std::vector<CByteImage> im1grad, std::vector<float> thetaV, std::vector<float> thetaC, int i, int gradContext);


SmoothnessCost *computeSmoothnessCost(std::vector<CByteImage> im1grad, std::vector<float> thetaV, std::vector<float> thetaZ, std::vector<float> thetaC, int i, int gradContext);
MRF::SmoothCostGeneralFn computeSmoothnessCostFunction(std::vector<CByteImage> im1grad, std::vector<float> thetaV, std::vector<float> thetaZ, std::vector<float> thetaC, int i, int gradContext);


/*
the parameters that determine the smoothness cost are stored in a
std::vector so they can be handled uniformly during learning:

  std::vector<float> theta;

  int nK = theta.size();
  //theta[0]: cost for |dp-dq| == 1
  //theta[k], k=1..nK-1: cost for |dp-dq| > 1 and gradient_pq < gradThresh[k]

OOOPS, I just realized that the MRF library doesn't support this, except if you set
up the smoothness cost as a general function, which is much slower.  If you
use the hCue and vCue parameters, then the smoothness cost (which only depends
on the two labels) simply gets multiplied by constants...

I can modify the graph cut code later, but for now let's stick with the old model:
  //theta[k], k=0..nK-1: cost for |dp-dq| > 0 and gradient_pq < gradThresh[k]
*/

// delete any arrays allocated for data and smoothness costs
void deleteGlobalArrays();


// main stereo matcher
// if WTAdisp is a valid image, it is used to initialize the labels
void crfmodel(int width, int height, int nD, DataCost *dcost, SmoothnessCost *scost,
	       CByteImage WTAdisp, CByteImage &disp, float closeEnoughPercent);

void crfmodel(int width, int height, int nD, DataCost *dcost, SmoothnessCost *scost,
               CByteImage WTAdisp, CByteImage &disp, float closeEnoughPercent,
               int inferencer, int innerIter = -1);

void crfmodel(int width, int height, CByteImage &disp, MRF *mrf,
               float closeEnoughPercent, int innerIter = 1,
               int maxIter = 50);

/*
MRF* crfmodel(int depth, int width, int height, int nD, DataCost *dcost, SmoothnessCost *scost,
	       std::vector<CByteImage> WTAdisp, std::vector<CByteImage> &disp, float closeEnoughPercent);
*/
MRF* crfmodel(int depth, int width, int height, int nD, DataCost *dcost, SmoothnessCost *scost,
               std::vector<CByteImage> &disp, float closeEnoughPercent,
               int inferencer, int innerIter, std::vector<CByteImage> interImage, std::vector<CByteImage> truedisp, int interactive, std::ofstream *logStream, int testv);

void crfmodel(int depth, int width, int height, std::vector<CByteImage> &disp, MRF *mrf,
               float closeEnoughPercent, int innerIter = 1,
               int maxIter = 50);





void writeDisparities(CByteImage disp, int outscale, char *outstem);
void writeDisparities(std::vector<CByteImage> disp, int outscale, char *outstem);

void setPairwise(bool v);
void setPairwiseInteraction(bool v);
void setLocal(bool v);
void updateThetaP(std::vector<float> theta);
int getNumStates();
float rangeNormailize(float rawMin, float rawMax, float min, float max, float real);

void updateCueIndex(unsigned int i);
bool getLocal();
void setRandom(bool rv);

double estimateLogLikelihood(std::vector<CByteImage> WTA, std::vector<float> thetaU, std::vector<float> thetaV, int index, int nD );


void logisticRegressionPairwise(std::vector <CImage2> im1, int nD, int numTrainingPats, std::vector <CByteImage> &wta, fvec thetaV, std::vector<int> gradThresh,
                                std::vector <CByteImage> im1grad, std::vector<int> startSliceNo, std::vector<int> endSliceNo, int index, double &loglikelihood, std::vector <CByteImage> gtIm);

void logisticRegressionPairwise(int nD, int numTrainingPats, std::vector <CByteImage> &wta, fvec thetaV, std::vector<int> gradThresh, std::vector<CByteImage> im1grad, std::vector<int> startSliceNo, std::vector<int> endSliceNo, int index, double &loglikelihood, std::vector <CByteImage> gtIm);
void logisticRegressionPairwise(int nD, int numTrainingPats, std::vector <CByteImage> &wta, fvec thetaV, std::vector<int> gradThresh , fvec thetaZ, std::vector<int> gradThreshZ , std::vector<CByteImage> im1grad, std::vector<int> startSliceNo, std::vector<int> endSliceNo, int index, double &loglikelihood, std::vector <CByteImage> gtIm);

void logisticRegressionPairwise(int nD, int numTrainingPats, std::vector <CByteImage> &wta, fvec thetaV, std::vector<int> gradThresh, fvec thetaZ, std::vector<int> gradThreshZ, std::vector<CByteImage> im1grad, std::vector<int> startSliceNo, std::vector<int> endSliceNo, int index, double &loglikelihood, std::vector <CByteImage> gtIm, std::vector <std::vector <CByteImage> > &dirDispXYZ);


void computeEmpirDistVlpd4(std::vector<CByteImage> validdisp, std::vector<CByteImage> disp, std::vector<CByteImage> im1grad, fvec &distV, int nD, int ignoreVal, int gradContext, int logregpaird4, std::vector<int> gradThreshVec, int empir, int index, int paramlreg); 

void logisticRegressionPairwiseD(int genparam, std::vector <CImage2> im1, std::vector <CByteImage>  hogIm,
      int nD, int numTrainingPats, std::vector <CByteImage> &wta,      // winner-take-all disparities
      fvec thetaU, fvec thetaA, fvec thetaH, fvec thetaL, int featureCode, std::vector <std::vector <CByteImage> > appdirImage,
      Appearance** appclass, LocationMV** LocationMVclass, AppearanceMV** AppMVclass,
      std::vector <std::vector  <matrixB<double> > >  appdirMVProb, intClass **intObj, std::vector <std::vector <CImage2> > locdirImage,
      std::vector <std::vector  <matrixB<double> > > hogdirMVProb, intensityNB* intNBclass, 
      std::vector<int> startSliceNo, std::vector<int> endSliceNo, int nbinsNB, int index, double &loglikelihood, 
                                std::vector <CByteImage> gtIm, fvec thetaV, std::vector<int> gradThresh, std::vector<CByteImage> im1grad);


void logisticRegressionPairwiseD4(int genparam, std::vector <CImage2> im1, std::vector <CByteImage>  hogIm,
      int nD, int numTrainingPats, std::vector <CByteImage> &wta,      // winner-take-all disparities
      fvec thetaU, fvec thetaA, fvec thetaH, fvec thetaL, int featureCode, std::vector <std::vector <CByteImage> > appdirImage,
      Appearance** appclass, LocationMV** LocationMVclass, AppearanceMV** AppMVclass,
      std::vector <std::vector  <matrixB<double> > >  appdirMVProb, intClass **intObj, std::vector <std::vector <CImage2> > locdirImage,
      std::vector <std::vector  <matrixB<double> > > hogdirMVProb, intensityNB* intNBclass, 
      std::vector<int> startSliceNo, std::vector<int> endSliceNo, int nbinsNB, int index, double &loglikelihood, 
                                  std::vector <CByteImage> gtIm, fvec thetaV, std::vector<int> gradThresh, std::vector<CByteImage> im1grad);


void reInitializeDataCostVectorPairwise(int numPatients, int genparam, int nD, int power);

void setGlobalValues(int nD, int nG, int nGz);

void copyImSet(std::vector <std::vector <CByteImage> > &dirDisp, std::vector <std::vector <CByteImage> > &dirWTA);


// void computeInterDataCost(DataCost *dcost, SmoothnessCost *scost, std::vector<CByteImage> interImage, std::vector<CByteImage> truedisp);

void checkReadImage(std::vector <std::vector <CByteImage> > interdirImage);

#endif
