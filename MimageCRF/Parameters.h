#ifndef PARAMETERS_FILE
#define PARAMETERS_FILE

#include <stdio.h>
#include <vector>
#include <map>
#include <fstream>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <ctype.h>
#include <climits>
#include <string.h>
#include <stdlib.h>
#include <string>
#include "common.h"

#include <libconfig.h++>
using namespace libconfig;


typedef std::vector<float> fvec;

// only for g++:
#define LOGBUG(verbose, debugfile, args...) { if (verbose) fprintf(debugfile, args); }
#define LOG(debugfile, args...) { fprintf(debugfile, args); }


//enum kernelType { RBF, CHI2, POLY, LINEAR, EUCLIDEAN, NONE };


class crfparams
{
public:

  int crfp;

  int hidden; // this is 1 if it is a hidden-crf. Note that unknown should be 1 as we need
  // some hidden nodes.

  // inference algorithm parameters and details
  // inference parameters (for crf)

  // USE_GC 0,  USE_BP 1,  USE_SBP 2,  USE_MF 3,  USE_SMF 4, USE_TRWS 5,  USE_BPS 6,  USE_MP 7

  int inferencer;
  int inferencerTest; 
  int random;			   // use to enable random in graphcut

  float damper;
  float msgTol;
//  MeanField::DiffType msgDiffFun = MeanField::ENERGY;

  int inferouteriter;
  int inferinneriter;

  crfparams();

};


class logregparams
{
public:

  int logreg;
  int logregpair;
  int logregpaird;
  int logregpaird4;
  int logregpl;

  logregparams();

};


class ioparams
{
public:

  // this will contain the main directory of input data
  std::string indirname;
  // this will contain the main directory of groundtruth data
  std::string gtdirname;
  // this will contain the main directory of interactive map
  std::string interdirname;
  // this will contain the main directory of output data
  std::string outdirname;
  // this name will be used to create two text files (same name as output directory above)
  // for dumping statistics of learning and accuracy results.
  std::string outstem;



  // to indicate test case
  unsigned int testDir;  // to register that test case exists
  //int testDirIndex = -1; // this will index the directory of test directory in our list of directories.
  // the test directory name has to be "test"
  std::vector<int> testDirIndexV;

  ioparams();

};

class klrparams
{
public:

  // klr-exp 0, ivm 1, svm 3
  int klrtype; // this will indicate whether the strong classifier is klr-exp, ivm, svm, etc

  std::string klrfName;

  int appKlr;
  int hogKlr;
  int intensityKlr;
  int locKlr;
  int motKlr;
  int rbmKlr;

  // parameters to control type of model
  int klr;

  // ivm parameters
  int hogdimklr;
  int motdimklr;
  int rbmdimklr;
  int appdimklr;
  int locdimklr;
  int biasklr;
  double lambda;
  std::string xfile;
  std::string wfile;
  std::string mmfile;
  std::string kernel;
  double p1;
  kernelType ker;

  std::map<std::string, kernelType> mapKernelType;

  int nW;

  std::string hogklrdescdirname;
  std::string appklrdescdirname;


  klrparams();

};




class featureparams
{

public:

  // common

  int generative;
  int featureCode; 
  int featureCodeKlr;
  float gaussSigma;      // Scale of gaussian regularizer (zero for none) - crf parameters only (not klr params)
  int bias;

  int nF;				 // Number of scaling factor. (vestigeal)
  fvec factorU;              // Scaling factor for each crf parameter (vestigeal) 

  std::string fileName;  // file for reading parameters

  int learngradonly;        // computeDataCost is not really rerun to save time
  int learncrfonly;             // really want this 


  // smoothness

  int gradOnly;  // to control whether to use non-pairwise or pairwise gradients or only gradients
  int context;
  int gradContext;
  int gradVZ;
  int csgrad;
  std::vector<int> gradThreshVec;  // 1 for only x-y, 2 for x-y-z (same parameters) and 3 for separate x-y and z parameters
  std::vector<int> gradThreshVecZ; // if gradVZ == 3, then we use this too
  fvec thetaV; // smoothness weights associated with gradient bins
  fvec thetaZ; // if gradVZ == 3, then we use this too along with thetaV
  fvec thetaC; // context weights
  int nV;
  int nZ;
  int nC;
  int nT;
  int nTz;
  int nG;
  int nGZ;

  // volume
  int nE;
  fvec thetaE;
  int volume;

  // intensity

  fvec thetaU; // data weights (at most one for now)
  int nU;
  int intensity;
  int nbinsintNB;
  int rangeI;
  int onlyUV; // if 1 , Y (intensity or luminance) is not used

  std::string intfilename;
  std::string intdirname;


  // location

  int loc; // to control whether location should be used in the generative-discriminative setting
  int locCubeSize;


  int loccubex, loccubey, loccubez;
  int skipXYZlocklr;


  int loc_params_movie;  // vestigeal
  fvec thetaL; //location weights
  std::string locdirname;
  int nL;

  std::string locfilename;  // for reading in parameters - this could be a directory too
//  std::string locdirname;   // usually for reading in per pixel data
  std::string imageSliceFile;


  // appearance

  int app;
  std::string appfilename;
  fvec thetaA; // appearance weights
  int nA;
  
//  std::string appfilename;
  std::string appdirname;


  // hog
  int hog;
  fvec thetaH; // HOG weights
  std::string hogdirname;
  int nH;
  std::string hogfilename;


  // motion

  int mot;
  fvec thetaM;
  std::string motdirname;
  int nM;
  std::string motfilename;

  // optical flow
  int opflowconnection;
  int opflowreverseconnection;

  fvec thetaO; // this feature is assuming for one frame
  int nO; // and for one volume only
  int opflowlocal; // this controls the use of these features
  std::vector<int> opflowlocalframes; // we assume only one folder i.e. no test folder allowed.
  int useopflowcache;
  std::string flowcachefile;
  std::string rflowcachefile;


  // bias
  fvec thetaB;
  int nB;

  // rbm
  int rbm;







  // klr

  klrparams klrpv;


  featureparams();

};



class graddescparams
{

public:

  float rate;
  int maxiter;          // maximum number of iterations for gradient descent
  float closeEnoughPercent; // when to stop GC or BP iterations

  // flags for optimization algorithm (choose either gradient descent or bfgs
  // for now, bfgs works only with logreg
  int grad_desc;
  int bfgs_flag; // when 1 call bfgs, when 2 call lbfgs

  int bfgs_outer_iter;
  int bfgs_inner_iter;
  // when using bfgs, set iterations provided by user to bfgs_outer_iter

  // need only maxiter or outer_iter
  // need only grad_desc or bfgs_flag

  int maxiter_out;

  int sgdskipV; // stochastic gradient descent skip volumes
  int sgdlast; // in last x iterations, all volumes will be used



  int m_bfgs; //limited memory version
  double beta;
  double beta_dash;
  double alpha;
  double stepmaxconst;

  graddescparams();

};







class Parameters 
{

public:

  // common

  int nD; // number of class labels
  int outscale8;         // scale factor for disparities; -1 means full range 255.0/(nD-1)
  int outscale16;        // (2^16-1)/(nD-1) 
  int verbose;               // print messages to stderr

  int ignoreVal;        // Value of ground truth disp image to ignore - (vestigeal)

  int parallelize;    // whether you want to run Par versions
  
  int timeseq;

  int interactive;

  int optClassAccuracy;

  int unknown; // This parameter dictates if you have some region in the image that are not to
  // be used for training. The values in images will lie inbetween valid values eg 127 for a
  // rgb image for two classes that have background 0 and foreground 255.

  // crf parameters

  crfparams crfpv;

  // logreg parameters

  logregparams logregpv;

  // gradient descent parameters

  graddescparams graddescpv;


  // feature parameters

  featureparams featurepv;


  // data
  // io files/directories parameters

  ioparams iopv;


  // constructor to put in default values
  Parameters();

  int readfromcfg(Config &cfg, const char *cfgfname);

  int printParameters();
  int readfromCommandline(int argc, char **argv);
  int updateDependParameters();
  int readfeaturefile();
  int updatePropertiesParameters();
  int returnConcatThetaVector(std::vector<float> *theta);
  int splitVectorToOriginal(std::vector<float> theta);

};


#endif
