/*
 * baseModel.h
 *
 *  Created on: Mar 15, 2011
 *      Author: bhole
 */

#ifndef BASEMODEL_H_
#define BASEMODEL_H_

#include "featuresM.h"
#include "IO.h"
#include "evaluation.h"
#include "bfgs.h"
#include "Parameters.h"
#include "common.h"



class baseModel: public features, public IO, public bfgs, public evaluation
{

private:
	void printfVecToLOG(fvec a, char *name);
	void sgdupdate(int &fullrun, int &startidx, int iter);
	int skipvolume(int i, int startidx, int fullrun);
	double getrmsdiffDist(fvec a);
	void weightClassValueTemp(double &dsiValue, int d); // hack - get rid of it in future
	void weightClassDistTemp(fvec &Dist);

	void weightClassValue(double &dsiValue, int d); // hack - get rid of it in future
	void weightClassDist(fvec &Dist);


protected:
  fvec totalEmpirDistU, totalModelDistU;
  fvec totalEmpirDistA, totalModelDistA;
  fvec totalEmpirDistH, totalModelDistH;
  fvec totalEmpirDistM, totalModelDistM;
  fvec totalEmpirDistL, totalModelDistL;
  fvec totalEmpirDistV, totalModelDistV;
  fvec totalEmpirDistZ, totalModelDistZ;
  fvec totalEmpirDistC, totalModelDistC;
  fvec totalEmpirDistW, totalModelDistW;
  fvec totalEmpirDistO, totalModelDistO;
  fvec totalEmpirDistB, totalModelDistB;
  fvec totalEmpirDistE, totalModelDistE;

  fvec classSize;
  fvec inverseClassSize;

  fvec oldTheta;
  fvec regTheta;

	Parameters *prm;

	std::vector<MRF::CostVal *> dCArray;  // change this better to have it as a two 2 array with nD as other dimension
	std::vector<MRF::CostVal *> dCArrayGen;

	double oldloglikelihood;
	double oldreldiffnorm;
	double oldrmsdiffDist;
	double loglikelihood;
	double loglikelihoodgen;

	// for deciding theta indices in pairwise terms
	void reorder(int *d1, int *d2);

	double getbadCost(int genparam);
	double getDataCost(int genparam, int width, int height, int depth, unsigned int j, unsigned int i, int y, int x, int d, int zslice);
	void getDataCostKlr(ublas::vector<DBL_TYPE> &f, ublas::vector<DBL_TYPE> &Xtest, unsigned int j, unsigned int i, int x, int y, int zslice, matrix<DBL_TYPE> &hogdescs, matrix<DBL_TYPE> &appdescs, int width, int height);

	void readklrdescfiles(matrix<DBL_TYPE> &hogdescs, matrix<DBL_TYPE> &appdescs, int index, int num);

	virtual void regularizeLL();
	void regularizationUpdate();

	int frameExistsInFrameList(const int &frame_num);

public:

	baseModel(Parameters* prmout, std::vector<float> theta);

	virtual ~baseModel();

	virtual void allocateDataCostSpace();
	virtual void deAllocateDataCostSpace();

	virtual void allocateSmoothCostGlobalSpace();
	virtual void deAllocateSmoothCostGlobalSpace();

	void checkSanity();

	void firstInference(int genparam, char* fnpart);

	void dumpThetaParameters(int verbose, int debug);

	void searchParameters();

	void initializeDistToZero();
	void initializeModelDistToZero();
	void initializeEmpirDistToZero();

	void initializeVarToZero();

	virtual void modelProcess(unsigned int index);
	virtual void hiddenModelProcess(unsigned int index, int iter, int iterout);

	virtual void evaluateResults(unsigned int index); // per volume computations
	void evaluateResults(); // for global computations

	void outputEvalutions(unsigned int index);
	void outputEvalutions();
	void outputFinalEvalutions();

	void outputConfMatrix(double **confMatrix, int debug);
	void outputDiceCoeff(std::vector<double> dcs, int debug);
	void outputPixelCounts(double* totalCounts, int debug);

	void outputImageResults(char* fnpart, int volno);
	void outputImageResults(unsigned int index, int iter);
	void outputImageResults(unsigned int index, int iter, const string prefix);

	int returnConcatEmpirDistVector(std::vector<float> *empirDist);
	int returnConcatModelDistVector(std::vector<float> *modelDist);

	int returnConcatThetaWWVector(std::vector<float> *theta);
	int splitThetaWWVectorToOriginal(std::vector<float> theta);

	virtual void updateDist(int index, int fiterflag);
	virtual void computeEmpirDist(int index);
	virtual void computeModelDist(int index);
	virtual void computeModelDist(int index, int z);

	void preComputeClassRatios();
	void preComputeInverseClassRatios();

	// virtual int preComputePairwiseFeaturesAndPutInGlobalSpace(std::vector <std::vector <CByteImage> > *inp) {}


/*
	virtual void computeEmpirDistU();
	virtual void computeEmpirDistA();
	virtual void computeEmpirDistH();
	virtual void computeEmpirDistM();
	virtual void computeEmpirDistL();
	virtual void computeEmpirDistV();
	virtual void computeEmpirDistW();

	virtual void computeModelDistU();
	virtual void computeModelDistA();
	virtual void computeModelDistH();
	virtual void computeModelDistM();
	virtual void computeModelDistL();
	virtual void computeModelDistV();
	virtual void computeModelDistW();
*/

	int makeNextThetaStep(int iterout, int iter);
	void stepBackSlowDownProgress(double reldiffnorm);
	void speedUpProgress(fvec gradTheta, double reldiffnorm, double rmsdiffDist);
	double getRelDiffNorm(fvec gradTheta, fvec empirDist);


	virtual void deletemrf(unsigned int index);

};



#endif /* BASEMODEL_H_ */
