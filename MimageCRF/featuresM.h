/*
 * featuresM.h
 *
 *  Created on: Mar 16, 2011
 *      Author: bhole
 */

#ifndef FEATURESM_H_
#define FEATURESM_H_



#include "locationMV.h"
#include "appearance.h"
#include "appearanceMV.h"
#include "intensityClass.h"
#include "intensityNB.h"

#include "Parameters.h"

// really need to chamge this
#include "features.h"

#include "Global.h"

//#include "opticalflow/project.h"
#include "opticalflow/Image.h"
#include "opticalflow/OpticalFlow.h"

class features
{

	int readIntensityFiles(int nD);
	int readLocationFiles(int nD, unsigned int testDir, int numP, int timeseq, std::vector<std::string> dirListNFP);
	int readHOGFiles(int nD, unsigned int testDir, int numP);
	int readMOTFiles(int nD, unsigned int testDir, int numP);
	int readAppearanceFiles(int timeseq, int nD);

	void read3DlocationSliceNumbers(int numP, std::vector<std::string> dirListNFP);
	void readLocationMVParameters(int nD);

	void computeLocation(std::vector <std::vector <CImage2> > *indirImage, int nD);
	int deAllocateLocationFiles(int timeseq);


	int computeIntensityNB(std::vector <std::vector <CImage2> > *inp, std::vector <std::vector <CByteImage> > *gtr, std::vector<int> testDirIndexV, int nD);
	int computeIntensityNB(std::vector <std::vector <CByteImage> > *inp, std::vector <std::vector <CByteImage> > *gtr, std::vector<int> testDirIndexV, int nD);


	void computeAppearanceVectors(std::vector <std::vector <CImage2> > * im1, int nD);
	void computeAppearanceVectors(std::vector <std::vector <CByteImage> > * im1, int nD);
	void computeAppearanceProb(std::vector <std::vector <CImage2> > * im1, int nD);
	int computeAppearance(std::vector <std::vector <CImage2> > *inp, std::vector <std::vector <CByteImage> > *gtr, int nD);
	int computeAppearance(std::vector <std::vector <CByteImage> > *inp, std::vector <std::vector <CByteImage> > *gtr, int nD);


	void computeQuantizedGradients(std::vector <std::vector <CImage2> > *inp);
	void computeQuantizedGradients(std::vector <std::vector <CByteImage> > *inp);
	void computeQuantizedGradientsOpflow(std::vector <std::vector <CByteImage> > *inp);

	void computeCSGradients(std::vector <std::vector <CImage2> > *inp);
	void computeCSGradients(std::vector <std::vector <CByteImage> > *inp);
	void computeCSGradientsOpflow(std::vector <std::vector <CByteImage> > *inp);


	int preComputeOpticalFlow(std::vector <string> &indirpath);
	void checkConsistencyCacheFiles(std::vector<string> &indirpath);
	void readOpflowCache();
	void writeOpflowCache();

	int readKlrFiles(int nD);
	void printklrParameters();




protected:

	featureparams* fpvp;

	LocationMV** LocationMVclass;
	Appearance** appclass;
	AppearanceMV** AppMVclass;
	intClass **intObj;
	intensityNB* intNBclass;


  std::vector<string> locdirListFP;
	std::vector <std::vector <CImage2> > locdirImage;
	std::vector <std::vector <IplImage*> > locdirImageClust;
	std::vector<int> startSliceNo, endSliceNo;


  std::vector<string> hogdirListFP;
  std::vector <std::vector <CImage2> > hogdirImage;
  std::vector <std::vector <std::vector  <matrixB<double> > > > hogdirMVProb;


  std::vector<string> motdirListFP;
  std::vector <std::vector <CImage2> > motdirImage;
  std::vector <std::vector <std::vector  <matrixB<double> > > > motdirMVProb;

  // last dimension is [0] for u or vx and [1] for v for vy
  std::vector<std::vector <std::vector <opflowns::DImage> > > flowvector;
  // last dimension is [0] for u or vx and [1] for v for vy
  std::vector<std::vector <std::vector <opflowns::DImage> > > reverseflowvector;


  std::vector <std::vector <std::vector <CByteImage> > > appdirImage;
  std::vector <std::vector <std::vector  <matrixB<double> > > > appdirMVProb;

  matrix<DBL_TYPE> Xtrain;
  matrix<DBL_TYPE> wparam;
  matrix<DBL_TYPE> minmaxM;


	std::vector <std::vector <CByteImage> > dirgrad;
	std::vector <std::vector <IplImage*> > dircsgrad;


	IVM* ivm;
	std::vector<string> hogklrdescdirListFP;
	std::vector<string> appklrdescdirListFP;


	double getCostLocationDiscGaussian3(int pno, unsigned int z, int x, int y, int d, int genparam, int nD, int &thetaid);


public:

	features(featureparams* featurepvp);
	int readFeatureFilesDirs(int timeseq, int nD, int testDir, int numP, std::vector<std::string> dirListNFP);
	int deleteFeatureAllocations(int nD, int timeseq);
	int deAllocatePairwiseAllocations();

	int preComputeLocalFeatures(std::vector <std::vector <CImage2> > *inp, std::vector <std::vector <CByteImage> > *gtr, std::vector<int> testDirIndexV, int nD);
	int preComputeLocalFeatures(std::vector <std::vector <CByteImage> > *inp, std::vector <std::vector <CByteImage> > *gtr, std::vector<int> testDirIndexV, int nD);

	virtual int preComputePairwiseFeatures(std::vector <std::vector <CImage2> > *inp);
	virtual int preComputePairwiseFeatures(std::vector <std::vector <CByteImage> > *inp);
	virtual int preComputePairwiseFeaturesAndPutInGlobalSpace(std::vector <std::vector <CByteImage> > *inp);

	void extractOpticalFlow(std::vector<string> &indirpath);



	void deleteLocationMVclass(int nD);
	void deleteAppMVclass(int nD);
	void deleteappclass(int nD);
	void deleteintObj(int nD);


	virtual ~features();

};



#endif /* FEATURESM_H_ */
