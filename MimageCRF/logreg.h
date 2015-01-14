/*
 * logreg.h
 *
 *  Created on: Mar 15, 2011
 *      Author: bhole
 */

#ifndef LOGREG_H_
#define LOGREG_H_

#include "Parameters.h"
#include "baseModel.h"


class logreg: public baseModel {

	void logisticRegression(int genparam, unsigned int index);

public:

	logreg(Parameters* prm, std::vector<float> theta);

	void allocateDataCostSpace();
	void deAllocateDataCostSpace();

	void allocateSmoothCostGlobalSpace();
	void deAllocateSmoothCostGlobalSpace();


	int preComputePairwiseFeatures(std::vector <std::vector <CImage2> > *inp);
	int preComputePairwiseFeatures(std::vector <std::vector <CByteImage> > *inp);

	virtual int preComputePairwiseFeaturesAndPutInGlobalSpace(std::vector <std::vector <CByteImage> > *inp) {}


	void modelProcess(unsigned int index);

	void updateDist(int index, int fiterflag);

	void computeEmpirDist(int index);
	void computeModelDist(int index, int z);
	void computeModelDist(int index, int z, matrix<DBL_TYPE> &hogdescs, matrix<DBL_TYPE> &appdescs);

	void computeEmpirDistU(int index);
	void computeEmpirDistA(int index);
	void computeEmpirDistH(int index);
	void computeEmpirDistM(int index);
	void computeEmpirDistL(int index);
	void computeEmpirDistW(int index);
	void computeEmpirDistO(int index);
	void computeEmpirDistB(int index);
	void computeEmpirDistE(int index);

	void computeModelDistU(int index, int z);
	void computeModelDistA(int index, int z);
	void computeModelDistH(int index, int z);
	void computeModelDistM(int index, int z);
	void computeModelDistL(int index, int z);
	void computeModelDistW(int index, int z, matrix<DBL_TYPE> &hogdescs, matrix<DBL_TYPE> &appdescs);
	void computeModelDistO(int index, int z);
	void computeModelDistB(int index, int z);
	void computeModelDistE(int index, int z);

};



#endif /* LOGREG_H_ */
