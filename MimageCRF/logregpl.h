/*
 * logregpl.h
 *
 *  Created on: Mar 15, 2011
 *      Author: bhole
 */

#ifndef LOGREGPL_H_
#define LOGREGPL_H_

#include "Parameters.h"
#include "baseModel.h"

extern GlobalPairwiseParam globalP;


class logregpl: public baseModel {

    MRF* mrf;
    MRF* mrfgen;

//	void lRpl(int genparam, unsigned int index);
    // this is used for training
    double getSmoothCost(unsigned int index, unsigned int z, int y, int x, int d);

    // these are used for testing
    DataCost* computeDataCost(int genparam, unsigned int index);
	SmoothnessCost* computeSmoothnessCost(int index);

public:

	logregpl(Parameters* prm, std::vector<float> theta);
	~logregpl();

	void deletemrf(unsigned int index);

	void allocateDataCostSpace();
	void deAllocateDataCostSpace();

	void allocateSmoothCostGlobalSpace();
	void deAllocateSmoothCostGlobalSpace();


	void modelProcess(unsigned int index);
	void modelProcessTrain(int genparam, unsigned int index);
	void modelProcessTest(int genparam, unsigned int index);

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
	void computeEmpirDistV(int index);


	void computeModelDistU(int index, int z);
	void computeModelDistA(int index, int z);
	void computeModelDistH(int index, int z);
	void computeModelDistM(int index, int z);
	void computeModelDistL(int index, int z);
	void computeModelDistW(int index, int z, matrix<DBL_TYPE> &hogdescs, matrix<DBL_TYPE> &appdescs);
	void computeModelDistV(int index, int z);

};



#endif /* LOGREGPL_H_ */
