/*
 * crf.h
 *
 *  Created on: Mar 15, 2011
 *      Author: bhole
 */

#ifndef CRF_H_
#define CRF_H_


#include "Parameters.h"
#include "baseModel.h"

extern GlobalPairwiseParam globalP;



class crf: public baseModel // might also inherit from class MRF
{

private:
  MRF* mrf;
  MRF* mrfgen;

  DataCost* computeDataCost(int genparam, unsigned int index);
  SmoothnessCost* computeSmoothnessCost(int index);

	void regularizeLL();

	void updateQuantizedPairwiseEmpirDistributionOpflow(int index);
	void updateQuantizedPairwiseEmpirDistribution(int index);

	void updateSmoothPairwiseEmpirDistributionOpflow(int index);
	void updateSmoothPairwiseEmpirDistribution(int index);

	void updateQuantizedPairwiseModelDistributionmpe(int index);
	void updateQuantizedPairwiseModelDistributionmpeOpflow(int index);

	void updateSmoothPairwiseModelDistributionmpe(int index);
	void updateSmoothPairwiseModelDistributionmpeOpflow(int index);

	void updateQuantizedPairwiseModelDistribution(int index);
	void updateQuantizedPairwiseModelDistributionOpflow(int index);

	void updateSmoothPairwiseModelDistribution(int index);
	void updateSmoothPairwiseModelDistributionOpflow(int index);

	// hidden
	void updateQuantizedPairwiseHiddenDistributionmpe(int index);
	void updateQuantizedPairwiseHiddenDistributionmpeOpflow(int index);

	void updateSmoothPairwiseHiddenDistributionmpe(int index);
	void updateSmoothPairwiseHiddenDistributionmpeOpflow(int index);

	void updateQuantizedPairwiseHiddenDistribution(int index);
	void updateQuantizedPairwiseHiddenDistributionOpflow(int index);

	void updateSmoothPairwiseHiddenDistribution(int index);
	void updateSmoothPairwiseHiddenDistributionOpflow(int index);


	void updateDataCostDueToInteractive(int index, DataCost *dcost, SmoothnessCost *scost);

public :
	crf(Parameters* prm, std::vector<float> theta);
	~crf();

	void allocateDataCostSpace();
	void deAllocateDataCostSpace();


	void allocateSmoothCostGlobalSpace();
	void deAllocateSmoothCostGlobalSpace();


//	int preComputePairwiseFeatures(std::vector <std::vector <CImage2> > *inp);
//	int preComputePairwiseFeatures(std::vector <std::vector <CByteImage> > *inp);


	void hiddenModelProcess(unsigned int index, int iter, int iterout);
	void modelProcess(unsigned int index);

	void updateDist(int index, int fiterflag);

	void computeEmpirDist(int index);
	void computeModelDist(int index);

	void computeHiddenDist(int index);

	void computeEmpirDistU(int index);
	void computeEmpirDistA(int index);
	void computeEmpirDistH(int index);
	void computeEmpirDistM(int index);
	void computeEmpirDistL(int index);
	void computeEmpirDistW(int index);
	void computeEmpirDistO(int index);
	void computeEmpirDistV(int index);
	void computeEmpirDistB(int index);
	void computeEmpirDistE(int index);

	void computeModelDistUmpe(int index);
	void computeModelDistAmpe(int index);
	void computeModelDistHmpe(int index);
	void computeModelDistMmpe(int index);
	void computeModelDistLmpe(int index);
	void computeModelDistWmpe(int index);
	void computeModelDistOmpe(int index);
	void computeModelDistVmpe(int index);
	void computeModelDistBmpe(int index);
	void computeModelDistEmpe(int index);


	void computeModelDistU(int index);
	void computeModelDistA(int index);
	void computeModelDistH(int index);
	void computeModelDistM(int index);
	void computeModelDistL(int index);
	void computeModelDistW(int index);
	void computeModelDistO(int index);
	void computeModelDistV(int index);
	void computeModelDistB(int index);
	void computeModelDistE(int index);


	void computeHiddenDistUmpe(int index);
	void computeHiddenDistAmpe(int index);
	void computeHiddenDistHmpe(int index);
	void computeHiddenDistMmpe(int index);
	void computeHiddenDistLmpe(int index);
	void computeHiddenDistWmpe(int index);
	void computeHiddenDistOmpe(int index);
	void computeHiddenDistVmpe(int index);
	void computeHiddenDistBmpe(int index);
	void computeHiddenDistEmpe(int index);


	void computeHiddenDistU(int index);
	void computeHiddenDistA(int index);
	void computeHiddenDistH(int index);
	void computeHiddenDistM(int index);
	void computeHiddenDistL(int index);
	void computeHiddenDistW(int index);
	void computeHiddenDistO(int index);
	void computeHiddenDistV(int index);
	void computeHiddenDistB(int index);
	void computeHiddenDistE(int index);


	void deletemrf(unsigned int index);


};




#endif /* CRF_H_ */
