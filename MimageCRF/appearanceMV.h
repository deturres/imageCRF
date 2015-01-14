/*
 * appearanceMV.h
 *
 *  Created on: Mar 9, 2010
 *      Author: bhole
 */

#ifndef APPEARANCEMV_H_
#define APPEARANCEMV_H_

class AppearanceMV
{
private:

	char *classname;
	int numClusters;
	int dimension;
	float *pii;
	float **mu;
	float ***Sigma;
	float ***SigmaInv;
	float *normConst;

public:

	AppearanceMV(int numC, int dim);
	~AppearanceMV();
	void setClassname(char cName[]);
	void readParameters(char filename[]);
	void printParameters();

	int getNClusters() { return numClusters; }

	int getDimension() {return dimension; }
	float* getPi() { return pii; };
	float** getMu() { return mu; };
	float*** getSigmaInv() { return SigmaInv; }
	float* getNormConst() { return normConst; }

};




#endif /* APPEARANCEMV_H_ */
