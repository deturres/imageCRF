/*
 * locationMV.h
 *
 *  Created on: Jan 4, 2010
 *      Author: bhole
 */

#ifndef LOCATIONMV_H_
#define LOCATIONMV_H_

#include<string>

class LocationMV
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

	LocationMV(int numC, int dim);
	~LocationMV();
	void setClassname(char cName[]);
	void readParameters(std::string filename);
	void printParameters();

	int getNClusters() { return numClusters; }

	float* getPi() { return pii; };
	float** getMu() { return mu; };
	float*** getSigmaInv() { return SigmaInv; }
	float* getNormConst() { return normConst; }

};

#endif /* LOCATIONMV_H_ */
