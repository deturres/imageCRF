/*
 * intensityClass.h
 *
 *  Created on: Feb 24, 2010
 *      Author: bhole
 */

#ifndef INTENSITYCLASS_H_
#define INTENSITYCLASS_H_

class intClass
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

	intClass(int numC, int dim);
	~intClass();
	void setClassname(char cName[]);
	void readParameters(char filename[]);
	void printParameters();

	int getNClusters() { return numClusters; }

	float* getPi() { return pii; };
	float** getMu() { return mu; };
	float*** getSigmaInv() { return SigmaInv; }
	float* getNormConst() { return normConst; }

};


#endif /* INTENSITYCLASS_H_ */
