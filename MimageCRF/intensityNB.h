/*
 * intensityNB.h
 *
 *  Created on: Mar 14, 2010
 *      Author: bhole
 */

#ifndef INTENSITYNB_H_
#define INTENSITYNB_H_

#include "imageLib.h"
#include <vector>

class intensityNB
{
private:

	std::vector<double> PofY;
	std::vector<double> PofX;

	std::vector<std::vector<double> > PofXcY;
	std::vector<std::vector<double> > PofYcX;

	int numBins;
	int nD;
	int alpha; //smoothing factor

public:

	intensityNB(int nB, int nDx);
	~intensityNB();

	void setSmoothFactor(double a) { alpha = a; }
	void readProbValues(std::vector <std::vector <CImage2> > im1, std::vector <std::vector <CByteImage> > gtIm, std::vector<int> testDirIndexV);
	void readProbValues(std::vector <std::vector <CImage2> >* im1, std::vector <std::vector <CByteImage> >* gtIm, std::vector<int> testDirIndexV, int rangeI);
	int computePofXcY();

	int getNBins() { return numBins; }
	int getND() { return nD; }

	std::vector<double> getPofY() { return PofY; }
	std::vector<double> getPofX() { return PofX; }

	std::vector<std::vector<double> > getPofXcY() { return PofXcY; }
	std::vector<std::vector<double> > getPofYcX() { return PofYcX; }

	double getPofXcYbd(int b, int d) { return PofXcY[b][d]; }

};

#endif /* INTENSITYNB_H_ */
