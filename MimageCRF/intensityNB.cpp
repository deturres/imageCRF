/*
 * intensityNB.cpp
 *
 *  Created on: Mar 14, 2010
 *      Author: bhole
 */

#include "intensityNB.h"
#include <stdio.h>
#include <iostream>
using namespace std;

intensityNB::~intensityNB()
{

}

intensityNB::intensityNB(int nB, int nDx)
{

	numBins = nB;
	nD = nDx;

	PofY.resize(nD);
	PofX.resize(numBins);

	PofXcY.resize(numBins);
	for (int i=0; i<numBins; i++)
	{
		PofXcY[i].resize(nD);
	}

	PofYcX.resize(nD);
	for (int i=0; i<nD; i++)
	{
		PofYcX[i].resize(numBins);
	}

}

void intensityNB::readProbValues(vector <vector <CImage2> > im1, vector <vector <CByteImage> > gtIm, std::vector<int> testDirIndexV)
{

	double binner = 1380.0/numBins;

	int numP = im1.size();
	int depth = im1[0].size();

	CShape sh = im1[0][0].Shape();
	int width = sh.width, height = sh.height;


	for (unsigned int j = 0; j < im1.size(); ++j)
	{

		// for each training patient else neglect test case
		if (testDirIndexV[j]==1) {
			continue;
		}

		for(unsigned int i = 0; i < im1[j].size(); ++i)
		{

			  for (int y = 0; y < height; y++)
			  {
				  for (int x = 0; x < width; x++)
				  {
					  unsigned short pix1 = im1[j][i].Pixel(x, y, 0);
					  uchar *gtpix = &gtIm[j][i].Pixel(x, y, 0);

					  int binid;
					  if (pix1 >= 1380)
						  binid = numBins-1;
					  else
						  binid = int(pix1/binner);

					  PofY[(int) *gtpix] ++;
					  PofX[binid] ++;
					  PofYcX[(int) *gtpix][binid] ++;

				  }
			  }
		}
	}


	// normalize probabilities

	double sumPofY=0.0;
	for (int d=0; d<nD; d++)
		sumPofY += PofY[d];

	for (int d=0; d<nD; d++)
		PofY[d] = PofY[d]/sumPofY;


	double sumPofX=0.0;
	for (int b=0; b<numBins; b++)
		sumPofX += PofX[b];

	for (int b=0; b<numBins; b++)
		PofX[b] = PofX[b]/sumPofX;

	for (int b=0; b<numBins; b++)
	{
		double sumPofYcX=0.0;
		for (int d=0; d<nD; d++)
			sumPofYcX += PofYcX[d][b];

		for (int d=0; d<nD; d++)
			PofYcX[d][b] = PofYcX[d][b]/sumPofYcX;

	}
}



void intensityNB::readProbValues(vector <vector <CImage2> >* im1, vector <vector <CByteImage> >* gtIm, std::vector<int> testDirIndexV, int rangeI)
{

	double binner = double(rangeI)/numBins;

	int numP = (*im1).size();
	int depth = (*im1)[0].size();

	CShape sh = (*im1)[0][0].Shape();
	int width = sh.width, height = sh.height;


	for (unsigned int j = 0; j < (*im1).size(); ++j)
	{

		// for each training patient else neglect test case
		if (testDirIndexV[j]==1) {
			continue;
		}

		for(unsigned int i = 0; i < (*im1)[j].size(); ++i)
		{

			  for (int y = 0; y < height; y++)
			  {
				  for (int x = 0; x < width; x++)
				  {
					  unsigned short pix1 = (*im1)[j][i].Pixel(x, y, 0);
					  uchar *gtpix = &(*gtIm)[j][i].Pixel(x, y, 0);

					  int binid;
					  if (pix1 >= rangeI)
						  binid = numBins-1;
					  else
						  binid = int(pix1/binner);

					  PofY[(int) *gtpix] ++;
					  PofX[binid] ++;
					  PofYcX[(int) *gtpix][binid] ++;

				  }
			  }
		}
	}


	// normalize probabilities

	double sumPofY=0.0;
	for (int d=0; d<nD; d++)
		sumPofY += PofY[d];

	for (int d=0; d<nD; d++)
		PofY[d] = PofY[d]/sumPofY;


	double sumPofX=0.0;
	for (int b=0; b<numBins; b++)
		sumPofX += PofX[b];

	for (int b=0; b<numBins; b++)
		PofX[b] = PofX[b]/sumPofX;

	for (int b=0; b<numBins; b++)
	{
		double sumPofYcX=0.0;
		for (int d=0; d<nD; d++)
			sumPofYcX += PofYcX[d][b];

		for (int d=0; d<nD; d++)
			PofYcX[d][b] = PofYcX[d][b]/sumPofYcX;

	}
}









int intensityNB::computePofXcY()
{

	for (int b=0; b<numBins; b++)
	{
		for (int d=0; d<nD; d++)
		{
			PofXcY[b][d] = PofYcX[d][b]*PofX[b]/(PofY[d] + alpha);
//			std::cout <<"  ["<<b<<"]["<<d<<"]  "<<PofXcY[b][d]<<std::endl;
		}
	}

	return 0;
}

