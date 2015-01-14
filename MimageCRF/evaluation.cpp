/*
 * evaluation.cpp
 *
 *  Created on: Apr 30, 2011
 *      Author: bhole
 */

#include "evaluation.h"


evaluation::evaluation(int nD)
{

	size = nD;

	// confusion matrix initialization
	cnterr = 0;
	totfull = 0;

	confMatrixFull = new double*[size];
	for (int i=0; i<size; i++) {
	  confMatrixFull[i] = new double[size];
	}

	totNDFull = new double[size];
	totresNDFull = new double[size];

	for (int i=0; i<size; i++) {
	  totNDFull[i] = 0;
	  totresNDFull[i] = 0;
	  for (int j=0; j<size; j++) {
		confMatrixFull[i][j] =  0;
	  }
	}

	// for test
	cnterrtest = 0;
	totfulltest = 0;

	confMatrixFullTest = new double*[size];
	for (int i=0; i<size; i++) {
	  confMatrixFullTest[i] = new double[size];
	}

	totNDFullTest = new double[size];
	totresNDFullTest = new double[size];

	for (int i=0; i<size; i++) {
	  totNDFullTest[i] = 0;
	  totresNDFullTest[i] = 0;
	  for (int j=0; j<size; j++) {
		confMatrixFullTest[i][j] =  0;
	  }
	}



	// for testgen
	cnterrtestgen = 0;
	totfulltestgen = 0;

	confMatrixFullTestgen = new double*[size];
	for (int i=0; i<size; i++) {
	  confMatrixFullTestgen[i] = new double[size];
	}

	totNDFullTestgen = new double[size];
	totresNDFullTestgen = new double[size];

	for (int i=0; i<size; i++) {
	  totNDFullTestgen[i] = 0;
	  totresNDFullTestgen[i] = 0;
	  for (int j=0; j<size; j++) {
		confMatrixFullTestgen[i][j] =  0;
	  }
	}



	// per volume

	cnterrPV = 0;
	totPV = 0;

	confMatrixPV = new double*[size];
	for (int i=0; i<size; i++) {
	  confMatrixPV[i] = new double[size];
	}

	totNDPV = new double[size];
	totresNDPV = new double[size];

	for (int i=0; i<size; i++) {
	  totNDPV[i] = 0;
	  totresNDPV[i] = 0;
	  for (int j=0; j<size; j++) {
		confMatrixPV[i][j] =  0;
	  }
	}


	dcs.resize(size, 0);
	dcsTr.resize(size, 0);
	dcsTe.resize(size, 0);
	dcsTgen.resize(size, 0);

	avgCaccTr = 0;
	pixAccTr = 0;

	avgCaccTe = 0;
	pixAccTe = 0;

	avgCaccTgen = 0;
	pixAccTgen = 0;

	avgCaccPV = 0;
	pixAccPV = 0;



}





evaluation::~evaluation()
{

	delete[] totNDFull;
	delete[] totresNDFull;
	for (int i=0; i<size; i++)
		delete[] confMatrixFull[i];
	delete[] confMatrixFull;


	delete[] totNDFullTest;
	delete[] totresNDFullTest;
	for (int i=0; i<size; i++)
		delete[] confMatrixFullTest[i];
	delete[] confMatrixFullTest;

	delete[] totNDFullTestgen;
	delete[] totresNDFullTestgen;
	for (int i=0; i<size; i++)
		delete[] confMatrixFullTestgen[i];
	delete[] confMatrixFullTestgen;


	delete[] totNDPV;
	delete[] totresNDPV;
	for (int i=0; i<size; i++)
		delete[] confMatrixPV[i];
	delete[] confMatrixPV;

}






void evaluation::initializeEvalVarToZero()
{

	initializeConfMatrixToZero();

	std::fill(dcs.begin(),dcs.end(),0);
	std::fill(dcsTr.begin(),dcsTr.end(),0);
	std::fill(dcsTe.begin(),dcsTe.end(),0);
	std::fill(dcsTgen.begin(),dcsTgen.end(),0);

	avgCaccTr = 0;
	pixAccTr = 0;

	avgCaccTe = 0;
	pixAccTe = 0;

	avgCaccTgen = 0;
	pixAccTgen = 0;

	avgCaccPV = 0;
	pixAccPV = 0;

}



void evaluation::initializeEvalVarPVToZero()
{

	// per volume
	// confusion matrix initialization
	initializeConfMatrixPVToZero();

	std::fill(dcs.begin(),dcs.end(),0);

	avgCaccPV = 0;
	pixAccPV = 0;

}



void evaluation::initializeConfMatrixToZero()
{

	// confusion matrix initialization
	cnterr = 0;
	totfull = 0;

	for (int i=0; i<size; i++) {
	  totNDFull[i] = 0;
	  totresNDFull[i] = 0;
	  for (int j=0; j<size; j++) {
		confMatrixFull[i][j] =  0;
	  }
	}

	// for test
	cnterrtest = 0;
	totfulltest = 0;

	for (int i=0; i<size; i++) {
	  totNDFullTest[i] = 0;
	  totresNDFullTest[i] = 0;
	  for (int j=0; j<size; j++) {
		confMatrixFullTest[i][j] =  0;
	  }
	}


	// for testgen
	cnterrtestgen = 0;
	totfulltestgen = 0;

	for (int i=0; i<size; i++) {
	  totNDFullTestgen[i] = 0;
	  totresNDFullTestgen[i] = 0;
	  for (int j=0; j<size; j++) {
		confMatrixFullTestgen[i][j] =  0;
	  }
	}


	// per volume
	// confusion matrix initialization
	cnterrPV = 0;
	totPV = 0;

	for (int i=0; i<size; i++) {
	  totNDPV[i] = 0;
	  totresNDPV[i] = 0;
	  for (int j=0; j<size; j++) {
		confMatrixPV[i][j] =  0;
	  }
	}

}



void evaluation::initializeConfMatrixPVToZero()
{

	// per volume
	// confusion matrix initialization
	cnterrPV = 0;
	totPV = 0;

	for (int i=0; i<size; i++)
	{
	  totNDPV[i] = 0;
	  totresNDPV[i] = 0;
	  for (int j=0; j<size; j++) {
		confMatrixPV[i][j] =  0;
	  }
	}

}



















