/*
 * evaluation.h
 *
 *  Created on: Apr 30, 2011
 *      Author: bhole
 */

#ifndef EVALUATION_H_
#define EVALUATION_H_

#include <vector>

class evaluation
{

	int size;

protected:


    // variables related to evaluation

    // confusion matrix initialization
    // for training
    double cnterr, totfull;
    double **confMatrixFull;
    double *totNDFull;  // gtruth
    double *totresNDFull;  // result

    // for test
    double cnterrtest, totfulltest;
    double **confMatrixFullTest;
    double *totNDFullTest; // gtruth
    double *totresNDFullTest; // result

    // for testgen
    double cnterrtestgen, totfulltestgen;
    double **confMatrixFullTestgen;
    double *totNDFullTestgen; // gtruth
    double *totresNDFullTestgen; // result

    // per volume
    double cnterrPV, totPV;
    double **confMatrixPV;
    double *totNDPV; // gtruth
    double *totresNDPV; // result

    // train
    std::vector<double> dcsTr;
    // double dcs_Xtr, dcs_Ytr, dcs_XinterYtr;
    // dcs_X registers segments from algorithm
    // dcs_Y registers groundtruth segmentation

    // test
    std::vector<double> dcsTe;
    // double dcs_Xte, dcs_Yte, dcs_XinterYte;

    // testgen
    std::vector<double> dcsTgen;
    // double dcs_Xtgen, dcs_Ytgen, dcs_XinterYtgen;

    // per volume
    std::vector<double> dcs;
    // double dcs_X, dcs_Y, dcs_XinterY;


    // training
    double avgCaccTr, pixAccTr;

    // testing
    double avgCaccTe, pixAccTe;

    // testing gen
    double avgCaccTgen, pixAccTgen;

    // per volume
    double avgCaccPV, pixAccPV;

public:

  evaluation(int nD);

  ~evaluation();


	void initializeConfMatrixToZero();
	void initializeConfMatrixPVToZero();

	void initializeEvalVarToZero();
	void initializeEvalVarPVToZero();



};








#endif /* EVALUATION_H_ */
