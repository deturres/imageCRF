/*
 * readParamFiles.h
 *
 *  Created on: May 22, 2010
 *      Author: bhole
 */

#ifndef READPARAMFILES_H_
#define READPARAMFILES_H_

#include "locationMV.h"
#include "appearance.h"
#include "appearanceMV.h"
#include "intensityClass.h"
#include "intensityNB.h"
#include "crfmodel.h"

#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;
using std::fstream;

int readCRFparamFile(char fileName[], string& paramstr);
int readLocationMVparams(LocationMV**& LocationMVclass , int nD);
int readAppearanceMVparams(AppearanceMV**& AppMVclass , int nD, std::string  &appfName);
int readAppearanceparams(Appearance**& appclass, int nD, std::string &appfName);
int readAppearanceparams1(Appearance**& appclass, int nD, std::string &appfName);
int readintClassparams(intClass**& intObj, int nD);


#endif /* READPARAMFILES_H_ */
