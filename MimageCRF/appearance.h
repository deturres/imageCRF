/*
 * appearance.h
 *
 *  Created on: Jan 6, 2010
 *      Author: bhole
 */

#ifndef APPEARANCE_H_
#define APPEARANCE_H_

#include<string>

class Appearance
{
private:

	char *classname;
	int numTextons;
	int patchSize;
	float **textons; //assuming they are 2D for now

public:

	Appearance(int numT, int ps, int colordims);
	~Appearance();
	void setClassname(char cName[]);
	void readParameters(std::string filename, int nD, int cno, int colordims);
	void printParameters(int colordims);

	int getNTextons() { return numTextons; }
	int getPatchSize() { return patchSize; }

	float** getTextons() { return textons; };

};


#endif /* APPEARANCE_H_ */
