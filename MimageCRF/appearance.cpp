/*
 * appearance.cpp
 *
 *  Created on: Jan 6, 2010
 *      Author: bhole
 */

#include "appearance.h"
#include<string.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>

using std::ifstream;

Appearance::Appearance(int numT, int ps, int colordims)
{
	classname = 0;
	numTextons = numT;
	patchSize = ps;

	textons = new float*[numT];
	for (int i=0; i<numT; i++)
	{
		textons[i] = new float[ps*ps*colordims];
	}
}


Appearance::~Appearance()
{
	if (classname)
		delete[] classname;

	for (int i=0; i<numTextons; i++)
		delete[] textons[i];
	delete[] textons;
}

void Appearance::setClassname(char cName[])
{
	classname = new char[strlen(cName)+1];
	strcpy(classname, cName);
}

void Appearance::printParameters(int colordims)
{
	printf("Class: %s\n\n", classname);

	printf(" Displaying Textons \n");

	for (int i=0; i<numTextons; i++)
	{
		for (int j=0; j<patchSize*patchSize*colordims; j++)
		{
			printf("%7.6f  ",textons[i][j]);
		}
		printf("\n");
	}
	printf("\n\n");
}


//fragile reader:: needs the numbers to be in exactly the same format
// with no additional lines except two header lines
// one for total number of textons (ie. textons per class * number of classes
// other line contains patch size
//sample file looks like this:

//numClusters = 60
//appBoxSize = 5
//12.2 1222.3 ....
//..

void Appearance::readParameters(std::string filename, int nD, int cno, int colordims)
{
	// cno indicates the class number of which Textons we want to read
	// hence we have to skip number of lines depending on that
	// cno ranges from 0 to nD-1

	char line[5000];

    ifstream myfile (filename.c_str());
	if (myfile.is_open())
	{
	  while (! myfile.eof() )
	  {
		  // skip first two lines
		  myfile.getline(line, 5000);

		  line[11]=0;
		  if (strcmp(line, "numClusters")==0)
		  {
			  // good to go, ignore this header line

			  // read second header line
			  myfile.getline(line, 5000);

			  // need to jump to correct class
			  int linestojump = cno*numTextons;
			  for (int jj=0; jj<linestojump; jj++) {
				  myfile.getline(line, 5000);
			  }

			  // update mu
			  for (int i=0; i<numTextons; i++)
			  {
				  for (int j=0; j<patchSize*patchSize*colordims; j++)
				  {
					  myfile >> line;
					  textons[i][j] = (float)atof(line);
				  }
			  }

			  break;

		  } else
		  {
			  std::cout<<" appearance Parameter file corrupt \n";
			  exit(1);
		  }
	  }
	  myfile.close();
	}
	else {
		std::cout << "Unable to open appearance parameter file";
		exit(1);
	}

}


