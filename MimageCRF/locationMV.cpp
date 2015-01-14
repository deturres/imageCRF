/*
 * locationMV.cpp
 *
 *  Created on: Jan 4, 2010
 *      Author: bhole
 */

#include "locationMV.h"
#include<string.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>

using std::ifstream;

LocationMV::LocationMV(int numC, int dim)
{
	classname = 0;
	numClusters = numC;
	dimension = dim;

	pii = new float[numC];
	mu = new float*[numC];
	for (int i=0; i<numC; i++)
	{
		mu[i] = new float[dim];
	}

	Sigma = new float**[numC];
	for (int i=0; i<numC; i++)
	{
		Sigma[i] = new float*[dim];
		for (int j=0; j<dim; j++)
		{
			Sigma[i][j] = new float[dim];
		}
	}

	SigmaInv = new float**[numC];
	for (int i=0; i<numC; i++)
	{
		SigmaInv[i] = new float*[dim];
		for (int j=0; j<dim; j++)
		{
			SigmaInv[i][j] = new float[dim];
		}
	}

	normConst = new float[numC];

}


LocationMV::~LocationMV()
{
	if (classname)
		delete[] classname;

	delete[] pii;
	delete[] normConst;

	for (int i=0; i<numClusters; i++)
		delete[] mu[i];
	delete[] mu;

	for (int i=0; i<numClusters; i++)
		for (int j=0; j<dimension; j++)
		{
			delete[] Sigma[i][j];
			delete[] SigmaInv[i][j];
		}

	for (int i=0; i<numClusters; i++)
	{
		delete[] Sigma[i];
		delete[] SigmaInv[i];
	}

	delete[] Sigma;
	delete[] SigmaInv;
}

void LocationMV::setClassname(char cName[])
{
	classname = new char[strlen(cName)+1];
	strcpy(classname, cName);
}

void LocationMV::printParameters()
{

	printf("Class: %s\n\n", classname);

	printf(" Displaying pi \n");
	for (int i=0; i<numClusters; i++)
	{
		printf("%7.6f  ",pii[i]);
	}
	printf("\n\n");

	printf(" Displaying mu \n");
	for (int i=0; i<numClusters; i++)
	{
		for (int j=0; j<dimension; j++)
		{
			printf("%7.6f  ",mu[i][j]);
		}
		printf("\n");
	}
	printf("\n\n");

	printf(" Displaying Sigma \n");
	for (int i=0; i<numClusters; i++)
	{
		for (int j=0; j<dimension; j++)
		{
			for (int k=0; k<dimension; k++)
			{
				printf("%7.6f  ",Sigma[i][j][k]);
			}
			printf("\n");
		}
		printf("\n\n");
	}
	printf("\n\n");

	printf(" Displaying SigmaInv \n");
	for (int i=0; i<numClusters; i++)
	{
		for (int j=0; j<dimension; j++)
		{
			for (int k=0; k<dimension; k++)
			{
				printf("%7.6f  ",SigmaInv[i][j][k]);
			}
			printf("\n");
		}
		printf("\n\n");
	}
	printf("\n\n");


	printf(" Displaying normConst \n");
	for (int i=0; i<numClusters; i++)
	{
		printf("%7.8f  ", normConst[i]);
	}
	printf("\n\n");





}


//fragile reader:: needs the numbers to be in exactly the same format
// with no additional lines
// the two numbers in header line indicate number of clusters(6) and dimension(3)
//sample file looks like this:
//###newclass bgnd 6 3
//pi1 pi2 pi3 pi4 pi5 pi6
//mu11 mu12 mu13
//mu21 ..
//..
//Sigma111 Sigma112 ..
//..
//SimgaInv ....
//normConst1 normConst2 normConst3

void LocationMV::readParameters(std::string filename)
{
	char line[100];
	int check = 0;

    ifstream myfile (filename.c_str());
	if (myfile.is_open())
	{
	  while (! myfile.eof() )
	  {
		  if (check==1) {
			  myfile >> line; //classname
			  myfile >> line; //numClusters
			  myfile >> line; //dimension
		  } else {
			  myfile >> line; // will have to read ###newclass or it will abort
		  }

		  if (strcmp(line, "###newclass")==0)
		  {
			  // good to go, ignore this header line
			  check = 1;
		  } else
		  {
			  if (check==0)
			  {
				  std::cout<<" MOG Parameter file corrupt \n";
				  exit(1);
			  }

			  // update pi
			  for (int i=0; i<numClusters; i++)
			  {
				  myfile >> line;
				  pii[i] = (float)atof(line);
			  }

			  // update mu
			  for (int i=0; i<numClusters; i++)
			  {
				  for (int j=0; j<dimension; j++)
				  {
					  myfile >> line;
					  mu[i][j] = (float)atof(line);
				  }
			  }

			  // update sigma
			  for (int i=0; i<numClusters; i++)
			  {
				  for (int j=0; j<dimension; j++)
				  {
					  for (int k=0; k<dimension; k++)
					  {
						  myfile >> line;
						  Sigma[i][j][k] = (float)atof(line);
					  }
				  }
			  }

			  // update sigmaInv
			  for (int i=0; i<numClusters; i++)
			  {
				  for (int j=0; j<dimension; j++)
				  {
					  for (int k=0; k<dimension; k++)
					  {
						  myfile >> line;
						  SigmaInv[i][j][k] = (float)atof(line);
					  }
				  }
			  }

			  // update normConst
			  for (int i=0; i<numClusters; i++)
			  {
				  myfile >> line;
				  normConst[i] = (float)atof(line);
			  }

			  break;
		  }
	  }
	  myfile.close();
	}
	else {
		std::cout << "Unable to open MOG parameter file : " << filename;
		exit(1);
	}

}

