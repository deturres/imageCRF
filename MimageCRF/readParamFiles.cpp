/*
 * readParamFiles..cpp
 *
 *  Created on: May 22, 2010
 *      Author: bhole
 */

#include "readParamFiles.h"


int readCRFparamFile(char fileName[], string& paramstr)
{
	string line;

    ifstream myfile (fileName);
    if (myfile.is_open())
    {
      while (! myfile.eof() )
      {
        getline (myfile,line);
        paramstr.append(" ");
        paramstr.append(line);

      }
      myfile.close();
    }
    else {
    	std::cout << "Unable to open parameter file";
    	exit(1);
    }

    return 0;

}





int readLocationMVparams(LocationMV**& LocationMVclass , int nD)
{

	LocationMVclass = new LocationMV*[nD];

	// read parameter files
	//read bgnd, liver, rk, lk, gb, spleen
	for (int kk=0; kk<nD; kk++)
	{
 	string line;
 	int check = 0;
		char cname[50], fname[50];
		int numC, dim;

		if (kk==0)
			strcpy(fname,"../parameters/mcovbgnd.txt");
		else if (kk==1)
			strcpy(fname,"../parameters/mcovliv.txt");
		else if (kk==2)
			strcpy(fname,"../parameters/mcovrk.txt");
		else if (kk==3)
			strcpy(fname,"../parameters/mcovlk.txt");
		else if (kk==4)
			strcpy(fname,"../parameters/mcovgb.txt");
		else if (kk==5)
			strcpy(fname,"../parameters/mcovspl.txt");


		ifstream myfile (fname);
		if (myfile.is_open())
		{
		  while (! myfile.eof() )
		  {
			  myfile >> line;
			  if (line.compare("###newclass")==0)
			  {
				  // good to go, ignore this header line
				  check = 1;
				  myfile >> cname;
				  myfile >> numC;
				  myfile >> dim;
				  break;

			  } else
			  {
				  if (check==0)
				  {
					  std::cout<<" LocationMV Parameter file corrupt \n";
					  exit(1);
				  }
			  }
		  }
		  myfile.close();
		}
		else {
			std::cout << "Unable to open LocationMV parameter file";
			exit(1);
		}

		DEBUG_OUT3(verbose, debugfile, "cname = %s numC = %d dim = %d\n", cname, numC, dim);

		LocationMVclass[kk] = new LocationMV(numC, dim);
		LocationMVclass[kk]->setClassname(cname);
		LocationMVclass[kk]->readParameters(fname);
//			LocationMVclass[kk]->printParameters();
	}

	return 0;

}




int readAppearanceMVparams(AppearanceMV**& AppMVclass , int nD, std::string &appfName)
{

	AppMVclass = new AppearanceMV*[nD];

	// read parameter files
	//read bgnd, liver, rk, lk, gb, spleen
	for (int kk=0; kk<nD; kk++)
	{
    	string line;
    	int check = 0;
		char cname[50], fname[200];
		int numC, dim;

		strcpy(fname,appfName.c_str());

		if (kk==0)
			strcat(fname,"bgnd.txt");
		else if (kk==1)
			strcat(fname,"liv.txt");
		else if (kk==2)
			strcat(fname,"rk.txt");
		else if (kk==3)
			strcat(fname,"lk.txt");
		else if (kk==4)
			strcat(fname,"gb.txt");
		else if (kk==5)
			strcat(fname,"spl.txt");


		ifstream myfile (fname);
		if (myfile.is_open())
		{
		  while (! myfile.eof() )
		  {
			  myfile >> line;
			  if (line.compare("###newclass")==0)
			  {
				  // good to go, ignore this header line
				  check = 1;
				  myfile >> cname;
				  myfile >> numC;
				  myfile >> dim;
				  break;

			  } else
			  {
				  if (check==0)
				  {
					  std::cout<<" AppearanceMV Parameter file corrupt \n";
					  exit(1);
				  }
			  }
		  }
		  myfile.close();
		}
		else {
			std::cout << "Unable to open AppearanceMV parameter file";
			exit(1);
		}

		DEBUG_OUT3(verbose, debugfile, "cname = %s numC = %d dim = %d\n", cname, numC, dim);

		AppMVclass[kk] = new AppearanceMV(numC, dim);
		AppMVclass[kk]->setClassname(cname);
		AppMVclass[kk]->readParameters(fname);
//			AppMVclass[kk]->printParameters();
	}

	return 0;
}




int readAppearanceparams(Appearance**& appclass, int nD, std::string  &appfName)
{

	appclass = new Appearance*[nD];

	// read parameter file: order of texton centers
	//read bgnd, liver, rk, lk, gb, spleen

	string line;
	int numT, patchsize;

	ifstream myfile (appfName.c_str());
	if (myfile.is_open())
	{
		  while (! myfile.eof() )
		  {
			  myfile >> line;
			  if (line.compare("numClusters")==0)
			  {
				  // good to go
				  myfile >> line; // take in "="
				  myfile >> numT;
				  myfile >> line; // take in appBoxSize
				  myfile >> line; // take in "="
				  myfile >> patchsize;
				  break;

			  } else
			  {
				  std::cout<<" appearance Parameter file corrupt \n";
				  exit(1);
			  }
		  }

		  myfile.close();

	} else {
		std::cout << "Unable to open appearance parameter file";
		exit(1);
	}

	for (int kk=0; kk<nD; kk++)
	{
		appclass[kk] = new Appearance(numT/nD, patchsize, 1);
		appclass[kk]->setClassname("general");

		appclass[kk]->readParameters(appfName, nD, kk, 1);
//			appclass[kk]->printParameters();
	}

	return 0;
}



int readAppearanceparams1(Appearance**& appclass, int nD, std::string &appfName)
{
	appclass = new Appearance*[1];

	// read parameter file: order of texton centers
	//read bgnd, liver, rk, lk, gb, spleen

	string line;
	int numT, patchsize;

	ifstream myfile (appfName.c_str());
	if (myfile.is_open())
	{
		  while (! myfile.eof() )
		  {
			  myfile >> line;
			  if (line.compare("numClusters")==0)
			  {
				  // good to go
				  myfile >> line; // take in "="
				  myfile >> numT;
				  myfile >> line; // take in appBoxSize
				  myfile >> line; // take in "="
				  myfile >> patchsize;
				  break;

			  } else
			  {
				  std::cout<<" appearance Parameter file corrupt \n";
				  exit(1);
			  }
		  }

		  myfile.close();

	} else {
		std::cout << "Unable to open appearance parameter file";
		exit(1);
	}

	for (int kk=0; kk<1; kk++)
	{
		appclass[kk] = new Appearance(numT, patchsize, 1);
		appclass[kk]->setClassname("general");

		appclass[kk]->readParameters(appfName, 1, kk, 1);
//			appclass[kk]->printParameters();
	}

	return 0;
}



int readintClassparams(intClass**& intObj, int nD)
{

    intObj = new intClass*[nD];

	// read parameter files
	//read bgnd, liver, rk, lk, gb, spleen
	for (int kk=0; kk<nD; kk++)
	{
    	string line;
    	int check = 0;
		char cname[50], fname[50];
		int numC, dim;

		if (kk==0)
			strcpy(fname,"../parameters/intmcovbgnd.txt");
		else if (kk==1)
			strcpy(fname,"../parameters/intmcovliv.txt");
		else if (kk==2)
			strcpy(fname,"../parameters/intmcovrk.txt");
		else if (kk==3)
			strcpy(fname,"../parameters/intmcovlk.txt");
		else if (kk==4)
			strcpy(fname,"../parameters/intmcovgb.txt");
		else if (kk==5)
			strcpy(fname,"../parameters/intmcovspl.txt");


		ifstream myfile (fname);
		if (myfile.is_open())
		{
		  while (! myfile.eof() )
		  {
			  myfile >> line;
			  if (line.compare("###newclass")==0)
			  {
				  // good to go, ignore this header line
				  check = 1;
				  myfile >> cname;
				  myfile >> numC;
				  myfile >> dim;
				  break;

			  } else
			  {
				  if (check==0)
				  {
					  std::cout<<" intClass Parameter file corrupt \n";
					  exit(1);
				  }
			  }
		  }
		  myfile.close();
		}
		else {
			std::cout << "Unable to open intClass parameter file";
			exit(1);
		}

		DEBUG_OUT3(verbose, debugfile, "cname = %s numC = %d dim = %d\n", cname, numC, dim);

		intObj[kk] = new intClass(numC, dim);
		intObj[kk]->setClassname(cname);
		intObj[kk]->readParameters(fname);
//			intObj[kk]->printParameters();
	}

	return 0;
}










