// evaldisps.cpp

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "crfmodel.h"


#define ERRLIMIT 4  /* highest absolute disparity error to count */
#define UNK 0     /* color for unknown true disparity and occlusion */
int errcol[ERRLIMIT+1] = {255, 200, 40, 20, 0};  // colors for err==0, 1, 2, 3, and >3

using std::cout;


void evaldisps(std::vector<CByteImage> disp, std::vector<CByteImage> truedisp, std::vector<CByteImage> &errormap,
               float &errpercent, int nD, double &cnterrfull, double &totfull, double **confMatrixFull, double *totNDFull, int paramlreg)
{

  CShape sh = disp[0].Shape();
  if (sh != truedisp[0].Shape())
	throw CError("image shapes don't match");
  int width = sh.width, height = sh.height;
  int depth = disp.size();

//  double confMatrix[nD][nD];
//  double totND[nD];

  double **confMatrix;
  confMatrix = new double*[nD];
  for (int i=0; i<nD; i++) {
	  confMatrix[i] = new double[nD];
  }

  double *totND = new double[nD];

  //initialize local confusion matrix
  // this is for the one patient only
  //while confMatrixFull (passed as a parameter) is the confusion matrix
  //for all patients ie all training data
  for (int i=0; i<nD; i++) {
	  totND[i] = 0;
	  for (int j=0; j<nD; j++) {
		  confMatrix[i][j] =  0;
	  }
  }

  for (int z=0; z<depth; z++)
	  errormap[z].ReAllocate(sh);

  int cnterr = 0;
  int total = depth*height*width;

  for (int z=0; z < depth; z++) {
	  for (int y = 0; y < height; y++) {
		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *tdi = &truedisp[z].Pixel(0, y, 0);
		uchar *errmap = &errormap[z].Pixel(0, y, 0);
		for (int x = 0; x < width; x++) {

	      if ((paramlreg == 2 || paramlreg == 3)  && x==width-1 && width%2==1)
	    	continue;

	      if (paramlreg == 4  && (y==0 || x==0 || y==height-1 ||  x==width-1))
	        continue;

          if (paramlreg == 5  && (y==0 || x==0 || (y==height-1 && height%2==0) || (y>=height-2 && height%2==1) ||  (x==width-1 && width%2==0) || (x>=width-2 && width%2==1)  ))
	        continue;


		  confMatrix[tdi[x]][di[x]] = confMatrix[tdi[x]][di[x]] + 1;
		  totND[tdi[x]] = totND[tdi[x]] + 1;

		  if (di[x]==tdi[x]) {
			  errmap[x] = 0;
		  } else {
			  cnterr++;
			  errmap[x] = 1;
		  }
		}
	  }
  }

  for (int i=0; i<nD; i++) {
	  totNDFull[i] += totND[i];
	  for (int j=0; j<nD; j++) {
		  // confMatrixFull[i*nD+j] += confMatrix[i][j];
		  confMatrixFull[i][j] += confMatrix[i][j];
	  }
  }

  cnterrfull += cnterr;
  totfull += total;



  std::vector<double> dcs;
  dcs.resize(nD);
  double dcs_X, dcs_Y, dcs_XinterY;
  // dcs_X registers segments from algorithm
  // dcs_Y registers groundtruth segmentation


  for (int i=0; i<nD; i++) {

	  dcs_X = 0;
	  dcs_Y = 0;
	  dcs_XinterY = 0;

	  dcs_XinterY = confMatrix[i][i];

	  for (int j=0; j<nD; j++)
		  dcs_X += confMatrix[i][j];

	  for (int j=0; j<nD; j++)
		  dcs_Y += confMatrix[j][i];

	  if (dcs_X + dcs_Y == 0)
		  dcs[i] = 1.0;
	  else
		  dcs[i] = 2*dcs_XinterY / (dcs_X + dcs_Y);

  }



  for (int i=0; i<nD; i++) {
	  for (int j=0; j<nD; j++) {
		  if (totND[i] == 0)
			  confMatrix[i][j] = 0.0;
		  else
			  confMatrix[i][j] =  confMatrix[i][j]/totND[i];
	  }
  }
  errpercent = (float)(100.0 * (float)cnterr/(float)total);

  DEBUG_OUT1(verbose, debugfile, "\n Training Case errpercent= %g \n", errpercent);

  DEBUG_OUT0(verbose, debugfile, "\n \t \t Training Confusion Matrix \n");

  string dumpstr="\t";

  for (int i = 0; i<nD; i++) {
	  dumpstr += "\t";

	  std::stringstream out;
	  out << i;
	  dumpstr += out.str();
  }
  dumpstr += "\n";
  DEBUG_OUT0(verbose, debugfile, dumpstr.c_str());

  for (int i = 0; i<nD; i++) {
	  dumpstr = "\t";
	  std::stringstream out;
	  out << i;
	  dumpstr += out.str();
	  dumpstr += "\t";
	  for (int j=0; j<nD; j++) {
		  std::stringstream out2;
		  out2 << fixed << setprecision(3) << (double)confMatrix[i][j];
		  dumpstr += out2.str();
		  dumpstr += "\t";
	  }
	  dumpstr += "\n";
	  DEBUG_OUT0(verbose, debugfile, dumpstr.c_str());
  }


  DEBUG_OUT0(verbose, debugfile, "\n \t \t Dice Coefficients \n");

  dumpstr = "\t";
  for (int i = 0; i<nD; i++) {
	  std::stringstream out2;
	  out2 << fixed << setprecision(3) << (double)dcs[i] << "\t";
	  dumpstr += out2.str();
  }
  dumpstr += "\n\n";
  DEBUG_OUT0(verbose, debugfile, dumpstr.c_str());

  double avgclassacc = 0.0;

  for (int i = 0; i<nD; i++) {
	  if (totND[i] == 0)
		  avgclassacc += 1.0;
	  else
		  avgclassacc += (double)confMatrix[i][i];
  }

  avgclassacc /= nD;

  DEBUG_OUT1(verbose, debugfile, "\n Average Class Accuracy = %g \n", avgclassacc);



  delete[] totND;
  for (int i=0; i<nD; i++)
	  delete[] confMatrix[i];
  delete[] confMatrix;


}

void evaldisps(std::vector<CByteImage> disp, std::vector<CByteImage> truedisp, std::vector<CByteImage> &errormap,
               float &errpercent, int nD, int testcase, int paramlreg, std::vector<CByteImage> interImage, int interactive)
{
	// this function is more specific for the test directory
	// so we don't add up the error counts of the training data
	// nor do we add up counts in the confusion matrix

  CShape sh = disp[0].Shape();
  if (sh != truedisp[0].Shape())
	throw CError("image shapes don't match");
  int width = sh.width, height = sh.height;
  int depth = disp.size();

//  double confMatrix[nD][nD];
//  double totND[nD];

  double **confMatrix;
  confMatrix = new double*[nD];
  for (int i=0; i<nD; i++) {
	  confMatrix[i] = new double[nD];
  }

  double *totND = new double[nD];


  for (int i=0; i<nD; i++) {
	  totND[i] = 0;
	  for (int j=0; j<nD; j++) {
		  confMatrix[i][j] =  0;
	  }
  }

  for (int z=0; z<depth; z++)
	  errormap[z].ReAllocate(sh);

  int cnterr = 0;
  int total = depth*height*width;

  int totalinter = 0;

  for (int z=0; z < depth; z++) {
	  for (int y = 0; y < height; y++) {
		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *tdi = &truedisp[z].Pixel(0, y, 0);
		uchar *errmap = &errormap[z].Pixel(0, y, 0);
		for (int x = 0; x < width; x++) {

	      if ((paramlreg == 2 || paramlreg == 3)  && x==width-1 && width%2==1)
	    	continue;

	      if (paramlreg == 4  && (y==0 || x==0 || y==height-1 ||  x==width-1))
	        continue;

          if (paramlreg == 5  && (y==0 || x==0 || (y==height-1 && height%2==0) || (y>=height-2 && height%2==1) ||  (x==width-1 && width%2==0) || (x>=width-2 && width%2==1)  ))
	        continue;

          if (interactive==0)
          {
        	  confMatrix[tdi[x]][di[x]] = confMatrix[tdi[x]][di[x]] + 1;
			  totND[tdi[x]] = totND[tdi[x]] + 1;
          } else
          {
			  uchar *inimp = &interImage[z].Pixel(x, y, 0);
			  if (*inimp==0)
			  {
				  confMatrix[tdi[x]][di[x]] = confMatrix[tdi[x]][di[x]] + 1;
				  totND[tdi[x]] = totND[tdi[x]] + 1;
				  totalinter++;
			  } else
			  {
				  // as if the pixels are not there, so no need of updating matrix
			  }
          }

		  if (di[x]==tdi[x]) {
			  errmap[x] = 0;
		  } else {
			  cnterr++;
			  errmap[x] = 1;
		  }
		}
	  }
  }

  if (interactive == 1)
	  total = totalinter;


  std::vector<double> dcs;
  dcs.resize(nD);
  double dcs_X, dcs_Y, dcs_XinterY;
  // dcs_X registers segments from algorithm
  // dcs_Y registers groundtruth segmentation


  for (int i=0; i<nD; i++) {

	  dcs_X = 0;
	  dcs_Y = 0;
	  dcs_XinterY = 0;

	  dcs_XinterY = confMatrix[i][i];

	  for (int j=0; j<nD; j++)
		  dcs_X += confMatrix[i][j];

	  for (int j=0; j<nD; j++)
		  dcs_Y += confMatrix[j][i];

	  if (dcs_X + dcs_Y == 0)
		  dcs[i] = 1.0;
	  else
		  dcs[i] = 2*dcs_XinterY / (dcs_X + dcs_Y);

  }


  for (int i=0; i<nD; i++) {
	  for (int j=0; j<nD; j++) {
		  if (totND[i]==0)
			  confMatrix[i][j] = 0.0;
		  else
			  confMatrix[i][j] =  confMatrix[i][j]/totND[i];
	  }
  }
  errpercent = (float)(100.0 * (float)cnterr/(float)total);

  if (testcase == 1) {
	  DEBUG_OUT1(verbose, debugfile, "\n Test Case errpercent= %g \n", errpercent);
	  DEBUG_OUT0(verbose, debugfile, "\n \t \t Test Confusion Matrix \n");
  } else {
	  DEBUG_OUT1(verbose, debugfile, "\n Training Case errpercent= %g \n", errpercent);
	  DEBUG_OUT0(verbose, debugfile, "\n \t \t Training Confusion Matrix \n");
  }

  string dumpstr="\t";

  for (int i = 0; i<nD; i++) {
	  dumpstr += "\t";

	  std::stringstream out;
	  out << i;
	  dumpstr += out.str();
  }
  dumpstr += "\n";
  DEBUG_OUT0(verbose, debugfile, dumpstr.c_str());

  for (int i = 0; i<nD; i++) {
	  dumpstr = "\t";
	  std::stringstream out;
	  out << i;
	  dumpstr += out.str();
	  dumpstr += "\t";
	  for (int j=0; j<nD; j++) {
		  std::stringstream out2;
		  out2 << fixed << setprecision(3) << (double)confMatrix[i][j];
		  dumpstr += out2.str();
		  dumpstr += "\t";
	  }
	  dumpstr += "\n";
	  DEBUG_OUT0(verbose, debugfile, dumpstr.c_str());
  }


  DEBUG_OUT0(verbose, debugfile, "\n \t \t Dice Coefficients \n");

  dumpstr = "\t";
  for (int i = 0; i<nD; i++) {
	  std::stringstream out2;
	  out2 << fixed << setprecision(3) << (double)dcs[i] << "\t";
	  dumpstr += out2.str();
  }
  dumpstr += "\n\n";
  DEBUG_OUT0(verbose, debugfile, dumpstr.c_str());


  double avgclassacc = 0.0;

  for (int i = 0; i<nD; i++) {
	  if (totND[i]==0)
		  avgclassacc += 1.0;
	  else
		  avgclassacc += (double)confMatrix[i][i];
  }

  avgclassacc /= nD;

  DEBUG_OUT1(verbose, debugfile, "\n Average Class Accuracy = %g \n", avgclassacc);


  delete[] totND;
  for (int i=0; i<nD; i++)
	  delete[] confMatrix[i];
  delete[] confMatrix;

}



void evaldisps(std::vector<CByteImage> disp, std::vector<CByteImage> truedisp, std::vector<CByteImage> &errormap,
               float &errpercent, int nD, int testcase, int paramlreg, std::vector<CByteImage> interImage, int interactive, double &cnterrfull, double &totfull, double **confMatrixFull, double *totNDFull)
{
	// this function is more specific for the test directory
	// so we don't add up the error counts of the training data
	// nor do we add up counts in the confusion matrix

  CShape sh = disp[0].Shape();
  if (sh != truedisp[0].Shape())
	throw CError("image shapes don't match");
  int width = sh.width, height = sh.height;
  int depth = disp.size();

//  double confMatrix[nD][nD];
//  double totND[nD];

  double **confMatrix;
  confMatrix = new double*[nD];
  for (int i=0; i<nD; i++) {
	  confMatrix[i] = new double[nD];
  }

  double *totND = new double[nD];


  for (int i=0; i<nD; i++) {
	  totND[i] = 0;
	  for (int j=0; j<nD; j++) {
		  confMatrix[i][j] =  0;
	  }
  }

  for (int z=0; z<depth; z++)
	  errormap[z].ReAllocate(sh);

  int cnterr = 0;
  int total = depth*height*width;

  int totalinter = 0;

  for (int z=0; z < depth; z++) {
	  for (int y = 0; y < height; y++) {
		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *tdi = &truedisp[z].Pixel(0, y, 0);
		uchar *errmap = &errormap[z].Pixel(0, y, 0);
		for (int x = 0; x < width; x++) {

	      if ((paramlreg == 2 || paramlreg == 3)  && x==width-1 && width%2==1)
	    	continue;

	      if (paramlreg == 4  && (y==0 || x==0 || y==height-1 ||  x==width-1))
	        continue;

          if (paramlreg == 5  && (y==0 || x==0 || (y==height-1 && height%2==0) || (y>=height-2 && height%2==1) ||  (x==width-1 && width%2==0) || (x>=width-2 && width%2==1)  ))
	        continue;

          if (interactive==0)
          {
        	  confMatrix[tdi[x]][di[x]] = confMatrix[tdi[x]][di[x]] + 1;
			  totND[tdi[x]] = totND[tdi[x]] + 1;
          } else
          {
			  uchar *inimp = &interImage[z].Pixel(x, y, 0);
			  if (*inimp==0)
			  {
				  confMatrix[tdi[x]][di[x]] = confMatrix[tdi[x]][di[x]] + 1;
				  totND[tdi[x]] = totND[tdi[x]] + 1;
				  totalinter++;
			  } else
			  {
				  // as if the pixels are not there, so no need of updating matrix
			  }
          }

		  if (di[x]==tdi[x]) {
			  errmap[x] = 0;
		  } else {
			  cnterr++;
			  errmap[x] = 1;
		  }
		}
	  }
  }

  if (interactive == 1)
	  total = totalinter;

  for (int i=0; i<nD; i++) {
	  totNDFull[i] += totND[i];
	  for (int j=0; j<nD; j++) {
		  // confMatrixFull[i*nD+j] += confMatrix[i][j];
		  confMatrixFull[i][j] += confMatrix[i][j];
	  }
  }

  cnterrfull += cnterr;
  totfull += total;





  std::vector<double> dcs;
  dcs.resize(nD);
  double dcs_X, dcs_Y, dcs_XinterY;
  // dcs_X registers segments from algorithm
  // dcs_Y registers groundtruth segmentation


  for (int i=0; i<nD; i++) {

	  dcs_X = 0;
	  dcs_Y = 0;
	  dcs_XinterY = 0;

	  dcs_XinterY = confMatrix[i][i];

	  for (int j=0; j<nD; j++)
		  dcs_X += confMatrix[i][j];

	  for (int j=0; j<nD; j++)
		  dcs_Y += confMatrix[j][i];

	  if (dcs_X + dcs_Y == 0)
		  dcs[i] = 1.0;
	  else
		  dcs[i] = 2*dcs_XinterY / (dcs_X + dcs_Y);

  }


  for (int i=0; i<nD; i++) {
	  for (int j=0; j<nD; j++) {
		  if (totND[i]==0)
			  confMatrix[i][j] = 0.0;
		  else
			  confMatrix[i][j] =  confMatrix[i][j]/totND[i];
	  }
  }
  errpercent = (float)(100.0 * (float)cnterr/(float)total);

  if (testcase == 1) {
	  DEBUG_OUT1(verbose, debugfile, "\n Test Case errpercent= %g \n", errpercent);
	  DEBUG_OUT0(verbose, debugfile, "\n \t \t Test Confusion Matrix \n");
  } else {
	  DEBUG_OUT1(verbose, debugfile, "\n Training Case errpercent= %g \n", errpercent);
	  DEBUG_OUT0(verbose, debugfile, "\n \t \t Training Confusion Matrix \n");
  }

  string dumpstr="\t";

  for (int i = 0; i<nD; i++) {
	  dumpstr += "\t";

	  std::stringstream out;
	  out << i;
	  dumpstr += out.str();
  }
  dumpstr += "\n";
  DEBUG_OUT0(verbose, debugfile, dumpstr.c_str());

  for (int i = 0; i<nD; i++) {
	  dumpstr = "\t";
	  std::stringstream out;
	  out << i;
	  dumpstr += out.str();
	  dumpstr += "\t";
	  for (int j=0; j<nD; j++) {
		  std::stringstream out2;
		  out2 << fixed << setprecision(3) << (double)confMatrix[i][j];
		  dumpstr += out2.str();
		  dumpstr += "\t";
	  }
	  dumpstr += "\n";
	  DEBUG_OUT0(verbose, debugfile, dumpstr.c_str());
  }


  DEBUG_OUT0(verbose, debugfile, "\n \t \t Dice Coefficients \n");

  dumpstr = "\t";
  for (int i = 0; i<nD; i++) {
	  std::stringstream out2;
	  out2 << fixed << setprecision(3) << (double)dcs[i] << "\t";
	  dumpstr += out2.str();
  }
  dumpstr += "\n\n";
  DEBUG_OUT0(verbose, debugfile, dumpstr.c_str());


  double avgclassacc = 0.0;

  for (int i = 0; i<nD; i++) {
	  if (totND[i]==0)
		  avgclassacc += 1.0;
	  else
		  avgclassacc += (double)confMatrix[i][i];
  }

  avgclassacc /= nD;

  DEBUG_OUT1(verbose, debugfile, "\n Average Class Accuracy = %g \n", avgclassacc);


  delete[] totND;
  for (int i=0; i<nD; i++)
	  delete[] confMatrix[i];
  delete[] confMatrix;

}




void userlabelfix(std::vector<CByteImage> &disp, std::vector<CByteImage> truedisp, std::vector<CByteImage> &errormap, std::vector<CByteImage> interImage)
{
	// this function is more specific for the test directory
	// so we don't add up the error counts of the training data
	// nor do we add up counts in the confusion matrix

  CShape sh = disp[0].Shape();
  if (sh != truedisp[0].Shape())
	throw CError("image shapes don't match");
  int width = sh.width, height = sh.height;
  int depth = disp.size();

  for (int z=0; z < depth; z++)
  {
	  for (int y = 0; y < height; y++)
	  {
		uchar *di = &disp[z].Pixel(0, y, 0);
		uchar *tdi = &truedisp[z].Pixel(0, y, 0);
		uchar *errmap = &errormap[z].Pixel(0, y, 0);

		for (int x = 0; x < width; x++)
		{
		  uchar *inimp = &interImage[z].Pixel(x, y, 0);
		  if (*inimp!=0)
		  {
			  di[x]=tdi[x];
			  errmap[x]=0;
		  }
		}
	  }
  }

}



void display_CM_avg(double **confMatrixFull, int nD, double* totNDFull, int train)
{


  std::vector<double> dcs;
  dcs.resize(nD);
  double dcs_X, dcs_Y, dcs_XinterY;
  // dcs_X registers segments from algorithm
  // dcs_Y registers groundtruth segmentation

  for (int i=0; i<nD; i++) 
  {
    dcs_X = 0;
    dcs_Y = 0;
    dcs_XinterY = 0;

    dcs_XinterY = confMatrixFull[i][i];

    for (int j=0; j<nD; j++)
      dcs_X += confMatrixFull[i][j];

    for (int j=0; j<nD; j++)
      dcs_Y += confMatrixFull[j][i];

    if (dcs_X + dcs_Y == 0)
    	dcs[i] = 1.0;
    else
    	dcs[i] = 2*dcs_XinterY / (dcs_X + dcs_Y);
  }


  for (int i=0; i<nD; i++) {
    for (int j=0; j<nD; j++) {
    	if (totNDFull[i] == 0)
    		confMatrixFull[i][j] = 0.0;
    	else
    		confMatrixFull[i][j] =  confMatrixFull[i][j]/totNDFull[i];
    }
  }


  DEBUG_OUT0(verbose, debugfile, "\n \t   Full Confusion Matrix \n");

  string dumpstr="\t";

  for (int i = 0; i<nD; i++) {
	  dumpstr += "\t";

	  std::stringstream out;
	  out << i;
	  dumpstr += out.str();
  }
  dumpstr += "\n";
  DEBUG_OUT0(verbose, debugfile, dumpstr.c_str());

  for (int i = 0; i<nD; i++) {
	  dumpstr = "\t";
	  std::stringstream out;
	  out << i;
	  dumpstr += out.str();
	  dumpstr += "\t";
	  for (int j=0; j<nD; j++) {
		  std::stringstream out2;
		  out2 << fixed << setprecision(3) << (double)confMatrixFull[i][j];
		  dumpstr += out2.str();
		  dumpstr += "\t";
	  }
	  dumpstr += "\n";
	  DEBUG_OUT0(verbose, debugfile, dumpstr.c_str());
  }


  DEBUG_OUT0(verbose, debugfile, "\n \t \t Dice Coefficients \n");

  dumpstr = "\t";
  for (int i = 0; i<nD; i++) {
	  std::stringstream out2;
	  out2 << fixed << setprecision(3) << (double)dcs[i] << "\t";
	  dumpstr += out2.str();
  }
  dumpstr += "\n\n";
  DEBUG_OUT0(verbose, debugfile, dumpstr.c_str());


  double avgclassacc = 0.0;

  for (int i = 0; i<nD; i++) {
	  if (totNDFull[i] == 0)
		  avgclassacc += 1.0;
	  else
		  avgclassacc += (double)confMatrixFull[i][i];
  }

  avgclassacc /= nD;

  if (train == 1)
  {
	  DEBUG_OUT1(verbose, debugfile, "\n Average Total Train Class Accuracy = %g \n", avgclassacc);
  } else
  {
	  if (train == 0) // test
	  DEBUG_OUT1(verbose, debugfile, "\n Average Total Test Class Accuracy = %g \n", avgclassacc);
  }
}







void confusionMatrix(CByteImage truedisp, CByteImage disp)
{
	int d, td;
	uchar *di, *tdi;
	int model_emp = 0;
	int model_notEmp = 0;
	int notMod_emp = 0;
	int notMod_notEmp = 0;
	int nState = getNumStates(); //I CAN'T GET THIS TO LINK!
	if(getLocal() == false) nState++;

	CShape sh = disp.Shape();
	if (sh != truedisp.Shape())
		throw CError("image shapes don't match");
	int width = sh.width, height = sh.height;

	 for (int y = 0; y < height; y++) {
			di = &disp.Pixel(0, y, 0);
			tdi = &truedisp.Pixel(0, y, 0);

			for (int x = 0; x < width; x++) {
				d = di[x];
				td = tdi[x];
				if(td == 0 & d == nState - 1)
					model_emp++;
				else if(td != 0 & d == nState - 1)
					model_notEmp++;
				else if(td == 0 & d != nState - 1)
					notMod_emp++;
				else
					notMod_notEmp++;
			}
	 }
	 DEBUG_OUT0(verbose, debugfile, "\n Confusion Matrix \n");
	 DEBUG_OUT0(verbose, debugfile, "Model \t O \t D \t\t O \t D \n");
	 DEBUG_OUT0(verbose, debugfile, "Truth \n");
	 DEBUG_OUT2(verbose, debugfile, "O \t %d \t %d \t\t", model_emp , notMod_emp );
	 DEBUG_OUT2(verbose, debugfile, "%f \t %f \n", model_emp/(float)(notMod_emp+model_emp) , notMod_emp/(float)(notMod_emp+model_emp) );
	 DEBUG_OUT2(verbose, debugfile, "D \t %d \t %d \t\t", model_notEmp , notMod_notEmp );
	 DEBUG_OUT2(verbose, debugfile, "%f \t %f \n", model_notEmp/(float)(model_notEmp+notMod_notEmp) , notMod_notEmp/(float)(model_notEmp+notMod_notEmp) );
}


// computes and prints disparity error statistics, creates errormap
// returns percentages of bad pixels (error > 1) and RMS disp error (both in nonoccluded areas)
// verbose==1: print lots, verbose==0: print single line, verbose<0: print nothing
void evaldispsSet(std::vector<CByteImage> disp, std::vector<CByteImage> truedisp, std::vector<CByteImage> &errormap,
               float &bad1, float &rms, int verbose)
{
	int cnt[ERRLIMIT+1];
	int tot = 0;
	int k;
	float sserr = 0;
	int nStates = getNumStates(); // I CAN'T GET THIS TO LINK! -JW
	for (k = 0; k <= ERRLIMIT; k++)
		cnt[k] = 0;


	for(unsigned int i = 0; i < disp.size(); ++i){
		CShape sh = disp[i].Shape();
		if (sh != truedisp[i].Shape())
			throw CError("image shapes don't match");
		int width = sh.width, height = sh.height;

		errormap[i].ReAllocate(sh);

		for (int y = 0; y < height; y++) {
			uchar *di = &(disp[i]).Pixel(0, y, 0);
			uchar *tdi = &(truedisp[i]).Pixel(0, y, 0);
			uchar *errmap = &(errormap[i]).Pixel(0, y, 0);
			for (int x = 0; x < width; x++) {
				int d = di[x];
				int td = tdi[x];
				if (td == 0 || d == (nStates -1) ) { // unknown true disp
					errmap[x] = UNK;
				} else {
					int err = td - d;
					float ferr = (float)err;
					if (err < 0) err = -err;
					if (err > ERRLIMIT) err = ERRLIMIT;
					tot++;
					cnt[err]++;
					sserr += ferr * ferr;
					errmap[x] = errcol[err];
				}
			}
		}
	}

	if (verbose > 0) {
		printf("err thresh:\t");
		for (k = 0; k < ERRLIMIT; k++)
			printf("%d\t", k);
		printf("\nbad disp %%:\t");
	}

	int ct = 0;
	for (k = 0; k < ERRLIMIT; k++) {
		ct += cnt[k];
		float errpercent = (float)(100.0 * (1.0 - (float)ct/(float)tot));
		if (k == 1)
			bad1 = errpercent;
		if (verbose > 0)
			printf("%.1f\t", errpercent);
		else if (verbose >= 0)
			printf("%d:%5.1f   ", k, errpercent);
	}

	rms = sqrt(sserr / tot);
	if (verbose > 0)
		printf("\nRMS disp error = %.2f\n", rms);
	else if (verbose >= 0)
		printf("rms: %.2f\n", rms);
}

void confusionMatrix(std::vector<CByteImage> &truedisp, std::vector<CByteImage> &disp)
{
	int d, td;
	uchar *di, *tdi;
	int model_emp = 0;
	int model_notEmp = 0;
	int notMod_emp = 0;
	int notMod_notEmp = 0;
	int nState = getNumStates();

	if(getLocal() == false) nState++;


	for(unsigned int i = 0; i < disp.size(); ++i){
		CShape sh = disp[i].Shape();
		if (sh != truedisp[i].Shape())
			throw CError("image shapes don't match");
		int width = sh.width, height = sh.height;

		for (int y = 0; y < height; y++) {
			di = &(disp[i]).Pixel(0, y, 0);
			tdi = &(truedisp[i]).Pixel(0, y, 0);

			for (int x = 0; x < width; x++) {
				d = di[x];
				td = tdi[x];
					notMod_notEmp++; // both are not occluded.
			}
		 }
	}
	 DEBUG_OUT0(verbose, debugfile, "\n Confusion Matrix \n");
	 DEBUG_OUT0(verbose, debugfile, "Model \t O \t D \t\t O \t D \n");
	 DEBUG_OUT0(verbose, debugfile, "Truth \n");
	 DEBUG_OUT2(verbose, debugfile, "O \t %d \t %d \t\t", model_emp , notMod_emp );
	 DEBUG_OUT2(verbose, debugfile, "%f \t %f \n", model_emp/(float)(notMod_emp+model_emp) , notMod_emp/(float)(notMod_emp+model_emp) );
	 DEBUG_OUT2(verbose, debugfile, "D \t %d \t %d \t\t", model_notEmp , notMod_notEmp );
	 DEBUG_OUT2(verbose, debugfile, "%f \t %f \n", model_notEmp/(float)(model_notEmp+notMod_notEmp) , notMod_notEmp/(float)(model_notEmp+notMod_notEmp) );
}
