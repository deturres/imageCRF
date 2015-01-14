
#include "project.h"
#include "Image.h"
#include "OpticalFlow.h"

#include <iostream>
#include <stdio.h>

int main(int argc, char* argv[])
{

  // default parameters
	double alpha= 1;
	double ratio=0.5;
	int minWidth= 40;
	int nOuterFPIterations = 3;
	int nInnerFPIterations = 1;
	int nSORIterations= 20; 

  // current settings
  alpha = 0.012;
  ratio = 0.75;
  minWidth = 20;
  nOuterFPIterations = 7;
  nInnerFPIterations = 1;
  nSORIterations = 30;


  opflowns::DImage Im1,Im2;

  Im1.imread("/Users/chetan/Documents/NoBackUp/opencv_samples/im00410.png");
  Im2.imread("/Users/chetan/Documents/NoBackUp/opencv_samples/im00411.png");

  printf("width %d   height %d   nchannels %d\n", Im1.width(),Im1.height(),Im1.nchannels());
  printf("width %d   height %d   nchannels %d\n", Im2.width(),Im2.height(),Im2.nchannels());

  if (Im1.matchDimension(Im2)==false)
  {
    printf("The two images don't match!");
    exit(1);
  }

  opflowns::DImage vx,vy,warpI2;
  opflowns::OpticalFlow::Coarse2FineFlow(vx,vy,warpI2,Im1,Im2,alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations);


  printf("width %d   height %d   nchannels %d\n", vx.width(), vx.height(), vx.nchannels());
  printf("width %d   height %d   nchannels %d\n", vy.width(), vy.height(), vy.nchannels());
  printf("width %d   height %d   nchannels %d\n", warpI2.width(), warpI2.height(), warpI2.nchannels());


  // channels is 1 for vx and vy
  for (int i=10; i < 20; i++)  // height
  {
    for (int j=20; j < 30; j++) // width
    {
      printf("%f ", vx[i*vx.width()+j]);
    }
    printf(";\n");
  }

  printf("\n\n");

  for (int i=10; i < 20; i++)  // height
  {
    for (int j=20; j < 30; j++) // width
    {
      printf("%f ", vy[i*vy.width()+j]);
    }
    printf(";\n");
  }

  // accessing elements.
  /*
  for(int i =  0; i <imHeight; i++)
    for(int j = 0; j<imWidth; j++)
    {
      int offset = (i*imWidth+j)*nChannels;
      for(int k = 0; k<nChannels; k++)
        pData[offset+k];
    }
  */



  printf(" vx min = %f,  vx max = %f \n", vx.immin(), vx.immax());
  printf(" vy min = %f,  vy max = %f \n", vy.immin(), vy.immax());







  return 0;
}
