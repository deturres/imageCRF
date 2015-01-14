/*
 * Global.cpp
 *
 *  Created on: May 2, 2011
 *      Author: bhole
 */

#include "Global.h"

GlobalPairwiseParam globalP;


GlobalPairwiseParam::GlobalPairwiseParam()
{

	horPairGlobalF.resize(0);
	verPairGlobalF.resize(0);
	depPairGlobalF.resize(0);

	edgeGlobalF.resize(0);

	thetaVGlobal.resize(0);
	thetaZGlobal.resize(0);
	thetaCGlobal.resize(0);

	nDG = 0;
	nGG = 0;
	nGzG = 0;

	gradVZ = 0;

	imageWidth = 0;
	imageHeight = 0;
	volumeDepth = 0;

	opflowconnection = 0;


}






// for quantized x-y, x-y-z for x-y and z mode
MRF::CostVal fnCostQuantizedGlobal(int pix1, int pix2, int di, int dj)
{


  int gradval = 0;
  int zflag = 0;

  // gradient cost
  if (globalP.getOpflowConnection() == 0)
  {
    if (pix1 < pix2) {
      if (pix2 - pix1 == 1)
          gradval = globalP.horPairGlobalF[0][pix1];
      else if (pix2 - pix1 == globalP.imageWidth) {
          gradval = globalP.verPairGlobalF[0][pix1];
      } else {
        zflag = 1;
        gradval = globalP.depPairGlobalF[0][pix1];
      }
    } else { // (pix2 < pix1)
      if (pix1 - pix2 == 1)
          gradval = globalP.horPairGlobalF[0][pix2];
      else if (pix1 - pix2 == globalP.imageWidth) {
          gradval = globalP.verPairGlobalF[0][pix2];
      } else {
        zflag = 1;
        gradval = globalP.depPairGlobalF[0][pix2];
      }
    }
  }
  else if (globalP.getOpflowConnection() == 1)
  {
    int pixi=pix1, pixj=pix2;
    if (pix1 < pix2) {}
    else
    {
      pixi=pix2;
      pixj=pix1;
    }

    int volidx = globalP.getOpflowVolumeIndex();
    std::pair <int, int> pv(pixi, pixj);
    std::map <std::pair <int, int>, std::pair <int, MRF::CostVal> >::iterator iter = globalP.edgeGlobalF[volidx].begin();
    iter = globalP.edgeGlobalF[volidx].find(pv);
    if (iter != globalP.edgeGlobalF[volidx].end())
    {
      int edgetype = iter->second.first;
      gradval = (int) iter->second.second;

      if (edgetype == 1)
        zflag = 1;
    }
    else
    {
      std::cout << " This should not be happening \n";
      exit(1);
    }
  }

  int d1V=di, d2V=dj;

  if (d1V>d2V)
  {
    d1V = dj;
    d2V = di;
  }


  if (zflag == 1 && globalP.gradVZ == 1)
	  return 0;

  if (zflag == 1 && globalP.gradVZ == 2)
	  return globalP.thetaVGlobal[(d1V*globalP.nDG + d2V - (d1V*(d1V+1))/2)*globalP.nGG + gradval];

  if (zflag == 1 && globalP.gradVZ == 3)
	  return globalP.thetaZGlobal[(d1V*globalP.nDG + d2V - (d1V*(d1V+1))/2)*globalP.nGzG + gradval];

  if (zflag == 0)
	  return globalP.thetaVGlobal[(d1V*globalP.nDG + d2V - (d1V*(d1V+1))/2)*globalP.nGG + gradval];


}






// for smooth x-y, x-y-z for x-y and z mode
MRF::CostVal fnCostSmoothGlobal(int pix1, int pix2, int di, int dj)
{


  double gradval = 0.0;
  int zflag = 0;

  // gradient cost
  if (globalP.getOpflowConnection() == 0)
  {
    if (pix1 < pix2) {
      if (pix2 - pix1 == 1)
          gradval = exp(globalP.horPairGlobalF[0][pix1]);
      else if (pix2 - pix1 == globalP.imageWidth)
          gradval = exp(globalP.verPairGlobalF[0][pix1]);
      else {
        zflag = 1;
        gradval = exp(globalP.depPairGlobalF[0][pix1]);
      }
    } else
    { // (pix2 < pix1)
      if (pix1 - pix2 == 1)
          gradval = exp(globalP.horPairGlobalF[0][pix2]);
      else if (pix1 - pix2 == globalP.imageWidth)
          gradval = exp(globalP.verPairGlobalF[0][pix2]);
      else {
        zflag = 1;
        gradval = exp(globalP.depPairGlobalF[0][pix2]);
      }
    }
  }
  else if (globalP.getOpflowConnection() == 1)
  {
    int pixi=pix1, pixj=pix2;
    if (pix1 < pix2) {}
    else
    {
      pixi=pix2;
      pixj=pix1;
    }

    int volidx = globalP.getOpflowVolumeIndex();
    std::pair <int, int> pv(pixi, pixj);
    std::map <std::pair <int, int>, std::pair <int, MRF::CostVal> >::iterator iter = globalP.edgeGlobalF[volidx].begin();
    iter = globalP.edgeGlobalF[volidx].find(pv);
    if (iter != globalP.edgeGlobalF[volidx].end())
    {
      int edgetype = iter->second.first;
      gradval = (double) exp(iter->second.second);

      if (edgetype == 1)
        zflag = 1;
    }
    else
    {
      std::cout << " This should not be happening \n";
      exit(1);
    }
  }

  int d1V=di, d2V=dj;

  if (d1V>d2V)
  {
    d1V = dj;
    d2V = di;
  }


  if (zflag == 1 && globalP.gradVZ == 1)
	  return 0; // set during initialization as 0 ie gradval will be 1 because of exp

  if (di == dj)
	  return 0; // indicator function constraint

  if (globalP.gradVZ == 1 || globalP.gradVZ == 2)
    gradval *= globalP.thetaVGlobal[0];

  if (globalP.gradVZ == 3 && zflag == 0)
    gradval *= globalP.thetaVGlobal[0];
  else if (globalP.gradVZ == 3 && zflag == 1)
    gradval *= globalP.thetaZGlobal[0];

  return gradval;


}






