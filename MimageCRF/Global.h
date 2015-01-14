/*
 * Global.h
 *
 *  Created on: May 2, 2011
 *      Author: bhole
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_

#include<vector>
#include<iostream>
#include "mrf.h"
#include<math.h>
#include<map>
#include <stdlib.h>

class GlobalPairwiseParam
{
public:

	GlobalPairwiseParam();

	std::vector<MRF::CostVal *> horPairGlobalF;
	std::vector<MRF::CostVal *> verPairGlobalF;
	std::vector<MRF::CostVal *> depPairGlobalF;

	// Two nodes that are neighbors are mapped to the spatial/opflow flag and gradient computation.
	// the spatial/opflow flag = 0 for spatial and 1 for opflow
	std::vector <std::map <std::pair <int, int>, std::pair<int, MRF::CostVal> > > edgeGlobalF;

	std::vector<float> thetaVGlobal, thetaZGlobal, thetaCGlobal;

	int nDG;
	int nGG;
	int nGzG;

	int gradVZ;

	int imageWidth, imageHeight, volumeDepth;

	void setOpflowConnection(int opflow)
	{
	  opflowconnection = opflow;
	}

	int getOpflowConnection()
  {
    return opflowconnection;
  }

	void setOpflowVolumeIndex(int opflowvi)
	{
	  opflowvolumeindex = opflowvi;
	}

	int getOpflowVolumeIndex()
	{
	  return opflowvolumeindex;
	}

private:
	int opflowconnection;
	int opflowvolumeindex;

};


MRF::CostVal fnCostQuantizedGlobal(int pix1, int pix2, int di, int dj);
MRF::CostVal fnCostSmoothGlobal(int pix1, int pix2, int di, int dj);


#endif /* GLOBAL_H_ */
