/*
 * crf.cpp
 *
 *  Created on: Mar 15, 2011
 *      Author: bhole
 */

#include "crf.h"

crf::crf(Parameters* prm, std::vector<float> theta) : baseModel(prm, theta)
{
	mrf=0;
  mrfgen=0;
}


crf::~crf()
{

}


void crf::deletemrf(unsigned int index)
{
	// both training and testing use mrf
	delete mrf;
//	if (fpvp->generative==1)
//		delete mrfgen;
}


void crf::allocateDataCostSpace()
{
	// find volume with largest number of slices and allocate that space once

	int numV = gtdirImage.size();
	int maxDepth = 0;

	for (int i=0; i<numV; i++)
	{
		if (gtdirImage[i].size() > maxDepth)
			maxDepth = gtdirImage[i].size();
	}

	// assuming all images have same sizes
	CShape sh = gtdirImage[0][0].Shape();
	int width = sh.width, height = sh.height;

	dCArray.push_back(new MRF::CostVal[maxDepth * width * height * prm->nD]);

	if (fpvp->generative==1)
		dCArrayGen.push_back(new MRF::CostVal[maxDepth * width * height * prm->nD]);


}


void crf::deAllocateDataCostSpace()
{

  for(std::vector<MRF::CostVal *>::reverse_iterator b = dCArray.rbegin(); b != dCArray.rend(); ++b)
	  delete [] *b;

  if (fpvp->generative==1)
  {
	  for(std::vector<MRF::CostVal *>::reverse_iterator b = dCArrayGen.rbegin(); b != dCArrayGen.rend(); ++b)
		  delete [] *b;
  }

}



void crf::allocateSmoothCostGlobalSpace()
{
	// find volume with largest number of slices and allocate that space once

	int numV = gtdirImage.size();
	int maxDepth = 0;

	for (int i=0; i<numV; i++)
	{
		if (gtdirImage[i].size() > maxDepth)
			maxDepth = gtdirImage[i].size();
	}

	// assuming all images have same sizes
	CShape sh = gtdirImage[0][0].Shape();
	int width = sh.width, height = sh.height;

	globalP.horPairGlobalF.push_back(new MRF::CostVal[maxDepth * width * height]);
	globalP.verPairGlobalF.push_back(new MRF::CostVal[maxDepth * width * height]);
	globalP.depPairGlobalF.push_back(new MRF::CostVal[maxDepth * width * height]);


}


void crf::deAllocateSmoothCostGlobalSpace()
{

  for(std::vector<MRF::CostVal *>::reverse_iterator b = globalP.horPairGlobalF.rbegin(); b != globalP.horPairGlobalF.rend(); ++b)
	  delete [] *b;

  for(std::vector<MRF::CostVal *>::reverse_iterator b = globalP.verPairGlobalF.rbegin(); b != globalP.verPairGlobalF.rend(); ++b)
	  delete [] *b;

  for(std::vector<MRF::CostVal *>::reverse_iterator b = globalP.depPairGlobalF.rbegin(); b != globalP.depPairGlobalF.rend(); ++b)
	  delete [] *b;

}






DataCost* crf::computeDataCost(int genparam, unsigned int index)
{

	MRF::CostVal badcost;

	badcost = getbadCost(genparam);

	// dCArray and dCArrayGen will be used to store the probability values in logistic regression instead
	// will need to see how to process dsiArrayVGen (derive equations etc)

	if (genparam==0)
	{
		if (dCArray.empty() == true)
			throw CError("call allocateDataCost first");
	} else if (genparam==1) {
		if (dCArrayGen.empty() == true)
			throw CError("call allocateDataCost first");
	}

	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	int zslice = 0;

	if (timeseq == 0 && (fpvp->loc != 0 || fpvp->klrpv.locKlr ==2))
		zslice = startSliceNo[index];


	// dirDisp.resize(depth);
	errormap[index].resize(depth);

	int dsiIndex = 0; // for crf , we use the full volum

	for(unsigned int i = 0; i < gtdirImage[index].size(); ++i)
	{

		// dsiIndex = 0; // since we are processing each image and not storing entire volume

		sh.nBands = 1;
		// dirDisp[i].ReAllocate(sh);
		errormap[index][i].ReAllocate(sh);

		matrix<DBL_TYPE> hogdescs, appdescs;

		if (fpvp->klrpv.klr==1)
			readklrdescfiles(hogdescs, appdescs, index, i+zslice);


		for (int y = 0; y < height; y++)
		{
			uchar *WTArow = &dirDisp[index][i].Pixel(0, y, 0);

			for (int x = 0; x < width; x++)
			{

				int numbest = 0;
				MRF::CostVal bestval = badcost;
				ublas::vector<DBL_TYPE> f;

				std::vector<double> dsiValueSum(prm->nD);

				if (fpvp->klrpv.klr==1)  // this may work only for intensity for timeseq=1 cases
				{
					ublas::vector<DBL_TYPE> Xtest(Xtrain.size2());
					getDataCostKlr(f, Xtest, index, i, x, y, zslice, hogdescs, appdescs, width, height);
				}


				for (int d = 0; d < prm->nD; d++)
				{

					MRF::CostVal dsiValue = 0;

					if (fpvp->klrpv.klr==1)
					{
					  dsiValue += -f(d);
					  dsiValue += fpvp->klrpv.biasklr;
					}

					dsiValue += getDataCost(genparam, width, height, depth, index, i , y, x, d, zslice);

					// The cost of pixel p and label l is stored at dCArray[p*nLabels+l]
					// Since in logistic regression, I am not using E_d ie e(-E_d), I negate the values of dCValue so i can work on prob()
					// directly.
					if (genparam==0)
						(dCArray[0])[dsiIndex++] =  (float)dsiValue;
					else if (genparam==1)
						(dCArrayGen[0])[dsiIndex++] =  (float)dsiValue;

				}

			}
		}

	}

	if (genparam==0)
		return new DataCost ((MRF::CostVal*)dCArray[0]);
	else if (genparam==1)
		return new DataCost ((MRF::CostVal*)dCArrayGen[0]);

}




SmoothnessCost* crf::computeSmoothnessCost(int index)
{

	// set global variables needed by MRF
	// - nDG, nGG, nGzG
	globalP.nDG = prm->nD;
	globalP.nGG = fpvp->gradThreshVec.size();
	globalP.nGzG = fpvp->gradThreshVecZ.size();

	// - globalP.gradVZ
	globalP.gradVZ = fpvp->gradVZ;

	// image properties for this volume
	globalP.volumeDepth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	globalP.imageWidth = sh.width;
	globalP.imageHeight = sh.height;

	int depth = globalP.volumeDepth;
	int width = globalP.imageWidth;
	int height = globalP.imageHeight;

	// - theta values
	globalP.thetaVGlobal.erase(globalP.thetaVGlobal.begin(), globalP.thetaVGlobal.end());
	globalP.thetaZGlobal.erase(globalP.thetaZGlobal.begin(), globalP.thetaZGlobal.end());
	globalP.thetaCGlobal.erase(globalP.thetaCGlobal.begin(), globalP.thetaCGlobal.end());


	int nK = (int)fpvp->thetaV.size();
	int nZ = (int)fpvp->thetaZ.size();
	int nC = (int)fpvp->thetaC.size();

	for (int k = 0; k < nK; k++)
		globalP.thetaVGlobal.push_back(fpvp->thetaV[k]);

	for (int k = 0; k < nZ; k++)
		globalP.thetaZGlobal.push_back(fpvp->thetaZ[k]);

	for (int k = 0; k < nC; k++)
		globalP.thetaCGlobal.push_back(fpvp->thetaC[k]);


	if (globalP.getOpflowConnection() == 1)
	{
	  // all the processing has occured in the initialization phase
	  // and so nothing else to do.
	} else {

    // - horPairGlobalF, verPairGlobalF, depPairGlobalF
    if (fpvp->csgrad == 0)
    {
      int n = 0;
      for (int z = 0; z < depth; z++)
      {
        for (int y = 0; y < height; y++)
        {
          for (int x = 0; x < width; x++)
          {
            uchar *grad   = &dirgrad[index][z].Pixel(x, y, 0);

            globalP.horPairGlobalF[0][n] = grad[0];
            globalP.verPairGlobalF[0][n] = grad[1];
            globalP.depPairGlobalF[0][n] = grad[2];

            n++;
          }
         }
       }
    }
    else if (fpvp->csgrad == 1)
    {

      CvScalar s;

      int n = 0;
      for (int z = 0; z < depth; z++)
      {
        for (int y = 0; y < height; y++)
        {
          for (int x = 0; x < width; x++)
          {
            s=cvGet2D(dircsgrad[index][z],y,x); // get the (i,j) pixel value
            globalP.horPairGlobalF[0][n] = s.val[0];
            globalP.verPairGlobalF[0][n] = s.val[1];
            globalP.depPairGlobalF[0][n] = s.val[2];

            n++;
          }
        }
      }
    }

	}

	// return function pointer
	if (fpvp->csgrad == 0)
		return new SmoothnessCost((MRF::SmoothCostGeneralFn) fnCostQuantizedGlobal);
	else if (fpvp->csgrad == 1)
		return new SmoothnessCost((MRF::SmoothCostGeneralFn) fnCostSmoothGlobal);

}


void crf::regularizeLL()
{
	// no need to change loglikelihood as it is not calculated correctly anyway
}


int checkNeighborPixel(int xp, int yp, int zp, int xq, int yq, int zq, int width, int height, int depth)
{
  if (xp==0 && xq<0)
    return -1;

  if (yp==0 && yq<0)
    return -1;

  if (zp==0 && zq<0)
    return -1;

  if (xp==width-1 && xq>width-1)
    return -1;

  if (yp==height-1 && yq>height-1)
    return -1;

  if (zp==depth-1 && zq>depth-1)
    return -1;

  return 0;
}


// using the 3D graph structure traverse thru all points.
// if the point is interactively labeled and all its neighbors are interactively labeled, don't update dcost
// if point is interactively labeled and its neighbor is not interactively labeled, update dcost of neighbor
// as dcost(d) + vcost(const, d) where const is the label of the interactively labeled point
// and d is the variable indexing possible labels for a point
void crf::updateDataCostDueToInteractive(int index, DataCost *dcost, SmoothnessCost *scost)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (globalP.getOpflowConnection() == 1 && timeseq == 0)
    assert(false);

  if (globalP.getOpflowConnection() == 1 && timeseq == 1)
  {
    std::map <std::pair <int, int>, std::pair <int, MRF::CostVal> >::iterator iter = globalP.edgeGlobalF[index].begin();
    int nD = prm->nD;

    for (; iter != globalP.edgeGlobalF[index].end(); iter++)
    {
      int node1 = iter->first.first;
      int node2 = iter->first.second;

      int z1 = node1 / (width*height);
      int z2 = node2 / (width*height);

      int rem1 = node1 - z1*(width*height);
      int rem2 = node2 - z2*(width*height);

      int y1 = rem1 / width;
      int y2 = rem2 / width;

      int x1 = rem1 - y1*width;
      int x2 = rem2 - y2*width;

      uchar *tdip = &gtdirImage[index][z1].Pixel(x1, y1, 0);
      uchar *tdiq = &gtdirImage[index][z2].Pixel(x2, y2, 0);

      uchar *inimp = &interdirImage[index][z1].Pixel(x1, y1, 0);
      uchar *inimq = &interdirImage[index][z2].Pixel(x2, y2, 0);

      if (*inimp == 255 && *inimq == 0) // here p is interactive, q is not so add to q
      {
        for (int d=0; d<prm->nD; d++)
        {
          dcost->updateDataCostRaw(node2*nD + d, dcost->getDataCostRaw(node2*nD + d) + scost->getSmoothnessCostRaw(node1, node2, tdip[0], d));
        }
      }

      if (*inimp == 0 && *inimq == 255) // here p is not interactive, q is so add to p
      {
        for (int d=0; d<prm->nD; d++)
        {
          dcost->updateDataCostRaw(node1*nD + d, dcost->getDataCostRaw(node1*nD + d) + scost->getSmoothnessCostRaw(node1, node2, d, tdiq[0]));
        }
      }
    }
  } else
  {
    int nD = prm->nD;
    int pIndex = 0;
    for(int z = 0; z < depth; z++)
    {
      for (int y = 0; y < height; y++)
      {
        uchar *tdi = &gtdirImage[index][z].Pixel(0, y, 0);
        uchar *inim = &interdirImage[index][z].Pixel(0, y, 0);

        for (int x = 0; x < width; x++)
        {
          if (inim[x]==0)
          {
            // skip because it (p) is a non-interactive pixel and cannot add smoothcosts to another pixel
          } else
          {
            // try all neighbors q of p(x,y,z) for the current interactive pixel
            // q(x+1, y, z)
            if (checkNeighborPixel(x, y, z, x+1, y, z, width, height, depth)== 0)
            { // fine
              int qIndex = pIndex + 1;
              if (inim[x+1] == 255)
              {
                // skip because it is labeled too ie q is interactive
              } else
              {
                // q is not interactive, so update q's datacost
                for (int d=0; d<nD; d++)
                {
                  dcost->updateDataCostRaw(qIndex*nD + d, dcost->getDataCostRaw(qIndex*nD + d) + scost->getSmoothnessCostRaw(pIndex, qIndex, tdi[x], d));
                }
              }
            }

            // q(x-1, y, z)
            if (checkNeighborPixel(x, y, z, x-1, y, z, width, height, depth)== 0)
            { // fine
              int qIndex = pIndex - 1;
              if (inim[x-1] == 255)
              {
                // skip because it is labeled too
              } else
              {
                for (int d=0; d<nD; d++)
                {
                  dcost->updateDataCostRaw(qIndex*nD + d, dcost->getDataCostRaw(qIndex*nD + d) + scost->getSmoothnessCostRaw(pIndex, qIndex, tdi[x], d));
                }
              }
            }


            // q(x, y+1, z)
            if (checkNeighborPixel(x, y, z, x, y+1, z, width, height, depth)== 0)
            { // fine

              int qIndex = pIndex + width;
              uchar *inimy = &interdirImage[index][z].Pixel(0, y+1, 0);
              if (inimy[x] == 255)
              {
                // skip because it is labeled too
              } else
              {
                for (int d=0; d<nD; d++)
                {
                  dcost->updateDataCostRaw(qIndex*nD + d, dcost->getDataCostRaw(qIndex*nD + d) + scost->getSmoothnessCostRaw(pIndex, qIndex, tdi[x], d));
                }
              }
            }

            // q(x, y-1, z)
            if (checkNeighborPixel(x, y, z, x, y-1, z, width, height, depth)== 0)
            { // fine

              int qIndex = pIndex - width;
              uchar *inimy = &interdirImage[index][z].Pixel(0, y-1, 0);
              if (inimy[x] == 255)
              {
                // skip because it is labeled too
              } else
              {
                for (int d=0; d<nD; d++)
                {
                  dcost->updateDataCostRaw(qIndex*nD + d, dcost->getDataCostRaw(qIndex*nD + d) + scost->getSmoothnessCostRaw(pIndex, qIndex, tdi[x], d));
                }
              }
            }

            // q(x, y, z+1)
            if (checkNeighborPixel(x, y, z, x, y, z+1, width, height, depth)== 0)
            { // fine

              int qIndex = pIndex + width*height;
              uchar *inimz = &interdirImage[index][z+1].Pixel(0, y, 0);
              if (inimz[x] == 255)
              {
                // skip because it is labeled too
              } else
              {
                for (int d=0; d<nD; d++)
                {
                  dcost->updateDataCostRaw(qIndex*nD + d, dcost->getDataCostRaw(qIndex*nD + d) + scost->getSmoothnessCostRaw(pIndex, qIndex, tdi[x], d));
                }
              }
            }

            // q(x, y, z-1)
            if (checkNeighborPixel(x, y, z, x, y, z-1, width, height, depth)== 0)
            { // fine

              int qIndex = pIndex - width*height;
              uchar *inimz = &interdirImage[index][z-1].Pixel(0, y, 0);
              if (inimz[x] == 255)
              {
                // skip because it is labeled too
              } else
              {
                for (int d=0; d<nD; d++)
                {
                  dcost->updateDataCostRaw(qIndex*nD + d, dcost->getDataCostRaw(qIndex*nD + d) + scost->getSmoothnessCostRaw(pIndex, qIndex, tdi[x], d));
                }
              }
            }
          }
          pIndex++;
        }
      }
    }

  }

}

// IMPORTANT NOTE: when we have a hidden-crf, the interactive flag is no use for a training example
// but the interactive images need to be in place because they are used by the method.
// the unknown and interactive work together in the hidden crf
//
// for training example
//    temporarily set prm->interactive = 1
//    call modelProcess
//    call empirical dist + hidden dist update (make sure picks all pixels)
//    delete mrf
//    temporarily set prm->interactive = 0
//    call modelProcess
//    evaluate results
//    outputEvaluations
//    call model dist (makes sure picks all pixels)
//    delete mrf
void crf::hiddenModelProcess(unsigned int index, int iter, int iterout)
{
  if (prm->crfpv.crfp == 1 && prm->crfpv.hidden == 1 && prm->unknown == 0)
    assert(0); // for a hidden crf, unknown has to be 1
  assert(interdirImage[index].size() > 0);

  int temp_interactive = prm->interactive;

  // using the crf to generate p(h|y,x)
  prm->interactive = 1;
  modelProcess(index);  // call model etc

  if (prm->graddescpv.bfgs_flag == 0)
  { // dump images for gradient descent methods
    outputImageResults(index, iterout*1000 + iter, "hidd"); // dump output images
  }

  computeEmpirDist(index);
  computeHiddenDist(index); // note that this updates the empir_dist variables
  deletemrf(index);

  // using the crf to generate p(h,y|x)
  prm->interactive = 0;
  modelProcess(index);  // call model etc
  evaluateResults(index); // for calling evaldisps
  outputEvalutions(index); // dump evaluations

  if (prm->graddescpv.bfgs_flag == 0)
  { // dump images for gradient descent methods
    outputImageResults(index, iterout*1000 + iter); // dump output images
  }

  computeModelDist(index);
  deletemrf(index);

  prm->interactive = temp_interactive;
}




void crf::modelProcess(unsigned int index)
{

	// image properties for this volume
	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width;
	int height = sh.height;

  DataCost *dcost;
  SmoothnessCost *scost;

  // datacost
  if (prm->featurepv.learngradonly==1)
  {
    std::cout << " Possibly not working !!! \n";
    exit(1);
    // dcost = computeDataCost(0, i); // the function in crfmodel.cpp
  } else
  {
    dcost = computeDataCost(0, index);
  }

  // smoothcost
  if (fpvp->context==0 && fpvp->gradContext==0)
  {
    std::cout << " Not defined (context or gradContext should be defined) \n";
    exit(1);
  }

  if (fpvp->context==1 && fpvp->gradContext==1)
  {
    std::cout << " Both context and gradContext should not be set \n";
    exit(1);
  }

  if (fpvp->context==1)
  {
      std::cout << " No longer functional - will need to use thetaC inplace of thetaV in code ...  \n";
      exit(1);
  }

  if (fpvp->gradContext==1)
  {
    scost = computeSmoothnessCost(index);
  }

  // assuming triangle of pairwise conditions
  if ((prm->iopv.testDirIndexV[index]==1 && prm->crfpv.inferencerTest == USE_GC) || (prm->iopv.testDirIndexV[index]==0 && prm->crfpv.inferencer == USE_GC))
  {
    // make sure submodularity is valid
    int cRC = 1;
    if (prm->featurepv.nV > 1)
      cRC = checkRegularityCondition(prm->featurepv.thetaV, prm->featurepv.gradThreshVec, prm->nD);

    int cRC2 = 1;
    if (prm->featurepv.nZ > 1)
      cRC2 = checkRegularityCondition(prm->featurepv.thetaZ, prm->featurepv.gradThreshVecZ, prm->nD);

    if (cRC == 0 || cRC2 == 0)
    {
      std::cout << " cRC failure \n";
      exit(1);
    }
  }

  if (prm->interactive == 1)
  {
    updateDataCostDueToInteractive(index, dcost, scost);
  }

  int infermac;
  int testv = 0;
  if (prm->iopv.testDirIndexV[index]==1)
  {
    infermac = prm->crfpv.inferencerTest;
    testv = 1;
  }
  if (prm->iopv.testDirIndexV[index]==0)
  {
    infermac = prm->crfpv.inferencer;
    testv = 0;
  }

  char logname[500];
  sprintf(logname, "%slog.txt", prm->iopv.outstem.c_str());
  std::ofstream logStream(logname);

	mrf=crfmodel(depth, width, height, prm->nD, dcost, scost, dirDisp[index], prm->graddescpv.closeEnoughPercent, infermac, 1, interdirImage[index], gtdirImage[index], prm->interactive, &logStream, testv);
	// we have a different mrf for each patient. Need to make sure we delete it before we allocate new one for new volume

	if (fpvp->generative==1)
	{
		std::cout << " generative in crf - not coded (look in older code) \n";
		exit(1);
	}

  delete scost;
  delete dcost;
}

void crf::updateDist(int index, int fiterflag)
{

	// if 1st iteration then do for empirical otherwise if not first, don't do for empirical
	if (fiterflag==0)
		computeEmpirDist(index);

	// do modelDist for all
	computeModelDist(index);

}


void crf::computeEmpirDist(int index)
{

  if (prm->optClassAccuracy == 1 || prm->optClassAccuracy == 2)
    assert(0); // error crf does not work with optimize class weighting

	if (fpvp->nU > 0 && fpvp->intensity==1)
		computeEmpirDistU(index);

	if (fpvp->nA > 0 && (fpvp->app==1 || fpvp->app==2))
	  computeEmpirDistA(index);

	if (fpvp->nH > 0 && fpvp->hog==1)
	  computeEmpirDistH(index);

	if (fpvp->nL > 0 && (fpvp->loc==1 || fpvp->loc==6 || fpvp->loc==8 || fpvp->loc==7 || fpvp->loc==9))
	  computeEmpirDistL(index);

	if (fpvp->nM > 0 && fpvp->mot==1)
	  computeEmpirDistM(index);

	if (fpvp->klrpv.klr==1)
	{
		computeEmpirDistW(index);
	}

	if (fpvp->nO > 0 && fpvp->opflowlocal==1)
	  computeEmpirDistO(index);

	// when nV == 1, we are not learning parameters (because constrast sensitive smoothing)
	// but in the function we still need to check again because nZ could be 1 , but not nV and so on
	if (fpvp->nV > 0 || fpvp->nC > 0 || fpvp-> nZ > 0)
		computeEmpirDistV(index);

	if (fpvp->nB > 0 && fpvp->bias==1)
	  computeEmpirDistB(index);

	if (fpvp->nE > 0 && fpvp->volume==1)
	  computeEmpirDistE(index);


}



void crf::computeModelDist(int index)
{

	int mpe = -1;

    if (prm->iopv.testDirIndexV[index]==1 && (prm->crfpv.inferencerTest == USE_GC || prm->crfpv.inferencerTest == USE_MP))
    {
    	mpe = 1;
    }

    if (prm->iopv.testDirIndexV[index]==0 && (prm->crfpv.inferencer == USE_GC || prm->crfpv.inferencer == USE_MP))
    {
    	mpe = 1;
    }

    if (prm->iopv.testDirIndexV[index]==1 && (prm->crfpv.inferencerTest == USE_MF ))
    {
    	mpe = 0;
    }

    if (prm->iopv.testDirIndexV[index]==0 && (prm->crfpv.inferencer == USE_MF) )
    {
    	mpe = 0;
    }


    if (mpe == -1)
    {
    	std::cout << " Need to initialize variable mpe" << std::endl;
    	exit(1);
    }





    // this uses point estimates
    if (mpe == 1)
    {

		if (fpvp->nU > 0 && fpvp->intensity==1)
			computeModelDistUmpe(index);

		if (fpvp->nA > 0 && (fpvp->app==1 || fpvp->app==2))
		  computeModelDistAmpe(index);

		if (fpvp->nH > 0 && fpvp->hog==1)
		  computeModelDistHmpe(index);

		if (fpvp->nL > 0 && (fpvp->loc==1 || fpvp->loc==6 || fpvp->loc==8 || fpvp->loc==7 || fpvp->loc==9))
		  computeModelDistLmpe(index);

		if (fpvp->nM > 0 && fpvp->mot==1)
		  computeModelDistMmpe(index);

		if (fpvp->klrpv.klr==1)
			computeModelDistWmpe(index);

	  if (fpvp->nO > 0 && fpvp->opflowlocal==1)
	    computeModelDistOmpe(index);

		// when nV == 1, we are not learning parameters (because constrast sensitive smoothing)
		// but in the function we still need to check again because nZ could be 1 , but not nV and so on
		if (fpvp->nV > 0 || fpvp->nC > 0 || fpvp-> nZ > 0)
			computeModelDistVmpe(index);

	  if (fpvp->nB > 0 && fpvp->bias==1)
	    computeModelDistBmpe(index);

	  if (fpvp->nE > 0 && fpvp->volume==1)
	      computeModelDistEmpe(index);

    }

    // this uses mrf probability values
	if (mpe == 0)
	{
		if (fpvp->nU > 0 && fpvp->intensity==1)
			computeModelDistU(index);

		if (fpvp->nA > 0 && (fpvp->app==1 || fpvp->app==2))
		  computeModelDistA(index);

		if (fpvp->nH > 0 && fpvp->hog==1)
		  computeModelDistH(index);

		if (fpvp->nL > 0 && (fpvp->loc==1 || fpvp->loc==6 || fpvp->loc==8 || fpvp->loc==7 || fpvp->loc==9))
		  computeModelDistL(index);

		if (fpvp->nM > 0 && fpvp->mot==1)
		  computeModelDistM(index);

		if (fpvp->klrpv.klr==1)
			computeModelDistW(index);

	  if (fpvp->nO > 0 && fpvp->opflowlocal==1)
	    computeModelDistO(index);

		// when nV == 1, we are not learning parameters (because constrast sensitive smoothing)
		// but in the function we still need to check again because nZ could be 1 , but not nV and so on
		if (fpvp->nV > 1 || fpvp->nC > 0 || fpvp-> nZ > 1)
			computeModelDistV(index);

	  if (fpvp->nB > 0 && fpvp->bias==1)
	    computeModelDistB(index);

	  if (fpvp->nE > 0 && fpvp->volume==1)
	      computeModelDistE(index);

	}

}


// grayscale and // RGB or Yuv or Lab
// does not current have interactive mode
void crf::computeEmpirDistU(int index)
{
	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nU == 0)
		return; // nothing to do

	int nbins = double(fpvp->nU)/prm->nD;
	double binner = fpvp->rangeI/nbins;

	int thetaid;
  int thetaidrgb[3]={0,0,0};
  int thetaidYuv[2]={0,0};


	uchar *pix1t = 0;
	unsigned short pix1 = 0;

	for (int z = 0; z < depth; z++)
	{
		for (int y = 0; y < height; y++)
		{
			uchar *di = &gtdirImage[index][z].Pixel(0, y, 0);

			for (int x = 0; x < width; x++)
			{
	      // unknown pixels do not count during learning parameters
	      if (prm->unknown == 1 && di[x]==prm->nD)
	        continue;

				if (prm->timeseq==1)
				{
					pix1t = &indirImaget[index][z].Pixel(x, y, 0);
			        // getCostIntensityDiscBins(pix1t, nbins, di[x], thetaU, thetaidrgb);
					if (prm->featurepv.onlyUV == 1)
			      getCostIntensityDiscBinsuv(pix1t, nbins, di[x], fpvp->thetaU, thetaidYuv);
					else
					  getCostIntensityDiscBinsYuv(pix1t, nbins, di[x], fpvp->thetaU, thetaidYuv);

					if (prm->featurepv.onlyUV == 0)
			      totalEmpirDistU[thetaidYuv[0]]++;
	        totalEmpirDistU[thetaidYuv[1]]++;

				}
				else if (prm->timeseq==0)
				{
					pix1 = indirImage[index][z].Pixel(x, y, 0);
					getCostIntensityDiscBins(pix1, nbins, di[x], fpvp->thetaU, thetaid);
					totalEmpirDistU[thetaid] ++;
				}

				// for us the costs are the theta param, but in differentiation, they vanish,
				// so its just whether the feature is enabled or not.
				// move costs pointer to next pixel


			}
		}
	}
}


void crf::computeEmpirDistA(int index)
{

	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nA == 0)
		return; // nothing to do

	int nT = fpvp->nA/prm->nD; //this gives centers per class
	// this will work for app==2 and app==1 because in app==2, we have nD*the number of total clusters as the number of parameters
	// while in app==1 it is nD*number of clusters of that class (only that clusters are also = number of params)

  int thetaid;


	for (int z=0; z < depth; z++)
	{
		for (int y = 0; y < height; y++)
		{
			uchar *di = &gtdirImage[index][z].Pixel(0, y, 0);

			for (int x = 0; x < width; x++)
			{
	      // unknown pixels do not count during learning parameters
	      if (prm->unknown == 1 && di[x]==prm->nD)
	        continue;

				int d = di[x];
				int dm = d;
        if (fpvp->app==2)
          dm=0;

            	if (((int)appdirImage[index][z][dm].Pixel(x, y, 0) >= nT && fpvp->app==1) || ((int)appdirImage[index][z][0].Pixel(x, y, 0) >= fpvp->nA && fpvp->app==2))
            	{
            		continue;
            	}
                getCostAppDiscPatch((int)appdirImage[index][z][dm].Pixel(x, y, 0), appclass[dm]->getPatchSize(), x, y, width, height, fpvp->thetaA, d, nT, thetaid);

                // assume the patch size is the same for all classes
                int patchSize = appclass[dm]->getPatchSize();

                if (x>(int)patchSize/2 && x<width-patchSize/2 && y>(int)patchSize/2 && y<height-patchSize/2)
                	totalEmpirDistA[thetaid] ++;

			}
		}
	}
}




void crf::computeEmpirDistH(int index)
{

	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nH == 0)
		return; // nothing to do

	int nHpC = fpvp->nH/prm->nD; // HoG vocabulary size

    int thetaid;

    for (int z = 0; z < depth; z++)
    {
    	for (int y = 0; y < height; y++)
    	{
    		uchar *di = &gtdirImage[index][z].Pixel(0, y, 0);

    		for (int x = 0; x < width; x++)
    		{
          // unknown pixels do not count during learning parameters
          if (prm->unknown == 1 && di[x]==prm->nD)
            continue;

    			unsigned short hogval = hogdirImage[index][z].Pixel(x, y, 0);
    			getCostHoGDiscBins((int) hogval, di[x], nHpC, fpvp->thetaH, thetaid);
    			//getCostHoGDiscBinsTemp((int) hogval, di[x], fpvp->thetaH, thetaid);

    			// for us the costs are the theta param, but in differentiation, they vanish,
    			// so its just whether the feature is enabled or not.
    			// move costs pointer to next pixel
    			if (thetaid >= 0)
    				totalEmpirDistH[thetaid] ++;

    		}
    	}
    }
}



void crf::computeEmpirDistM(int index)
{
	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nM == 0)
		return; // nothing to do

	int nMpC = fpvp->nM/prm->nD; // Mot vocabulary size

    int thetaid;

    for (int z = 0; z < depth; z++)
    {
    	for (int y = 0; y < height; y++)
    	{
    		uchar *di = &gtdirImage[index][z].Pixel(0, y, 0);

    		for (int x = 0; x < width; x++)
    		{
          // unknown pixels do not count during learning parameters
          if (prm->unknown == 1 && di[x]==prm->nD)
            continue;

    			unsigned short motval = motdirImage[index][z].Pixel(x, y, 0);
    			getCostmotDiscBins((int) motval, di[x], nMpC, fpvp->thetaM, thetaid);

    			totalEmpirDistM[thetaid] ++;

    		}
    	}
    }
}


void crf::computeEmpirDistB(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nB == 0)
    return; // nothing to do

  int thetaid;

  for (int z = 0; z < depth; z++)
  {
      for (int y = 0; y < height; y++)
      {
        uchar *di = &gtdirImage[index][z].Pixel(0, y, 0);

        for (int x = 0; x < width; x++)
        {
          // unknown pixels do not count during learning parameters
          if (prm->unknown == 1 && di[x]==prm->nD)
            continue;

          getCostBias(di[x], fpvp->thetaB, thetaid);

          // for us the costs are the theta param, but in differentiation, they vanish,
          // so its just whether the feature is enabled or not.
          // move costs pointer to next pixel
          if (thetaid >= 0)
          {
          if (prm->optClassAccuracy == 1)
          {
            uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
            totalEmpirDistB[thetaid] += 1.0/classSize[gtpix[0]];
          } else
            totalEmpirDistB[thetaid] ++;
          }
        }
      }
  }
}




void crf::computeEmpirDistE(int index)
{

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nE == 0)
    return; // nothing to do

  int thetaid;

  for (int z = 0; z < depth; z++)
  {
      for (int y = 0; y < height; y++)
      {
        uchar *di = &gtdirImage[index][z].Pixel(0, y, 0);

        for (int x = 0; x < width; x++)
        {
          // unknown pixels do not count during learning parameters
          if (prm->unknown == 1 && di[x]==prm->nD)
            continue;

          getCostInverseClassSize(di[x], inverseClassSize, fpvp->thetaE, thetaid);

          // for us the costs are the theta param, but in differentiation, they vanish,
          // so its just whether the feature is enabled or not.
          // move costs pointer to next pixel
          if (thetaid >= 0)
          {
          if (prm->optClassAccuracy == 1)
          {
            uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
            totalEmpirDistE[thetaid] += 1.0/classSize[gtpix[0]];
          } else
            totalEmpirDistE[thetaid] += inverseClassSize[di[x]];
          }

        }
      }
  }
}





void crf::computeEmpirDistL(int index)
{

	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nL == 0)
		return; // nothing to do

	int nLpC = fpvp->nL/prm->nD;  // this variable will be useful only for cube setting (loc==1) or loc==9

    int thetaid;

    int zslice = 0;
	if (timeseq == 0 && (fpvp->loc != 0 || fpvp->klrpv.locKlr ==2))
		zslice = startSliceNo[index];

    for (int z = 0; z < depth; z++)
    {

    	for (int y = 0; y < height; y++)
    	{
    		uchar *di = &gtdirImage[index][z].Pixel(0, y, 0);

    		for (int x = 0; x < width; x++)
    		{
          // unknown pixels do not count during learning parameters
          if (prm->unknown == 1 && di[x]==prm->nD)
            continue;

    			if (fpvp->loc==1)
    			{
    				/*
    				unsigned short locpix = locdirImage[index][di[x]][z].Pixel(x,y,0);
    				getCostLocationDiscCube(locpix, x, y, z, fpvp->loccubex, fpvp->loccubey, fpvp->loccubez, numTrainingPats, nLpC, fpvp->thetaL, width, height, depth, di[x], thetaid);
    				double locpixd = (double)locpix/ (fpvp->loccubex * fpvp->loccubey * fpvp->loccubez * numTrainingPats);

    				totalEmpirDistL[thetaid] += (1-locpixd);
    				*/

    				std::cout << " Need to fix this - location " << std::endl;
    				exit(1);


    			} else if (fpvp->loc==6 || fpvp->loc==8) { // CRF hard assignment and soft

    				getCostLocationDiscGaussian(LocationMVclass, z, zslice, x, y, fpvp->thetaL, di[x], fpvp->loc, fpvp->generative, thetaid);
    				totalEmpirDistL[thetaid]++;

    			} else if (fpvp->loc==9 || (fpvp->loc==7 && timeseq == 1)) {  //CRF soft or hard assignment redefined

    				getCostLocationDiscGaussian2(LocationMVclass, z, zslice, x, y, fpvp->thetaL, di[x], fpvp->loc, fpvp->generative, prm->nD, thetaid);
    				totalEmpirDistL[thetaid]++;

    			} else if (fpvp->loc==7 && timeseq == 0) {  //CRF hard assignment redefined

    				getCostLocationDiscGaussian3(index, z, x, y, di[x], fpvp->generative, prm->nD, thetaid);
    				totalEmpirDistL[thetaid]++;
    			}


    		}
    	}
    }
}


void crf::computeEmpirDistW(int index)
{
	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->klrpv.nW == 0)
		return; // nothing to do

	int dimW = wparam.size2();
	int Ntr = wparam.size1();


	int zslice = 0;

	if (timeseq == 0 && (fpvp->loc != 0 || fpvp->klrpv.locKlr ==2))
		zslice = startSliceNo[index];


	uchar *pix1t = 0;
	unsigned short pix1 = 0;

	ublas::vector<DBL_TYPE> Xtest(Xtrain.size2());
	ublas::vector<DBL_TYPE> f;


	for (int z = 0; z < depth; z++)
	{

		matrix<DBL_TYPE> hogdescs, appdescs;

		if (fpvp->klrpv.klr==1)
			readklrdescfiles(hogdescs, appdescs, index, z+zslice);

		for (int y = 0; y < height; y++)
		{
			uchar *di = &gtdirImage[index][z].Pixel(0, y, 0);

			for (int x = 0; x < width; x++)
			{
        // unknown pixels do not count during learning parameters
        if (prm->unknown == 1 && di[x]==prm->nD)
          continue;

				getDataCostKlr(f, Xtest, index, z, x, y, zslice, hogdescs, appdescs, width, height);

				for (int j=0; j<Ntr; j++)
				{
					double Kval = ivm->classify_klrexp(&Xtest, j);

					int thetaid = j*dimW + di[x];
					totalEmpirDistW[thetaid] += Kval;
				}

			}
		}
	}
}


// optical flow
void crf::computeEmpirDistO(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nO == 0)
    return; // nothing to do

  int nOpC = 0;
  if (fpvp->opflowlocal == 1)
  {
    nOpC = fpvp->nO/(prm->nD * fpvp->opflowlocalframes.size()); // there are (u,v) variables per frame except last frame for continuous features
  // but is more complicated for bin features (modeled similar to color uv.
  }

  int thetaid;
  int thetaido[2]={0,0};
  float indcost[2] = {0.0, 0.0};

  // the last frame should not be counted towards the cost
  for (int z = 0; z < depth-1; z++)
  {
    for (int y = 0; y < height; y++)
    {
      uchar *di = &gtdirImage[index][z].Pixel(0, y, 0);

      for (int x = 0; x < width; x++)
      {
        // unknown pixels do not count during learning parameters
        if (prm->unknown == 1 && di[x]==prm->nD)
          continue;

        if (prm->timeseq==1)
        {
          int idx = frameExistsInFrameList(z);
          if (idx != -1)
          {
            float u = flowvector[index][z][0][y*width+x];
            float v = flowvector[index][z][1][y*width+x];

            if (0) // optical flow continuous
            {
              getIndicatorCostOpticalFlow(depth, z, di[x], fpvp->thetaO, thetaido, indcost, fpvp->opflowlocalframes.size(), idx);
              totalEmpirDistO[thetaido[0]] += u;
              totalEmpirDistO[thetaido[1]] += v;
            }
            else // binned optical flow features
            {
              getCostOpflowBinsuv(u, v, nOpC, di[x], fpvp->thetaO, thetaid, fpvp->opflowlocalframes.size(), idx);
              totalEmpirDistO[thetaid] ++;
            }
          }
        }
        // for us the costs are the optical flow values
        // since the theta value svanish after differentiation
      }
    }
  }
}



// this automatically handles interactive
void crf::updateQuantizedPairwiseEmpirDistribution(int index)
{

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  int depthrun = 0;

  /*
  if (fpvp->gradVZ == 1)
    depthrun = depth;
  else if (fpvp->gradVZ == 2 || fpvp->gradVZ == 3)
    depthrun = depth-1;
   */

  // changed so that x-y gradients on last slice can be used.
  depthrun = depth;

  uchar *di = 0;
  uchar *di2 = 0;
  uchar *di3 = 0;


  for (int z=0; z < depthrun; z++)
  {
    for (int y = 0; y < height-1; y++)
    {
      di = &gtdirImage[index][z].Pixel(0, y, 0);
      di2 = &gtdirImage[index][z].Pixel(0, y+1, 0);

      if ((fpvp->gradVZ == 2 || fpvp->gradVZ == 3)&&(z<depthrun-1))
        di3 = &gtdirImage[index][z+1].Pixel(0, y, 0);

      for (int x = 0; x < width-1; x++)
      {
        int gh, gv, gd;
        gh  = dirgrad[index][z].Pixel(x, y, 0);
        gv  = dirgrad[index][z].Pixel(x, y, 1);
        gd  = dirgrad[index][z].Pixel(x, y, 2);


        //CHECK need to make sure it has the gradContext option too
        // i am not fixing this gradOnly option
        if (fpvp->gradContext==0)
        {

          std::cout << " Not supposed to be working really check this !!! \n";
          exit(1);

          // Horizontal Edge
          if (di[x] != di[x+1])   // this is the border of the parts
            totalEmpirDistC[gh]++;

          // Vertical Edge
          if (di[x] != di2[x]) // this is the border of the parts
            totalEmpirDistC[gv]++;

          if (fpvp->gradVZ == 2)
          {
            // depth edge
            if (di[x] != di3[x]) // this is the border of the parts
              totalEmpirDistC[gd]++;
          }

        } else { //gradContext option

          if (fpvp->nV > 1)  // this makes sure it is not smooth constrast sensitive
          {

            // Horizontal Edge
            int d1 = (int) di[x];
            int d2 = (int) di[x+1];

            if (prm->unknown == 1 && (d1==prm->nD || d2==prm->nD)) {}
            else
            {
              if (d1>d2)
              {
                int swap = d2;
                d2 = d1;
                d1 = swap;
              }

              totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gh] ++;

            }


            // Vertical Edge
            d1 = (int) di[x];
            d2 = (int) di2[x];

            if (prm->unknown == 1 && (d1==prm->nD || d2==prm->nD)) {}
            else
            {
              if (d1>d2)
              {
                int swap = d2;
                d2 = d1;
                d1 = swap;
              }

              totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gv] ++;

            }

          }


          if ((fpvp->gradVZ == 2 && fpvp->nV > 1) && (z<depthrun-1))
          {

            // Depth Edge
            int d1 = (int) di[x];
            int d2 = (int) di3[x];

            if (prm->unknown == 1 && (d1==prm->nD || d2==prm->nD)) {}
            else
            {
              if (d1>d2)
              {
                int swap = d2;
                d2 = d1;
                d1 = swap;
              }

              totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gd] ++;
            }

          }


          if ((fpvp->gradVZ == 3 && fpvp->nZ > 1) &&(z<depthrun-1))
          {

            // Depth Edge
            int d1 = (int) di[x];
            int d2 = (int) di3[x];

            if (prm->unknown == 1 && (d1==prm->nD || d2==prm->nD)) {}
            else
            {
              if (d1>d2)
              {
                int swap = d2;
                d2 = d1;
                d1 = swap;
              }

              totalEmpirDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nGZ + gd] ++;

            }
          }
        }
      }
    }
  }
}


void crf::updateQuantizedPairwiseEmpirDistributionOpflow(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  assert(fpvp->gradVZ != 1);
  assert(globalP.getOpflowVolumeIndex() == index);

  if (fpvp->gradContext == 0)
  {
    std::cout << " Not implemented gradcontext = 0 \n";
  }

  uchar *di = 0;
  uchar *di2 = 0;


  std::map <std::pair <int, int>, std::pair <int, MRF::CostVal> >::iterator iter = globalP.edgeGlobalF[index].begin();

  for (; iter != globalP.edgeGlobalF[index].end(); iter++)
  {
    int node1 = iter->first.first;
    int node2 = iter->first.second;

    int optf = iter->second.first;
    int gradval = (int) iter->second.second;

    int z1 = node1 / (width*height);
    int z2 = node2 / (width*height);

    int rem1 = node1 - z1*(width*height);
    int rem2 = node2 - z2*(width*height);

    int y1 = rem1 / width;
    int y2 = rem2 / width;

    int x1 = rem1 - y1*width;
    int x2 = rem2 - y2*width;

    di = &gtdirImage[index][z1].Pixel(x1, y1, 0);
    di2 = &gtdirImage[index][z2].Pixel(x2, y2, 0);


    if (fpvp->nV > 1)  // this makes sure it is not smooth constrast sensitive
    {
      int d1 = (int) di[0];
      int d2 = (int) di2[0];

      if (prm->unknown == 1 && (d1==prm->nD || d2==prm->nD)) {}
      else
      {
        if (d1>d2)
        {
          int swap = d2;
          d2 = d1;
          d1 = swap;
        }

        if (optf == 0) // spatial ie x or y
        {
          totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] ++;
        } else { // optical flow connection
          if (fpvp->gradVZ == 2)
          {
            totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] ++;
          }
          if (fpvp->gradVZ == 3)
          {
            totalEmpirDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] ++;
          }
        }

      }
    }
  }
}



// this automatically handles interactive
void crf::updateSmoothPairwiseEmpirDistribution(int index)
{

  assert(fpvp->nV <= 1);
  assert(fpvp->nZ <= 1);


  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  int depthrun = 0;

  /*
  if (fpvp->gradVZ == 1)
    depthrun = depth;
  else if (fpvp->gradVZ == 2 || fpvp->gradVZ == 3)
    depthrun = depth-1;
  */

  // changed so that x-y gradients on last slice can be used.
  depthrun = depth;

  uchar *di = 0;
  uchar *di2 = 0;
  uchar *di3 = 0;


  for (int z=0; z < depthrun; z++)
  {
    for (int y = 0; y < height-1; y++)
    {
      di = &gtdirImage[index][z].Pixel(0, y, 0);
      di2 = &gtdirImage[index][z].Pixel(0, y+1, 0);

      if ((fpvp->gradVZ == 2 || fpvp->gradVZ == 3) && (z<depth-1))
        di3 = &gtdirImage[index][z+1].Pixel(0, y, 0);

      for (int x = 0; x < width-1; x++)
      {
        double valh, valv, vald;
        CvScalar s;
        s=cvGet2D(dircsgrad[index][z],y,x); // get the (i,j) pixel value
        valh = exp(s.val[0]);
        valv = exp(s.val[1]);
        vald = exp(s.val[2]);

        //CHECK need to make sure it has the gradContext option too
        // i am not fixing this gradOnly option
        if (fpvp->gradContext==0)
        {
          std::cout << " Not supposed to be working really check this !!! \n";
          exit(1);

        } else { //gradContext option

          if (fpvp->nV == 1)  // this makes sure it is smooth contrast sensitive
          {

            // Horizontal Edge
            int d1 = (int) di[x];
            int d2 = (int) di[x+1];

            if (prm->unknown == 1 && (d1==prm->nD || d2==prm->nD)) {}
            else if (d1 == d2)
            {} // since [yp != yq] in Vpq
            else
            {
              totalEmpirDistV[0] += valh;
            }


            // Vertical Edge
            d1 = (int) di[x];
            d2 = (int) di2[x];

            if (prm->unknown == 1 && (d1==prm->nD || d2==prm->nD)) {}
            else if (d1 == d2)
            {} // since [yp != yq] in Vpq
            else
            {
              totalEmpirDistV[0] += valv;
            }

          }


          if ((fpvp->gradVZ == 2 && fpvp->nV == 1)  && (z<depth-1))
          {

            // Depth Edge
            int d1 = (int) di[x];
            int d2 = (int) di3[x];

            if (prm->unknown == 1 && (d1==prm->nD || d2==prm->nD)) {}
            else if (d1 == d2)
            {} // since [yp != yq] in Vpq
            else
            {
              totalEmpirDistV[0] += vald;
            }

          }


          if ((fpvp->gradVZ == 3 && fpvp->nZ == 1)  && (z<depth-1))
          {

            // Depth Edge
            int d1 = (int) di[x];
            int d2 = (int) di3[x];

            if (prm->unknown == 1 && (d1==prm->nD || d2==prm->nD)) {}
            else if (d1 == d2)
            {} // since [yp != yq] in Vpq
            else
            {
              totalEmpirDistZ[0] += vald;
            }
          }
        }
      }
    }
  }
}


void crf::updateSmoothPairwiseEmpirDistributionOpflow(int index)
{
  // we currently don't update V and Z if they are single
  if ((fpvp->nV == 1 && fpvp->nZ == 1) || (fpvp->nV + fpvp->nZ == 1))
    return;


  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  assert(fpvp->gradVZ != 1);
  assert(globalP.getOpflowVolumeIndex() == index);

  if (fpvp->gradContext == 0)
  {
    std::cout << " Not implemented gradcontext = 0 \n";
  }

  uchar *di = 0;
  uchar *di2 = 0;


  std::map <std::pair <int, int>, std::pair <int, MRF::CostVal> >::iterator iter = globalP.edgeGlobalF[index].begin();

  for (; iter != globalP.edgeGlobalF[index].end(); iter++)
  {
    int node1 = iter->first.first;
    int node2 = iter->first.second;

    int optf = iter->second.first;
    int gradval = (int) iter->second.second;

    int z1 = node1 / (width*height);
    int z2 = node2 / (width*height);

    int rem1 = node1 - z1*(width*height);
    int rem2 = node2 - z2*(width*height);

    int y1 = rem1 / width;
    int y2 = rem2 / width;

    int x1 = rem1 - y1*width;
    int x2 = rem2 - y2*width;

    di = &gtdirImage[index][z1].Pixel(x1, y1, 0);
    di2 = &gtdirImage[index][z2].Pixel(x2, y2, 0);


    if (fpvp->nV > 1 || fpvp->nZ > 1)  // this makes sure it is not smooth constrast sensitive
    {
      int d1 = (int) di[0];
      int d2 = (int) di2[0];

      if (prm->unknown == 1 && (d1==prm->nD || d2==prm->nD)) {}
      else
      {
        if (d1>d2)
        {
          int swap = d2;
          d2 = d1;
          d1 = swap;
        }

        if (optf == 0) // spatial ie x or y
        {
          totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] ++;
        } else { // optical flow connection
          if (fpvp->gradVZ == 2)
          {
            totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] ++;
          }
          if (fpvp->gradVZ == 3)
          {
            totalEmpirDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] ++;
          }
        }

      }
    }
  }
}










// Note unknown region can currently only be used with fixed pairwise parameter
// ie when there is only one -v value.
// So currently there are no changes to computeEmpirDistV

// this will not learn parameters related to smoothed contrast sensitive gradients
void crf::computeEmpirDistV(int index)
{

  if (fpvp->nV == 0 && fpvp->nC == 0 && fpvp->nZ == 0)
    return; // nothing to do

  if (fpvp->csgrad == 0)
  {
    if (globalP.getOpflowConnection() == 1)
    {
      updateQuantizedPairwiseEmpirDistributionOpflow(index);
    } else {
      updateQuantizedPairwiseEmpirDistribution(index);
    }
  }
  else if (fpvp->csgrad == 1)
  {
    if (globalP.getOpflowConnection() == 1)
    {
      updateSmoothPairwiseEmpirDistributionOpflow(index);
    } else {
      updateSmoothPairwiseEmpirDistribution(index);
    }
  }

}



// mpe functions

// grayscale and // RGB or Yuv or Lab
// does not current have interactive mode
void crf::computeModelDistUmpe(int index)
{
	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nU == 0)
		return; // nothing to do

	int nbins = double(fpvp->nU)/prm->nD;
	double binner = fpvp->rangeI/nbins;

	int thetaid;
  int thetaidrgb[3]={0,0,0};
  int thetaidYuv[2]={0,0};


	uchar *pix1t = 0;
	unsigned short pix1 = 0;

	for (int z = 0; z < depth; z++)
	{
		for (int y = 0; y < height; y++)
		{
			uchar *di = &dirDisp[index][z].Pixel(0, y, 0);

      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

			for (int x = 0; x < width; x++)
			{
        if (prm->crfpv.hidden == 1) {} // use all pixels
        else if (prm->unknown == 1 && dig[x]==prm->nD) // unknown pixels do not count during learning parameters
          continue;

			  if (prm->timeseq==1)
				{
					pix1t = &indirImaget[index][z].Pixel(x, y, 0);
			        // getCostIntensityDiscBins(pix1t, nbins, di[x], thetaU, thetaidrgb);
					if (prm->featurepv.onlyUV == 1)
			      getCostIntensityDiscBinsuv(pix1t, nbins, di[x], fpvp->thetaU, thetaidYuv);
					else
					  getCostIntensityDiscBinsYuv(pix1t, nbins, di[x], fpvp->thetaU, thetaidYuv);

					if (prm->featurepv.onlyUV == 0)
			      totalModelDistU[thetaidYuv[0]]++;
			    totalModelDistU[thetaidYuv[1]]++;

				}
				else if (prm->timeseq==0)
				{
					pix1 = indirImage[index][z].Pixel(x, y, 0);
					getCostIntensityDiscBins(pix1, nbins, di[x], fpvp->thetaU, thetaid);
					totalModelDistU[thetaid] ++;
				}

				// for us the costs are the theta param, but in differentiation, they vanish,
				// so its just whether the feature is enabled or not.
				// move costs pointer to next pixel
			}
		}
	}
}


void crf::computeModelDistAmpe(int index)
{

	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nA == 0)
		return; // nothing to do

	int nT = fpvp->nA/prm->nD; //this gives centers per class
	// this will work for app==2 and app==1 because in app==2, we have nD*the number of total clusters as the number of parameters
	// while in app==1 it is nD*number of clusters of that class (only that clusters are also = number of params)

    int thetaid;


	for (int z=0; z < depth; z++)
	{
		for (int y = 0; y < height; y++)
		{
			uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
			uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

			for (int x = 0; x < width; x++)
			{
        if (prm->crfpv.hidden == 1) {} // use all pixels
        else if (prm->unknown == 1 && dig[x]==prm->nD)         // unknown pixels do not count during learning parameters
          continue;

				int d = di[x];
				int dm = d;
                if (fpvp->app==2)
                  dm=0;

            	if (((int)appdirImage[index][z][dm].Pixel(x, y, 0) >= nT && fpvp->app==1) || ((int)appdirImage[index][z][0].Pixel(x, y, 0) >= fpvp->nA && fpvp->app==2))
            	{
            		continue;
            	}
                getCostAppDiscPatch((int)appdirImage[index][z][dm].Pixel(x, y, 0), appclass[dm]->getPatchSize(), x, y, width, height, fpvp->thetaA, d, nT, thetaid);

                // assume the patch size is the same for all classes
                int patchSize = appclass[dm]->getPatchSize();

                if (x>(int)patchSize/2 && x<width-patchSize/2 && y>(int)patchSize/2 && y<height-patchSize/2)
                	totalModelDistA[thetaid] ++;

			}
		}
	}
}




void crf::computeModelDistHmpe(int index)
{

	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nH == 0)
		return; // nothing to do

	int nHpC = fpvp->nH/prm->nD; // HoG vocabulary size

    int thetaid;

    for (int z = 0; z < depth; z++)
    {
    	for (int y = 0; y < height; y++)
    	{
    		uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
    		uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

    		for (int x = 0; x < width; x++)
    		{
          if (prm->crfpv.hidden == 1) {} // use all pixels
          else if (prm->unknown == 1 && dig[x]==prm->nD)         // unknown pixels do not count during learning parameters
            continue;

    			unsigned short hogval = hogdirImage[index][z].Pixel(x, y, 0);
    			getCostHoGDiscBins((int) hogval, di[x], nHpC, fpvp->thetaH, thetaid);
    			//getCostHoGDiscBinsTemp((int) hogval, di[x], fpvp->thetaH, thetaid);

    			// for us the costs are the theta param, but in differentiation, they vanish,
    			// so its just whether the feature is enabled or not.
    			// move costs pointer to next pixel
    			if (thetaid >= 0)
    				totalModelDistH[thetaid] ++;

    		}
    	}
    }
}

void crf::computeModelDistBmpe(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nB == 0)
    return; // nothing to do

  int thetaid;

  for (int z = 0; z < depth; z++)
  {
    for (int y = 0; y < height; y++)
    {
      uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

      for (int x = 0; x < width; x++)
      {
        if (prm->crfpv.hidden == 1) {} // use all pixels
        else if (prm->unknown == 1 && dig[x]==prm->nD)         // unknown pixels do not count during learning parameters
          continue;

        getCostBias(di[x], fpvp->thetaB, thetaid);

        // for us the costs are the theta param, but in differentiation, they vanish,
        // so its just whether the feature is enabled or not.
        // move costs pointer to next pixel
        if (thetaid >= 0)
          totalModelDistB[thetaid] ++;

      }
    }
  }
}

void crf::computeModelDistEmpe(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nE == 0)
    return; // nothing to do

  int thetaid;

  for (int z = 0; z < depth; z++)
  {
    for (int y = 0; y < height; y++)
    {
      uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

      for (int x = 0; x < width; x++)
      {
        if (prm->crfpv.hidden == 1) {} // use all pixels
        else if (prm->unknown == 1 && dig[x]==prm->nD)         // unknown pixels do not count during learning parameters
          continue;

        getCostInverseClassSize(di[x], inverseClassSize, fpvp->thetaE, thetaid);

        // for us the costs are the theta param, but in differentiation, they vanish,
        // so its just whether the feature is enabled or not.
        // move costs pointer to next pixel
        if (thetaid >= 0)
          totalModelDistE[thetaid] += inverseClassSize[di[x]];

      }
    }
  }
}

void crf::computeModelDistMmpe(int index)
{
	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nM == 0)
		return; // nothing to do

	int nMpC = fpvp->nM/prm->nD; // Mot vocabulary size

    int thetaid;

    for (int z = 0; z < depth; z++)
    {
    	for (int y = 0; y < height; y++)
    	{
    		uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
    		uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

    		for (int x = 0; x < width; x++)
    		{
          if (prm->crfpv.hidden == 1) {} // use all pixels
          else if (prm->unknown == 1 && dig[x]==prm->nD)         // unknown pixels do not count during learning parameters
            continue;

    			unsigned short motval = motdirImage[index][z].Pixel(x, y, 0);
    			getCostmotDiscBins((int) motval, di[x], nMpC, fpvp->thetaM, thetaid);

    			totalModelDistM[thetaid] ++;

    		}
    	}
    }
}


void crf::computeModelDistLmpe(int index)
{

	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nL == 0)
		return; // nothing to do

	int nLpC = fpvp->nL/prm->nD;  // this variable will be useful only for cube setting (loc==1) or loc==9

  int thetaid;

  int zslice = 0;
	if (timeseq == 0 && (fpvp->loc != 0 || fpvp->klrpv.locKlr ==2))
		zslice = startSliceNo[index];

    for (int z = 0; z < depth; z++)
    {

    	for (int y = 0; y < height; y++)
    	{
    		uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
    		uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

    		for (int x = 0; x < width; x++)
    		{
          if (prm->crfpv.hidden == 1) {} // use all pixels
          else if (prm->unknown == 1 && dig[x]==prm->nD)         // unknown pixels do not count during learning parameters
            continue;

    			if (fpvp->loc==1)
    			{
    				/*
    				unsigned short locpix = locdirImage[index][di[x]][z].Pixel(x,y,0);
    				getCostLocationDiscCube(locpix, x, y, z, fpvp->loccubex, fpvp->loccubey, fpvp->loccubez, numTrainingPats, nLpC, fpvp->thetaL, width, height, depth, di[x], thetaid);
    				double locpixd = (double)locpix/ (fpvp->loccubex * fpvp->loccubey * fpvp->loccubez * numTrainingPats);

    				totalEmpirDistL[thetaid] += (1-locpixd);
    				*/

    				std::cout << " Need to fix this - location " << std::endl;
    				exit(1);


    			} else if (fpvp->loc==6 || fpvp->loc==8) { // CRF hard assignment and soft

    				getCostLocationDiscGaussian(LocationMVclass, z, zslice, x, y, fpvp->thetaL, di[x], fpvp->loc, fpvp->generative, thetaid);
    				totalModelDistL[thetaid]++;

    			} else if (fpvp->loc==9 ||  (fpvp->loc==7 && timeseq == 1)) {  //CRF soft or hard assignment redefined

    				getCostLocationDiscGaussian2(LocationMVclass, z, zslice, x, y, fpvp->thetaL, di[x], fpvp->loc, fpvp->generative, prm->nD, thetaid);
    				totalModelDistL[thetaid]++;

    			} else if  (fpvp->loc==7 && timeseq == 0) {  //CRF soft or hard assignment redefined

    				getCostLocationDiscGaussian3(index, z, x, y, di[x], fpvp->generative, prm->nD, thetaid);
    				totalModelDistL[thetaid]++;
    			}



    		}
    	}
    }
}


void crf::computeModelDistWmpe(int index)
{
	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->klrpv.nW == 0)
		return; // nothing to do

	int dimW = wparam.size2();
	int Ntr = wparam.size1();


	int zslice = 0;

	if (timeseq == 0 && (fpvp->loc != 0 || fpvp->klrpv.locKlr ==2))
		zslice = startSliceNo[index];


	uchar *pix1t = 0;
	unsigned short pix1 = 0;

	ublas::vector<DBL_TYPE> Xtest(Xtrain.size2());
	ublas::vector<DBL_TYPE> f;


	for (int z = 0; z < depth; z++)
	{

		matrix<DBL_TYPE> hogdescs, appdescs;

		if (fpvp->klrpv.klr==1)
			readklrdescfiles(hogdescs, appdescs, index, z+zslice);

		for (int y = 0; y < height; y++)
		{
			uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
			uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

			for (int x = 0; x < width; x++)
			{
        if (prm->crfpv.hidden == 1) {} // use all pixels
        else if (prm->unknown == 1 && dig[x]==prm->nD)         // unknown pixels do not count during learning parameters
          continue;

				getDataCostKlr(f, Xtest, index, z, x, y, zslice, hogdescs, appdescs, width, height);

				for (int j=0; j<Ntr; j++)
				{
					double Kval = ivm->classify_klrexp(&Xtest, j);

					int thetaid = j*dimW + di[x];
					totalModelDistW[thetaid] += Kval;
				}

			}
		}
	}
}



void crf::computeModelDistOmpe(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (prm->timeseq!=1)
    return;

  if (fpvp->nO == 0)
    return; // nothing to do

  int nOpC = 0;
  if (fpvp->opflowlocal == 1)
  {
    nOpC = fpvp->nO/(prm->nD * fpvp->opflowlocalframes.size()); // there are (u,v) variables per frame except last frame for continuous features
  // but is more complicated for bin features (modeled similar to color uv.
  }

  int thetaid;
  int thetaido[2]={0,0};
  float indcost[2] = {0.0, 0.0};

  // the last frame should not be counted towards the cost
  for (int z = 0; z < depth-1; z++)
  {
    for (int y = 0; y < height; y++)
    {
      uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);
      for (int x = 0; x < width; x++)
      {
        if (prm->crfpv.hidden == 1) {} // use all pixels
        else if (prm->unknown == 1 && dig[x]==prm->nD)         // unknown pixels do not count during learning parameters
          continue;

        int idx = frameExistsInFrameList(z);
        if (idx != -1)
        {
          float u = flowvector[index][z][0][y*width+x];
          float v = flowvector[index][z][1][y*width+x];
          if (0) // optical flow continuous
          {
            getIndicatorCostOpticalFlow(depth, z, di[x], fpvp->thetaO, thetaido, indcost, fpvp->opflowlocalframes.size(), idx);
            totalModelDistO[thetaido[0]] +=  u;
            totalModelDistO[thetaido[1]] +=  v;
          }
          else // binned optical flow features
          {
            getCostOpflowBinsuv(u, v, nOpC, di[x], fpvp->thetaO, thetaid, fpvp->opflowlocalframes.size(), idx);
            totalModelDistO[thetaid] ++;
          }
        }
        // consider the costs for the theta param,
        // since they exist after differentiation
        // no use of the indcost since they are the theta values

      }
    }

  }
}





void crf::updateQuantizedPairwiseModelDistributionmpeOpflow(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  assert(fpvp->gradVZ != 1);
  assert(globalP.getOpflowVolumeIndex() == index);

  if (fpvp->gradContext == 0)
  {
    std::cout << " Not implemented gradcontext = 0 \n";
  }

  uchar *di = 0;
  uchar *di2 = 0;

  uchar *dig = 0;
  uchar *dig2 = 0;


  std::map <std::pair <int, int>, std::pair <int, MRF::CostVal> >::iterator iter = globalP.edgeGlobalF[index].begin();

  for (; iter != globalP.edgeGlobalF[index].end(); iter++)
  {
    int node1 = iter->first.first;
    int node2 = iter->first.second;

    int optf = iter->second.first;
    int gradval = (int) iter->second.second;

    int z1 = node1 / (width*height);
    int z2 = node2 / (width*height);

    int rem1 = node1 - z1*(width*height);
    int rem2 = node2 - z2*(width*height);

    int y1 = rem1 / width;
    int y2 = rem2 / width;

    int x1 = rem1 - y1*width;
    int x2 = rem2 - y2*width;

    di = &dirDisp[index][z1].Pixel(x1, y1, 0);
    di2 = &dirDisp[index][z2].Pixel(x2, y2, 0);

    dig = &gtdirImage[index][z1].Pixel(x1, y1, 0);
    dig2 = &gtdirImage[index][z2].Pixel(x2, y2, 0);

    if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
    else
    {

      if (fpvp->nV > 1)  // this makes sure it is not smooth constrast sensitive
      {
        int d1 = (int) di[0];
        int d2 = (int) di2[0];

        if (d1>d2)
        {
          int swap = d2;
          d2 = d1;
          d1 = swap;
        }

        if (optf == 0) // spatial ie x or y
        {
          totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] ++;
        } else { // optical flow connection
          if (fpvp->gradVZ == 2)
          {
            totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] ++;
          }
          if (fpvp->gradVZ == 3)
          {
            totalModelDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] ++;
          }
        }

      }
    }
  }
}


void crf::updateSmoothPairwiseModelDistributionmpeOpflow(int index)
{
  // we currently don't update V and Z if they are single
  if ((fpvp->nV == 1 && fpvp->nZ == 1) || (fpvp->nV + fpvp->nZ == 1))
    return;

  assert(0); // might need to double check validity of code in this function

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  assert(fpvp->gradVZ != 1);
  assert(globalP.getOpflowVolumeIndex() == index);

  if (fpvp->gradContext == 0)
  {
    std::cout << " Not implemented gradcontext = 0 \n";
  }

  uchar *di = 0;
  uchar *di2 = 0;

  uchar *dig = 0;
  uchar *dig2 = 0;


  std::map <std::pair <int, int>, std::pair <int, MRF::CostVal> >::iterator iter = globalP.edgeGlobalF[index].begin();

  for (; iter != globalP.edgeGlobalF[index].end(); iter++)
  {
    int node1 = iter->first.first;
    int node2 = iter->first.second;

    int optf = iter->second.first;
    int gradval = (int) iter->second.second;

    int z1 = node1 / (width*height);
    int z2 = node2 / (width*height);

    int rem1 = node1 - z1*(width*height);
    int rem2 = node2 - z2*(width*height);

    int y1 = rem1 / width;
    int y2 = rem2 / width;

    int x1 = rem1 - y1*width;
    int x2 = rem2 - y2*width;

    di = &dirDisp[index][z1].Pixel(x1, y1, 0);
    di2 = &dirDisp[index][z2].Pixel(x2, y2, 0);

    dig = &gtdirImage[index][z1].Pixel(x1, y1, 0);
    dig2 = &gtdirImage[index][z2].Pixel(x2, y2, 0);

    if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
    else
    {

      if (fpvp->nV > 1 || fpvp->nZ > 1)  // this makes sure it is not smooth constrast sensitive
      {
        int d1 = (int) di[0];
        int d2 = (int) di2[0];

        if (d1>d2)
        {
          int swap = d2;
          d2 = d1;
          d1 = swap;
        }

        if (optf == 0) // spatial ie x or y
        {
          totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] ++;
        } else { // optical flow connection
          if (fpvp->gradVZ == 2)
          {
            totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] ++;
          }
          if (fpvp->gradVZ == 3)
          {
            totalModelDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] ++;
          }
        }

      }
    }
  }
}






// Note unknown region can currently only be used with fixed pairwise parameter
// ie when there is only one -v value.
// So currently there are no changes to computeModelDistVmpe

// this will not learn parameters related to smoothed contrast sensitive gradients
void crf::computeModelDistVmpe(int index)
{

  if (fpvp->nV == 0 && fpvp->nC == 0 && fpvp->nZ == 0)
    return; // nothing to do

  if (fpvp->csgrad == 0)
  {
    if (globalP.getOpflowConnection() == 1)
    {
      updateQuantizedPairwiseModelDistributionmpeOpflow(index);
    } else {
      updateQuantizedPairwiseModelDistributionmpe(index);
    }
  }
  else if (fpvp->csgrad == 1)
  {
    if (globalP.getOpflowConnection() == 1)
    {
      updateSmoothPairwiseModelDistributionmpeOpflow(index);
    } else {
      updateSmoothPairwiseModelDistributionmpe(index);
    }
  }

}



void crf::updateQuantizedPairwiseModelDistributionmpe(int index)
{

	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	int depthrun = 0;

	/*
	if (fpvp->gradVZ == 1)
		depthrun = depth;
	else if (fpvp->gradVZ == 2 || fpvp->gradVZ == 3)
		depthrun = depth-1;
	 */
	depthrun = depth;

	uchar *di = 0;
	uchar *di2 = 0;
	uchar *di3 = 0;

  uchar *dig = 0;
  uchar *dig2 = 0;
  uchar *dig3 = 0;


	for (int z=0; z < depthrun; z++)
	{
		for (int y = 0; y < height-1; y++)
		{
			di = &dirDisp[index][z].Pixel(0, y, 0);
			di2 = &dirDisp[index][z].Pixel(0, y+1, 0);

			if ((fpvp->gradVZ == 2 || fpvp->gradVZ == 3)&&(z<depthrun-1))
				di3 = &dirDisp[index][z+1].Pixel(0, y, 0);


			for (int x = 0; x < width-1; x++)
			{
				int gh, gv, gd;
				gh	= dirgrad[index][z].Pixel(x, y, 0);
				gv	= dirgrad[index][z].Pixel(x, y, 1);
				gd	= dirgrad[index][z].Pixel(x, y, 2);



				//CHECK need to make sure it has the gradContext option too
				// i am not fixing this gradOnly option
				if (fpvp->gradContext==0)
				{

					std::cout << " Not supposed to be working really check this !!! \n";
					exit(1);

					// Horizontal Edge
					if (di[x] != di[x+1])   // this is the border of the parts
						totalModelDistC[gh]++;

					// Vertical Edge
					if (di[x] != di2[x]) // this is the border of the parts
						totalModelDistC[gv]++;

					if (fpvp->gradVZ == 2)
					{
						// depth edge
						if (di[x] != di3[x]) // this is the border of the parts
							totalModelDistC[gd]++;
					}

				} else { //gradContext option

					if (fpvp->nV > 1)  // this makes sure it is not smooth constrast sensitive
					{

						// Horizontal Edge
						int d1 = (int) di[x];
						int d2 = (int) di[x+1];

				    dig = &gtdirImage[index][z].Pixel(x, y, 0);
				    dig2 = &gtdirImage[index][z].Pixel(x+1, y, 0);

				    if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
				    else
				    {
				      if (d1>d2)
              {
                int swap = d2;
                d2 = d1;
                d1 = swap;
              }

				      totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gh] ++;
				    }



						// Vertical Edge
						d1 = (int) di[x];
						d2 = (int) di2[x];

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z].Pixel(x, y+1, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
            else
            {
              if (d1>d2)
              {
                int swap = d2;
                d2 = d1;
                d1 = swap;
              }

              totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gv] ++;

            }

					}


					if (fpvp->gradVZ == 2 && fpvp->nV > 1 && z < depthrun-1)
					{

						// Depth Edge
						int d1 = (int) di[x];
						int d2 = (int) di3[x];

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z+1].Pixel(x, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
            else
            {
              if (d1>d2)
              {
                int swap = d2;
                d2 = d1;
                d1 = swap;
              }

              totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gd] ++;
            }

					}


					if (fpvp->gradVZ == 3 && fpvp->nZ > 1 && z < depthrun-1)
					{

						// Depth Edge
						int d1 = (int) di[x];
						int d2 = (int) di3[x];

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z+1].Pixel(x, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
            else
            {
              if (d1>d2)
              {
                int swap = d2;
                d2 = d1;
                d1 = swap;
              }
              totalModelDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nGZ + gd] ++;

            }
					}
				}
			}
		}
	}
}



void crf::updateSmoothPairwiseModelDistributionmpe(int index)
{

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  int depthrun = 0;

  /*
  if (fpvp->gradVZ == 1)
    depthrun = depth;
  else if (fpvp->gradVZ == 2 || fpvp->gradVZ == 3)
    depthrun = depth-1;
   */

  depthrun = depth;

  uchar *di = 0;
  uchar *di2 = 0;
  uchar *di3 = 0;

  uchar *dig = 0;
  uchar *dig2 = 0;
  uchar *dig3 = 0;


  for (int z=0; z < depthrun; z++)
  {
    for (int y = 0; y < height-1; y++)
    {
      di = &dirDisp[index][z].Pixel(0, y, 0);
      di2 = &dirDisp[index][z].Pixel(0, y+1, 0);

      if ((fpvp->gradVZ == 2 || fpvp->gradVZ == 3)  && (z<depth-1))
        di3 = &dirDisp[index][z+1].Pixel(0, y, 0);


      for (int x = 0; x < width-1; x++)
      {
        double valh, valv, vald;
        CvScalar s;
        s=cvGet2D(dircsgrad[index][z],y,x); // get the (i,j) pixel value
        valh = exp(s.val[0]);
        valv = exp(s.val[1]);
        vald = exp(s.val[2]);

        //CHECK need to make sure it has the gradContext option too
        // i am not fixing this gradOnly option
        if (fpvp->gradContext==0)
        {
          std::cout << " Not supposed to be working really check this !!! \n";
          exit(1);

        } else { //gradContext option

          if (fpvp->nV == 1)  // this makes sure it is smooth constrast sensitive
          {

            // Horizontal Edge
            int d1 = (int) di[x];
            int d2 = (int) di[x+1];

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z].Pixel(x+1, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
            else if (d1 == d2)
            {} // since [yp != yq] in Vpq
            else
            {
              totalModelDistV[0] += valh;
            }

            // Vertical Edge
            d1 = (int) di[x];
            d2 = (int) di2[x];

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z].Pixel(x, y+1, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
            else if (d1 == d2)
            {} // since [yp != yq] in Vpq
            else
            {
              totalModelDistV[0] += valv;
            }

          }


          if (fpvp->gradVZ == 2 && fpvp->nV == 1 && z < depthrun-1)
          {

            // Depth Edge
            int d1 = (int) di[x];
            int d2 = (int) di3[x];

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z+1].Pixel(x, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
            else if (d1 == d2)
            {} // since [yp != yq] in Vpq
            else
            {
              totalModelDistV[0] += vald;
            }
          }


          if (fpvp->gradVZ == 3 && fpvp->nZ == 1 && z < depthrun-1)
          {

            // Depth Edge
            int d1 = (int) di[x];
            int d2 = (int) di3[x];

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z+1].Pixel(x, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
            else if (d1 == d2)
            {} // since [yp != yq] in Vpq
            else
            {
              totalModelDistZ[0] += vald;
            }
          }
        }
      }
    }
  }
}




// grayscale and // RGB or Yuv or Lab
// does not current have interactive mode
void crf::computeModelDistU(int index)
{
	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nU == 0)
		return; // nothing to do

	int nbins = double(fpvp->nU)/prm->nD;
	double binner = fpvp->rangeI/nbins;

	int thetaid;
    int thetaidrgb[3]={0,0,0};
    int thetaidYuv[2]={0,0};


	uchar *pix1t = 0;
	unsigned short pix1 = 0;

	MRFEnergy* mrfe = (MRFEnergy*)mrf;

  float *belief = new float[prm->nD];
  int pixel=0;

	for (int z = 0; z < depth; z++)
	{
		for (int y = 0; y < height; y++)
		{
			uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
			uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

			for (int x = 0; x < width; x++)
			{
        if (prm->crfpv.hidden == 1) {} // use all pixels
        else if (prm->unknown == 1 && dig[x]==prm->nD)         // unknown pixels do not count during learning parameters
        {
          pixel++;
          continue;
        }

				mrfe->getBelief(pixel,belief);

				for (int d=0; d<prm->nD; d++)
				{

					if (prm->timeseq==1)
					{
						pix1t = &indirImaget[index][z].Pixel(x, y, 0);
						// getCostIntensityDiscBins(pix1t, nbins, d, thetaU, thetaidrgb);
						if (prm->featurepv.onlyUV == 1)
						  getCostIntensityDiscBinsuv(pix1t, nbins, d, fpvp->thetaU, thetaidYuv);
						else
						  getCostIntensityDiscBinsYuv(pix1t, nbins, d, fpvp->thetaU, thetaidYuv);

						if (prm->featurepv.onlyUV == 0)
						  totalModelDistU[thetaidYuv[0]] += exp(belief[d]);
						totalModelDistU[thetaidYuv[1]] += exp(belief[d]);

					}
					else if (prm->timeseq==0)
					{
						pix1 = indirImage[index][z].Pixel(x, y, 0);
						getCostIntensityDiscBins(pix1, nbins, d, fpvp->thetaU, thetaid);
						totalModelDistU[thetaid] += exp(belief[d]);
					}

				}

				// for us the costs are the theta param, but in differentiation, they vanish,
				// so its just whether the feature is enabled or not.
				// move costs pointer to next pixel

				pixel++;

			}
		}

	}

	delete [] belief;

}


void crf::computeModelDistA(int index)
{

	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nA == 0)
		return; // nothing to do

	int nT = fpvp->nA/prm->nD; //this gives centers per class
	// this will work for app==2 and app==1 because in app==2, we have nD*the number of total clusters as the number of parameters
	// while in app==1 it is nD*number of clusters of that class (only that clusters are also = number of params)

    int thetaid;

    float *belief = new float[prm->nD];
    int pixel=0;

    MRFEnergy* mrfe = (MRFEnergy*)mrf;

	for (int z=0; z < depth; z++)
	{
		for (int y = 0; y < height; y++)
		{
			uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
			uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

			for (int x = 0; x < width; x++)
			{
        if (prm->crfpv.hidden == 1) {} // use all pixels
        else if (prm->unknown == 1 && dig[x]==prm->nD)         // unknown pixels do not count during learning parameters
        {
          pixel++;
          continue;
        }

				mrfe->getBelief(pixel,belief);

				for (int d=0; d<prm->nD; d++)
				{
					int dm = d;
					if (fpvp->app==2)
					  dm=0;

					if (((int)appdirImage[index][z][dm].Pixel(x, y, 0) >= nT && fpvp->app==1) || ((int)appdirImage[index][z][0].Pixel(x, y, 0) >= fpvp->nA && fpvp->app==2))
					{
						continue;
					}
					getCostAppDiscPatch((int)appdirImage[index][z][dm].Pixel(x, y, 0), appclass[dm]->getPatchSize(), x, y, width, height, fpvp->thetaA, d, nT, thetaid);

					// assume the patch size is the same for all classes
					int patchSize = appclass[dm]->getPatchSize();

					if (x>(int)patchSize/2 && x<width-patchSize/2 && y>(int)patchSize/2 && y<height-patchSize/2)
						totalModelDistA[thetaid] ++;
				}

        		pixel++;

			}
		}
	}

	delete [] belief;

}




void crf::computeModelDistH(int index)
{

	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nH == 0)
		return; // nothing to do

	int nHpC = fpvp->nH/prm->nD; // HoG vocabulary size

    int thetaid;

    MRFEnergy* mrfe = (MRFEnergy*)mrf;

    float *belief = new float[prm->nD];
    int pixel=0;

    for (int z = 0; z < depth; z++)
    {
    	for (int y = 0; y < height; y++)
    	{
    		uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
    		uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

    		for (int x = 0; x < width; x++)
    		{
          if (prm->crfpv.hidden == 1) {} // use all pixels
          else if (prm->unknown == 1 && dig[x]==prm->nD)         // unknown pixels do not count during learning parameters
          {
            pixel++;
            continue;
          }

				mrfe->getBelief(pixel,belief);

				for (int d=0; d<prm->nD; d++)
				{
					unsigned short hogval = hogdirImage[index][z].Pixel(x, y, 0);
					getCostHoGDiscBins((int) hogval, d, nHpC, fpvp->thetaH, thetaid);
					//getCostHoGDiscBinsTemp((int) hogval, d, fpvp->thetaH, thetaid);

					// for us the costs are the theta param, but in differentiation, they vanish,
					// so its just whether the feature is enabled or not.
					// move costs pointer to next pixel
					if (thetaid >= 0)
						totalModelDistH[thetaid] += exp(belief[d]);

				}

    			pixel++;

    		}
    	}
    }

    delete [] belief;

}


void crf::computeModelDistB(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nB == 0)
    return; // nothing to do

  int thetaid;

  MRFEnergy* mrfe = (MRFEnergy*)mrf;

  float *belief = new float[prm->nD];
  int pixel=0;

  for (int z = 0; z < depth; z++)
  {
    for (int y = 0; y < height; y++)
    {
      uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

      for (int x = 0; x < width; x++)
      {
        if (prm->crfpv.hidden == 1) {} // use all pixels
        else if (prm->unknown == 1 && dig[x]==prm->nD)         // unknown pixels do not count during learning parameters
        {
          pixel++;
          continue;
        }

        mrfe->getBelief(pixel,belief);

        for (int d=0; d<prm->nD; d++)
        {
          getCostBias(d, fpvp->thetaB, thetaid);

          // for us the costs are the theta param, but in differentiation, they vanish,
          // so its just whether the feature is enabled or not.
          // move costs pointer to next pixel
          if (thetaid >= 0)
            totalModelDistB[thetaid] += exp(belief[d]);

        }

        pixel++;

      }
    }
  }

  delete [] belief;
}

void crf::computeModelDistE(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nE == 0)
    return; // nothing to do

  int thetaid;

  MRFEnergy* mrfe = (MRFEnergy*)mrf;

  float *belief = new float[prm->nD];
  int pixel=0;

  for (int z = 0; z < depth; z++)
  {
    for (int y = 0; y < height; y++)
    {
      uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

      for (int x = 0; x < width; x++)
      {
        if (prm->crfpv.hidden == 1) {} // use all pixels
        else if (prm->unknown == 1 && dig[x]==prm->nD)         // unknown pixels do not count during learning parameters
        {
          pixel++;
          continue;
        }

        mrfe->getBelief(pixel,belief);

        for (int d=0; d<prm->nD; d++)
        {
          getCostInverseClassSize(d, inverseClassSize, fpvp->thetaE, thetaid);

          // for us the costs are the theta param, but in differentiation, they vanish,
          // so its just whether the feature is enabled or not.
          // move costs pointer to next pixel
          // assert(0); // do i need to multiply by classinverseSize here ? above comment is false
          if (thetaid >= 0)
            totalModelDistE[thetaid] += exp(belief[d]) * inverseClassSize[d];

        }

        pixel++;

      }
    }
  }

  delete [] belief;
}


void crf::computeModelDistM(int index)
{
	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nM == 0)
		return; // nothing to do

	int nMpC = fpvp->nM/prm->nD; // Mot vocabulary size

    int thetaid;

    MRFEnergy* mrfe = (MRFEnergy*)mrf;

    float *belief = new float[prm->nD];
    int pixel=0;

    for (int z = 0; z < depth; z++)
    {
    	for (int y = 0; y < height; y++)
    	{
    		uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
    		uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

    		for (int x = 0; x < width; x++)
    		{
          if (prm->crfpv.hidden == 1) {} // use all pixels
          else if (prm->unknown == 1 && dig[x]==prm->nD)         // unknown pixels do not count during learning parameters
          {
            pixel++;
            continue;
          }

				mrfe->getBelief(pixel,belief);

				for (int d=0; d<prm->nD; d++)
				{
					unsigned short motval = motdirImage[index][z].Pixel(x, y, 0);
					getCostmotDiscBins((int) motval, d, nMpC, fpvp->thetaM, thetaid);

					totalModelDistM[thetaid] += exp(belief[d]);
				}

    			pixel++;

    		}
    	}
    }

    delete [] belief;

}


void crf::computeModelDistL(int index)
{

	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nL == 0)
		return; // nothing to do

	int nLpC = fpvp->nL/prm->nD;  // this variable will be useful only for cube setting (loc==1) or loc==9

    int thetaid;

    MRFEnergy* mrfe = (MRFEnergy*)mrf;

    int zslice = 0;
	if (timeseq == 0 && (fpvp->loc != 0 || fpvp->klrpv.locKlr ==2))
		zslice = startSliceNo[index];

    float *belief = new float[prm->nD];
    int pixel=0;

    for (int z = 0; z < depth; z++)
    {

    	for (int y = 0; y < height; y++)
    	{
    		uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
    		uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

    		for (int x = 0; x < width; x++)
    		{
          if (prm->crfpv.hidden == 1) {} // use all pixels
          else if (prm->unknown == 1 && dig[x]==prm->nD)         // unknown pixels do not count during learning parameters
          {
            pixel++;
            continue;
          }

				mrfe->getBelief(pixel,belief);

				for (int d=0; d<prm->nD; d++)
				{
					if (fpvp->loc==1)
					{
						/*
						unsigned short locpix = locdirImage[index][di[x]][z].Pixel(x,y,0);
						getCostLocationDiscCube(locpix, x, y, z, fpvp->loccubex, fpvp->loccubey, fpvp->loccubez, numTrainingPats, nLpC, fpvp->thetaL, width, height, depth, di[x], thetaid);
						double locpixd = (double)locpix/ (fpvp->loccubex * fpvp->loccubey * fpvp->loccubez * numTrainingPats);

						totalModelDistL[thetaid] += += exp(belief[d])*(1-locpixd);
						*/

						std::cout << " Need to fix this - location " << std::endl;
						exit(1);


					} else if (fpvp->loc==6 || fpvp->loc==8) { // CRF hard assignment and soft

						getCostLocationDiscGaussian(LocationMVclass, z, zslice, x, y, fpvp->thetaL, d, fpvp->loc, fpvp->generative, thetaid);
						totalModelDistL[thetaid] += exp(belief[d]);

					} else if (fpvp->loc==9 ||  (fpvp->loc==7 && timeseq == 1)) {  //CRF soft or hard assignment redefined

						getCostLocationDiscGaussian2(LocationMVclass, z, zslice, x, y, fpvp->thetaL, d, fpvp->loc, fpvp->generative, prm->nD, thetaid);
						totalModelDistL[thetaid] += exp(belief[d]);

					} else if  (fpvp->loc==7 && timeseq == 0) {  //CRF hard assignment redefined

						getCostLocationDiscGaussian3(index, z, x, y, d, fpvp->generative, prm->nD, thetaid);
						totalModelDistL[thetaid] += exp(belief[d]);
					}
				}

    			pixel++;

    		}
    	}
    }

    delete [] belief;

}


void crf::computeModelDistW(int index)
{
	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->klrpv.nW == 0)
		return; // nothing to do

	int dimW = wparam.size2();
	int Ntr = wparam.size1();

	MRFEnergy* mrfe = (MRFEnergy*)mrf;

	int zslice = 0;

	if (timeseq == 0 && (fpvp->loc != 0 || fpvp->klrpv.locKlr ==2))
		zslice = startSliceNo[index];


	uchar *pix1t = 0;
	unsigned short pix1 = 0;

	ublas::vector<DBL_TYPE> Xtest(Xtrain.size2());
	ublas::vector<DBL_TYPE> f;


    float *belief = new float[prm->nD];
    int pixel=0;

	for (int z = 0; z < depth; z++)
	{

		matrix<DBL_TYPE> hogdescs, appdescs;

		if (fpvp->klrpv.klr==1)
			readklrdescfiles(hogdescs, appdescs, index, z+zslice);

		for (int y = 0; y < height; y++)
		{
			uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
			uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

			for (int x = 0; x < width; x++)
			{
        if (prm->crfpv.hidden == 1) {} // use all pixels
        else if (prm->unknown == 1 && dig[x]==prm->nD)         // unknown pixels do not count during learning parameters
        {
          pixel++;
          continue;
        }

				mrfe->getBelief(pixel,belief);


				getDataCostKlr(f, Xtest, index, z, x, y, zslice, hogdescs, appdescs, width, height);

				for (int j=0; j<Ntr; j++)
				{
					double Kval = ivm->classify_klrexp(&Xtest, j);

					for (int d=0; d<prm->nD; d++)
					{
						int thetaid = j*dimW + d;
						totalModelDistW[thetaid] += exp(belief[d]) * Kval;
					}
				}

				pixel++;

			}
		}
	}

	delete [] belief;

}



void crf::computeModelDistO(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (prm->timeseq!=1)
    return;

  if (fpvp->nO == 0)
    return; // nothing to do

  int nOpC = 0;
  if (fpvp->opflowlocal == 1)
  {
    nOpC = fpvp->nO/(prm->nD * fpvp->opflowlocalframes.size()); // there are (u,v) variables per frame except last frame for continuous features
  // but is more complicated for bin features (modeled similar to color uv.
  }

  int thetaid;
  int thetaido[2]={0,0};
  float indcost[2] = {0.0, 0.0};

  MRFEnergy* mrfe = (MRFEnergy*)mrf;

  float *belief = new float[prm->nD];
  int pixel = 0;

  // the last frame should not be counted towards the cost
  for (int z = 0; z < depth-1; z++)
  {
    for (int y = 0; y < height; y++)
    {
      uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);
      for (int x = 0; x < width; x++)
      {
        if (prm->crfpv.hidden == 1) {} // use all pixels
        else if (prm->unknown == 1 && dig[x]==prm->nD)         // unknown pixels do not count during learning parameters
        {
          pixel++;
          continue;
        }

        mrfe->getBelief(pixel,belief);

        for (int d=0; d < prm->nD ; d++)
        {
          int idx = frameExistsInFrameList(z);
          if (idx != -1)
          {
            float u = flowvector[index][z][0][y*width+x];
            float v = flowvector[index][z][1][y*width+x];
            if (0) // optical flow continuous
            {
              getIndicatorCostOpticalFlow(depth, z, d, fpvp->thetaO, thetaido, indcost, fpvp->opflowlocalframes.size(), idx);
              totalModelDistO[thetaido[0]] += exp(belief[d]) * u;
              totalModelDistO[thetaido[1]] += exp(belief[d]) * v;
            }
            else // binned optical flow features
            {
              getCostOpflowBinsuv(u, v, nOpC, d, fpvp->thetaO, thetaid, fpvp->opflowlocalframes.size(), idx);
              totalModelDistO[thetaid] += exp(belief[d]);
            }
          }
          // consider the costs for the theta param,
          // since they exist after differentiation
          // no use of the indcost since they are the theta values
        }
        pixel++;
      }
    }

  }

  delete [] belief;
}





// Note unknown region can currently only be used with fixed pairwise parameter
// ie when there is only one -v value.
// So currently there are no changes to computeModelDistV


// this will not learn parameters related to smoothed contrast sensitive gradients
void crf::computeModelDistV(int index)
{

  if (fpvp->nV == 0 && fpvp->nC == 0 && fpvp->nZ == 0)
    return; // nothing to do

  if (fpvp->csgrad == 0)
  {
    if (globalP.getOpflowConnection() == 1)
    {
      updateQuantizedPairwiseModelDistributionOpflow(index);
    } else {
      updateQuantizedPairwiseModelDistribution(index);
    }
  }
  else if (fpvp->csgrad == 1)
  {
    if (globalP.getOpflowConnection() == 1)
    {
      updateSmoothPairwiseModelDistributionOpflow(index);
    } else {
      updateSmoothPairwiseModelDistribution(index);
    }
  }
}


void crf::updateQuantizedPairwiseModelDistributionOpflow(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  assert(fpvp->gradVZ != 1);
  assert(globalP.getOpflowVolumeIndex() == index);

  if (fpvp->gradContext == 0)
  {
    std::cout << " Not implemented gradcontext = 0 \n";
  }

  MRFEnergy* mrfe = (MRFEnergy*)mrf;

  uchar *dig = 0;
  uchar *dig2 = 0;


  float totBelEqDisp, belNotEq;
  float *belEqDisp;
  float *belief = new float[prm->nD*prm->nD];

  std::map <std::pair <int, int>, std::pair <int, MRF::CostVal> >::iterator iter = globalP.edgeGlobalF[index].begin();

  for (; iter != globalP.edgeGlobalF[index].end(); iter++)
  {
    int node1 = iter->first.first;
    int node2 = iter->first.second;

    int optf = iter->second.first;
    int gradval = (int) iter->second.second;

    mrfe->getEdgeBelief( node1, node2, belief);

    int z1 = node1 / (width*height);
    int z2 = node2 / (width*height);

    int rem1 = node1 - z1*(width*height);
    int rem2 = node2 - z2*(width*height);

    int y1 = rem1 / width;
    int y2 = rem2 / width;

    int x1 = rem1 - y1*width;
    int x2 = rem2 - y2*width;

    dig = &gtdirImage[index][z1].Pixel(x1, y1, 0);
    dig2 = &gtdirImage[index][z2].Pixel(x2, y2, 0);

    if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
    else
    {


      for (int d1v=0; d1v<prm->nD; d1v++)
      {
        for (int d2v=0; d2v<prm->nD; d2v++)
        {
          if (fpvp->nV > 1)  // this makes sure it is not smooth constrast sensitive
          {
            int d1 = d1v;
            int d2 = d2v;

            if (d1>d2)
            {
              int swap = d2;
              d2 = d1;
              d1 = swap;
            }

            if (optf == 0) // spatial ie x or y
            {
              totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] +=  exp(belief[d1*prm->nD + d2]);
            } else { // optical flow connection
              if (fpvp->gradVZ == 2)
              {
                totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] +=  exp(belief[d1*prm->nD + d2]);
              }
              if (fpvp->gradVZ == 3)
              {
                totalModelDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] +=  exp(belief[d1*prm->nD + d2]);
              }
            }

          }
        }
      }

    }
  }

  delete[] belief;

}







void crf::updateQuantizedPairwiseModelDistribution(int index)
{

	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	MRFEnergy* mrfe = (MRFEnergy*)mrf;

	int depthrun = 0;

	/*
	if (fpvp->gradVZ == 1)
		depthrun = depth;
	else if (fpvp->gradVZ == 2 || fpvp->gradVZ == 3)
		depthrun = depth-1;
	 */

	depthrun = depth;

	uchar *di = 0;
	uchar *di2 = 0;
	uchar *di3 = 0;

  uchar *dig = 0;
  uchar *dig2 = 0;
  uchar *dig3 = 0;

	float totBelEqDisp, belNotEq;
	float *belEqDisp;
	float *belief = new float[prm->nD*prm->nD];


	for (int z=0; z < depthrun; z++)
	{
		for (int y = 0; y < height-1; y++)
		{
			di = &dirDisp[index][z].Pixel(0, y, 0);
			di2 = &dirDisp[index][z].Pixel(0, y+1, 0);

			if ((fpvp->gradVZ == 2 || fpvp->gradVZ == 3)&&(z<depthrun-1))
				di3 = &dirDisp[index][z+1].Pixel(0, y, 0);

			for (int x = 0; x < width-1; x++)
			{
				int gh, gv, gd;
				gh	= dirgrad[index][z].Pixel(x, y, 0);
				gv	= dirgrad[index][z].Pixel(x, y, 1);
				gd	= dirgrad[index][z].Pixel(x, y, 2);


				//CHECK need to make sure it has the gradContext option too
				// i am not fixing this gradOnly option
				if (fpvp->gradContext==0)
				{

					std::cout << " Not supposed to be working really check this !!! \n";
					exit(1);


				} else { //gradContext option


					if (fpvp->nV > 1)  // this makes sure it is not smooth constrast sensitive
					{

						// Horizontal Edge

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z].Pixel(x+1, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
            else
            {
              int p1, p2;

              p1 = z*(width*height) + y*width + x;
              p2 = z*(width*height) + y*width + x+1;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrfe->getEdgeBelief( p1, p2, belief);

              for (int d1v=0; d1v<prm->nD; d1v++)
              {
                for (int d2v=0; d2v<prm->nD; d2v++)
                {
                  int d1 = d1v;
                  int d2 = d2v;

                  if (d1>d2)
                  {
                    int swap = d2;
                    d2 = d1;
                    d1 = swap;
                  }

                  totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gh] +=  exp(belief[d1*prm->nD + d2]);

                }
              }
            }



						// Vertical Edge
            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z].Pixel(x, y+1, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
            else
            {

              int p1 = z*(width*height) + y*width + x;
              int p2 = z*(width*height) + (y+1)*width + x;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrfe->getEdgeBelief( p1, p2, belief);

              for (int d1v=0; d1v<prm->nD; d1v++)
              {
                for (int d2v=0; d2v<prm->nD; d2v++)
                {
                  int d1 = d1v;
                  int d2 = d2v;

                  if (d1>d2)
                  {
                    int swap = d2;
                    d2 = d1;
                    d1 = swap;
                  }
                  totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gv] +=  exp(belief[d1*prm->nD + d2]);

                }
              }
            }
					}


					if (fpvp->gradVZ == 2 && fpvp->nV > 1 && z < depthrun-1)
					{

						// Depth Edge
            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z+1].Pixel(x, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
            else
            {
              int p1 = z*(width*height) + y*width + x;
              int p2 = (z+1)*(width*height) + y*width + x;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrfe->getEdgeBelief( p1, p2, belief);

              for (int d1v=0; d1v<prm->nD; d1v++)
              {
                for (int d2v=0; d2v<prm->nD; d2v++)
                {
                  int d1 = d1v;
                  int d2 = d2v;

                  if (d1>d2)
                  {
                    int swap = d2;
                    d2 = d1;
                    d1 = swap;
                  }

                  totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gd] +=  exp(belief[d1*prm->nD + d2]);

                }
			        }
            }

					}


					if (fpvp->gradVZ == 3 && fpvp->nZ > 1 && z < depthrun-1)
					{

						// Depth Edge
            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z+1].Pixel(x, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
            else
            {

              int p1 = z*(width*height) + y*width + x;
              int p2 = (z+1)*(width*height) + y*width + x;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrfe->getEdgeBelief( p1, p2, belief);


              for (int d1v=0; d1v<prm->nD; d1v++)
              {
                for (int d2v=0; d2v<prm->nD; d2v++)
                {
                  int d1 = d1v;
                  int d2 = d2v;

                  if (d1>d2)
                  {
                    int swap = d2;
                    d2 = d1;
                    d1 = swap;
                  }

                  totalModelDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nGZ + gd] +=  exp(belief[d1*prm->nD + d2]);

                }
              }
			      }
					}
				}
			}
		}
	}

	delete[] belief;

}




void crf::updateSmoothPairwiseModelDistributionOpflow(int index)
{
  assert(0); // not coded

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  assert(fpvp->gradVZ != 1);
  assert(globalP.getOpflowVolumeIndex() == index);

  if (fpvp->gradContext == 0)
  {
    std::cout << " Not implemented gradcontext = 0 \n";
  }

  MRFEnergy* mrfe = (MRFEnergy*)mrf;

  uchar *dig = 0;
  uchar *dig2 = 0;


  float totBelEqDisp, belNotEq;
  float *belEqDisp;
  float *belief = new float[prm->nD*prm->nD];

  std::map <std::pair <int, int>, std::pair <int, MRF::CostVal> >::iterator iter = globalP.edgeGlobalF[index].begin();

  for (; iter != globalP.edgeGlobalF[index].end(); iter++)
  {
    int node1 = iter->first.first;
    int node2 = iter->first.second;

    int optf = iter->second.first;
    int gradval = (int) iter->second.second;

    mrfe->getEdgeBelief( node1, node2, belief);

    int z1 = node1 / (width*height);
    int z2 = node2 / (width*height);

    int rem1 = node1 - z1*(width*height);
    int rem2 = node2 - z2*(width*height);

    int y1 = rem1 / width;
    int y2 = rem2 / width;

    int x1 = rem1 - y1*width;
    int x2 = rem2 - y2*width;

    dig = &gtdirImage[index][z1].Pixel(x1, y1, 0);
    dig2 = &gtdirImage[index][z2].Pixel(x2, y2, 0);

    if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
    else
    {


      for (int d1v=0; d1v<prm->nD; d1v++)
      {
        for (int d2v=0; d2v<prm->nD; d2v++)
        {
          if (fpvp->nV > 1)  // this makes sure it is not smooth constrast sensitive
          {
            int d1 = d1v;
            int d2 = d2v;

            if (d1>d2)
            {
              int swap = d2;
              d2 = d1;
              d1 = swap;
            }

            if (optf == 0) // spatial ie x or y
            {
              totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] +=  exp(belief[d1*prm->nD + d2]);
            } else { // optical flow connection
              if (fpvp->gradVZ == 2)
              {
                totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] +=  exp(belief[d1*prm->nD + d2]);
              }
              if (fpvp->gradVZ == 3)
              {
                totalModelDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] +=  exp(belief[d1*prm->nD + d2]);
              }
            }

          }
        }
      }

    }
  }

  delete[] belief;

}







void crf::updateSmoothPairwiseModelDistribution(int index)
{

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  MRFEnergy* mrfe = (MRFEnergy*)mrf;

  int depthrun = 0;

  /*
  if (fpvp->gradVZ == 1)
    depthrun = depth;
  else if (fpvp->gradVZ == 2 || fpvp->gradVZ == 3)
    depthrun = depth-1;
   */

  depthrun = depth;

  uchar *di = 0;
  uchar *di2 = 0;
  uchar *di3 = 0;

  uchar *dig = 0;
  uchar *dig2 = 0;
  uchar *dig3 = 0;

  float totBelEqDisp, belNotEq;
  float *belEqDisp;
  float *belief = new float[prm->nD*prm->nD];


  for (int z=0; z < depthrun; z++)
  {
    for (int y = 0; y < height-1; y++)
    {
      di = &dirDisp[index][z].Pixel(0, y, 0);
      di2 = &dirDisp[index][z].Pixel(0, y+1, 0);

      if ((fpvp->gradVZ == 2 || fpvp->gradVZ == 3)&&(z<depthrun-1))
        di3 = &dirDisp[index][z+1].Pixel(0, y, 0);

      for (int x = 0; x < width-1; x++)
      {
        double valh, valv, vald;
        CvScalar s;
        s=cvGet2D(dircsgrad[index][z],y,x); // get the (i,j) pixel value
        valh = exp(s.val[0]);
        valv = exp(s.val[1]);
        vald = exp(s.val[2]);


        //CHECK need to make sure it has the gradContext option too
        // i am not fixing this gradOnly option
        if (fpvp->gradContext==0)
        {
          std::cout << " Not supposed to be working really check this !!! \n";
          exit(1);

        } else { //gradContext option

          if (fpvp->nV == 1)  // this makes sure it is not smooth constrast sensitive
          {

            // Horizontal Edge

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z].Pixel(x+1, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
            else
            {
              int p1, p2;

              p1 = z*(width*height) + y*width + x;
              p2 = z*(width*height) + y*width + x+1;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrfe->getEdgeBelief( p1, p2, belief);

              for (int d1v=0; d1v<prm->nD; d1v++)
              {
                for (int d2v=0; d2v<prm->nD; d2v++)
                {
                  int d1 = d1v;
                  int d2 = d2v;

                  if (prm->unknown == 1 && (d1==prm->nD || d2==prm->nD) && prm->crfpv.hidden == 0) {}
                  else if (d1 == d2)
                  {} // since [yp != yq] in Vpq
                  else
                  {
                    totalModelDistV[0] += valh * exp(belief[d1*prm->nD + d2]);
                  }

                }
              }
            }



            // Vertical Edge
            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z].Pixel(x, y+1, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
            else
            {

              int p1 = z*(width*height) + y*width + x;
              int p2 = z*(width*height) + (y+1)*width + x;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrfe->getEdgeBelief( p1, p2, belief);

              for (int d1v=0; d1v<prm->nD; d1v++)
              {
                for (int d2v=0; d2v<prm->nD; d2v++)
                {
                  int d1 = d1v;
                  int d2 = d2v;

                  if (prm->unknown == 1 && (d1==prm->nD || d2==prm->nD) && prm->crfpv.hidden == 0) {}
                  else if (d1 == d2)
                  {} // since [yp != yq] in Vpq
                  else
                  {
                    totalModelDistV[0] += valv * exp(belief[d1*prm->nD + d2]);
                  }

                }
              }
            }
          }


          if (fpvp->gradVZ == 2 && fpvp->nV == 1 && z < depthrun-1)
          {

            // Depth Edge
            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z+1].Pixel(x, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
            else
            {
              int p1 = z*(width*height) + y*width + x;
              int p2 = (z+1)*(width*height) + y*width + x;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrfe->getEdgeBelief( p1, p2, belief);

              for (int d1v=0; d1v<prm->nD; d1v++)
              {
                for (int d2v=0; d2v<prm->nD; d2v++)
                {
                  int d1 = d1v;
                  int d2 = d2v;

                  if (prm->unknown == 1 && (d1==prm->nD || d2==prm->nD) && prm->crfpv.hidden == 0) {}
                  else if (d1 == d2)
                  {} // since [yp != yq] in Vpq
                  else
                  {
                    totalModelDistV[0] += vald * exp(belief[d1*prm->nD + d2]);
                  }
                }
              }
            }

          }


          if (fpvp->gradVZ == 3 && fpvp->nZ == 1 && z < depthrun-1)
          {

            // Depth Edge
            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z+1].Pixel(x, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD) && prm->crfpv.hidden == 0) {}
            else
            {

              int p1 = z*(width*height) + y*width + x;
              int p2 = (z+1)*(width*height) + y*width + x;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrfe->getEdgeBelief( p1, p2, belief);


              for (int d1v=0; d1v<prm->nD; d1v++)
              {
                for (int d2v=0; d2v<prm->nD; d2v++)
                {
                  int d1 = d1v;
                  int d2 = d2v;

                  if (prm->unknown == 1 && (d1==prm->nD || d2==prm->nD) && prm->crfpv.hidden == 0) {}
                  else if (d1 == d2)
                  {} // since [yp != yq] in Vpq
                  else
                  {
                    totalModelDistZ[0] += vald * exp(belief[d1*prm->nD + d2]);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  delete[] belief;

}



// A lot of this code is similar to the model distribution code
// except that the pixels to be used are the hidden or unknown nodes.
// and that we update the empirical distribution instead of the model distribution
// it might be better to have a flag in the model_distribution functions to select to use the unkwown pixels
// and empirical dist when we are using hidden crf

void crf::computeHiddenDist(int index)
{

  int mpe = -1;

    if (prm->iopv.testDirIndexV[index]==1 && (prm->crfpv.inferencerTest == USE_GC || prm->crfpv.inferencerTest == USE_MP))
    {
      mpe = 1;
    }

    if (prm->iopv.testDirIndexV[index]==0 && (prm->crfpv.inferencer == USE_GC || prm->crfpv.inferencer == USE_MP))
    {
      mpe = 1;
    }

    if (prm->iopv.testDirIndexV[index]==1 && (prm->crfpv.inferencerTest == USE_MF ))
    {
      mpe = 0;
    }

    if (prm->iopv.testDirIndexV[index]==0 && (prm->crfpv.inferencer == USE_MF) )
    {
      mpe = 0;
    }


    if (mpe == -1)
    {
      std::cout << " Need to initialize variable mpe" << std::endl;
      exit(1);
    }





    // this uses point estimates
    if (mpe == 1)
    {

    if (fpvp->nU > 0 && fpvp->intensity==1)
      computeHiddenDistUmpe(index);

    if (fpvp->nA > 0 && (fpvp->app==1 || fpvp->app==2))
      computeHiddenDistAmpe(index);

    if (fpvp->nH > 0 && fpvp->hog==1)
      computeHiddenDistHmpe(index);

    if (fpvp->nL > 0 && (fpvp->loc==1 || fpvp->loc==6 || fpvp->loc==8 || fpvp->loc==7 || fpvp->loc==9))
      computeHiddenDistLmpe(index);

    if (fpvp->nM > 0 && fpvp->mot==1)
      computeHiddenDistMmpe(index);

    if (fpvp->klrpv.klr==1)
      computeHiddenDistWmpe(index);

    if (fpvp->nO > 0 && fpvp->opflowlocal==1)
      computeHiddenDistOmpe(index);

    // when nV == 1, we are not learning parameters (because constrast sensitive smoothing)
    // but in the function we still need to check again because nZ could be 1 , but not nV and so on
    if (fpvp->nV > 0 || fpvp->nC > 0 || fpvp-> nZ > 0)
      computeHiddenDistVmpe(index);

    if (fpvp->nB > 0 && fpvp->bias==1)
      computeHiddenDistBmpe(index);

    if (fpvp->nE > 0 && fpvp->volume==1)
        computeHiddenDistEmpe(index);

    }

    // this uses mrf probability values
  if (mpe == 0)
  {
    if (fpvp->nU > 0 && fpvp->intensity==1)
      computeHiddenDistU(index);

    if (fpvp->nA > 0 && (fpvp->app==1 || fpvp->app==2))
      computeHiddenDistA(index);

    if (fpvp->nH > 0 && fpvp->hog==1)
      computeHiddenDistH(index);

    if (fpvp->nL > 0 && (fpvp->loc==1 || fpvp->loc==6 || fpvp->loc==8 || fpvp->loc==7 || fpvp->loc==9))
      computeHiddenDistL(index);

    if (fpvp->nM > 0 && fpvp->mot==1)
      computeHiddenDistM(index);

    if (fpvp->klrpv.klr==1)
      computeHiddenDistW(index);

    if (fpvp->nO > 0 && fpvp->opflowlocal==1)
      computeHiddenDistO(index);

    // when nV == 1, we are not learning parameters (because constrast sensitive smoothing)
    // but in the function we still need to check again because nZ could be 1 , but not nV and so on
    if (fpvp->nV > 1 || fpvp->nC > 0 || fpvp-> nZ > 1)
      computeHiddenDistV(index);

    if (fpvp->nB > 0 && fpvp->bias==1)
      computeHiddenDistB(index);

    if (fpvp->nE > 0 && fpvp->volume==1)
        computeHiddenDistE(index);

  }

}












// mpe functions

// grayscale and // RGB or Yuv or Lab
// does not current have interactive mode
void crf::computeHiddenDistUmpe(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nU == 0)
    return; // nothing to do

  int nbins = double(fpvp->nU)/prm->nD;
  double binner = fpvp->rangeI/nbins;

  int thetaid;
  int thetaidrgb[3]={0,0,0};
  int thetaidYuv[2]={0,0};


  uchar *pix1t = 0;
  unsigned short pix1 = 0;

  for (int z = 0; z < depth; z++)
  {
    for (int y = 0; y < height; y++)
    {
      uchar *di = &dirDisp[index][z].Pixel(0, y, 0);

      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

      for (int x = 0; x < width; x++)
      {
        if (!(prm->unknown == 1 && dig[x]==prm->nD)) // unknown pixels do not count during learning parameters
          continue;

        if (prm->timeseq==1)
        {
          pix1t = &indirImaget[index][z].Pixel(x, y, 0);
              // getCostIntensityDiscBins(pix1t, nbins, di[x], thetaU, thetaidrgb);
          if (prm->featurepv.onlyUV == 1)
            getCostIntensityDiscBinsuv(pix1t, nbins, di[x], fpvp->thetaU, thetaidYuv);
          else
            getCostIntensityDiscBinsYuv(pix1t, nbins, di[x], fpvp->thetaU, thetaidYuv);

          if (prm->featurepv.onlyUV == 0)
            totalEmpirDistU[thetaidYuv[0]]++;
          totalEmpirDistU[thetaidYuv[1]]++;

        }
        else if (prm->timeseq==0)
        {
          pix1 = indirImage[index][z].Pixel(x, y, 0);
          getCostIntensityDiscBins(pix1, nbins, di[x], fpvp->thetaU, thetaid);
          totalEmpirDistU[thetaid] ++;
        }

        // for us the costs are the theta param, but in differentiation, they vanish,
        // so its just whether the feature is enabled or not.
        // move costs pointer to next pixel
      }
    }
  }
}


void crf::computeHiddenDistAmpe(int index)
{

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nA == 0)
    return; // nothing to do

  int nT = fpvp->nA/prm->nD; //this gives centers per class
  // this will work for app==2 and app==1 because in app==2, we have nD*the number of total clusters as the number of parameters
  // while in app==1 it is nD*number of clusters of that class (only that clusters are also = number of params)

    int thetaid;


  for (int z=0; z < depth; z++)
  {
    for (int y = 0; y < height; y++)
    {
      uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

      for (int x = 0; x < width; x++)
      {
        if (!(prm->unknown == 1 && dig[x]==prm->nD))         // unknown pixels do not count during learning parameters
          continue;

        int d = di[x];
        int dm = d;
                if (fpvp->app==2)
                  dm=0;

              if (((int)appdirImage[index][z][dm].Pixel(x, y, 0) >= nT && fpvp->app==1) || ((int)appdirImage[index][z][0].Pixel(x, y, 0) >= fpvp->nA && fpvp->app==2))
              {
                continue;
              }
                getCostAppDiscPatch((int)appdirImage[index][z][dm].Pixel(x, y, 0), appclass[dm]->getPatchSize(), x, y, width, height, fpvp->thetaA, d, nT, thetaid);

                // assume the patch size is the same for all classes
                int patchSize = appclass[dm]->getPatchSize();

                if (x>(int)patchSize/2 && x<width-patchSize/2 && y>(int)patchSize/2 && y<height-patchSize/2)
                  totalEmpirDistA[thetaid] ++;

      }
    }
  }
}




void crf::computeHiddenDistHmpe(int index)
{

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nH == 0)
    return; // nothing to do

  int nHpC = fpvp->nH/prm->nD; // HoG vocabulary size

    int thetaid;

    for (int z = 0; z < depth; z++)
    {
      for (int y = 0; y < height; y++)
      {
        uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
        uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

        for (int x = 0; x < width; x++)
        {
          if (!(prm->unknown == 1 && dig[x]==prm->nD))         // unknown pixels do not count during learning parameters
            continue;

          unsigned short hogval = hogdirImage[index][z].Pixel(x, y, 0);
          getCostHoGDiscBins((int) hogval, di[x], nHpC, fpvp->thetaH, thetaid);
          //getCostHoGDiscBinsTemp((int) hogval, di[x], fpvp->thetaH, thetaid);

          // for us the costs are the theta param, but in differentiation, they vanish,
          // so its just whether the feature is enabled or not.
          // move costs pointer to next pixel
          if (thetaid >= 0)
            totalEmpirDistH[thetaid] ++;

        }
      }
    }
}

void crf::computeHiddenDistBmpe(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nB == 0)
    return; // nothing to do

  int thetaid;

  for (int z = 0; z < depth; z++)
  {
    for (int y = 0; y < height; y++)
    {
      uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

      for (int x = 0; x < width; x++)
      {
        if (!(prm->unknown == 1 && dig[x]==prm->nD))         // unknown pixels do not count during learning parameters
          continue;

        getCostBias(di[x], fpvp->thetaB, thetaid);

        // for us the costs are the theta param, but in differentiation, they vanish,
        // so its just whether the feature is enabled or not.
        // move costs pointer to next pixel
        if (thetaid >= 0)
          totalEmpirDistB[thetaid] ++;

      }
    }
  }
}

void crf::computeHiddenDistEmpe(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nE == 0)
    return; // nothing to do

  int thetaid;

  for (int z = 0; z < depth; z++)
  {
    for (int y = 0; y < height; y++)
    {
      uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

      for (int x = 0; x < width; x++)
      {
        if (!(prm->unknown == 1 && dig[x]==prm->nD))         // unknown pixels do not count during learning parameters
          continue;

        getCostInverseClassSize(di[x], inverseClassSize, fpvp->thetaE, thetaid);

        // for us the costs are the theta param, but in differentiation, they vanish,
        // so its just whether the feature is enabled or not.
        // move costs pointer to next pixel
        if (thetaid >= 0)
          totalEmpirDistE[thetaid] += inverseClassSize[di[x]];

      }
    }
  }
}

void crf::computeHiddenDistMmpe(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nM == 0)
    return; // nothing to do

  int nMpC = fpvp->nM/prm->nD; // Mot vocabulary size

    int thetaid;

    for (int z = 0; z < depth; z++)
    {
      for (int y = 0; y < height; y++)
      {
        uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
        uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

        for (int x = 0; x < width; x++)
        {
          if (!(prm->unknown == 1 && dig[x]==prm->nD))         // unknown pixels do not count during learning parameters
            continue;

          unsigned short motval = motdirImage[index][z].Pixel(x, y, 0);
          getCostmotDiscBins((int) motval, di[x], nMpC, fpvp->thetaM, thetaid);

          totalEmpirDistM[thetaid] ++;

        }
      }
    }
}


void crf::computeHiddenDistLmpe(int index)
{

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nL == 0)
    return; // nothing to do

  int nLpC = fpvp->nL/prm->nD;  // this variable will be useful only for cube setting (loc==1) or loc==9

  int thetaid;

  int zslice = 0;
  if (timeseq == 0 && (fpvp->loc != 0 || fpvp->klrpv.locKlr ==2))
    zslice = startSliceNo[index];

    for (int z = 0; z < depth; z++)
    {

      for (int y = 0; y < height; y++)
      {
        uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
        uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

        for (int x = 0; x < width; x++)
        {
          if (!(prm->unknown == 1 && dig[x]==prm->nD))         // unknown pixels do not count during learning parameters
            continue;

          if (fpvp->loc==1)
          {
            /*
            unsigned short locpix = locdirImage[index][di[x]][z].Pixel(x,y,0);
            getCostLocationDiscCube(locpix, x, y, z, fpvp->loccubex, fpvp->loccubey, fpvp->loccubez, numTrainingPats, nLpC, fpvp->thetaL, width, height, depth, di[x], thetaid);
            double locpixd = (double)locpix/ (fpvp->loccubex * fpvp->loccubey * fpvp->loccubez * numTrainingPats);

            totalEmpirDistL[thetaid] += (1-locpixd);
            */

            std::cout << " Need to fix this - location " << std::endl;
            exit(1);


          } else if (fpvp->loc==6 || fpvp->loc==8) { // CRF hard assignment and soft

            getCostLocationDiscGaussian(LocationMVclass, z, zslice, x, y, fpvp->thetaL, di[x], fpvp->loc, fpvp->generative, thetaid);
            totalEmpirDistL[thetaid]++;

          } else if (fpvp->loc==9 ||  (fpvp->loc==7 && timeseq == 1)) {  //CRF soft or hard assignment redefined

            getCostLocationDiscGaussian2(LocationMVclass, z, zslice, x, y, fpvp->thetaL, di[x], fpvp->loc, fpvp->generative, prm->nD, thetaid);
            totalEmpirDistL[thetaid]++;

          } else if  (fpvp->loc==7 && timeseq == 0) {  //CRF soft or hard assignment redefined

            getCostLocationDiscGaussian3(index, z, x, y, di[x], fpvp->generative, prm->nD, thetaid);
            totalEmpirDistL[thetaid]++;
          }



        }
      }
    }
}


void crf::computeHiddenDistWmpe(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->klrpv.nW == 0)
    return; // nothing to do

  int dimW = wparam.size2();
  int Ntr = wparam.size1();


  int zslice = 0;

  if (timeseq == 0 && (fpvp->loc != 0 || fpvp->klrpv.locKlr ==2))
    zslice = startSliceNo[index];


  uchar *pix1t = 0;
  unsigned short pix1 = 0;

  ublas::vector<DBL_TYPE> Xtest(Xtrain.size2());
  ublas::vector<DBL_TYPE> f;


  for (int z = 0; z < depth; z++)
  {

    matrix<DBL_TYPE> hogdescs, appdescs;

    if (fpvp->klrpv.klr==1)
      readklrdescfiles(hogdescs, appdescs, index, z+zslice);

    for (int y = 0; y < height; y++)
    {
      uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

      for (int x = 0; x < width; x++)
      {
        if (!(prm->unknown == 1 && dig[x]==prm->nD))         // unknown pixels do not count during learning parameters
          continue;

        getDataCostKlr(f, Xtest, index, z, x, y, zslice, hogdescs, appdescs, width, height);

        for (int j=0; j<Ntr; j++)
        {
          double Kval = ivm->classify_klrexp(&Xtest, j);

          int thetaid = j*dimW + di[x];
          totalEmpirDistW[thetaid] += Kval;
        }

      }
    }
  }
}



void crf::computeHiddenDistOmpe(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (prm->timeseq!=1)
    return;

  if (fpvp->nO == 0)
    return; // nothing to do

  int nOpC = 0;
  if (fpvp->opflowlocal == 1)
  {
    nOpC = fpvp->nO/(prm->nD * fpvp->opflowlocalframes.size()); // there are (u,v) variables per frame except last frame for continuous features
  // but is more complicated for bin features (Hiddened similar to color uv.
  }

  int thetaid;
  int thetaido[2]={0,0};
  float indcost[2] = {0.0, 0.0};

  // the last frame should not be counted towards the cost
  for (int z = 0; z < depth-1; z++)
  {
    for (int y = 0; y < height; y++)
    {
      uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);
      for (int x = 0; x < width; x++)
      {
        if (!(prm->unknown == 1 && dig[x]==prm->nD))         // unknown pixels do not count during learning parameters
          continue;

        int idx = frameExistsInFrameList(z);
        if (idx != -1)
        {
          float u = flowvector[index][z][0][y*width+x];
          float v = flowvector[index][z][1][y*width+x];
          if (0) // optical flow continuous
          {
            getIndicatorCostOpticalFlow(depth, z, di[x], fpvp->thetaO, thetaido, indcost, fpvp->opflowlocalframes.size(), idx);
            totalEmpirDistO[thetaido[0]] +=  u;
            totalEmpirDistO[thetaido[1]] +=  v;
          }
          else // binned optical flow features
          {
            getCostOpflowBinsuv(u, v, nOpC, di[x], fpvp->thetaO, thetaid, fpvp->opflowlocalframes.size(), idx);
            totalEmpirDistO[thetaid] ++;
          }
        }
        // consider the costs for the theta param,
        // since they exist after differentiation
        // no use of the indcost since they are the theta values

      }
    }

  }
}





void crf::updateQuantizedPairwiseHiddenDistributionmpeOpflow(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  assert(fpvp->gradVZ != 1);
  assert(globalP.getOpflowVolumeIndex() == index);

  if (fpvp->gradContext == 0)
  {
    std::cout << " Not implemented gradcontext = 0 \n";
  }

  uchar *di = 0;
  uchar *di2 = 0;

  uchar *dig = 0;
  uchar *dig2 = 0;


  std::map <std::pair <int, int>, std::pair <int, MRF::CostVal> >::iterator iter = globalP.edgeGlobalF[index].begin();

  for (; iter != globalP.edgeGlobalF[index].end(); iter++)
  {
    int node1 = iter->first.first;
    int node2 = iter->first.second;

    int optf = iter->second.first;
    int gradval = (int) iter->second.second;

    int z1 = node1 / (width*height);
    int z2 = node2 / (width*height);

    int rem1 = node1 - z1*(width*height);
    int rem2 = node2 - z2*(width*height);

    int y1 = rem1 / width;
    int y2 = rem2 / width;

    int x1 = rem1 - y1*width;
    int x2 = rem2 - y2*width;

    di = &dirDisp[index][z1].Pixel(x1, y1, 0);
    di2 = &dirDisp[index][z2].Pixel(x2, y2, 0);

    dig = &gtdirImage[index][z1].Pixel(x1, y1, 0);
    dig2 = &gtdirImage[index][z2].Pixel(x2, y2, 0);

    if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
    {

      if (fpvp->nV > 1)  // this makes sure it is not smooth constrast sensitive
      {
        int d1 = (int) di[0];
        int d2 = (int) di2[0];

        if (d1>d2)
        {
          int swap = d2;
          d2 = d1;
          d1 = swap;
        }

        if (optf == 0) // spatial ie x or y
        {
          totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] ++;
        } else { // optical flow connection
          if (fpvp->gradVZ == 2)
          {
            totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] ++;
          }
          if (fpvp->gradVZ == 3)
          {
            totalEmpirDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] ++;
          }
        }

      }
    }
  }
}


void crf::updateSmoothPairwiseHiddenDistributionmpeOpflow(int index)
{
  // we currently don't update V and Z if they are single
  if ((fpvp->nV == 1 && fpvp->nZ == 1) || (fpvp->nV + fpvp->nZ == 1))
    return;

  assert(0); // might need to double check validity of code in this function

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  assert(fpvp->gradVZ != 1);
  assert(globalP.getOpflowVolumeIndex() == index);

  if (fpvp->gradContext == 0)
  {
    std::cout << " Not implemented gradcontext = 0 \n";
  }

  uchar *di = 0;
  uchar *di2 = 0;

  uchar *dig = 0;
  uchar *dig2 = 0;


  std::map <std::pair <int, int>, std::pair <int, MRF::CostVal> >::iterator iter = globalP.edgeGlobalF[index].begin();

  for (; iter != globalP.edgeGlobalF[index].end(); iter++)
  {
    int node1 = iter->first.first;
    int node2 = iter->first.second;

    int optf = iter->second.first;
    int gradval = (int) iter->second.second;

    int z1 = node1 / (width*height);
    int z2 = node2 / (width*height);

    int rem1 = node1 - z1*(width*height);
    int rem2 = node2 - z2*(width*height);

    int y1 = rem1 / width;
    int y2 = rem2 / width;

    int x1 = rem1 - y1*width;
    int x2 = rem2 - y2*width;

    di = &dirDisp[index][z1].Pixel(x1, y1, 0);
    di2 = &dirDisp[index][z2].Pixel(x2, y2, 0);

    dig = &gtdirImage[index][z1].Pixel(x1, y1, 0);
    dig2 = &gtdirImage[index][z2].Pixel(x2, y2, 0);

    if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
    {

      if (fpvp->nV > 1 || fpvp->nZ > 1)  // this makes sure it is not smooth constrast sensitive
      {
        int d1 = (int) di[0];
        int d2 = (int) di2[0];

        if (d1>d2)
        {
          int swap = d2;
          d2 = d1;
          d1 = swap;
        }

        if (optf == 0) // spatial ie x or y
        {
          totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] ++;
        } else { // optical flow connection
          if (fpvp->gradVZ == 2)
          {
            totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] ++;
          }
          if (fpvp->gradVZ == 3)
          {
            totalEmpirDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] ++;
          }
        }

      }
    }
  }
}






// Note unknown region can currently only be used with fixed pairwise parameter
// ie when there is only one -v value.
// So currently there are no changes to computeHiddenDistVmpe

// this will not learn parameters related to smoothed contrast sensitive gradients
void crf::computeHiddenDistVmpe(int index)
{

  if (fpvp->nV == 0 && fpvp->nC == 0 && fpvp->nZ == 0)
    return; // nothing to do

  if (fpvp->csgrad == 0)
  {
    if (globalP.getOpflowConnection() == 1)
    {
      updateQuantizedPairwiseHiddenDistributionmpeOpflow(index);
    } else {
      updateQuantizedPairwiseHiddenDistributionmpe(index);
    }
  }
  else if (fpvp->csgrad == 1)
  {
    if (globalP.getOpflowConnection() == 1)
    {
      updateSmoothPairwiseHiddenDistributionmpeOpflow(index);
    } else {
      updateSmoothPairwiseHiddenDistributionmpe(index);
    }
  }

}



void crf::updateQuantizedPairwiseHiddenDistributionmpe(int index)
{

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  int depthrun = 0;

  /*
  if (fpvp->gradVZ == 1)
    depthrun = depth;
  else if (fpvp->gradVZ == 2 || fpvp->gradVZ == 3)
    depthrun = depth-1;
   */
  depthrun = depth;

  uchar *di = 0;
  uchar *di2 = 0;
  uchar *di3 = 0;

  uchar *dig = 0;
  uchar *dig2 = 0;
  uchar *dig3 = 0;


  for (int z=0; z < depthrun; z++)
  {
    for (int y = 0; y < height-1; y++)
    {
      di = &dirDisp[index][z].Pixel(0, y, 0);
      di2 = &dirDisp[index][z].Pixel(0, y+1, 0);

      if ((fpvp->gradVZ == 2 || fpvp->gradVZ == 3)&&(z<depthrun-1))
        di3 = &dirDisp[index][z+1].Pixel(0, y, 0);


      for (int x = 0; x < width-1; x++)
      {
        int gh, gv, gd;
        gh  = dirgrad[index][z].Pixel(x, y, 0);
        gv  = dirgrad[index][z].Pixel(x, y, 1);
        gd  = dirgrad[index][z].Pixel(x, y, 2);



        //CHECK need to make sure it has the gradContext option too
        // i am not fixing this gradOnly option
        if (fpvp->gradContext==0)
        {

          std::cout << " Not supposed to be working really check this !!! \n";
          exit(1);

          // Horizontal Edge
          if (di[x] != di[x+1])   // this is the border of the parts
            totalEmpirDistC[gh]++;

          // Vertical Edge
          if (di[x] != di2[x]) // this is the border of the parts
            totalEmpirDistC[gv]++;

          if (fpvp->gradVZ == 2)
          {
            // depth edge
            if (di[x] != di3[x]) // this is the border of the parts
              totalEmpirDistC[gd]++;
          }

        } else { //gradContext option

          if (fpvp->nV > 1)  // this makes sure it is not smooth constrast sensitive
          {

            // Horizontal Edge
            int d1 = (int) di[x];
            int d2 = (int) di[x+1];

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z].Pixel(x+1, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
            {
              if (d1>d2)
              {
                int swap = d2;
                d2 = d1;
                d1 = swap;
              }

              totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gh] ++;
            }



            // Vertical Edge
            d1 = (int) di[x];
            d2 = (int) di2[x];

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z].Pixel(x, y+1, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
            {
              if (d1>d2)
              {
                int swap = d2;
                d2 = d1;
                d1 = swap;
              }

              totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gv] ++;

            }

          }


          if (fpvp->gradVZ == 2 && fpvp->nV > 1 && z < depthrun-1)
          {

            // Depth Edge
            int d1 = (int) di[x];
            int d2 = (int) di3[x];

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z+1].Pixel(x, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
            {
              if (d1>d2)
              {
                int swap = d2;
                d2 = d1;
                d1 = swap;
              }

              totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gd] ++;
            }

          }


          if (fpvp->gradVZ == 3 && fpvp->nZ > 1 && z < depthrun-1)
          {

            // Depth Edge
            int d1 = (int) di[x];
            int d2 = (int) di3[x];

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z+1].Pixel(x, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
            {
              if (d1>d2)
              {
                int swap = d2;
                d2 = d1;
                d1 = swap;
              }
              totalEmpirDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nGZ + gd] ++;

            }
          }
        }
      }
    }
  }
}



void crf::updateSmoothPairwiseHiddenDistributionmpe(int index)
{

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  int depthrun = 0;

  /*
  if (fpvp->gradVZ == 1)
    depthrun = depth;
  else if (fpvp->gradVZ == 2 || fpvp->gradVZ == 3)
    depthrun = depth-1;
   */

  depthrun = depth;

  uchar *di = 0;
  uchar *di2 = 0;
  uchar *di3 = 0;

  uchar *dig = 0;
  uchar *dig2 = 0;
  uchar *dig3 = 0;


  for (int z=0; z < depthrun; z++)
  {
    for (int y = 0; y < height-1; y++)
    {
      di = &dirDisp[index][z].Pixel(0, y, 0);
      di2 = &dirDisp[index][z].Pixel(0, y+1, 0);

      if ((fpvp->gradVZ == 2 || fpvp->gradVZ == 3)  && (z<depth-1))
        di3 = &dirDisp[index][z+1].Pixel(0, y, 0);


      for (int x = 0; x < width-1; x++)
      {
        double valh, valv, vald;
        CvScalar s;
        s=cvGet2D(dircsgrad[index][z],y,x); // get the (i,j) pixel value
        valh = exp(s.val[0]);
        valv = exp(s.val[1]);
        vald = exp(s.val[2]);

        //CHECK need to make sure it has the gradContext option too
        // i am not fixing this gradOnly option
        if (fpvp->gradContext==0)
        {
          std::cout << " Not supposed to be working really check this !!! \n";
          exit(1);

        } else { //gradContext option

          if (fpvp->nV == 1)  // this makes sure it is smooth constrast sensitive
          {

            // Horizontal Edge
            int d1 = (int) di[x];
            int d2 = (int) di[x+1];

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z].Pixel(x+1, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
            {
              if (d1 == d2)
              {} // since [yp != yq] in Vpq
              else
              {
                totalEmpirDistV[0] += valh;
              }
            }

            // Vertical Edge
            d1 = (int) di[x];
            d2 = (int) di2[x];

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z].Pixel(x, y+1, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
            {
              if (d1 == d2)
              {} // since [yp != yq] in Vpq
              else
              {
                totalEmpirDistV[0] += valv;
              }
            }

          }


          if (fpvp->gradVZ == 2 && fpvp->nV == 1 && z < depthrun-1)
          {

            // Depth Edge
            int d1 = (int) di[x];
            int d2 = (int) di3[x];

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z+1].Pixel(x, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
            {
              if (d1 == d2)
              {} // since [yp != yq] in Vpq
              else
              {
                totalEmpirDistV[0] += vald;
              }
            }
          }


          if (fpvp->gradVZ == 3 && fpvp->nZ == 1 && z < depthrun-1)
          {

            // Depth Edge
            int d1 = (int) di[x];
            int d2 = (int) di3[x];

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z+1].Pixel(x, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
            {
              if (d1 == d2)
              {} // since [yp != yq] in Vpq
              else
              {
                totalEmpirDistZ[0] += vald;
              }
            }
          }
        }
      }
    }
  }
}




// grayscale and // RGB or Yuv or Lab
// does not current have interactive mode
void crf::computeHiddenDistU(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nU == 0)
    return; // nothing to do

  int nbins = double(fpvp->nU)/prm->nD;
  double binner = fpvp->rangeI/nbins;

  int thetaid;
    int thetaidrgb[3]={0,0,0};
    int thetaidYuv[2]={0,0};


  uchar *pix1t = 0;
  unsigned short pix1 = 0;

  MRFEnergy* mrfe = (MRFEnergy*)mrf;

    float *belief = new float[prm->nD];
    int pixel=0;

  for (int z = 0; z < depth; z++)
  {
    for (int y = 0; y < height; y++)
    {
      uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

      for (int x = 0; x < width; x++)
      {
        if (!(prm->unknown == 1 && dig[x]==prm->nD))         // unknown pixels do not count during learning parameters
          continue;

        mrfe->getBelief(pixel,belief);

        for (int d=0; d<prm->nD; d++)
        {

          if (prm->timeseq==1)
          {
            pix1t = &indirImaget[index][z].Pixel(x, y, 0);
            // getCostIntensityDiscBins(pix1t, nbins, d, thetaU, thetaidrgb);
            if (prm->featurepv.onlyUV == 1)
              getCostIntensityDiscBinsuv(pix1t, nbins, d, fpvp->thetaU, thetaidYuv);
            else
              getCostIntensityDiscBinsYuv(pix1t, nbins, d, fpvp->thetaU, thetaidYuv);

            if (prm->featurepv.onlyUV == 0)
              totalEmpirDistU[thetaidYuv[0]] += exp(belief[d]);
            totalEmpirDistU[thetaidYuv[1]] += exp(belief[d]);

          }
          else if (prm->timeseq==0)
          {
            pix1 = indirImage[index][z].Pixel(x, y, 0);
            getCostIntensityDiscBins(pix1, nbins, d, fpvp->thetaU, thetaid);
            totalEmpirDistU[thetaid] += exp(belief[d]);
          }

        }

        // for us the costs are the theta param, but in differentiation, they vanish,
        // so its just whether the feature is enabled or not.
        // move costs pointer to next pixel

        pixel++;

      }
    }

  }

  delete [] belief;

}


void crf::computeHiddenDistA(int index)
{

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nA == 0)
    return; // nothing to do

  int nT = fpvp->nA/prm->nD; //this gives centers per class
  // this will work for app==2 and app==1 because in app==2, we have nD*the number of total clusters as the number of parameters
  // while in app==1 it is nD*number of clusters of that class (only that clusters are also = number of params)

    int thetaid;

    float *belief = new float[prm->nD];
    int pixel=0;

    MRFEnergy* mrfe = (MRFEnergy*)mrf;

  for (int z=0; z < depth; z++)
  {
    for (int y = 0; y < height; y++)
    {
      uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

      for (int x = 0; x < width; x++)
      {
        if (!(prm->unknown == 1 && dig[x]==prm->nD))         // unknown pixels do not count during learning parameters
          continue;

        mrfe->getBelief(pixel,belief);

        for (int d=0; d<prm->nD; d++)
        {
          int dm = d;
          if (fpvp->app==2)
            dm=0;

          if (((int)appdirImage[index][z][dm].Pixel(x, y, 0) >= nT && fpvp->app==1) || ((int)appdirImage[index][z][0].Pixel(x, y, 0) >= fpvp->nA && fpvp->app==2))
          {
            continue;
          }
          getCostAppDiscPatch((int)appdirImage[index][z][dm].Pixel(x, y, 0), appclass[dm]->getPatchSize(), x, y, width, height, fpvp->thetaA, d, nT, thetaid);

          // assume the patch size is the same for all classes
          int patchSize = appclass[dm]->getPatchSize();

          if (x>(int)patchSize/2 && x<width-patchSize/2 && y>(int)patchSize/2 && y<height-patchSize/2)
            totalEmpirDistA[thetaid] ++;
        }

            pixel++;

      }
    }
  }

  delete [] belief;

}




void crf::computeHiddenDistH(int index)
{

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nH == 0)
    return; // nothing to do

  int nHpC = fpvp->nH/prm->nD; // HoG vocabulary size

    int thetaid;

    MRFEnergy* mrfe = (MRFEnergy*)mrf;

    float *belief = new float[prm->nD];
    int pixel=0;

    for (int z = 0; z < depth; z++)
    {
      for (int y = 0; y < height; y++)
      {
        uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
        uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

        for (int x = 0; x < width; x++)
        {
          if (!(prm->unknown == 1 && dig[x]==prm->nD))         // unknown pixels do not count during learning parameters
            continue;

        mrfe->getBelief(pixel,belief);

        for (int d=0; d<prm->nD; d++)
        {
          unsigned short hogval = hogdirImage[index][z].Pixel(x, y, 0);
          getCostHoGDiscBins((int) hogval, d, nHpC, fpvp->thetaH, thetaid);
          //getCostHoGDiscBinsTemp((int) hogval, d, fpvp->thetaH, thetaid);

          // for us the costs are the theta param, but in differentiation, they vanish,
          // so its just whether the feature is enabled or not.
          // move costs pointer to next pixel
          if (thetaid >= 0)
            totalEmpirDistH[thetaid] += exp(belief[d]);

        }

          pixel++;

        }
      }
    }

    delete [] belief;

}


void crf::computeHiddenDistB(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nB == 0)
    return; // nothing to do

  int thetaid;

  MRFEnergy* mrfe = (MRFEnergy*)mrf;

  float *belief = new float[prm->nD];
  int pixel=0;

  for (int z = 0; z < depth; z++)
  {
    for (int y = 0; y < height; y++)
    {
      uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

      for (int x = 0; x < width; x++)
      {
        if (!(prm->unknown == 1 && dig[x]==prm->nD))         // unknown pixels do not count during learning parameters
          continue;

        mrfe->getBelief(pixel,belief);

        for (int d=0; d<prm->nD; d++)
        {
          getCostBias(d, fpvp->thetaB, thetaid);

          // for us the costs are the theta param, but in differentiation, they vanish,
          // so its just whether the feature is enabled or not.
          // move costs pointer to next pixel
          if (thetaid >= 0)
            totalEmpirDistB[thetaid] += exp(belief[d]);

        }

        pixel++;

      }
    }
  }

  delete [] belief;
}

void crf::computeHiddenDistE(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nE == 0)
    return; // nothing to do

  int thetaid;

  MRFEnergy* mrfe = (MRFEnergy*)mrf;

  float *belief = new float[prm->nD];
  int pixel=0;

  for (int z = 0; z < depth; z++)
  {
    for (int y = 0; y < height; y++)
    {
      uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

      for (int x = 0; x < width; x++)
      {
        if (!(prm->unknown == 1 && dig[x]==prm->nD))         // unknown pixels do not count during learning parameters
          continue;

        mrfe->getBelief(pixel,belief);

        for (int d=0; d<prm->nD; d++)
        {
          getCostInverseClassSize(d, inverseClassSize, fpvp->thetaE, thetaid);

          // for us the costs are the theta param, but in differentiation, they vanish,
          // so its just whether the feature is enabled or not.
          // move costs pointer to next pixel
          // assert(0); // do i need to multiply by classinverseSize here ? above comment is false
          if (thetaid >= 0)
            totalEmpirDistE[thetaid] += exp(belief[d]) * inverseClassSize[d];

        }

        pixel++;

      }
    }
  }

  delete [] belief;
}


void crf::computeHiddenDistM(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nM == 0)
    return; // nothing to do

  int nMpC = fpvp->nM/prm->nD; // Mot vocabulary size

    int thetaid;

    MRFEnergy* mrfe = (MRFEnergy*)mrf;

    float *belief = new float[prm->nD];
    int pixel=0;

    for (int z = 0; z < depth; z++)
    {
      for (int y = 0; y < height; y++)
      {
        uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
        uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

        for (int x = 0; x < width; x++)
        {
          if (!(prm->unknown == 1 && dig[x]==prm->nD))         // unknown pixels do not count during learning parameters
            continue;

        mrfe->getBelief(pixel,belief);

        for (int d=0; d<prm->nD; d++)
        {
          unsigned short motval = motdirImage[index][z].Pixel(x, y, 0);
          getCostmotDiscBins((int) motval, d, nMpC, fpvp->thetaM, thetaid);

          totalEmpirDistM[thetaid] += exp(belief[d]);
        }

          pixel++;

        }
      }
    }

    delete [] belief;

}


void crf::computeHiddenDistL(int index)
{

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nL == 0)
    return; // nothing to do

  int nLpC = fpvp->nL/prm->nD;  // this variable will be useful only for cube setting (loc==1) or loc==9

    int thetaid;

    MRFEnergy* mrfe = (MRFEnergy*)mrf;

    int zslice = 0;
  if (timeseq == 0 && (fpvp->loc != 0 || fpvp->klrpv.locKlr ==2))
    zslice = startSliceNo[index];

    float *belief = new float[prm->nD];
    int pixel=0;

    for (int z = 0; z < depth; z++)
    {

      for (int y = 0; y < height; y++)
      {
        uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
        uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

        for (int x = 0; x < width; x++)
        {
          if (!(prm->unknown == 1 && dig[x]==prm->nD))         // unknown pixels do not count during learning parameters
            continue;

        mrfe->getBelief(pixel,belief);

        for (int d=0; d<prm->nD; d++)
        {
          if (fpvp->loc==1)
          {
            /*
            unsigned short locpix = locdirImage[index][di[x]][z].Pixel(x,y,0);
            getCostLocationDiscCube(locpix, x, y, z, fpvp->loccubex, fpvp->loccubey, fpvp->loccubez, numTrainingPats, nLpC, fpvp->thetaL, width, height, depth, di[x], thetaid);
            double locpixd = (double)locpix/ (fpvp->loccubex * fpvp->loccubey * fpvp->loccubez * numTrainingPats);

            totalEmpirDistL[thetaid] += += exp(belief[d])*(1-locpixd);
            */

            std::cout << " Need to fix this - location " << std::endl;
            exit(1);


          } else if (fpvp->loc==6 || fpvp->loc==8) { // CRF hard assignment and soft

            getCostLocationDiscGaussian(LocationMVclass, z, zslice, x, y, fpvp->thetaL, d, fpvp->loc, fpvp->generative, thetaid);
            totalEmpirDistL[thetaid] += exp(belief[d]);

          } else if (fpvp->loc==9 ||  (fpvp->loc==7 && timeseq == 1)) {  //CRF soft or hard assignment redefined

            getCostLocationDiscGaussian2(LocationMVclass, z, zslice, x, y, fpvp->thetaL, d, fpvp->loc, fpvp->generative, prm->nD, thetaid);
            totalEmpirDistL[thetaid] += exp(belief[d]);

          } else if  (fpvp->loc==7 && timeseq == 0) {  //CRF hard assignment redefined

            getCostLocationDiscGaussian3(index, z, x, y, d, fpvp->generative, prm->nD, thetaid);
            totalEmpirDistL[thetaid] += exp(belief[d]);
          }
        }

          pixel++;

        }
      }
    }

    delete [] belief;

}


void crf::computeHiddenDistW(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->klrpv.nW == 0)
    return; // nothing to do

  int dimW = wparam.size2();
  int Ntr = wparam.size1();

  MRFEnergy* mrfe = (MRFEnergy*)mrf;

  int zslice = 0;

  if (timeseq == 0 && (fpvp->loc != 0 || fpvp->klrpv.locKlr ==2))
    zslice = startSliceNo[index];


  uchar *pix1t = 0;
  unsigned short pix1 = 0;

  ublas::vector<DBL_TYPE> Xtest(Xtrain.size2());
  ublas::vector<DBL_TYPE> f;


    float *belief = new float[prm->nD];
    int pixel=0;

  for (int z = 0; z < depth; z++)
  {

    matrix<DBL_TYPE> hogdescs, appdescs;

    if (fpvp->klrpv.klr==1)
      readklrdescfiles(hogdescs, appdescs, index, z+zslice);

    for (int y = 0; y < height; y++)
    {
      uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);

      for (int x = 0; x < width; x++)
      {
        if (!(prm->unknown == 1 && dig[x]==prm->nD))         // unknown pixels do not count during learning parameters
          continue;

        mrfe->getBelief(pixel,belief);


        getDataCostKlr(f, Xtest, index, z, x, y, zslice, hogdescs, appdescs, width, height);

        for (int j=0; j<Ntr; j++)
        {
          double Kval = ivm->classify_klrexp(&Xtest, j);

          for (int d=0; d<prm->nD; d++)
          {
            int thetaid = j*dimW + d;
            totalEmpirDistW[thetaid] += exp(belief[d]) * Kval;
          }
        }

        pixel++;

      }
    }
  }

  delete [] belief;

}



void crf::computeHiddenDistO(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (prm->timeseq!=1)
    return;

  if (fpvp->nO == 0)
    return; // nothing to do

  int nOpC = 0;
  if (fpvp->opflowlocal == 1)
  {
    nOpC = fpvp->nO/(prm->nD * fpvp->opflowlocalframes.size()); // there are (u,v) variables per frame except last frame for continuous features
  // but is more complicated for bin features (Hiddened similar to color uv.
  }

  int thetaid;
  int thetaido[2]={0,0};
  float indcost[2] = {0.0, 0.0};

  MRFEnergy* mrfe = (MRFEnergy*)mrf;

  float *belief = new float[prm->nD];
  int pixel = 0;

  // the last frame should not be counted towards the cost
  for (int z = 0; z < depth-1; z++)
  {
    for (int y = 0; y < height; y++)
    {
      uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
      uchar *dig = &gtdirImage[index][z].Pixel(0, y, 0);
      for (int x = 0; x < width; x++)
      {
        if (!(prm->unknown == 1 && dig[x]==prm->nD))         // unknown pixels do not count during learning parameters
          continue;

        mrfe->getBelief(pixel,belief);

        for (int d=0; d < prm->nD ; d++)
        {
          int idx = frameExistsInFrameList(z);
          if (idx != -1)
          {
            float u = flowvector[index][z][0][y*width+x];
            float v = flowvector[index][z][1][y*width+x];
            if (0) // optical flow continuous
            {
              getIndicatorCostOpticalFlow(depth, z, d, fpvp->thetaO, thetaido, indcost, fpvp->opflowlocalframes.size(), idx);
              totalEmpirDistO[thetaido[0]] += exp(belief[d]) * u;
              totalEmpirDistO[thetaido[1]] += exp(belief[d]) * v;
            }
            else // binned optical flow features
            {
              getCostOpflowBinsuv(u, v, nOpC, d, fpvp->thetaO, thetaid, fpvp->opflowlocalframes.size(), idx);
              totalEmpirDistO[thetaid] += exp(belief[d]);
            }
          }
          // consider the costs for the theta param,
          // since they exist after differentiation
          // no use of the indcost since they are the theta values
        }
        pixel++;
      }
    }

  }

  delete [] belief;
}





// Note unknown region can currently only be used with fixed pairwise parameter
// ie when there is only one -v value.
// So currently there are no changes to computeHiddenDistV


// this will not learn parameters related to smoothed contrast sensitive gradients
void crf::computeHiddenDistV(int index)
{

  if (fpvp->nV == 0 && fpvp->nC == 0 && fpvp->nZ == 0)
    return; // nothing to do

  if (fpvp->csgrad == 0)
  {
    if (globalP.getOpflowConnection() == 1)
    {
      updateQuantizedPairwiseHiddenDistributionOpflow(index);
    } else {
      updateQuantizedPairwiseHiddenDistribution(index);
    }
  }
  else if (fpvp->csgrad == 1)
  {
    if (globalP.getOpflowConnection() == 1)
    {
      updateSmoothPairwiseHiddenDistributionOpflow(index);
    } else {
      updateSmoothPairwiseHiddenDistribution(index);
    }
  }
}


void crf::updateQuantizedPairwiseHiddenDistributionOpflow(int index)
{
  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  assert(fpvp->gradVZ != 1);
  assert(globalP.getOpflowVolumeIndex() == index);

  if (fpvp->gradContext == 0)
  {
    std::cout << " Not implemented gradcontext = 0 \n";
  }

  MRFEnergy* mrfe = (MRFEnergy*)mrf;

  uchar *dig = 0;
  uchar *dig2 = 0;


  float totBelEqDisp, belNotEq;
  float *belEqDisp;
  float *belief = new float[prm->nD*prm->nD];

  std::map <std::pair <int, int>, std::pair <int, MRF::CostVal> >::iterator iter = globalP.edgeGlobalF[index].begin();

  for (; iter != globalP.edgeGlobalF[index].end(); iter++)
  {
    int node1 = iter->first.first;
    int node2 = iter->first.second;

    int optf = iter->second.first;
    int gradval = (int) iter->second.second;

    mrfe->getEdgeBelief( node1, node2, belief);

    int z1 = node1 / (width*height);
    int z2 = node2 / (width*height);

    int rem1 = node1 - z1*(width*height);
    int rem2 = node2 - z2*(width*height);

    int y1 = rem1 / width;
    int y2 = rem2 / width;

    int x1 = rem1 - y1*width;
    int x2 = rem2 - y2*width;

    dig = &gtdirImage[index][z1].Pixel(x1, y1, 0);
    dig2 = &gtdirImage[index][z2].Pixel(x2, y2, 0);

    if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
    {
      for (int d1v=0; d1v<prm->nD; d1v++)
      {
        for (int d2v=0; d2v<prm->nD; d2v++)
        {
          if (fpvp->nV > 1)  // this makes sure it is not smooth constrast sensitive
          {
            int d1 = d1v;
            int d2 = d2v;

            if (d1>d2)
            {
              int swap = d2;
              d2 = d1;
              d1 = swap;
            }

            if (optf == 0) // spatial ie x or y
            {
              totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] +=  exp(belief[d1*prm->nD + d2]);
            } else { // optical flow connection
              if (fpvp->gradVZ == 2)
              {
                totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] +=  exp(belief[d1*prm->nD + d2]);
              }
              if (fpvp->gradVZ == 3)
              {
                totalEmpirDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] +=  exp(belief[d1*prm->nD + d2]);
              }
            }

          }
        }
      }

    }
  }

  delete[] belief;

}







void crf::updateQuantizedPairwiseHiddenDistribution(int index)
{

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  MRFEnergy* mrfe = (MRFEnergy*)mrf;

  int depthrun = 0;

  /*
  if (fpvp->gradVZ == 1)
    depthrun = depth;
  else if (fpvp->gradVZ == 2 || fpvp->gradVZ == 3)
    depthrun = depth-1;
   */

  depthrun = depth;

  uchar *di = 0;
  uchar *di2 = 0;
  uchar *di3 = 0;

  uchar *dig = 0;
  uchar *dig2 = 0;
  uchar *dig3 = 0;

  float totBelEqDisp, belNotEq;
  float *belEqDisp;
  float *belief = new float[prm->nD*prm->nD];


  for (int z=0; z < depthrun; z++)
  {
    for (int y = 0; y < height-1; y++)
    {
      di = &dirDisp[index][z].Pixel(0, y, 0);
      di2 = &dirDisp[index][z].Pixel(0, y+1, 0);

      if ((fpvp->gradVZ == 2 || fpvp->gradVZ == 3)&&(z<depthrun-1))
        di3 = &dirDisp[index][z+1].Pixel(0, y, 0);

      for (int x = 0; x < width-1; x++)
      {
        int gh, gv, gd;
        gh  = dirgrad[index][z].Pixel(x, y, 0);
        gv  = dirgrad[index][z].Pixel(x, y, 1);
        gd  = dirgrad[index][z].Pixel(x, y, 2);


        //CHECK need to make sure it has the gradContext option too
        // i am not fixing this gradOnly option
        if (fpvp->gradContext==0)
        {

          std::cout << " Not supposed to be working really check this !!! \n";
          exit(1);


        } else { //gradContext option


          if (fpvp->nV > 1)  // this makes sure it is not smooth constrast sensitive
          {

            // Horizontal Edge

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z].Pixel(x+1, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
            {
              int p1, p2;

              p1 = z*(width*height) + y*width + x;
              p2 = z*(width*height) + y*width + x+1;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrfe->getEdgeBelief( p1, p2, belief);

              for (int d1v=0; d1v<prm->nD; d1v++)
              {
                for (int d2v=0; d2v<prm->nD; d2v++)
                {
                  int d1 = d1v;
                  int d2 = d2v;

                  if (d1>d2)
                  {
                    int swap = d2;
                    d2 = d1;
                    d1 = swap;
                  }

                  totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gh] +=  exp(belief[d1*prm->nD + d2]);

                }
              }
            }



            // Vertical Edge
            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z].Pixel(x, y+1, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
            {

              int p1 = z*(width*height) + y*width + x;
              int p2 = z*(width*height) + (y+1)*width + x;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrfe->getEdgeBelief( p1, p2, belief);

              for (int d1v=0; d1v<prm->nD; d1v++)
              {
                for (int d2v=0; d2v<prm->nD; d2v++)
                {
                  int d1 = d1v;
                  int d2 = d2v;

                  if (d1>d2)
                  {
                    int swap = d2;
                    d2 = d1;
                    d1 = swap;
                  }
                  totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gv] +=  exp(belief[d1*prm->nD + d2]);

                }
              }
            }
          }


          if (fpvp->gradVZ == 2 && fpvp->nV > 1 && z < depthrun-1)
          {

            // Depth Edge
            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z+1].Pixel(x, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
            {
              int p1 = z*(width*height) + y*width + x;
              int p2 = (z+1)*(width*height) + y*width + x;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrfe->getEdgeBelief( p1, p2, belief);

              for (int d1v=0; d1v<prm->nD; d1v++)
              {
                for (int d2v=0; d2v<prm->nD; d2v++)
                {
                  int d1 = d1v;
                  int d2 = d2v;

                  if (d1>d2)
                  {
                    int swap = d2;
                    d2 = d1;
                    d1 = swap;
                  }

                  totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gd] +=  exp(belief[d1*prm->nD + d2]);

                }
              }
            }

          }


          if (fpvp->gradVZ == 3 && fpvp->nZ > 1 && z < depthrun-1)
          {

            // Depth Edge
            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z+1].Pixel(x, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
            {

              int p1 = z*(width*height) + y*width + x;
              int p2 = (z+1)*(width*height) + y*width + x;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrfe->getEdgeBelief( p1, p2, belief);


              for (int d1v=0; d1v<prm->nD; d1v++)
              {
                for (int d2v=0; d2v<prm->nD; d2v++)
                {
                  int d1 = d1v;
                  int d2 = d2v;

                  if (d1>d2)
                  {
                    int swap = d2;
                    d2 = d1;
                    d1 = swap;
                  }

                  totalEmpirDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nGZ + gd] +=  exp(belief[d1*prm->nD + d2]);

                }
              }
            }
          }
        }
      }
    }
  }

  delete[] belief;

}




void crf::updateSmoothPairwiseHiddenDistributionOpflow(int index)
{
  assert(0); // not coded

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  assert(fpvp->gradVZ != 1);
  assert(globalP.getOpflowVolumeIndex() == index);

  if (fpvp->gradContext == 0)
  {
    std::cout << " Not implemented gradcontext = 0 \n";
  }

  MRFEnergy* mrfe = (MRFEnergy*)mrf;

  uchar *dig = 0;
  uchar *dig2 = 0;


  float totBelEqDisp, belNotEq;
  float *belEqDisp;
  float *belief = new float[prm->nD*prm->nD];

  std::map <std::pair <int, int>, std::pair <int, MRF::CostVal> >::iterator iter = globalP.edgeGlobalF[index].begin();

  for (; iter != globalP.edgeGlobalF[index].end(); iter++)
  {
    int node1 = iter->first.first;
    int node2 = iter->first.second;

    int optf = iter->second.first;
    int gradval = (int) iter->second.second;

    mrfe->getEdgeBelief( node1, node2, belief);

    int z1 = node1 / (width*height);
    int z2 = node2 / (width*height);

    int rem1 = node1 - z1*(width*height);
    int rem2 = node2 - z2*(width*height);

    int y1 = rem1 / width;
    int y2 = rem2 / width;

    int x1 = rem1 - y1*width;
    int x2 = rem2 - y2*width;

    dig = &gtdirImage[index][z1].Pixel(x1, y1, 0);
    dig2 = &gtdirImage[index][z2].Pixel(x2, y2, 0);

    if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
    {

      for (int d1v=0; d1v<prm->nD; d1v++)
      {
        for (int d2v=0; d2v<prm->nD; d2v++)
        {
          if (fpvp->nV > 1)  // this makes sure it is not smooth constrast sensitive
          {
            int d1 = d1v;
            int d2 = d2v;

            if (d1>d2)
            {
              int swap = d2;
              d2 = d1;
              d1 = swap;
            }

            if (optf == 0) // spatial ie x or y
            {
              totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] +=  exp(belief[d1*prm->nD + d2]);
            } else { // optical flow connection
              if (fpvp->gradVZ == 2)
              {
                totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] +=  exp(belief[d1*prm->nD + d2]);
              }
              if (fpvp->gradVZ == 3)
              {
                totalEmpirDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gradval] +=  exp(belief[d1*prm->nD + d2]);
              }
            }

          }
        }
      }

    }
  }

  delete[] belief;

}







void crf::updateSmoothPairwiseHiddenDistribution(int index)
{

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  MRFEnergy* mrfe = (MRFEnergy*)mrf;

  int depthrun = 0;

  /*
  if (fpvp->gradVZ == 1)
    depthrun = depth;
  else if (fpvp->gradVZ == 2 || fpvp->gradVZ == 3)
    depthrun = depth-1;
   */

  depthrun = depth;

  uchar *di = 0;
  uchar *di2 = 0;
  uchar *di3 = 0;

  uchar *dig = 0;
  uchar *dig2 = 0;
  uchar *dig3 = 0;

  float totBelEqDisp, belNotEq;
  float *belEqDisp;
  float *belief = new float[prm->nD*prm->nD];


  for (int z=0; z < depthrun; z++)
  {
    for (int y = 0; y < height-1; y++)
    {
      di = &dirDisp[index][z].Pixel(0, y, 0);
      di2 = &dirDisp[index][z].Pixel(0, y+1, 0);

      if ((fpvp->gradVZ == 2 || fpvp->gradVZ == 3)&&(z<depthrun-1))
        di3 = &dirDisp[index][z+1].Pixel(0, y, 0);

      for (int x = 0; x < width-1; x++)
      {
        double valh, valv, vald;
        CvScalar s;
        s=cvGet2D(dircsgrad[index][z],y,x); // get the (i,j) pixel value
        valh = exp(s.val[0]);
        valv = exp(s.val[1]);
        vald = exp(s.val[2]);


        //CHECK need to make sure it has the gradContext option too
        // i am not fixing this gradOnly option
        if (fpvp->gradContext==0)
        {
          std::cout << " Not supposed to be working really check this !!! \n";
          exit(1);

        } else { //gradContext option

          if (fpvp->nV == 1)  // this makes sure it is not smooth constrast sensitive
          {

            // Horizontal Edge

            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z].Pixel(x+1, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
            {
              int p1, p2;

              p1 = z*(width*height) + y*width + x;
              p2 = z*(width*height) + y*width + x+1;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrfe->getEdgeBelief( p1, p2, belief);

              for (int d1v=0; d1v<prm->nD; d1v++)
              {
                for (int d2v=0; d2v<prm->nD; d2v++)
                {
                  int d1 = d1v;
                  int d2 = d2v;

                  if (prm->unknown == 1 && (d1==prm->nD || d2==prm->nD))
                  {
                    if (d1 == d2)
                    {} // since [yp != yq] in Vpq
                    else
                    {
                      totalEmpirDistV[0] += valh * exp(belief[d1*prm->nD + d2]);
                    }
                  }

                }
              }
            }



            // Vertical Edge
            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z].Pixel(x, y+1, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
            {

              int p1 = z*(width*height) + y*width + x;
              int p2 = z*(width*height) + (y+1)*width + x;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrfe->getEdgeBelief( p1, p2, belief);

              for (int d1v=0; d1v<prm->nD; d1v++)
              {
                for (int d2v=0; d2v<prm->nD; d2v++)
                {
                  int d1 = d1v;
                  int d2 = d2v;

                  if (prm->unknown == 1 && (d1==prm->nD || d2==prm->nD))
                  {
                    if (d1 == d2)
                    {} // since [yp != yq] in Vpq
                    else
                    {
                      totalEmpirDistV[0] += valv * exp(belief[d1*prm->nD + d2]);
                    }
                  }

                }
              }
            }
          }


          if (fpvp->gradVZ == 2 && fpvp->nV == 1 && z < depthrun-1)
          {

            // Depth Edge
            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z+1].Pixel(x, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
            {
              int p1 = z*(width*height) + y*width + x;
              int p2 = (z+1)*(width*height) + y*width + x;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrfe->getEdgeBelief( p1, p2, belief);

              for (int d1v=0; d1v<prm->nD; d1v++)
              {
                for (int d2v=0; d2v<prm->nD; d2v++)
                {
                  int d1 = d1v;
                  int d2 = d2v;

                  if (prm->unknown == 1 && (d1==prm->nD || d2==prm->nD))
                  {
                    if (d1 == d2)
                    {} // since [yp != yq] in Vpq
                    else
                    {
                      totalEmpirDistV[0] += vald * exp(belief[d1*prm->nD + d2]);
                    }
                  }
                }
              }
            }

          }


          if (fpvp->gradVZ == 3 && fpvp->nZ == 1 && z < depthrun-1)
          {

            // Depth Edge
            dig = &gtdirImage[index][z].Pixel(x, y, 0);
            dig2 = &gtdirImage[index][z+1].Pixel(x, y, 0);

            if (prm->unknown == 1 && (dig[0]==prm->nD || dig2[0]==prm->nD))
            {

              int p1 = z*(width*height) + y*width + x;
              int p2 = (z+1)*(width*height) + y*width + x;
              // don't i need to check the border conditions? as well as the original place
              // i copied from: CHECK
              mrfe->getEdgeBelief( p1, p2, belief);


              for (int d1v=0; d1v<prm->nD; d1v++)
              {
                for (int d2v=0; d2v<prm->nD; d2v++)
                {
                  int d1 = d1v;
                  int d2 = d2v;

                  if (prm->unknown == 1 && (d1==prm->nD || d2==prm->nD))
                  {
                    if (d1 == d2)
                    {} // since [yp != yq] in Vpq
                    else
                    {
                      totalEmpirDistZ[0] += vald * exp(belief[d1*prm->nD + d2]);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  delete[] belief;
}
