/*
 * baseModel.cpp
 *
 *  Created on: Mar 15, 2011
 *      Author: bhole
 */

#include "baseModel.h"

extern GlobalPairwiseParam globalP;


void baseModel::allocateDataCostSpace()
{
	// in the base model you only need to bother about the dataCost for one image
	// for model like crf you will need space for the entire volume/patient
	// for logregpd4 you will need two or 3 slices at a time. So override this function

	std::cout << " allocateDataCostSpace() This should be model specific " << std::endl;
	exit(1);
}


void baseModel::deAllocateDataCostSpace()
{
	std::cout << " deAllocateDataCostSpace() This should be model specific " << std::endl;
	exit(1);
}




void baseModel::allocateSmoothCostGlobalSpace()
{
	// in the base model you only need to bother about the dataCost for one image
	// for model like crf you will need space for the entire volume/patient
	// for logregpd4 you will need two or 3 slices at a time. So override this function

	std::cout << " allocateSmoothCostGlobalSpace() This should be model specific " << std::endl;
	exit(1);
}


void baseModel::deAllocateSmoothCostGlobalSpace()
{
	std::cout << " deAllocateSmoothCostGlobalSpace() This should be model specific " << std::endl;
	exit(1);
}




baseModel::baseModel(Parameters* prmout, std::vector<float> theta) :
    features(&prmout->featurepv), IO(prmout->iopv.outstem, prmout->timeseq),
    bfgs(prmout->graddescpv.bfgs_flag, theta, prmout->graddescpv.beta_dash, prmout->graddescpv.beta, prmout->graddescpv.alpha, prmout->graddescpv.m_bfgs, prmout->graddescpv.stepmaxconst),
    evaluation(prmout->nD)
{
	prm = prmout;
	oldreldiffnorm = -1;
	oldloglikelihood = 0.0;
	oldrmsdiffDist = -1;
}

baseModel::~baseModel()
{

	prm = 0;  // we don't free because prm has not been dynamically allocated.

}



void baseModel::initializeVarToZero()
{

	initializeEvalVarToZero();

	// oldloglikelihood = 0.0;
	loglikelihood = 0.0;
	loglikelihoodgen = 0.0;

}


// d2 should be larger than d1
void baseModel::reorder(int *d1, int *d2)
{
	if (*d1> *d2)
	{
		int temp = *d2;
		*d2 = *d1;
		*d1 = temp;
	}

}

int baseModel::frameExistsInFrameList(const int &frame_num)
{
  for (int i=0; i<fpvp->opflowlocalframes.size(); i++)
  {
    if (frame_num == fpvp->opflowlocalframes[i])
      return i;
  }
  return -1;
}

void baseModel::checkSanity()
{
  // check test directory
  int testpresent = 0;
  int trainpresent = 0;

  for (int i=0; i < prm->iopv.testDirIndexV.size(); i++)
  {
    if (prm->iopv.testDirIndexV[i] == 1)
      testpresent = 1;
    if (prm->iopv.testDirIndexV[i] == 0)
      trainpresent = 1;
  }
  if (testpresent == 1 && prm->iopv.testDir == 1)
  {}
  else if (testpresent == 0 && prm->iopv.testDir == 1)
  {
    std::cout << " No test directory present though test flag is set " <<std::endl;
    exit(1);
  }
  else if (testpresent == 1 && prm->iopv.testDir == 0)
  {
    std::cout << " test directory present though test flag is not set " <<std::endl;
    exit(1);
  }


  if (prm->crfpv.crfp == 1 && prm->crfpv.hidden == 1 && prm->unknown == 0)
  {
    std::cout << " unknown has to be 1 for hcrf " << std::endl;
  }

  if (prm->crfpv.crfp == 1 && prm->crfpv.hidden == 1 && prm->interactive == 1)
  {
    std::cout << " interactive=1 for hcrf never go together" << std::endl;
  }

  if (prm->crfpv.crfp == 1)
  {
    if (prm->crfpv.hidden == 1 && prm->unknown == 0 && prm->interactive == 1 && testpresent == 1)
    {
      std::cout << " Set hidden to zero to avoid confusion " <<std::endl;
      exit(1);
    }

    if (prm->interactive == 1 && trainpresent == 1)
    {
      std::cout << " interactive will be ignored for training folders " <<std::endl;
    }

  }



}



double baseModel::getDataCost(int genparam, int width, int height, int depth, unsigned int j, unsigned int i, int y, int x, int d, int zslice)
{

	int nbins = double(fpvp->nU)/prm->nD;
	double binner = 0;
	if (fpvp->nU>0)
		binner = fpvp->rangeI/nbins;

	int nT = fpvp->nA/prm->nD; //this gives centers per class
	// this will work for app==2 and app==1 because in app==2, we have nD*the number of total clusters as the number of parameters
	// while in app==1 it is nD*number of clusters of that class (only that clusters are also = number of params)

	int nHpC = fpvp->nH/prm->nD; // HoG vocabulary size

	int nMpC = fpvp->nM/prm->nD; // MOT vocabulary size

	int nLpC = fpvp->nL/prm->nD;  // this variable will be useful only for cube setting (loc==1) or loc==9

	int nOpC = 0;
	if (fpvp->opflowlocal == 1)
	{
	  nOpC = fpvp->nO/(prm->nD * fpvp->opflowlocalframes.size()); // there are (u,v) variables per frame except last frame for continuous features
	// but is more complicated for bin features (modeled similar to color uv.
	}


	if (prm->timeseq==1 && fpvp->loc==5 && prm->nD > 2)
	{
	  std::cout << " Currently the datacost terms for location support only 2 labels \n";
	  std::cout << " To extend to multiple labels, you need an input location image for each label \n";
	  std::cout << " like that done (previous revisions) for location cube \n";
	  exit(1);
	}

	int thetaid;
  int thetaidrgb[3]={0,0,0};
  int thetaidYuv[2]={0,0};

  // for optical flow
  int thetaido[2] = {0,0};
  float indcost[2] = {0.0, 0.0};

	double dsiValue = 0;

	unsigned short hogval;
	if (fpvp->nH > 0)
		hogval = hogdirImage[j][i].Pixel(x, y, 0);

	unsigned short motval;
	if (fpvp->nM > 0)
		motval = motdirImage[j][i].Pixel(x, y, 0);



	uchar *pix1t = 0;
	unsigned short pix1 = 0;

	if (prm->timeseq==1)
		pix1t = &indirImaget[j][i].Pixel(x, y, 0);
	else if (prm->timeseq==0)
		pix1 = indirImage[j][i].Pixel(x, y, 0);


	if (prm->timeseq==0)
	{
		if (fpvp->nU>0 && fpvp->intensity==1)  // intensity as histogram bins
		{
		  dsiValue += getCostIntensityDiscBins(pix1, nbins, binner, fpvp->rangeI,  d, fpvp->thetaU, thetaid);

		} else if (fpvp->intensity==3 && genparam==1) // generative setting with Gaussians
		{
		  dsiValue += getCostIntensityGenGaussian(intObj[d], pix1);

		} else if (fpvp->intensity == 6 && genparam==1) // using NBayes and bins
		{
		  dsiValue += getCostIntensityGenNB(pix1, fpvp->nbinsintNB, d, intNBclass);
		}
	} else if (prm->timeseq==1)
	{

		if (fpvp->nU>0 && fpvp->intensity==1)  // intensity as histogram bins
		{
		  if (prm->featurepv.onlyUV == 1)
		    dsiValue += getCostIntensityDiscBinsuv(pix1t, nbins, d, fpvp->thetaU, thetaidYuv);
		  else
		    dsiValue += getCostIntensityDiscBinsYuv(pix1t, nbins, d, fpvp->thetaU, thetaidYuv);

		} else if (fpvp->intensity==3 && genparam==1) // generative setting with Gaussians
		{
			//dsiValue += getCostIntensityGenGaussian(intObj[d], pix1);
		  std::cout << " Not yet implemented ! " << std::endl;
		  exit(1);

		} else if (fpvp->intensity == 6 && genparam==1) // using NBayes and bins
		{
			//dsiValue += getCostIntensityGenNB(pix1, nbins, d, intNBclass);
		  std::cout << " Not yet implemented ! " << std::endl;
		  exit(1);
		}

	}


	if (fpvp->nA>0 && fpvp->app==1) // we have appearance features
	{
	  // different features for different classes
	  if ((int)appdirImage[j][i][d].Pixel(x, y, 0) >= nT)
	  {
		  // this is to handle the possibility that some of the pixels in the image borders might have invalid values
	  } else
	  {
		  dsiValue += getCostAppDiscPatch((int)appdirImage[j][i][d].Pixel(x, y, 0), appclass[d]->getPatchSize(), x, y, width, height, fpvp->thetaA, d, nT, thetaid);
	  }
	} else if (fpvp->nA>0 && fpvp->app==2) // we have appearance features
	{
	  // common features for all classes
	  if ((int)appdirImage[j][i][0].Pixel(x, y, 0) >= fpvp->nA)
	  {
		  // this is to handle the possibility that some of the pixels in the image borders might have invalid values
	  } else
	  {
		  dsiValue += getCostAppDiscPatch((int)appdirImage[j][i][0].Pixel(x, y, 0), appclass[0]->getPatchSize(), x, y, width, height, fpvp->thetaA, d, nT, thetaid);
	  }
	}
	else if (fpvp->app==3 && genparam==1) // generative setting
	{
	  dsiValue += getCostAppGenPatch((appdirMVProb[j][i][d])(y,x), AppMVclass[d]->getDimension(), x, y, width, height);
	}



	if (fpvp->nH > 0 && fpvp->hog==1)
	{
	  dsiValue += getCostHoGDiscBins((int) hogval, d, nHpC, fpvp->thetaH, thetaid);
	  //dsiValue += getCostHoGDiscBinsTemp((int) hogval, d, fpvp->thetaH, thetaid);

	} else if (fpvp->hog==3 && genparam==1) // generative setting
	{

	  int startPatch = 8; // index starting from 0
	  // this variable indicates where the Hog descriptors are defined
	  // since i have generated them using the entire patient body, don't need to
	  // worry about the z axis and can use all slices here.
	  //WORK-THIS need to change above stuff for flawless working

	  dsiValue += getCostHoGGenBins((hogdirMVProb[j][i][d])(y,x), x, y, width, height, startPatch);

	}


	if (fpvp->nM > 0 && fpvp->mot==1)
	{
	  dsiValue += getCostmotDiscBins((int) motval, d, nMpC, fpvp->thetaM, thetaid);

	} else if (fpvp->mot==3 && genparam==1) // generative setting
	{

		std::cout << " Not yet implemented !!! ";
		exit(1);

	}

	if (fpvp->nL>0 && fpvp->loc==1)
	{ // CRF cube code

	  std::cout << " out dated \n";
	  exit(1);
	  // this function is completely out-dated
	  // and for older series of data, will need to modify significantly for use

	  //dsiValue += getCostLocationDiscCube((unsigned short) locdirImage[j][d][i].Pixel(x,y,0), x, y, i, fpvp->loccubex, fpvp->loccubey, fpvp->loccubez, numTrainingPats, nLpC, fpvp->thetaL, width, height, depth, d, thetaid);

	} else if (fpvp->nL>0 && (fpvp->loc==6 || fpvp->loc==8) )  //CRF hard or soft assignment
	{

	  dsiValue += getCostLocationDiscGaussian(LocationMVclass, i, zslice, x, y, fpvp->thetaL, d, fpvp->loc, genparam, thetaid);

	} else if (fpvp->nL>0 && (fpvp->loc==9) )  //CRF soft assignment redefined
	{
		dsiValue += getCostLocationDiscGaussian2(LocationMVclass, i, zslice, x, y, fpvp->thetaL, d, fpvp->loc, genparam, prm->nD, thetaid);

	} else if (fpvp->nL>0 && (fpvp->loc==7) )  //CRF hard assignment redefined - using features class function
	{
		dsiValue += getCostLocationDiscGaussian3(j , i, x, y, d, genparam, prm->nD, thetaid);

	} else if (fpvp->loc==3 && genparam==1) // generative setting
	{
		dsiValue += getCostLocationGenGaussian(LocationMVclass[d], i, zslice, x, y);
	} else if (prm->timeseq==1 && fpvp->loc==5 && fpvp->nL>0)
	{
	  dsiValue += getCostLocationDT((unsigned short) locdirImage[j][i].Pixel(x,y,0), fpvp->thetaL, d, thetaid);
	} else if (prm->timeseq==1 && fpvp->loc==4 && fpvp->nL>0 && genparam==1)
  {
    dsiValue += getCostLocationDTgen((unsigned short) locdirImage[j][i].Pixel(x,y,0), fpvp->thetaL, d, thetaid);
  }


	// optical flow features
	if (prm->timeseq==1 && fpvp->nO > 0 && fpvp->opflowlocal == 1 )
	{
	  // the last frame should not be counted towards the cost
	  if (i < depth-1)
	  {
	    int idx = frameExistsInFrameList(i);
	    if (idx != -1)
	    {
        float u = flowvector[j][i][0][y*width+x];
        float v = flowvector[j][i][1][y*width+x];
        if (0) // for continuous flow local optical flow feature
        {
          getIndicatorCostOpticalFlow(depth, i, d, fpvp->thetaO, thetaido, indcost, fpvp->opflowlocalframes.size(), idx);
          dsiValue += indcost[0]*u + indcost[1]*v;
        } else // binned approach for optical flow feature
        {
          dsiValue += getCostOpflowBinsuv(u, v, nOpC, d, fpvp->thetaO, thetaid, fpvp->opflowlocalframes.size(), idx);
        }
	    }
	  }
	}

	if (fpvp->nB > 0 && fpvp->bias==1)
	{
	  dsiValue += getCostBias(d, fpvp->thetaB, thetaid);
	}

	if (fpvp->nE > 0 && fpvp->volume==1)
	{
	  dsiValue += getCostInverseClassSize(d, inverseClassSize, fpvp->thetaE, thetaid);
	}



	if (prm->optClassAccuracy == 2)  // type II weighting scheme
		weightClassValue(dsiValue, d);

	return (dsiValue);

}



void baseModel::preComputeClassRatios()
{
  classSize.clear();
	classSize.resize(prm->nD, 0);

	if (prm->optClassAccuracy == 0)
	{
		for (int d=0; d<prm->nD; d++)
			classSize[d] = 0;

	} else
	{

		int numP = gtdirImage.size();
		double totalPoints = 0.0;

		for(int j = 0; j < numP; ++j)
		{

			if (prm->iopv.testDirIndexV[j] == 0) // training
			{
				int depth = gtdirImage[j].size();
				CShape sh = gtdirImage[j][0].Shape();
				int width = sh.width, height = sh.height;

				for(unsigned int i = 0; i < depth; ++i)
				{
					for (int y = 0; y < height; y++)
					{

						for (int x = 0; x < width; x++)
						{
							uchar *gtpix = &gtdirImage[j][i].Pixel(x, y, 0);

							classSize[gtpix[0]] = classSize[gtpix[0]] + 1;
							totalPoints++;
						}
					}
				}
			}

		}

		for (int d=0; d<prm->nD; d++)
			classSize[d] /= totalPoints;
	}

}



void baseModel::preComputeInverseClassRatios()
{
  inverseClassSize.clear();
  inverseClassSize.resize(prm->nD, 0);
  std::vector<int> perClass(prm->nD, 0);

  if (prm->featurepv.volume == 0)
  {
    for (int d=0; d<prm->nD; d++)
      inverseClassSize[d] = 0;

  } else
  {

    int numP = gtdirImage.size();
    double totalPoints = 0.0;

    for(int j = 0; j < numP; ++j)
    {

      if (prm->iopv.testDirIndexV[j] == 0) // training
      {
        int depth = gtdirImage[j].size();
        CShape sh = gtdirImage[j][0].Shape();
        int width = sh.width, height = sh.height;

        for(unsigned int i = 0; i < depth; ++i)
        {
          for (int y = 0; y < height; y++)
          {

            for (int x = 0; x < width; x++)
            {
              uchar *gtpix = &gtdirImage[j][i].Pixel(x, y, 0);

              perClass[gtpix[0]] = perClass[gtpix[0]] + 1;
              totalPoints++;
            }
          }
        }
      }

    }

    double parNorm = 0.0;
    for (int d=0; d<prm->nD; d++)
    {
      assert(perClass[d] > 0); // we should have atleast one pixel of the organ
      parNorm += 1.0/double(perClass[d]);
    }
    for (int d=0; d<prm->nD; d++)
      inverseClassSize[d] =  1.0 / (perClass[d] * parNorm);
  }

}





void baseModel::weightClassValueTemp(double &dsiValue, int d)
{
	if (d==0)
		dsiValue = dsiValue/263.65;
	if (d==1)
		dsiValue = dsiValue/47.24;
	if (d==2)
		dsiValue = dsiValue/4.8;
	if (d==3)
		dsiValue = dsiValue/4.5;
	if (d==4)
		dsiValue = dsiValue/1.0;
	if (d==5)
		dsiValue = dsiValue/7.9;

	if (d>=6 || d<0)
	{
		std::cout << " Shouldn't be reaching here - this is hack code \n";
		exit(0);
	}

}

void baseModel::weightClassValue(double &dsiValue, int d)
{
	dsiValue = dsiValue/classSize[d];
}




double baseModel::getbadCost(int genparam)
{

	double badcost = 1e20;

	if (fpvp->nU > 0)
	  for(int j=0; j<fpvp->nU; j++)
		  badcost += fpvp->thetaU[j];

	if (fpvp->nA > 0)
	for(int j=0; j<fpvp->nA; j++)
	  badcost += fpvp->thetaA[j];

	if (fpvp->nH > 0)
	for(int j=0; j<fpvp->nH; j++)
	  badcost += fpvp->thetaH[j];

	if (fpvp->nM > 0)
	for(int j=0; j<fpvp->nM; j++)
	  badcost += fpvp->thetaM[j];

	if (fpvp->nL > 0)
	for(int j=0; j<fpvp->nL; j++)
	  badcost += fpvp->thetaL[j];

	if (fpvp->nO > 0)
	  for(int j=0; j<fpvp->nO; j++)
	    badcost += fpvp->thetaO[j];

	 if (fpvp->nB > 0)
	    for(int j=0; j<fpvp->nB; j++)
	      badcost += fpvp->thetaB[j];

   if (fpvp->nE > 0)
      for(int j=0; j<fpvp->nE; j++)
        badcost += fpvp->thetaE[j];


	if (genparam==1)
	{

		int numGens=0;
		if (fpvp->intensity==3 || fpvp->intensity==6)
		  numGens++;
		if (fpvp->loc==3 || fpvp->loc==8 || fpvp->loc==9)
		  numGens++;
		if (fpvp->app==3)
		  numGens++;
		if (fpvp->hog==3)
		  numGens++;
		if (fpvp->mot==3)
		  numGens++;

		// assuming worst cost after taking log is -log(10^-99)
		badcost += numGens*(99*3);  //that 3 accounts for 2.7 value of e (natural log)
	}

	return (badcost);

}


void baseModel::getDataCostKlr(ublas::vector<DBL_TYPE> &f, ublas::vector<DBL_TYPE> &Xtest, unsigned int j, unsigned int i, int x, int y, int zslice, matrix<DBL_TYPE> &hogdescs, matrix<DBL_TYPE> &appdescs, int width, int height)
{

//	int dim = Xtrain.size2();

	uchar *pix1t = 0;
	unsigned short pix1 = 0;

	if (prm->timeseq==1)
		pix1t = &indirImaget[j][i].Pixel(x, y, 0);
	else if (prm->timeseq==0)
		pix1 = indirImage[j][i].Pixel(x, y, 0);

	int idxx = 0;

	if (fpvp->klrpv.intensityKlr==2) // for now only intensity case
	{
	  if (prm->timeseq==0)
		  getCostIntensityKlr(Xtest, (double)pix1, minmaxM, idxx);
	  else if (prm->timeseq==1)
		  getCostIntensityKlrYuv(Xtest, pix1t, minmaxM, idxx);
	}

	if (fpvp->klrpv.locKlr==2)
	{
	  // 25 values. // this should be a parameter  // location features
	  // x y z <prob of belonging to each cluster>*16 <prob of belong to each class>*6
	  getCostLocationKlr(Xtest, i, zslice, x, y, prm->nD, minmaxM, LocationMVclass, idxx, fpvp->skipXYZlocklr);

	}

	if (fpvp->klrpv.appKlr==2)
	{
	  getCostAppKlr(Xtest, x, y, width, minmaxM, idxx, fpvp->klrpv.appdimklr, &appdescs);
	}


	if (fpvp->klrpv.hogKlr==2)
	{
	  getCostHoGKlr(Xtest, x, y, width, minmaxM, idxx, fpvp->klrpv.hogdimklr, &hogdescs);
	}

	f = ivm->classify_klrexp(&Xtest);


}


void baseModel::readklrdescfiles(matrix<DBL_TYPE> &hogdescs, matrix<DBL_TYPE> &appdescs, int index, int num)
{

	if (fpvp->klrpv.hogKlr==2)
	{
	  // read file
	  char hogfile[500];
	  sprintf(hogfile, "%s/hog2ddata_%d.csv", hogklrdescdirListFP.at(index).c_str(), num );

	  // read this file in matrix
	  if(!load_one_data(hogfile, hogdescs))
	  {
		  std::cout << "csv read error: Bad data file filename : "<< hogfile << std::endl;
		  exit(1);
	  }

	}

	// need to load appfile if appklr==2

	if (fpvp->klrpv.appKlr==2)
	{
	  // read file
	  char appfile[500];
	  sprintf(appfile, "%s/app2ddatan_%d.csv", appklrdescdirListFP.at(index).c_str(), num );

	  // read this file in matrix
	  if(!load_one_data(appfile, appdescs))
	  {
		  std::cout << "csv read error: Bad data file filename : "<< appfile << endl;
		  exit(1);
	  }

	}

}


void baseModel::firstInference(int genparam, char* fnpart)
{

	MRF::CostVal badcost;

	badcost = getbadCost(genparam);

	for (unsigned int j = 0; j < gtdirImage.size(); ++j)
	{

		int depth = gtdirImage[j].size();
		CShape sh = gtdirImage[j][0].Shape();
		int width = sh.width, height = sh.height, nB = sh.nBands;

		int zslice = 0;

		if (timeseq == 0 && (fpvp->loc != 0  || fpvp->klrpv.locKlr ==2))
			zslice = startSliceNo[j];

		int dsiIndex = 0;

		// dirDisp[j].resize(gtdirImage[j].size());

		for(unsigned int i = 0; i < gtdirImage[j].size(); ++i)
		{
			// sh.nBands = 1;
			// dirDisp[j][i].ReAllocate(sh);

			matrix<DBL_TYPE> hogdescs, appdescs;

			if (fpvp->klrpv.klr==1)
				readklrdescfiles(hogdescs, appdescs, j, i+zslice);


			for (int y = 0; y < height; y++)
			{
				uchar *WTArow = &dirDisp[j][i].Pixel(0, y, 0);

				for (int x = 0; x < width; x++)
				{

					WTArow[x] = 0;

					ublas::vector<DBL_TYPE> f;
					MRF::CostVal bestval = badcost;
					int* bestd = new int [prm->nD];
					int numbest = 0;

					for (int d = 0; d < prm->nD; d++)
					{
					  bestd[d] = d;
					}

					if (fpvp->klrpv.klr==1)  // this may work only for intensity for timeseq=1 cases
					{
						ublas::vector<DBL_TYPE> Xtest(Xtrain.size2());
						getDataCostKlr(f, Xtest, j, i, x, y, zslice, hogdescs, appdescs, width, height);
					}


					for (int d = 0; d < prm->nD; d++)
					{

						MRF::CostVal dsiValue = 0;

						if (fpvp->klrpv.klr==1)
						{
						  dsiValue += -f(d);
						  dsiValue += fpvp->klrpv.biasklr;
						}

						dsiValue += getDataCost(genparam, width, height, depth, j, i , y, x, d, zslice);


						if (dsiValue < bestval)
						{
						  bestval = dsiValue;
						  bestd[0] = d;
						  numbest = 1;
						} else if (dsiValue == bestval)
						{
						  bestd[numbest] = d;
						  numbest++;
						}

						if (prm->timeseq==0)
						{
						  unsigned short pix1 = indirImage[j][i].Pixel(x, y, 0);
						  if (pix1 >= fpvp->rangeI)
						  {
							  bestd[0] = 0; // assume 0 is background
							  numbest = 1;
						  }
						}
					}

					if (numbest == 0)
					  numbest = prm->nD;

					int curr = rand() % numbest;
					WTArow[x] = bestd[curr];

					delete [] bestd;
				}
			}
		}

		outputImageResults(fnpart, j);

	}
}


void baseModel::outputImageResults(char* fnpart, int volno)
{

    char outnameW[500];
    int cnter = 0;
    string fullPath;
    for(std::vector<string>::iterator b = outdirListFP.begin(); b < outdirListFP.end(); ++b)
    {
    	if (volno != cnter)
    	{
    		cnter++;
    		continue;
    	}

    	fullPath = *b;
		if( fullPath[fullPath.length() -1] != dirs)
			fullPath += dirs;

		for(unsigned int j = 0; j < gtdirImage[cnter].size(); ++j)
		{
			sprintf(outnameW, "%s%s%d", fullPath.c_str(), fnpart, j);

			// since the values are in form 0, 1, 2, .. need to scale to 255/nd-1 just so that output is visible in png viewer
			writeDisparities(dirDisp[cnter][j], prm->outscale8, outnameW);

		}
		cnter++;
	}

}



void baseModel::initializeEmpirDistToZero()
{
  totalEmpirDistU.resize(fpvp->nU);
  for (int k = 0; k < fpvp->nU; k++)
  {
	  totalEmpirDistU[k] = 0;
  }

  totalEmpirDistV.resize(fpvp->nV);
  for (int k = 0; k < fpvp->nV; k++)
  {
	  totalEmpirDistV[k] = 0;
  }

  totalEmpirDistC.resize(fpvp->nC);
  for (int k = 0; k < fpvp->nC; k++)
  {
	  totalEmpirDistC[k] = 0;
  }

  totalEmpirDistZ.resize(fpvp->nZ);
  for (int k = 0; k < fpvp->nZ; k++)
  {
	  totalEmpirDistZ[k] = 0;
  }

  totalEmpirDistA.resize(fpvp->nA);
  for (int k = 0; k < fpvp->nA; k++)
  {
	  totalEmpirDistA[k] = 0;
  }

  totalEmpirDistH.resize(fpvp->nH);
  for (int k = 0; k < fpvp->nH; k++)
  {
	  totalEmpirDistH[k] = 0;
  }

  totalEmpirDistM.resize(fpvp->nM);
  for (int k = 0; k < fpvp->nM; k++)
  {
	  totalEmpirDistM[k] = 0;
  }

  totalEmpirDistL.resize(fpvp->nL);
  for (int k = 0; k < fpvp->nL; k++)
  {
	  totalEmpirDistL[k] = 0;
  }

  totalEmpirDistO.resize(fpvp->nO);
  for (int k = 0; k < fpvp->nO; k++)
  {
    totalEmpirDistO[k] = 0;
  }

  totalEmpirDistB.resize(fpvp->nB);
  for (int k = 0; k < fpvp->nB; k++)
  {
    totalEmpirDistB[k] = 0;
  }

  totalEmpirDistE.resize(fpvp->nE);
  for (int k = 0; k < fpvp->nE; k++)
  {
    totalEmpirDistE[k] = 0;
  }


  totalEmpirDistW.resize(fpvp->klrpv.nW);
  for (int k = 0; k < fpvp->klrpv.nW; k++)
  {
	  totalEmpirDistW[k] = 0;
  }
}





void baseModel::initializeModelDistToZero()
{
  totalModelDistU.resize(fpvp->nU);
  for (int k = 0; k < fpvp->nU; k++)
  {
	  totalModelDistU[k] = 0;
  }

  totalModelDistV.resize(fpvp->nV);
  for (int k = 0; k < fpvp->nV; k++)
  {
	  totalModelDistV[k] = 0;
  }

  totalModelDistC.resize(fpvp->nC);
  for (int k = 0; k < fpvp->nC; k++)
  {
	  totalModelDistC[k] = 0;
  }

  totalModelDistZ.resize(fpvp->nZ);
  for (int k = 0; k < fpvp->nZ; k++)
  {
	  totalModelDistZ[k] = 0;
  }

  totalModelDistA.resize(fpvp->nA);
  for (int k = 0; k < fpvp->nA; k++)
  {
	  totalModelDistA[k] = 0;
  }

  totalModelDistH.resize(fpvp->nH);
  for (int k = 0; k < fpvp->nH; k++)
  {
	  totalModelDistH[k] = 0;
  }

  totalModelDistM.resize(fpvp->nM);
  for (int k = 0; k < fpvp->nM; k++)
  {
	  totalModelDistM[k] = 0;
  }


  totalModelDistL.resize(fpvp->nL);
  for (int k = 0; k < fpvp->nL; k++)
  {
	  totalModelDistL[k] = 0;
  }

  totalModelDistO.resize(fpvp->nO);
  for (int k = 0; k < fpvp->nO; k++)
  {
    totalModelDistO[k] = 0;
  }

  totalModelDistB.resize(fpvp->nB);
  for (int k = 0; k < fpvp->nB; k++)
  {
    totalModelDistB[k] = 0;
  }

  totalModelDistE.resize(fpvp->nE);
  for (int k = 0; k < fpvp->nE; k++)
  {
    totalModelDistE[k] = 0;
  }


  totalModelDistW.resize(fpvp->klrpv.nW);
  for (int k = 0; k < fpvp->klrpv.nW; k++)
  {
	  totalModelDistW[k] = 0;
  }
}




void baseModel::initializeDistToZero()
{
  totalEmpirDistU.resize(fpvp->nU);
  totalModelDistU.resize(fpvp->nU);
  for (int k = 0; k < fpvp->nU; k++)
  {
	  totalEmpirDistU[k] = 0;
	  totalModelDistU[k] = 0;
  }

  totalEmpirDistV.resize(fpvp->nV);
  totalModelDistV.resize(fpvp->nV);
  for (int k = 0; k < fpvp->nV; k++)
  {
	  totalEmpirDistV[k] = 0;
	  totalModelDistV[k] = 0;
  }

  totalEmpirDistC.resize(fpvp->nC);
  totalModelDistC.resize(fpvp->nC);
  for (int k = 0; k < fpvp->nC; k++)
  {
	  totalEmpirDistC[k] = 0;
	  totalModelDistC[k] = 0;
  }

  totalEmpirDistZ.resize(fpvp->nZ);
  totalModelDistZ.resize(fpvp->nZ);
  for (int k = 0; k < fpvp->nZ; k++)
  {
	  totalEmpirDistZ[k] = 0;
	  totalModelDistZ[k] = 0;
  }

  totalEmpirDistA.resize(fpvp->nA);
  totalModelDistA.resize(fpvp->nA);
  for (int k = 0; k < fpvp->nA; k++)
  {
	  totalEmpirDistA[k] = 0;
	  totalModelDistA[k] = 0;
  }

  totalEmpirDistH.resize(fpvp->nH);
  totalModelDistH.resize(fpvp->nH);
  for (int k = 0; k < fpvp->nH; k++)
  {
	  totalEmpirDistH[k] = 0;
	  totalModelDistH[k] = 0;
  }

  totalEmpirDistM.resize(fpvp->nM);
  totalModelDistM.resize(fpvp->nM);
  for (int k = 0; k < fpvp->nM; k++)
  {
	  totalEmpirDistM[k] = 0;
	  totalModelDistM[k] = 0;
  }


  totalEmpirDistL.resize(fpvp->nL);
  totalModelDistL.resize(fpvp->nL);
  for (int k = 0; k < fpvp->nL; k++)
  {
	  totalEmpirDistL[k] = 0;
	  totalModelDistL[k] = 0;
  }

  totalEmpirDistO.resize(fpvp->nO);
  totalModelDistO.resize(fpvp->nO);
  for (int k = 0; k < fpvp->nO; k++)
  {
    totalEmpirDistO[k] = 0;
    totalModelDistO[k] = 0;
  }

  totalEmpirDistB.resize(fpvp->nB);
  totalModelDistB.resize(fpvp->nB);
  for (int k = 0; k < fpvp->nB; k++)
  {
    totalEmpirDistB[k] = 0;
    totalModelDistB[k] = 0;
  }

  totalEmpirDistE.resize(fpvp->nB);
  totalModelDistE.resize(fpvp->nB);
  for (int k = 0; k < fpvp->nB; k++)
  {
    totalEmpirDistE[k] = 0;
    totalModelDistE[k] = 0;
  }


  totalEmpirDistW.resize(fpvp->klrpv.nW);
  totalModelDistW.resize(fpvp->klrpv.nW);
  for (int k = 0; k < fpvp->klrpv.nW; k++)
  {
	  totalEmpirDistW[k] = 0;
	  totalModelDistW[k] = 0;
  }

}






void baseModel::modelProcess(unsigned int index)
{
	std::cout << " baseModel modelProcess Shouldn't Reach here " <<std::endl;
	exit(1);
}

void baseModel::hiddenModelProcess(unsigned int index, int iter, int iterout)
{
  std::cout << " baseModel hiddenModelProcess Shouldn't Reach here " <<std::endl;
  exit(1);
}

void baseModel::evaluateResults(unsigned int index)
{
	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;


	for(int z = 0; z < depth; ++z)
	{
		for (int y = 0; y < height; y++)
		{
			uchar *di = &dirDisp[index][z].Pixel(0, y, 0);
			uchar *errmapp = &errormap[index][z].Pixel(0, y, 0);
			uchar *tdi = &gtdirImage[index][z].Pixel(0, y, 0);

			for (int x = 0; x < width; x++)
			{

			  if (prm->unknown == 1)
			  {

			    if (tdi[x] == prm->nD) // unknown region has value nD
			    {
			      continue;
			    }
			  }

				confMatrixPV[tdi[x]][di[x]] = confMatrixPV[tdi[x]][di[x]] + 1;
				totNDPV[tdi[x]] = totNDPV[tdi[x]] + 1;
				totresNDPV[di[x]] = totresNDPV[di[x]] + 1;

				if (di[x]==tdi[x])
				{
					pixAccPV++;
					errmapp[x] = 0;
				}
				else
				{
					cnterrPV++;
					errmapp[x] = 1;
				}

				if (prm->iopv.testDirIndexV[index] == 0) // training
				{
					confMatrixFull[tdi[x]][di[x]] = confMatrixFull[tdi[x]][di[x]] + 1;
					totNDFull[tdi[x]] = totNDFull[tdi[x]] + 1;
					totresNDFull[di[x]] = totresNDFull[di[x]] + 1;

					if (di[x]==tdi[x])
						pixAccTr++;
					else
						cnterr++;

				} else // Test
				{
					confMatrixFullTest[tdi[x]][di[x]] = confMatrixFullTest[tdi[x]][di[x]] + 1;
					totNDFullTest[tdi[x]] = totNDFullTest[tdi[x]] + 1;
					totresNDFullTest[di[x]] = totresNDFullTest[di[x]] + 1;

					if (di[x]==tdi[x])
						pixAccTe++;
					else
						cnterrtest++;
				}
			}
		}
	}

	totPV += double(pixAccPV)+double(cnterrPV++);

	if (prm->iopv.testDirIndexV[index] == 0) // training
	{
	  totfull = totfull + double(pixAccPV)+double(cnterrPV++);
	} else
	{
	  totfulltest = totfulltest + double(pixAccPV)+double(cnterrPV++);
	}


	double dcs_X, dcs_Y, dcs_XinterY;
	// dcs_X registers segments from algorithm
	// dcs_Y registers groundtruth segmentation

	for (int i=0; i<prm->nD; i++)
	{

		  dcs_X = 0;
		  dcs_Y = 0;
		  dcs_XinterY = 0;

		  dcs_XinterY = confMatrixPV[i][i];

		  for (int j=0; j<prm->nD; j++)
			  dcs_X += confMatrixPV[i][j];

		  for (int j=0; j<prm->nD; j++)
			  dcs_Y += confMatrixPV[j][i];

		  if (dcs_X + dcs_Y == 0)
			  dcs[i] = 1.0;
		  else
			  dcs[i] = 2*dcs_XinterY / (dcs_X + dcs_Y);

	}


	for (int i=0; i<prm->nD; i++)
	{
		for (int j=0; j<prm->nD; j++)
		{
			if (totNDPV[i] == 0)
				confMatrixPV[i][j] = 0.0;
			else
				confMatrixPV[i][j] =  confMatrixPV[i][j]/totNDPV[i];
		}
	}

	pixAccPV = pixAccPV / totPV;
	cnterrPV = cnterrPV / totPV;

	for (int i = 0; i<prm->nD; i++)
	{
		  if (totNDPV[i] == 0)
			  avgCaccPV += 1.0;
		  else
			  avgCaccPV += (double)confMatrixPV[i][i];
	}

	avgCaccPV /= prm->nD;

}





// for all volumes (training and testing separate of course)
void baseModel::evaluateResults()
{

	double dcs_X, dcs_Y, dcs_XinterY;
	// dcs_X registers segments from algorithm
	// dcs_Y registers groundtruth segmentation

	// if there are no training examples, then don't count (will result in divide by zero)

	if (numTrainingPats != 0)
	{
		for (int i=0; i<prm->nD; i++)
		{

			  dcs_X = 0;
			  dcs_Y = 0;
			  dcs_XinterY = 0;

			  dcs_XinterY = confMatrixFull[i][i];

			  for (int j=0; j<prm->nD; j++)
				  dcs_X += confMatrixFull[i][j];

			  for (int j=0; j<prm->nD; j++)
				  dcs_Y += confMatrixFull[j][i];

			  if (dcs_X + dcs_Y == 0)
				  dcsTr[i] = 1.0;
			  else
				  dcsTr[i] = 2*dcs_XinterY / (dcs_X + dcs_Y);

		}
	}



	if (numTestPats != 0)
	{
		for (int i=0; i<prm->nD; i++)
		{

			  dcs_X = 0;
			  dcs_Y = 0;
			  dcs_XinterY = 0;

			  dcs_XinterY = confMatrixFullTest[i][i];

			  for (int j=0; j<prm->nD; j++)
				  dcs_X += confMatrixFullTest[i][j];

			  for (int j=0; j<prm->nD; j++)
				  dcs_Y += confMatrixFullTest[j][i];

			  if (dcs_X + dcs_Y == 0)
				  dcsTe[i] = 1.0;
			  else
				  dcsTe[i] = 2*dcs_XinterY / (dcs_X + dcs_Y);

		}
	}


	for (int i=0; i<prm->nD; i++)
	{
		for (int j=0; j<prm->nD; j++)
		{
			if (totNDFull[i] == 0)
				confMatrixFull[i][j] = 0.0;
			else
				confMatrixFull[i][j] =  confMatrixFull[i][j]/totNDFull[i];
		}
	}

	for (int i=0; i<prm->nD; i++)
	{
		for (int j=0; j<prm->nD; j++)
		{
			if (totNDFullTest[i] == 0)
				confMatrixFullTest[i][j] = 0.0;
			else
				confMatrixFullTest[i][j] =  confMatrixFullTest[i][j]/totNDFullTest[i];
		}
	}


	if (numTrainingPats != 0)
	{
		pixAccTr = pixAccTr / totfull;
		cnterr = cnterr / totfull;
	}

	if (numTestPats != 0)
	{
		pixAccTe = pixAccTe / totfulltest;
		cnterrtest = cnterrtest / totfulltest;
	}


	for (int i = 0; i<prm->nD; i++)
	{
		  if (totNDFull[i] == 0)
			  avgCaccTr += 1.0;
		  else
			  avgCaccTr += (double)confMatrixFull[i][i];
	}

	avgCaccTr /= prm->nD;


	for (int i = 0; i<prm->nD; i++)
	{
		  if (totNDFullTest[i] == 0)
			  avgCaccTe += 1.0;
		  else
			  avgCaccTe += (double)confMatrixFullTest[i][i];
	}

	avgCaccTe /= prm->nD;



}



void baseModel::outputConfMatrix(double **confMatrix, int debug)
{
  FILE* filep;
  if (debug == 1)
    filep = debugfile;
  else
    filep = resultoutfile;

	string dumpstr="\t";

	for (int i = 0; i<prm->nD; i++) {
		  dumpstr += "\t";

		  std::stringstream out;
		  out << i;
		  dumpstr += out.str();
	}
	dumpstr += "\n";
	LOG(filep, dumpstr.c_str());

	for (int i = 0; i<prm->nD; i++) {
		  dumpstr = "\t";
		  std::stringstream out;
		  out << i;
		  dumpstr += out.str();
		  dumpstr += "\t";
		  for (int j=0; j<prm->nD; j++) {
			  std::stringstream out2;
			  out2 << fixed << setprecision(3) << (double)confMatrix[i][j];
			  dumpstr += out2.str();
			  dumpstr += "\t";
		  }
		  dumpstr += "\n";
		  LOG(filep, dumpstr.c_str());
	}

}


void baseModel::outputDiceCoeff(std::vector<double> dcs, int debug)
{

  FILE* filep;
  if (debug == 1)
    filep = debugfile;
  else
    filep = resultoutfile;


	string dumpstr="\t";

	for (int i = 0; i<prm->nD; i++) {
		  std::stringstream out2;
		  out2 << fixed << setprecision(3) << (double)dcs[i] << "\t";
		  dumpstr += out2.str();
	}
	dumpstr += "\n\n";
	LOG(filep, dumpstr.c_str());

}



void baseModel::outputPixelCounts(double* totalCounts, int debug)
{

  FILE* filep;
  if (debug == 1)
    filep = debugfile;
  else
    filep = resultoutfile;


  string dumpstr="\t";

  for (int i = 0; i<prm->nD; i++) {
      std::stringstream out2;
      out2 << fixed << setprecision(1) << (double)totalCounts[i] << "\t";
      dumpstr += out2.str();
  }
  dumpstr += "\n\n";
  LOG(filep, dumpstr.c_str());

}






void baseModel::outputEvalutions(unsigned int index)
{

	if (prm->verbose)
	{
		// LOGBUG(prm->verbose, debugfile, "\n Training Case errpercent= %g \n", errpercent);

		if (prm->iopv.testDirIndexV[index]==0)
			LOGBUG(prm->verbose, debugfile, "\n \t \t Training Confusion Matrix \n");
		if (prm->iopv.testDirIndexV[index]==1)
			LOGBUG(prm->verbose, debugfile, "\n \t \t Test Confusion Matrix \n");
		outputConfMatrix(confMatrixPV, 1);

		LOGBUG(prm->verbose, debugfile, "\n \t \t Dice Coefficients \n");
		outputDiceCoeff(dcs, 1);

    LOGBUG(prm->verbose, debugfile, "\n \t \t Ground Truth Pixel Counts \n");
    outputPixelCounts(totNDPV, 1);

    LOGBUG(prm->verbose, debugfile, "\n \t \t Result Pixel Counts \n");
    outputPixelCounts(totresNDPV, 1);

		LOGBUG(prm->verbose, debugfile, "\n Pixel Accuracy : %g \n", pixAccPV);
		LOGBUG(prm->verbose, debugfile, "\n Average Class Accuracy : %g \n", avgCaccPV);
	}

}



void baseModel::outputEvalutions()
{

	LOG(debugfile, "\n  Log Likelihood : %f \n", loglikelihood);
	LOG(debugfile, "  Log Likelihood Gen : %f \n", loglikelihoodgen);

	// cannot print this here since it has not been evaluated yet!!! and variables are temporay
	// printfVec(empirDist, "empirDist (GT)");
	// printfVec(modelDist, "modelDist     ");
	// printfVec(diffDist,  "diffDist      ");




	LOG(debugfile, "\n \t \t Total Training Confusion Matrix \n");
	outputConfMatrix(confMatrixFull, 1);
	LOG(debugfile, "\n \t \t Total Training Dice Coefficients \n");
	outputDiceCoeff(dcsTr, 1);

  LOGBUG(prm->verbose, debugfile, "\n \t \t Ground Truth Pixel Counts \n");
  outputPixelCounts(totNDFull, 1);
  LOGBUG(prm->verbose, debugfile, "\n \t \t Result Pixel Counts \n");
  outputPixelCounts(totresNDFull, 1);

	LOG(debugfile, "\n Total Training Pixel Accuracy : %g \n", pixAccTr);
	LOG(debugfile, "\n Total Training Average Class Accuracy : %g \n", avgCaccTr);


	LOG(debugfile, "\n \t \t Total Test Confusion Matrix \n");
	outputConfMatrix(confMatrixFullTest, 1);
	LOG(debugfile, "\n \t \t Total Test Dice Coefficients \n");
	outputDiceCoeff(dcsTe, 1);

  LOGBUG(prm->verbose, debugfile, "\n \t \t Ground Truth Pixel Counts \n");
  outputPixelCounts(totNDFullTest, 1);
  LOGBUG(prm->verbose, debugfile, "\n \t \t Result Pixel Counts \n");
  outputPixelCounts(totresNDFullTest, 1);

	LOG(debugfile, "\n Total Testing Pixel Accuracy : %g \n", pixAccTe);
	LOG(debugfile, "\n Total Testing Average Class Accuracy : %g \n", avgCaccTe);


	// if (prm.featurepv.generative==1) {
	//	DEBUG_OUT1(verbose, debugfile, "   Total Hybrid Test Error Percent = %g \n", (100*cnterrtestgen)/totfulltestgen);
	//	display_CM_avg(confMatrixFullTestgen, prm.nD, totNDFullTestgen, 0);
	// }

	// DEBUG_OUT0(verbose, debugfile, "Theta Parameters \n");
	// dumpParameters(prm);

}

void baseModel::outputFinalEvalutions()
{

  LOG(resultoutfile, "\n  Log Likelihood : %f \n", loglikelihood);
  LOG(resultoutfile, "  Log Likelihood Gen : %f \n", loglikelihoodgen);

  LOG(resultoutfile, "\n \t \t Total Training Confusion Matrix \n");
  outputConfMatrix(confMatrixFull, 0);
  LOG(resultoutfile, "\n \t \t Total Training Dice Coefficients \n");
  outputDiceCoeff(dcsTr, 0);

  LOGBUG(prm->verbose, debugfile, "\n \t \t Ground Truth Pixel Counts \n");
  outputPixelCounts(totNDFull, 0);
  LOGBUG(prm->verbose, debugfile, "\n \t \t Result Pixel Counts \n");
  outputPixelCounts(totresNDFull, 0);

  LOG(resultoutfile, "\n Total Training Pixel Accuracy : %g \n", pixAccTr);
  LOG(resultoutfile, "\n Total Training Average Class Accuracy : %g \n", avgCaccTr);


  LOG(resultoutfile, "\n \t \t Total Test Confusion Matrix \n");
  outputConfMatrix(confMatrixFullTest, 0);
  LOG(resultoutfile, "\n \t \t Total Test Dice Coefficients \n");
  outputDiceCoeff(dcsTe, 0);

  LOGBUG(prm->verbose, debugfile, "\n \t \t Ground Truth Pixel Counts \n");
  outputPixelCounts(totNDFullTest, 0);
  LOGBUG(prm->verbose, debugfile, "\n \t \t Result Pixel Counts \n");
  outputPixelCounts(totresNDFullTest, 0);

  LOG(resultoutfile, "\n Total Testing Pixel Accuracy : %g \n", pixAccTe);
  LOG(resultoutfile, "\n Total Testing Average Class Accuracy : %g \n", avgCaccTe);

}


void baseModel::outputImageResults(unsigned int index, int iter)
{
  outputImageResults(index, iter, "");
}



void baseModel::outputImageResults(unsigned int index, int iter, const string prefix)
{
    // output inference (argmax label) results
    char outnameW[500];
    string fullPath = outdirListFP[index];
    if( fullPath[fullPath.length() -1] != dirs)
      fullPath += dirs;

    int dumpimages = 0;

    if (prm->graddescpv.maxiter < 30 || iter % 10 == 0 || iter == 1 || iter ==1001)
    {
      dumpimages = 1;
    }

    if (prm->graddescpv.bfgs_flag == 1)
    {
      dumpimages = 1;
    }

    if (dumpimages == 1)
    {
    for(unsigned int j = 0; j < gtdirImage[index].size(); ++j)
    {
        sprintf(outnameW, "%s%sout-%d-iter-%05d", fullPath.c_str(), prefix.c_str(), j, iter);

        // since the values are in form 0, 1, 2, .. need to scale to 255/nd-1 just so that output is visible in png viewer
        writeDisparities(dirDisp[index][j], prm->outscale8, outnameW);

        sprintf(outnameW, "%s%serr-%d-iter-%05d", fullPath.c_str(), prefix.c_str(), j, iter);
        writeDisparities(errormap[index][j], 255, outnameW); // err are 0 or 1 so scaled it to 255

        // generative just puts info in dirDisp and errormap
/*
        if (fpvp->generative==1)
        {
          sprintf(outnameW, "%soutGen-%d-iter-%05d", fullPath.c_str(), j, iter);

          // since the values are in form 0, 1, 2, .. need to scale to 255/nd-1 just so that output is visible in png viewer
          writeDisparities(dirDispGen[index][j], prm->outscale8, outnameW);

          sprintf(outnameW, "%serrGen-%d-iter-%05d", fullPath.c_str(), j, iter);
          writeDisparities(errormapGen[index][j], 255, outnameW); // err are 0 or 1 so scaled it to 255
        }
*/
    }

    }
}


void baseModel::updateDist(int index, int fiterflag)
{
	std::cout << " baseModel Shouldn't Reach here " <<std::endl;
}


void baseModel::computeEmpirDist(int index)
{
	std::cout << " baseModel Shouldn't Reach here " <<std::endl;
}

void baseModel::computeModelDist(int index)
{
	std::cout << " baseModel Shouldn't Reach here " <<std::endl;
}

void baseModel::computeModelDist(int index, int z)
{
	std::cout << " baseModel Shouldn't Reach here " <<std::endl;
}

/*
void baseModel::computeEmpirDistU()
{
	std::cout << " baseModel Shouldn't Reach here " <<std::endl;
}
void baseModel::computeModelDistU()
{
	std::cout << " baseModel Shouldn't Reach here " <<std::endl;
}

void baseModel::computeEmpirDistA()
{
	std::cout << " baseModel Shouldn't Reach here " <<std::endl;
}
void baseModel::computeModelDistA()
{
	std::cout << " baseModel Shouldn't Reach here " <<std::endl;
}

void baseModel::computeEmpirDistH()
{
	std::cout << " baseModel Shouldn't Reach here " <<std::endl;
}
void baseModel::computeModelDistH()
{
	std::cout << " baseModel Shouldn't Reach here " <<std::endl;
}

void baseModel::computeEmpirDistM()
{
	std::cout << " baseModel Shouldn't Reach here " <<std::endl;
}
void baseModel::computeModelDistM()
{
	std::cout << " baseModel Shouldn't Reach here " <<std::endl;
}

void baseModel::computeEmpirDistL()
{
	std::cout << " baseModel Shouldn't Reach here " <<std::endl;
}
void baseModel::computeModelDistL()
{
	std::cout << " baseModel Shouldn't Reach here " <<std::endl;
}


void baseModel::computeEmpirDistV()
{
	std::cout << " baseModel Shouldn't Reach here " <<std::endl;
}
void baseModel::computeModelDistV()
{
	std::cout << " baseModel Shouldn't Reach here " <<std::endl;
}


void baseModel::computeEmpirDistW()
{
	std::cout << " baseModel Shouldn't Reach here " <<std::endl;
}
void baseModel::computeModelDistW()
{
	std::cout << " baseModel Shouldn't Reach here " <<std::endl;
}


*/


void baseModel::dumpThetaParameters(int verbose, int debug)
{
  FILE* filep = 0;
  if (debug == 1)
    filep = debugfile;
  else
    filep = paramoutfile;


	if (fpvp->nU>0 && fpvp->intensity==1)
	{
		for (int ii=0; ii<fpvp->nU; ii++)
		{
			LOGBUG(verbose, filep, "-u %g ", fpvp->thetaU[ii]);
			if ((ii+1)%20 == 0)
				LOGBUG(verbose, filep, "\n");
		}
		LOGBUG(verbose, filep, "\n");

	}


  if (fpvp->nH>0 && fpvp->hog==1)
  {
    for (int ii=0; ii<fpvp->nH; ii++)
    {
    	LOGBUG(verbose, filep, "-h %g ", fpvp->thetaH[ii]);
    	if ((ii+1)%20 == 0)
    	  LOGBUG(verbose, filep, "\n");
    }
    LOGBUG(verbose, filep, "\n");
  }

  if (fpvp->nM>0 && fpvp->mot==1)
  {
    for (int ii=0; ii<fpvp->nM; ii++)
    {
    	LOGBUG(verbose, filep, "-q %g ", fpvp->thetaM[ii]);
    	if ((ii+1)%20 == 0)
    	  LOGBUG(verbose, filep, "\n");
    }
    LOGBUG(verbose, filep, "\n");
  }


  if (fpvp->nA>0 && fpvp->app==1)
  {
    for (int ii=0; ii<fpvp->nA; ii++)
    {
    	LOGBUG(verbose, filep, "-s %g ", fpvp->thetaA[ii]);
    	if ((ii+1)%20 == 0)
    	  LOGBUG(verbose, filep, "\n");
    }
    LOGBUG(verbose, filep, "\n");
  }

  if (fpvp->nL>0 && (fpvp->loc==5 || fpvp->loc==6 || fpvp->loc==7))
  {
    for (int ii=0; ii<fpvp->nL; ii++)
    {
    	LOGBUG(verbose, filep, "-l %g ", fpvp->thetaL[ii]);
    	if ((ii+1)%20 == 0)
    	  LOGBUG(verbose, filep, "\n");
    }
    LOGBUG(verbose, filep, "\n");
  }

  if (fpvp->nO>0 && fpvp->opflowlocal==1)
  {
    for (int ii=0; ii<fpvp->nO; ii++)
    {
      LOGBUG(verbose, filep, "-f %g ", fpvp->thetaO[ii]);
      if ((ii+1)%20 == 0)
        LOGBUG(verbose, filep, "\n");
    }
    LOGBUG(verbose, filep, "\n");
  }


  if (fpvp->nC>0 && fpvp->context==1)
  {
    for (int ii=0; ii<fpvp->nC; ii++)
    {
    	LOGBUG(verbose, filep, "-c %g ", fpvp->thetaC[ii]);
    	if ((ii+1)%20 == 0)
    	  LOGBUG(verbose, filep, "\n");
    }
    LOGBUG(verbose, filep, "\n");
  }


  if (fpvp->nV>0)
  {
    for (int ii=0; ii<fpvp->nV; ii++)
    {
    	LOGBUG(verbose, filep, "-v %g ", fpvp->thetaV[ii]);
    	if ((ii+1)%20 == 0)
    	  LOGBUG(verbose, filep, "\n");
    }
    LOGBUG(verbose, filep, "\n");
  }

  if (fpvp->nZ>0)
  {
    for (int ii=0; ii<fpvp->nZ; ii++)
    {
    	LOGBUG(verbose, filep, "-d %g ", fpvp->thetaZ[ii]);
    	if ((ii+1)%20 == 0)
    	  LOGBUG(verbose, filep, "\n");
      }
    LOGBUG(verbose, filep, "\n");
  }

  if (fpvp->nB>0 && fpvp->bias==1)
  {
    for (int ii=0; ii<fpvp->nB; ii++)
    {
      LOGBUG(verbose, filep, "-j %g ", fpvp->thetaB[ii]);
      if ((ii+1)%20 == 0)
        LOGBUG(verbose, filep, "\n");
    }
    LOGBUG(verbose, filep, "\n");
  }

  if (fpvp->nE>0 && fpvp->volume==1)
  {
    for (int ii=0; ii<fpvp->nE; ii++)
    {
      LOGBUG(verbose, filep, "-p %g ", fpvp->thetaE[ii]);
      if ((ii+1)%20 == 0)
        LOGBUG(verbose, filep, "\n");
    }
    LOGBUG(verbose, filep, "\n");
  }


}




int baseModel::returnConcatThetaWWVector(std::vector<float> *theta)
{
	// follow order thetaU, thetaV, thetaZ, thetaC, thetaA, thetaH, thetaL, thetaM, thetaO, thetaB, thetaE, wparams

	int total = fpvp->nU + fpvp->nV + fpvp->nZ + fpvp->nC + fpvp->nA + fpvp->nH + fpvp->nL + fpvp->nM + fpvp->nO + fpvp->nB + fpvp->nE + fpvp->klrpv.nW ;
	(*theta).resize(total);

	int k, i = 0;

	for (k = 0; k < fpvp->nU; k++)
		(*theta)[i++] = fpvp->thetaU[k];
	for (k = 0; k < fpvp->nV; k++)
		(*theta)[i++] = fpvp->thetaV[k];
	for (k = 0; k < fpvp->nZ; k++)
		(*theta)[i++] = fpvp->thetaZ[k];
	for (k = 0; k < fpvp->nC; k++)
		(*theta)[i++] = fpvp->thetaC[k];
	for (k = 0; k < fpvp->nA; k++)
		(*theta)[i++] = fpvp->thetaA[k];
	for (k = 0; k < fpvp->nH; k++)
		(*theta)[i++] = fpvp->thetaH[k];
	for (k = 0; k < fpvp->nL; k++)
		(*theta)[i++] = fpvp->thetaL[k];
	for (k = 0; k < fpvp->nM; k++)
		(*theta)[i++] = fpvp->thetaM[k];
	for (k = 0; k < fpvp->nO; k++)
	  (*theta)[i++] = fpvp->thetaO[k];
	for (k = 0; k < fpvp->nB; k++)
	  (*theta)[i++] = fpvp->thetaB[k];
	for (k = 0; k < fpvp->nE; k++)
	  (*theta)[i++] = fpvp->thetaE[k];

	fvec thetaW;
	matrixToVectorCopy(wparam, thetaW);
	for (k = 0; k < fpvp->klrpv.nW; k++)
		(*theta)[i++] = thetaW[k];


}




int baseModel::splitThetaWWVectorToOriginal(std::vector<float> theta)
{
	// follow order thetaU, thetaV, thetaZ, thetaC, thetaA, thetaH, thetaL, thetaM, thetaB, thetaE, wparams

	int total = fpvp->nU + fpvp->nV + fpvp->nZ + fpvp->nC + fpvp->nA + fpvp->nH + fpvp->nL + fpvp->nM + fpvp->nO + fpvp->nB + fpvp->nE;

	int k, i = 0;

	for (k = 0; k < fpvp->nU; k++)
		fpvp->thetaU[k] = theta[i++];
	for (k = 0; k < fpvp->nV; k++)
		fpvp->thetaV[k] = theta[i++];
	for (k = 0; k < fpvp->nZ; k++)
		fpvp->thetaZ[k] = theta[i++];
	for (k = 0; k < fpvp->nC; k++)
		fpvp->thetaC[k] = theta[i++];
	for (k = 0; k < fpvp->nA; k++)
		fpvp->thetaA[k] = theta[i++];
	for (k = 0; k < fpvp->nH; k++)
		fpvp->thetaH[k] = theta[i++];
	for (k = 0; k < fpvp->nL; k++)
		fpvp->thetaL[k] = theta[i++];
	for (k = 0; k < fpvp->nM; k++)
		fpvp->thetaM[k] = theta[i++];
	 for (k = 0; k < fpvp->nO; k++)
	    fpvp->thetaO[k] = theta[i++];
   for (k = 0; k < fpvp->nB; k++)
      fpvp->thetaB[k] = theta[i++];
   for (k = 0; k < fpvp->nE; k++)
      fpvp->thetaE[k] = theta[i++];

	fvec thetaW;
	thetaW.resize(fpvp->klrpv.nW);
	for (k = 0; k < fpvp->klrpv.nW; k++)
		thetaW[k] = theta[i++];

	overwriteVectorToMatrix(thetaW, wparam); // ie dim of matrix are known


}


















int baseModel::returnConcatEmpirDistVector(std::vector<float> *empirDist)
{
	// follow order thetaU, thetaV, thetaZ, thetaC, thetaA, thetaH, thetaL, thetaM, thetaO, thetaB, thetaE, nW

	int total = fpvp->nU + fpvp->nV + fpvp->nZ + fpvp->nC + fpvp->nA + fpvp->nH + fpvp->nL + fpvp->nM + fpvp->nO + fpvp->nB + fpvp->nE + fpvp->klrpv.nW;
	(*empirDist).resize(total);

	int k, i = 0;

	for (k = 0; k < fpvp->nU; k++)
		(*empirDist)[i++] = totalEmpirDistU[k];
	for (k = 0; k < fpvp->nV; k++)
		(*empirDist)[i++] = totalEmpirDistV[k];
	for (k = 0; k < fpvp->nZ; k++)
		(*empirDist)[i++] = totalEmpirDistZ[k];
	for (k = 0; k < fpvp->nC; k++)
		(*empirDist)[i++] = totalEmpirDistC[k];
	for (k = 0; k < fpvp->nA; k++)
		(*empirDist)[i++] = totalEmpirDistA[k];
	for (k = 0; k < fpvp->nH; k++)
		(*empirDist)[i++] = totalEmpirDistH[k];
	for (k = 0; k < fpvp->nL; k++)
		(*empirDist)[i++] = totalEmpirDistL[k];
	for (k = 0; k < fpvp->nM; k++)
		(*empirDist)[i++] = totalEmpirDistM[k];
	for (k = 0; k < fpvp->nO; k++)
	  (*empirDist)[i++] = totalEmpirDistO[k];
	for (k = 0; k < fpvp->nB; k++)
	    (*empirDist)[i++] = totalEmpirDistB[k];
	for (k = 0; k < fpvp->nE; k++)
	      (*empirDist)[i++] = totalEmpirDistE[k];
	for (k = 0; k < fpvp->klrpv.nW; k++)
		(*empirDist)[i++] = -totalEmpirDistW[k];  // note W's have to be negated


}




int baseModel::returnConcatModelDistVector(std::vector<float> *modelDist)
{
	// follow order thetaU, thetaV, thetaZ, thetaC, thetaA, thetaH, thetaL, thetaM, thetaO, thetaB, thetaE, nW

	int total = fpvp->nU + fpvp->nV + fpvp->nZ + fpvp->nC + fpvp->nA + fpvp->nH + fpvp->nL + fpvp->nM + fpvp->nO + fpvp->nB + fpvp->nE + fpvp->klrpv.nW;
	(*modelDist).resize(total);

	int k, i = 0;

	for (k = 0; k < fpvp->nU; k++)
		(*modelDist)[i++] = totalModelDistU[k];
	for (k = 0; k < fpvp->nV; k++)
		(*modelDist)[i++] = totalModelDistV[k];
	for (k = 0; k < fpvp->nZ; k++)
		(*modelDist)[i++] = totalModelDistZ[k];
	for (k = 0; k < fpvp->nC; k++)
		(*modelDist)[i++] = totalModelDistC[k];
	for (k = 0; k < fpvp->nA; k++)
		(*modelDist)[i++] = totalModelDistA[k];
	for (k = 0; k < fpvp->nH; k++)
		(*modelDist)[i++] = totalModelDistH[k];
	for (k = 0; k < fpvp->nL; k++)
		(*modelDist)[i++] = totalModelDistL[k];
	for (k = 0; k < fpvp->nM; k++)
		(*modelDist)[i++] = totalModelDistM[k];
	for (k = 0; k < fpvp->nO; k++)
	  (*modelDist)[i++] = totalModelDistO[k];
	for (k = 0; k < fpvp->nB; k++)
	    (*modelDist)[i++] = totalModelDistB[k];
	for (k = 0; k < fpvp->nE; k++)
	      (*modelDist)[i++] = totalModelDistE[k];
	for (k = 0; k < fpvp->klrpv.nW; k++)
		(*modelDist)[i++] = -totalModelDistW[k];

}


double baseModel::getRelDiffNorm(fvec gradTheta, fvec empirDist)
{
	fvec relDiff;
    // terminate if gradient is small enough percentage of empirdist
    vecEltWiseQuot(gradTheta, empirDist, relDiff);
    // printfVec(relDiff,  " relDiff  : ");

    double norm = 0;

    //WORK-THIS will need to generalize code below
    if (fpvp->learngradonly == 1 && fpvp->klrpv.klr==1)
    {
    	norm = vecNorm(relDiff, 0, fpvp->nV-1);
    } else
    {
    	norm = vecNorm(relDiff);
    }
    // DEBUG_OUT1(verbose, debugfile, "norm of gradient = %.10f%%\n", 100 * norm);
    // DEBUG_OUT1(verbose, debugfile, "old norm of gradient = %.10f%%\n", 100 * oldnorm);
    // DEBUG_OUT1(verbose, debugfile, "diff old norm - norm  = %.20f%%\n", 100 * oldnorm - 100 * norm);

    return (norm);
}


void baseModel::speedUpProgress(fvec gradTheta, double reldiffnorm, double rmsdiffDist)
{
	fvec theta;

	if (prm->featurepv.generative == 0)
	  oldloglikelihood = loglikelihood;
	else
	  oldloglikelihood = loglikelihoodgen;

	// if current solution is better than previous one,
	// we're fine, keep going with bigger and bigger steps in same direction
	returnConcatThetaWWVector(&theta);
	vecCopy(theta, oldTheta); // remember current point

//	if (prm.featurepv.klrpv.klr==1)
//		matrixCopy(wparam, oldwparam);

	// update theta
	vecScale(gradTheta, prm->graddescpv.rate);

	// need to see how to handle klr
	vecSum(theta, gradTheta, theta); // theta = theta + step * gradTheta
	splitThetaWWVectorToOriginal(theta);

	prm->graddescpv.rate *= 1.05f;

	LOGBUG(prm->verbose, debugfile, "\n reldiff = %g %g\n", oldreldiffnorm, reldiffnorm);

	oldreldiffnorm = reldiffnorm;
	oldrmsdiffDist = rmsdiffDist;

	LOGBUG(prm->verbose, debugfile, "\n speed up learning rate = %g \n", prm->graddescpv.rate);
	// DEBUG_OUT1(verbose, debugfile, "learning prm.graddescpv.rate = %g\n", prm.graddescpv.rate);

	// copy new result(dirDisp) to old temporary output variable(dirWTA)
	// for use in next iteration
	// copy from to to
	//copyImSet(dirDisp, dirWTA);
}



void baseModel::stepBackSlowDownProgress(double reldiffnorm)
{
	// if norm gets much bigger, go back half-way and try again with much smaller step size
	//BKUPrateflag = 1;
	// in case prm.featurepv.gradOnly=1, since we are not updating prm.featurepv.thetaU in other if condition,
	// the following two operatings will maintain same value of prm.featurepv.thetaU

	//			vecSum(theta, oldTheta, theta);  //theta = (theta + oldTheta)/2
	//			vecScale(theta, 0.5);

	fvec thetatemp, theta;

	returnConcatThetaWWVector(&theta);
	vecScale(theta, 0.001);
	vecCopy(oldTheta, thetatemp);
	vecScale(thetatemp, 0.999);
	vecSum(theta, thetatemp, theta);  //theta = 0.001*theta + 0.999*oldTheta
	splitThetaWWVectorToOriginal(theta);


	if (reldiffnorm>2*oldreldiffnorm)
		prm->graddescpv.rate *= 0.5f;
	else
		prm->graddescpv.rate *= 0.85f;


	if (prm->graddescpv.rate < 1e-12)
	  prm->graddescpv.rate = 1e-12;

	LOGBUG(prm->verbose, debugfile, "\n backing up learning rate = %g\n", prm->graddescpv.rate);
	// DEBUG_OUT1(verbose, debugfile, "*** BACKING UP, learning prm.graddescpv.rate = %g\n", prm.graddescpv.rate);

}


void baseModel::printfVecToLOG(fvec a, char *name)
{
  int nK = (int)a.size();
  LOGBUG(prm->verbose, debugfile, "%s: (%s", name, nK==0? ")\n" : "");
  for (int k = 0; k < nK; k++)
	  LOGBUG(prm->verbose, debugfile, "%g%s", a[k], k==nK-1? ")\n" : ",");
}


double baseModel::getrmsdiffDist(fvec a)
{
	double rms = 0.0;
	int nK = (int)a.size();
	for (int k = 0; k < nK; k++)
		rms += a[k]*a[k];
	rms = sqrt(rms);

	return(rms);
}


void baseModel::weightClassDistTemp(fvec &Dist)
{
	// order  fpvp->nU + fpvp->nV + fpvp->nZ + fpvp->nC + fpvp->nA + fpvp->nH + fpvp->nL + fpvp->nM + fpvp->klrpv.nW;
	// assume nW  not there


	// assume the results are not fractional so some bins/parameters are not lost
	int nbins = double(fpvp->nU)/prm->nD;
	int nT = fpvp->nA/prm->nD; //this gives centers per class
	int nHpC = fpvp->nH/prm->nD; // HoG vocabulary size
	int nMpC = fpvp->nM/prm->nD; // MOT vocabulary size
	int nLpC = fpvp->nL/prm->nD;  // this variable will be useful only for cube setting (loc==1) or loc==9

	double divi[6] = {263.65, 47.25, 4.8, 4.5, 1.0, 7.9};

	int center=0;

	if (fpvp->nU > 0)
	{
		for (int d=0; d<prm->nD; d++)
		{
			for (int b=0; b<nbins; b++)
			{
				Dist[center] = Dist[center]/divi[d];
				center++;
			}
		}
	}


	if (fpvp->nA > 0)
	{
		for (int d=0; d<prm->nD; d++)
		{
			for (int b=0; b<nT; b++)
			{
				Dist[center] = Dist[center]/divi[d];
				center++;
			}
		}
	}

	if (fpvp->nH > 0)
	{
		for (int d=0; d<prm->nD; d++)
		{
			for (int b=0; b<nHpC; b++)
			{
				Dist[center] = Dist[center]/divi[d];
				center++;
			}
		}
	}


	if (fpvp->nL > 0)
	{
		for (int d=0; d<prm->nD; d++)
		{
			for (int b=0; b<nLpC; b++)
			{
				Dist[center] = Dist[center]/divi[d];
				center++;
			}
		}
	}


	if (fpvp->nM > 0)
	{
		for (int d=0; d<prm->nD; d++)
		{
			for (int b=0; b<nMpC; b++)
			{
				Dist[center] = Dist[center]/divi[d];
				center++;
			}
		}
	}

	if (fpvp->nO > 0)
	  assert(false);

	 if (fpvp->nB > 0)
	    assert(false);

   if (fpvp->nE > 0)
      assert(false);

}



void baseModel::weightClassDist(fvec &Dist)
{
	// order  fpvp->nU + fpvp->nV + fpvp->nZ + fpvp->nC + fpvp->nA + fpvp->nH + fpvp->nL + fpvp->nM + fpvp->klrpv.nW;
	// assume nW  not there

	if (fpvp->klrpv.nW != 0)
	{
		std::cout << " klr not handled !\n";
		exit(1);
	}


	// assume the results are not fractional so some bins/parameters are not lost
	int nbins = double(fpvp->nU)/prm->nD;
	int nT = fpvp->nA/prm->nD; //this gives centers per class
	int nHpC = fpvp->nH/prm->nD; // HoG vocabulary size
	int nMpC = fpvp->nM/prm->nD; // MOT vocabulary size
	int nLpC = fpvp->nL/prm->nD;  // this variable will be useful only for cube setting (loc==1) or loc==9

	int cnter=0;

	if (fpvp->nU > 0)
	{
		for (int d=0; d<prm->nD; d++)
		{
			for (int b=0; b<nbins; b++)
			{
				Dist[cnter] = Dist[cnter]/classSize[d];
				cnter++;
			}
		}
	}

	cnter += fpvp->nV + fpvp->nZ + fpvp->nC;


	if (fpvp->nA > 0)
	{
		for (int d=0; d<prm->nD; d++)
		{
			for (int b=0; b<nT; b++)
			{
				Dist[cnter] = Dist[cnter]/classSize[d];
				cnter++;
			}
		}
	}

	if (fpvp->nH > 0)
	{
		for (int d=0; d<prm->nD; d++)
		{
			for (int b=0; b<nHpC; b++)
			{
				Dist[cnter] = Dist[cnter]/classSize[d];
				cnter++;
			}
		}
	}


	if (fpvp->nL > 0)
	{
		for (int d=0; d<prm->nD; d++)
		{
			for (int b=0; b<nLpC; b++)
			{
				Dist[cnter] = Dist[cnter]/classSize[d];
				cnter++;
			}
		}
	}


	if (fpvp->nM > 0)
	{
		for (int d=0; d<prm->nD; d++)
		{
			for (int b=0; b<nMpC; b++)
			{
				Dist[cnter] = Dist[cnter]/classSize[d];
				cnter++;
			}
		}
	}

	if (fpvp->nO > 0)
	  assert(false);

	 if (fpvp->nB > 0)
	    assert(false);

   if (fpvp->nE > 0)
      assert(false);


}

















int baseModel::makeNextThetaStep(int iterout, int iter)
{

	int breakflag = 0;

	double loglikelihoodt = loglikelihood;

	if (prm->featurepv.generative == 1)
	  loglikelihoodt = loglikelihoodgen;

  // temporary variables
  fvec empirDist, modelDist, diffDist, gradTheta, theta;

  returnConcatEmpirDistVector(&empirDist);
  returnConcatModelDistVector(&modelDist);


  if (prm->optClassAccuracy == 2)
	{
		weightClassDist(empirDist);
		weightClassDist(modelDist);
	}

  vecDiff(modelDist, empirDist, diffDist);
  vecSum(diffDist, regTheta, diffDist); // regularization addition

  double rmsdiffDist = getrmsdiffDist(diffDist);

  printfVecToLOG(empirDist, "empirDist (GT)");
  printfVecToLOG(modelDist, "modelDist     ");
  printfVecToLOG(diffDist,  "diffDist      ");
  LOGBUG(prm->verbose, debugfile, "rms diffDist    %g\n", rmsdiffDist);


  if (prm->graddescpv.bfgs_flag != 0) // for both bfgs and lbfgs
  {
    returnConcatThetaWWVector(&theta);
    if (iterout==0 && iter==0)
      nrfirstUpdate(loglikelihoodt, diffDist, theta);
    else  // otherwise you are now doing linesearch and follow procedure
    {
      int another_flag = nrlinesearchUpdate(loglikelihoodt, diffDist, theta); // theta gets updated
      if (another_flag == 0)
        breakflag = 1;
      splitThetaWWVectorToOriginal(theta);
    }
  }

  if (prm->graddescpv.grad_desc==1)
  {
    vecCopy(diffDist, gradTheta);
    double reldiffnorm = getRelDiffNorm(gradTheta, empirDist);


/*		if (100*reldiffnorm < 0.5) // terminate when norm of diff is less than 0.5 percent
			breakflag = 1;
		else if (oldreldiffnorm < 0)
		{
			speedUpProgress(gradTheta, reldiffnorm);
			oldrmsdiffDist = rmsdiffDist;
		}
		else if ((iter > 0 || iterout > 0) && (reldiffnorm <= oldreldiffnorm) && (oldrmsdiffDist <= rmsdiffDist))
		{
			speedUpProgress(gradTheta, reldiffnorm);
			oldrmsdiffDist = rmsdiffDist;
		}
//		else if (oldreldiffnorm < 0 || reldiffnorm <= oldreldiffnorm || 100*abs(oldreldiffnorm-reldiffnorm) < 0.1 )
//			speedUpProgress(gradTheta, reldiffnorm);
		else
			stepBackSlowDownProgress(reldiffnorm);
*/

    if (prm->logregpv.logreg == 1 || prm->logregpv.logregpl == 1)
    {
      if (100*reldiffnorm < 0.5) // terminate when norm of diff is less than 0.5 percent
        breakflag = 1;
      else if (iter==0 || oldloglikelihood <= loglikelihoodt)
      	speedUpProgress(gradTheta, reldiffnorm, rmsdiffDist);
      else
      	stepBackSlowDownProgress(reldiffnorm);
    } else
    {
			if (100*reldiffnorm < 0.5) // terminate when norm of diff is less than 0.5 percent
				breakflag = 1;
			else if (oldreldiffnorm < 0 || reldiffnorm <= oldreldiffnorm || 100*abs(oldreldiffnorm-reldiffnorm) < 0.1 )
				speedUpProgress(gradTheta, reldiffnorm, rmsdiffDist);
			else
				stepBackSlowDownProgress(reldiffnorm);
    }

	} else if (prm->graddescpv.grad_desc == 2)
  {
	  vecCopy(diffDist, gradTheta);
    double reldiffnorm = getRelDiffNorm(gradTheta, empirDist);

    if (rmsdiffDist < 0.001)
      breakflag = 1;
    else if (oldrmsdiffDist  < 0 || rmsdiffDist <= oldrmsdiffDist || abs(oldrmsdiffDist-rmsdiffDist) < 1)
      speedUpProgress(gradTheta, reldiffnorm, rmsdiffDist);
    else
      stepBackSlowDownProgress(reldiffnorm);

  }

  return (breakflag);

}

void baseModel::deletemrf(unsigned int index)
{
	// do nothing
}

void baseModel::regularizeLL()
{
	fvec theta;
	returnConcatThetaWWVector(&theta);

	if (fpvp->nO > 0)
	  assert(false);

	 if (fpvp->nB > 0)
	    assert(false);

	int numRegParams = 0;
	numRegParams = 	fpvp->nU + fpvp->nV + fpvp->nZ + fpvp->nC + fpvp->nA + fpvp->nH + fpvp->nL + fpvp->nM;

	double regLLpart = 0.0;

	for (int i=0; i<numRegParams; i++)
	{
		regLLpart += theta[i]*theta[i];
	}
	regLLpart /= (2*(fpvp->gaussSigma)*(fpvp->gaussSigma));

	loglikelihood -= regLLpart;

	if (prm->featurepv.generative == 1)
	  loglikelihoodgen -= regLLpart;

}

// gaussian regularizer can be used (assume klr cannot be used with regularizer for now)
void baseModel::regularizationUpdate()
{

    if (fpvp->gaussSigma!=0)
    {
    	if (fpvp->klrpv.klr == 1)
    		std::cout << " Warning : KLR parameters are not regularized !!! " << std::endl;

    	 if (fpvp->nO > 0)
    	    assert(false);

		// update loglikelihood
		regularizeLL();

		// update regTheta
		fvec theta;
		returnConcatThetaWWVector(&theta);
		regTheta.resize(theta.size());

		int numRegParams = 0;
		numRegParams = 	fpvp->nU + fpvp->nV + fpvp->nZ + fpvp->nC + fpvp->nA + fpvp->nH + fpvp->nL + fpvp->nM + fpvp->nO + fpvp->nB;

		for (int i=0; i<numRegParams; i++)
			regTheta[i] = theta[i]/(fpvp->gaussSigma*fpvp->gaussSigma);

		for (int i=numRegParams; i<theta.size(); i++) // for nW if there are any at the end
			regTheta[i] = 0.0;

		LOG(debugfile, "\n  Regularized Log Likelihood : %f \n", loglikelihood);
		LOG(debugfile, "  Regularized Log Likelihood Gen : %f \n", loglikelihoodgen);


    } else
    {
      fvec theta;
      returnConcatThetaWWVector(&theta);
      regTheta.resize(theta.size());
      for (int i=0; i<theta.size(); i++)
        regTheta[i] = 0.0;
    }


}


// the function tells us whether next iteration is to be run using all data or not (fullrun = 1, full data)
// and what the start index for the next iteration should be (startidx)
// remember anything between sgdlast and maxiter is run using full data
// and once it goes through the data in a online / stochastic manner, it will do one iteration over the full data
// and then continue in a stochastic manner.
void baseModel::sgdupdate(int &fullrun, int &startidx, int iter)
{

	if (fullrun == 1)
		startidx = -1;

    fullrun = 0;
    if (prm->graddescpv.sgdlast == 0 || prm->graddescpv.sgdskipV == 0)
    {
    	fullrun = 1;
    	startidx = 0;
    }

	if (iter >= (prm->graddescpv.maxiter  - prm->graddescpv.sgdlast))
	{
		fullrun = 1;
		startidx = 0;
	}

	if (prm->graddescpv.sgdskipV == startidx)
	{
		fullrun = 1;
		startidx = 0;
	}
	else
		startidx++;

}

int baseModel::skipvolume(int i, int startidx, int fullrun)
{
	int skip = 0;

	int shifti = i - startidx;
	if (fullrun == 1) {}
	else if ((shifti % (prm->graddescpv.sgdskipV + 1)) == 0) {}
	else
		skip = 1; // skip this volume

	return skip;
}


///// define optimize in baseModel

  // select gradient method
  // for each iteration
  // for each training instance
  // for each image
  // (if 2 model, then do image+1 slice as well)
  //   computeDataCost(array)
  //   computeSmoothnessCost(function)  {specific to crfs}
  //   (a) inference
  //   (b) computeEmpirDist (only in iteration 1)
  //   (c) computeModelDist


/// a, b and c are in the same function? for logreg models unlike crf model



void baseModel::searchParameters()
{

  time_t start, end1, end2;
  double timediff=0.0;

	std::vector<float> theta;

	srand(time(0));
	time(&start);

  loglikelihood = 0.0;
  loglikelihoodgen = 0.0;

  int another_flag = 1;

  initializeEmpirDistToZero();

  int startidx = 0;
  int fullrun = 0;
  if (prm->graddescpv.sgdlast == 0 || prm->graddescpv.sgdskipV == 0)
  	fullrun = 1;

    for (int iterout=0; iterout < prm->graddescpv.maxiter_out && another_flag==1; iterout++)
    {
    	resetalphaTo1();
    	resetflagaTo1();

    	for (int iter=0; iter < prm->graddescpv.maxiter && getflaga() == 1; iter++)
    	{

    		initializeModelDistToZero();
    		initializeVarToZero();

    		if (prm->crfpv.crfp == 1 && prm->crfpv.hidden == 1) // we need to update empir dist every time
    		  initializeEmpirDistToZero();

    		// for each training instance
        // i indexes 3D volume

    		if (iter >= (prm->graddescpv.maxiter  - prm->graddescpv.sgdlast))
    			fullrun = 1;

        for(unsigned int i = startidx; i < gtdirImage.size(); ++i)
        {
          // sgd
          if (skipvolume(i, startidx, fullrun)==1)
            continue;

          globalP.setOpflowVolumeIndex(i);

          //if (iter > (prm->graddescpv.maxiter  - 5) && prm->iopv.testDirIndexV[i]==1) {}
          //else if ((iter % 10 != 0) && prm->iopv.testDirIndexV[i]==1)
          //	continue;
          if (prm->graddescpv.bfgs_flag != 0)
          {
            LOG(debugfile, "\n\n Outer Iter : %d Inner Iter : %d   Image No: %d \n\n",iterout, iter, i );
          }
          else
          {
            LOG(debugfile, "\n\n Iteration : %d   Image No: %d \n\n",iter, i );
          }
          // updateCueIndex(i); -- need to do something about this

          // reset variables per volume
          initializeEvalVarPVToZero();

          if (prm->crfpv.crfp == 1 && prm->crfpv.hidden == 1 && prm->iopv.testDirIndexV[i]==0)
          {
            // this function includes modelProcess, evaluateResults, outputing results and updating
            // distributions
            hiddenModelProcess(i, iter, iterout);
          }
          else
          {

            modelProcess(i);  // allocated space for mrf and calls model_dist
            evaluateResults(i); // for calling evaldisps
            outputEvalutions(i); // dump evaluations

            if (prm->graddescpv.bfgs_flag == 0)
            { // dump images for gradient descent methods
              outputImageResults(i, iterout*1000 + iter); // dump output images
            }

            if (prm->iopv.testDirIndexV[i]==0) // training
              updateDist(i, iterout || iter );

            deletemrf(i);

          }
        } // for each training example

        // sgd
        sgdupdate(fullrun, startidx, iter);

        evaluateResults(); // global evaluation of all volumes
        outputEvalutions(); // dump global evaluations

    	dumpThetaParameters(1, 1);

        regularizationUpdate();


        int breakflag = 0;
        if (iter == prm->graddescpv.maxiter-1 && iterout == prm->graddescpv.maxiter_out-1)
        {
          // we don't update thetas in the last iteration.
        }
        else
          breakflag = makeNextThetaStep(iterout, iter); // thetas and wparams get updated

        if (prm->graddescpv.bfgs_flag != 0 && breakflag==1)
        { // dump images at the last inner iteration for bfgs or lbfgs
          for(unsigned int i = startidx; i < gtdirImage.size(); ++i)
          {
            outputImageResults(i, iterout*1000 + iter); // dump output images
          }
        }

    		time(&end2);
    		timediff = difftime(end2, start);
    		LOG(debugfile, " Iteration : %d  Time : %f \n", iter, timediff);

        if (breakflag == 1)
          break;

  	    // start next iteration with current solution
    	} // inner iteration loop


		// bfgs processing
    if (getflaga() == 100)
    {
      std::cout << "reached small alpha value " <<std::endl;
    	break;
    }

    if (iterout == prm->graddescpv.maxiter_out-1)
    {
      // we don't update thetas in the last iteration
      // so the output images are the same as the parameter values
      break;
    }
      // we don't update thetas in the last iteration.
    if (prm->graddescpv.bfgs_flag != 0)
    {
      // call concatenate again
      returnConcatThetaWWVector(&theta);

      if (prm->graddescpv.bfgs_flag == 1 && iterout==0)
      {
        nrStartNewDescent(theta); // update theta
      }
      else if (prm->graddescpv.bfgs_flag == 1)
      {
        another_flag = nrouteriterUpdate(theta);  // update theta
      }
      else if (prm->graddescpv.bfgs_flag == 2 && iterout==0)
      {
        nrupdatedklimited(iterout);
        nrStartNewDescentlimited(theta); // update theta
      } else if (prm->graddescpv.bfgs_flag == 2)
      {
        another_flag = nrouteriterUpdatelimited(theta, iterout);  // update theta
      }
      splitThetaWWVectorToOriginal(theta);
    }


    } // outer iteration loop


    // output final results


}


