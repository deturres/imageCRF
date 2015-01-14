/*
 * logreg.cpp
 *
 *  Created on: Mar 15, 2011
 *      Author: bhole
 */

#include "logreg.h"


logreg::logreg(Parameters* prm, std::vector<float> theta) : baseModel(prm, theta)
{


}


int logreg::preComputePairwiseFeatures(std::vector <std::vector <CImage2> > *inp)
{
  dirgrad.erase(dirgrad.begin(), dirgrad.end());

  int numP = (*inp).size();

  for (int i=0; i< numP; i++)
  {
	  dirgrad.push_back((std::vector<CByteImage>) 0);
  }

}


int logreg::preComputePairwiseFeatures(std::vector <std::vector <CByteImage> > *inp)
{
  dirgrad.erase(dirgrad.begin(), dirgrad.end());

  int numP = (*inp).size();

  for (int i=0; i< numP; i++)
  {
	  dirgrad.push_back((std::vector<CByteImage>) 0);
  }
}

void logreg::allocateDataCostSpace()
{
	int depth = gtdirImage[0].size();
	CShape sh = gtdirImage[0][0].Shape();
	int width = sh.width, height = sh.height;

	dCArray.push_back(new MRF::CostVal[width * height * prm->nD]);

	if (fpvp->generative==1)
	    dCArrayGen.push_back(new MRF::CostVal[width * height * prm->nD]);

}


void logreg::deAllocateDataCostSpace()
{
  // we expect to have only one!
  for(std::vector<MRF::CostVal *>::reverse_iterator b = dCArray.rbegin(); b != dCArray.rend(); ++b)
	  delete [] *b;

  if (fpvp->generative==1)
  {
	  for(std::vector<MRF::CostVal *>::reverse_iterator b = dCArrayGen.rbegin(); b != dCArrayGen.rend(); ++b)
		  delete [] *b;
  }

}



void logreg::allocateSmoothCostGlobalSpace()
{
	// no need for space
}


void logreg::deAllocateSmoothCostGlobalSpace()
{
	// no need for space
}




void logreg::modelProcess(unsigned int index)
{
  //std::cout << " \n Dis \n";
  logisticRegression(	0, index);
  //std::cout << " \n Gen \n";
  if (fpvp->generative == 1)
    logisticRegression( 1, index);

  /*
    if (fpvp->generative == 1)
    {
    	std::cout << " Needs to be implemented yet !!" << std::endl;
    	exit(1);
      // will have to look into log reg code as mentioned above.
 //     logisticRegression(	1, indirImage[i], hogdirImage[i], prm.nD, numTrainingPats, dirDispGen[i], prm.featurepv.thetaU, prm.featurepv.thetaA, prm.featurepv.thetaH, prm.featurepv.thetaL, prm.featurepv.featureCode, appdirImage[i], appclass, LocationMVclass, AppMVclass, appdirMVProb[i], intObj, locdirImage[i], hogdirMVProb[i], intNBclass, startSliceNo, endSliceNo, prm.featurepv.nbinsintNB, i, loglikelihoodgen, gtdirImage[i]);
    }
*/



}


void logreg::logisticRegression(int genparam, unsigned int index)
{

	MRF::CostVal badcost;

	badcost = getbadCost(genparam);

	// dCArray and dCArrayGen will be used to store the probability values in logistic regression instead
	// will need to see how to process dsiArrayVGen (derive equations etc)

	if (genparam==0) {
		if (dCArray.empty() == true)
			throw CError("call allocateDataCost first");
	} else if (genparam==1) {
		if (dCArrayGen.empty() == true)
			throw CError("call allocateDataCost first");
	}

	int* bestd = new int[prm->nD];

	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	int zslice = 0;

	if (timeseq == 0 && (fpvp->loc != 0 || fpvp->klrpv.locKlr ==2))
		zslice = startSliceNo[index];


	// dirDisp.resize(depth);
	errormap[index].resize(depth);

	int dsiIndex;

	for(unsigned int i = 0; i < gtdirImage[index].size(); ++i)
	{

		dsiIndex = 0; // since we are processing each image and not storing entire volume

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
				uchar *gtpix = &gtdirImage[index][i].Pixel(x, y, 0);

				int numbest = 0;
				MRF::CostVal bestval = badcost;
				ublas::vector<DBL_TYPE> f;

				std::vector<double> dsiValueSum(prm->nD);

				for (int d = 0; d < prm->nD; d++)
				{
					bestd[d] = d;
				}

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
					// if classAccu is on, then klr parameters must also be scaled.
					dsiValue += getDataCost(genparam, width, height, depth, index, i , y, x, d, zslice);


					// The cost of pixel p and label l is stored at dCArray[p*nLabels+l]
					// Since in logistic regression, I am not using E_d ie e(-E_d), I negate the values of dCValue so i can work on prob()
					// directly.
					if (genparam==0)
						(dCArray[0])[dsiIndex++] =  -(float)dsiValue;
					else if (genparam==1)
						(dCArrayGen[0])[dsiIndex++] =  -(float)dsiValue;

					dsiValueSum[d] = - (double) dsiValue;


					// dsiValue here is the cost, so we want the minimum
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

/*          if ((x == 6 && y == 10))// || (x == 10 && y == 10))
          {
            std::cout << "d" << d << " bestd " << bestd[0]  ;
            if (genparam == 0)
              std::cout << " datacost : " << (dCArray[0])[dsiIndex-d]  << "\n";
            if (genparam == 1)
              std::cout << " datagen : " << (dCArrayGen[0])[dsiIndex-d] << "\n";
          }
*/

					if (prm->timeseq==0)
					{
						unsigned short pix1 = indirImage[index][i].Pixel(x, y, 0);
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

				double dsiValueLogSum  = getLogSumExp(dsiValueSum);



				for (int d = prm->nD; d > 0; d--)
				{

					if (genparam==0)
						(dCArray[0])[dsiIndex-d] = (double)(dCArray[0])[dsiIndex-d] - dsiValueLogSum;
					else if (genparam==1)
						(dCArrayGen[0])[dsiIndex-d] = (double)(dCArrayGen[0])[dsiIndex-d] - dsiValueLogSum;

				}

				// dsiIndex is already pointing to the new pixel and output level index and so need to rewind
				int tempid = dsiIndex-prm->nD+gtpix[0]; //bestd[curr];

				double factor = 1.0;
				if (prm->optClassAccuracy == 1)
				{
					factor = 1.0/classSize[gtpix[0]];
				}

        // unknown pixels do not count during loglikelihood value
        if (prm->unknown == 1 && gtpix[0]==prm->nD)
          continue;

				if (genparam == 0)
				{
					loglikelihood += (dCArray[0])[tempid] * factor;
				}
				else if (genparam==1)
				{
					loglikelihoodgen += (dCArrayGen[0])[tempid] * factor;
				}

			}
		}

		// updating model distribution here only if it is discriminative and it is a training example
		if (genparam == 0)
		  if (prm->iopv.testDirIndexV[index] == 0)
		    computeModelDist(index, i, hogdescs, appdescs);

	}

  delete [] bestd;

}



void logreg::updateDist(int index, int fiterflag)
{

	// if 1st iteration then do for empirical otherwise if not first, don't do for empirical
	if (fiterflag==0)
		computeEmpirDist(index);

	// modelDist already updated in modelProcess code

}


void logreg::computeEmpirDist(int index)
{

	if (fpvp->nU > 0 && fpvp->intensity==1)
		computeEmpirDistU(index);

	if (fpvp->nA > 0 && (fpvp->app==1 || fpvp->app==2))
	  computeEmpirDistA(index);

	if (fpvp->nH > 0 && fpvp->hog==1)
	  computeEmpirDistH(index);

	// no learning currently if loc = 4 or 5
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

	if (fpvp->nB > 0 && fpvp->bias==1)
	  computeEmpirDistB(index);

	if (fpvp->nE > 0 && fpvp->volume==1)
	    computeEmpirDistE(index);

}


void logreg::computeModelDist(int index, int z)
{
	// for optimal memory usage, this computation is done in modelprocess() for logreg model
}


void logreg::computeModelDist(int index, int z, matrix<DBL_TYPE> &hogdescs, matrix<DBL_TYPE> &appdescs)
{
	if (fpvp->nU > 0 && fpvp->intensity==1)
		computeModelDistU(index, z);

	if (fpvp->nA > 0 && (fpvp->app==1 || fpvp->app==2))
	  computeModelDistA(index, z);

	if (fpvp->nH > 0 && fpvp->hog==1)
	  computeModelDistH(index, z);

	if (fpvp->nL > 0 && (fpvp->loc==1 || fpvp->loc==6 || fpvp->loc==8 || fpvp->loc==7 || fpvp->loc==9))
	  computeModelDistL(index, z);

	if (fpvp->nM > 0 && fpvp->mot==1)
	  computeModelDistM(index, z);

	if (fpvp->klrpv.klr==1)
		computeModelDistW(index, z, hogdescs, appdescs);

	if (fpvp->nO > 0 && fpvp->opflowlocal==1)
	  computeModelDistO(index, z);

	if (fpvp->nB > 0 && fpvp->bias==1)
	  computeModelDistB(index, z);

  if (fpvp->nE > 0 && fpvp->volume==1)
      computeModelDistE(index, z);

}




// grayscale and // RGB or Yuv or Lab
// does not current have interactive mode
void logreg::computeEmpirDistU(int index)
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

					if (prm->optClassAccuracy == 1)
					{
						uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
						if (prm->featurepv.onlyUV == 0)
						  totalEmpirDistU[thetaidYuv[0]] += 1.0/classSize[gtpix[0]];
						totalEmpirDistU[thetaidYuv[1]] += 1.0/classSize[gtpix[0]];
					} else
					{
					  if (prm->featurepv.onlyUV == 0)
					    totalEmpirDistU[thetaidYuv[0]]++;
						totalEmpirDistU[thetaidYuv[1]]++;
					}

				}
				else if (prm->timeseq==0)
				{
					pix1 = indirImage[index][z].Pixel(x, y, 0);
					// getCostIntensityDiscBins(pix1, nbins, di[x], fpvp->thetaU, thetaid);
					getCostIntensityDiscBins(pix1, nbins, binner, fpvp->rangeI,  di[x], fpvp->thetaU, thetaid);
					if (prm->optClassAccuracy == 1)
					{
						uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
						totalEmpirDistU[thetaid] += 1.0/classSize[gtpix[0]];
					} else
						totalEmpirDistU[thetaid] ++;
				}

				// for us the costs are the theta param, but in differentiation, they vanish,
				// so its just whether the feature is enabled or not.
				// move costs pointer to next pixel


			}
		}
	}

}



void logreg::computeEmpirDistA(int index)
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
                {
					if (prm->optClassAccuracy == 1)
					{
						uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
						totalEmpirDistA[thetaid] += 1.0/classSize[gtpix[0]];
					} else
						totalEmpirDistA[thetaid] ++;
                }

			}
		}
	}
}




void logreg::computeEmpirDistH(int index)
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
    			{
					if (prm->optClassAccuracy == 1)
					{
						uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
						totalEmpirDistH[thetaid] += 1.0/classSize[gtpix[0]];
					} else
						totalEmpirDistH[thetaid] ++;
    			}

    		}
    	}
    }
}


void logreg::computeEmpirDistM(int index)
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

				if (prm->optClassAccuracy == 1)
				{
					uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
					totalEmpirDistM[thetaid] += 1.0/classSize[gtpix[0]];
				} else
					totalEmpirDistM[thetaid] ++;

    		}
    	}
    }
}


void logreg::computeEmpirDistB(int index)
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




void logreg::computeEmpirDistE(int index)
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







void logreg::computeEmpirDistL(int index)
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
					if (prm->optClassAccuracy == 1)
					{
						uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
						totalEmpirDistL[thetaid] += 1.0/classSize[gtpix[0]];
					} else
						totalEmpirDistL[thetaid]++;

    			} else if (fpvp->loc==9 || (fpvp->loc==7 && timeseq == 1)) {  //CRF soft or hard assignment redefined

    				getCostLocationDiscGaussian2(LocationMVclass, z, zslice, x, y, fpvp->thetaL, di[x], fpvp->loc, fpvp->generative, prm->nD, thetaid);
					if (prm->optClassAccuracy == 1)
					{
						uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
						totalEmpirDistL[thetaid] += 1.0/classSize[gtpix[0]];
					}
						totalEmpirDistL[thetaid]++;
    			} else if ((fpvp->loc==7 && timeseq == 0)) {  //CRF hard assignment redefined

    				getCostLocationDiscGaussian3(index, z, x, y, di[x], fpvp->generative, prm->nD, thetaid);
					if (prm->optClassAccuracy == 1)
					{
						uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
						totalEmpirDistL[thetaid] += 1.0/classSize[gtpix[0]];
					} else
						totalEmpirDistL[thetaid]++;
    			}

    		}
    	}
    }
}


void logreg::computeEmpirDistW(int index)
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
void logreg::computeEmpirDistO(int index)
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








// grayscale and // RGB or Yuv or Lab
// does not current have interactive mode
void logreg::computeModelDistU(int index, int z)
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

	int pixelnD = 0;

	for (int y = 0; y < height; y++)
	{
	  uchar *di = &gtdirImage[index][z].Pixel(0, y, 0);
		for (int x = 0; x < width; x++)
		{
      // unknown pixels do not count during learning parameters
      if (prm->unknown == 1 && di[x]==prm->nD)
      {
        pixelnD += prm->nD;
        continue;
      }

			for (int d=0; d < prm->nD ; d++)
			{

				if (prm->timeseq==1)
				{
					pix1t = &indirImaget[index][z].Pixel(x, y, 0);
					// getCostIntensityDiscBins(pix1t, nbins, d, thetaU, thetaidrgb);
	         if (prm->featurepv.onlyUV == 1)
	           getCostIntensityDiscBinsuv(pix1t, nbins, d, fpvp->thetaU, thetaidYuv);
	         else
	           getCostIntensityDiscBinsYuv(pix1t, nbins, d, fpvp->thetaU, thetaidYuv);


					if (prm->optClassAccuracy == 1)
					{
						uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
						if (prm->featurepv.onlyUV == 0)
						  totalModelDistU[thetaidYuv[0]] += exp(-log(classSize[gtpix[0]]) + (dCArray[0])[pixelnD]);
						totalModelDistU[thetaidYuv[1]] += exp(-log(classSize[gtpix[0]]) + (dCArray[0])[pixelnD]);

					} else
					{
					  if (prm->featurepv.onlyUV == 0)
					    totalModelDistU[thetaidYuv[0]] += exp((dCArray[0])[pixelnD]);
						totalModelDistU[thetaidYuv[1]] += exp((dCArray[0])[pixelnD]);
					}

				}
				else if (prm->timeseq==0)
				{
					pix1 = indirImage[index][z].Pixel(x, y, 0);
					// getCostIntensityDiscBins(pix1, nbins, d, fpvp->thetaU, thetaid);
					getCostIntensityDiscBins(pix1, nbins, binner, fpvp->rangeI,  d, fpvp->thetaU, thetaid);

					if (prm->optClassAccuracy == 1)
					{
						uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
						totalModelDistU[thetaid] += exp(-log(classSize[gtpix[0]]) + (dCArray[0])[pixelnD]);

					} else
					{
						totalModelDistU[thetaid] += exp((dCArray[0])[pixelnD]);
					}
				}

				pixelnD ++;

				// for us the costs are the theta param, but in differentiation, they vanish,
				// so its just whether the feature is enabled or not.
				// move costs pointer to next pixel
			}

		}
	}

}


void logreg::computeModelDistA(int index, int z)
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

    int pixelnD = 0;

	for (int y = 0; y < height; y++)
	{
	  uchar *di = &gtdirImage[index][z].Pixel(0, y, 0);
		for (int x = 0; x < width; x++)
		{
      // unknown pixels do not count during learning parameters
      if (prm->unknown == 1 && di[x]==prm->nD)
      {
        pixelnD += prm->nD;
        continue;
      }

			for (int d=0; d < prm->nD; d++)
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
				{
					if (prm->optClassAccuracy == 1)
					{
						uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
						totalModelDistA[thetaid] += exp(-log(classSize[gtpix[0]]) + (dCArray[0])[pixelnD]);
					} else
						totalModelDistA[thetaid] += exp((dCArray[0])[pixelnD]);
				}

			    pixelnD ++;

			}
		}
	}
}




void logreg::computeModelDistH(int index, int z)
{

	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nH == 0)
		return; // nothing to do

	int nHpC = fpvp->nH/prm->nD; // HoG vocabulary size

    int thetaid;

    int pixelnD = 0;


	for (int y = 0; y < height; y++)
	{
	  uchar *di = &gtdirImage[index][z].Pixel(0, y, 0);
		for (int x = 0; x < width; x++)
		{
      // unknown pixels do not count during learning parameters
      if (prm->unknown == 1 && di[x]==prm->nD)
      {
        pixelnD += prm->nD;
        continue;
      }

			for (int d=0; d < prm->nD; d++)
			{
				unsigned short hogval = hogdirImage[index][z].Pixel(x, y, 0);
    			getCostHoGDiscBins((int) hogval, d, nHpC, fpvp->thetaH, thetaid);
				//getCostHoGDiscBinsTemp((int) hogval, d, fpvp->thetaH, thetaid);

				if (thetaid>=0)
				{
					if (prm->optClassAccuracy == 1)
					{
						uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
						totalModelDistH[thetaid] += exp(-log(classSize[gtpix[0]]) + (dCArray[0])[pixelnD]);
					} else
						totalModelDistH[thetaid] += exp((dCArray[0])[pixelnD]);

				}
    			pixelnD ++;

    		}
    	}
    }
}


void logreg::computeModelDistM(int index, int z)
{
	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nM == 0)
		return; // nothing to do

	int nMpC = fpvp->nM/prm->nD; // Mot vocabulary size

    int thetaid;

    int pixelnD = 0;

	for (int y = 0; y < height; y++)
	{
	  uchar *di = &gtdirImage[index][z].Pixel(0, y, 0);
		for (int x = 0; x < width; x++)
		{
      // unknown pixels do not count during learning parameters
      if (prm->unknown == 1 && di[x]==prm->nD)
      {
        pixelnD += prm->nD;
        continue;
      }

			for (int d=0; d < prm->nD; d++)
			{
				unsigned short motval = motdirImage[index][z].Pixel(x, y, 0);
    			getCostmotDiscBins((int) motval, d, nMpC, fpvp->thetaM, thetaid);

				if (prm->optClassAccuracy == 1)
				{
					uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
					totalModelDistM[thetaid] += exp(-log(classSize[gtpix[0]]) + (dCArray[0])[pixelnD]);
				} else
					totalModelDistM[thetaid] += exp((dCArray[0])[pixelnD]);


    			pixelnD ++;

    		}
    	}
    }
}


void logreg::computeModelDistB(int index, int z)
{

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nB == 0)
    return; // nothing to do

  int thetaid;
  int pixelnD = 0;

  for (int y = 0; y < height; y++)
  {
    uchar *di = &gtdirImage[index][z].Pixel(0, y, 0);
    for (int x = 0; x < width; x++)
    {
      // unknown pixels do not count during learning parameters
      if (prm->unknown == 1 && di[x]==prm->nD)
      {
        pixelnD += prm->nD;
        continue;
      }

      for (int d=0; d < prm->nD; d++)
      {
        getCostBias(d, fpvp->thetaB, thetaid);

        if (thetaid>=0)
        {
          if (prm->optClassAccuracy == 1)
          {
            uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
            totalModelDistB[thetaid] += exp(-log(classSize[gtpix[0]]) + (dCArray[0])[pixelnD]);
          } else
            totalModelDistB[thetaid] += exp((dCArray[0])[pixelnD]);

        }
        pixelnD ++;

      }
    }
  }
}




void logreg::computeModelDistE(int index, int z)
{

  int depth = gtdirImage[index].size();
  CShape sh = gtdirImage[index][0].Shape();
  int width = sh.width, height = sh.height;

  if (fpvp->nE == 0)
    return; // nothing to do

  int thetaid;
  int pixelnD = 0;

  for (int y = 0; y < height; y++)
  {
    uchar *di = &gtdirImage[index][z].Pixel(0, y, 0);
    for (int x = 0; x < width; x++)
    {
      // unknown pixels do not count during learning parameters
      if (prm->unknown == 1 && di[x]==prm->nD)
      {
        pixelnD += prm->nD;
        continue;
      }

      for (int d=0; d < prm->nD; d++)
      {
        getCostInverseClassSize(d, inverseClassSize, fpvp->thetaE, thetaid);

        if (thetaid>=0)
        {
          if (prm->optClassAccuracy == 1)
          {
            uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
            totalModelDistE[thetaid] += exp(-log(classSize[gtpix[0]]) + (dCArray[0])[pixelnD]);
          } else
            totalModelDistE[thetaid] += exp((dCArray[0])[pixelnD]) * inverseClassSize[d];

        }
        pixelnD ++;

      }
    }
  }
}





void logreg::computeModelDistL(int index, int z)
{

	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nL == 0)
		return; // nothing to do

	int nLpC = fpvp->nL/prm->nD;  // this variable will be useful only for cube setting (loc==1) or loc==9

    int thetaid;

    int zslice=0;
    //zslice = startSliceNo[index]; // - WORK on this

    int pixelnD = 0;

	for (int y = 0; y < height; y++)
	{
	  uchar *di = &gtdirImage[index][z].Pixel(0, y, 0);
		for (int x = 0; x < width; x++)
		{
      // unknown pixels do not count during learning parameters
      if (prm->unknown == 1 && di[x]==prm->nD)
      {
        pixelnD += prm->nD;
        continue;
      }

			for (int d=0; d < prm->nD; d++)
			{
    			if (fpvp->loc==1)
    			{
    				/*
    				unsigned short locpix = locdirImage[index][di[x]][z].Pixel(x,y,0);
    				getCostLocationDiscCube(locpix, x, y, z, fpvp->loccubex, fpvp->loccubey, fpvp->loccubez, numTrainingPats, nLpC, fpvp->thetaL, width, height, depth, d, thetaid);
    				double locpixd = (double)locpix/ (fpvp->loccubex * fpvp->loccubey * fpvp->loccubez * numTrainingPats);

    				totalModelDistL[thetaid] += exp((dsiArrayV[index])[pixelnD]) * (1-locpixd);

    				// need to check this for boundary conditions especially when loccube? not multiple of size
					// or when it is actually multiple of size of image volume

					// int thetaid = nLpC*d + hno*maxwidth*maxheight + rowno*maxwidth + colno;

    				*/

    				std::cout << " Need to fix this - location " << std::endl;
    				exit(1);


    			} else if (fpvp->loc==6 || fpvp->loc==8) { // CRF hard assignment and soft

    				getCostLocationDiscGaussian(LocationMVclass, z, zslice, x, y, fpvp->thetaL, d, fpvp->loc, fpvp->generative, thetaid);
					if (prm->optClassAccuracy == 1)
					{
						uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
						totalModelDistL[thetaid] += exp(-log(classSize[gtpix[0]]) + (dCArray[0])[pixelnD]);
					} else
						totalModelDistL[thetaid] += exp((dCArray[0])[pixelnD]);


    			} else if (fpvp->loc==9 || (fpvp->loc==7 && timeseq == 1)) {  //CRF soft or hard assignment redefined

    				getCostLocationDiscGaussian2(LocationMVclass, z, zslice, x, y, fpvp->thetaL, d, fpvp->loc, fpvp->generative, prm->nD, thetaid);
					if (prm->optClassAccuracy == 1)
					{
						uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
						totalModelDistL[thetaid] += exp(-log(classSize[gtpix[0]]) + (dCArray[0])[pixelnD]);
					} else
						totalModelDistL[thetaid] += exp((dCArray[0])[pixelnD]);


    			} else if ((fpvp->loc==7 && timeseq == 0)) {  //CRF hard assignment redefined

    				getCostLocationDiscGaussian3(index, z, x, y, d, fpvp->generative, prm->nD, thetaid);
					if (prm->optClassAccuracy == 1)
					{
						uchar *gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
						totalModelDistL[thetaid] += exp(-log(classSize[gtpix[0]]) + (dCArray[0])[pixelnD]);
					} else
						totalModelDistL[thetaid] += exp((dCArray[0])[pixelnD]);

    			}

    			pixelnD ++;

    		}
    	}
    }
}


void logreg::computeModelDistW(int index, int z, matrix<DBL_TYPE> &hogdescs, matrix<DBL_TYPE> &appdescs)
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

	ublas::vector<DBL_TYPE> Xtest(Xtrain.size2());
	ublas::vector<DBL_TYPE> f;

	int pixelnD = 0;

	uchar *pix1t = 0;
	unsigned short pix1 = 0;

	for (int y = 0; y < height; y++)
	{
		uchar *di = &gtdirImage[index][z].Pixel(0, y, 0);

		for (int x = 0; x < width; x++)
		{

      // unknown pixels do not count during learning parameters
      if (prm->unknown == 1 && di[x]==prm->nD)
      {
        pixelnD += prm->nD;
        continue;
      }

			getDataCostKlr(f, Xtest, index, z, x, y, zslice, hogdescs, appdescs, width, height);

			for (int j=0; j<Ntr; j++)
			{
				double Kval = ivm->classify_klrexp(&Xtest, j);

				for (int d=0; d<prm->nD; d++)
				{
					int thetaid = j*dimW + d;

					totalModelDistW[thetaid] += exp((dCArray[0])[pixelnD])*Kval;

					pixelnD ++;
				}
			}
		}
	}

}









void logreg::computeModelDistO(int index, int z)
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
  if (z < depth-1)
  {
    int pixelnD = 0;

    for (int y = 0; y < height; y++)
    {
      uchar *di = &gtdirImage[index][z].Pixel(0, y, 0);
      for (int x = 0; x < width; x++)
      {
        // unknown pixels do not count during learning parameters
        if (prm->unknown == 1 && di[x]==prm->nD)
        {
          pixelnD += prm->nD;
          continue;
        }

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
              totalModelDistO[thetaido[0]] += exp((dCArray[0])[pixelnD]) * u;
              totalModelDistO[thetaido[1]] += exp((dCArray[0])[pixelnD]) * v;
            }
            else // binned optical flow features
            {
              getCostOpflowBinsuv(u, v, nOpC, d, fpvp->thetaO, thetaid, fpvp->opflowlocalframes.size(), idx);
              totalModelDistO[thetaid] += exp((dCArray[0])[pixelnD]);
            }
          }

          pixelnD ++;

          // consider the costs for the theta param,
          // since they exist after differentiation
          // no use of the indcost since they are the theta values
        }

      }
    }

  }
}






