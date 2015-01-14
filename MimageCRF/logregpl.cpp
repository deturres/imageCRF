/*
 * logregpl.cpp
 *
 *  Created on: Mar 15, 2011
 *      Author: bhole
 */

#include "logregpl.h"


logregpl::logregpl(Parameters* prm, std::vector<float> theta) : baseModel(prm, theta)
{
	mrf=0;
    mrfgen=0;
}


logregpl::~logregpl()
{

}

void logregpl::deletemrf(unsigned int index)
{
	if (prm->iopv.testDirIndexV[index]==1)
		delete mrf;
//	if (fpvp->generative==1)
//		delete mrfgen;
}


void logregpl::allocateDataCostSpace()
{
	int depth = gtdirImage[0].size();
	CShape sh = gtdirImage[0][0].Shape();
	int width = sh.width, height = sh.height;

	dCArray.push_back(new MRF::CostVal[width * height * prm->nD]);

	if (fpvp->generative==1)
	    dCArrayGen.push_back(new MRF::CostVal[width * height * prm->nD]);

}


void logregpl::deAllocateDataCostSpace()
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



void logregpl::allocateSmoothCostGlobalSpace()
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


void logregpl::deAllocateSmoothCostGlobalSpace()
{

  for(std::vector<MRF::CostVal *>::reverse_iterator b = globalP.horPairGlobalF.rbegin(); b != globalP.horPairGlobalF.rend(); ++b)
	  delete [] *b;

  for(std::vector<MRF::CostVal *>::reverse_iterator b = globalP.verPairGlobalF.rbegin(); b != globalP.verPairGlobalF.rend(); ++b)
	  delete [] *b;

  for(std::vector<MRF::CostVal *>::reverse_iterator b = globalP.depPairGlobalF.rbegin(); b != globalP.depPairGlobalF.rend(); ++b)
	  delete [] *b;

}



void logregpl::modelProcess(unsigned int index)
{

    if (prm->iopv.testDirIndexV[index]==1)
    	modelProcessTest(0, index);
    if (prm->iopv.testDirIndexV[index]==0)
    	modelProcessTrain(0, index);


    if (fpvp->generative == 1)
    {
    	std::cout << " Needs to be implemented yet !!" << std::endl;
    	exit(1);
    }
}




double logregpl::getSmoothCost(unsigned int index, unsigned int z, int y, int x, int d)
{

	double dsiValuepq = 0;

	int d2, d1;

	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;


    // l,r,u, d == left, right, up, down = 2D connections
	// t, b == top, bottom = 3D connections

    uchar *gtpixql = 0;
    uchar *gtpixqr = 0;
    uchar *gtpixqu = 0;
    uchar *gtpixqd = 0;
    uchar *gtpixqt = 0;
    uchar *gtpixqb = 0;

    uchar *gradl = 0;
    uchar *gradr = 0;
    uchar *gradd = 0;
    uchar *gradu = 0;
    uchar *gradt = 0;
    uchar *gradb = 0;


    if (x>0)
    	gtpixql = &gtdirImage[index][z].Pixel(x-1, y, 0);

    if (x<width-1)
    	gtpixqr = &gtdirImage[index][z].Pixel(x+1, y, 0);

    if (y>0)
    	gtpixqu = &gtdirImage[index][z].Pixel(x, y-1, 0);

    if (y<height-1)
    	gtpixqd = &gtdirImage[index][z].Pixel(x, y+1, 0);

    if (z>0)
    	gtpixqt = &gtdirImage[index][z-1].Pixel(x, y, 0);

    if (z>depth-1)
    	gtpixqb = &gtdirImage[index][z+1].Pixel(x, y, 0);


    if (fpvp->csgrad == 0)
    {

        int gr = 0;
        int gd = 0;
        int gl = 0;
        int gu = 0;
        int gt = 0;
        int gb = 0;

        // process left
    	if (x>0)
    	{
    		gradl = &dirgrad[index][z].Pixel(x-1, y, 0);
    		gl = gradl[0];

    		d1 = d;
    		d2 = gtpixql[0];
    		reorder(&d1, &d2);
    		dsiValuepq +=  fpvp->thetaV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gl];
    	}

    	// process right
    	if (x<width-1)
    	{
    		gradr = &dirgrad[index][z].Pixel(x, y, 0);
    		gr = gradr[0];

    		d1 = d;
    		d2 = gtpixqr[0];
    		reorder(&d1, &d2);
    		dsiValuepq +=  fpvp->thetaV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gr];
    	}

    	// process up
        if (y>0)
        {
        	gradu = &dirgrad[index][z].Pixel(x, y-1, 1);
        	gu = gradu[0];

    		d1 = d;
    		d2 = gtpixqu[0];
    		reorder(&d1, &d2);
    		dsiValuepq +=  fpvp->thetaV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gu];
        }

        // process down
    	if (y<height-1)
    	{
    		gradd = &dirgrad[index][z].Pixel(x, y, 1);
    		gd = gradd[0];

    		d1 = d;
    		d2 = gtpixqd[0];
    		reorder(&d1, &d2);
    		dsiValuepq +=  fpvp->thetaV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gd];
    	}

    	// process top
    	if (fpvp->gradVZ == 2 || fpvp->gradVZ == 3)
        if (z>0)
        {
        	gradt = &dirgrad[index][z-1].Pixel(x, y, 2);
        	gt = gradt[0];

    		d1 = d;
    		d2 = gtpixqt[0];
    		reorder(&d1, &d2);
    		if (fpvp->gradVZ == 1)
    		{}
    		else if (fpvp->gradVZ == 2)
    			dsiValuepq +=  fpvp->thetaV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gt];
    		else if (fpvp->gradVZ == 3)
    			dsiValuepq +=  fpvp->thetaZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nGZ + gt];

        }

        // process bottom
    	if (fpvp->gradVZ == 2 || fpvp->gradVZ == 3)
    	if (z<depth-1)
    	{
    		gradb = &dirgrad[index][z].Pixel(x, y, 2);
    		gb = gradb[0];

    		d1 = d;
    		d2 = gtpixqb[0];
    		reorder(&d1, &d2);
    		if (fpvp->gradVZ == 1)
    		{}
    		else if (fpvp->gradVZ == 2)
    			dsiValuepq +=  fpvp->thetaV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gb];
    		else if (fpvp->gradVZ == 3)
    			dsiValuepq +=  fpvp->thetaZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nGZ + gb];
    	}

    }
    else if (fpvp->csgrad == 1)
    {
    	std::cout << " This needs to be done !" << std::endl;
    	exit(1);
    }


    return dsiValuepq;

}


void logregpl::modelProcessTrain(int genparam, unsigned int index)
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

					dsiValue += getDataCost(genparam, width, height, depth, index, i , y, x, d, zslice);

					dsiValue += getSmoothCost(index, i , y, x, d);

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


				if (genparam == 0)
				{
					loglikelihood += (dCArray[0])[tempid] *factor;
				}
				else if (genparam==1)
				{
					loglikelihoodgen += (dCArrayGen[0])[tempid] * factor;
				}

			}
		}

		// updating model distribution here
		if (prm->iopv.testDirIndexV[index] == 0)
			computeModelDist(index, i, hogdescs, appdescs);

	}

  delete [] bestd;

}




void logregpl::modelProcessTest(int genparam, unsigned int index)
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
    if (prm->iopv.testDirIndexV[index]==1 && prm->crfpv.inferencerTest == USE_GC)
    {
    	// make sure submodularity is valid
    	int cRC = checkRegularityCondition(prm->featurepv.thetaV, prm->featurepv.gradThreshVec, prm->nD);
    	int cRC2 = checkRegularityCondition(prm->featurepv.thetaZ, prm->featurepv.gradThreshVecZ, prm->nD);
    	if (cRC == 0 || cRC2 == 0)
    	{
    		std::cout << " cRC failure \n";
    		exit(1);
    	}

    }

    // interactive mode - not coded (look in older code)
    if (fpvp->gradContext==1 && prm->interactive==1)
    {
		std::cout << " interactive failure - not coded (look in older code) \n";
		exit(1);
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
		std::cout << " generative in logregpl - not coded (look in older code) \n";
		exit(1);
	}

    delete scost;
    delete dcost;


}



DataCost* logregpl::computeDataCost(int genparam, unsigned int index)
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

	int dsiIndex = 0;

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










SmoothnessCost* logregpl::computeSmoothnessCost(int index)
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


	// return function pointer
	if (fpvp->csgrad == 0)
		return new SmoothnessCost((MRF::SmoothCostGeneralFn) fnCostQuantizedGlobal);
	else if (fpvp->csgrad == 1)
		return new SmoothnessCost((MRF::SmoothCostGeneralFn) fnCostSmoothGlobal);

}

void logregpl::updateDist(int index, int fiterflag)
{

	// if 1st iteration then do for empirical otherwise if not first, don't do for empirical
	if (fiterflag==0)
		computeEmpirDist(index);

	// modelDist already updated in modelProcess code

}


void logregpl::computeEmpirDist(int index)
{

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

	// when nV == 1, we are not learning parameters (because constrast sensitive smoothing)
	// but in the function we still need to check again because nZ could be 1 , but not nV and so on
	if (fpvp->nV > 1 || fpvp->nC > 0 || fpvp-> nZ > 1)
		computeEmpirDistV(index);

}






void logregpl::computeModelDist(int index, int z)
{
	// for optimal memory usage, this computation is done in modelprocess() for logregpl model
}


void logregpl::computeModelDist(int index, int z, matrix<DBL_TYPE> &hogdescs, matrix<DBL_TYPE> &appdescs)
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

	// when nV == 1, we are not learning parameters (because constrast sensitive smoothing)
	// but in the function we still need to check again because nZ could be 1 , but not nV and so on
	if (fpvp->nV > 1 || fpvp->nC > 0 || fpvp-> nZ > 1)
		computeModelDistV(index, z);

}





// grayscale and // RGB or Yuv or Lab
// does not current have interactive mode
void logregpl::computeEmpirDistU(int index)
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
					getCostIntensityDiscBins(pix1, nbins, di[x], fpvp->thetaU, thetaid);
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



void logregpl::computeEmpirDistA(int index)
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




void logregpl::computeEmpirDistH(int index)
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


void logregpl::computeEmpirDistM(int index)
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




void logregpl::computeEmpirDistL(int index)
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

    				getCostLocationDiscGaussian3(index, z,  x, y,  di[x],  fpvp->generative, prm->nD,  thetaid);
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


void logregpl::computeEmpirDistW(int index)
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





// this will not learn parameters related to smoothed contrast sensitive gradients
void logregpl::computeEmpirDistV(int index)
{
	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nV == 0 && fpvp->nC == 0 && fpvp->nZ == 0)
		return; // nothing to do

	if (prm->interactive==1)
	{
		std::cout << " interactive not implemented yet" << std::endl;
		exit(1);
	}

	if (fpvp->gradContext==0)
	{
		std::cout << " Not  working really check this !!! \n";
		exit(1);
	}


	int depthrun = 0;

	if (fpvp->gradVZ == 1)
		depthrun = depth;
	else if (fpvp->gradVZ == 2 || fpvp->gradVZ == 3)
		depthrun = depth-1;


	depthrun = depth;


	int d1, d2;

	for (int z=0; z < depthrun; z++)
	{

		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{

				// l,r,u, d == left, right, up, down = 2D connections
				// t, b == top, bottom = 3D connections

				uchar *gtpixql = 0;
				uchar *gtpixqr = 0;
				uchar *gtpixqu = 0;
				uchar *gtpixqd = 0;
				uchar *gtpixqt = 0;
				uchar *gtpixqb = 0;

				uchar *gradl = 0;
				uchar *gradr = 0;
				uchar *gradd = 0;
				uchar *gradu = 0;
				uchar *gradt = 0;
				uchar *gradb = 0;


				if (x>0)
					gtpixql = &gtdirImage[index][z].Pixel(x-1, y, 0);

				if (x<width-1)
					gtpixqr = &gtdirImage[index][z].Pixel(x+1, y, 0);

				if (y>0)
					gtpixqu = &gtdirImage[index][z].Pixel(x, y-1, 0);

				if (y<height-1)
					gtpixqd = &gtdirImage[index][z].Pixel(x, y+1, 0);

				if (z>0)
					gtpixqt = &gtdirImage[index][z-1].Pixel(x, y, 0);

				if (z<depth-1)
					gtpixqb = &gtdirImage[index][z+1].Pixel(x, y, 0);



				uchar *gtpix = 0;
				gtpix = &gtdirImage[index][z].Pixel(x, y, 0);
				int d = gtpix[0];

        if (prm->unknown == 1 && d==prm->nD)
          continue;


				double factor = 1.0;
				if (prm->optClassAccuracy == 1)
					factor = 1.0/classSize[d];


				if (fpvp->csgrad == 0)
				{

					int gr = 0;
					int gd = 0;
					int gl = 0;
					int gu = 0;
					int gt = 0;
					int gb = 0;

					// process left
					if (x>0)
					{
						gradl = &dirgrad[index][z].Pixel(x-1, y, 0);
						gl = gradl[0];

						d1 = d;
						d2 = gtpixql[0];

		        if (prm->unknown == 1 && d2==prm->nD) {}
		        else
		        {
              reorder(&d1, &d2);
              totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gl] += factor;
		        }
					}

					// process right
					if (x<width-1)
					{
						gradr = &dirgrad[index][z].Pixel(x, y, 0);
						gr = gradr[0];

						d1 = d;
						d2 = gtpixqr[0];
            if (prm->unknown == 1 && d2==prm->nD) {}
            else
            {
              reorder(&d1, &d2);
              totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gr] += factor;
            }
					}

					// process up
					if (y>0)
					{
						gradu = &dirgrad[index][z].Pixel(x, y-1, 1);
						gu = gradu[0];

						d1 = d;
						d2 = gtpixqu[0];
            if (prm->unknown == 1 && d2==prm->nD) {}
            else
            {
              reorder(&d1, &d2);
              totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gu] += factor;
            }
					}

					// process down
					if (y<height-1)
					{
						gradd = &dirgrad[index][z].Pixel(x, y, 1);
						gd = gradd[0];

						d1 = d;
						d2 = gtpixqd[0];
            if (prm->unknown == 1 && d2==prm->nD) {}
            else
            {
              reorder(&d1, &d2);
              totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gd] += factor;
            }
					}

					// process top
					if (fpvp->gradVZ == 2 || fpvp->gradVZ == 3)
					if (z>0 && z<depth)
					{
						gradt = &dirgrad[index][z-1].Pixel(x, y, 2);
						gt = gradt[0];

						d1 = d;
						d2 = gtpixqt[0];

            if (prm->unknown == 1 && d2==prm->nD) {}
            else
            {
              reorder(&d1, &d2);
              if (fpvp->gradVZ == 1)
              {}
              else if (fpvp->gradVZ == 2 && fpvp->nZ > 1)
                totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gt] += factor;
              else if (fpvp->gradVZ == 3 && fpvp->nZ > 1)
                totalEmpirDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nGZ + gt]  += factor;
            }

					}

					// process bottom
					if (fpvp->gradVZ == 2 || fpvp->gradVZ == 3)
					if (z>= 0 && z<depth-1)
					{
						gradb = &dirgrad[index][z].Pixel(x, y, 2);
						gb = gradb[0];

						d1 = d;
						d2 = gtpixqb[0];
            if (prm->unknown == 1 && d2==prm->nD) {}
            else
            {
              reorder(&d1, &d2);
              if (fpvp->gradVZ == 1)
              {}
              else if (fpvp->gradVZ == 2 && fpvp->nZ > 1)
                totalEmpirDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gb] += factor;
              else if (fpvp->gradVZ == 3 && fpvp->nZ > 1)
                totalEmpirDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nGZ + gb]  += factor;
            }
					}

				}
				else if (fpvp->csgrad == 1)
				{
					std::cout << " This needs to be done !" << std::endl;
					exit(1);
				}
			}
		}
	}
}







// grayscale and // RGB or Yuv or Lab
// does not current have interactive mode
void logregpl::computeModelDistU(int index, int z)
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
        continue;

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
					getCostIntensityDiscBins(pix1, nbins, d, fpvp->thetaU, thetaid);

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


void logregpl::computeModelDistA(int index, int z)
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
        continue;

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




void logregpl::computeModelDistH(int index, int z)
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
        continue;

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


void logregpl::computeModelDistM(int index, int z)
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
        continue;

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




void logregpl::computeModelDistL(int index, int z)
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
        continue;

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


void logregpl::computeModelDistW(int index, int z, matrix<DBL_TYPE> &hogdescs, matrix<DBL_TYPE> &appdescs)
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
        continue;

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

// this will not learn parameters related to smoothed contrast sensitive gradients
void logregpl::computeModelDistV(int index, int z)
{
	int depth = gtdirImage[index].size();
	CShape sh = gtdirImage[index][0].Shape();
	int width = sh.width, height = sh.height;

	if (fpvp->nV == 0 && fpvp->nC == 0 && fpvp->nZ == 0)
		return; // nothing to do

	if (prm->interactive==1)
	{
		std::cout << " interactive not implemented yet" << std::endl;
		exit(1);
	}

	if (fpvp->gradContext==0)
	{
		std::cout << " Not  working really check this !!! \n";
		exit(1);
	}

	int pixelnD = 0;
	int d1, d2;

	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{

			// l,r,u, d == left, right, up, down = 2D connections
			// t, b == top, bottom = 3D connections

			uchar *gtpixql = 0;
			uchar *gtpixqr = 0;
			uchar *gtpixqu = 0;
			uchar *gtpixqd = 0;
			uchar *gtpixqt = 0;
			uchar *gtpixqb = 0;

			uchar *gradl = 0;
			uchar *gradr = 0;
			uchar *gradd = 0;
			uchar *gradu = 0;
			uchar *gradt = 0;
			uchar *gradb = 0;

			uchar *gtpixc = &gtdirImage[index][z].Pixel(x, y, 0);
			int dc = gtpixc[0];

      if (prm->unknown == 1 && dc==prm->nD)
        continue;

			double factor = 1.0;
			if (prm->optClassAccuracy == 1)
				factor = 1.0/classSize[dc];


			if (x>0)
				gtpixql = &gtdirImage[index][z].Pixel(x-1, y, 0);

			if (x<width-1)
				gtpixqr = &gtdirImage[index][z].Pixel(x+1, y, 0);

			if (y>0)
				gtpixqu = &gtdirImage[index][z].Pixel(x, y-1, 0);

			if (y<height-1)
				gtpixqd = &gtdirImage[index][z].Pixel(x, y+1, 0);

			if (z>0)
				gtpixqt = &gtdirImage[index][z-1].Pixel(x, y, 0);

			if (z>depth-1)
				gtpixqb = &gtdirImage[index][z+1].Pixel(x, y, 0);


			for (int d=0; d < prm->nD; d++)
			{

				if (fpvp->csgrad == 0)
				{

					int gr = 0;
					int gd = 0;
					int gl = 0;
					int gu = 0;
					int gt = 0;
					int gb = 0;

					// process left
					if (x>0)
					{
						gradl = &dirgrad[index][z].Pixel(x-1, y, 0);
						gl = gradl[0];

						d1 = d;
						d2 = gtpixql[0];

            if (prm->unknown == 1 && d2==prm->nD) {}
            else
            {
              reorder(&d1, &d2);
              totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gl] += exp((dCArray[0])[pixelnD]) * factor ;
            }
					}

					// process right
					if (x<width-1)
					{
						gradr = &dirgrad[index][z].Pixel(x, y, 0);
						gr = gradr[0];

						d1 = d;
						d2 = gtpixqr[0];
            if (prm->unknown == 1 && d2==prm->nD) {}
            else
            {
              reorder(&d1, &d2);
              totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gr] += exp((dCArray[0])[pixelnD]) * factor;
            }
					}


					// process up
					if (y>0)
					{
						gradu = &dirgrad[index][z].Pixel(x, y-1, 1);
						gu = gradu[0];

						d1 = d;
						d2 = gtpixqu[0];
            if (prm->unknown == 1 && d2==prm->nD) {}
            else
            {
              reorder(&d1, &d2);
              totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gu] += exp((dCArray[0])[pixelnD]) * factor;
            }
					}

					// process down
					if (y<height-1)
					{
						gradd = &dirgrad[index][z].Pixel(x, y, 1);
						gd = gradd[0];

						d1 = d;
						d2 = gtpixqd[0];
            if (prm->unknown == 1 && d2==prm->nD) {}
            else
            {
              reorder(&d1, &d2);
              totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gd] += exp((dCArray[0])[pixelnD]) * factor;
            }
					}

					// process top
					if (fpvp->gradVZ == 2 || fpvp->gradVZ == 3)
					if (z>0)
					{
						gradt = &dirgrad[index][z-1].Pixel(x, y, 2);
						gt = gradt[0];

						d1 = d;
						d2 = gtpixqt[0];
            if (prm->unknown == 1 && d2==prm->nD) {}
            else
            {
              reorder(&d1, &d2);
              if (fpvp->gradVZ == 1)
              {}
              else if (fpvp->gradVZ == 2 && fpvp->nZ > 1)
                totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gt] += exp((dCArray[0])[pixelnD]) * factor;
              else if (fpvp->gradVZ == 3 && fpvp->nZ > 1)
                totalModelDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nGZ + gt] += exp((dCArray[0])[pixelnD]) * factor;
            }
					}

					// process bottom
					if (fpvp->gradVZ == 2 || fpvp->gradVZ == 3)
					if (z<depth-1)
					{
						gradb = &dirgrad[index][z].Pixel(x, y, 2);
						gb = gradb[0];

						d1 = d;
						d2 = gtpixqb[0];
            if (prm->unknown == 1 && d2==prm->nD) {}
            else
            {
              reorder(&d1, &d2);
              if (fpvp->gradVZ == 1)
              {}
              else if (fpvp->gradVZ == 2 && fpvp->nZ > 1)
                totalModelDistV[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nG + gb] += exp((dCArray[0])[pixelnD]) * factor;
              else if (fpvp->gradVZ == 3 && fpvp->nZ > 1)
                totalModelDistZ[(d1*prm->nD + d2 - (d1*(d1+1))/2)*fpvp->nGZ + gb] += exp((dCArray[0])[pixelnD]) * factor;
            }
          }
				}
				else if (fpvp->csgrad == 1)
				{
					std::cout << " This needs to be done !" << std::endl;
					exit(1);
				}

				pixelnD ++;

			}
		}
	}

}



