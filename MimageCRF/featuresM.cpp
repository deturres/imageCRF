/*
 * featuresM.cpp
 *
 *  Created on: Mar 16, 2011
 *      Author: bhole
 */

#include "featuresM.h"

extern GlobalPairwiseParam globalP;


features::features(featureparams* featurepvp)
{
	  LocationMVclass=0;
	  appclass = 0;
	  AppMVclass=0;
	  intObj=0;
	  intNBclass = 0;

	  fpvp = featurepvp;

//	  dirgrad.erase(dirgrad.begin(), dirgrad.end());
//	  dircsgrad.erase(dircsgrad.begin(), dircsgrad.end());

}


features::~features()
{

	if (fpvp->klrpv.klr==1)
		delete ivm; // pointers in this point to static objects so no need to delete

	fpvp = 0; // we don't free since it is allocated statically

}



int features::readIntensityFiles(int nD)
{
	  if (fpvp->intensity==3)
	    readintClassParameter(intObj, nD);
}


int features::readLocationFiles(int nD, unsigned int testDir, int numP, int timeseq, std::vector<std::string> dirListNFP)
{

	// for location with mean and covariances
	if (fpvp->loc==3 || fpvp->loc==6 || fpvp->loc==8 || fpvp->loc==7 || fpvp->loc==9 || fpvp->klrpv.locKlr==2)
		readLocationMVParameters(nD);

	if (fpvp->loc==1)
	{ // location as cubes
	  std::cout << " Location as cubes has been disabled \n";
	  exit(1);
	}


	std::vector<CImage2> setLocDisp;

	locdirImage.erase(locdirImage.begin(), locdirImage.end());

	if (fpvp->loc==4 || fpvp->loc==5)
	{
		std::vector<string> dirList; // temporary variable: stores the subdirectory names (patient folder names)
		dirList.erase(dirList.begin(), dirList.end());
		readDirectories(fpvp->locdirname, dirList, locdirListFP, testDir);

		for(std::vector<string>::iterator b = locdirListFP.begin(); b < locdirListFP.end(); ++b)
		{
			cout << *b <<"\n";
			setLocDisp.erase(setLocDisp.begin(), setLocDisp.end());
			readImages(*b, setLocDisp);
			locdirImage.push_back(setLocDisp);
		}
	} else {

	  //for(std::vector<string>::iterator b = indirListFP.begin(); b < indirListFP.end(); ++b) {
	  for (int i=0; i<numP; i++)
	  {
	    locdirImage.push_back((std::vector<CImage2>) 0);
	  }
	}



	if (timeseq == 0 && (fpvp->loc != 0 || fpvp->klrpv.locKlr != 0))
	{
		read3DlocationSliceNumbers(numP, dirListNFP);

	}
}




void features::read3DlocationSliceNumbers(int numP, std::vector<std::string> dirListNFP)
{

    startSliceNo.resize(numP);
    endSliceNo.resize(numP);

    char fname[500];
    std::string locd = fpvp->imageSliceFile;
    sprintf(fname, "%s", locd.c_str());

    std::string lineS, foldname;
    int startno, endno;

    ifstream myfileS (fname);
    if (myfileS.is_open())
      {

        myfileS >> foldname ;

        while (! myfileS.eof() )
          {

            myfileS >> startno >> endno;

            int jj=0;
            for(std::vector<string>::iterator b = dirListNFP.begin(); b < dirListNFP.end(); ++b)
            {
              if (foldname.compare(*b)==0)
              {
                startSliceNo[jj] = startno;
                endSliceNo[jj] = endno;
              }
              jj++;
            }

            myfileS >> foldname ;
            if ( myfileS.eof() )
              break;

          }
        myfileS.close();
      }
    else {
      std::cout << "Unable to open parameter file for imageSliceInfo ";
      exit(1);
    }

}


// reading parameter files

void features::readLocationMVParameters(int nD)
{

  LocationMVclass = new LocationMV*[nD];

  // read parameter files
  //read bgnd, liver, rk, lk, gb, spleen === older stuff
  for (int kk=0; kk<nD; kk++) {
    std::string line;
    char cname[500], fname[500];
    int numC, dim;
    int check = 0;

    std::string locd = fpvp->locdirname;

    sprintf(fname, "%s/class%03d.txt", locd.c_str(), kk);

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

    // DEBUG_OUT3(verbose, debugfile, "cname = %s numC = %d dim = %d\n", cname, numC, dim);

    LocationMVclass[kk] = new LocationMV(numC, dim);
    LocationMVclass[kk]->setClassname(cname);
    LocationMVclass[kk]->readParameters(fname);
    //			LocationMVclass[kk]->printParameters();

  }

}






















int features::readHOGFiles(int nD, unsigned int testDir, int numP)
{


    if (fpvp->hog==1)
    {
    	std::vector<string> dirList; // temporary variable: stores the subdirectory names (patient folder names)
		dirList.erase(dirList.begin(), dirList.end());
		readDirectories(fpvp->hogdirname, dirList, hogdirListFP, testDir);
    }


    std::vector<CImage2> setHoGDisp;

    hogdirImage.erase(hogdirImage.begin(), hogdirImage.end());
    if (fpvp->hog==1)
    {
    	for(std::vector<string>::iterator b = hogdirListFP.begin(); b < hogdirListFP.end(); ++b)
    	{
    		cout << *b <<"\n";
    		setHoGDisp.erase(setHoGDisp.begin(), setHoGDisp.end());
    		readImages(*b, setHoGDisp);
    		hogdirImage.push_back(setHoGDisp);
    	}
    } else {
      for (int i=0; i< numP; i++) {
        hogdirImage.push_back((std::vector<CImage2>) 0);
      }
    }

    hogdirMVProb.erase(hogdirMVProb.begin(), hogdirMVProb.end());

    if (fpvp->hog==3)
    {

    	std::cout << " This option has been disabled due to memory requirements " <<std::endl;
    	exit(1);
/*
        hogdirMVProb.resize(numP);
        for (int mm=0; mm<numP; mm++)
          {
            int numSlices = gtdirImage[mm].size();
            hogdirMVProb[mm].resize(numSlices);
            for (int nn=0; nn<numSlices; nn++ )
              {
                hogdirMVProb[mm][nn].erase(hogdirMVProb[mm][nn].begin(), hogdirMVProb[mm][nn].end());
              }
          }
        // will need to change the working on prm.iopv.testDirIndexV here
        readHoGProb(indirImage, nD, hogdirMVProb, testDirIndexV);
*/

    } else
    {
      for (int i=0; i< numP; i++) {
        hogdirMVProb.push_back((std::vector <std::vector <matrixB<double> > >)0); // runtime bug in mac os x
      }

    }

}







int features::readMOTFiles(int nD, unsigned int testDir, int numP)
{


    if (fpvp->mot==1)
    {
    	std::vector<string> dirList; // temporary variable: stores the subdirectory names (patient folder names)
		dirList.erase(dirList.begin(), dirList.end());
		readDirectories(fpvp->motdirname, dirList, motdirListFP, testDir);
    }


    std::vector<CImage2> setMotDisp;

    motdirImage.erase(motdirImage.begin(), motdirImage.end());
    if (fpvp->mot==1)
    {
    	for(std::vector<string>::iterator b = motdirListFP.begin(); b < motdirListFP.end(); ++b)
    	{
    		cout << *b <<"\n";
    		setMotDisp.erase(setMotDisp.begin(), setMotDisp.end());
    		readImages(*b, setMotDisp);
    		motdirImage.push_back(setMotDisp);
    	}
    } else {
      for (int i=0; i< numP; i++) {
        motdirImage.push_back((std::vector<CImage2>) 0);
      }
    }

    motdirMVProb.erase(motdirMVProb.begin(), motdirMVProb.end());

    if (fpvp->mot==3)
    {

    	std::cout << " This option has been disabled due to memory requirements " <<std::endl;
    	exit(1);
/*
        hogdirMVProb.resize(numP);
        for (int mm=0; mm<numP; mm++)
          {
            int numSlices = gtdirImage[mm].size();
            hogdirMVProb[mm].resize(numSlices);
            for (int nn=0; nn<numSlices; nn++ )
              {
                hogdirMVProb[mm][nn].erase(hogdirMVProb[mm][nn].begin(), hogdirMVProb[mm][nn].end());
              }
          }
        // will need to change the working on prm.iopv.testDirIndexV here
        readHoGProb(indirImage, nD, hogdirMVProb, testDirIndexV);
*/

    } else
    {
      for (int i=0; i< numP; i++)
        motdirMVProb.push_back((std::vector <std::vector <matrixB<double> > >)0);
    }

}












int features::readAppearanceFiles(int timeseq, int nD)
{
	  if (fpvp->app==1 && timeseq==0)
	    readApprearanceParametersCPC(appclass, nD, fpvp->appfilename, 1);
	  else if (fpvp->app==1 && timeseq==1)
	  {
		  // assume uv stores patches as patchsize for u and then patch size for v
		  readApprearanceParametersCPC(appclass, nD, fpvp->appfilename, 2);
	  }

	  if (fpvp->app==2)  // we must remember that there aren't clusters per class.
	    readApprearanceParameters(appclass, nD, fpvp->appfilename);


	  if (fpvp->app==3)
	    readAppearanceMVParameter(AppMVclass, nD, fpvp->appfilename);

}

int features::readKlrFiles(int nD)
{
	if (fpvp->klrpv.klr==1)  // read main file,  weight parameter file and points or indices
	{
		// we do not read a klr parameter file since the main parameter file takes care of it
		// in the future the most variable arguments could be passed using the command line

		printklrParameters();

		// readWeightParams
		if(!load_one_data(fpvp->klrpv.wfile.c_str(), wparam)){
		std::cout << "ivm file (error): Bad data file filename : "<< fpvp->klrpv.wfile << endl;
		exit(1);
		}

		fpvp->klrpv.nW = wparam.size1() * wparam.size2();

		// readWeightPoints and neglect 1st column
		if(!load_one_data(fpvp->klrpv.xfile.c_str(), Xtrain, 1)){
		std::cout << "ivm file (error): Bad data file filename : "<< fpvp->klrpv.xfile << endl;
		exit(1);
		}

		// readMinMaxMatrix - this is needed for the normalization of the test data point
		// and neglect 1st column
		if(!load_one_data(fpvp->klrpv.mmfile.c_str(), minmaxM, 1)){
		std::cout << "ivm file (error): Bad data file filename : "<< fpvp->klrpv.mmfile << endl;
		exit(1);
		}


		fpvp->klrpv.nW = 0;
	    if (fpvp->klrpv.klr==1) {
	    	fpvp->klrpv.nW = wparam.size1() * wparam.size2();
	    }

	    // create ivm
	    // this statement is throwing errors
		// ivm = new IVM::IVM(&Xtrain, &wparam, fpvp->klrpv.ker, fpvp->klrpv.lambda, nD, fpvp->klrpv.p1);
	    ivm = 0;

		// read desc directory listing

		std::vector<string> dirList; // temporary variable: stores the subdirectory names (patient folder names)

		if (fpvp->klrpv.hogKlr==2)
		{
			dirList.erase(dirList.begin(), dirList.end());
			readDirectories(fpvp->klrpv.hogklrdescdirname, dirList, hogklrdescdirListFP, 0);
		}

		if (fpvp->klrpv.appKlr==2)
		{
			dirList.erase(dirList.begin(), dirList.end());
			readDirectories(fpvp->klrpv.appklrdescdirname, dirList, appklrdescdirListFP, 0);
		}


	}
}


void features::printklrParameters()
{

	std::cout << " filename : " << fpvp->klrpv.klrfName << std::endl;
	std::cout << " hog dims : " << fpvp->klrpv.hogdimklr << std::endl;
	std::cout << " mot dims : " << fpvp->klrpv.motdimklr << std::endl;
	std::cout << " app dims : " << fpvp->klrpv.appdimklr << std::endl;
	std::cout << " loc dims : " << fpvp->klrpv.locdimklr << std::endl;
	std::cout << " bias  : " << fpvp->klrpv.biasklr << std::endl;
	std::cout << " lambda : " << fpvp->klrpv.lambda << std::endl;
	std::cout << " xfile : " << fpvp->klrpv.xfile << std::endl;
	std::cout << " wfile : " << fpvp->klrpv.wfile << std::endl;
	std::cout << " mmfile : " << fpvp->klrpv.mmfile << std::endl;
	std::cout << " kernel : " << fpvp->klrpv.kernel << std::endl;
	std::cout << " p1 : " << fpvp->klrpv.p1 << std::endl;
	std::cout << " ker : " << fpvp->klrpv.ker << std::endl;

}



int features::readFeatureFilesDirs(int timeseq, int nD, int testDir, int numP, std::vector<std::string> dirListNFP)
{

	// intensity
	readIntensityFiles(nD);

	// location
	readLocationFiles(nD, testDir, numP, timeseq, dirListNFP);

	// hog
	readHOGFiles(nD, testDir, numP);

	// mot
	readMOTFiles(nD, testDir, numP);


	// appearance
	readAppearanceFiles(timeseq, nD);

	// klr
	readKlrFiles(nD);

}



int features::deleteFeatureAllocations(int nD, int timeseq)
{
	// location
	deleteLocationMVclass(nD);
	deAllocateLocationFiles(timeseq);

	// appearance
	deleteAppMVclass(nD);
	deleteappclass(nD);

	//intensity
	deleteintObj(nD);

}


int features::deAllocatePairwiseAllocations()
{

	// if quantized it is not malloced but as static objects
	// if smooth then we have used opencv and need to free images

	if (fpvp->csgrad == 1)
	{
		int numP = dircsgrad.size();

		// deallocating space
		for (int i=0; i< numP; i++)
		{
			int depth = dircsgrad[i].size();

			for (int j=0; j<depth; j++)
			{
				cvReleaseImage(&dircsgrad[i][j]);
			}

		}

	}
}




int features::computeIntensityNB(std::vector <std::vector <CImage2> > *inp, std::vector <std::vector <CByteImage> > *gtr, std::vector<int> testDirIndexV, int nD)
{

	// this will run only for timeseq==0
	if (fpvp->intensity==6)
	{
		intNBclass = new intensityNB(fpvp->nbinsintNB, nD);
		intNBclass->setSmoothFactor(0.0001);
		intNBclass->readProbValues(inp, gtr, testDirIndexV, fpvp->rangeI);
		intNBclass->computePofXcY();
	}
}


int features::computeIntensityNB(std::vector <std::vector <CByteImage> > *inp, std::vector <std::vector <CByteImage> > *gtr, std::vector<int> testDirIndexV, int nD)
{

	// this will run only for timeseq==0
	if (fpvp->intensity==6)
	{
		std::cout << " This has not been implemented yet !!!" <<std::endl;
		exit(1);
	}
}




void features::computeAppearanceVectors(std::vector <std::vector <CImage2> > * im1, int nD)
{
  // checkForFloatCosts();

  int nA = fpvp->nA;

  int nT;
  if (fpvp->app==1)
    nT = nA/nD; //this gives centers per class
  else
    nT = appclass[0]->getNTextons();  // since nD = 1 when we pass app==2

  if (nA<=0)
    return;



  for (unsigned int j = 0; j < (*im1).size(); ++j)
  {

    int depth = (*im1)[j].size();

    CShape sh = (*im1)[j][0].Shape();
    int width = sh.width, height = sh.height, nB = sh.nBands;
    int tempglobalNpixels = depth * width * height;
    int nColors = __min(3, nB);

    for(unsigned int i = 0; i < (*im1)[j].size(); ++i)
    {

    	// DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);

    	sh.nBands = 1;

		for (int ss=0; ss<nD; ss++)
			appdirImage[j][i][ss].ReAllocate(sh);

		int dsiIndex = 0;

		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				for (int d = 0; d < nD; d++) {

					uchar *appIm = &appdirImage[j][i][d].Pixel(x, y, 0);

					// assume the patch size is the same for all classes
					int patchSize = appclass[d]->getPatchSize();
					int nTextons = appclass[d]->getNTextons();

					if (nTextons!=nT)
					{
						printf(" Mismatch in number of textons per class \n");
						exit(1);
					}

					float** textons = appclass[d]->getTextons();

					if (x>(int)patchSize/2 && x<width-patchSize/2 && y>(int)patchSize/2 && y<height-patchSize/2)
					{
						// select best texton center
						int bestT=-1;
						double bestDist=9e99;

						for (int mm=0; mm<nTextons; mm++)
						{
							double distT=0;
							int cnter=0;
							for (int qq=-(int)patchSize/2 + y; qq<=(int)patchSize/2 + y; qq++)
								for (int pp=-(int)patchSize/2 + x; pp<=(int)patchSize/2 + x; pp++)
								{
									distT += (textons[mm][cnter]-(float)(*im1)[j][i].Pixel(pp, qq, 0)) * (textons[mm][cnter]-(float)(*im1)[j][i].Pixel(pp, qq, 0));
									cnter++;
								}

							if (bestDist>distT)
							{
								bestDist = distT;
								bestT = mm;
							}

						}

						if (bestT==-1)
						{
							printf(" some error in appearance model \n");
							exit(1);
						}

						*appIm = (unsigned char)bestT;
					}
				}
			}
		}
	}
  }
}






void features::computeAppearanceProb(std::vector <std::vector <CImage2> > * im1, int nD)
{

  for (unsigned int j = 0; j < (*im1).size(); ++j)
  {

    int depth = (*im1)[j].size();

    CShape sh = (*im1)[j][0].Shape();
    int width = sh.width, height = sh.height, nB = sh.nBands;

    matrixB<double> improb(height, width);

    for(unsigned int i = 0; i < (*im1)[j].size(); ++i)
    {

      // DEBUG_OUT4(verbose, debugfile, "Appearance: Image [%d][%d] \n", j, i);

    	for (int d = 0; d < nD; d++)
    	{
    		for (int y = 0; y < height; y++)
    		{
    			for (int x = 0; x < width; x++)
    			{

    				float appcost;
    				float *pii;
    				float **mu;
    				float ***SigmaInv;
    				float *normConst;
    				int dim;

    				dim = AppMVclass[d]->getDimension();
    				pii = AppMVclass[d]->getPi();
    				mu = AppMVclass[d]->getMu();
    				SigmaInv = AppMVclass[d]->getSigmaInv();
    				normConst = AppMVclass[d]->getNormConst();


            //					  AppMVclass[d]->printParameters();

    				double appval = 0.0;

    				int patchSize = sqrt(dim);

    				if (x>(int)patchSize/2 && x<width-patchSize/2 && y>(int)patchSize/2 && y<height-patchSize/2)
    				{
    					for (int c=0; c<AppMVclass[d]->getNClusters(); c++)
    					{
    						int nnn = 0;
    						double tempval=0;
    						for (int xn=x-(int)patchSize/2; xn<=x+(int)patchSize/2; xn++)
    						{
    							for (int yn=y-(int)patchSize/2; yn<=y+(int)patchSize/2; yn++)
    							{
    								unsigned short pixt = (*im1)[j][i].Pixel(xn, yn, 0);

    								double xs = (double)pixt-mu[c][nnn];

    								tempval += xs*xs*SigmaInv[c][nnn][nnn];

    							}
    						}
    						appval += exp(log(pii[c])+(-tempval/2)-log(normConst[c]));
    					}

    					if (appval==0.0)
    						appval = pow(10.0, -300.0);  // tiny value greater than zero

    					appcost = -log(appval);

    					(improb)(y, x) = appcost;

    				} else {
    					(improb)(y,x) = -log(1/nD); //uniform distribution
    				}
    			}
    		} //finished filling matrix(y,x)

    		appdirMVProb[j][i].push_back(improb);

      }
    }
  }

}






// this uses Yuv (for now only uv) and shrinks image to half for computation reasons
// remember k-means and/or EM are trained using this half images to find the centers too
void features::computeAppearanceVectors(std::vector <std::vector <CByteImage> > *im1, int nD)
{
  //checkForFloatCosts();

  int nA = fpvp->nA;

  int nT;
  if (fpvp->app==1)
    nT = nA/nD; //this gives centers per class
  else
    nT = appclass[0]->getNTextons();  // since nD = 1 when we pass app==2

  if (nA<=0)
    return;

  for (unsigned int j = 0; j < (*im1).size(); ++j)
  {

	  int depth = (*im1)[j].size();

	  CShape sh = (*im1)[j][0].Shape();
	  int width = sh.width, height = sh.height, nB = sh.nBands;
	  int tempglobalNpixels = depth * width * height;
	  int nColors = __min(3, nB);

	  for(unsigned int i = 0; i < (*im1)[j].size(); ++i){

		  // DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);

		  sh.nBands = 1;

		  for (int ss=0; ss<nD; ss++)
			  appdirImage[j][i][ss].ReAllocate(sh);

      // create image buffer of half size
      std::vector < std::vector <std::vector <unsigned char> > > halfImage;
      halfImage.resize(height/2);
      for (int kk=0; kk< height/2; kk++)
      {
    	  halfImage[kk].resize(width/2);
    	  for (int mm=0; mm<width/2; mm++)
    	  {
    		  halfImage[kk][mm].resize(3);
    	  }
      }

      // copy image to make half and convert to Yuv
      for (int kk=0; kk< height/2; kk++)
      {
    	  for (int mm=0; mm<width/2; mm++)
    	  {
    		  double tempR=0, tempG=0, tempB=0;
    		  int rgb[3];
    		  for (int nn=kk*2; nn<=kk*2+1; nn++)
    		  {
    			  for (int ss=mm*2; ss<=mm*2+1; ss++)
    			  {
    				  if (nn>=height || ss>=width)
    					  continue;

					  tempB += (float)(*im1)[j][i].Pixel(ss, nn, 0);
					  tempG += (float)(*im1)[j][i].Pixel(ss, nn, 1);
					  tempR += (float)(*im1)[j][i].Pixel(ss, nn, 2);
    			  }
    		  }

    		  tempR = tempR/4;
    		  tempG = tempG/4;
    		  tempB = tempB/4;

    		  rgb[0]=round(tempB);
    		  rgb[1]=round(tempG);
    		  rgb[2]=round(tempR);
    		  int Yuv[3];
    		  rgb2Yuv(rgb, Yuv);

    		  int r=rgb[2], g=rgb[1], b=rgb[0];

    		  halfImage[kk][mm][0] = Yuv[0];  //v
    		  halfImage[kk][mm][1] = Yuv[1];  //u
    		  halfImage[kk][mm][2] = Yuv[2];  //Y

    	  }
      }


      int dsiIndex = 0;
      for (int y = 0; y < height/2; y++)
      {
    	  for (int x = 0; x < width/2; x++)
    	  {
    		  for (int d = 0; d < nD; d++)
    		  {

    			  // assume the patch size is the same for all classes
    			  int patchSize = appclass[d]->getPatchSize();
    			  int nTextons = appclass[d]->getNTextons();

    			  if (nTextons!=nT)
    			  {
    				  printf(" Mismatch in number of textons per class \n");
    				  exit(1);
    			  }

    			  float** textons = appclass[d]->getTextons();
    			  // textons are arranged all u then v for uv for a given patch
    			  // or all Y and then all u and then v for a given patch for Yuv

    			  // here using only uv for now

    			  if (x>(int)patchSize/2 && x<width/2-patchSize/2 && y>(int)patchSize/2 && y<height/2-patchSize/2)
    			  {
    				  // select best texton center
    				  int bestT=-1;
    				  double bestDist=9e99;

    				  for (int mm=0; mm<nTextons; mm++)
    				  {
    					  double distT=0;
    					  int cnter=0;

    					  // u => note texton is stored as [u v]
    					  for (int qq=-(int)patchSize/2 + y; qq<=(int)patchSize/2 + y; qq++)
    						  for (int pp=-(int)patchSize/2 + x; pp<=(int)patchSize/2 + x; pp++)
    						  {
    							  distT += (textons[mm][cnter]-halfImage[qq][pp][1]) * (textons[mm][cnter]-halfImage[qq][pp][1]);
    							  cnter++;
    						  }

    					  // v
    					  for (int qq=-(int)patchSize/2 + y; qq<=(int)patchSize/2 + y; qq++)
    						  for (int pp=-(int)patchSize/2 + x; pp<=(int)patchSize/2 + x; pp++)
    						  {
    							  distT += (textons[mm][cnter]-halfImage[qq][pp][0]) * (textons[mm][cnter]-halfImage[qq][pp][0]);
    							  cnter++;
    						  }

    					  if (bestDist>distT)
    					  {
    						  bestDist = distT;
    						  bestT = mm;
    					  }

    				  }

    				  if (bestT==-1)
    				  {
    					  printf(" some error in appearance model \n");
    					  exit(1);
    				  }

    				  for (int yy=2*y; yy<=2*y+1; yy++)
    				  {
    					  for(int xx=2*x;xx<=2*x+1; xx++)
    					  {
    						  if (xx>=width || yy>=height)
    							  continue;
    						  uchar *appIm = &appdirImage[j][i][d].Pixel(xx, yy, 0);
    						  *appIm = (unsigned char)bestT;
    					  }
    				  }

    			  }
    		  }
    	  }
      }
    }
  }
}








int features::computeAppearance(std::vector <std::vector <CImage2> > *indirImage, std::vector <std::vector <CByteImage> > *gtdirImage, int nD)
{

	appdirImage.erase(appdirImage.begin(), appdirImage.end());
	appdirMVProb.erase(appdirMVProb.begin(), appdirMVProb.end());

	int numPatients = (*gtdirImage).size();

	if (fpvp->app==1 || fpvp->app==2)
	{
		CByteImage tempImage;

		appdirImage.resize(numPatients);
		for (int mm=0; mm<numPatients; mm++)
		{
			int numSlices = (*gtdirImage)[mm].size();
			appdirImage[mm].resize(numSlices);
			for (int nn=0; nn<numSlices; nn++ )
			{
				appdirImage[mm][nn].resize(nD);
			}
		}

		if (fpvp->app==1)
		  computeAppearanceVectors(indirImage, nD);
		else if (fpvp->app==2)
		  computeAppearanceVectors(indirImage, 1);

	} else {
	  for (int i=0; i< numPatients; i++)
		appdirImage.push_back((std::vector <std::vector <CByteImage> >)0);
	}


	if (fpvp->app==3)
	{
		appdirMVProb.resize(numPatients);
		for (int mm=0; mm<numPatients; mm++)
		{
			int numSlices = (*gtdirImage)[mm].size();
			appdirMVProb[mm].resize(numSlices);
			for (int nn=0; nn<numSlices; nn++ )
			{
				appdirMVProb[mm][nn].erase(appdirMVProb[mm][nn].begin(), appdirMVProb[mm][nn].end());
			}
		}
		computeAppearanceProb(indirImage, nD);

	  } else {
		  for (int i=0; i< numPatients; i++)
			appdirMVProb.push_back((std::vector <std::vector <matrixB<double> > >)0);
	  }


}





int features::computeAppearance(std::vector <std::vector <CByteImage> > *indirImaget, std::vector <std::vector <CByteImage> > *gtdirImage, int nD)
{

	appdirImage.erase(appdirImage.begin(), appdirImage.end());
	appdirMVProb.erase(appdirMVProb.begin(), appdirMVProb.end());

	int numPatients = (*gtdirImage).size();

	if (fpvp->app==1 || fpvp->app==2)
	{
		CByteImage tempImage;

		appdirImage.resize(numPatients);
		for (int mm=0; mm<numPatients; mm++)
		{
			int numSlices = (*gtdirImage)[mm].size();
			appdirImage[mm].resize(numSlices);
			for (int nn=0; nn<numSlices; nn++ )
			{
				appdirImage[mm][nn].resize(nD);
			}
		}

		// these are Yuv specific
		if (fpvp->app==1)
		  computeAppearanceVectors(indirImaget, nD);
		else if (fpvp->app==2)
		  computeAppearanceVectors(indirImaget, 1);

	} else {
	  for (int i=0; i< numPatients; i++)
		appdirImage.push_back((std::vector <std::vector <CByteImage> >)0);
	}


	if (fpvp->app==3)
	{
		std::cout << " This has not been implemented yet !!!" <<std::endl;
		exit(1);

	} else {
		  for (int i=0; i< numPatients; i++)
			appdirMVProb.push_back((std::vector <std::vector <matrixB<double> > >)0);
	}


}




// grayscale (for medical images)
void features::computeCSGradients(std::vector <std::vector <CImage2> > *inp)
{
  // simple2
  double simple2[2] = {-1, +1};

  int numP = (*inp).size();

  std::vector<double> expectedGradDiff;
  expectedGradDiff.resize(numP);

  std::vector<double> expectedGradDiffZ;
  expectedGradDiffZ.resize(numP);


  // allocating space first
  for (int i=0; i< numP; i++)
  {

	  int depth = (*inp)[i].size();
	  CShape sh = (*inp)[i][0].Shape();
	  int width = sh.width, height = sh.height, nB = sh.nBands;
	  nB = 3;  // band 0: x gradient, band 1: y gradient, band 2: z gradient

	  std::vector <IplImage*> setImgrad;
	  setImgrad.erase(setImgrad.begin(), setImgrad.end());

	  for (int j=0; j<depth; j++)
	  {
		  IplImage* img = cvCreateImage( cvSize(width,height), IPL_DEPTH_32F, nB );
		  setImgrad.push_back(img);
	  }

	  dircsgrad.push_back(setImgrad);


	  // computing the expected gradient difference
	  expectedGradDiff[i] = 0.0;
	  expectedGradDiffZ[i] = 0.0;

	  int edges = 0;
	  int edgesZ = 0;

	  for (int z = 0; z < depth; z++)
	  {
		  for (int y = 0; y < height; y++)
		  {
			  for (int x = 0; x < width; x++)
			  {
				  int pix1 = (*inp)[i][z].Pixel(x, y, 0);
				  int pix2 = 0;
				  int pix3 = 0;
				  int pix4 = 0;

				  double dx=0.0, dy=0.0, dz=0.0;

				  if (x==width-1)
					  dx = 0;
				  else {
					  pix2 = (*inp)[i][z].Pixel(x+1, y, 0);
					  dx = fabs(simple2[0]*pix2 + simple2[1]*pix1);
					  edges++;
				  }

				  if (y==height-1)
					  dy = 0;
				  else {
					  pix3 = (*inp)[i][z].Pixel(x, y+1, 0);
					  dy = fabs(simple2[0]*pix3 + simple2[1]*pix1);
					  edges++;
				  }

				  if (z==depth-1)
					  dz = 0;
				  else {
					  pix4 = (*inp)[i][z+1].Pixel(x, y, 0);
					  dz = fabs(simple2[0]*pix4 + simple2[1]*pix1);
	          if (fpvp->gradVZ == 2)
	            edges++;
	          else if (fpvp->gradVZ == 3)
	            edgesZ++;
				  }

				  if (fpvp->gradVZ == 1) // (x-y)
				  {
					  expectedGradDiff[i] += dx*dx + dy*dy;
				  }
				  else if (fpvp->gradVZ == 2) // (x-y-z)
				  {
					  expectedGradDiff[i] += dx*dx + dy*dy + dz*dz;
				  }
				  else if (fpvp->gradVZ == 3) // (x-y and z)
				  {
					  expectedGradDiff[i] += dx*dx + dy*dy;
					  expectedGradDiffZ[i] += dz*dz;
				  }

			  }
		  }
	  }

	   if (edges != 0)
	      expectedGradDiff[i] /= edges;
	   if (edgesZ != 0 && fpvp->gradVZ == 3)
	      expectedGradDiffZ[i] /= edgesZ;

  }


  // then computing gradients
  for (int i=0; i< numP; i++)
  {
	  int depth = (*inp)[i].size();
	  CShape sh = (*inp)[i][0].Shape();
	  int width = sh.width, height = sh.height, nB = sh.nBands;
	  nB = 3;  // band 0: x gradient, band 1: y gradient, band 2: z gradient


	  for (int z = 0; z < depth; z++)
	  {
		  for (int y = 0; y < height; y++)
		  {
			  for (int x = 0; x < width; x++)
			  {

				  int pix1 = (*inp)[i][z].Pixel(x, y, 0);
				  int pix2 = 0;
				  int pix3 = 0;
				  int pix4 = 0;

				  double dx=0.0, dy=0.0, dz=0.0;

				  if (x==width-1)
					  dx = 0;
				  else {
					  pix2 = (*inp)[i][z].Pixel(x+1, y, 0);
					  dx = fabs(simple2[0]*pix2 + simple2[1]*pix1);
				  }

				  if (y==height-1)
					  dy = 0;
				  else {
					  pix3 = (*inp)[i][z].Pixel(x, y+1, 0);
					  dy = fabs(simple2[0]*pix3 + simple2[1]*pix1);
				  }

				  if (z==depth-1)
					  dz = 0;
				  else {
					  pix4 = (*inp)[i][z+1].Pixel(x, y, 0);
					  dz = fabs(simple2[0]*pix4 + simple2[1]*pix1);
				  }

				  if (fpvp->gradVZ == 1)
				  {
					  CvScalar s;
					  s=cvGet2D(dircsgrad[i][z],y,x); // get the (i,j) pixel value
					  s.val[0] = log(fpvp->thetaV[0]) - dx*dx/(4*expectedGradDiff[i]);
					  s.val[1] = log(fpvp->thetaV[0]) - dy*dy/(4*expectedGradDiff[i]);
					  s.val[2] = 0;
					  cvSet2D(dircsgrad[i][z],y,x,s); // set the (i,j) pixel value
				  }

				  if (fpvp->gradVZ == 2)
				  {
					  CvScalar s;
					  s=cvGet2D(dircsgrad[i][z],y,x); // get the (i,j) pixel value
					  s.val[0] = log(fpvp->thetaV[0]) - dx*dx/(4*expectedGradDiff[i]);
					  s.val[1] = log(fpvp->thetaV[0]) - dy*dy/(4*expectedGradDiff[i]);
					  s.val[2] = log(fpvp->thetaV[0]) - dz*dz/(4*expectedGradDiff[i]);
					  cvSet2D(dircsgrad[i][z],y,x,s); // set the (i,j) pixel value
				  }

				  if (fpvp->gradVZ == 3)
				  {
					  CvScalar s;
					  s=cvGet2D(dircsgrad[i][z],y,x); // get the (i,j) pixel value
					  s.val[0] = log(fpvp->thetaV[0]) - dx*dx/(4*expectedGradDiff[i]);
					  s.val[1] = log(fpvp->thetaV[0]) - dy*dy/(4*expectedGradDiff[i]);
					  s.val[2] = log(fpvp->thetaZ[0]) - dz*dz/(4*expectedGradDiffZ[i]);
					  cvSet2D(dircsgrad[i][z],y,x,s); // set the (i,j) pixel value
				  }

			  }
		  }
	  }

  }


}













// grayscale (for medical images)
void features::computeQuantizedGradients(std::vector <std::vector <CImage2> > *inp)
  /*
    use nK = gradThresh.size() + 1 levels of quantization.
    quantization is done as follows:
    compute absolute gradient (of color/grayscale bands) in x and y (and z)
    given gradient g, assign largest k such that gradThresh[k] <= g
    the last entry (gradThresh[nK-1]) is ignored and is assumed to be infinity
    in im1grad, use band 0 for x gradient, and band 1 for y gradients and band 2 for z gradients
  */
{

  // simple2
  double simple2[2] = {-1, +1};

  int k;

  int numP = (*inp).size();

  // allocating space first
  std::vector<CByteImage> setImgrad;

  for (int i=0; i< numP; i++)
  {

	  setImgrad.erase(setImgrad.begin(), setImgrad.end());

	  int depth = (*inp)[i].size();
	  setImgrad.resize(depth);
	  dirgrad.push_back(setImgrad);

  }


  // then computing gradients
  for (int i=0; i< numP; i++)
  {

	  int depth = (*inp)[i].size();
	  CShape sh = (*inp)[i][0].Shape();
	  int width = sh.width, height = sh.height, nB = sh.nBands;
	  sh.nBands = 3;  // band 0: x gradient, band 1: y gradient, band 2: z gradient

	  for (int z = 0; z<depth; z++) {
		  dirgrad[i][z].ReAllocate(sh);
	  }

	  for (int z = 0; z < depth; z++) {
		  for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {

			  uchar *grad   = &dirgrad[i][z].Pixel(x, y, 0);

			  int pix1 = (*inp)[i][z].Pixel(x, y, 0);
			  int pix2 = 0;
			  int pix3 = 0;
			  int pix4 = 0;

			  double dx=0.0, dy=0.0, dz=0.0;

			  if (x==width-1)
				  dx = 0;
			  else {
				  pix2 = (*inp)[i][z].Pixel(x+1, y, 0);
				  dx = fabs(simple2[0]*pix2 + simple2[1]*pix1);
			  }

			  if (y==height-1)
				  dy = 0;
			  else {
				  pix3 = (*inp)[i][z].Pixel(x, y+1, 0);
				  dy = fabs(simple2[0]*pix3 + simple2[1]*pix1);
			  }

			  if (z==depth-1)
				  dz = 0;
			  else {
				  pix4 = (*inp)[i][z+1].Pixel(x, y, 0);
				  dz = fabs(simple2[0]*pix4 + simple2[1]*pix1);
			  }


			  for (k = 0; fpvp->gradThreshVec[k] <= dx; k++)
				; // find lowest bin for x gradient
			  assert(k < fpvp->nG);
			  grad[0] = k;

			  for (k = 0; fpvp->gradThreshVec[k] <= dy; k++)
				; // find lowest bin for y gradient
			  assert(k < fpvp->nG);
			  grad[1] = k;

			  if (fpvp->gradVZ == 1) // (x-y)
			  {

			  } else if (fpvp->gradVZ == 2) // (x-y-z)
			  {
				  for (k = 0; fpvp->gradThreshVec[k] <= dz; k++)
					; // find lowest bin for z gradient
				  assert(k < fpvp->nG);
				  grad[2] = k;

			  } else if (fpvp->gradVZ == 3) // (x-y and z)
			  {
				  for (k = 0; fpvp->gradThreshVecZ[k] <= dz; k++)
					; // find lowest bin for z gradient
				  assert(k < fpvp->nGZ);
				  grad[2] = k;
			  }

			}
		  }
	  }

  }
}










// RGB converted to grayscale (x-y or x-y-z or (x-y and z) )
void features::computeCSGradients(std::vector <std::vector <CByteImage> > *inp)
{

  // simple2
  double simple2[2] = {-1, +1};

  int numP = (*inp).size();

  std::vector<double> expectedGradDiff;
  expectedGradDiff.resize(numP);

  std::vector<double> expectedGradDiffZ;
  expectedGradDiffZ.resize(numP);


  // allocating space first
  for (int i=0; i< numP; i++)
  {

	  int depth = (*inp)[i].size();
	  CShape sh = (*inp)[i][0].Shape();
	  int width = sh.width, height = sh.height, nB = sh.nBands;
	  nB = 3;  // band 0: x gradient, band 1: y gradient, band 2: z gradient

	  std::vector <IplImage*> setImgrad;
	  setImgrad.erase(setImgrad.begin(), setImgrad.end());

	  for (int j=0; j<depth; j++)
	  {
		  IplImage* img = cvCreateImage( cvSize(width,height), IPL_DEPTH_32F, nB );
		  setImgrad.push_back(img);
	  }

	  dircsgrad.push_back(setImgrad);


	  // computing the expected gradient difference
	  expectedGradDiff[i] = 0.0;
	  expectedGradDiffZ[i] = 0.0;

	  int edges = 0;
	  int edgesZ = 0;

	  for (int z = 0; z < depth; z++)
	  {
		  for (int y = 0; y < height; y++)
		  {
			  for (int x = 0; x < width; x++)
			  {
				  unsigned char *pix1 = &(*inp)[i][z].Pixel(x, y, 0);
				  unsigned char *pix2 = 0;
				  unsigned char *pix3 = 0;
				  unsigned char *pix4 = 0;

				  double dx=0.0, dy=0.0, dz=0.0;
				  double pix1g, pix2g, pix3g, pix4g;

				  pix1g = 0.11*pix1[0] + 0.59*pix1[1] + 0.3*pix1[2];

				  if (x==width-1)
					  dx = 0;
				  else {
					  pix2 = &(*inp)[i][z].Pixel(x+1, y, 0);
					  pix2g = 0.11*pix2[0] + 0.59*pix2[1] + 0.3*pix2[2];
					  dx = fabs(simple2[0]*pix2g + simple2[1]*pix1g);
					  edges++;
				  }

				  if (y==height-1)
					  dy = 0;
				  else {
					  pix3 = &(*inp)[i][z].Pixel(x, y+1, 0);
					  pix3g = 0.11*pix3[0] + 0.59*pix3[1] + 0.3*pix3[2];
					  dy = fabs(simple2[0]*pix3g + simple2[1]*pix1g);
					  edges++;
				  }

				  if (z==depth-1)
					  dz = 0;
				  else {
					  pix4 = &(*inp)[i][z+1].Pixel(x, y, 0);
					  pix4g = 0.11*pix4[0] + 0.59*pix4[1] + 0.3*pix4[2];
					  dz = fabs(simple2[0]*pix4g + simple2[1]*pix1g);
					  if (fpvp->gradVZ == 2)
					    edges++;
					  else if (fpvp->gradVZ == 3)
					    edgesZ++;
				  }

				  if (fpvp->gradVZ == 1) // (x-y)
				  {
					  expectedGradDiff[i] += dx*dx + dy*dy;
				  }
				  else if (fpvp->gradVZ == 2) // (x-y-z)
				  {
					  expectedGradDiff[i] += dx*dx + dy*dy + dz*dz;
				  }
				  else if (fpvp->gradVZ == 3) // (x-y and z)
				  {
					  expectedGradDiff[i] += dx*dx + dy*dy;
					  expectedGradDiffZ[i] += dz*dz;
				  }

			  }
		  }
	  }

	  if (edges != 0)
	    expectedGradDiff[i] /= edges;
	  if (edgesZ != 0 && fpvp->gradVZ == 3)
	    expectedGradDiffZ[i] /= edgesZ;

  } // i


  // then computing gradients
  for (int i=0; i< numP; i++)
  {
	  int depth = (*inp)[i].size();
	  CShape sh = (*inp)[i][0].Shape();
	  int width = sh.width, height = sh.height, nB = sh.nBands;
	  nB = 3;  // band 0: x gradient, band 1: y gradient, band 2: z gradient

	  assert(expectedGradDiff[i]!=0);
	  if (fpvp->gradVZ == 3 && depth > 1)
	    assert(expectedGradDiffZ[i]!=0);

	  assert(fpvp->thetaV[0] > 0);
	  if (fpvp->gradVZ == 3)
	    assert(fpvp->thetaZ[0] > 0);

	  for (int z = 0; z < depth; z++)
	  {
		  for (int y = 0; y < height; y++)
		  {
			  for (int x = 0; x < width; x++)
			  {

				  unsigned char *pix1 = &(*inp)[i][z].Pixel(x, y, 0);
				  unsigned char *pix2 = 0;
				  unsigned char *pix3 = 0;
				  unsigned char *pix4 = 0;

				  double dx=0.0, dy=0.0, dz=0.0;
				  double pix1g, pix2g, pix3g, pix4g;

				  pix1g = 0.11*pix1[0] + 0.59*pix1[1] + 0.3*pix1[2];

				  if (x==width-1)
					  dx = 0;
				  else {
					  pix2 = &(*inp)[i][z].Pixel(x+1, y, 0);
					  pix2g = 0.11*pix2[0] + 0.59*pix2[1] + 0.3*pix2[2];
					  dx = fabs(simple2[0]*pix2g + simple2[1]*pix1g);
				  }

				  if (y==height-1)
					  dy = 0;
				  else {
					  pix3 = &(*inp)[i][z].Pixel(x, y+1, 0);
					  pix3g = 0.11*pix3[0] + 0.59*pix3[1] + 0.3*pix3[2];
					  dy = fabs(simple2[0]*pix3g + simple2[1]*pix1g);
				  }

				  if (z==depth-1)
					  dz = 0;
				  else {
					  pix4 = &(*inp)[i][z+1].Pixel(x, y, 0);
					  pix4g = 0.11*pix4[0] + 0.59*pix4[1] + 0.3*pix4[2];
					  dz = fabs(simple2[0]*pix4g + simple2[1]*pix1g);
				  }

				  if (fpvp->gradVZ == 1)
				  {
					  CvScalar s;
					  s=cvGet2D(dircsgrad[i][z],y,x); // get the (i,j) pixel value
					  // s.val[0] = log(fpvp->thetaV[0]) - dx*dx/(4*expectedGradDiff[i]);
					  // s.val[1] = log(fpvp->thetaV[0]) - dy*dy/(4*expectedGradDiff[i]);
	          s.val[0] = - dx*dx/(4*expectedGradDiff[i]);
	          s.val[1] = - dy*dy/(4*expectedGradDiff[i]);
					  s.val[2] = 0;
					  cvSet2D(dircsgrad[i][z],y,x,s); // set the (i,j) pixel value
				  }

				  if (fpvp->gradVZ == 2)
				  {
					  CvScalar s;
					  s=cvGet2D(dircsgrad[i][z], y,x); // get the (i,j) pixel value
					  // s.val[0] = log(fpvp->thetaV[0]) - dx*dx/(4*expectedGradDiff[i]);
					  // s.val[1] = log(fpvp->thetaV[0]) - dy*dy/(4*expectedGradDiff[i]);
					  // s.val[2] = log(fpvp->thetaV[0]) - dz*dz/(4*expectedGradDiff[i]);
					  s.val[0] = - dx*dx/(4*expectedGradDiff[i]);
	          s.val[1] = - dy*dy/(4*expectedGradDiff[i]);
	          s.val[2] = - dz*dz/(4*expectedGradDiff[i]);
					  cvSet2D(dircsgrad[i][z],y,x,s); // set the (i,j) pixel value
				  }

				  if (fpvp->gradVZ == 3)
				  {
					  CvScalar s;
					  s=cvGet2D(dircsgrad[i][z],y,x); // get the (i,j) pixel value
					  // s.val[0] = log(fpvp->thetaV[0]) - dx*dx/(4*expectedGradDiff[i]);
					  // s.val[1] = log(fpvp->thetaV[0]) - dy*dy/(4*expectedGradDiff[i]);
					  // s.val[2] = log(fpvp->thetaZ[0]) - dz*dz/(4*expectedGradDiffZ[i]);
					  s.val[0] = - dx*dx/(4*expectedGradDiff[i]);
	          s.val[1] = - dy*dy/(4*expectedGradDiff[i]);
	          s.val[2] = - dz*dz/(4*expectedGradDiffZ[i]);
					  cvSet2D(dircsgrad[i][z],y,x,s); // set the (i,j) pixel value
				  }

			  }
		  }
	  }

  }


}

void features::computeCSGradientsOpflow(std::vector <std::vector <CByteImage> > *inp)
{
  // simple2
  double simple2[2] = {-1, +1};

  int numP = (*inp).size();

  assert(fpvp->gradVZ != 1); // it cannot be (x-y) type for optical flow connections

  std::vector<double> expectedGradDiff;
  expectedGradDiff.resize(numP);

  std::vector<double> expectedGradDiffZ;
  expectedGradDiffZ.resize(numP);


  for (int i=0; i< numP; i++)
  {

    int depth = (*inp)[i].size();
    CShape sh = (*inp)[i][0].Shape();
    int width = sh.width, height = sh.height, nB = sh.nBands;
    nB = 3;  // band 0: x gradient, band 1: y gradient, band 2: z gradient

    // computing the expected gradient difference
    expectedGradDiff[i] = 0.0;
    expectedGradDiffZ[i] = 0.0;

    int edges = 0;
    int edgesZ = 0;

    for (int z = 0; z < depth; z++)
    {
      for (int y = 0; y < height; y++)
      {
        for (int x = 0; x < width; x++)
        {
          unsigned char *pix1 = &(*inp)[i][z].Pixel(x, y, 0);
          unsigned char *pix2 = 0;
          unsigned char *pix3 = 0;
          unsigned char *pix4 = 0;
          unsigned char *pix4R = 0;


          double dx=0.0, dy=0.0, dz=0.0, dzR=0.0; // dzR stands for reverse
          double pix1g, pix2g, pix3g, pix4g, pix4Rg;

          pix1g = 0.11*pix1[0] + 0.59*pix1[1] + 0.3*pix1[2];

          if (x==width-1)
            dx = 0;
          else {
            pix2 = &(*inp)[i][z].Pixel(x+1, y, 0);
            pix2g = 0.11*pix2[0] + 0.59*pix2[1] + 0.3*pix2[2];
            dx = fabs(simple2[0]*pix2g + simple2[1]*pix1g);
            edges++;
          }

          if (y==height-1)
            dy = 0;
          else {
            pix3 = &(*inp)[i][z].Pixel(x, y+1, 0);
            pix3g = 0.11*pix3[0] + 0.59*pix3[1] + 0.3*pix3[2];
            dy = fabs(simple2[0]*pix3g + simple2[1]*pix1g);
            edges++;
          }

          if (z==depth-1)
          {
            dz = 0;
            dzR = 0;
          }
          else {

            // use optical flow to get neighbor of [z](x,y)
            float u = flowvector[i][z][0][y*width+x];
            float v = flowvector[i][z][1][y*width+x];

            int new_x = x + int(round(u));
            int new_y = y + int(round(v));

            // check to see that new node has correct connection within image bounds
            // otherwise set exists to 0.
            if (new_x < width && new_y < height && new_x >= 0 && new_y >= 0)
            {
              pix4 = &(*inp)[i][z+1].Pixel(new_x, new_y, 0);
              pix4g = 0.11*pix4[0] + 0.59*pix4[1] + 0.3*pix4[2];
              dz = fabs(simple2[0]*pix4g + simple2[1]*pix1g);

              if (fpvp->gradVZ == 2)
                edges++;
              else if (fpvp->gradVZ == 3)
                edgesZ++;

            } else
            {
              dz = 0;// since point is out of image
            }

            // this is actually using (x,y) points of z+1 and not z
            if (fpvp->opflowreverseconnection == 1)
            {
              // use optical flow to get neighbor of [z](x,y)
              float u = reverseflowvector[i][z][0][y*width+x];
              float v = reverseflowvector[i][z][1][y*width+x];

              int new_x = x + int(round(u));
              int new_y = y + int(round(v));

              // x,y are values on [z+1] slice while new_x, new_y are values on [z] slice
              // because we are in reverse flow now.

              // check to see that new node has correct connection within image bounds
              // otherwise set exists to 0.
              if (new_x < width && new_y < height && new_x >= 0 && new_y >= 0)
              {
                pix4 = &(*inp)[i][z+1].Pixel(x, y, 0);
                pix4R = &(*inp)[i][z].Pixel(new_x, new_y, 0);
                pix4g = 0.11*pix4[0] + 0.59*pix4[1] + 0.3*pix4[2];
                pix4Rg = 0.11*pix4R[0] + 0.59*pix4R[1] + 0.3*pix4R[2];
                dzR = fabs(simple2[0]*pix4g + simple2[1]*pix4Rg);

                if (fpvp->gradVZ == 2)
                  edges++;
                else if (fpvp->gradVZ == 3)
                  edgesZ++;

              } else
              {
                dzR = 0;// since point is out of image
              }
            }


          }


          if (fpvp->gradVZ == 2) // (x-y-z)
          {
            expectedGradDiff[i] += dx*dx + dy*dy + dz*dz + dzR*dzR;
          }
          else if (fpvp->gradVZ == 3) // (x-y and z)
          {
            expectedGradDiff[i] += dx*dx + dy*dy;
            expectedGradDiffZ[i] += dz*dz + dzR*dzR;
          }

        }
      }
    }

    if (edges != 0)
      expectedGradDiff[i] /= edges;
    if (edgesZ != 0 && fpvp->gradVZ == 3)
      expectedGradDiffZ[i] /= edgesZ;

  } // i

  // Important thing to note is that the expected values in the flow direction may lie between 1* to 2*
  // the real expectation values. This is because while computing the reverse flows, they might already
  // be the same as the forward flow and so it will be counted twice.

  // However, during graph contrucntion iw storing in the map below (edgeFlowPV) we check to see so that
  // we don't add duplicates



  // then computing gradients
  for (int i=0; i< numP; i++)
  {
    std::map <std::pair <int, int>, std::pair<int, MRF::CostVal> > edgeFlowPV; // per volume

    int depth = (*inp)[i].size();
    CShape sh = (*inp)[i][0].Shape();
    int width = sh.width, height = sh.height, nB = sh.nBands;
    nB = 3;  // band 0: x gradient, band 1: y gradient, band 2: z gradient

    assert(expectedGradDiff[i]!=0);
    if (fpvp->gradVZ == 3 && depth > 1)
      assert(expectedGradDiffZ[i]!=0);

    assert(fpvp->thetaV[0] > 0);
    if (fpvp->gradVZ == 3)
      assert(fpvp->thetaZ[0] > 0);

    int noden = 0;
    for (int z = 0; z < depth; z++)
    {
      for (int y = 0; y < height; y++)
      {
        for (int x = 0; x < width; x++)
        {

          int exists = 0;
          std::pair <int, int> pnodes;
          std::pair <int, MRF::CostVal> pvals;
          int new_noden = 0;

          unsigned char *pix1 = &(*inp)[i][z].Pixel(x, y, 0);
          unsigned char *pix2 = 0;
          unsigned char *pix3 = 0;
          unsigned char *pix4 = 0;
          unsigned char *pix4R = 0;

          double dx=0.0, dy=0.0, dz=0.0, dzR = 0.0;
          double pix1g, pix2g, pix3g, pix4g, pix4Rg;

          pix1g = 0.11*pix1[0] + 0.59*pix1[1] + 0.3*pix1[2];

          if (x==width-1)
            dx = 0;
          else {
            pix2 = &(*inp)[i][z].Pixel(x+1, y, 0);
            pix2g = 0.11*pix2[0] + 0.59*pix2[1] + 0.3*pix2[2];
            dx = fabs(simple2[0]*pix2g + simple2[1]*pix1g);

            pnodes.first = noden;
            pnodes.second = noden + 1;
            pvals.first = 0; //spatial
            // pvals.second = log(fpvp->thetaV[0]) - dx*dx/(4*expectedGradDiff[i]);
            pvals.second = - dx*dx/(4*expectedGradDiff[i]);
            edgeFlowPV.insert(std::pair<std::pair <int, int>, std::pair <int, MRF::CostVal> >(pnodes, pvals));
          }

          if (y==height-1)
            dy = 0;
          else {
            pix3 = &(*inp)[i][z].Pixel(x, y+1, 0);
            pix3g = 0.11*pix3[0] + 0.59*pix3[1] + 0.3*pix3[2];
            dy = fabs(simple2[0]*pix3g + simple2[1]*pix1g);

            pnodes.first = noden;
            pnodes.second = noden + width;
            pvals.first = 0; //spatial
            // pvals.second = log(fpvp->thetaV[0]) - dy*dy/(4*expectedGradDiff[i]);
            pvals.second = - dy*dy/(4*expectedGradDiff[i]);
            edgeFlowPV.insert(std::pair<std::pair <int, int>, std::pair <int, MRF::CostVal> >(pnodes, pvals));
          }

          if (z==depth-1)
          {
            dz = 0;
            dzR = 0;
          }
          else {

            // use optical flow to get neighbor of [z](x,y)
            float u = flowvector[i][z][0][y*width+x];
            float v = flowvector[i][z][1][y*width+x];

            int new_x = x + int(round(u));
            int new_y = y + int(round(v));


            // check to see that new node has correct connection within image bounds
            // otherwise set exists to 0.
            if (new_x < width && new_y < height && new_x >= 0 && new_y >= 0)
            {// valid
              pix4 = &(*inp)[i][z+1].Pixel(new_x, new_y, 0);
              pix4g = 0.11*pix4[0] + 0.59*pix4[1] + 0.3*pix4[2];
              dz = fabs(simple2[0]*pix4g + simple2[1]*pix1g);
              new_noden = (width*height*(z+1)) + (new_y*width + new_x);

              pnodes.first = noden;
              pnodes.second = new_noden;

              int notthere = 1;

              // if the reverse connection is not enabled we are guaranteed that edgeFlowPV does not
              // have the pnodes key and save computation time by not searching in the map.
              // and so check only if reverse connection is enabled
              if (fpvp->opflowreverseconnection == 1)
              {
                std::pair <int, int> pnodesR;
                pnodesR.first = new_noden;
                pnodesR.second = noden;

                if (edgeFlowPV.count(pnodes) > 0)
                  notthere = 0;
                if (edgeFlowPV.count(pnodesR) > 0)
                  notthere = 0;
              }

              if (notthere == 1)
              {

                pvals.first = 1; //optical

                if (fpvp->gradVZ == 2)
                {
                  // pvals.second = log(fpvp->thetaV[0]) - dz*dz/(4*expectedGradDiff[i]);
                  pvals.second = - dz*dz/(4*expectedGradDiff[i]);
                }
                if (fpvp->gradVZ == 3)
                {
                  // pvals.second = log(fpvp->thetaZ[0]) - dz*dz/(4*expectedGradDiffZ[i]);
                  pvals.second = - dz*dz/(4*expectedGradDiffZ[i]);
                }
                edgeFlowPV.insert(std::pair<std::pair <int, int>, std::pair <int, MRF::CostVal> >(pnodes, pvals));

              }
            }

            if (fpvp->opflowreverseconnection == 1)
            {

              // use optical flow to get neighbor of [z](x,y)
              float u = reverseflowvector[i][z][0][y*width+x];
              float v = reverseflowvector[i][z][1][y*width+x];

              int new_x = x + int(round(u));
              int new_y = y + int(round(v));

              // x,y belong to slice z+1
              // new_x, new_y belong to slice z


              // check to see that new node has correct connection within image bounds
              // otherwise set exists to 0.
              if (new_x < width && new_y < height && new_x >= 0 && new_y >= 0)
              {// valid

                pix4 = &(*inp)[i][z+1].Pixel(x, y, 0);
                pix4R = &(*inp)[i][z].Pixel(new_x, new_y, 0);
                pix4g = 0.11*pix4[0] + 0.59*pix4[1] + 0.3*pix4[2];
                pix4Rg = 0.11*pix4R[0] + 0.59*pix4R[1] + 0.3*pix4R[2];
                dzR = fabs(simple2[0]*pix4g + simple2[1]*pix4Rg);

                int curr_node = noden + width*height; // because its on z+1
                new_noden = (width*height*(z)) + (new_y*width + new_x); // on frame z

                pnodes.first = new_noden;
                pnodes.second = curr_node;

                std::pair <int, int> pnodesR;
                pnodesR.first = curr_node;
                pnodesR.second = new_noden;

                if (edgeFlowPV.count(pnodes) == 0 && edgeFlowPV.count(pnodesR) == 0)
                {
                  pvals.first = 1; //optical

                  if (fpvp->gradVZ == 2)
                  {
                    // pvals.second = log(fpvp->thetaV[0]) - dzR*dzR/(4*expectedGradDiff[i]);
                    pvals.second = - dzR*dzR/(4*expectedGradDiff[i]);
                  }
                  if (fpvp->gradVZ == 3)
                  {
                    // pvals.second = log(fpvp->thetaZ[0]) - dzR*dzR/(4*expectedGradDiffZ[i]);
                    pvals.second = - dzR*dzR/(4*expectedGradDiffZ[i]);
                  }
                  edgeFlowPV.insert(std::pair<std::pair <int, int>, std::pair <int, MRF::CostVal> >(pnodes, pvals));
                }
              }


            }


          }

          noden++;
        }
      }
    }

    globalP.edgeGlobalF.push_back(edgeFlowPV);
  }

}





// RGB converted to grayscale (x-y or x-y-z or (x-y and z) )
void features::computeQuantizedGradients(std::vector <std::vector <CByteImage> > *inp)
{

  // simple2
  double simple2[2] = {-1, +1};

  int k, kz;

  int numP = (*inp).size();

  // allocating space first
  std::vector<CByteImage> setImgrad;

  for (int i=0; i< numP; i++)
  {

	  setImgrad.erase(setImgrad.begin(), setImgrad.end());

	  int depth = (*inp)[i].size();
	  setImgrad.resize(depth);
	  dirgrad.push_back(setImgrad);

  }



  // then computing gradients
  for (int i=0; i< numP; i++)
  {

	  int depth = (*inp)[i].size();
	  CShape sh = (*inp)[i][0].Shape();
	  int width = sh.width, height = sh.height, nB = sh.nBands;
	  sh.nBands = 3;  // band 0: x gradient, band 1: y gradient, band 2: z gradient

	  for (int z = 0; z<depth; z++) {
		  dirgrad[i][z].ReAllocate(sh);
	  }


	  for (int z = 0; z < depth; z++)
	  {
		  for (int y = 0; y < height; y++)
		  {
			  for (int x = 0; x < width; x++)
			  {

				  uchar *grad   = &dirgrad[i][z].Pixel(x, y, 0);

				  unsigned char *pix1 = &(*inp)[i][z].Pixel(x, y, 0);
				  unsigned char *pix2 = 0;
				  unsigned char *pix3 = 0;
				  unsigned char *pix4 = 0;

				  double dx=0.0, dy=0.0, dz=0.0;
				  double pix1g, pix2g, pix3g, pix4g;

				  pix1g = 0.11*pix1[0] + 0.59*pix1[1] + 0.3*pix1[2];

				  if (x==width-1)
					  dx = 0;
				  else {
					  pix2 = &(*inp)[i][z].Pixel(x+1, y, 0);
					  pix2g = 0.11*pix2[0] + 0.59*pix2[1] + 0.3*pix2[2];
					  dx = fabs(simple2[0]*pix2g + simple2[1]*pix1g);
				  }

				  if (y==height-1)
					  dy = 0;
				  else {
					  pix3 = &(*inp)[i][z].Pixel(x, y+1, 0);
					  pix3g = 0.11*pix3[0] + 0.59*pix3[1] + 0.3*pix3[2];
					  dy = fabs(simple2[0]*pix3g + simple2[1]*pix1g);
				  }

				  if (z==depth-1)
					  dz = 0;
				  else {
					  pix4 = &(*inp)[i][z+1].Pixel(x, y, 0);
					  pix4g = 0.11*pix4[0] + 0.59*pix4[1] + 0.3*pix4[2];
					  dz = fabs(simple2[0]*pix4g + simple2[1]*pix1g);
				  }


				  for (k = 0; fpvp->gradThreshVec[k] <= dx; k++)
					  ; // find lowest bin for x gradient
				  assert(k < fpvp->nG);
				  grad[0] = k;

				  for (k = 0; fpvp->gradThreshVec[k] <= dy; k++)
					  ; // find lowest bin for y gradient
				  assert(k < fpvp->nG);
				  grad[1] = k;

				  if (fpvp->gradVZ == 1) // (x-y)
				  {

				  } else if (fpvp->gradVZ == 2) // (x-y-z)
				  {
					  for (k = 0; fpvp->gradThreshVec[k] <= dz; k++)
						  ; // find lowest bin for z gradient
					  assert(k < fpvp->nG);
					  grad[2] = k;
				  } else if (fpvp->gradVZ == 3) // (x-y and z)
				  {
					  for (kz = 0; fpvp->gradThreshVecZ[k] <= dz; kz++)
						  ; // find lowest bin for z gradient
					  assert(kz < fpvp->nGZ);
					  grad[2] = kz;
				  }

			  }
		  }
	  }
  }
}



void features::computeQuantizedGradientsOpflow(std::vector <std::vector <CByteImage> > *inp)
{
  // simple2
  double simple2[2] = {-1, +1};
  int numP = (*inp).size();

  int k, kz;

  // then computing gradients
  for (int i=0; i< numP; i++)
  {
    std::map <std::pair <int, int>, std::pair<int, MRF::CostVal> > edgeFlowPV; // per volume
    int depth = (*inp)[i].size();
    CShape sh = (*inp)[i][0].Shape();
    int width = sh.width, height = sh.height, nB = sh.nBands;
    sh.nBands = 3;  // band 0: x gradient, band 1: y gradient, band 2: z gradient

    int noden = 0;

    for (int z = 0; z < depth; z++)
    {
      for (int y = 0; y < height; y++)
      {
        for (int x = 0; x < width; x++)
        {

          int exists = 0;
          std::pair <int, int> pnodes;
          std::pair <int, MRF::CostVal> pvals;
          int new_noden = 0;

          unsigned char *pix1 = &(*inp)[i][z].Pixel(x, y, 0);
          unsigned char *pix2 = 0;
          unsigned char *pix3 = 0;
          unsigned char *pix4 = 0;

          double dx=0.0, dy=0.0, dz=0.0;
          double pix1g, pix2g, pix3g, pix4g;

          pix1g = 0.11*pix1[0] + 0.59*pix1[1] + 0.3*pix1[2];

          if (x==width-1)
          {
            exists = 0;
            dx = 0;
          }
          else {
            exists = 1;
            pix2 = &(*inp)[i][z].Pixel(x+1, y, 0);
            pix2g = 0.11*pix2[0] + 0.59*pix2[1] + 0.3*pix2[2];
            dx = fabs(simple2[0]*pix2g + simple2[1]*pix1g);
          }

          if (exists == 1)
          {
            for (k = 0; fpvp->gradThreshVec[k] <= dx; k++)
              ; // find lowest bin for x gradient
            assert(k < fpvp->nG);
            pnodes.first = noden;
            pnodes.second = noden + 1;
            pvals.first = 0; //spatial
            pvals.second = k; // grad bin
            edgeFlowPV.insert(std::pair<std::pair <int, int>, std::pair <int, MRF::CostVal> >(pnodes, pvals));
          }


          if (y==height-1)
          {
            exists = 0;
            dy = 0;
          }
          else {
            exists = 1;
            pix3 = &(*inp)[i][z].Pixel(x, y+1, 0);
            pix3g = 0.11*pix3[0] + 0.59*pix3[1] + 0.3*pix3[2];
            dy = fabs(simple2[0]*pix3g + simple2[1]*pix1g);
          }


          if (exists == 1)
          {
            for (k = 0; fpvp->gradThreshVec[k] <= dy; k++)
              ; // find lowest bin for y gradient
            assert(k < fpvp->nG);
            pnodes.first = noden;
            pnodes.second = noden + width;
            pvals.first = 0; //spatial
            pvals.second = k; // grad bin
            edgeFlowPV.insert(std::pair<std::pair <int, int>, std::pair <int, MRF::CostVal> >(pnodes, pvals));
          }


          if (z==depth-1)
          {
            exists = 0;
            dz = 0;
          }
          else {
            exists = 1;

            // use optical flow to get neighbor of [z](x,y)
            float u = flowvector[i][z][0][y*width+x];
            float v = flowvector[i][z][1][y*width+x];

            int new_x = x + int(round(u));
            int new_y = y + int(round(v));

            // check to see that new node has correct connection within image bounds
            // otherwise set exists to 0.
            if (new_x < width && new_y < height && new_x >= 0 && new_y >= 0)
            {// valid
              pix4 = &(*inp)[i][z+1].Pixel(new_x, new_y, 0);
              pix4g = 0.11*pix4[0] + 0.59*pix4[1] + 0.3*pix4[2];
              dz = fabs(simple2[0]*pix4g + simple2[1]*pix1g);

              new_noden = (width*height*(z+1)) + (new_y*width + new_x);

            } else {
              exists = 0;
            }
          }

          assert(fpvp->gradVZ != 1); // it cannot be (x-y) type

          if (exists == 1)
          {
            pnodes.first = noden;
            pnodes.second = new_noden;
            pvals.first = 1; //optical

            if (fpvp->gradVZ == 2) // (x-y-z)
            {
              for (k = 0; fpvp->gradThreshVec[k] <= dz; k++)
                ; // find lowest bin for z gradient
              assert(k < fpvp->nG);
              pvals.second = k; // grad bin
              edgeFlowPV.insert(std::pair<std::pair <int, int>, std::pair <int, MRF::CostVal> >(pnodes, pvals));
            }
            else if (fpvp->gradVZ == 3) // (x-y and z)
            {
              for (kz = 0; fpvp->gradThreshVecZ[k] <= dz; kz++)
                ; // find lowest bin for z gradient
              assert(kz < fpvp->nGZ);
              pvals.second = kz; // grad bin
              edgeFlowPV.insert(std::pair<std::pair <int, int>, std::pair <int, MRF::CostVal> >(pnodes, pvals));
            }
          }
          noden++;
        } // x
      } // y
    } // z
    globalP.edgeGlobalF.push_back(edgeFlowPV);
  }

}



// we have parameters for each cluster,class pair.
// ie each class will have the parameters for each cluster (even if some clusters don't belong to the class)
// this will try to make the parameters of the clusters belonging to the class (more so the one cluster closest) with a higher impact value
// than the rest of the cluster.
double features::getCostLocationDiscGaussian3(int pno, unsigned int z, int x, int y, int d, int genparam, int nD, int &thetaid)
{

	int nL = fpvp->thetaL.size();
	int nLpC = nL/nD;

	double loccost = 0.0;

	CvScalar s;
	s=cvGet2D(locdirImageClust[pno][z],y,x); // get the (i,j) pixel value

	// find the theta corresponding to c
	thetaid = nLpC*d + s.val[0];
	loccost = fpvp->thetaL[thetaid];

	return loccost;

}











// this function is called only when timeseq == 0
void features::computeLocation(std::vector <std::vector <CImage2> > *indirImage, int nD)
{

	if (fpvp->loc == 7) // this is a hard assignment used from older getCostLocationDiscGaussian2
		// for the soft version loc==9, you will need to run older getCostLocationDiscGaussian2
		// because we need the value maxLocationMVval
	{

		// allocate space
		locdirImageClust.erase(locdirImageClust.begin(), locdirImageClust.end());
		int numP = (*indirImage).size();

		for (int i=0; i< numP; i++)
		{
			int depth = (*indirImage)[i].size();
			CShape sh = (*indirImage)[i][0].Shape();
			int width = sh.width, height = sh.height, nB = sh.nBands;
			nB = 1;

			std::vector <IplImage*> setImLoc;
			setImLoc.erase(setImLoc.begin(), setImLoc.end());

			for (int j=0; j<depth; j++)
			{
				IplImage* img = cvCreateImage( cvSize(width,height), IPL_DEPTH_16U, nB );
				setImLoc.push_back(img);
			}

			locdirImageClust.push_back(setImLoc);

		}


		/*
		locdirImage.resize(numPatients);
		for (int mm=0; mm<numPatients; mm++)
		{
			int numSlices = (*indirImage)[mm].size();
			locdirImage[mm].resize(numSlices);
			CShape sh = (*indirImage)[mm][0].Shape();
			for (int z = 0; z<numSlices; z++) {
				//locdirImage[mm][z].ReAllocate(sh);
			}
		}
		*/


		int nL = fpvp->thetaL.size();
		int nLpC = nL/nD;

		float loccost=0;
		float *pii;
		float **mu;
		float ***SigmaInv;
		float *normConst;


		for (int i=0; i< numP; i++)
		{
			int depth = (*indirImage)[i].size();
			CShape sh = (*indirImage)[i][0].Shape();
			int width = sh.width, height = sh.height;

			int zslice = startSliceNo[i];


			for (int z = 0; z < depth; z++)
			{
				for (int y = 0; y < height; y++)
				{
					for (int x = 0; x < width; x++)
					{

						double LocationMVval = 0.0;
						int maxclass = 0;
						double maxLocationMVval = 0.0;

						int totClusters = 0;

						for (int ff=0; ff<nD; ff++)
						{

							pii = LocationMVclass[ff]->getPi();
							mu = LocationMVclass[ff]->getMu();
							SigmaInv = LocationMVclass[ff]->getSigmaInv();
							normConst = LocationMVclass[ff]->getNormConst();


							for (int c=0; c<LocationMVclass[ff]->getNClusters(); c++)
							{

								float zs = z + zslice - mu[c][2];
								float xs = y-mu[c][0];
								float ys = x-mu[c][1];

								double tempval = xs*xs*SigmaInv[c][0][0] + xs*ys*(SigmaInv[c][1][0] + SigmaInv[c][0][1]) + xs*zs*(SigmaInv[c][2][0] + SigmaInv[c][0][2]) +
								  ys*ys*SigmaInv[c][1][1] + ys*zs*(SigmaInv[c][2][1] + SigmaInv[c][1][2]) + zs*zs*SigmaInv[c][2][2];

								LocationMVval = ((-tempval/2)-log(normConst[c])); // no need to take exp
								if (c==0 && ff==0) {
									maxclass = 0;
									maxLocationMVval = LocationMVval;
								} else {
									if (maxLocationMVval<LocationMVval) {
										maxclass = totClusters;
										maxLocationMVval = LocationMVval;
									}
								}
								totClusters++;

							}

						}

						//unsigned short *pix1 = 0;
						//pix1 = &locdirImage[i][z].Pixel(x, y, 0);
						//*pix1 = maxclass;

						CvScalar s;
						s=cvGet2D(locdirImageClust[i][z],y,x); // get the (i,j) pixel value
						s.val[0] = maxclass;
						cvSet2D(locdirImageClust[i][z],y,x,s); // set the (i,j) pixel value

					}
				}
			}
		}

//		if (totClusters != nLpC) // sanity error check
//		{
//			  std::cout<<"Error in number of clusters for location"<<endl;
//			  exit(1);
//		}

	}

}

int features::deAllocateLocationFiles(int timeseq)
{
	if (timeseq == 0 && fpvp->loc == 7)
	{
/*		int numD = locdirImage.size();

		for (int i=0; i<numD; i++)
		{
			int numD2 = locdirImage[i].size();
			for (int j=0; j<numD2; j++)
			{
				// std::cout << "Deallocating " <<i <<"  " << j <<"\n";
				locdirImage[i][j].DeAllocate();
			}
		}
*/

		int numP = locdirImageClust.size();

		// deallocating space
		for (int i=0; i< numP; i++)
		{
			int depth = locdirImageClust[i].size();

			for (int j=0; j<depth; j++)
			{
					cvReleaseImage(&locdirImageClust[i][j]);
			}

		}

	}

}



int features::preComputeLocalFeatures(std::vector <std::vector <CImage2> > *inp, std::vector <std::vector <CByteImage> > *gtr, std::vector<int> testDirIndexV, int nD)
{

	// intensity
	computeIntensityNB(inp, gtr, testDirIndexV, nD);

	// appearance
	computeAppearance(inp, gtr, nD);

	// location
	computeLocation(inp, nD);

}


int features::preComputeLocalFeatures(std::vector <std::vector <CByteImage> > *inp, std::vector <std::vector <CByteImage> > *gtr, std::vector<int> testDirIndexV, int nD)
{

	// intensity
	computeIntensityNB(inp, gtr, testDirIndexV, nD);

	// appearance
	computeAppearance(inp, gtr, nD);

}




int features::preComputePairwiseFeatures(std::vector <std::vector <CImage2> > *inp)
{
	if (fpvp->csgrad == 0)
		computeQuantizedGradients(inp);
	else if (fpvp->csgrad == 1)
		computeCSGradients(inp);
}


int features::preComputePairwiseFeatures(std::vector <std::vector <CByteImage> > *inp)
{
	if (fpvp->csgrad == 0)
		computeQuantizedGradients(inp);
	else if (fpvp->csgrad == 1)
		computeCSGradients(inp);
}

int features::preComputePairwiseFeaturesAndPutInGlobalSpace(std::vector <std::vector <CByteImage> > *inp)
{
  globalP.setOpflowConnection(1);

  if (fpvp->csgrad == 0)
    computeQuantizedGradientsOpflow(inp);
  else if (fpvp->csgrad == 1)
    computeCSGradientsOpflow(inp);
}


void features::writeOpflowCache()
{
  // forward flow
  int numP = flowvector.size();
  std::vector<int> numf;
  for (int i = 0; i <numP; i++)
    numf.push_back(flowvector[i].size());

  // reverse flow
  int numPr =reverseflowvector.size();
  std::vector<int> numfr;
  for (int i = 0; i <numPr; i++)
    numfr.push_back(reverseflowvector[i].size());

  ofstream myfileo(fpvp->flowcachefile.c_str(), ios::out | ios::binary);
  if(myfileo.is_open())
  {
    myfileo.write((char *)&numP, sizeof(int));
    for (int i = 0; i <numP; i++)
      myfileo.write((char *)&numf[i], sizeof(int));

    for (int i=0; i < numP; i++)
    {
      for (int j = 0; j < numf[i]; j++)
      {
        bool foou = flowvector[i][j][0].saveImageData(myfileo);
        bool foov = flowvector[i][j][1].saveImageData(myfileo);
        assert(foou && foov);
      }
    }
    myfileo.close();
  }
  else
  {
    std::cout << " Error writing opflow cache file \n";
    exit(1);
  }


  ofstream myfileor(fpvp->rflowcachefile.c_str(), ios::out | ios::binary);
  if(myfileor.is_open())
  {
    myfileor.write((char *)&numPr, sizeof(int));
    for (int i = 0; i <numPr; i++)
      myfileor.write((char *)&numfr[i], sizeof(int));

    for (int i=0; i < numPr; i++)
    {
      for (int j = 0; j < numfr[i]; j++)
      {
        bool foou = reverseflowvector[i][j][0].saveImageData(myfileor);
        bool foov = reverseflowvector[i][j][1].saveImageData(myfileor);
        assert(foou && foov);
      }
    }
    myfileor.close();
  }
  else
  {
    std::cout << " Error writing opflow cache file \n";
    exit(1);
  }

}


void features::readOpflowCache()
{
  flowvector.clear();
  reverseflowvector.clear();

  int numP = 0, temp;
  std::vector<int> numf;
  int numPr = 0;
  std::vector<int> numfr;

  ifstream myfilei(fpvp->flowcachefile.c_str(), ios::in | ios::binary);
  if(myfilei.is_open())
  {
    myfilei.read((char *)&numP, sizeof(int));
    for (int i = 0; i <numP; i++)
    {
      myfilei.read((char *)&temp, sizeof(int));
      numf.push_back(temp);
    }

    for (int i=0; i < numP; i++)
    {
      std::vector<std::vector <opflowns::DImage> > tempfvarr(0);
      for (int j = 0; j < numf[i]; j++)
      {
        opflowns::DImage vx,vy;
        vx.loadImageData(myfilei);
        vy.loadImageData(myfilei);
        std::vector <opflowns::DImage> tempuv(0);
        tempuv.push_back(vx);
        tempuv.push_back(vy);
        tempfvarr.push_back(tempuv);
      }
      flowvector.push_back(tempfvarr);
    }

    myfilei.close();
  }
  else
  {
    std::cout << " Error writing opflow cache file \n";
    exit(1);
  }


  ifstream myfileir(fpvp->rflowcachefile.c_str(), ios::in | ios::binary);
  if(myfileir.is_open())
  {
    myfileir.read((char *)&numPr, sizeof(int));
    for (int i = 0; i <numPr; i++)
    {
      myfileir.read((char *)&temp, sizeof(int));
      numfr.push_back(temp);
    }

    for (int i=0; i < numPr; i++)
    {
      std::vector<std::vector <opflowns::DImage> > tempfvarr(0);
      for (int j = 0; j < numfr[i]; j++)
      {
        opflowns::DImage vx,vy;
        vx.loadImageData(myfileir);
        vy.loadImageData(myfileir);
        std::vector <opflowns::DImage> tempuv(0);
        tempuv.push_back(vx);
        tempuv.push_back(vy);
        tempfvarr.push_back(tempuv);
      }
      reverseflowvector.push_back(tempfvarr);
    }


    myfileir.close();
  }
  else
  {
    std::cout << " Error reading opflow cache file \n";
    exit(1);
  }

}


void features::extractOpticalFlow(std::vector<string> &indirpath)
{
  if (fpvp->useopflowcache == 0)
  {
    preComputeOpticalFlow(indirpath);
  } else {
    // if files exist then no need to compute optical flow.
    int need_to_compute = 0;
    ifstream myfilei(fpvp->flowcachefile.c_str(), ios::in | ios::binary);
    if(myfilei.is_open())
      myfilei.close();
    else
      need_to_compute = 1;

    ifstream myfileir(fpvp->rflowcachefile.c_str(), ios::in | ios::binary);
    if(myfileir.is_open())
      myfileir.close();
    else
      need_to_compute = 1;

    if (need_to_compute == 1)
    {
      preComputeOpticalFlow(indirpath);
      // write to files
      std::cout << " Writing optical flow cache file " << std::endl;
      writeOpflowCache();
    } else
    {
      std::cout << " Reading optical flow cache file " << std::endl;
      // read files
      readOpflowCache();
      // consistency check regarding input files and cache files
      checkConsistencyCacheFiles(indirpath);
    }

  }
}

void features::checkConsistencyCacheFiles(std::vector<string> &indirpath)
{

  int vid = 0;
  for(std::vector<string>::iterator b = indirpath.begin(); b < indirpath.end(); ++b)
  {
    std::vector<string> filepaths(0);
    readImagesPaths(*b, filepaths);
    int numf = filepaths.size();
    if (flowvector[vid].size() != numf - 1)
    {
      std::cout << " flow cache file consistency error !! \n";
      assert(false);
    }

    if (reverseflowvector[vid].size() != numf - 1)
    {
      std::cout << " reverse cache file consistency error !! \n";
      assert(false);
    }

    vid++;
  }

}


int features::preComputeOpticalFlow(std::vector<string> &indirpath)
{

  std::cout << " Computing optical flow ... \n";

  flowvector.clear();
  reverseflowvector.clear();

  // current settings
  double alpha = 0.012;
  double ratio = 0.75;
  int minWidth = 20;
  int nOuterFPIterations = 7;
  int nInnerFPIterations = 1;
  int nSORIterations = 30;

  for(std::vector<string>::iterator b = indirpath.begin(); b < indirpath.end(); ++b)
  {
    std::vector<string> filepaths(0);
    readImagesPaths(*b, filepaths);
    std::vector<opflowns::DImage> opinput(0);

    for(std::vector<string>::iterator bf = filepaths.begin(); bf < filepaths.end(); ++bf)
    {
      cout << *bf << "\n";
      opflowns::DImage Im1;
      Im1.imread(bf->c_str());
      opinput.push_back(Im1);
    }

    std::vector<std::vector <opflowns::DImage> > tempfvarr(0);
    std::vector<std::vector <opflowns::DImage> > tempfRvarr(0);

    int depth = opinput.size();
    for (int j=0; j<depth-1; j++)
    {
      opflowns::DImage vx,vy,warpI2;
      opflowns::OpticalFlow::Coarse2FineFlow(vx,vy,warpI2,opinput[j],opinput[j+1],alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations);
      std::vector <opflowns::DImage> tempuv(0);
      tempuv.push_back(vx);
      tempuv.push_back(vy);
      tempfvarr.push_back(tempuv);
    }
    flowvector.push_back(tempfvarr);

    if (fpvp->opflowreverseconnection == 1 || fpvp->useopflowcache == 1)
    {
      for (int j=1; j<depth; j++)
      {
        opflowns::DImage vx,vy,warpI2;
        opflowns::OpticalFlow::Coarse2FineFlow(vx,vy,warpI2,opinput[j],opinput[j-1],alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations);
        std::vector <opflowns::DImage> tempuv(0);
        tempuv.push_back(vx);
        tempuv.push_back(vy);
        tempfRvarr.push_back(tempuv);
      }
      reverseflowvector.push_back(tempfRvarr);
    }

  }

}




void features::deleteLocationMVclass(int nD)
{

	if (fpvp->loc == 3 || fpvp->loc==6 || fpvp->loc==8 || fpvp->loc==7 || fpvp->loc==9)
	{
		for (int i=0; i<nD; i++)
		{
			delete LocationMVclass[i];
		}
		delete [] LocationMVclass;
	}

}

void features::deleteAppMVclass(int nD)
{

	if (fpvp->app==3)
	{
		for (int i=0; i<nD; i++)
		{
			delete AppMVclass[i];
		}
		delete AppMVclass;
	}

}

void features::deleteappclass(int nD)
{

	if (fpvp->app==1)
	{
		for (int i=0; i<nD; i++)
		{
			delete appclass[i];
		}
		delete appclass;
	}

	if (fpvp->app==2)
	{
		for (int i=0; i<1; i++)
		{
			delete appclass[i];
		}
		delete appclass;

	}

}

void features::deleteintObj(int nD)
{

	if (fpvp->intensity==3)
	{
		for (int i=0; i<nD; i++)
		{
			delete intObj[i];
		}
		delete intObj;
	}

}






/*

void readHoGProb(std::vector <std::vector <CImage2> > im1, int nD, std::vector <std::vector <std::vector <matrixB<double> > > > &hogdirMVProb, std::vector<int> testDirIndexV)
{


  for (unsigned int j = 0; j < im1.size(); ++j) {

    int depth = im1[j].size();

    CShape sh = im1[j][0].Shape();
    int width = sh.width, height = sh.height, nB = sh.nBands;

    //WORK-THIS - how to take into account multiple patients etc in an elegant way

    //WORK-THIS assume the patient is pa13, will need to modify this code
    // to account for multiple patients

    char cname[50];
    // need to work on testDirIndex
    if (testDirIndexV[j]==1)
      strcpy(cname, "test/");
    else {

      if (j==0)
        strcpy(cname, "pa13/");
      else if (j==1)
        strcpy(cname, "test/"); //WORK-THIS will need to change depend on the data

    }

    for (int d = 0; d < nD; d++) {

      char fname[200];
      strcpy(fname, "../art/HoGgen30/" );
      strcat(fname, cname);

      if (d==0)
        strcat(fname,"hogprobbgnd.txt");
      else if (d==1)
        strcat(fname,"hogprobliv.txt");
      else if (d==2)
        strcat(fname,"hogprobrk.txt");
      else if (d==3)
        strcat(fname,"hogproblk.txt");
      else if (d==4)
        strcat(fname,"hogprobgb.txt");
      else if (d==5)
        strcat(fname,"hogprobspl.txt");

      std::cout << "filename : "<<fname<<endl;

      ifstream myfile (fname);

      if (myfile.is_open())
        {
          string line;

          int z=0, y=0, x=0;  //y is for height or row, x is for width or column
          matrixB<double> improb(height, width);

          while( getline( myfile, line ) )
            {
              tokenizer<escaped_list_separator<char> > tok(line);
              for(tokenizer<escaped_list_separator<char> >::iterator beg=tok.begin(); beg!=tok.end(); ++beg){
                improb(y,x) = boost::lexical_cast<double>(*beg);
                x++;
                if (x>width) {
                  std::cout <<" Error in the data in hog param file \n";
                  exit(1);
                }
              }
              y++; x=0;
              if (y>=height) {
                hogdirMVProb[j][z].push_back(improb);
                z++; y=0; x=0;
              }
            }
          myfile.close();
        } else {
        std::cout << "Unable to open hogprob parameter file";
        exit(1);
      }

    }

  }
}


 */
