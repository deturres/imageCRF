/*
 * modelmain.cpp
 *
 *  Created on: Mar 15, 2011
 *      Author: bhole
 */

/*
*
* Copyright 2013 Chetan Bhole
*
* This file is part of imageCRF.
*
*    imageCRF is free software: you can redistribute it and/or modify
*    it under the terms of the GNU Lesser General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    imageCRF is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public License
*    along with imageCRF.  If not, see <http://www.gnu.org/licenses/>.
*    
*    
* Please make sure to see the individual licences in the packages used for this code :
* boost, eigen, imageLib, libconfig, MRF2.0float, opticalflow, randgen and zlib.
* They all do not follow LGPL. So you might need to alternates depending on your needs. 
*
*/ 



#include "modelmain.h"

// global variables for debugging output
extern int verbose;
extern FILE *debugfile;

int main(int argc, char **argv)
{

	// temporary fix for verbose and debugfile needed in crfmodel.cpp
  verbose = 0;

  time_t start, end1, end2;
  double timediff=0.0;

  srand(time(0));
  time(&start);

  Parameters prm;
  Config cfg;

  if (argc < 3)
  {
    std::cerr << " You need to pass the parameter config file as the first parameter and debug file as the second parameter " <<std::endl;
    std::cerr << " ./modelmain path_to_cfg_file path_to_debug_file " << std::endl;
    std::cerr << " The debug file usage will be deprecated but currently must be provided." << std::endl;
    exit(0);
  }

  debugfile = fopen(argv[2], "w");
  prm.readfromcfg(cfg, argv[1]);

  prm.readfromCommandline(argc, argv); // no need to call update
  prm.readfeaturefile();  // read parameters of crf
  prm.updatePropertiesParameters();

  prm.printParameters();

  std::vector<float> theta;
  prm.returnConcatThetaVector(&theta);

  baseModel* bM;

  if (prm.logregpv.logreg == 1)
	  bM = new logreg(&prm, theta);
  else if (prm.crfpv.crfp == 1)
	  bM = new crf(&prm, theta);
  else if (prm.logregpv.logregpl == 1)
	  bM = new logregpl(&prm, theta);
  else
  {
	  std::cout << " Not yet implemented ! \n";
	  exit(1);
  }



  // IO (input, gtruth, interactive, output)
  bM->readIOfolderStructures(prm.iopv, prm.interactive, prm.crfpv.hidden);
  bM->readIOfiles(prm.interactive, prm.crfpv.hidden);
  bM->allocateOutputSpace(prm.featurepv.generative);


  // gtruth is currently spread over 0 and 255 so scaling is needed
  if (prm.unknown == 1)
    bM->scaleGtruthWithUnknown(prm.outscale8, prm.nD);
  else
    bM->scaleGtruth(prm.outscale8);


  bM->checkSanity();

  // features (need to check consistency)
  int numP = bM->getNumPats();
  bM->readFeatureFilesDirs(prm.timeseq, prm.nD, prm.iopv.testDir, numP, bM->getVolumeNames());


  if (prm.timeseq==0)
  {
	  std::vector <std::vector <CImage2> > *inp = bM->getInputImageDirRef();
	  std::vector <std::vector <CByteImage> > *gtr = bM->getGtImageDirRef();
	  bM->preComputeLocalFeatures(inp, gtr ,prm.iopv.testDirIndexV, prm.nD);
	  bM->preComputePairwiseFeatures(inp);

  } else if (prm.timeseq==1)
  {
	  // read comments in files regarding shrinking image into half for computational reasons
	  std::vector <std::vector <CByteImage> > *inp = bM->getInputImagetDirRef();
	  std::vector <std::vector <CByteImage> > *gtr = bM->getGtImageDirRef();
	  bM->preComputeLocalFeatures(inp, gtr ,prm.iopv.testDirIndexV, prm.nD);


	  if (prm.featurepv.opflowlocal == 1)
	  {
	    if (inp->size() == 1)
	    {}
	    else if (inp->size() == 2 && prm.crfpv.crfp == 1 && prm.iopv.testDir == 1)
	    {}
	    else
	    {
	      std::cout << " This feature is supported only when one folder exists ";
	      assert(false);
	    }
	  }

	  if (prm.featurepv.opflowconnection == 1 || prm.featurepv.opflowreverseconnection == 1 || prm.featurepv.opflowlocal == 1) // flowlocal not working for crf
	  {
	    std::vector <string> indirpath = bM->getInputImagetDirPath();
	    bM->extractOpticalFlow(indirpath);
	    bM->preComputePairwiseFeaturesAndPutInGlobalSpace(inp);
	  } else
	  {
	    bM->preComputePairwiseFeatures(inp);
	  }
  }

  bM->preComputeClassRatios();
  bM->preComputeInverseClassRatios();

  // should be called only after images are loaded
  bM->allocateDataCostSpace(); // assumes all images have same width and height
  if (prm.featurepv.opflowconnection == 1 || prm.featurepv.opflowreverseconnection == 1)
  {
    // do nothing as the smooth cost space was allocated in the preComputePairwise functions
    // this was done to avoid having an extra copy of dirGradient and instead only have the
    // global structure that will be needed during inference.
  } else
  {
    bM->allocateSmoothCostGlobalSpace();
  }

  time(&end2);
  timediff = difftime(end2, start);
  LOG(bM->debugfile, " Time : For loading and initialization : %f \n", timediff);


  // logreg output for initial parameters

  time(&start);

  bM->firstInference(0, "out-WTA-");

  if (prm.featurepv.generative == 1)
  {
	  bM->firstInference(1, "out-WTA-Gen-");
  }


  time(&end2);
  timediff = difftime(end2, start);
  LOG(bM->debugfile, " Time : For first inference : %f \n", timediff);

  std::cout << " \n\nStarting ... \n\n";

  // optimize (for train data) or // inference or average-inference (for test data)
  bM->searchParameters();

  std::cout << " Output evaluations! \n";
  bM->outputFinalEvalutions();

  // dump parameters
  bM->dumpThetaParameters(1, 0);

  // remember to free any data structures created earlier and close files
  bM->deAllocatePairwiseAllocations();
  bM->deAllocateSmoothCostGlobalSpace();
  bM->deAllocateDataCostSpace();

  bM->deleteFeatureAllocations(prm.nD, prm.timeseq);

  bM->deAllocateIOfiles();

  delete bM;

  fclose(debugfile);

  return 0;
}
