#include "Parameters.h"
#include <sstream>


Parameters::Parameters()
{

  nD = 2; // number of class labels
  outscale8 = -1;         // scale factor for disparities; -1 means full range 255.0/(nD-1)
  outscale16 = -1;        // (2^16-1)/(nD-1) 
  verbose = 1;               // print messages to stderr

  ignoreVal = -1;        // Value of ground truth disp image to ignore - (vestigeal)

  parallelize = 0;    // whether you want to run Par versions (openMP for klr code)
  timeseq = 0;

  interactive = 0;

  optClassAccuracy = 0;

}



crfparams::crfparams()
{

  crfp = 1;

  hidden = 0;

  // USE_GC 0,  USE_BP 1,  USE_SBP 2,  USE_MF 3,  USE_SMF 4, USE_TRWS 5,  USE_BPS 6,  USE_MP 7

  inferencer = USE_GC;
  inferencerTest = USE_GC; 
  random = 0;			   // use to enable random in graphcut

  damper = 1;
  msgTol = 0.001;
//  MeanField::DiffType msgDiffFun = MeanField::ENERGY;

  inferouteriter = 1;
  inferinneriter = 10;

}


logregparams::logregparams()
{

  logreg = 0;
  logregpair = 0;
  logregpaird = 0;
  logregpaird4 = 0;
  logregpl = 0;

}


ioparams::ioparams()
{

  testDir = 1;  // to register that test case exists
//  std::vector<int> testDirIndexV;

}



klrparams::klrparams()
{
  // klr-exp 0, ivm 1, svm 3
  klrtype = 0; // this will indicate whether the strong classifier is klr-exp, ivm, svm, etc

//  char klrfName[100];

  appKlr = 0;
  hogKlr = 0;
  motKlr = 0;
  rbmKlr = 0;
  intensityKlr = 0;
  locKlr = 0;

  // parameters to control type of model
  klr = 0;

  // ivm parameters
  hogdimklr = 0;
  motdimklr = 0;
  rbmdimklr = 0;
  appdimklr = 0;
  locdimklr = 0;
  biasklr = 0;
  lambda = 0;
//  std::string xfile;
//  std::string wfile;
//  std::string mmfile;
//  std::string kernel;
  p1 = 0;
  ker = RBF;

  mapKernelType["RBF"] = RBF;
  mapKernelType["CHI2"] = CHI2;
  mapKernelType["POLY"] = POLY;
  mapKernelType["LINEAR"] = LINEAR;
  mapKernelType["EUCLIDEAN"] = EUCLIDEAN;
  mapKernelType["NONE"] = NONE;
  mapKernelType[""] = NONE;

  nW = 0;

}




featureparams::featureparams()
{

  // common

  generative = 0;
  featureCode = 0;  
  featureCodeKlr = 0;
  gaussSigma = 0.0;      // Scale of gaussian regularizer (zero for none) - crf parameters only (not klr params)
  bias = 0;

  nF = 0;				 // Number of scaling factor. (vestigeal)

//  char fileName[100];

  learngradonly = 0;        // computeDataCost is not really rerun to save time
  learncrfonly = 0;             // really want this 


  // smoothness

  gradOnly = 0;  // to control whether to use non-pairwise or pairwise gradients or only gradients
  context = 0;
  gradContext = 0;
  gradVZ = 0;
  csgrad = 0;
  nV = 0;
  nC = 0;
  nT = 0;
  nZ = 0;

  // volume
  nE = 0;
  volume = 0;

  // intensity

  intensity = 0;
  nbinsintNB = 50;
  rangeI = 256;
  nU = 0;
  onlyUV = 0;


  // location

  loc = 0; // to control whether location should be used in the generative-discriminative setting
  locCubeSize = 0; // 4*4*4;
  loccubex = 0;
  loccubey = 0;
  loccubez = 0;
  skipXYZlocklr = 0;
  loc_params_movie = 0;  // vestigeal
//  char *locdirname;
  nL = 0;

  // appearance

  app = 0;
//  char appfilename[100];
  nA = 0;

  // hog
  hog = 0;
//  char *hogdirname;
  nH = 0;


  // motion
  mot = 0;
  nM = 0;

  opflowconnection = 0;
  opflowreverseconnection = 0;
  opflowlocal = 0;
  nO = 0;
  opflowlocalframes.clear();
  useopflowcache = 0;

  // rbm
  rbm = 0;


}


graddescparams::graddescparams()
{

  rate = 0.001f;
  maxiter = 100;          // maximum number of iterations for gradient descent
  closeEnoughPercent = 2.0f; // when to stop GC or BP iterations

  // flags for optimization algorithm (choose either gradient descent or bfgs
  // for now, bfgs works only with logreg
  grad_desc = 1;
  bfgs_flag = 0;

  sgdskipV = 0;
  sgdlast = 0;


  bfgs_outer_iter = 20;
  bfgs_inner_iter = 100;

  // when using bfgs, set iterations provided by user to bfgs_outer_iter

  // need only maxiter or outer_iter
  // need only grad_desc or bfgs_flag

  maxiter_out = 1;


  m_bfgs = 6; //limited memory version, not yet implemented
  beta = 0.75;
  beta_dash = 0.25;
  alpha = 1;
  stepmaxconst = 1;
}


int Parameters::readfromcfg(Config &cfg, const char* cfgfname)
{

  try
  {
    cfg.readFile(cfgfname);
  }
  catch(const FileIOException &fioex)
  {
    std::cerr << "I/O error while reading file." << std::endl;
    return(EXIT_FAILURE);
  }
  catch(const ParseException &pex)
  {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    return(EXIT_FAILURE);
  }

  // Get the common parameter values
  try
  {
    cfg.lookupValue("nD", nD);
	  
    cfg.lookupValue("outscale8", outscale8);
    cfg.lookupValue("outscale16", outscale16);
    cfg.lookupValue("verbose", verbose);
    cfg.lookupValue("ignoreVal", ignoreVal);
    cfg.lookupValue("parallelize", parallelize);
    cfg.lookupValue("timeseq", timeseq);
    cfg.lookupValue("interactive", interactive);
    cfg.lookupValue("optClassAccuracy", optClassAccuracy);
    cfg.lookupValue("unknown", unknown);

    cfg.lookupValue("crfparams.crfp", crfpv.crfp);
    cfg.lookupValue("crfparams.hidden", crfpv.hidden);
    cfg.lookupValue("crfparams.inferencer", crfpv.inferencer);
    cfg.lookupValue("crfparams.inferencerTest", crfpv.inferencerTest); 
    cfg.lookupValue("crfparams.random", crfpv.random);			   
		  
    cfg.lookupValue("crfparams.damper", crfpv.damper);
    cfg.lookupValue("crfparams.msgTol", crfpv.msgTol);
		  //  MeanField::DiffType msgDiffFun = MeanField::ENERGY");
		  
    cfg.lookupValue("crfparams.inferouteriter", crfpv.inferouteriter);
    cfg.lookupValue("crfparams.inferinneriter", crfpv.inferinneriter);
		  


    cfg.lookupValue("logregparams.logreg", logregpv.logreg);
    cfg.lookupValue("logregparams.logregpair", logregpv.logregpair);
    cfg.lookupValue("logregparams.logregpaird", logregpv.logregpaird);
    cfg.lookupValue("logregparams.logregpaird4", logregpv.logregpaird4);
    cfg.lookupValue("logregparams.logregpl", logregpv.logregpl );

    cfg.lookupValue("ioparams.indirname", iopv.indirname);
    cfg.lookupValue("ioparams.gtdirname", iopv.gtdirname);
    cfg.lookupValue("ioparams.interdirname", iopv.interdirname);
    cfg.lookupValue("ioparams.outdirname", iopv.outdirname);
    cfg.lookupValue("ioparams.outstem", iopv.outstem);
		  
		  
    cfg.lookupValue("ioparams.testDir", iopv.testDir);

		  // klr-exp 0, ivm 1, svm 3
    cfg.lookupValue("klrparams.klrtype", featurepv.klrpv.klrtype); // this will indicate whether the strong classifier is klr-exp, ivm, svm, etc
		  
    cfg.lookupValue("klrparams.klrfName", featurepv.klrpv.klrfName);
		  
    cfg.lookupValue("klrparams.appKlr", featurepv.klrpv.appKlr);
    cfg.lookupValue("klrparams.hogKlr", featurepv.klrpv.hogKlr);
    cfg.lookupValue("klrparams.motKlr", featurepv.klrpv.motKlr);
    cfg.lookupValue("klrparams.rbmKlr", featurepv.klrpv.rbmKlr);
    cfg.lookupValue("klrparams.intensityKlr", featurepv.klrpv.intensityKlr);
    cfg.lookupValue("klrparams.locKlr", featurepv.klrpv.locKlr);
		  
		  // parameters to control type of model
    cfg.lookupValue("klrparams.klr", featurepv.klrpv.klr);
		  
		  // ivm parameters
    cfg.lookupValue("klrparams.hogdimklr", featurepv.klrpv.hogdimklr);
    cfg.lookupValue("klrparams.motdimklr", featurepv.klrpv.motdimklr);
    cfg.lookupValue("klrparams.rbmdimklr", featurepv.klrpv.rbmdimklr);
    cfg.lookupValue("klrparams.appdimklr", featurepv.klrpv.appdimklr);
    cfg.lookupValue("klrparams.locdimklr", featurepv.klrpv.locdimklr);
    cfg.lookupValue("klrparams.biasklr", featurepv.klrpv.biasklr);
    cfg.lookupValue("klrparams.lambda", featurepv.klrpv.lambda);
    cfg.lookupValue("klrparams.xfile", featurepv.klrpv.xfile);
    cfg.lookupValue("klrparams.wfile", featurepv.klrpv.wfile);
    cfg.lookupValue("klrparams.mmfile", featurepv.klrpv.mmfile);
    cfg.lookupValue("klrparams.kernel", featurepv.klrpv.kernel);
    cfg.lookupValue("klrparams.p1", featurepv.klrpv.p1);
//need to fix this
    std::string tempkernType;
    cfg.lookupValue("klrparams.ker", tempkernType);
    featurepv.klrpv.ker = featurepv.klrpv.mapKernelType[tempkernType];
    cfg.lookupValue("klrparams.appklrdescdirname", featurepv.klrpv.appklrdescdirname);
    cfg.lookupValue("klrparams.hogklrdescdirname", featurepv.klrpv.hogklrdescdirname);


    cfg.lookupValue("featureparams.generative", featurepv.generative);
    cfg.lookupValue("featureparams.featureCode", featurepv.featureCode); 
    cfg.lookupValue("featureparams.featureCodeKlr", featurepv.featureCodeKlr);
    cfg.lookupValue("featureparams.gaussSigma", featurepv.gaussSigma);
    cfg.lookupValue("featureparams.bias", featurepv.bias);
		  
    cfg.lookupValue("featureparams.nF", featurepv.nF);
		  //  	factorU"); 
		  
    cfg.lookupValue("featureparams.fileName", featurepv.fileName);
		  
    cfg.lookupValue("featureparams.learngradonly", featurepv.learngradonly);
    cfg.lookupValue("featureparams.learncrfonly", featurepv.learncrfonly);
		  
		  // smoothness
		  
    cfg.lookupValue("featureparams.gradOnly", featurepv.gradOnly);  // to control whether to use non-pairwise or pairwise gradients or only gradients
    cfg.lookupValue("featureparams.context", featurepv.context);
    cfg.lookupValue("featureparams.gradContext", featurepv.gradContext);
    cfg.lookupValue("featureparams.gradVZ", featurepv.gradVZ);
		  //  	gradThreshVec");
		  //  	thetaV"); // smoothness weights associated with gradient bins
		  //  	thetaC"); // context weights
    cfg.lookupValue("featureparams.nV", featurepv.nV);
    cfg.lookupValue("featureparams.nZ", featurepv.nZ);
    cfg.lookupValue("featureparams.nC", featurepv.nC);
    cfg.lookupValue("featureparams.nT", featurepv.nT);
    cfg.lookupValue("featureparams.nTz", featurepv.nTz);

    // volume
    cfg.lookupValue("featureparams.nE", featurepv.nE);
    cfg.lookupValue("featureparams.volume", featurepv.volume);
		  
		  // intensity
		  
		  //  	thetaU"); // data weights (at most one for now)
    cfg.lookupValue("featureparams.nU", featurepv.nU);
    cfg.lookupValue("featureparams.intensity", featurepv.intensity);
    cfg.lookupValue("featureparams.nbinsintNB", featurepv.nbinsintNB);
    cfg.lookupValue("featureparams.rangeI", featurepv.rangeI);
    cfg.lookupValue("featureparams.onlyUV", featurepv.onlyUV);
		  
    cfg.lookupValue("featureparams.intfilename", featurepv.intfilename);
    cfg.lookupValue("featureparams.intdirname", featurepv.intdirname);
		  
		  
		  // location
		  
    cfg.lookupValue("featureparams.loc", featurepv.loc); // to control whether location should be used in the generative-discriminative setting
	cfg.lookupValue("featureparams.locCubeSize",featurepv.locCubeSize);

	cfg.lookupValue("featureparams.loccubex",featurepv.loccubex);
	cfg.lookupValue("featureparams.loccubey",featurepv.loccubey);
	cfg.lookupValue("featureparams.loccubez",featurepv.loccubez);
	cfg.lookupValue("featureparams.skipXYZlocklr",featurepv.skipXYZlocklr);

	cfg.lookupValue("featureparams.loc_params_movie", featurepv.loc_params_movie);  // vestigeal
		  //  	thetaL",); //location weights
  cfg.lookupValue("featureparams.nL", featurepv.nL);
		  
	cfg.lookupValue("featureparams.locfilename", featurepv.locfilename);  // for reading in parameters - this could be a directory too
	cfg.lookupValue("featureparams.locdirname", featurepv.locdirname);
	cfg.lookupValue("featureparams.imageSliceFile", featurepv.imageSliceFile);
		  		  
		  // appearance
		  
	cfg.lookupValue("featureparams.app", featurepv.app);
		  //  	thetaA",); // appearance weights
	cfg.lookupValue("featureparams.nA", featurepv.nA);
		  
	cfg.lookupValue("featureparams.appfilename", featurepv.appfilename);
	cfg.lookupValue("featureparams.appdirname", featurepv.appdirname);
		  
		  
		  // hog
	cfg.lookupValue("featureparams.hog", featurepv.hog);
  cfg.lookupValue("featureparams.nH", featurepv.nH);
	cfg.lookupValue("featureparams.hogfilename", featurepv.hogfilename);
	cfg.lookupValue("featureparams.hogdirname", featurepv.hogdirname);


	// motion
	cfg.lookupValue("featureparams.mot", featurepv.mot);
  cfg.lookupValue("featureparams.nM", featurepv.nM);
	cfg.lookupValue("featureparams.motfilename", featurepv.motfilename);
	cfg.lookupValue("featureparams.motdirname", featurepv.motdirname);

	// optical flow
	cfg.lookupValue("featureparams.opflowconnection", featurepv.opflowconnection);
	cfg.lookupValue("featureparams.opflowreverseconnection", featurepv.opflowreverseconnection);
	cfg.lookupValue("featureparams.opflowlocal", featurepv.opflowlocal);

	if (featurepv.opflowlocal == 1)
	{
    const Setting& root = cfg.getRoot();
    try
    {
      const Setting &frames = root["featureparams"]["frames"];
      int count = frames.getLength();
      for(int i = 0; i < count; ++i)
      {
        const Setting &frame = frames[i];
        int number;

        if(!frame.lookupValue("num", number))
          continue;
        featurepv.opflowlocalframes.push_back(number);
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
        // Ignore.
    }

    if (featurepv.opflowlocalframes.size() == 0)
    {
      std::cout << " If opflowlocal is enabled at least one frame number should be provided. " <<std::endl;
      assert(false);
    }
	}

	cfg.lookupValue("featureparams.useopflowcache", featurepv.useopflowcache);
	cfg.lookupValue("featureparams.flowcachefile", featurepv.flowcachefile);
	cfg.lookupValue("featureparams.rflowcachefile", featurepv.rflowcachefile);

	cfg.lookupValue("graddescparams.rate", graddescpv.rate);  
	cfg.lookupValue("graddescparams.maxiter", graddescpv.maxiter);          // maximum number of iterations for gradient descent
	cfg.lookupValue("graddescparams.closeEnoughPercent", graddescpv.closeEnoughPercent); // when to stop GC or BP iterations
		  
	cfg.lookupValue("graddescparams.grad_desc", graddescpv.grad_desc);
	cfg.lookupValue("graddescparams.bfgs_flag", graddescpv.bfgs_flag);

	cfg.lookupValue("graddescparams.sgdskipV", graddescpv.sgdskipV);
	cfg.lookupValue("graddescparams.sgdlast", graddescpv.sgdlast);

		  
	cfg.lookupValue("graddescparams.bfgs_outer_iter", graddescpv.bfgs_outer_iter);
	cfg.lookupValue("graddescparams.bfgs_inner_iter", graddescpv.bfgs_inner_iter);
		  
	cfg.lookupValue("graddescparams.maxiter_out", graddescpv.maxiter_out);
		  
	cfg.lookupValue("graddescparams.m_bfgs", graddescpv.m_bfgs); 
	cfg.lookupValue("graddescparams.beta", graddescpv.beta);
	cfg.lookupValue("graddescparams.beta_dash", graddescpv.beta_dash);
	cfg.lookupValue("graddescparams.alpha", graddescpv.alpha);
	cfg.lookupValue("graddescparams.stepmaxconst", graddescpv.stepmaxconst);
	  

  }
  catch(const SettingNotFoundException &nfex)
  {
    std::cerr << "Variable/setting  missing in configuration file." << std::endl;
  }

}



int Parameters::printParameters()
{

  std::cout << " Parameters : \n";
  std::cout<<"\n General  : "<< std::endl;
  std::cout<<"nD  : "<< nD << std::endl;
  std::cout<<"outscale8  : "<< outscale8 << std::endl;
  std::cout<<"outscale16  : "<< outscale16 << std::endl;
  std::cout<<"verbose  : "<< verbose << std::endl;
  std::cout<<"ignoreVal  : "<< ignoreVal << std::endl;
  std::cout<<"parallelize  : "<< parallelize << std::endl;
  std::cout<<"time sequence : "<< timeseq << std::endl;
  std::cout<<"interactive  : "<< interactive << std::endl;
  std::cout<<"optClassAccuracy  : "<< optClassAccuracy << std::endl;
  std::cout<<"unknown  : "<< unknown << std::endl;

  std::cout<<"\n CRF  : "<< std::endl;
  std::cout<<"crfparams.crfp  : "<< crfpv.crfp << std::endl;
  std::cout<<"crfparams.hidden  : "<< crfpv.hidden << std::endl;
  std::cout<<"crfparams.inferencer  : "<< crfpv.inferencer << std::endl;
  std::cout<<"crfparams.inferencerTest  : "<< crfpv.inferencerTest << std::endl; 
  std::cout<<"crfparams.random  : "<< crfpv.random << std::endl;
  std::cout<<"crfparams.damper  : "<< crfpv.damper << std::endl;
  std::cout<<"crfparams.msgTol  : "<< crfpv.msgTol << std::endl;
//  MeanField::DiffType msgDiffFun = MeanField::ENERGY" << std::endl;
  std::cout<<"crfparams.inferouteriter  : "<< crfpv.inferouteriter << std::endl;
  std::cout<<"crfparams.inferinneriter  : "<< crfpv.inferinneriter << std::endl;

  std::cout<<"\nlogreg type  : "<< std::endl;
  std::cout<<"logregparams.logreg  : "<< logregpv.logreg << std::endl;
  std::cout<<"logregparams.logregpair  : "<< logregpv.logregpair << std::endl;
  std::cout<<"logregparams.logregpaird  : "<< logregpv.logregpaird << std::endl;
  std::cout<<"logregparams.logregpaird4  : "<< logregpv.logregpaird4 << std::endl;
  std::cout<<"logregparams.logregpl  : "<< logregpv.logregpl  << std::endl;

  std::cout<<"\nIO  : "<< std::endl;
  std::cout<<"ioparams.indirname  : "<< iopv.indirname << std::endl;
  std::cout<<"ioparams.gtdirname  : "<< iopv.gtdirname << std::endl;
  std::cout<<"ioparams.interdirname  : "<< iopv.interdirname << std::endl;
  std::cout<<"ioparams.outdirname  : "<< iopv.outdirname << std::endl;
  std::cout<<"ioparams.outstem  : "<< iopv.outstem << std::endl;
		  

		  
  std::cout<<"ioparams.testDir  : "<< iopv.testDir << std::endl;

  std::cout<<"\nklr  : "<< std::endl;
// klr-exp 0, ivm 1, svm 3
  std::cout<<"klrparams.klrtype  : "<< featurepv.klrpv.klrtype << std::endl; // this will indicate whether the strong classifier is klr-exp, ivm, svm, etc
		  
  std::cout<<"klrparams.klrfName  : "<< featurepv.klrpv.klrfName << std::endl;
		  
  std::cout<<"klrparams.appKlr  : "<< featurepv.klrpv.appKlr << std::endl;
  std::cout<<"klrparams.hogKlr  : "<< featurepv.klrpv.hogKlr << std::endl;
  std::cout<<"klrparams.motKlr  : "<< featurepv.klrpv.motKlr << std::endl;
  std::cout<<"klrparams.rbmKlr  : "<< featurepv.klrpv.rbmKlr << std::endl;
  std::cout<<"klrparams.intensityKlr  : "<< featurepv.klrpv.intensityKlr << std::endl;
  std::cout<<"klrparams.locKlr  : "<< featurepv.klrpv.locKlr << std::endl;
		  
// parameters to control type of model
  std::cout<<"klrparams.klr  : "<< featurepv.klrpv.klr << std::endl;
		  
// ivm parameters
  std::cout<<"klrparams.hogdimklr  : "<< featurepv.klrpv.hogdimklr << std::endl;
  std::cout<<"klrparams.motdimklr  : "<< featurepv.klrpv.motdimklr << std::endl;
  std::cout<<"klrparams.rbmdimklr  : "<< featurepv.klrpv.rbmdimklr << std::endl;
  std::cout<<"klrparams.appdimklr  : "<< featurepv.klrpv.appdimklr << std::endl;
  std::cout<<"klrparams.locdimklr  : "<< featurepv.klrpv.locdimklr << std::endl;
  std::cout<<"klrparams.biasklr  : "<< featurepv.klrpv.biasklr << std::endl;
  std::cout<<"klrparams.lambda  : "<< featurepv.klrpv.lambda << std::endl;
  std::cout<<"klrparams.xfile  : "<< featurepv.klrpv.xfile << std::endl;
  std::cout<<"klrparams.wfile  : "<< featurepv.klrpv.wfile << std::endl;
  std::cout<<"klrparams.mmfile  : "<< featurepv.klrpv.mmfile << std::endl;
  std::cout<<"klrparams.kernel  : "<< featurepv.klrpv.kernel << std::endl;
  std::cout<<"klrparams.p1  : "<< featurepv.klrpv.p1 << std::endl;
//need to fix this
  std::cout<<"klrparams.ker  : "<< (int)featurepv.klrpv.ker << std::endl;
  std::cout<<"klrparams.hogklrdescdirname  : "<< featurepv.klrpv.hogklrdescdirname << std::endl;
  std::cout<<"klrparams.appklrdescdirname  : "<< featurepv.klrpv.appklrdescdirname << std::endl;

  std::cout<<"\nfeatures  : "<< std::endl;
  std::cout<<"featureparams.generative  : "<< featurepv.generative << std::endl;
  std::cout<<"featureparams.featureCode  : "<< featurepv.featureCode << std::endl; 
  std::cout<<"featureparams.featureCodeKlr  : "<< featurepv.featureCodeKlr << std::endl;
  std::cout<<"featureparams.gaussSigma  : "<< featurepv.gaussSigma << std::endl;
  std::cout<<"featureparams.bias  : "<< featurepv.bias << std::endl;
		  
  std::cout<<"featureparams.nF  : "<< featurepv.nF << std::endl;
//  	factorU" << std::endl; 
		  
  std::cout<<"featureparams.fileName  : "<< featurepv.fileName << std::endl;
		  
  std::cout<<"featureparams.learngradonly  : "<< featurepv.learngradonly << std::endl;
  std::cout<<"featureparams.learncrfonly  : "<< featurepv.learncrfonly << std::endl;
		  
// smoothness
		  
  std::cout<<"featureparams.gradOnly  : "<< featurepv.gradOnly << std::endl;  // to control whether to use non-pairwise or pairwise gradients or only gradients
  std::cout<<"featureparams.context  : "<< featurepv.context << std::endl;
  std::cout<<"featureparams.gradContext  : "<< featurepv.gradContext << std::endl;
  std::cout<<"featureparams.gradVZ  : "<< featurepv.gradVZ << std::endl;
  std::cout<<"featureparams.csgrad  : "<< featurepv.csgrad << std::endl;
//  	gradThreshVec" << std::endl;
//  	thetaV" << std::endl; // smoothness weights associated with gradient bins
//  	thetaC" << std::endl; // context weights
//  std::cout<<"featureparams.nV  : "<< featurepv.nV << std::endl;
//  std::cout<<"featureparams.nZ  : "<< featurepv.nZ << std::endl;
//  std::cout<<"featureparams.nC  : "<< featurepv.nC << std::endl;
//  std::cout<<"featureparams.nT  : "<< featurepv.nT << std::endl;
//  std::cout<<"featureparams.nTz  : "<< featurepv.nTz << std::endl;
		  
  // volume
  std::cout<<"featureparams.volume  : "<< featurepv.volume << std::endl;

// intensity
		  
//  	thetaU" << std::endl; // data weights (at most one for now)
//  std::cout<<"featureparams.nU  : "<< featurepv.nU << std::endl;
  std::cout<<"featureparams.intensity  : "<< featurepv.intensity << std::endl;
  std::cout<<"featureparams.nbinsintNB  : "<< featurepv.nbinsintNB << std::endl;
  std::cout<<"featureparams.rangeI  : "<< featurepv.rangeI << std::endl;
  std::cout<<"featureparams.onlyUV  : "<< featurepv.onlyUV << std::endl;
		  
  std::cout<<"featureparams.intfilename  : "<< featurepv.intfilename << std::endl;
  std::cout<<"featureparams.intdirname  : "<< featurepv.intdirname << std::endl;
		  
		  
// location
		  
  std::cout<<"featureparams.loc  : "<< featurepv.loc << std::endl; // to control whether location should be used in the generative-discriminative setting
  std::cout<<"featureparams.locCubeSize  : " << featurepv.locCubeSize << std::endl;

  std::cout<<"featureparams.loccubex  : " << featurepv.loccubex << std::endl;
  std::cout<<"featureparams.loccubey  : " << featurepv.loccubey << std::endl;
  std::cout<<"featureparams.loccubez  : " << featurepv.loccubez << std::endl;

  std::cout<<"featureparams.skipXYZlocklr  : " << featurepv.skipXYZlocklr << std::endl;


  std::cout<<"featureparams.loc_params_movie  : "<< featurepv.loc_params_movie << std::endl;  // vestigeal
//  	thetaL", << std::endl; //location weights
//  std::cout<<"featureparams.nL  : "<< featurepv.nL << std::endl;
		  
  std::cout<<"featureparams.locfilename  : "<< featurepv.locfilename << std::endl;  // for reading in parameters - this could be a directory too
  std::cout<<"featureparams.locdirname  : "<< featurepv.locdirname << std::endl;
  std::cout<<"featureparams.imageSliceFile  : "<< featurepv.imageSliceFile << std::endl;

// appearance
		  
  std::cout<<"featureparams.app  : "<< featurepv.app << std::endl;
//  	thetaA", << std::endl; // appearance weights
//  std::cout<<"featureparams.nA  : "<< featurepv.nA << std::endl;
		  
  std::cout<<"featureparams.appfilename  : "<< featurepv.appfilename << std::endl;
  std::cout<<"featureparams.appdirname  : "<< featurepv.appdirname << std::endl;
		  
		  
// hog
  std::cout<<"featureparams.hog  : "<< featurepv.hog << std::endl;
//  	thetaH", << std::endl; // HOG weights
//  std::cout<<"featureparams.nH  : "<< featurepv.nH << std::endl;
		  
  std::cout<<"featureparams.hogfilename  : "<< featurepv.hogfilename << std::endl;
  std::cout<<"featureparams.hogdirname  : "<< featurepv.hogdirname << std::endl;


// motion
  std::cout<<"featureparams.mot  : "<< featurepv.mot << std::endl;
  std::cout<<"featureparams.motfilename  : "<< featurepv.motfilename << std::endl;
  std::cout<<"featureparams.motdirname  : "<< featurepv.motdirname << std::endl;

  // optical flow
  std::cout<<"featureparams.opflowconnection  : "<< featurepv.opflowconnection << std::endl;
  std::cout<<"featureparams.opflowreverseconnection  : "<< featurepv.opflowreverseconnection << std::endl;
  std::cout<<"featureparams.opflowlocal  : "<< featurepv.opflowlocal << std::endl;
  std::cout<<"featureparams.opflowlocalframes  : " << std::endl;
  for (int i = 0; i < featurepv.opflowlocalframes.size(); i++)
  {
    std::cout << "  " << featurepv.opflowlocalframes[i];
  }
  std::cout<< std::endl;
  std::cout << "featureparams.useopflowcache : " << featurepv.useopflowcache << std::endl;
  std::cout << "featureparams.flowcachefile : " << featurepv.flowcachefile << std::endl;
  std::cout << "featureparams.rflowcachefile : " <<  featurepv.rflowcachefile << std::endl;


  std::cout<<"\ngradient learning  : "<< std::endl;
  std::cout<<"graddescparams.rate  : "<< graddescpv.rate << std::endl;  
  std::cout<<"graddescparams.maxiter  : "<< graddescpv.maxiter << std::endl;          // maximum number of iterations for gradient descent
  std::cout<<"graddescparams.closeEnoughPercent  : "<< graddescpv.closeEnoughPercent << std::endl; // when to stop GC or BP iterations
		  
  std::cout<<"graddescparams.grad_desc  : "<< graddescpv.grad_desc << std::endl;
  std::cout<<"graddescparams.bfgs_flag  : "<< graddescpv.bfgs_flag << std::endl;

  std::cout<<"graddescparams.sgdskipV  : "<< graddescpv.sgdskipV << std::endl;
  std::cout<<"graddescparams.sgdlast  : "<< graddescpv.sgdlast << std::endl;

		  
  std::cout<<"graddescparams.bfgs_outer_iter  : "<< graddescpv.bfgs_outer_iter << std::endl;
  std::cout<<"graddescparams.bfgs_inner_iter  : "<< graddescpv.bfgs_inner_iter << std::endl;
		  
  std::cout<<"graddescparams.maxiter_out  : "<< graddescpv.maxiter_out << std::endl;
		  
  std::cout<<"graddescparams.m_bfgs  : "<< graddescpv.m_bfgs << std::endl; 
  std::cout<<"graddescparams.beta  : "<< graddescpv.beta << std::endl;
  std::cout<<"graddescparams.beta_dash  : "<< graddescpv.beta_dash << std::endl;
  std::cout<<"graddescparams.alpha  : "<< graddescpv.alpha << std::endl;
  std::cout<<"graddescparams.stepmaxconst  : "<< graddescpv.stepmaxconst << std::endl;

}



int Parameters::updateDependParameters()
{

	iopv.outstem = iopv.outdirname;

  featurepv.intensity = (int)featurepv.featureCode/1000000;
  featurepv.hog = (int) (featurepv.featureCode/100000)%10;
  featurepv.app = (int) (featurepv.featureCode/10000)%10;
  featurepv.loc = (int) (featurepv.featureCode/1000)%10;
  featurepv.mot = (int) (featurepv.featureCode/100)%10;
  featurepv.rbm = (int) (featurepv.featureCode/10)%10;

  int tempGr = (int) featurepv.featureCode%10;

  if (tempGr == 0)
    logregpv.logreg = 1;
  else if (tempGr == 1)
    featurepv.gradOnly = 1;  //CRF or KLR
  else if (tempGr == 2)
    featurepv.gradContext = 1; //CRF or KLR
  //	else if (tempGr == 4)
  //	prm.featurepv.gradOnly = 2;   // KLR
  //		else if (tempGr == 5)
  //	prm.featurepv.gradContext = 2; //KLR
  else if (tempGr == 5)
    logregpv.logregpaird4 = 1;
  else if (tempGr == 6)
    featurepv.context = 1; //CRF
  //else if (tempGr == 7)
  //	prm.featurepv.context = 2; //KLR
  else if (tempGr == 7)
    logregpv.logregpl = 1;
  else if (tempGr == 8)  // this is by default prm.featurepv.gradContext
    logregpv.logregpair = 1;
  else if (tempGr == 9)
    logregpv.logregpaird = 1;


  if (featurepv.intensity==3 || featurepv.intensity == 6 || featurepv.app==3 || featurepv.hog==3
      || featurepv.loc==3 || featurepv.loc==4 || featurepv.loc==8 || featurepv.loc==9)
    featurepv.generative = 1;



  featurepv.klrpv.intensityKlr = (int)featurepv.featureCodeKlr/1000000;
  featurepv.klrpv.hogKlr = (int) (featurepv.featureCodeKlr/100000)%10;
  featurepv.klrpv.appKlr = (int) (featurepv.featureCodeKlr/10000)%10;
  featurepv.klrpv.locKlr = (int) (featurepv.featureCodeKlr/1000)%10;
  featurepv.klrpv.motKlr = (int) (featurepv.featureCodeKlr/100)%10;
  featurepv.klrpv.rbmKlr = (int) (featurepv.featureCodeKlr/10)%10;

  if (featurepv.klrpv.intensityKlr==2 || featurepv.klrpv.appKlr==2 || featurepv.klrpv.hogKlr==2 || featurepv.klrpv.locKlr==2 || featurepv.klrpv.motKlr==2 || featurepv.klrpv.rbmKlr==2)
    featurepv.klrpv.klr = 1;



  if (featurepv.klrpv.klr==1)
  {
    if (featurepv.klrpv.intensityKlr>0)
      featurepv.intensity=0;
    if (featurepv.klrpv.appKlr>0)
      featurepv.app=0;
    if (featurepv.klrpv.locKlr>0)
      featurepv.loc=0;
    if (featurepv.klrpv.hogKlr>0)
      featurepv.hog=0;
    if (featurepv.klrpv.motKlr>0)
      featurepv.mot=0;
    if (featurepv.klrpv.rbmKlr>0)
      featurepv.rbm=0;
  }


  if (featurepv.useopflowcache == 0)
  {
    featurepv.flowcachefile.clear();
    featurepv.rflowcachefile.clear();
  }

}

// we no longer can read feature value parameters on command line
int Parameters::readfromCommandline(int argc, char **argv)
{

	  // parse command-line arguments using the "getopt" utility
	  int o;
	  while ((o = getopt(argc, argv, "a:b:e:i:k:m:n:o:r:t:w:x:z:")) != -1)
		switch (o) {
		case 'n': nD = atoi(optarg); break;
		case 'r': graddescpv.rate = (float)atof(optarg); break;
		case 'i': graddescpv.maxiter = atoi(optarg); break;
    case 'm': {

      int temp_grtype = atoi(optarg);
      if (temp_grtype==0)
      {
        graddescpv.grad_desc = 1;
        graddescpv.bfgs_flag = 0;
      } else if (temp_grtype==1)
      {
        graddescpv.grad_desc = 0;
        graddescpv.bfgs_flag = 1;
      } else if (temp_grtype==2)
      {
        graddescpv.grad_desc = 0;
        graddescpv.bfgs_flag = 2;
      }

    }
		case 'x': graddescpv.closeEnoughPercent = (float)atof(optarg); break;
		case 'o': outscale8 = atoi(optarg); break;
		case 'e': interactive = atoi(optarg); break;
    case 'z': featurepv.gaussSigma = (float)atof(optarg); break;
		case 'w': iopv.testDir = (int) atoi(optarg); break;  // used if a test directory of images exist
		case 'a': featurepv.fileName = optarg; break;
		case 't': crfpv.inferencer = (int) atoi(optarg); break;
		case 'b': {
	      featurepv.featureCode = (int) atoi(optarg);

	      featurepv.intensity = (int)featurepv.featureCode/1000000;
	      featurepv.hog = (int) (featurepv.featureCode/100000)%10;
	      featurepv.app = (int) (featurepv.featureCode/10000)%10;
	      featurepv.loc = (int) (featurepv.featureCode/1000)%10;
	      featurepv.mot = (int) (featurepv.featureCode/100)%10;
	      featurepv.rbm = (int) (featurepv.featureCode/10)%10;

	      int tempGr = (int) featurepv.featureCode%10;

	      if (tempGr == 0)
	        logregpv.logreg = 1;
	      else if (tempGr == 1)
	        featurepv.gradOnly = 1;  //CRF or KLR
	      else if (tempGr == 2)
	        featurepv.gradContext = 1; //CRF or KLR
	      else if (tempGr == 5)
	        logregpv.logregpaird4 = 1;
	      else if (tempGr == 6)
	        featurepv.context = 1; //CRF
	      else if (tempGr == 7)
	        logregpv.logregpl = 1;
	      else if (tempGr == 8)  // this is by default prm.featurepv.gradContext
	        logregpv.logregpair = 1;
	      else if (tempGr == 9)
	        logregpv.logregpaird = 1;

	      break;
		}
		case 'k': {
	      featurepv.featureCodeKlr = (int) atoi(optarg);

	      featurepv.klrpv.intensityKlr = (int)featurepv.featureCodeKlr/1000000;
	      featurepv.klrpv.hogKlr = (int) (featurepv.featureCodeKlr/100000)%10;
	      featurepv.klrpv.appKlr = (int) (featurepv.featureCodeKlr/10000)%10;
	      featurepv.klrpv.locKlr = (int) (featurepv.featureCodeKlr/1000)%10;
	      featurepv.klrpv.motKlr = (int) (featurepv.featureCodeKlr/100)%10;
	      featurepv.klrpv.rbmKlr = (int) (featurepv.featureCodeKlr/10)%10;

	      if (featurepv.klrpv.intensityKlr==2 || featurepv.klrpv.appKlr==2 || featurepv.klrpv.hogKlr==2 || featurepv.klrpv.locKlr==2 || featurepv.klrpv.motKlr==2 || featurepv.klrpv.rbmKlr==2)
	        featurepv.klrpv.klr = 1;
	      break;
		}
		default:
	      fprintf(stderr, "Ignoring unrecognized option %s\n", argv[optind-1]);
		}





}



int Parameters::readfeaturefile()
{
  std::string line, paramstr;

  std::cout<<featurepv.fileName.c_str()<< std::endl;

  std::ifstream myfile (featurepv.fileName.c_str());
  if (myfile.is_open())
	{
	  while (! myfile.eof() )
		{
		  getline (myfile,line);
		  paramstr.append(" ");
		  paramstr.append(line);

		}
	  std::cout<<" " <<line<<std::endl;
	  myfile.close();
	}
  else {
	std::cout << "Unable to open feature parameter file";
	exit(1);
  }

  std::stringstream ostr(paramstr);
  std::string temp;

  int nump = 0;
  while (ostr >> temp) {
	if (temp.compare("-u")==0) {
	  ostr >> temp;
	  featurepv.thetaU.push_back((float)atof(temp.c_str()));
	} else if (temp.compare("-v")==0) {
	  ostr >> temp;
	  featurepv.thetaV.push_back((float)atof(temp.c_str()));
	} else if (temp.compare("-d")==0) {
		ostr >> temp;
		featurepv.thetaZ.push_back((float)atof(temp.c_str()));
	} else if (temp.compare("-g")==0) {
	  ostr >> temp;
	  featurepv.gradThreshVec.push_back(atoi(temp.c_str()));
	} else if (temp.compare("-y")==0) {
		ostr >> temp;
		featurepv.gradThreshVecZ.push_back(atoi(temp.c_str()));
	} else if (temp.compare("-c")==0) {
	  ostr >> temp;
	  featurepv.thetaC.push_back((float)atof(temp.c_str()));
	} else if (temp.compare("-s")==0) {
	  ostr >> temp;
	  featurepv.thetaA.push_back((float)atof(temp.c_str()));
	} else if (temp.compare("-l")==0) {
	  ostr >> temp;
	  featurepv.thetaL.push_back((float)atof(temp.c_str()));
	} else if (temp.compare("-q")==0) {
        ostr >> temp;
      featurepv.thetaM.push_back((float)atof(temp.c_str()));
	} else if (temp.compare("-f")==0) {
	  ostr >> temp;
	  featurepv.thetaO.push_back((float)atof(temp.c_str()));
	} else if (temp.compare("-j")==0) {
    ostr >> temp;
    featurepv.thetaB.push_back((float)atof(temp.c_str()));
	} else if (temp.compare("-p")==0) {
	    ostr >> temp;
	    featurepv.thetaE.push_back((float)atof(temp.c_str()));
	} else if (temp.compare("-h")==0) {
	  ostr >> temp;
	  featurepv.thetaH.push_back((float)atof(temp.c_str()));
	}else {
	  std::cout<< " Error in arguments -u, -v, -g, -d, -y, -j, -q, -c, -l, -s, -f or -h ::: " << temp.c_str() << " \n";
	  exit(1);
	}
	nump++;
  }

}



int Parameters::updatePropertiesParameters()
{

	featurepv.intensity = (int)featurepv.featureCode/1000000;
	featurepv.hog = (int) (featurepv.featureCode/100000)%10;
	featurepv.app = (int) (featurepv.featureCode/10000)%10;
	featurepv.loc = (int) (featurepv.featureCode/1000)%10;
  featurepv.mot = (int) (featurepv.featureCode/100)%10;
  featurepv.rbm = (int) (featurepv.featureCode/10)%10;

	int tempGr = (int) featurepv.featureCode%10;

	if (tempGr == 0)
		logregpv.logreg = 1;
	else if (tempGr == 1)
		featurepv.gradOnly = 1;  //CRF or KLR
	else if (tempGr == 2)
		featurepv.gradContext = 1; //CRF or KLR
	//	else if (tempGr == 4)
	//	prm.featurepv.gradOnly = 2;   // KLR
	//		else if (tempGr == 5)
	//	prm.featurepv.gradContext = 2; //KLR
	else if (tempGr == 5)
		logregpv.logregpaird4 = 1;
	else if (tempGr == 6)
		featurepv.context = 1; //CRF
	//else if (tempGr == 7)
	//	prm.featurepv.context = 2; //KLR
	else if (tempGr == 7)
		logregpv.logregpl = 1;
	else if (tempGr == 8)  // this is by default prm.featurepv.gradContext
		logregpv.logregpair = 1;
	else if (tempGr == 9)
		logregpv.logregpaird = 1;


	featurepv.klrpv.intensityKlr = (int)featurepv.featureCodeKlr/1000000;
	featurepv.klrpv.hogKlr = (int) (featurepv.featureCodeKlr/100000)%10;
	featurepv.klrpv.appKlr = (int) (featurepv.featureCodeKlr/10000)%10;
	featurepv.klrpv.locKlr = (int) (featurepv.featureCodeKlr/1000)%10;
	featurepv.klrpv.motKlr = (int) (featurepv.featureCodeKlr/100)%10;
	featurepv.klrpv.rbmKlr = (int) (featurepv.featureCodeKlr/10)%10;

	if (featurepv.klrpv.intensityKlr==2 || featurepv.klrpv.appKlr==2 || featurepv.klrpv.hogKlr==2 || featurepv.klrpv.locKlr==2 || featurepv.klrpv.motKlr==2 || featurepv.klrpv.rbmKlr==2)
		featurepv.klrpv.klr = 1;


	iopv.outstem = iopv.outdirname;

	if (featurepv.klrpv.klr==1)
	{
		if (featurepv.klrpv.intensityKlr>0)
		  featurepv.intensity=0;
		if (featurepv.klrpv.appKlr>0)
		  featurepv.app=0;
		if (featurepv.klrpv.locKlr>0)
		  featurepv.loc=0;
		if (featurepv.klrpv.hogKlr>0)
		  featurepv.hog=0;
		if (featurepv.klrpv.motKlr>0)
		  featurepv.mot=0;
		if (featurepv.klrpv.rbmKlr>0)
		  featurepv.rbm=0;

	}

	if (featurepv.intensity==3 || featurepv.intensity == 6 || featurepv.app==3 || featurepv.hog==3 || featurepv.loc==3 || featurepv.loc==8 || featurepv.loc==9)
		featurepv.generative = 1;

	// if the image file is in 8bit format, we should use this variable
	if (outscale8 < 0)
		outscale8 = 255 / (nD-1);

	// if the image file is in 16bit format, we should use this variable
	if (outscale16 < 0)
		outscale16 = 65535 / (nD-1);


	featurepv.gradThreshVec.push_back(INT_MAX);
	featurepv.gradThreshVecZ.push_back(INT_MAX);

   // these lines make sure that n? correspond to the theta values in the theta parameter file

    if (featurepv.intensity != 1)
      featurepv.thetaU.erase(featurepv.thetaU.begin(), featurepv.thetaU.end());

    featurepv.nU = (int)featurepv.thetaU.size(); // number of data parameters

    // number of texton parameters = number of texton centers in parameter file
    if (featurepv.app != 1 && featurepv.app != 2)
      featurepv.thetaA.erase(featurepv.thetaA.begin(), featurepv.thetaA.end());

    featurepv.nA = (int)featurepv.thetaA.size();

    if (featurepv.hog != 1)
      featurepv.thetaH.erase(featurepv.thetaH.begin(), featurepv.thetaH.end());

    featurepv.nH = (int)featurepv.thetaH.size();

    if (featurepv.bias != 1)
      featurepv.thetaB.erase(featurepv.thetaB.begin(), featurepv.thetaB.end());

    featurepv.nB = (int)featurepv.thetaB.size();

    if (featurepv.volume != 1)
      featurepv.thetaE.erase(featurepv.thetaE.begin(), featurepv.thetaE.end());

    featurepv.nE = (int)featurepv.thetaE.size();



    if (featurepv.loc != 1 && featurepv.loc != 4 && featurepv.loc != 5 && featurepv.loc!=6 && featurepv.loc!=8 && featurepv.loc!=7 && featurepv.loc!=9)
      featurepv.thetaL.erase(featurepv.thetaL.begin(), featurepv.thetaL.end());

    featurepv.nL = (int) featurepv.thetaL.size();


    if (featurepv.mot != 1)
      featurepv.thetaM.erase(featurepv.thetaM.begin(), featurepv.thetaM.end());

    featurepv.nM = (int)featurepv.thetaM.size();

    if (featurepv.opflowlocal != 1)
      featurepv.thetaO.erase(featurepv.thetaO.begin(), featurepv.thetaO.end());

    featurepv.nO = (int)featurepv.thetaO.size();


    if (logregpv.logreg==1) {
      featurepv.thetaV.erase(featurepv.thetaV.begin(), featurepv.thetaV.end());
      featurepv.gradThreshVec.erase(featurepv.gradThreshVec.begin(), featurepv.gradThreshVec.end());
      featurepv.thetaZ.erase(featurepv.thetaZ.begin(), featurepv.thetaZ.end());
      featurepv.gradThreshVecZ.erase(featurepv.gradThreshVecZ.begin(), featurepv.gradThreshVecZ.end());
    }

    featurepv.nV = (int)featurepv.thetaV.size(); // number of smoothness parameters
    if (featurepv.nV == 1)
    	featurepv.csgrad = 1;
    featurepv.nZ = (int)featurepv.thetaZ.size();
    featurepv.nG = (int)featurepv.gradThreshVec.size();
    featurepv.nGZ = (int)featurepv.gradThreshVecZ.size();


    if (featurepv.context == 0)
      featurepv.thetaC.erase(featurepv.thetaC.begin(), featurepv.thetaC.end());

    featurepv.nC = (int)featurepv.thetaC.size();

    // bfgs
    if (graddescpv.bfgs_flag != 0)
    {
    	graddescpv.maxiter = graddescpv.bfgs_inner_iter;
    	graddescpv.maxiter_out = graddescpv.bfgs_outer_iter;
    }


}



int Parameters::returnConcatThetaVector(std::vector<float> *theta)
{
	// follow order thetaU, thetaV, thetaZ, thetaC, thetaA, thetaH, thetaL, thetaM, thetaO, thetaB, thetaE

	int total = featurepv.nU + featurepv.nV + featurepv.nZ + featurepv.nC + featurepv.nA + featurepv.nH + featurepv.nL + featurepv.nM + featurepv.nO + featurepv.nB + featurepv.nE;
	(*theta).resize(total);

	int k, i = 0;

	for (k = 0; k < featurepv.nU; k++)
		(*theta)[i++] = featurepv.thetaU[k];
	for (k = 0; k < featurepv.nV; k++)
		(*theta)[i++] = featurepv.thetaV[k];
	for (k = 0; k < featurepv.nZ; k++)
		(*theta)[i++] = featurepv.thetaZ[k];
	for (k = 0; k < featurepv.nC; k++)
		(*theta)[i++] = featurepv.thetaC[k];
	for (k = 0; k < featurepv.nA; k++)
		(*theta)[i++] = featurepv.thetaA[k];
	for (k = 0; k < featurepv.nH; k++)
		(*theta)[i++] = featurepv.thetaH[k];
	for (k = 0; k < featurepv.nL; k++)
		(*theta)[i++] = featurepv.thetaL[k];
	for (k = 0; k < featurepv.nM; k++)
		(*theta)[i++] = featurepv.thetaM[k];
	 for (k = 0; k < featurepv.nO; k++)
	    (*theta)[i++] = featurepv.thetaO[k];
   for (k = 0; k < featurepv.nB; k++)
      (*theta)[i++] = featurepv.thetaB[k];
   for (k = 0; k < featurepv.nE; k++)
      (*theta)[i++] = featurepv.thetaE[k];

}




int Parameters::splitVectorToOriginal(std::vector<float> theta)
{
	// follow order thetaU, thetaV, thetaZ, thetaC, thetaA, thetaH, thetaL, thetaM, thetaB, thetaE

	int total = featurepv.nU + featurepv.nV + featurepv.nZ + featurepv.nC + featurepv.nA + featurepv.nH + featurepv.nL + featurepv.nM + featurepv.nO + featurepv.nB + featurepv.nE;

	int k, i = 0;

	for (k = 0; k < featurepv.nU; k++)
		 featurepv.thetaU[k] = theta[i++];
	for (k = 0; k < featurepv.nV; k++)
		featurepv.thetaV[k] = theta[i++];
	for (k = 0; k < featurepv.nZ; k++)
		featurepv.thetaZ[k] = theta[i++];
	for (k = 0; k < featurepv.nC; k++)
		featurepv.thetaC[k] = theta[i++];
	for (k = 0; k < featurepv.nA; k++)
		featurepv.thetaA[k] = theta[i++];
	for (k = 0; k < featurepv.nH; k++)
		featurepv.thetaH[k] = theta[i++];
	for (k = 0; k < featurepv.nL; k++)
		featurepv.thetaL[k] = theta[i++];
	for (k = 0; k < featurepv.nM; k++)
		featurepv.thetaM[k] = theta[i++];
	for (k = 0; k < featurepv.nO; k++)
	  featurepv.thetaO[k] = theta[i++];
	for (k = 0; k < featurepv.nB; k++)
	  featurepv.thetaB[k] = theta[i++];
  for (k = 0; k < featurepv.nE; k++)
    featurepv.thetaE[k] = theta[i++];

}
