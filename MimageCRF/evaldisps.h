
// computes and prints disparity error statistics, creates errormap
// returns percentages of bad pixels (error > 1) and RMS disp error (both in nonoccluded areas)
// verbose==1: print lots, verbose==0: print single line, verbose<0: print nothing

//void evaldisps(CByteImage disp, CByteImage truedisp, CByteImage &errormap,
//	       float &bad1, float &rms, int verbose);

void evaldisps(std::vector<CByteImage> disp, std::vector<CByteImage> truedisp, std::vector<CByteImage> &errormap,
	       float &errpercent, int nD, double &cnterrfull, double &totfull, double **confMatrixFull, double *totNDFull, int paramlreg);

void evaldisps(std::vector<CByteImage> disp, std::vector<CByteImage> truedisp, std::vector<CByteImage> &errormap,
	       float &errpercent, int nD, int testcase, int paramlreg, std::vector<CByteImage> interImage, int interactive);

void evaldisps(std::vector<CByteImage> disp, std::vector<CByteImage> truedisp, std::vector<CByteImage> &errormap,
               float &errpercent, int nD, int testcase, int paramlreg, std::vector<CByteImage> interImage, int interactive, double &cnterrfull, double &totfull, double **confMatrixFull, double *totNDFull);


void userlabelfix(std::vector<CByteImage> &disp, std::vector<CByteImage> truedisp, std::vector<CByteImage> &errormap, std::vector<CByteImage> interImage);

void confusionMatrix(CByteImage truedisp, CByteImage disp);

void evaldispsSet(std::vector<CByteImage> disp, std::vector<CByteImage> truedisp, std::vector<CByteImage> &errormap,
               float &bad1, float &rms, int verbose);

void confusionMatrix(std::vector<CByteImage> &truedisp, std::vector<CByteImage> &disp);

void display_CM_avg(double **confMatrixFull, int nD, double* totNDFull, int train);

