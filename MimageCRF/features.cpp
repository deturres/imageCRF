/*
 * features.cpp
 *
 *  Created on: Jul 16, 2010
 *      Author: bhole
 */

#include "features.h"

// code modified from http://www.rawness.es/cielab/?lang=en
// To the extent possible under law, Manuel Llorens <manuelllorens@gmail.com>
// has waived all copyright and related or neighboring rights to this work.
// This code is licensed under CC0 v1.0, see license information at
// http://creativecommons.org/publicdomain/zero/1.0/

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define LIM(x,min,max) MAX(min,MIN(x,max))
#define CLIP(x) LIM(x,0,65535)

double  ep=8.85645167903563110267661784291703952476382255554199e-03; //   ep=216.0/24389.0;
double  ka=9.03296296296296304717543534934520721435546875000000e+02; //   ka=24389.0/27.0;

float d50_white[3]={0.964220,1,0.825211};
static const double rgb_xyz[3][3] =
       { { 0.7976748465, 0.1351917082, 0.0313534088 },
         { 0.2880402025, 0.7118741325, 0.0000856651 },
         { 0.0000000000, 0.0000000000, 0.8252114389 } }; // From Jacques Desmis

// using pseudoinverse function below
static const double xyz_rgb_2[3][3] =
			{{ 1.3459434662015494765796575e+00, -5.4459884249024337332656387e-01, -6.9388939039072283776476979e-18 },
			{-2.5560754075639424698351831e-01,  1.5081672430343562307797356e+00, 3.9573379295720911841272027e-18},
			{-5.1111772187973379677483621e-02, 2.0535140503506361941976621e-02, 1.2118106376869808293861297e+00}};

// using matlab pinv (this is the transposed matrix of that obtained by the pseudoinverse function below)
//static const double xyz_rgb_2[3][3] =
//			{{ 1.3459434662015494765796575e+00, -2.5560754075639424698351831e-01 , -5.1111772187973379677483621e-02},
//			{-5.4459884249024337332656387e-01,  1.5081672430343562307797356e+00,  2.0535140503506361941976621e-02},
//			{-6.9388939039072283776476979e-18 , 3.9573379295720911841272027e-18, 1.2118106376869808293861297e+00}};


/*
void pseudoinverse (double (*in)[3], double (*out)[3], int size)
{
   double work[3][6], num;
   int i, j, k;

   for (i=0; i < 3; i++) {
     for (j=0; j < 6; j++)
       work[i][j] = j == i+3;
     for (j=0; j < 3; j++)
       for (k=0; k < size; k++)
    	   work[i][j] += in[k][i] * in[k][j];
   }

   for (i=0; i < 3; i++) {
     num = work[i][i];
     for (j=0; j < 6; j++)
       work[i][j] /= num;
     for (k=0; k < 3; k++) {
       if (k==i) continue;
       num = work[k][i];
       for (j=0; j < 6; j++)
    	   work[k][j] -= work[i][j] * num;
     }
   }

   for (i=0; i < size; i++)
     for (j=0; j < 3; j++)
       for (out[i][j]=k=0; k < 3; k++)
    	   out[i][j] += work[j][k+3] * in[i][k];
 }
*/


double f_cbrt(double r)
{
     r/=65535.0;
     return r>ep?pow(r,1/3.0):(ka*r+16)/116.0;
}

int rgb2Lab(int rgb[], double Lab[])
{
    int c;
    double xyz[3]={0,0,0};

	// Convert RGB to XYZ
	for (c=0; c < 3; c++)
	{
		xyz[c]+=rgb_xyz[c][0]*(double)rgb[0];
		xyz[c]+=rgb_xyz[c][1]*(double)rgb[1];
		xyz[c]+=rgb_xyz[c][2]*(double)rgb[2];
	}

	// Convert XYZ to L*a*b*
	for (c=0; c < 3; c++)
		xyz[c]=f_cbrt(xyz[c]/d50_white[c]);

	Lab[0]=116.0*xyz[1]-16.0;
	Lab[1]=500.0*(xyz[0]-xyz[1]);;
	Lab[2]=200.0*(xyz[1]-xyz[2]);

}


int Lab2rgb(int rgb[], double Lab[])
{

    int c;
    double L;
    double xyz[3],f[3], rgbt[3]={0,0,0};


	// Convert L*a*b* to XYZ
	L=(double)Lab[0];
	f[1]=(L+16.0)/116.0;    // fy
	f[0]=f[1]+(double)Lab[1]/500.0; // fx
	f[2]=f[1]-(double)Lab[2]/200.0; // fz

	xyz[0]=65535.0*d50_white[0]*(f[0]*f[0]*f[0]>ep?f[0]*f[0]*f[0]:(116.0*f[0]-16.0)/ka);
	xyz[1]=65535.0*d50_white[1]*(L>ka*ep?pow(f[1],3.0):L/ka);
	xyz[2]=65535.0*d50_white[2]*(f[2]*f[2]*f[2]>ep?f[2]*f[2]*f[2]:(116.0*f[2]-16.0)/ka);

	// Convert XYZ to RGB
	for (c=0; c < 3; c++)
	{
		rgbt[0]+=xyz[c]*xyz_rgb_2[c][0];
		rgbt[1]+=xyz[c]*xyz_rgb_2[c][1];
		rgbt[2]+=xyz[c]*xyz_rgb_2[c][2];
	}

	for (c=0; c < 3; c++)
	 rgb[c]=(unsigned short)CLIP(rgbt[c]);

}


int rgb2YCbCr(int rgb[], int YCbCr[])
{

}

int YCbCr2rgb(int rgb[], int YCbCr[])
{

}

// color conversion formula taken from http://msdn.microsoft.com/en-us/library/ms893078.aspx

int rgb2Yuv(int rgb[], int Yuv[])
{

	double b,g,r, Y, u, v;

	b = rgb[0];
	g = rgb[1];
	r = rgb[2];

	int sign2 = 1, sign1 = 1, sign0 = 1;
	int temp2 =   66 * r + 129 * g +  25 * b + 128;
	int temp1 =  -38 * r -  74 * g + 112 * b + 128;
	int temp0 =  112 * r -  94 * g -  18 * b + 128;

	if (temp2<0)
	{
		sign2 = -1;
		temp2 = -temp2;
	}

	if (temp1<0)
	{
		sign1 = -1;
		temp1 = -temp1;
	}

	if (temp0<0)
	{
		sign0 = -1;
		temp0 = -temp0;
	}

	Yuv[2] = sign2*(temp2 >> 8) +  16; //Y
	Yuv[1] = sign1*(temp1 >> 8) + 128; //u
	Yuv[0] = sign0*(temp0 >> 8) + 128; //v

}

int Yuv2rgb(int rgb[], int Yuv[])
{

	double C, D, E;

	C = Yuv[2] - 16;
	D = Yuv[1] - 128;
	E = Yuv[0] - 128;

	rgb[2] = (int( 298 * C           + 409 * E + 128) >> 8);
	rgb[1] = (int( 298 * C - 100 * D - 208 * E + 128) >> 8);
	rgb[0] = (int( 298 * C + 516 * D           + 128) >> 8);

}





// intensity features
// for grayscale
double getCostIntensityDiscBins(unsigned int pix1, int nbins, double binner, int range, int d, std::vector<float> thetaU, int &thetaid)
{
	int binid;

	if (pix1 >= range)
		  binid = nbins-1;
	else
		  binid = int(pix1/binner);

//	if (binid < 5 || binid > 110 || d>=2)
//	  std::cout << " " << (int)pix1;

	if (binid >= nbins)
		binid = nbins - 1;



	thetaid = d*nbins + binid;

	return thetaU[thetaid];
}



double getCostIntensityDiscBins(unsigned int pix1, int nbins, int d, std::vector<float> thetaU, int &thetaid)
{
	int binid;
	double binner = 256.0/nbins;

	if (pix1 >= 256.0)
		  binid = nbins-1;
	else
		  binid = int(pix1/binner);

	thetaid = d*nbins + binid;

	return thetaU[thetaid];
}

// for RGB
// this uses cubes as bins
double getCostIntensityDiscBins(unsigned char* pix1, int nbins, int d, std::vector<float> thetaU, int &thetaid)
{
	int r = (int)pix1[0];
	int g = (int)pix1[1];
	int b = (int)pix1[2];

	double nbinsCR3 = pow(abs(nbins), 1.0/3.0);

	int binidr, binidg, binidb;
	double binner = 256.0/nbinsCR3;

	binidr = int(r/binner);
	binidg = int(g/binner);
	binidb = int(b/binner);

	thetaid = d*nbins + binidr*64*16 + binidg*8*4 + binidb;

	return thetaU[thetaid];
}

// RGB - this uses linear combination of R, G and B
double getCostIntensityDiscBins(unsigned char* pix1, int nbins, int d, std::vector<float> thetaU, int thetaid[])
{
	int r = (int)pix1[0];
	int g = (int)pix1[1];
	int b = (int)pix1[2];

	double nbinsby3 = nbins/3.0;

	int binidr, binidg, binidb;
	double binner = 256.0/nbinsby3;

	binidr = int(r/binner);
	binidg = int(g/binner);
	binidb = int(b/binner);

	thetaid[0] = d*nbins + binidb;
	thetaid[1] = d*nbins + binidg + nbinsby3;
	thetaid[2] = d*nbins + binidr + nbinsby3*2;

/*	int v1=thetaid[0];
	int v2=thetaid[1];
	int v3=thetaid[2];

	int sz=thetaU.size();
*/
	return thetaU[thetaid[0]] + thetaU[thetaid[1]] + thetaU[thetaid[2]];
}

// Yuv Y has its histogram and u,v have a joint space
// range of Yuv == Y(16, 235) u(17, 240) v(17, 240)

double getCostIntensityDiscBinsYuv(unsigned char* pix1, int nbins, int d, std::vector<float> thetaU, int thetaid[])
{
	int rgb[3], Yuv[3];

	rgb[0] = (int)pix1[0];  // b
	rgb[1] = (int)pix1[1];  // g
	rgb[2] = (int)pix1[2];  // r

	rgb2Yuv(rgb, Yuv);

	// this assumes that the bins are divided in the following manner
	// (nb*nb) + (nb) ==> first term is area of (u,v) divided into bins and Y divided in nb bins

	// nb^2 + nb - nbins = 0, solve and use positive integer value, if real number it is advisable to change parameters listed in parameter file

	double nbinsY = int((-1 + sqrt(1+4*nbins))/2);
	double nbinsuv = nbinsY * nbinsY;

	double bwY = 220/nbinsY;
	double bwu = 224/nbinsY;
	double bwv = 224/nbinsY;

	double idxY = (Yuv[2] - 15)/bwY;
	if (idxY>=nbinsY)
		idxY = nbinsY - 1;

	double idxu = (Yuv[1] - 16)/bwY;
	if (idxu>=nbinsY)
		idxu = nbinsY - 1;

	double idxv = (Yuv[0] - 16)/bwY;
	if (idxv>=nbinsY)
		idxv = nbinsY - 1;

	int idxuv = idxu*nbinsY + idxv;
	if (idxuv>=nbinsuv)
		idxuv = nbinsuv - 1;

	thetaid[0] = d*nbins + idxY;
	thetaid[1] = d*nbins + nbinsY + idxuv;

	return thetaU[thetaid[0]] + thetaU[thetaid[1]];
}


double getCostIntensityDiscBinsuv(unsigned char* pix1, int nbins, int d, std::vector<float> thetaU, int thetaid[])
{
  int rgb[3], Yuv[3];

  rgb[0] = (int)pix1[0];  // b
  rgb[1] = (int)pix1[1];  // g
  rgb[2] = (int)pix1[2];  // r

  rgb2Yuv(rgb, Yuv);

  // this assumes that the bins are divided in the following manner
  // (nb*nb) + (nb) ==> first term is area of (u,v) divided into bins and Y divided in nb bins

  // nb^2 + nb - nbins = 0, solve and use positive integer value, if real number it is advisable to change parameters listed in parameter file

  double nbinsY = int((-1 + sqrt(1+4*nbins))/2);
  double nbinsuv = nbinsY * nbinsY;

  double bwY = 220/nbinsY;
  double bwu = 224/nbinsY;
  double bwv = 224/nbinsY;

  double idxY = (Yuv[2] - 15)/bwY;
  if (idxY>=nbinsY)
    idxY = nbinsY - 1;

  double idxu = (Yuv[1] - 16)/bwY;
  if (idxu>=nbinsY)
    idxu = nbinsY - 1;

  double idxv = (Yuv[0] - 16)/bwY;
  if (idxv>=nbinsY)
    idxv = nbinsY - 1;

  int idxuv = idxu*nbinsY + idxv;
  if (idxuv>=nbinsuv)
    idxuv = nbinsuv - 1;

  thetaid[0] = d*nbins + idxY;
  thetaid[1] = d*nbins + nbinsY + idxuv;

  return thetaU[thetaid[1]];
}

// nbins is the number of parameters per class per opflowlocal frame.

double getCostOpflowBinsuv(double u, double v, int nbins, int d, std::vector<float> thetaO, int &thetaid, int numframes, int frameidx)
{
  // this assumes that the bins are divided in the following manner
  // (nb*nb) ==> first term is area of (u,v) divided into bins

  // nb^2 + nb - nbins = 0, solve and use positive integer value, if real number it is advisable to change parameters listed in parameter file

  double nbinss = int(sqrt(nbins));
  double nbinsuv = nbins;

  double bwu = 400/nbinss;  // range is -200 to 200 ie jumps of 1/4th screen for 800*800 image
  double bwv = 400/nbinss;  // range is -200 to 200

  double idxu = (u + 200);
  if (idxu < 0)
    idxu = 0;
  if (idxu > 400)
    idxu = 400;

  idxu = idxu/bwu;

  if (idxu>=nbinss)
    idxu = nbinss - 1;

  double idxv = (v + 200);
  if (idxv < 0)
    idxv = 0;
  if (idxv > 400)
    idxv = 400;

  idxv = idxv/bwv;

  if (idxv>=nbinss)
    idxv = nbinss - 1;

  int idxuv = idxu*nbinss + idxv;
  if (idxuv>=nbinsuv)
    idxuv = nbinsuv - 1;

  thetaid = d*numframes*nbins + frameidx*nbins + idxuv;

  return thetaO[thetaid];
}




// Lab L has its histogram and a,b have a joint space
// range of Lab == L(0, 100) a(-80, 94) b(-113, 94)

double getCostIntensityDiscBinsLab(unsigned char* pix1, int nbins, int d, std::vector<float> thetaU, int thetaid[])
{
	int rgb[3];
	double Lab[3];

	rgb[0] = (int)pix1[0];  // b
	rgb[1] = (int)pix1[1];  // g
	rgb[2] = (int)pix1[2];  // r

	rgb2Lab(rgb, Lab);

	// this assumes that the bins are divided in the following manner
	// (nb*nb) + (nb) ==> first term is area of (u,v) divided into bins and Y divided in nb bins

	// nb^2 + nb - nbins = 0, solve and use positive integer value, if real number it is advisable to change parameters listed in parameter file

	double nbinsL = int((-1 + sqrt(1+4*nbins))/2);
	double nbinsab = nbinsL * nbinsL;

	double bwL = 101/nbinsL;
	double bwa = 175/nbinsL;
	double bwb = 208/nbinsL;

	double idxL = (Lab[2] + 1)/bwL;
	if (idxL>=nbinsL)
		idxL = nbinsL - 1;

	double idxa = (Lab[1] + 81)/bwL;
	if (idxa>=nbinsL)
		idxa = nbinsL - 1;

	double idxb = (Lab[0] + 114)/bwL;
	if (idxb>=nbinsL)
		idxb = nbinsL - 1;

	int idxab = idxa*nbinsL + idxb;
	if (idxab>=nbinsab)
		idxab = nbinsab - 1;

	thetaid[0] = d*nbins + idxL;
	thetaid[1] = d*nbins + nbinsL + idxab;

	return thetaU[thetaid[0]] + thetaU[thetaid[1]];
}




double getCostIntensityGenGaussian(intClass *intObj, unsigned short pix1)
{
	  float intcost;
	  float *pii;
	  float **mu;
	  float ***SigmaInv;
	  float *normConst;

	  pii = intObj->getPi();
	  mu = intObj->getMu();
	  SigmaInv = intObj->getSigmaInv();
	  normConst = intObj->getNormConst();

	  double intval = 0.0;

	  for (int c=0; c<intObj->getNClusters(); c++)
	  {

		  float xs = (float)pix1-mu[c][0];

		  double tempval = xs*xs*SigmaInv[c][0][0];

		  intval += exp(log(pii[c])+(-tempval/2)-log(normConst[c]));
	  }

	  if (intval==0.0)
		  intval = pow(10.0, -99.0);

	  intcost = -log(intval);

	  return intcost;

}


double getCostIntensityGenNB(unsigned int pix1, int nbinsNB, int d, intensityNB* intNBclass)
{
	int binid;
	double binner = 1380.0/nbinsNB;

	if (pix1 >= 1380)
		  binid = nbinsNB-1;
	else
		  binid = int(pix1/binner);

	double intcost= intNBclass->getPofXcYbd(binid, d);

	if (intcost==0.0)
		  intcost = pow(10.0,-99.0);

	intcost = -log(intcost);

	return intcost;
}









// appearance features

double getCostAppDiscPatch(int clustNo , int patchSize, int x, int y, int width, int height, std::vector<float> thetaA, int d, int nT, int &thetaid)
{

	double appcost = 0;
	thetaid = -1;

	if (x>(int)patchSize/2 && x<width-(int)patchSize/2 && y>(int)patchSize/2 && y<height-(int)patchSize/2)
	{
	  thetaid = clustNo + d*nT;
	  appcost =  thetaA[thetaid];
	}

	return appcost;

}


double getCostAppGenPatch(double appcostp, int dim, int x, int y, int width, int height)
{

	double appcost = 0;

	int patchSize = sqrt(dim);

	if (x>(int)patchSize/2 && x<width-patchSize/2 && y>(int)patchSize/2 && y<height-patchSize/2)
	{
		  appcost = appcostp;
	}

	return appcost;

}


// hog features
double getCostHoGDiscBinsTemp(int hogbin, int d, std::vector<float> thetaH, int &thetaid)
{
	int classv = -1;

	if (hogbin<350)
		classv = 0;
	else if (hogbin < 560)
		classv = 1;
	else if (hogbin < 610)
		classv = 2;
	else if (hogbin < 660)
		classv = 3;
	else if (hogbin < 700)
		classv = 4;
	else
		classv = 5;

	thetaid = -1;
	if (classv == d)
	{
		thetaid = hogbin;
		return thetaH[thetaid];
	}

	return 0;

}


double getCostHoGDiscBins(int hogbin, int d, int nHpC, std::vector<float> thetaH, int &thetaid)
{

	thetaid = d*nHpC + hogbin;
	// std::cout << hogbin << "  " << thetaid << "  " << thetaH.size() << std::endl;
	return thetaH[thetaid];

}


double getCostHoGGenBins(double hogcostP, int x, int y, int width, int height, int startPatch)
{

	//assuming i already have the data in log space
	// and that log(0) = max-val

	double hogcost = 0;

	if (x>(int)startPatch && x<width-startPatch && y>(int)startPatch && y<height-startPatch)
	{
		  hogcost = -hogcostP; //remember it is -log()
	}

	return hogcost;

}

// bias features
double getCostBias(int d, std::vector<float> thetaB, int &thetaid)
{
  thetaid = d;
  return thetaB[thetaid];
}

// volume features
double getCostInverseClassSize(int d, std::vector<float> inverseClassSize, std::vector<float> thetaE, int &thetaid)
{
  thetaid = d;
  return thetaE[thetaid] * inverseClassSize[d];
}



// motion features

double getCostmotDiscBins(int motbin, int d, int nMpC, std::vector<float> thetaM, int &thetaid)
{

	thetaid = d*nMpC + motbin;

	return thetaM[thetaid];

}


double getCostmotGenBins(double motcostP, int x, int y, int width, int height, int startPatch)
{

	//assuming i already have the data in log space
	// and that log(0) = max-val

	double motcost = 0;

	if (x>(int)startPatch && x<width-startPatch && y>(int)startPatch && y<height-startPatch)
	{
		  motcost = -motcostP; //remember it is -log()
	}

	return motcost;

}





// location features


// we have parameters for each cluster (and so is weaker)
double getCostLocationDiscGaussian(LocationMV** LocationMVclass, unsigned int i, int zslice, int x, int y, std::vector<float> thetaL, int d, int loc, int genparam, int &thetaid)
{

	float loccost=0;
	float *pii;
	float **mu;
	float ***SigmaInv;
	float *normConst;

	pii = LocationMVclass[d]->getPi();
	mu = LocationMVclass[d]->getMu();
	SigmaInv = LocationMVclass[d]->getSigmaInv();
	normConst = LocationMVclass[d]->getNormConst();

	//					  LocationMVclass[d]->printParameters();

	double LocationMVval = 0.0;
	int maxclass = 0;
	double maxLocationMVval = 0.0;

	for (int c=0; c<LocationMVclass[d]->getNClusters(); c++)
	{

		  float zs = i + zslice - mu[c][2];
		  float xs = y-mu[c][0];
		  float ys = x-mu[c][1];

		  double tempval = xs*xs*SigmaInv[c][0][0] + xs*ys*(SigmaInv[c][1][0] + SigmaInv[c][0][1]) + xs*zs*(SigmaInv[c][2][0] + SigmaInv[c][0][2]) +
						  ys*ys*SigmaInv[c][1][1] + ys*zs*(SigmaInv[c][2][1] + SigmaInv[c][1][2]) + zs*zs*SigmaInv[c][2][2];

		  LocationMVval = ((-tempval/2)-log(normConst[c])); // no need to take exp
		  if (c==0) {
			  maxclass = 0;
			  maxLocationMVval = LocationMVval;
		  } else {
			  if (maxLocationMVval<LocationMVval) {
				  maxclass = c;
				  maxLocationMVval = LocationMVval;
			  }
		  }

	}
	// find the theta corresponding to c
	thetaid = 0;
	for (int km=0; km<d; km++) {
		  thetaid += LocationMVclass[km]->getNClusters();
	}
	thetaid += maxclass;
	loccost = thetaL[thetaid];

	if (loc==8 && genparam==1) {
		  loccost -= maxLocationMVval; //remember it is -Ed , hence the -ve sign for normal distribution
	}

	return loccost;

}


// we have parameters for each cluster,class pair.
// ie each class will have the parameters for each cluster (even if some clusters don't belong to the class)
// this will try to make the parameters of the clusters belonging to the class (more so the one cluster closest) with a higher impact value
// than the rest of the cluster.
double getCostLocationDiscGaussian2(LocationMV** LocationMVclass, unsigned int z, int zslice, int x, int y, std::vector<float> thetaL, int d, int loc, int genparam, int nD, int &thetaid)
{

  int globalNdisps = nD;

	int nL = thetaL.size();
	int nLpC = nL/globalNdisps;

	float loccost=0;
	float *pii;
	float **mu;
	float ***SigmaInv;
	float *normConst;

	double LocationMVval = 0.0;
	int maxclass = 0;
	double maxLocationMVval = 0.0;

	int totClusters = 0;

	for (int ff=0; ff<nD; ff++) {

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

	if (totClusters != nLpC) {
		  std::cout<<"Error in number of clusters for location"<<endl;
		  exit(1);
	}

	// find the theta corresponding to c
	thetaid = nLpC*d + maxclass;
	loccost = thetaL[thetaid];

	if (loc==9 && genparam==1) {
		  loccost -= maxLocationMVval; //remember it is -Ed , hence the -ve sign for normal distribution
	}

	return loccost;

}


double getCostLocationGenGaussian(LocationMV* LocationMVclass, unsigned int z, int zslice, int x, int y)
{

	float loccost;
	float *pii;
	float **mu;
	float ***SigmaInv;
	float *normConst;

	pii = LocationMVclass->getPi();
	mu = LocationMVclass->getMu();
	SigmaInv = LocationMVclass->getSigmaInv();
	normConst = LocationMVclass->getNormConst();

	float LocationMVval = 0.0;

	for (int c=0; c<LocationMVclass->getNClusters(); c++)
	{

		  float zs = z + zslice - mu[c][2];
		  float xs = y-mu[c][0];
		  float ys = x-mu[c][1];

		  float tempval = xs*xs*SigmaInv[c][0][0] + xs*ys*(SigmaInv[c][1][0] + SigmaInv[c][0][1]) + xs*zs*(SigmaInv[c][2][0] + SigmaInv[c][0][2]) +
						  ys*ys*SigmaInv[c][1][1] + ys*zs*(SigmaInv[c][2][1] + SigmaInv[c][1][2]) + zs*zs*SigmaInv[c][2][2];

		  LocationMVval += exp(log(pii[c])+(-tempval/2)-log(normConst[c]));
	}

	if (LocationMVval==0.0)
		  LocationMVval = pow(10.0, -99.0);  //WORK-THIS

	loccost = -log(LocationMVval);

	return loccost;

}





double getCostLocationDiscCube(unsigned short locpix, int x, int y, int i, int loccubex, int loccubey, int loccubez, int numTrainingPats, int nLpC, std::vector<float> thetaL, int width, int height, int depth, int d, int &thetaid)
{

	double locpixd = (double)locpix/ (loccubex * loccubey * loccubez * numTrainingPats);

	int colno = x/loccubex;
	int rowno = y/loccubey;
	int hno = i/loccubez;

	int maxwidth = width/loccubex;
	if (width%loccubex!=0)
		  maxwidth++;
	int maxheight = height/loccubey;
	if (height%loccubey!=0)
		  maxheight++;
	int maxdepth = depth/loccubez;
	if (depth%loccubez!=0)
		  maxdepth++;

	if (maxwidth*maxheight*maxdepth != nLpC) {
		  std::cout << " Location parameters don't match up!! "<<endl;
		  exit(1);
	}

	// need to check this for boundary conditions
	// especially when loccube? not multiple of size
	// or when it is actually multiple of size of image volume
	thetaid = nLpC*d + hno*maxwidth*maxheight + rowno*maxwidth + colno;

	return (thetaL[thetaid] * (1-locpixd));

}

// distance transform location
double getCostLocationDT(unsigned short locpixval, std::vector<float> thetaL, int d, int &thetaid)
{

  thetaid = 0;

  if (locpixval == 65535)
    return 0; // 65535 is the ignore pixel

  double val = double(locpixval)/65536.0; // this is the foreground map value
  // and so is true if d == 1

  /*
  if (val != 0)
    val = log(val);
  else
    val = -1e10;

  if (d == 0)
    val = log(1) - val;
   */

  val = exp(val)-1;
  if (d == 0)
    val = (exp(1)-1) - val;

  return (thetaL[thetaid] * val);
}


// distance transform location
double getCostLocationDTgen(unsigned short locpixval, std::vector<float> thetaL, int d, int &thetaid)
{

  thetaid = d;

  if (locpixval == 65535)
    return 0; // 65535 is the ignore pixel

  double val = double(locpixval)/65536.0; // this is the foreground map value
  // and so is true if d == 1

  //val = exp(val)-1;
  //if (d == 0)
  //  val = (exp(1)-1) - val;

  if (d == 0)
    val = 1.0 - val;

  // 0.399 => 0.5log(2*sigma^2), note that sigma = 1/sqrt(2*thetaL[0])
  val = val*val*thetaL[d] + 0.399 - 0.5*log(2*thetaL[d]);
  // should above line be negated or not
  // and should val be exponentiated further as done in discriminative version with a bias of -1?

  return (val);
}




// rbm features


// daisy features



// klr features

void readklrParameters(std::string &klrfName, int &hogdimklr, int &appdimklr, int &locdimklr, int &biasklr, double &lambda, std::string &xfile, std::string &wfile, std::string &mmfile, std::string &kernel, double &p1, kernelType &ker)
{

  string line;
  string temps;

  ifstream myfile (klrfName.c_str());
  if (myfile.is_open()) {

    while (! myfile.eof() ) {
      myfile >> line;
      if (line.compare("hogdimklr=")==0)
      {
    	  myfile >> hogdimklr;

      } else if (line.compare("appdimklr=")==0)
      {
    	  myfile >> appdimklr;
      } else if (line.compare("locdimklr=")==0)
      {
    	  myfile >> locdimklr;

      } else if (line.compare("biasklr=")==0)
      {
    	  myfile >> biasklr;

      } else if (line.compare("lambda=")==0)
      {
    	  myfile >> lambda;
      } else if (line.compare("xfile=")==0)
      {
    	  myfile >> xfile;

      }else if (line.compare("wfile=")==0)
      {
    	  myfile >> wfile;
      }else if (line.compare("mmfile=")==0)
      {
    	  myfile >> mmfile;
      } else if (line.compare("kernel=")==0)
      {
    	  myfile >> kernel;
      } else if (line.compare("p1=")==0)
      {
    	  myfile >> p1;
      }else if (line.compare("ker=")==0)
      {
    	  myfile >> temps;
    	  //chet work: need to fix this to read enum
      }
      else
      {
        std::cout<<" klr Parameter file corrupt \n";
        exit(1);
      }
    }

    myfile.close();

  } else {
    std::cout << "Unable to open klr parameter file";
    exit(1);
  }

}


void readklrParameters(klrparams &klrpv)
{

  string line;
  string temps;

  ifstream myfile (klrpv.klrfName.c_str());
  if (myfile.is_open()) {

    while (! myfile.eof() ) {
      myfile >> line;
      if (line.compare("hogdimklr=")==0)
      {
    	  myfile >> klrpv.hogdimklr;

      } else if (line.compare("appdimklr=")==0)
      {
    	  myfile >> klrpv.appdimklr;
      } else if (line.compare("locdimklr=")==0)
      {
    	  myfile >> klrpv.locdimklr;

      } else if (line.compare("biasklr=")==0)
      {
    	  myfile >> klrpv.biasklr;

      } else if (line.compare("lambda=")==0)
      {
    	  myfile >> klrpv.lambda;
      } else if (line.compare("xfile=")==0)
      {
    	  myfile >> klrpv.xfile;

      }else if (line.compare("wfile=")==0)
      {
    	  myfile >> klrpv.wfile;
      }else if (line.compare("mmfile=")==0)
      {
    	  myfile >> klrpv.mmfile;
      } else if (line.compare("kernel=")==0)
      {
    	  myfile >> klrpv.kernel;
      } else if (line.compare("p1=")==0)
      {
    	  myfile >> klrpv.p1;
      }else if (line.compare("ker=")==0)
      {
    	  myfile >> temps;
    	  //chet work: need to fix this to read enum
      }
      else
      {
        std::cout<<" klr Parameter file corrupt \n";
        exit(1);
      }
    }

    myfile.close();

  } else {
    std::cout << "Unable to open klr parameter file";
    exit(1);
  }

}





void printklrParameters(std::string &klrfName, int hogdimklr, int appdimklr, int locdimklr, int biasklr, double lambda, std::string xfile, std::string wfile, std::string mmfile, std::string kernel, double p1, kernelType ker)
{

	std::cout << " filename : " << klrfName << std::endl;
	std::cout << " hog dims : " << hogdimklr << std::endl;
	std::cout << " app dims : " << appdimklr << std::endl;
	std::cout << " loc dims : " << locdimklr << std::endl;
	std::cout << " bias  : " << biasklr << std::endl;
	std::cout << " lambda : " << lambda << std::endl;
	std::cout << " xfile : " << xfile << std::endl;
	std::cout << " wfile : " << wfile << std::endl;
	std::cout << " mmfile : " << mmfile << std::endl;
	std::cout << " kernel : " << kernel << std::endl;
	std::cout << " p1 : " << p1 << std::endl;
	std::cout << " ker : " << ker << std::endl;

}


void printklrParameters(klrparams klrpv)
{

	std::cout << " filename : " << klrpv.klrfName << std::endl;
	std::cout << " hog dims : " << klrpv.hogdimklr << std::endl;
	std::cout << " app dims : " << klrpv.appdimklr << std::endl;
	std::cout << " loc dims : " << klrpv.locdimklr << std::endl;
	std::cout << " bias  : " << klrpv.biasklr << std::endl;
	std::cout << " lambda : " << klrpv.lambda << std::endl;
	std::cout << " xfile : " << klrpv.xfile << std::endl;
	std::cout << " wfile : " << klrpv.wfile << std::endl;
	std::cout << " mmfile : " << klrpv.mmfile << std::endl;
	std::cout << " kernel : " << klrpv.kernel << std::endl;
	std::cout << " p1 : " << klrpv.p1 << std::endl;
	std::cout << " ker : " << klrpv.ker << std::endl;

}






/*
int skipXYZlocklr;

int setskipXYZlocklr(int val)
{
	skipXYZlocklr = val;
}
*/


void getCostLocationKlr(ublas::vector<DBL_TYPE> &Xloctest, unsigned int i, int zslice, int x, int y, int nD, matrix<DBL_TYPE> minmaxM, LocationMV** LocationMVclass, int &idxx, int skipXYZlocklr)
{

	float zs1 = i + zslice;
	float xs1 = y + 1;
	float ys1 = x + 1;

	if (skipXYZlocklr==0)
	{
		Xloctest(idxx) = (xs1-minmaxM(0,idxx))/(minmaxM(1,idxx));
		idxx++;
		Xloctest(idxx) = (ys1-minmaxM(0,idxx))/(minmaxM(1,idxx));;
		idxx++;
		Xloctest(idxx) = (zs1-minmaxM(0,idxx))/(minmaxM(1,idxx));;
		idxx++;
	}

	float loccost;
	float *pii;
	float **mu;
	float ***SigmaInv;
	float *normConst;


	for (int d = 0; d < nD; d++) {

		  pii = LocationMVclass[d]->getPi();
		  mu = LocationMVclass[d]->getMu();
		  SigmaInv = LocationMVclass[d]->getSigmaInv();
		  normConst = LocationMVclass[d]->getNormConst();

		  for (int c=0; c<LocationMVclass[d]->getNClusters(); c++)
		  {
			  float zs = i + zslice - mu[c][2] - 1;
			  float xs = y + 1 - mu[c][0];
			  float ys = x + 1 - mu[c][1];

			  double tempval = xs*xs*SigmaInv[c][0][0] + xs*ys*(SigmaInv[c][1][0] + SigmaInv[c][0][1]) + xs*zs*(SigmaInv[c][2][0] + SigmaInv[c][0][2]) +
							  ys*ys*SigmaInv[c][1][1] + ys*zs*(SigmaInv[c][2][1] + SigmaInv[c][1][2]) + zs*zs*SigmaInv[c][2][2];

			  double lv = ((-tempval/2)-log(normConst[c]));

			  Xloctest(idxx) = (lv-minmaxM(0,idxx))/(minmaxM(1,idxx));

			  idxx++;
		  }

	}


	for (int d = 0; d < nD; d++) {

		  pii = LocationMVclass[d]->getPi();
		  mu = LocationMVclass[d]->getMu();
		  SigmaInv = LocationMVclass[d]->getSigmaInv();
		  normConst = LocationMVclass[d]->getNormConst();

		  double sumprob = 0.0;

		  for (int c=0; c<LocationMVclass[d]->getNClusters(); c++)
		  {
			  float zs = i + zslice - mu[c][2] - 1;
			  float xs = y + 1 - mu[c][0];
			  float ys = x + 1 - mu[c][1];

			  double tempval = xs*xs*SigmaInv[c][0][0] + xs*ys*(SigmaInv[c][1][0] + SigmaInv[c][0][1]) + xs*zs*(SigmaInv[c][2][0] + SigmaInv[c][0][2]) +
							  ys*ys*SigmaInv[c][1][1] + ys*zs*(SigmaInv[c][2][1] + SigmaInv[c][1][2]) + zs*zs*SigmaInv[c][2][2];

			  sumprob += (log(pii[c])+(-tempval/2)-log(normConst[c]));

		  }

		  Xloctest(idxx) = (sumprob-minmaxM(0,idxx))/(minmaxM(1,idxx));

		  idxx++;

	}

}



// grayscale
void getCostIntensityKlr(ublas::vector<DBL_TYPE> &Xtest, double pix1, matrix<DBL_TYPE> minmaxM, int &idxx)
{
	DBL_TYPE pixv = pix1;

	if (pix1>400)
		  pixv = 0;

	DBL_TYPE normpix = (pixv-minmaxM(0,idxx))/(minmaxM(1,idxx));
	Xtest(idxx) = normpix;
	idxx++;
}

// Yuv
void getCostIntensityKlrYuv(ublas::vector<DBL_TYPE> &Xtest, unsigned char* pix1, matrix<DBL_TYPE> minmaxM, int &idxx)
{

	int rgb[3], Yuv[3];

	rgb[0] = (int)pix1[0];  // b
	rgb[1] = (int)pix1[1];  // g
	rgb[2] = (int)pix1[2];  // r

	rgb2Yuv(rgb, Yuv);

	// Y is Yuv[2] etc
	DBL_TYPE normpix = (Yuv[2]-minmaxM(0,idxx))/(minmaxM(1,idxx));
	Xtest(idxx) = normpix;
	idxx++;
	normpix = (Yuv[1]-minmaxM(0,idxx))/(minmaxM(1,idxx));
	Xtest(idxx) = normpix;
	idxx++;
	normpix = (Yuv[0]-minmaxM(0,idxx))/(minmaxM(1,idxx));
	Xtest(idxx) = normpix;
	idxx++;
}


void getCostAppKlr(ublas::vector<DBL_TYPE> &Xtest, int x, int y, int width, matrix<DBL_TYPE> minmaxM, int &idxx, int appdim, matrix<DBL_TYPE> *appdescs)
{

	for (int jk=0; jk<appdim; jk++)
	{
		  Xtest(idxx) = (double((*appdescs)(y*width+x,jk))-minmaxM(0,idxx))/(minmaxM(1,idxx));
		  idxx++;
	}

}

void getCostHoGKlr(ublas::vector<DBL_TYPE> &Xtest, int x, int y, int width, matrix<DBL_TYPE> minmaxM, int &idxx, int hogdim, matrix<DBL_TYPE> *hogdescs)
{

	for (int jk=0; jk<hogdim; jk++)
	{
		  Xtest(idxx) = (double((*hogdescs)(y*width+x,jk))-minmaxM(0,idxx))/(minmaxM(1,idxx));
		  idxx++;
	}

}





// optical flow features
// this function returns the correct theta value and not theta*opflow_u or theta*opticalflow_v
void getIndicatorCostOpticalFlow(int depth, int z, int d, std::vector<float> thetaO, int thetaid[], float indcost[], int numframes, int frameidx)
{
  int pC = numframes * 2; // 2 because of u and v, pC is more like per class

  thetaid[0] = d*pC + frameidx*2;
  thetaid[1] = d*pC + frameidx*2 + 1;

  // indicator cost ie theta values
  indcost[0] = thetaO[thetaid[0]];
  indcost[1] = thetaO[thetaid[1]];
}









// reading parameter files

void readLocationMVParameters(LocationMV** &LocationMVclass, int nD)
{

  LocationMVclass = new LocationMV*[nD];
  //need to write the destructor for this
  // and anything related to this

  // read parameter files
  //read bgnd, liver, rk, lk, gb, spleen
  for (int kk=0; kk<nD; kk++) {
    string line;
    int check = 0;
    char cname[50], fname[100];
    int numC, dim;
/*
    if (kk==0)
      strcpy(fname,"/p/imageseg/medicalSeg/parameters/location/mcovbgnd_7_14n.txt");
    else if (kk==1)
      strcpy(fname,"/p/imageseg/medicalSeg/parameters/location/mcovliv_7_14n.txt");
    else if (kk==2)
      strcpy(fname,"/p/imageseg/medicalSeg/parameters/location/mcovrk_7_14n.txt");
    else if (kk==3)
      strcpy(fname,"/p/imageseg/medicalSeg/parameters/location/mcovlk_7_14n.txt");
    else if (kk==4)
      strcpy(fname,"/p/imageseg/medicalSeg/parameters/location/mcovgb_7_14n.txt");
    else if (kk==5)
      strcpy(fname,"/p/imageseg/medicalSeg/parameters/location/mcovspl_7_14n.txt");
*/

    if (kk==0)
      strcpy(fname,"/p/imageseg/medicalSeg/parameters/location/mcovbgndTrial001.txt");
    else if (kk==1)
      strcpy(fname,"/p/imageseg/medicalSeg/parameters/location/mcovlivTrial001.txt");
    else if (kk==2)
      strcpy(fname,"/p/imageseg/medicalSeg/parameters/location/mcovrkTrial001.txt");
    else if (kk==3)
      strcpy(fname,"/p/imageseg/medicalSeg/parameters/location/mcovlkTrial001.txt");
    else if (kk==4)
      strcpy(fname,"/p/imageseg/medicalSeg/parameters/location/mcovgbTrial001.txt");
    else if (kk==5)
      strcpy(fname,"/p/imageseg/medicalSeg/parameters/location/mcovspTrial001.txt");

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

    DEBUG_OUT3(verbose, debugfile, "cname = %s numC = %d dim = %d\n", cname, numC, dim);

    LocationMVclass[kk] = new LocationMV(numC, dim);
    LocationMVclass[kk]->setClassname(cname);
    LocationMVclass[kk]->readParameters(fname);
    //			LocationMVclass[kk]->printParameters();
  }

}



// clusters per class
void readApprearanceParametersCPC(Appearance** &appclass, int nD, std::string &appfName, int colordims)
{

  appclass = new Appearance*[nD];
  // need to write the destructor of this and what follows

  // read parameter file: order of texton centers
  //read bgnd, liver, rk, lk, gb, spleen

  string line;
  int numT, patchsize;

  ifstream myfile (appfName.c_str());
  if (myfile.is_open()) {
 
    while (! myfile.eof() ) {
      myfile >> line;
      if (line.compare("numClusters")==0)
        {
          // good to go
          myfile >> line; // take in "="
          myfile >> numT;
          myfile >> line; // take in appBoxSize
          myfile >> line; // take in "="
          myfile >> patchsize;
          break;

        } else
        {
          std::cout<<" appearance Parameter file corrupt \n";
          exit(1);
        }
    }

    myfile.close();

  } else {
    std::cout << "Unable to open appearance parameter file";
    exit(1);
  }

  for (int kk=0; kk<nD; kk++)
    {
      appclass[kk] = new Appearance(numT/nD, patchsize, colordims);
      appclass[kk]->setClassname("general");

      appclass[kk]->readParameters(appfName, nD, kk, colordims);
      //			appclass[kk]->printParameters(colordims);
    }

}




// common clusters for all classes
void readApprearanceParameters(Appearance** &appclass, int nD, std::string &appfName)
{

  appclass = new Appearance*[1];
  // need to write the destructor of this and what follows

  // read parameter file: order of texton centers
  //read bgnd, liver, rk, lk, gb, spleen

  string line;
  int numT, patchsize;

  ifstream myfile (appfName.c_str());
  if (myfile.is_open())
    {
      while (! myfile.eof() )
        {
          myfile >> line;
          if (line.compare("numClusters")==0)
            {
              // good to go
              myfile >> line; // take in "="
              myfile >> numT;
              myfile >> line; // take in appBoxSize
              myfile >> line; // take in "="
              myfile >> patchsize;
              break;

            } else
            {
              std::cout<<" appearance Parameter file corrupt \n";
              exit(1);
            }
        }

      myfile.close();

    } else {
    std::cout << "Unable to open appearance parameter file";
    exit(1);
  }

  for (int kk=0; kk<1; kk++)
    {
      appclass[kk] = new Appearance(numT, patchsize, 1);
      appclass[kk]->setClassname("general");

      appclass[kk]->readParameters(appfName, 1, kk, 1);
      //			appclass[kk]->printParameters();
    }
}


void readAppearanceMVParameter(AppearanceMV** &AppMVclass, int nD, std::string &appfName)
{

  AppMVclass = new AppearanceMV*[nD];
  //need to write the destructor for this
  // and anything related to this

  // read parameter files
  //read bgnd, liver, rk, lk, gb, spleen
  for (int kk=0; kk<nD; kk++)
    {
      string line;
      int check = 0;
      char cname[50], fname[200];
      int numC, dim;

      strcpy(fname,appfName.c_str());

      if (kk==0)
        strcat(fname,"bgnd.txt");
      else if (kk==1)
        strcat(fname,"liv.txt");
      else if (kk==2)
        strcat(fname,"rk.txt");
      else if (kk==3)
        strcat(fname,"lk.txt");
      else if (kk==4)
        strcat(fname,"gb.txt");
      else if (kk==5)
        strcat(fname,"spl.txt");


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
                      std::cout<<" AppearanceMV Parameter file corrupt \n";
                      exit(1);
                    }
                }
            }
          myfile.close();
        }
      else {
        std::cout << "Unable to open AppearanceMV parameter file";
        exit(1);
      }

      DEBUG_OUT3(verbose, debugfile, "cname = %s numC = %d dim = %d\n", cname, numC, dim);

      AppMVclass[kk] = new AppearanceMV(numC, dim);
      AppMVclass[kk]->setClassname(cname);
      AppMVclass[kk]->readParameters(fname);
      //			AppMVclass[kk]->printParameters();
    }

}




void readintClassParameter(intClass** &intObj, int nD) 
{
  intObj = new intClass*[nD];
  //need to write the destructor for this
  // and anything related to this

  // read parameter files
  //read bgnd, liver, rk, lk, gb, spleen
  for (int kk=0; kk<nD; kk++)
    {
      string line;
      int check = 0;
      char cname[50], fname[50];
      int numC, dim;

      if (kk==0)
        strcpy(fname,"/p/imageseg/medicalSeg/parameters/intensity/intmcovbgnd.txt");
      else if (kk==1)
        strcpy(fname,"/p/imageseg/medicalSeg/parameters/intensity/intmcovliv.txt");
      else if (kk==2)
        strcpy(fname,"/p/imageseg/medicalSeg/parameters/intensity/intmcovrk.txt");
      else if (kk==3)
        strcpy(fname,"/p/imageseg/medicalSeg/parameters/intensity/intmcovlk.txt");
      else if (kk==4)
        strcpy(fname,"/p/imageseg/medicalSeg/parameters/intensity/intmcovgb.txt");
      else if (kk==5)
        strcpy(fname,"/p/imageseg/medicalSeg/parameters/intensity/intmcovspl.txt");


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
                      std::cout<<" intClass Parameter file corrupt \n";
                      exit(1);
                    }
                }
            }
          myfile.close();
        }
      else {
        std::cout << "Unable to open intClass parameter file";
        exit(1);
      }

      DEBUG_OUT3(verbose, debugfile, "cname = %s numC = %d dim = %d\n", cname, numC, dim);

      intObj[kk] = new intClass(numC, dim);
      intObj[kk]->setClassname(cname);
      intObj[kk]->readParameters(fname);
      //			intObj[kk]->printParameters();
    }

}








void deleteLocationMVclass(int loc, int nD, LocationMV** &LocationMVclass)
{

	if (loc == 3 || loc==6 || loc==8 || loc==7 || loc==9)
	{
		for (int i=0; i<nD; i++)
		{
			delete LocationMVclass[i];
		}
		delete [] LocationMVclass;
	}

}

void deleteAppMVclass(int app, int nD, AppearanceMV** &AppMVclass)
{

	if (app==3)
	{
		for (int i=0; i<nD; i++)
		{
			delete AppMVclass[i];
		}
		delete AppMVclass;
	}

}


void deleteappclass(int app, int nD, Appearance** &appclass)
{

	if (app==1)
	{
		for (int i=0; i<nD; i++)
		{
			delete appclass[i];
		}
		delete appclass;
	}

	if (app==2)
	{
		for (int i=0; i<1; i++)
		{
			delete appclass[i];
		}
		delete appclass;

	}

}


void deleteintObj(int intensity, int nD, intClass** &intObj)
{

	if (intensity==3)
	{
		for (int i=0; i<nD; i++)
		{
			delete intObj[i];
		}
		delete intObj;
	}

}
