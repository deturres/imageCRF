/*
 * bfgs.h
 *
 *  Created on: Mar 26, 2011
 *      Author: bhole
 */

#ifndef BFGS_H_
#define BFGS_H_

#include <vector>
#include <deque>
#include "eigen/Eigen/Dense"

class bfgs
{
  int bfgs_type; //1 bfgs, 2 lbfgs

  int vec_size;

	double beta_dash;
	double beta;
	double alpha;
	double alphastart;
	double stepmaxconst;

	// needed for the numerical recipes linesearch
  double alamin; // used in numerical recipe line search
  double slope;
  double f2; // stores previous failed fk value
  double stpmax;

  Eigen::VectorXd xk_bfgs;
  Eigen::VectorXd gk_bfgs;
  Eigen::VectorXd dk_bfgs;

  double fk_bfgs;

  Eigen::MatrixXd eye;
  Eigen::MatrixXd Hk_bfgs;

  Eigen::VectorXd xknew_bfgs;
  Eigen::VectorXd gknew_bfgs;
  double fknew_bfgs;
  Eigen::MatrixXd Hkp1_bfgs;

  int flag_a;


  // lets have lbfgs variables
  int lmemory; // M in Nocedals paper (1980, Updating Quasi Newton ...)

  Eigen::VectorXd lalpham; // temp variables to compute optimized H*d
  // this is the alpha vector in Nocedals paper
  Eigen::VectorXd lbetam;  // temp variables to compute optimized H*d
  // this is the beta vector in Nocedals paper

  // these variables have the same name as in Nocedal's paper
  int incr;
  int bound;
  int iterout;
  Eigen::MatrixXd q;
  Eigen::MatrixXd r; // r_bound is your new d_k
  std::deque<Eigen::VectorXd> s;
  std::deque<Eigen::VectorXd> y;
  std::deque<double> rho;

public:

    bfgs(int bfgs_flag, std::vector<float> theta, double betad, double beta1, double alp, int m, double stpmaxconst);
    virtual ~bfgs();

    void firstUpdate(double loglikelihood, std::vector<float> diffDist, std::vector<float> theta);
    void linesearchUpdate(double loglikelihood, std::vector<float> diffDist, std::vector<float> &theta);
    void firstOuteriterUpdate(std::vector<float> &theta);
    int outeriterUpdate(std::vector<float> &theta);

    // the nr prefix stands for the code when the linesearch from numerical recipes in C++ was used.
    void nrfirstUpdate(double loglikelihood, std::vector<float> diffDist, std::vector<float> theta);
    int nrlinesearchUpdate(double loglikelihood, std::vector<float> diffDist, std::vector<float> &theta);
    void nrStartNewDescent(std::vector<float> &theta);
    int nrouteriterUpdate(std::vector<float> &theta);

    void resetalphaTo1();
    void resetflagaTo1();
    int getflaga();
    void setbetadash(double);

    // additional functions needed by lbfgs
    void nrupdatedklimited(int iterout);
    void nrprodHkgklimited(int iterout);
    void nrStartNewDescentlimited(std::vector<float> &theta);
    int nrouteriterUpdatelimited(std::vector<float> &theta, int iterout);

};

// Note:
// if refering to Numerical Recipes in C++ (2nd edition), the following variables are changed.
// dk_bfgs is p in Numerical Recipes
// alpha is alam
// beta_dash is ALF



#endif /* BFGS_H_ */
