/*
 * bfgs.cpp
 *
 *  Created on: Mar 26, 2011
 *      Author: bhole
 */

#include "bfgs.h"
#include <iostream>

double dmaxfn(double a, double b)
{
  if (a > b)
    return a;
  else
    return b;
}


bfgs::bfgs(int bfgs_flag, std::vector<float> theta, double betad, double beta1, double alp, int m, double stpmaxconst)
{

  if (bfgs_flag == 1)
  {
    bfgs_type = 1;

    // try to set beta_dash to 1e-4 when using the numerical recipes  linesearch
    beta_dash = betad; // Ensures sufficient decrease in functionvalue.
    beta = beta1;
    alpha = alp;
    alphastart = alp;
    lmemory = 0;

    f2 = 0.0;
    slope = 0.0;
    alamin = 0.0;
    stpmax = 0.0;

    fk_bfgs = 0.0;
    fknew_bfgs = 0.0;

    vec_size = theta.size();

    xk_bfgs.resize(vec_size);
    gk_bfgs.resize(vec_size);
    xknew_bfgs.resize(vec_size);
    gknew_bfgs.resize(vec_size);

    dk_bfgs.resize(vec_size);

    for (int jin=0; jin<vec_size; jin++)
      xk_bfgs(jin) = theta[jin];

    eye.resize(vec_size, vec_size);
    eye.setIdentity();

    Hk_bfgs.resize(vec_size, vec_size);
    Hk_bfgs = eye;

    Hkp1_bfgs.resize(vec_size, vec_size);
  } else if (bfgs_flag == 2)
  {
    bfgs_type = 2;

    beta_dash = betad; // Ensures sufficient decrease in functionvalue.
    beta = beta1;
    alpha = alp;
    alphastart = alp;

    lmemory = m;
    stepmaxconst = stpmaxconst;

    f2 = 0.0;
    slope = 0.0;
    alamin = 0.0;
    stpmax = 0.0;

    fk_bfgs = 0.0;
    fknew_bfgs = 0.0;

    vec_size = theta.size();

    xk_bfgs.resize(vec_size);
    gk_bfgs.resize(vec_size);
    xknew_bfgs.resize(vec_size);
    gknew_bfgs.resize(vec_size);

    dk_bfgs.resize(vec_size);

    for (int jin=0; jin<vec_size; jin++)
      xk_bfgs(jin) = theta[jin];

    eye.resize(vec_size, vec_size);
    eye.setIdentity();

    lalpham.resize(lmemory);
    lbetam.resize(lmemory);

    incr = 0;
    bound = 0;
    iterout = 0;

    q.resize(lmemory+1, vec_size);
    // size is lmemory+1, vec_size
    r.resize(lmemory+1, vec_size);
    // size is lmemory+1, vec_size
    s.clear();
    // size is lmemory+1, vec_size
    y.clear();
    // size is lmemory
    rho.clear();

  }

}


bfgs::~bfgs()
{
}


void bfgs::firstUpdate(double loglikelihood, std::vector<float> diffDist, std::vector<float> theta)
{

    fk_bfgs = - loglikelihood;

    for (int jin=0; jin<theta.size(); jin++)
      gk_bfgs(jin) = - diffDist[jin];
    // the reason the above quantity is -ve is because
    // 1. Log likelihood is concave, but these equations are for convex function
    // log p(y|x,theta) = log(exp(-theta*func(x))) - log (Z(x, theta))
    // dL/dtheta = -empir_dist + model_dist
    // -dL/dtheta = empir_dist - model_dist
    // note that diffDist = model_dist - empir_dist, hence the -ve sign;



    for (int jin=0; jin<theta.size(); jin++)
      xk_bfgs(jin) = theta[jin];

    // break;  // break out of inner loop so i can compute dk
    flag_a = 0;
}


void bfgs::linesearchUpdate(double loglikelihood, std::vector<float> diffDist, std::vector<float> &theta)
{
    flag_a = 0;

    fknew_bfgs = - loglikelihood;

    for (int jin=0; jin<theta.size(); jin++)
      gknew_bfgs(jin) = - diffDist[jin];


    //  % if first condition not satisfied

    double tempv = gk_bfgs.transpose() * dk_bfgs;
    double tempvnew = gknew_bfgs.transpose() * dk_bfgs;

    // std::cout << " gd = " << tempv << "  gdnew = " << tempvnew << std::endl;
    // std::cout << " fnew = " << fknew_bfgs << "  f = " << fk_bfgs << std::endl;

    if (fknew_bfgs > fk_bfgs + beta_dash*alpha*tempv) {
      flag_a = 1;
      alpha = alpha * 0.8;
    } else if ((tempvnew < beta*tempv)  && flag_a == 0) {
      //  % if second condition not satisfied
      flag_a = 1;
      alpha = alpha * 0.8;
    }


    if (alpha<=1e-20) {
      flag_a = 100;  //% no step size
      // break; // break out of inner loop, no further progress possible

    } else if (flag_a == 1) {

      xknew_bfgs = xk_bfgs + alpha*dk_bfgs;

      for (int jin=0; jin<theta.size(); jin++)
        theta[jin] = xknew_bfgs(jin);

    } else if (flag_a == 0) {
      // break; // we can do an outer loop to update H now.
      // flag_a is in the inner loop and will cause it to break.
    }

}


void bfgs::firstOuteriterUpdate(std::vector<float> &theta)
{

  dk_bfgs = -1 * (Hk_bfgs*gk_bfgs);

  xknew_bfgs = xk_bfgs + alpha* dk_bfgs;

  for (int jin=0; jin<theta.size(); jin++)
    theta[jin] = xknew_bfgs(jin);

}


int bfgs::outeriterUpdate(std::vector<float> &theta)
{

	int another_flag = 1;

	vec_size = theta.size();

	Eigen::VectorXd sk_bfgs(vec_size), yk_bfgs(vec_size);
  double rhok_bfgs;
  Eigen::MatrixXd Vk_bfgs(vec_size, vec_size), Hkp1_bfgs(vec_size, vec_size);
  Eigen::MatrixXd Hkp1_bfgs_temp(vec_size, vec_size);

  // % find new Hkp1
  sk_bfgs = xknew_bfgs - xk_bfgs;
  yk_bfgs = gknew_bfgs - gk_bfgs;

  rhok_bfgs = 1/(yk_bfgs.transpose() * sk_bfgs);
  Vk_bfgs = eye - rhok_bfgs * (yk_bfgs * sk_bfgs.transpose());
  Hkp1_bfgs_temp = Vk_bfgs * Hk_bfgs;
  Hkp1_bfgs = Hkp1_bfgs_temp * Vk_bfgs + rhok_bfgs*(sk_bfgs * sk_bfgs.transpose());

  xk_bfgs = xknew_bfgs;
  gk_bfgs = gknew_bfgs;
  Hk_bfgs = Hkp1_bfgs;
  fk_bfgs = fknew_bfgs;

  // starting new search
  resetalphaTo1();
  resetflagaTo1();

  dk_bfgs = - Hk_bfgs * gk_bfgs;

  xknew_bfgs = xk_bfgs + alpha* dk_bfgs;

  for (int jin=0; jin<theta.size(); jin++)
    theta[jin] = xknew_bfgs(jin);

  //      for (int jin=0; jin<theta.size(); jin++)
  //  theta[jin] = xknew_bfgs(jin);


  double normxk = xk_bfgs.squaredNorm();
  double normgk = gk_bfgs.squaredNorm();

  double max1 = 1;

  if (normxk > max1)
    max1 = normxk;
  another_flag = 0;
  if (normgk >= 1e-5 * max1)
    another_flag = 1;

  return another_flag;
}

void bfgs::nrfirstUpdate(double loglikelihood, std::vector<float> diffDist, std::vector<float> theta)
{

    fk_bfgs =  -loglikelihood;

    for (int jin=0; jin<theta.size(); jin++)
      gk_bfgs(jin) =  -diffDist[jin];

    for (int jin=0; jin<theta.size(); jin++)
      xk_bfgs(jin) = theta[jin];

    // break;  // break out of inner loop so i can compute dk
    flag_a = 0;

    vec_size = theta.size();

    stpmax = stepmaxconst*dmaxfn(sqrt(xk_bfgs.squaredNorm()), (double)vec_size);
}


void bfgs::nrStartNewDescent(std::vector<float> &theta)
{

  double sum = 0.0;
  double test = 0.0, temp = 0.0;
  int i;


  const double TOLX = 1.0e-7; // Convergence criterion on delta x.

  dk_bfgs = -1 * (Hk_bfgs * gk_bfgs);

  vec_size = theta.size();

  sum = sqrt(dk_bfgs.squaredNorm());

  if (sum > stpmax)
    for (i=0; i < vec_size; i++)
      dk_bfgs(i) *= stpmax/sum; // Scale if attempted step is too big.

  slope = gk_bfgs.transpose() * dk_bfgs;

  if (slope > 0.0)
  {
    std::cout << "Round off problem in lnsrch.";
    flag_a = 100;
    return;
  } else if (slope == 0.0) // found solution because g = 0 hence slope = g*d = 0
  {
    flag_a = 100;
    return;
  }

  test = 0.0; // Compute λmin.

  for (i=0; i<vec_size; i++)
  {
    temp = fabs(dk_bfgs(i)) / dmaxfn(fabs(xk_bfgs(i)), 1.0);
    if (temp > test)
      test = temp;
  }

  alamin = TOLX / test;
  resetalphaTo1(); // Always try full Newton step ﬁrst.


  xknew_bfgs = xk_bfgs + alpha* dk_bfgs;

  for (int jin=0; jin<theta.size(); jin++)
    theta[jin] = xknew_bfgs(jin);

}




int bfgs::nrlinesearchUpdate(double loglikelihood, std::vector<float> diffDist, std::vector<float> &theta)
{
  double alam2, tmplam, rhs1, rhs2;

  flag_a = 0;

  fknew_bfgs =  - loglikelihood;

  for (int jin=0; jin<theta.size(); jin++)
    gknew_bfgs(jin) =  -diffDist[jin];

  if (alpha < alamin)
  { //Convergence on delta x. For zero ﬁnding, the calling program should verify the convergence.
    for (int jin=0; jin<theta.size(); jin++)
      theta[jin] = xk_bfgs(jin);
    // this will break out of inner loop
    // in addition we want to break out of outer loop
    flag_a = 100;
    return (0);
  }
  else if (fknew_bfgs <= fk_bfgs + beta_dash*alpha*slope)
    return (0); // Sufficient function decrease.
  else {
    // Backtrack. Find new value of alpha
    flag_a = 1;

    if(alpha==1.0)
      tmplam = -slope/(2.0*(fknew_bfgs-fk_bfgs-slope)); // First time.
    else
    { // Subsequent backtracks.
      rhs1 = fknew_bfgs - fk_bfgs - alpha*slope;
      rhs2 = f2 - fk_bfgs - alam2*slope;
      double a = (rhs1/(alpha*alpha) - rhs2/(alam2*alam2))/(alpha-alam2);
      double b = (-alam2*rhs1/(alpha*alpha) + alpha*rhs2/(alam2*alam2))/(alpha-alam2);
      if (a == 0.0)
        tmplam=-slope/(2.0*b);
      else {
        double disc = b*b - 3.0*a*slope;
        if (disc<0.0)
          tmplam = 0.5*alpha;
        else if (b<=0.0)
          tmplam = (-b+sqrt(disc))/(3.0*a);
        else
          tmplam = -slope/(b+sqrt(disc));
      }
      if(tmplam > 0.5*alpha)
        tmplam = 0.5*alpha; // alpha <= 0.5 alpha.
    }
  }
  alam2 = alpha;
  f2=fknew_bfgs;
  alpha = dmaxfn(tmplam, 0.1*alpha); // alpha >= 0.1*alpha.

  if (alpha < 1e-40)
    flag_a = 100;
  else if (flag_a == 1)
  {
    xknew_bfgs = xk_bfgs + alpha*dk_bfgs;

    for (int jin=0; jin<theta.size(); jin++)
      theta[jin] = xknew_bfgs(jin);

  }

  return (1);

}




int bfgs::nrouteriterUpdate(std::vector<float> &theta)
{

  int another_flag = 1;

  vec_size = theta.size();

  Eigen::VectorXd sk_bfgs(vec_size), yk_bfgs(vec_size);
  double rhok_bfgs;
  Eigen::MatrixXd Vk_bfgs(vec_size, vec_size), Hkp1_bfgs(vec_size, vec_size);
  Eigen::MatrixXd Hkp1_bfgs_temp(vec_size, vec_size);

  // % find new Hkp1
  sk_bfgs = xknew_bfgs - xk_bfgs;
  yk_bfgs = gknew_bfgs - gk_bfgs;

  rhok_bfgs = 1/(yk_bfgs.transpose() * sk_bfgs);
  Vk_bfgs = eye - rhok_bfgs * (yk_bfgs * sk_bfgs.transpose());
  Hkp1_bfgs_temp = Vk_bfgs * Hk_bfgs;
  Hkp1_bfgs = Hkp1_bfgs_temp * Vk_bfgs + rhok_bfgs*(sk_bfgs * sk_bfgs.transpose());

  xk_bfgs = xknew_bfgs;
  gk_bfgs = gknew_bfgs;
  Hk_bfgs = Hkp1_bfgs;
  fk_bfgs = fknew_bfgs;

  // starting new search
  resetalphaTo1();
  resetflagaTo1();

  nrStartNewDescent(theta);

  /// this stop criteria might not be necesssary
  double normxk = xk_bfgs.squaredNorm();
  double normgk = gk_bfgs.squaredNorm();

  double max1 = 1;

  if (normxk > max1)
    max1 = normxk;
  another_flag = 0;
  if (normgk >= 1e-10 * max1)
    another_flag = 1;

  return another_flag;
}




void bfgs::nrprodHkgklimited(int iterout)
{
    if (iterout <= lmemory)
    {
      incr = 0;
      bound = iterout;
    } else {
      incr = iterout - lmemory;
      bound = lmemory;
    }


    if (iterout == 0)
    {
      for (int jin=0; jin<vec_size; jin++)
        q(0, jin) = gk_bfgs(jin);

      for (int jin=0; jin<vec_size; jin++)
        r(0,jin) = q(0, jin);
    } else {

      for (int jin=0; jin<vec_size; jin++)
        q(bound, jin) = gknew_bfgs(jin);

      // std::cout <<" s size " << s.size() << "\n";
      // std::cout <<" y size " << y.size() << "\n";
      // std::cout <<" rho size " << rho.size() << "\n";
      // std::cout <<" lalpham " << lalpham.size() << "\n";

      // since i am using a queue to store s, y, rho we don't use j but instead use a variable jj
      // only for clarity (though the value is the same as i)
      for (int i=bound-1; i>=0; i--)
      {
        int j = i + incr;
        int jj = i;
        // alpha_i = rho_j * s_j^T * q_{i+1}
        double stimesq = 0.0;
        for (int t=0; t<vec_size; t++)
          stimesq += s[jj](t) * q(i+1,t);

        lalpham(i) = rho[jj]*stimesq;

        // q_i = q_{i+1} - alpha_i y_j
        for (int jin=0; jin<vec_size; jin++)
          q(i, jin) = q(i+1, jin) - lalpham(i)*y[jj](jin);
      }

      for (int jin=0; jin<vec_size; jin++)
        r(0,jin) = q(0, jin);  // since H0 is the identity matrix

      for (int i=0; i <= bound-1; i++)
      {
        int j = i + incr;
        int jj = i;
        // beta_j = rho_j * y_j^T * r_i
        double ytimesr = 0.0;
        for (int t=0; t<vec_size; t++)
          ytimesr += y[jj](t) * r(i,t);

        // the paper says lbetam(j) which is confusing and should probably be indexed by i in the following equation
        lbetam(jj) = rho[jj] * ytimesr;

        // r_{i+1} = r_i + s_j*(alpha_i - beta_i)
        for (int t=0; t<vec_size; t++)
          r(i+1,t) = r(i,t) + s[jj](t)*(lalpham(i) - lbetam(i));
      }
    }
}


void bfgs::nrupdatedklimited(int iterout)
{
  nrprodHkgklimited(iterout); // this also sets bound useful in call later

  for (int i=0; i<vec_size; i++)
    dk_bfgs(i) =-1 * r(bound,i);
}


void bfgs::nrStartNewDescentlimited(std::vector<float> &theta)
{
  // dk should have been updated before this function call.

  double sum = 0.0;
  double test = 0.0, temp = 0.0;
  int i;

  const double TOLX = 1.0e-7; // Convergence criterion on delta x.

  sum = sqrt(dk_bfgs.squaredNorm());

  if (sum > stpmax)
    for (i=0; i < vec_size; i++)
      dk_bfgs(i) *= stpmax/sum; // Scale if attempted step is too big.

  slope = gk_bfgs.transpose() * dk_bfgs;

  if (slope >= 0.0)
  {
    std::cout << "Round off problem in lnsrch. \n";
    flag_a = 100;
    return;
  } else if (slope == 0.0) // found solution because g = 0 hence slope = g*d = 0, but what is d = 0 ??
  {
    flag_a = 100;
    return;
  }

  test = 0.0; // Compute λmin.

  for (i=0; i<vec_size; i++)
  {
    temp = fabs(dk_bfgs(i)) / dmaxfn(fabs(xk_bfgs(i)), 1.0);
    if (temp > test)
      test = temp;
  }

  alamin = TOLX / test;
  resetalphaTo1(); // Always try full Newton step ﬁrst.

  xknew_bfgs = xk_bfgs + alpha* dk_bfgs;

  for (int jin=0; jin<theta.size(); jin++)
    theta[jin] = xknew_bfgs(jin);

}


int bfgs::nrouteriterUpdatelimited(std::vector<float> &theta, int iterout)
{

  int another_flag = 1;

  Eigen::VectorXd sk_bfgs(vec_size), yk_bfgs(vec_size);
  sk_bfgs = xknew_bfgs - xk_bfgs;
  yk_bfgs = gknew_bfgs - gk_bfgs;
  double rhok_bfgs = 1/(yk_bfgs.transpose() * sk_bfgs);

  // still need to update rho, s, y

  if (iterout > lmemory)
  {
    if (y.size() == 0)
    {
      std::cout << "deque cannot be empty \n";
      exit(1);
    }
    y.pop_front();
    s.pop_front();
    rho.pop_front();
  }
  y.push_back(yk_bfgs);
  s.push_back(sk_bfgs);
  rho.push_back(rhok_bfgs);


  // computing d_1 = - H_1 * g_1
  nrprodHkgklimited(iterout);

  // std::cout << "bound = " << bound <<" \n ";

  for (int i=0; i<vec_size; i++)
    dk_bfgs(i) = -1 * r(bound,i);


  xk_bfgs = xknew_bfgs;
  gk_bfgs = gknew_bfgs;
  fk_bfgs = fknew_bfgs;

  //stpmax = stepmaxconst*dmaxfn(sqrt(xk_bfgs.squaredNorm()), (double)vec_size);

  // starting new search
  resetalphaTo1();
  resetflagaTo1();

  nrStartNewDescentlimited(theta);

  /// this stop criteria might not be necesssary
  double normxk = xk_bfgs.squaredNorm();
  double normgk = gk_bfgs.squaredNorm();

  double max1 = 1;

  if (normxk > max1)
    max1 = normxk;
  another_flag = 0;
  if (normgk >= 1e-10 * max1)
    another_flag = 1;

  return another_flag;
}












void bfgs::resetalphaTo1()
{
	alpha = alphastart;
}

void bfgs::resetflagaTo1()
{
	flag_a = 1;
}

int bfgs::getflaga()
{
  return(flag_a);
}

void bfgs::setbetadash(double value)
{
  beta_dash = value;
}

// Explanation from numerical recipes in C++ (second edition) page 385
// Given an n-dimensional point xold[1..n], the value of the function and gradient there, fold and
// g[1..n], and a direction p[1..n], ﬁnds a newpoint x[1..n] along the direction p from xold where the
// function funch as decreased “sufficiently.” The new function value is returned inf. stpmax is an
// input quantity that limits the length of the steps so that you do not try to evaluate the function in
// regions where it is undeﬁned or subject to overﬂow. p is usually the Newton direction. The output
// quantity check is false(0) on a normal exit. It is true(1) when x is too close to xold. In a
// minimization algorithm, this usually signals convergence and can be ignored. However, in a zero-ﬁnding
// algorithm the calling program should check whether the convergence is spurious.
// Some “difficult” problems may require double precision in this routine.




