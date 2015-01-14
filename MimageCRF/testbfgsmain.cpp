/*
 *  main.cpp
 *  
 *
 *  Created by Chetan Bhole on 2/22/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "bfgs.h"
#include <vector>
#include <iostream>

void paraboloid(const std::vector<float>& x, float &v, std::vector<float> &g)
{
  float A[2] = {1, 10};
  float B[2] = {3, 7};
  v = 0.0;
  for (int i=0; i<2; i++)
    v += -A[i]*(x[i]-B[i])*(x[i]-B[i]);

  for (int i=0; i<2; i++)
    g[i] = -2*A[i]*(x[i]-B[i]);

}

void testfunction()
{
  std::vector<float> theta(2);
  std::vector<float> g(2);
  
  float fx = 0.0;
  
  theta[0] = 3.4;
  theta[1] = 5.6;
  paraboloid(theta, fx, g);
  
  std::cout << " fval = " << fx <<std::endl;
  for (int i=0; i<2; i++)
    std::cout << " " << g[i];
  std::cout << std::endl;
}


int testnrbfgs()
{
  // testfunction();
  std::vector<float> theta(2);
  std::vector<float> g(2);
  float fx = 0.0;
  
  // initialize theta
  theta[0] = 22.4;
  theta[1] = 34.6;
  
  bfgs bfgsobj(1, theta, 0.25, 0.75, 1.0);
  
  int another_flag = 1;
  int flag_a = 0;
  
  for (int iterout=0; iterout < 50 && another_flag==1; iterout++)
  {
    bfgsobj.resetalphaTo1();
    bfgsobj.resetflagaTo1();
    
    for (int iter=0; iter < 100 && bfgsobj.getflaga() == 1; iter++)
    {
      std::cout << " Outer: " << iterout << " Inner: " <<iter <<std::endl;
      paraboloid(theta, fx, g);
      if (iterout==0 && iter==0)
        bfgsobj.nrfirstUpdate(fx, g, theta);
      else  // otherwise you are now doing linesearch and follow procedure
      {
        another_flag = bfgsobj.nrlinesearchUpdate(fx, g, theta); // theta gets updated
      }
    }
    
    // bfgs processing
    if (bfgsobj.getflaga() == 100)
    {
      std::cout << "reached small alpha value " << std::endl;
      break;
    }
    
    if (iterout==0)
		{
			bfgsobj.nrStartNewDescent(theta); // update theta			
		} else
		{
			another_flag = bfgsobj.nrouteriterUpdate(theta);  // update theta
		}  
    
  }
  
  paraboloid(theta, fx, g);
  std::cout << theta[0] << "  " << theta[1] << "  " << fx << std::endl;

  return (0);
}	


int testbfgs()
{
  // testfunction();
  std::vector<float> theta(2);
  std::vector<float> g(2);
  float fx = 0.0;
  
  // initialize theta
  theta[0] = 22.4;
  theta[1] = 34.6;
  
  bfgs bfgsobj(1, theta, 0.25, 0.75, 1.0);
  
  int another_flag = 1;
  int flag_a = 0;
  
  for (int iterout=0; iterout < 50 && another_flag==1; iterout++)
  {
    bfgsobj.resetalphaTo1();
    bfgsobj.resetflagaTo1();
    
    for (int iter=0; iter < 100 && bfgsobj.getflaga() == 1; iter++)
    {
      std::cout << " Outer: " << iterout << " Inner: " <<iter <<std::endl;
      paraboloid(theta, fx, g);
      if (iterout==0 && iter==0)
        bfgsobj.firstUpdate(fx, g, theta);
      else  // otherwise you are now doing linesearch and follow procedure
      {
        bfgsobj.linesearchUpdate(fx, g, theta); // theta gets updated
      }
    }
    
    // bfgs processing
    if (bfgsobj.getflaga() == 100)
    {
      std::cout << "reached small alpha value " << std::endl;
      break;
    }
    
    if (iterout==0)
		{
			bfgsobj.firstOuteriterUpdate(theta); // update theta			
		} else
		{
			another_flag = bfgsobj.outeriterUpdate(theta);  // update theta
		}  
    
  }
  
  paraboloid(theta, fx, g);
  std::cout << theta[0] << "  " << theta[1] << "  " << fx << std::endl;
    
  return 0;
}	




int main(int argc, char* argv[])
{

  std::cout << " Testing bfgs where linesearch has alpha decreasing by some constant amount \n"; 
  testbfgs();
  
  std::cout << " Testing bfgs where linesearch is from numerical recipes\n";
  testnrbfgs();

  return 0;
}
		
