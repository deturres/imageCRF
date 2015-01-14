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


int main(int argc, char* argv[])
{
  // testfunction();
  std::vector<float> theta(2);
  std::vector<float> g(2);
  float fx = 0.0;
  
  int bfgs_type = 2; // 1 for bfgs and 2 for lbfgs
  int memory = 2; 
  
  // initialize theta
  theta[0] = 22.4;
  theta[1] = 34.6;
  
  double stepconst = 100.0;
  bfgs bfgsobj(bfgs_type, theta, 0.1, 0.75, 1.0, memory, stepconst);
  
  int another_flag = 1;
  int flag_a = 0;
  
  for (int iterout=0; iterout < 200 && another_flag==1; iterout++)
  {
    bfgsobj.resetalphaTo1();
    bfgsobj.resetflagaTo1();
    
    for (int iter=0; iter < 100 && bfgsobj.getflaga() == 1; iter++)
    {
      std::cout << " Outer: " << iterout << " Inner: " <<iter <<std::endl;
      paraboloid(theta, fx, g);
      if (iterout==0 && iter==0)
      {
        bfgsobj.nrfirstUpdate(fx, g, theta);
      }
      else  // otherwise you are now doing linesearch and follow procedure
      {
        another_flag = bfgsobj.nrlinesearchUpdate(fx, g, theta); // theta gets updated
      }

      std::cout << theta[0] << "  " << theta[1] << "  " << fx << std::endl;
    }
    
    // bfgs processing
    if (bfgsobj.getflaga() == 100)
    {
      std::cout << "reached small alpha value " << std::endl;
      break;
    }
    
    if (bfgs_type == 2) // lbfgs
    {
      if (iterout==0) 
		  {
        bfgsobj.nrupdatedklimited(iterout);
			  bfgsobj.nrStartNewDescentlimited(theta); // update theta			
		  } else
		  {
			  another_flag = bfgsobj.nrouteriterUpdatelimited(theta, iterout);  // update theta
        // the above call calls nrStartNewDescentlimited for the correct iterout
		  }
    } else if (bfgs_type == 1) // bfgs
    {
      if (iterout==0)
		  {
			  bfgsobj.nrStartNewDescent(theta); // update theta			
		  } else
		  {
			  another_flag = bfgsobj.nrouteriterUpdate(theta);  // update theta
		  }
    }

    // std::cout << " fx value = " << fx << "\n";
    std::cout << theta[0] << "  " << theta[1] << "  " << fx << std::endl;
  }
  
  paraboloid(theta, fx, g);
  std::cout << theta[0] << "  " << theta[1] << "  " << fx << std::endl;
  
  
  return 0;
}	
		
