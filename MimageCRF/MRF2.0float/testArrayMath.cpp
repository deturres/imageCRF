#include <iostream>
#include "ArrayMath.h"

int main(int argc, char ** argv)
{

  FloatType array[2] = {0,0};

  std::cout << "logSum(array) = " << ArrayMath::logSum(array,2) << "\n";
}
