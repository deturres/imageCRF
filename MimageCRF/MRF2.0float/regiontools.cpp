#include "regiontools.h"
#include <math.h>

FLOATTYPE sumNegLogProb(FLOATTYPE a, FLOATTYPE b)
		{
			if (a == 9999 && b == 9999)
				return 9999;
			else if (a > b)
				return b - log (1 + exp(b-a));
			else
				return a - log (1 + exp(a-b));
		}

FLOATTYPE sumLogProb(FLOATTYPE a, FLOATTYPE b)
		{
			if (a == 9999 && b == 9999)
				return 9999;
			else if (b > a)
				return b + log (1 + exp(a-b));
			else
				return a + log (1 + exp(b-a));
		}