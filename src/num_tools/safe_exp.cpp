#include <math.h>
#include <algorithm>

#include "num_tools.hpp"

// Cut-off parameter for exponential evaluation
// Warning: this value was tuned on the previous version of the code
//          (with C++ math functions) 
const double fExpUppLim = 120.;          // possible alternative choice: const*DBL_MAX_EXP
const double fExpLowLim = -1.0E+10; //-fExpUppLim;   //                              const*DBL_MIN_EXP

// Safe exp function to avoid underflow/overflow
double safe_FDexp(double x) {
	//double arg = std::min(std::max(x,fExpLowLim),fExpUppLim);
    double arg = (x <= 120.) ? x : 120.;
	if (arg < 20.) {
		return 1. - exp(arg) / (1. + exp(arg));	
	} else {
		return 1. / (1. + exp(arg));
	}
}




