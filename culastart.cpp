#include "mex.h"
#include "cula_lapack.hpp"
#include "culamex.hpp"
#include "iostream"
 
// MATLAB Gateway Function
void mexFunction(int nlhs,              /* number of expected outputs */
                 mxArray* plhs[],       /* output pointer array */
                 int nrhs,              /* number of inputs */
                 const mxArray* prhs[]  /* input pointer array */ )
{
    // Initialize CULA
    culaStatus status = culaInitialize();
    checkStatus(status, "culaInitialize");
}

