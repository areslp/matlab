// culamex.hpp
// Common helper functions for integrating CULA into MATLAB
 
#ifndef __CULAMEX_HPP__
#define __CULAMEX_HPP__

#include "string.h"
 
// Used to find real type associated with complex type
template<class T> struct ToReal { typedef T type; };
template<>        struct ToReal<culaFloatComplex> { typedef float type; };
template<>        struct ToReal<culaDoubleComplex> { typedef double type; };
 
// Convert from MATLAB complex format to CULA complex format
template <typename FloatType>
void MatToCula(FloatType* buf, const mxArray* src)
{
    typedef typename ToReal<FloatType>::type RealType;
 
    const RealType* r0 = (const RealType*) mxGetPr(src);
    const RealType* i0 = (const RealType*) mxGetPi(src);
 
    int M = (int)mxGetM(src);
    int N = (int)mxGetN(src);
 
    for(int j = 0; j < N; ++j)
    {
        for(int i = 0; i < M; ++i)
        {
            int p  = j*M+i;
            buf[p].x = r0[p];
            buf[p].y = i0[p];
        }
    }
}
 
// Convert from CULA complex format to MATLAB complex format
template <typename FloatType>
void CulaToMat(mxArray* src, FloatType* buf)
{
    typedef typename ToReal<FloatType>::type RealType;
 
    RealType* r0 = (RealType*) mxGetPr(src);
    RealType* i0 = (RealType*) mxGetPi(src);
 
    int M = mxGetM(src);
    int N = mxGetN(src);
 
    for(int j = 0; j < N; ++j)
    {
        for(int i = 0; i < M; ++i)
        {
            int p = j*M+i;
            r0[p] = buf[p].x;
            i0[p] = buf[p].y;
        }
    }
}
 
// Do nothing for non-complex cases
template <> void MatToCula(float* buf, const mxArray* src) {}
template <> void MatToCula(double* buf, const mxArray* src) {}
template <> void CulaToMat(mxArray* src, float* buf) {}
template <> void CulaToMat(mxArray* src, double* buf) {}
 
void checkStatus(culaStatus status, const char* funcname)
{
    if(!status)
        return;
 
    culaShutdown();
 
    char id[128];
    sprintf(id, "CULA:%s:", funcname);
 
    if(status == culaArgumentError)
    {
        strcat(id, "culaArgumentError");
        mexErrMsgIdAndTxt(id, "%s: Invalid value for parameter %d\n", funcname, culaGetErrorInfo());
    }
    else if(status == culaDataError)
    {
        strcat(id, "culaDataError");
        mexErrMsgIdAndTxt(id, "%s: Data error (%d)\n", funcname, culaGetErrorInfo());
    }
    else if(status == culaBlasError)
    {
        strcat(id, "culaBlasError");
        mexErrMsgIdAndTxt(id, "%s: Blas error (%d)\n", funcname, culaGetErrorInfo());
    }
    else if(status == culaRuntimeError)
    {
        strcat(id, "culaRuntimeError");
        mexErrMsgIdAndTxt(id, "%s: Runtime error (%d)\n", funcname, culaGetErrorInfo());
    }
    else if(status == culaNotInitialized)
        strcat(id, "culaNotInitialized");
    else if(status == culaNoHardware)
        strcat(id, "culaNoHardware");
    else if(status == culaInsufficientRuntime)
        strcat(id, "culaInsufficientRuntime");
    else if(status == culaInsufficientComputeCapability)
        strcat(id, "culaInsufficientComputeCapability");
    else if(status == culaInsufficientMemory)
        strcat(id, "culaInsufficientMemory");
    else if(status == culaFeatureNotImplemented)
        strcat(id, "culaFeatureNotImplemented");
    else
        strcat(id, "unknown");
 
    // Message that don't have error info fall through to here
    mexErrMsgIdAndTxt(id, "%s: %s\n", funcname, culaGetStatusString(status));
}


#endif //__CULAMEX_HPP__

