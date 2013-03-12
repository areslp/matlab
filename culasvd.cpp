#include "mex.h"
#include "cula_lapack.hpp"
#include "culamex.hpp"
#include "iostream"
#include "windows.h"

using std::max;
using std::min;
using std::cout;
using std::endl;

// Complex conjugation of complex data
template<class T> void Conjugate(T* a) { /* Do nothing */ };
template<> void Conjugate(culaFloatComplex* a) { a->y = -(a->y); }
template<> void Conjugate(culaDoubleComplex* a) { a->y = -(a->y); }

 
// Templated wrapper
template <typename T>
void mexCulaGesvd(int nlhs,           /* number of expected outputs */
              mxArray* plhs[],        /* output pointer array */
              int nrhs,               /* number of inputs */
              const mxArray* prhs[],  /* input pointer array */
              mxClassID id,
              mxComplexity complexity)
{
    unsigned long t1,t2;
    // Function core, details in "Using CULA in MATLAB, Part 2"
    // culaGesvd(...)
    // Initialize flags and types
    typedef typename ToReal<T>::type RealType;
    bool isReal = (complexity == mxREAL);
    bool isComplex = (complexity == mxCOMPLEX);
 
    // Initialize sizes
    int M = (int) mxGetM(prhs[0]);
    int N = (int) mxGetN(prhs[0]);
    int K = min(M,N);
    int L = max(M,N);
 
    t1=GetTickCount();
    // Allocate a temporary to not destroy input data
    T* A = (T*) mxMalloc( M * N * sizeof(T) );
 
    if (isReal)
    {
        // Copy input data directly into temporary
        memcpy( A, mxGetPr( prhs[0] ), M * N * sizeof(T) );
    }
    else if (isComplex)
    {
        // If complex, convert from MATLAB format into CULA format
        MatToCula( A, prhs[0] );
    }
 
    // Create MATLAB output matrices
    plhs[0] = mxCreateNumericMatrix(M, M, id, complexity);  // U (M x M)
    plhs[1] = mxCreateNumericMatrix(M, N, id, 0);           // S (M x N, Real)
    plhs[2] = mxCreateNumericMatrix(N, N, id, complexity);  // V (N x N)
 
    // Allocate CULA intermediate
    RealType* SVEC = (RealType*) mxMalloc( K * sizeof(RealType) );
 
    // CULA Memory Pointers
    T* U;
    T* VT;
 
    if (isReal)
    {
        // Get CULA memory pointers from allocated MATLAB matrices
        U = (T*) mxGetPr( plhs[0] );
        VT = (T*) mxGetPr( plhs[2] );
    }
    else if (isComplex)
    {
        // If complex, allocate an AoS complex buffer for CULA
        U = (T*) mxMalloc( M * M * sizeof(T) );
        VT = (T*) mxMalloc( N * N * sizeof(T) );
    }
    t2=GetTickCount();
    cout<<"memory preparation takes:"<<t2-t1<<endl;
    // culaStatus status;
    // Initialize CULA
    t1=GetTickCount();
    culaStatus status = culaInitialize();
    checkStatus(status, "culaInitialize");
    t2=GetTickCount();
    cout<<"culaInitialize takes:"<<t2-t1<<endl;

 
    // CULA SVD Factorization
    t1=GetTickCount();
    status = culaGesvd('A', 'A', M, N, A, M, SVEC, U, M, VT, N);
    checkStatus(status, "culaGesvd");
    t2=GetTickCount();
    cout<<"culaGesvd takes:"<<t2-t1<<endl;
 
    // Shutdown CULA
    t1=GetTickCount();
    culaShutdown();
    t2=GetTickCount();
    cout<<"culaShutdown takes:"<<t2-t1<<endl;
 
    // Get pointer to output matrix, S
    t1=GetTickCount();
    RealType* S = (RealType*) mxGetPr( plhs[1] );
 
    // Copy SVEC to diagonal of S
    for (int i=0; i<K; i++)
        S[i*M+i] = SVEC[i];
 
    // Inplace transpose of VT
    for (int i=0; i<N; i++)
    {
        for (int j=i; j<N; j++)
        {
            T temp = VT[j+i*N];
            VT[j+i*N] = VT[i+j*N];
            VT[i+j*N] = temp;
        }
    }
 
    // If complex, conjugate VT
    if (isComplex)
        for (int i=0; i<N; i++)
            for (int j=0; j<N; j++)
                Conjugate( &VT[j+i*N] );
 
    if (isComplex)
    {
        // If complex, convert from CULA format into MATLAB format
        CulaToMat( plhs[0], U );
        CulaToMat( plhs[2], VT );
 
        // Free MATLAB buffers
        mxFree(U);
        mxFree(VT);
    }
 
    // Free allocate data
    mxFree(A);
    mxFree(SVEC);
    t2=GetTickCount();
    cout<<"output and clear takes:"<<t2-t1<<endl;

}
 
// MATLAB Gateway Function
void mexFunction(int nlhs,              /* number of expected outputs */
                 mxArray* plhs[],       /* output pointer array */
                 int nrhs,              /* number of inputs */
                 const mxArray* prhs[]  /* input pointer array */ )
{
    // We only support a full SVD in this example
    if(nrhs != 1)
        mexErrMsgTxt("culasvd: Must have 1 input argument [X]");
    if(nlhs != 3)
        mexErrMsgTxt("culasvd: Must have 3 output arguments [U,S,V]");
 
    // Get precision (single or double)
    mxClassID classID = mxGetClassID(prhs[0]);
 
    // Get complexity (real or complex)
    mxComplexity complexity = mxIsComplex(prhs[0]) ? mxCOMPLEX : mxREAL;
 
    // Switch based on data type
    if ( classID == mxSINGLE_CLASS && complexity == mxREAL)
        mexCulaGesvd<culaFloat>(nlhs, plhs, nrhs, prhs, classID, complexity);
    else if (classID == mxDOUBLE_CLASS && complexity == mxREAL )
        mexCulaGesvd<culaDouble>(nlhs, plhs, nrhs, prhs, classID, complexity);
    else if ( classID == mxSINGLE_CLASS && complexity == mxCOMPLEX )
        mexCulaGesvd<culaFloatComplex>(nlhs, plhs, nrhs, prhs, classID, complexity);
    else if ( classID == mxDOUBLE_CLASS && complexity == mxCOMPLEX )
        mexCulaGesvd<culaDoubleComplex>(nlhs, plhs, nrhs, prhs, classID, complexity);
    else
        mexErrMsgTxt("culasvd: Unknown or unsupported data type");
}

