/*
 *
 * mexCCACollectdata.c
 *      prepare data for cca.m to form SDP
 *
 * by feisha@cis.upenn.edu
 */
 
 #include "mex.h"
 #include "matrix.h"
#include <stdlib.h>
#include <float.h>
#include <string.h> 
#include <math.h>

 /* the computation engine */
void collectdata(double *x, double *y, int* edgerow, int *edgecol, int *relative, 
int *nv, int *vidx, int D, int d, int n, int ks, double *a, double *b, double *g);

/* auxiliary functions */
void sanity_check(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/* the gateway */ 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int n, ks=1, D, d;
    double *x, *y, *a, *b, *g;
    int *neighbors;
    
    /* sanity_check(nlhs, plhs, nrhs, prhs); */
    
    n = mxGetN(prhs[0]);
    D = mxGetM(prhs[0]);
    d = mxGetM(prhs[1]);
    /*printf("%d data points, reducing from dimension %d ---> dimension %d, using \
%d nearest neighbors\n", n, D, d, ks); */
    
    plhs[0] = mxCreateDoubleMatrix(d*d, d*d, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(d*d, n, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n,1, mxREAL);
    
    collectdata( (double*)mxGetData(prhs[0]),
                 (double*)mxGetData(prhs[1]),
                 (int*)mxGetData(prhs[2]),
                 (int*)mxGetData(prhs[3]),
                 (int*)mxGetData(prhs[4]),
                 (int*)mxGetData(prhs[5]),
                 (int*)mxGetData(prhs[6]),
                 D, d, n, ks, 
                 (double*)mxGetData(plhs[0]),
                 (double*)mxGetData(plhs[1]),
                 (double*)mxGetData(plhs[2]));
    return;
}


/* check input arguments */
void sanity_check(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    int n1, n2, n3;
    /* Check for proper number of input and output arguments. */    
    if (nrhs != 4) {
        mexErrMsgTxt("Four input argument required.");
    } 
    if (nlhs != 3) {
        mexErrMsgTxt("Three output arguments required.");
    }

    /* Check data type of input argument. */
    /*if (!(mxIsDouble(prhs[0])) || !(mxIsDouble(prhs[1])) || !(mxIsDouble(prhs[2]))) {
        mexErrMsgTxt("Input array must be of type double.");
    }*/
    
    /* Check dimensions */
    if(mxGetNumberOfDimensions(prhs[0])!=2 || mxGetNumberOfDimensions(prhs[1])!=2 || mxGetNumberOfDimensions(prhs[2])!=2) {
        mexErrMsgTxt("Input arrays must be two-dimension arrays");
    }

    /* matching dimensions */
    n1 = mxGetN(prhs[0]);
    n2 = mxGetN(prhs[1]);
    n3 = mxGetN(prhs[2]);
    if (n1 != n2 || n1 !=n3) {
        mexErrMsgTxt("Input arrays must have same number of columns!");
    }
    return;
}

/* computing */
double vecdot(double* v1, int n) {  /* v1'*v2 */
    double x = 0;
    int i;
    for(i = 0; i < n; i++) {
        x += v1[i]*v1[i];
    }
    return x;
}
void vecsub(double* v1,double* v2, double* v3, int n) { /* v3=v1-v2 */
    int incx=1;
    double none = -1;
    /*       dcopy_(&n, v1, &incx,v3, &incx);
	     daxpy_(&n, &none, v2, &incx, v3, &incx);*/
    int i;
     for(i=0;i<n;i++) v3[i]=v1[i]-v2[i]; 
    
}
void vecout(double* v1, double* v3, int d) { /* v3 = v1*v2' */
  double temp,one = 1.0;
    int inca = 1;
    int i,j,c;
 
    /*   dger_(&d , &d, &one, v1, &inca,v1, &inca, v3, &d); */
    c=0;
    for(i=0;i<d;i++){ 
      temp=v1[i];
      for(j=0;j<d;j++){
   	v3[c]=temp*v1[j];
      c++;
	}    
    }
}
void makevecsym(double *v, int d) { /* assume V is a dxd matrix, return (V+V')/2 */
    int row, col, i;
    double temp;
    for(col=0; col < d; col++) {
        for(row=col; row < d; row++) {
            temp = (v[row*d+col] + v[col*d+row])/2;
            v[row*d+col] = v[col*d+row] = temp;
        }
    }
}
void updateA(double *a, double *v, double alpha, int d) {  /* update A */
    int inca = 1;
    char *uplo= "U";
    double temp;
    int i,j,c;

    /*            dspr_(uplo, &d , &alpha, v, &inca, a);*/
 
    if(alpha==1.0){ /*optimize the common case */
     c=0;
     for(i=0;i<d;i++){
      temp=v[i];
      for(j=0;j<=i;j++)  {
	a[c]+=temp*v[j];
       c++;
      }
      }
    } else{
     c=0;
     for(i=0;i<d;i++){
      temp=v[i]*alpha;
      for(j=0;j<=i;j++)  {
	a[c]+=temp*v[j];
       c++;
      }
      }
    }
}
void updateB(double * b, double alpha, double *v, int d) { /* update B(:,i) */
    int i;
    int incx = 1;
    /*   daxpy_(&d, &alpha, v, &incx, b, &incx);*/
    for(i=0;i<d;i++) b[i]+=alpha*v[i];
}
/* transform upper triangular matrix stored in A as full matrix */
void recoverA(int D, double *A, double *tempA)
{
  int l = 0, idx, jdx,s, incx=1, incy=1;
  int n = D*(D+1)/2;
  for(jdx=0; jdx <D; jdx++) {
     for(idx=0; idx <=jdx; idx++) {
    	A[jdx*D+idx] = tempA[l];
        l++;    
      }
  }
  for(jdx=0; jdx <D; jdx++) 
    for(idx=jdx; idx <D; idx++) {
        A[jdx*D+idx] = A[idx*D+jdx];
    }
 

}



void collectdata(double *x, double *y, int* edgesrow, int *edgescol, int *relative, int *nv, 
int *vidx, int D, int d, int n, int ks, double *a, double *b, double *g) 
{
    double *diffx, *diffy, *yij, *nid;
    double gij;
    int i, j, k, iEdge, iVertex,itemp, ii;
    int nn_start, nn_end, nn;
    double *diffxjk, *diffyjk;
    double gjk, *yjk, *tempA;
    int relflg = *relative;
    /* clear data area */
    memset(a, 0, sizeof(double)*d*d*d*d);
    memset(b, 0, sizeof(double)*d*d*n);
    memset(g, 0, sizeof(double)*n);

    /* temporary working space */
    diffx = (double*)calloc(D, sizeof(double));
    diffy = (double*)calloc(d, sizeof(double));
    yij = (double*)calloc(d*d, sizeof(double));
    
    
    tempA = (double*)calloc((d*d)*(d*d+1)/2, sizeof(double));
    
    if(!diffx || !diffy || !yij || !tempA) {
        mexErrMsgTxt("Out of memory..cannot allocate working space..");
        return;
    }
    iEdge = 0;
    iVertex = 0;
    if(relflg==0) {
    
        for(i=0; i <n ; i++) { 
            /* figure out the nearest neighbors */
            nn_start = edgescol[i]; nn_end = edgescol[i+1]-1;
            if (nn_end >= nn_start) {
                for(nn=nn_start; nn<=nn_end; nn++) {
                    j = edgesrow[nn];
                    
                    /* compute diffx and diffy */
                    vecsub(x+i*D, x+j*D, diffx, D);
                    gij = vecdot(diffx, D);
                
                    vecsub(y+i*d, y+j*d, diffy, d);
                    vecout(diffy,  yij, d);
            
                    /* update A, b */
                    updateA(tempA, yij, 1.0*nv[iEdge],d*d);
                    for(itemp=0; itemp < nv[iEdge]; itemp++) {
                        ii = vidx[iVertex]-1;
                        updateB(b+ii*d*d, gij, yij, d*d);
                        g[ii] = g[ii]+gij*gij;
                        iVertex++;
                    }
                    
                    iEdge++;
                }   
            }
        }
        printf("Iedge :%d, ivertex: %d\n", iEdge, iVertex);
    } else {
        /* reproducing the code above to some degree so that we don't have to 
         do if statement inside cache-intensive loop */
        for(i=0; i <n ; i++) { 
            /* figure out the nearest neighbors */
            nn_start = edgescol[i]; nn_end = edgescol[i+1]-1;
            if (nn_end >= nn_start) {
                for(nn=nn_start; nn<=nn_end; nn++) {
                    j = edgesrow[nn];
            
                    /* compute diffx and diffy */
                    vecsub(x+i*D, x+j*D, diffx, D);
                    gij = vecdot(diffx, D);
                
                    vecsub(y+i*d, y+j*d, diffy, d);
                    vecout(diffy,  yij, d);
            
                    /* update A, b */
                    updateA(tempA, yij, 1*nv[iEdge]/(gij*gij), d*d);
                    for(itemp=0; itemp < nv[iEdge]; itemp++) {
                        ii = vidx[iVertex]-1;
                        updateB(b+ii*d*d, 1/gij, yij, d*d);
                        g[ii] = g[ii]+1;
                        iVertex++;
                    }
                    iEdge++;
                }   
            }
        }
    }
    /* make A symmetric */
    recoverA(d*d, a, tempA); 
}   


  