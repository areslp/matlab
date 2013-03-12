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
void collectdata(double *y, int* edgerow, int *edgecol, double *edgesdist, int *relative,  int d, int n, int ks, double *a, double *b, double *g);

/* auxiliary functions */
void sanity_check(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/* the gateway */ 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int n, ks=1,  d;
    double *x, *y, *a, *b, *g;
    int *neighbors;
    
    /* sanity_check(nlhs, plhs, nrhs, prhs); */
    

    d = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    /*printf("%d data points, reducing from dimension %d ---> dimension %d, using \
%d nearest neighbors\n", n, D, d, ks); */
    plhs[0] = mxCreateDoubleMatrix(d*d, d*d, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(d*d, n, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n,1, mxREAL);
    
    collectdata( (double*)mxGetData(prhs[0]),
                 (int*)mxGetData(prhs[1]),
                 (int*)mxGetData(prhs[2]),
                 (double*)mxGetData(prhs[3]),
                 (int*)mxGetData(prhs[4]),
                  d, n, ks, 
                 (double*)mxGetData(plhs[0]),
                 (double*)mxGetData(plhs[1]),
                 (double*)mxGetData(plhs[2]));
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
    int i;
    /*       dcopy_(&n, v1, &incx,v3, &incx);
	     daxpy_(&n, &none, v2, &incx, v3, &incx);*/
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



void collectdata(double *y, int* edgesrow, int *edgescol, double * edgesdist, int *relative,  int d, int n, int ks, double *a, double *b, double *g) 
{
    double  *diffy, *yij, *nid;
    double gij;
    int i, j, k;
    int nn_start, nn_end, nn;
    double *diffxjk, *diffyjk;
    double gjk, *yjk, *tempA;
    int relflg = *relative;

    /* clear data area */
    memset(a, 0, sizeof(double)*d*d*d*d);
    memset(b, 0, sizeof(double)*d*d*n);
    memset(g, 0, sizeof(double)*n);

    /* temporary working space */
    diffy = (double*)calloc(d, sizeof(double));
    yij = (double*)calloc(d*d, sizeof(double));
    
    tempA = (double*)calloc((d*d)*(d*d+1)/2, sizeof(double));
    
    if(!diffy || !yij || !tempA) {
        mexErrMsgTxt("Out of memory..cannot allocate working space..");
        return;
    }
    
    if(relflg==0) {
    
        for(i=0; i <n ; i++) { 
            /* figure out the nearest neighbors */
            nn_start = edgescol[i]; nn_end = edgescol[i+1]-1;
            if (nn_end > nn_start) {
                for(nn=nn_start; nn<=nn_end; nn++) {
                    j = edgesrow[nn];
            
                    /* compute diffx and diffy */
		    /*                    vecsub(x+i*D, x+j*D, diffx, D);*/
                    gij = edgesdist[nn];
                
                    vecsub(y+i*d, y+j*d, diffy, d);
                    vecout(diffy,  yij, d);
            
                    /* update A, b */
                    updateA(tempA, yij, 1.0,d*d);
            
                    updateB(b+i*d*d, gij, yij, d*d);
                    g[i] = g[i]+gij*gij;
                }   
            }
        }
    } else {
        /* reproducing the code above to some degree so that we don't have to 
         do if statement inside cache-intensive loop */
        for(i=0; i <n ; i++) { 
            /* figure out the nearest neighbors */
            nn_start = edgescol[i]; nn_end = edgescol[i+1]-1;
            if (nn_end > nn_start) {
                for(nn=nn_start; nn<=nn_end; nn++) {
                    j = edgesrow[nn];
            
		    gij=edgesdist[nn];
                
                    vecsub(y+i*d, y+j*d, diffy, d);
                    vecout(diffy,  yij, d);
		    /*		    printf("diff: %f %f %f\n",yij[0],gij,edgesdist[nn]);*/
                    /* update A, b */
                    updateA(tempA, yij, 1.0/(gij*gij), d*d);
            
                    updateB(b+i*d*d, 1/gij, yij, d*d);
                    g[i] = g[i]+1;
                }   
            }
        }
    }
    /* make A symmetric */
    recoverA(d*d, a, tempA); 
}   

    
