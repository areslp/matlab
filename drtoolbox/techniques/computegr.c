/*
 * =============================================================
 * computegr.c 
 *  
 * input: x,i,j,d
 *    x    : matrix DxN
 *    i,j  : indices of neighbors 1xC
 *    d    : distances 1xC
 *
 * output: 
 *    dx : matrix DxN
 * =============================================================
 */

/* $Revision: 1.1 $ */

#include "mex.h"

/* If you are using a compiler that equates NaN to zero, you must
 * compile this example using the flag -DNAN_EQUALS_ZERO. For 
 * example:ç
 *
 *     mex -DNAN_EQUALS_ZERO findnz.c  
 *
 * This will correctly define the IsNonZero macro for your
   compiler. */

#if NAN_EQUALS_ZERO
#define IsNonZero(d) ((d) != 0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d) != 0.0)
#endif

#define Xp prhs[0]
#define ip prhs[1]
#define jp prhs[2]
#define distp prhs[3]

double square(double x) { return(x*x);}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  /* Declare variables. */ 

  double *obj,objective,*dX, *X, *dist,weight,curdist, *ind1,*ind2,*v1,*v2,*vo, *dummy;
  double pev=1.0,mev=1.0,*ev,*av,*bv;
  int m,p,n,inds;
  int oldi1=-1,i1,i2,j,i,r,c;
  int ione=1;
  char *chu="U"; 
  char *chl="L";
  char *chn2="T";
  char *chn="N";
  double minusone=-1.0,one=1.0, zero=0.0;



  /* Check for proper number of input and output arguments. */    
  if (nrhs!=4) mexErrMsgTxt("Exactly one or two input arguments required.");
 

  if (nlhs > 2) {
    mexErrMsgTxt("Too many output arguments.");
  }

  /* Check data type of input argument. */
  if (!(mxIsDouble(prhs[0]))) {
   mexErrMsgTxt("Input array must be of type double.");
  }
      
  n = mxGetN(Xp);
  m = mxGetM(Xp);
  c = mxGetNumberOfElements(distp);

  /* Get the data. */
  X  = mxGetPr(Xp);
  ind1  = mxGetPr(ip);
  ind2  = mxGetPr(jp);
  dist  = mxGetPr(distp);

  /* Create output matrix */
  plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
  dX=mxGetPr(plhs[0]);
  plhs[1]=mxCreateDoubleMatrix(1,1,mxREAL);
  obj=mxGetPr(plhs[1]);

  /* create space for dummy vectors */
  dummy=malloc(m*sizeof(double));

  weight=0;
  oldi1=-1;
  v1=&X[0];v2=&X[0];vo=&dX[0];
  objective=0;
  for(i=0;i<c;i++){
    i1=(int) ind1[i]-1;
    i2=(int) ind2[i]-1;
    v2 =&X[(int) (i2*m)];

    if(i1!=oldi1){
      for(j=0;j<m;j++) {vo[j]=weight*v1[j]-dummy[j];dummy[j]=0;}
      weight=0;
     v1 =&X[(int) (i1*m)];
     vo=&dX[(int) (i1*m)];
     oldi1=i1;
    };
    curdist=0;
    for(j=0;j<m;j++) curdist=curdist+square(v1[j]-v2[j]);
    /*     curdist=(curdist-dist[i])/dist[i];objective+=square(curdist);*/
    curdist=(curdist-dist[i]);objective+=square(curdist);
    for(j=0;j<m;j++) dummy[j]+=curdist*v2[j];
    weight=weight+curdist;

  }

  for(j=0;j<m;j++) {vo[j]=weight*v1[j]-dummy[j];}

  obj[0]=objective;

  free(dummy);
}



