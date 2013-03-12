#include "mex.h"
#include <math.h>


void mexFunction(
        int nlhs,       mxArray *plhs[],
        int nrhs, const mxArray *prhs[]
        )
{
  
    /* Declare variables. */
    int *irs, *jcs, *k_ind, i, j, k, l, m, n, p, cc, ck, nzmax, min_ind, *ind;
    double *d, *k_d, *pr, *sr, dist, tmp, min_d;
    mxArray *rhs[1], *lhs[2];

    /* Check for proper number of input and output arguments. */    
    if (nrhs < 1) {
        mexErrMsgTxt("At least one input argument required.");
    } 
    if (nrhs > 2) {
        mexErrMsgTxt("No more than two input arguments allowed.");
    } 
    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }  
    if (!(mxIsDouble(prhs[0]))) {
        mexErrMsgTxt("Input argument must be of type double.");
    } 

    /* Get all input data. */
    m  = mxGetM(prhs[0]);                       // height
    n  = mxGetN(prhs[0]);                       // width
    pr = mxGetPr(prhs[0]);                      // pointer to data array (double)
    if (nrhs > 1) {
        k = (int)mxGetScalar(prhs[1]);          // number of neighbors
    }
    else {
        k = 12;                                 // number of neighbors
    }

    /* Allocate space for sparse matrix (don't allow complex values). */
    nzmax = m * k;
    plhs[0] = mxCreateSparse(m, n, nzmax, 0);
    sr  = mxGetPr(plhs[0]);                     // pointer to sparse values
    irs = mxGetIr(plhs[0]);                     // pointer to sparse row indices
    jcs = mxGetJc(plhs[0]);                     // pointer to sparse column indices
    cc = 0;                                     // current row position in sparse matrix
    ck = 0;                                     // current column position in sparse matrix

    /* Allocate space for matrix containing the distances between one datapoint and all other datapoints. */
    d     = malloc(m * sizeof(double));         // stores distances with n datapoints
    k_d   = malloc(k * sizeof(double));         // stores distance values to k nearest neighbors
    k_ind = malloc(k * sizeof(int));            // stores indices of k nearest neighbors
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            /* Compute Euclidean distance between datapoint i and j. */
            dist = 0;
            for (l = 0; l < n; l++) {
                tmp = pr[i * n + l] - pr[j * n + l];
                dist += (tmp * tmp);
            }
            d[j] = sqrt(dist);
        }
            
        /* Find the k nearest neighbors. */
/*        for (l = 0; l < k; l++) {
            min_ind = 0;
            min_d = 0xffffff;
            for (p = 0; p < m; p++) {
                if(d[p] < min_d) {
                    min_d   = d[p];
                    min_ind = p;
                }
            }
            k_d[l]     = min_d;
            k_ind[l]   = min_ind;
            d[min_ind] = 0xffffff;
        }     */
        
        /* Find the k nearest neighbors */
        rhs[0] = (mxArray*)d;
        mexCallMATLAB(2, lhs, 1, rhs, 'sort');
        d   = lhs[0];
        ind = lhs[1];
        for (l = 0; l < k; l++) {
            k_d[l]   = d[l];
            k_ind[l] = (int)ind[l];
        }
               
        /* Store the k nearest neighbors in the sparse dissimilarity matrix. */
        for (l = 0; l < k; l++) {
            if(cc < nzmax) {
                sr[cc]  = k_d[l];           // store sparse value
                irs[cc] = k_ind[l];         // store row index
                cc++;
            }
        }
        jcs[i] = ck;                        // store column index
        ck += k;        
    }
    jcs[m] = ck;

    /* Clean up memory. */
    free(d);
    free(k_d);
    free(k_ind);

    return;
}
      
      
