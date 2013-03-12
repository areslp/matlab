#include "mex.h"
#include "math.h"
#include "string.h"


void computeKernelRow(double* X, int n, int d, int index, const char* function, int param1, int param2, double* row);
void computeColumnSums(double* X, int n, int d, const char* function, int param1, int param2, double* column_sums, double *total_sum);
void centerKernelRow(double* row, int index, int n, double* column_sums, double *total_sum);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // Initialize variables
    int i, j, n, d;
    double *v, *X, *column_sums, *total_sum, *row, *res;
    bool center;
    char *function, *type;
    double param1, param2;
    
    // Check and process inputs
    if(nrhs < 2) {
        mexErrMsgTxt("Not enough inputs.");
    }
    if(mxIsClass(prhs[0], "double")) {
        v = (double*) mxGetPr(prhs[0]);
    }
    else {
         mexErrMsgTxt("First input should be of type double.");
    }
    if(mxIsClass(prhs[1], "double")) {
        X = (double*) mxGetPr(prhs[1]);
    }
    else {
        mexErrMsgTxt("Second input should be of type double.");
    }
    if(nrhs < 3) {
        center = true;
    }
    else {
        if(mxIsClass(prhs[2], "logical")) {
            center = *((bool*) mxGetPr(prhs[2]));
        }
        else if(mxIsClass(prhs[2], "uint8")) {
            if (*((char*) mxGetPr(prhs[2])) == 0) center = false;
            else center = true;
        }
        else if(mxIsClass(prhs[2], "double")) {
            if (*((double*) mxGetPr(prhs[2])) == 0.0) center = false;
            else center = true;
        }
        else {
           mexErrMsgTxt("Third input should be of type logical, uint8, or double.");
        }
    }
    if(nrhs < 4) {
        function = "kernel";
    }
    else {
        if(mxIsClass(prhs[3], "char")) {
            function = (char*) mxGetPr(prhs[3]);
        }
        else {
            mexErrMsgTxt("Fourth input should be of type char.");
        }
    }
    if(nrhs < 5) {
        param1 = 1.0;
    }
    else {
        if(mxIsClass(prhs[4], "double")) {
            param1 = *((double*) mxGetPr(prhs[4]));
        }
        else {
            mexErrMsgTxt("Fifth input should be of type double.");
        }
    }
    if(nrhs < 6) {
        param2 = 3.0;
    }
    else {
        if(mxIsClass(prhs[5], "double")) {
            param2 = *((double*) mxGetPr(prhs[5]));
        }
        else {
            mexErrMsgTxt("Fifth input should be of type double.");
        }
    }
    if(nrhs < 7) {
        type = "Normal";
    }
    else {
        if(mxIsClass(prhs[6], "char")) {
            type = (char*) mxGetPr(prhs[6]);
        }
        else {
            mexErrMsgTxt("Sixth input should be of type char.");
        }
    }
    if(mxGetN(prhs[1]) != mxGetM(prhs[0]) && mxGetM(prhs[0]) != 0) {
        mexErrMsgTxt("Number of instances does not equal length of vector v.");
    }
        
    // Allocate some memory and set some variables
    n = mxGetN(prhs[1]);
    d = mxGetM(prhs[1]);
    plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL);
    res = mxGetPr(plhs[0]);
    
    // Compute column sums and total sum
    if(center || strcmp(type, "C") == 0) {
        column_sums = (double*) malloc(n * sizeof(double));
        total_sum   = (double*) malloc(sizeof(double));
        for(i = 0; i < n; i++) *(column_sums + i) = 0.0;
        *total_sum = 0.0;
        computeColumnSums(X, n, d, function, param1, param2, column_sums, total_sum);
        
        // Return only the column sums if type set to "ColumnSums"
        if(strcmp(type, "C") == 0) {
            for(i = 0; i < n; i++) {
                *(res + i) = *(column_sums + i);
            }
            free(column_sums);
            free(total_sum);
            return;
        }
    }
    
    // Compute all rows and multiply by v
    row = (double*) malloc(n * sizeof(double));
    for(i = 0; i < n; i++) {
        
        // Compute rows
        computeKernelRow(X, n, d, i, function, param1, param2, row);
        if(center) {
            centerKernelRow(row, i, n, column_sums, total_sum);
        }
        
        // Multiply row by v
        *(res + i) = 0.0;
        for(j = 0; j < n; j++) {
            *(res + i) += (*(row + j) * *(v + j));
        }
    }
    mexPrintf(".");
    
    // Clean up some memory
    free(row);
    if(center) {
        free(column_sums);
        free(total_sum);
    }
}


/**
 * 
 * Computes a single row of the kernel matrix, viz. row index.
 *
 */
void computeKernelRow(double* X, int n, int d, int index, const char* function, int param1, int param2, double* row) {

    // Initialize variables
    int i, j;
    
    // Compute linear kernel
    if(strcmp(function, "l") == 0) {
        for(i = 0; i < n; i++) {
            *(row + i) = 0.0;
            for(j = 0; j < d; j++) {
                *(row + i) += *(X + (index * d) + j) * *(X + (i * d) + j);
            }
        }        
    }
    
    // Compute polynomial kernel
    else if(strcmp(function, "p") == 0) {
        for(i = 0; i < n; i++) {
            *(row + i) = 0.0;
            for(j = 0; j < d; j++) {
                *(row + i) += *(X + (index * d) + j) * *(X + (i * d) + j);
            }
            *(row + i) = pow(*(row + i) + param1, param2);
        } 
    }
    
    // Compute Gaussian kernel
    else if(strcmp(function, "g") == 0) {
        for(i = 0; i < n; i++) {
            *(row + i) = 0.0;
            for(j = 0; j < d; j++) {
                *(row + i) += pow(*(X + (index * d) + j) - *(X + (i * d) + j), 2);
            }
            *(row + i) = exp(-(*(row + i) / (2 * pow(param1, 2))));
        } 
    }
    
    // Unknown kernel function
    else {
        mexErrMsgTxt("Unknown kernel function.");
    }
}


/**
 *
 * Compute all column sums and the total sum of a kernel matrix.
 *
 */
void computeColumnSums(double* X, int n, int d, const char* function, int param1, int param2, double* column_sums, double* total_sum) {
 
    // Initialize variables
    int i, j;
    double *row;
    
    // Preallocate memory
    row = (double*) malloc(n * sizeof(double));
    for(i = 0; i < n; i++) *(row + i) = 0.0;
    *total_sum = 0.0;
    
    // Compute column sums and total sum
    for(i = 0; i < n; i++) {
        
        // Compute kernel row
        computeKernelRow(X, n, d, i, function, param1, param2, row);
        
        // Update sums
        for(j = 0; j < n; j++) {
            *(column_sums + j) += *(row + j);
            *total_sum += *(row + j);
        }
    }
    *total_sum = *total_sum / pow(n, 2);
    
    // Clean up memory
    free(row);
}


/**
 *
 * Centers a kernel row as computed by computeKernelRow() using precomputed sums (e.g., sums as given by computeColumnSums(...)).
 *
 */
void centerKernelRow(double* row, int index, int n, double* column_sums, double* total_sum) {
    
    // Initialize variables
    int i;
    
    // Update kernel row
    for(i = 0; i < n; i++) {
        *(row + i) = *(row + i) - (*(column_sums + i) / (double) n) - (*(column_sums + index) / (double) n) + *total_sum;
    }
}
