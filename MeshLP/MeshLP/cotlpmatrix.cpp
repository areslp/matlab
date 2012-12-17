#include <mex.h>
#include <math.h>
#include <vector>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <assert.h>

#include "meshlpmatrix.h"


using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
	/*  check for proper number of arguments */
	/* NOTE: You do not need an else statement when using mexErrMsgTxt
	 *       within an if statement, because it will never get to the else
	 *       statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
	 *       the MEX-file) 
	*/
	if(nrhs != 1){ 
		mexErrMsgTxt("One inputs required.");
	}
	if(nlhs != 4){
		mexErrMsgTxt("Four output required.");
	}

	/* check to make sure the first input argument is a string */
	if( mxGetClassID(prhs[0]) != mxCHAR_CLASS){
		mexErrMsgTxt("Input x must be a string.");
	}
	char* filename;
	filename = mxArrayToString(prhs[0]);
	//mexPrintf("%s!!!!", filename);

	//cout<<"nn: "<<nn<<" hs: "<<hs<<" rho: "<<rho<<endl;
	
	double *II, *JJ, *SS, *DD;
	vector<double> SSV, DDV;
	vector<unsigned int> IIV, JJV;
	
	generate_Xu_Meyer_laplace_matrix_matlab(filename, IIV, JJV, SSV, DDV);

	unsigned int dim = DDV.size();
	unsigned int nelem = IIV.size();
	if( nelem != JJV.size() ||  nelem != SSV.size() ){
		mexErrMsgTxt("Dimension of II, JJ, SS has to be the same");
	}
	plhs[0] = mxCreateDoubleMatrix(nelem, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(nelem, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(nelem, 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(dim, 1, mxREAL);
	II = mxGetPr(plhs[0]);
	JJ = mxGetPr(plhs[1]);
	SS = mxGetPr(plhs[2]);
	DD = mxGetPr(plhs[3]);

	for(mwSize i = 0; i < nelem; i ++){
		II[i] = IIV[i];
		JJ[i] = JJV[i];
		SS[i] = SSV[i];
	}
	for(mwSize i = 0; i < dim; i ++){
		DD[i] = DDV[i];
	}

	mxFree(filename);
}


