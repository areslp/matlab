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
	if(nrhs < 1){ 
		mexErrMsgTxt("At least one inputs required.");
	}
	if(nrhs > 2){
		mexErrMsgTxt("At most two inputs required.");
	}
	if(nlhs != 5){
		mexErrMsgTxt("Five output required.");
	}

	/* check to make sure the first input argument is a string */
	if( mxGetClassID(prhs[0]) != mxCHAR_CLASS){
		mexErrMsgTxt("Input x must be a string.");
	}
	char* filename;
	filename = mxArrayToString(prhs[0]);
	//mexPrintf("%s!!!!", filename);

	unsigned int htype = 0, dtype = 0;
	double hs = 2, rho = 3;
	
	if(nrhs == 2){
		const mxArray *opt = prhs[1];
		mxClassID category = mxGetClassID(opt);
		if(category != mxSTRUCT_CLASS){
			mexErrMsgTxt("Third input must be a structure");
		}
	 	mwSize total_num_of_elements;
  		mwIndex index;
  		int number_of_fields, field_index;
  		const char  *field_name;
  		const mxArray *field_array_ptr;
  		total_num_of_elements = mxGetNumberOfElements(opt); 
  		number_of_fields = mxGetNumberOfFields(opt);
  
  		/* Walk through each structure element. */
  		for (index=0; index<total_num_of_elements; index++)  {
    
    		/* For the given index, walk through each field. */ 
    		for (field_index=0; field_index<number_of_fields; field_index++)  {
         	field_name = mxGetFieldNameByNumber(opt, field_index);
				field_array_ptr = mxGetFieldByNumber(opt, index, field_index);
		      if (field_array_ptr == NULL) {
					continue;
				}

				if(strcmp(field_name, "hs") == 0){
					mxClassID   cat = mxGetClassID(field_array_ptr);
					if(cat == mxDOUBLE_CLASS){
						//mexPrintf("hs: %f", hs);
						hs = (*mxGetPr(field_array_ptr));
					}
				}
				else if(strcmp(field_name, "rho") == 0){
					mxClassID   cat = mxGetClassID(field_array_ptr);
					if(cat == mxDOUBLE_CLASS){
						//mexPrintf("rho: %f", rho);
						rho = (*mxGetPr(field_array_ptr));
					}
				}
				else if(strcmp(field_name, "htype") == 0){
					mxClassID   cat = mxGetClassID(field_array_ptr);
					if(cat == mxCHAR_CLASS){
						char buf[256];
    					mwSize buflen;

		    			/* Allocate enough memory to hold the converted string. */
				    	buflen = mxGetNumberOfElements(field_array_ptr) + 1;
						if(buflen > 256){
							buflen = 256;
						}

					  	/* Copy the string data from string_array_ptr and place it into buf. */
				  	  	if (mxGetString(field_array_ptr, buf, buflen) == 0){
					  		if( strcmp(buf, "ddr") == 0 ){
								htype = 0;	
							}
							else if( strcmp(buf, "psp") == 0 ){
								htype = 1;
							}
					  	}
					}
				}
				else if(strcmp(field_name, "dtype") == 0){
					mxClassID   cat = mxGetClassID(field_array_ptr);
					if(cat == mxCHAR_CLASS){
						char buf[256];
    					mwSize buflen;

		    			/* Allocate enough memory to hold the converted string. */
				    	buflen = mxGetNumberOfElements(field_array_ptr) + 1;
						if(buflen > 256){
							buflen = 256;
						}

					  	/* Copy the string data from string_array_ptr and place it into buf. */
				  	  	if (mxGetString(field_array_ptr, buf, buflen) == 0){
					  		if( strcmp(buf, "euclidean") == 0 ){
								dtype = 0;	
							}
							else if( strcmp(buf, "geodesic") == 0 ){
								dtype = 1;
							}
							//else if( strcmp(buf, "fastmarching") == 0 ){
							//	dtype = 2;
							//}

					  	}
					}
				}
			}// for field_index
      }//for index
   }//if(nrhs == 3)
	mexPrintf("dtype: %d, htype: %d, hs: %.2f, rho: %.2f\n", dtype, htype, hs, rho);

	//cout<<"nn: "<<nn<<" hs: "<<hs<<" rho: "<<rho<<endl;
	
	double *II, *JJ, *SS, *AA, *pt_h, h;
	vector<double> SSV, AAV;
	vector<unsigned int> IIV, JJV;
	
	if(dtype == 0){
		generate_sym_meshlp_matrix_matlab(filename, htype, hs, rho, h, IIV, JJV, SSV, AAV);
	}
	else{ // if (dtype == 1 )
		generate_sym_meshlp_matrix_geod_matlab(filename, htype, hs, rho, h, IIV, JJV, SSV, AAV);
	}
	//else{
	//	generate_sym_meshlp_matrix_geod_matlab_fastmarching(filename, htype, hs, rho, IIV, JJV, SSV, AAV);
	//}
	

	unsigned int dim = AAV.size();
	unsigned int nelem = IIV.size();
	if( nelem != JJV.size() ||  nelem != SSV.size() ){
		mexErrMsgTxt("Dimension of II, JJ, SS has to be the same");
	}
	plhs[0] = mxCreateDoubleMatrix(nelem, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(nelem, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(nelem, 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(dim, 1, mxREAL);
	plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);

	
	II = mxGetPr(plhs[0]);
	JJ = mxGetPr(plhs[1]);
	SS = mxGetPr(plhs[2]);
	AA = mxGetPr(plhs[3]);
	pt_h = mxGetPr(plhs[4]);

	cout<<"h: "<<h<<endl;
		
	*pt_h = h;

	for(mwSize i = 0; i < nelem; i ++){
		II[i] = IIV[i];
		JJ[i] = JJV[i];
		SS[i] = SSV[i];
	}
	for(mwSize i = 0; i < dim; i ++){
		AA[i] = AAV[i];
	}

	mxFree(filename);
}


