#include "mex.h"
#include "matrix.h"
#include <stdio.h>

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){

	const size_t dims[4]={100,100,100,100};
  
	unsigned short int *out_array= (unsigned short int*) mxCalloc(dims[0]*dims[1]*dims[2]*dims[3], sizeof(unsigned short int)); 			  
	plhs[0]=mxCreateNumericArray(4, dims, mxUINT16_CLASS, mxREAL);
	mxSetData(plhs[0], out_array);
	
}