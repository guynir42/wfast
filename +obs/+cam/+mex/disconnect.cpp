#include "mex.h"
#include "matrix.h"

#include <stdio.h>

#include "atcore.h"
// #include "atutility.h"

#include <iostream> 

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){

	if(nrhs<1) mexErrMsgIdAndTxt( "MATLAB:obs:cam:disconnect:noArguments", "Must give a handle to a camera in order to disconnect!");
	
	if(mxIsNumeric(prhs[0])==0) mexErrMsgIdAndTxt("MATLAB:obs:cam:disconnect:inputNotNumeric", "Input 1 to mexCutout is not numeric...");
	
	int rc=0; // return code! 
	
	AT_H hndl = mxGetScalar(prhs[0]);
	
	rc = AT_Close(hndl); 
	if(rc) mexErrMsgIdAndTxt("MATLAB:obs:cam:disconnect:cannotClose", "Unable to close handl %ld. Received error code %d", hndl, rc);
	
	rc = AT_FinaliseLibrary();
	if(rc) mexErrMsgIdAndTxt("MATLAB:obs:cam:disconnect:cannotClose", "Unable to finalize library. Received error code %d", rc);
	
}