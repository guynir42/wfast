#include "mex.h"
#include "matrix.h"

#include <stdio.h>

#include "atcore.h"
// #include "atutility.h"

#include <iostream> 

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){

	printf("Trying to connect to Andor camera!\n"); 
				  
	int rc=0; // return code! 
	
	long int *output = (long int *) mxCalloc(1, sizeof(long int));
	output[0] = AT_HANDLE_SYSTEM;
	
	rc = AT_InitialiseLibrary();
	if(rc) mexErrMsgIdAndTxt( "MATLAB:obs:cam:mex:connect:initializeLibrary", "Got error code %d when calling AT_InitialiseLibrary!", rc);
	
	AT_64 num_devices=0;
	rc=AT_GetInt(AT_HANDLE_SYSTEM, L"DeviceCount", &num_devices);
	if(rc) mexErrMsgIdAndTxt( "MATLAB:obs:cam:mex:connect:deviceCount", "Got error code %d when getting DeviceCount!", rc);
	if(num_devices<=0) mexErrMsgIdAndTxt( "MATLAB:obs:cam:mex:connect:numDevices", "Found %d devices...", num_devices);
	
	mexPrintf("sizeof(AT_H)= %d | sizeof(long int)= %d\n", sizeof(AT_H), sizeof(long int));
	
	AT_H *hndl=(AT_H*)mxCalloc(1,sizeof(AT_H));
	
	rc = AT_Open(0, hndl);
	if(rc) mexErrMsgIdAndTxt( "MATLAB:obs:cam:mex:connect:open", "Got error code %d when using AT_Open!", rc);
	
	mexPrintf("hndl= %ld\n", *hndl);
	
	rc=AT_SetBool(*hndl, L"SensorCooling", (bool) 1); 
	rc=AT_SetEnumString(*hndl, L"FanSpeed", L"On");
	
	mwSize dims[2] = {1,1};
	if(nlhs>0){
		plhs[0]=mxCreateNumericArray(1,dims, mxUINT64_CLASS, mxREAL);
		mxSetData(plhs[0], hndl);
	}
	else{ 
	
		mexPrintf("Successfully connected to camera with handle %d\n", *hndl);
		
		mexPrintf("Disconnecting from camera!\n");
		
		rc = AT_Close(*hndl); 
		if(rc) mexErrMsgIdAndTxt("MATLAB:obs:cam:disconnect:cannotClose", "Unable to close handl %ld. Received error code %d", *hndl, rc);
		
		rc = AT_FinaliseLibrary();
		if(rc) mexErrMsgIdAndTxt("MATLAB:obs:cam:disconnect:cannotClose", "Unable to finalize library. Received error code %d", rc);
		
	}
	
}

