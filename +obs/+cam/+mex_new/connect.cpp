#include "mex.h"
#include "matrix.h"

#include <stdio.h>

#include "atcore.h"
// #include "atutility.h"

#include <iostream> 

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){

	printf("Connecting to Andor camera!\n"); 
				  
	int rc=0; // return code! 
	
	long int *output = (long int *) mxCalloc(1, sizeof(long int));
	output[0] = AT_HANDLE_SYSTEM;
	
	rc = AT_InitialiseLibrary();
	if(rc) mexErrMsgIdAndTxt( "MATLAB:obs:cam:mex:connect:initializeLibrary", "Got error code %d when calling AT_InitialiseLibrary!", rc);
	
	AT_64 num_devices=0;
	rc=AT_GetInt(AT_HANDLE_SYSTEM, L"DeviceCount", &num_devices);
	if(rc) mexErrMsgIdAndTxt( "MATLAB:obs:cam:mex:connect:deviceCount", "Got error code %d when getting DeviceCount!", rc);
	if(num_devices<=0) mexErrMsgIdAndTxt( "MATLAB:obs:cam:mex:connect:numDevices", "Found %d devices...", num_devices);
	
	// printf("sizeof(AT_H)= %d | sizeof(long int)= %d\n", sizeof(AT_H), sizeof(long int));
	
	AT_H *hndl=(AT_H*)mxCalloc(1,sizeof(AT_H));
	
	rc = AT_Open(0, hndl);
	if(rc) mexErrMsgIdAndTxt( "MATLAB:obs:cam:mex:connect:open", "Got error code %d when using AT_Open!", rc);
	
	printf("hndl= %ld\n", *hndl);
	
	rc=AT_Command(*hndl, L"AcquisitionStop");
	if(rc) mexErrMsgIdAndTxt( "MATLAB:obs:cam:mex:connect:stopAcquistionFailed", "Got error code %d when trying to stop acquisition!", rc);

	rc=AT_SetBool(*hndl, L"SensorCooling", (bool) 1); 
	if(rc) mexErrMsgIdAndTxt( "MATLAB:obs:cam:mex:connect:setSensorCooling", "Got error code %d when setting sensor cooling!", rc);

	rc=AT_SetEnumString(*hndl, L"FanSpeed", L"On");
	if(rc) mexErrMsgIdAndTxt( "MATLAB:obs:cam:mex:connect:setFanSpeed", "Got error code %d when setting fan speed!", rc);

	AT_BOOL bool_value=0;
	rc=AT_IsWritable(*hndl, L"TargetSensorTemperature", &bool_value);
	if(rc) mexErrMsgIdAndTxt( "MATLAB:obs:cam:mex:connect:checkTempIsWritable", "Got error code %d when checking if target temperature is writable!", rc);

	if(bool_value){ 
		rc=AT_SetFloat(*hndl, L"TargetSensorTemperature", 10); // set target temperature to 10. If it was set too low it causes a camera overheat! 
		if(rc) mexErrMsgIdAndTxt( "MATLAB:obs:cam:mex:connect:setTargetTemperature", "Got error code %d when setting target temperature!", rc);
	}
	
	mwSize dims[2] = {1,1};
	if(nlhs>0){
		plhs[0]=mxCreateNumericArray(1,dims, mxUINT64_CLASS, mxREAL);
		mxSetData(plhs[0], hndl);
	}
	else{ 
	
		printf("Successfully connected to camera with handle %d\n", *hndl);
		
		printf("Disconnecting from camera!\n");
		
		rc = AT_Close(*hndl); 
		if(rc) mexErrMsgIdAndTxt("MATLAB:obs:cam:disconnect:cannotClose", "Unable to close handl %ld. Received error code %d", *hndl, rc);
		
		rc = AT_FinaliseLibrary();
		if(rc) mexErrMsgIdAndTxt("MATLAB:obs:cam:disconnect:cannotClose", "Unable to finalize library. Received error code %d", rc);
		
	}
	
}

