#include "mex.h"
#include "matrix.h"

#include <stdio.h>

#include "atcore.h"
// #include "atutility.h"

#include <iostream> 
using namespace std;

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){

	mexPrintf("trying to initialize...\n");
				  
	int i_retCode;
  cout << "Initialising ..." << endl << endl;
  i_retCode = AT_InitialiseLibrary();
  if (i_retCode != AT_SUCCESS) {
    mexPrintf("Error initialising library!\n");
  }
  else {
    AT_64 iNumberDevices = 0;
    AT_GetInt(AT_HANDLE_SYSTEM, L"Device Count", &iNumberDevices);
    if (iNumberDevices <= 0) {
      mexPrintf("No cameras detected! iNumberDevices= %d\n", iNumberDevices);
    }
    else {
      AT_H Hndl;
      i_retCode = AT_Open(0, &Hndl);
      if (i_retCode != AT_SUCCESS) {
        cout << "Error condition, could not initialise camera" << endl << endl;
      }
      else {
        cout << "Successfully initialised camera" << endl << endl;
      }
      AT_WC szValue[64];
      i_retCode= AT_GetString(Hndl, L"Serial Number", szValue, 64);
      if (i_retCode == AT_SUCCESS) {
        //The serial number of the camera is szValue
        mexPrintf("The serial number is %ld\n",szValue);
      }
      else {
        cout << "Error obtaining Serial number" << endl << endl;
      }
      AT_Close(Hndl);
    }
    AT_FinaliseLibrary();
  }

}

void func(){
	
	int rc=0; // return code! 
	
	long int *output = (long int *) mxCalloc(1, sizeof(long int));
	output[0] = AT_HANDLE_SYSTEM;
	
	rc = AT_InitialiseLibrary();
	if(rc) mexErrMsgIdAndTxt( "MATLAB:obs:cam:initialize:initialize_library", "Got error code %d when calling AT_InitializeLibrary!", rc);
	
	AT_64 num_devices=0;
	rc=AT_GetInt(AT_HANDLE_SYSTEM, L"DeviceCount", &num_devices);
	if(rc) mexErrMsgIdAndTxt( "MATLAB:obs:cam:initialize:device_count", "Got error code %d when getting DeviceCount!", rc);
	if(num_devices<=0) mexErrMsgIdAndTxt( "MATLAB:obs:cam:initialize:num_devices", "Found %d devices...", num_devices);
	
	
	
	
	mwSize dims[2] = {1,1};
	//plhs[0]=mxCreateNumericArray(1,dims, mxUINT64_CLASS, mxREAL);
	//mxSetData(plhs[0], output);
	
	}