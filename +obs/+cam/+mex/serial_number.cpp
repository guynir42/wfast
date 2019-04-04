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
	
	cout << "Disconnected from camera!" << endl << endl;
	
  }

}
