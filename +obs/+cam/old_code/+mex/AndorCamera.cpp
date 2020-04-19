#include "AndorCamera.h"

AndorCamera::AndorCamera(){

	_err_code=0;
	
	printf("Initializing library!\n"); 
	
	_err_code = AT_InitialiseLibrary(); // add error check? 
	
	if(_err_code){ printf("Failed to initialize library (%d)!\n", _err_code); return; }
	
	AT_64 num_devices=0;
	_err_code=AT_GetInt(AT_HANDLE_SYSTEM, L"DeviceCount", &num_devices);
	
	if(_err_code || num_devices<=0){ printf("Cannot find any devices (%d)!\n", _err_code); return; }
	else printf("Found %d devices!\n", (int) num_devices); 
	
	// at this stage the library is initialized but no camera is yet connected!
	connect(); 
	
}

AndorCamera::~AndorCamera(){

	printf("Deleting camera object!\n"); 

	// get rid of all allocated resources! 
	
	// close camera handle
	disconnect(); 
	
	// finalize library
	_err_code = AT_FinaliseLibrary();
	if(_err_code){ printf("Problem finalizing library (%d)!\n", _err_code); return; }
	
}

void AndorCamera::connect(){
	
	_err_code = AT_Open(0, &_hndl);
	if(_err_code){ printf("Problem opening the camera handle (%d)!\n", _err_code); return; }
	else printf("_hndl= %ld\n", _hndl);
	
}

void AndorCamera::disconnect(){

	_err_code = AT_Close(_hndl); 

	if(_err_code){ printf("Problem closing camera handle (%d)!\n", _err_code); return; }
	else printf("Camera handle closed!\n"); 
	
	_hndl=0;
	
}

double AndorCamera::getImageSize(){

	AT_64 value=0;
	
	_err_code=AT_GetInt(_hndl, L"ImageSizeBytes", &value);
	printf("value= %d\n", value); 
	if(_err_code){ printf("Problem getting image size (%d)!\n", _err_code); return 0; }
	
	return (double) value; 

}