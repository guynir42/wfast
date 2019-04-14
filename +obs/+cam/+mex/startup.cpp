#include "mex.h"
#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <thread>
#include <chrono>
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "CameraControl.h"
#include "SimCameraControl.h"
#include "ZylaCameraControl.h"

#define STRLN 256

int cs(const char *keyword, const char *compare_str, int num_letters=3);
int cs(const char *keyword, const char *str1, const char *str2, int num_letters=3);

#define INDEX_CAM 0
#define INDEX_FLAG 1
#define INDEX_BUF 2
#define INDEX_REC 3
#define INDEX_NUM 4
#define INDEX_SIZE 5
#define INDEX_LOG 6

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){
	
	if(nrhs<1){ // usage is detailed when function is called without arguments!
		
		mexPrintf("Usage: capture(cam_object, mex_flag, buffer_struct_array, index_rec_vector, num_batches=1, [batch_size], error_log=[])\n");
		mexPrintf("Activates camera and starts to capture images into the buffer structure.\n");
		mexPrintf("Will stop capturing when 'num_batches' is reached or when 'mex_flag' second element is set to 1\n");
		mexPrintf("\nInputs: \n");
		mexPrintf("-cam_object is a CameraControl object.\n");		
		mexPrintf("-mex_flag should be input as [0,0,0,0]. \n   mex_flag[0] means camera is now recording, \n   mex_flag[1] is used to stop recording. \n   mex_flag[2] is an error flag, \n   mex_flag[3] is counter\n");
		mexPrintf("-buffer_struct_array is the struct array inside the BufferWheel.\n");
		mexPrintf("-index_rec_vector is a vector indicating which struct to write to.\n");
		mexPrintf("-num_batches specifies how many batches to capture / buffers to fill before stopping (default=1).\n");
		mexPrintf("-batch_size is an optional argument that overrides the batch_size in cam_object (number of frames per batch).\n");
		mexPrintf("-error log is used to record non-critical errors from camera.\n");

		return;
	}
	
	////////// check minimal inputs ///////////
	if(nrhs<INDEX_NUM) mexErrMsgIdAndTxt("MATLAB:obs:mexWrite:invalidNumInputs", "At least %d arguments are required!", INDEX_NUM);
	if(mxIsClass(prhs[INDEX_CAM], "obs.cam.CameraControl")!=1)
		mexErrMsgIdAndTxt( "MATLAB:obs:cam:capture:inputWrongClass", "Input %d must be a CameraControl object.", INDEX_CAM+1);
	if(mxGetN(prhs[INDEX_FLAG])!=4 || mxGetM(prhs[INDEX_FLAG])!=1 || mxIsNumeric(prhs[INDEX_FLAG])==0) 
		mexErrMsgIdAndTxt( "MATLAB:obs:cam:capture:inputSizeMismatch", "Input %d must be a numeric 4 element vector.", INDEX_FLAG+1);
	if(mxIsStruct(prhs[INDEX_BUF])!=1)
		mexErrMsgIdAndTxt( "MATLAB:obs:cam:capture:inputNotStruct", "Input %d must be a BufferWheel 'buf' structure array.", INDEX_BUF+1);
	if(mxGetN(prhs[INDEX_REC])!=2 || mxGetM(prhs[INDEX_REC])!=1 || mxIsNumeric(prhs[INDEX_REC])==0) 
		mexErrMsgIdAndTxt( "MATLAB:obs:cam:capture:inputWrongSize", "Input %d must be a numeric 2 element vector.", INDEX_REC+1); 
	if(nrhs>=INDEX_NUM && (mxIsEmpty(prhs[INDEX_NUM]) || mxIsScalar(prhs[INDEX_NUM])!=1 || mxIsNumeric(prhs[INDEX_NUM])!=1) ) 
		mexErrMsgIdAndTxt( "MATLAB:obs:cam:capture:inputNotScalar", "Input %d must be a numeric scalar.", INDEX_NUM+1); 
	if(nrhs>=INDEX_SIZE && mxIsEmpty(prhs[INDEX_SIZE])!=1 && (mxIsScalar(prhs[INDEX_SIZE])!=1 || mxIsNumeric(prhs[INDEX_SIZE])!=1) )
		mexErrMsgIdAndTxt( "MATLAB:obs:cam:capture:inputNotScalar", "Input %d must be a numeric scalar.", INDEX_SIZE+1); 
	
	const mxArray *camera=prhs[INDEX_CAM]; // MATLAB object
	mxArray *buffers=(mxArray*)prhs[INDEX_BUF]; // struct array	
	mxArray *prop=(mxArray *)mxGetProperty(camera, 0, "cam_name");
	if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:capture:cannotFindProperty", "Cannot find property 'cam_name' inside CameraControl!");
	const char *cam_name=mxArrayToString(prop);
	
	CameraControl *cc=0;
	//	choose which camera to use
	if (cs("sim", cam_name)){ 
		cc=new SimCameraControl;
	}
	else if (cs("zyla", cam_name)){ 
		cc=new ZylaCameraControl;
	}
	else mexErrMsgIdAndTxt( "MATLAB:obs:cam:capture:unknownCameraName", "Unknown camera name '%s', use 'sim' or 'zyla'...", cam_name);

	cc->mex_flag_cam=mxGetPr(prhs[INDEX_FLAG]);
	
	if(cc->mex_flag_cam[0]){
		if(cc->debug_bit) printf("Camera is already recording... \n"); 
		return; 
	}
	cc->num_batches=1;	
	if(nrhs>INDEX_NUM) cc->num_batches=(unsigned int) mxGetScalar(prhs[INDEX_NUM]);
	
	
	cc->loadFromCamera(camera);
	cc->index_rec_vector=mxGetPr(prhs[INDEX_REC]);
	
	if(nrhs>INDEX_SIZE && mxIsEmpty(prhs[INDEX_SIZE])!=1) 
		cc->batch_size=mxGetScalar(prhs[INDEX_SIZE]);
	
	if(nrhs>INDEX_LOG && mxIsEmpty(prhs[INDEX_LOG])!=1) 
		cc->error_log=mxGetPr(prhs[INDEX_LOG]);
	
	cc->loadFromBuffers(buffers);
	
	if(cc->debug_bit>1) cc->printout();
	
	///////// reading buffer wheel strut //////////
	
	std::thread mythread(&CameraControl::loop, cc);
	mythread.detach();
	
}

int cs(const char *keyword, const char *compare_str, int num_letters){
	
	char str1[STRLN]={0};
	char str2[STRLN]={0};
	
	// clean up string 1 (keyword)
	size_t N=STRLN;
	if(strlen(keyword)<N) N=strlen(keyword);
	
	int j=0;
	for(int i=0; i<N; i++){
		
		if(keyword[i]=='_' || keyword[i]==' ') continue;
		
		str1[j]=tolower(keyword[i]);
		j++;
		
	}
	
	// clean up string 2 (compare_str)
	N = STRLN;
	if(strlen(compare_str)<N) N=strlen(compare_str);
	
	j=0;
	for(int i=0; i<N; i++){
		
		if(compare_str[i]=='_' || compare_str[i]==' ') continue;
		
		str2[j]=tolower(compare_str[i]);
		j++;
		
	}
	
	// compare the strings
	int success=1;
	
	N=num_letters;
	if(strlen(str1)>N) N=strlen(str1); // number of letters to compare (minimum 3, or length of keyword).
	
	for(int i=0;i<N;i++){
		
		if(str1[i]!=str2[i]){
			success=0;
			break;
		}
		
	}
	
	return success;
	
}

int cs(const char *keyword, const char *str1, const char *str2, int num_letters){

	return cs(keyword, str1, num_letters) || cs(keyword, str2, num_letters);

}
