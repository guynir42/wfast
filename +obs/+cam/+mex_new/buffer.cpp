#include "mex.h"
#include "matrix.h"
#include <math.h>
#include "atcore.h"
#include "atutility.h"

#define STRLN 256

// usage: buffer(hndl, command_string, parameter)

// this is used to report errors back to matlab (definition at the end)
void throw_error(const char *description, const char *par_name=0, int error_code=0); 

// utility to get enumerated value from camera
int get_enum(AT_H hndl, const AT_WC *par_name, AT_WC *str);

// utility functions to compare strings
bool cs(const char *keyword, const char *compare_str, int num_letters=3);
bool cs(const char *keyword, const char *str1, const char *str2, int num_letters=3);
bool cs(const char *keyword, const char *str1, const char *str2, const char *str3, int num_letters=3);

class BufferQueue {

	public:

	int index=0;
	int num_arrays=0;
	unsigned char **arrays=0;
	AT_64 im_size_bytes=0;
	unsigned char *latest_image=0;
	AT_H hndl;
	
	BufferQueue();
	~BufferQueue();
	
	void printout();
	
	void allocate(AT_H camera_handle, int num_images=10);
	void release(); 
	void queue(); 
	void wait(long int timeout=10000); 
	mxArray *BufferQueue::getImageMatrix();
	
} buffer; // create a global member in memory, which is deleted when "clear mex" is called or when compiling

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){

	// if no arguments are given, print the help section
	if (nrhs==0){
		const char *string[1]={"obs.cam.mex_new.buffer"};
		mxArray *array[1]={mxCreateCharMatrixFromStrings(1, string)};
		mexCallMATLAB(0,0,1,array,"help"); 
		return;
	}
	
	// check 1st argument is the camera handle
	if(mxIsScalar(prhs[0])==0 || mxIsNumeric(prhs[0])==0) throw_error("First input to buffer() must be numeric scalar (camera handle)"); 
	AT_H hndl = (AT_H) mxGetScalar(prhs[0]);
	
	// get the 2nd argument as the command to the buffers
	char key[STRLN];
	if(mxIsChar(prhs[1])==0) throw_error("Second input to buffer() must be character array (name of command)"); 
	mxGetString(prhs[1], key, STRLN); // copy the string data
	
	// optional 3rd argument is a scalar parameter
	double par=NAN;
	if(nrhs>2 && mxIsEmpty(prhs[2])==0) par=mxGetScalar(prhs[2]); 
	
	if(cs(key, "printout")){
		buffer.printout(); 
	}
	if(cs(key, "allocate")){
		
		if(isnan(par)) buffer.allocate(hndl); // use the default number of images
		else buffer.allocate(hndl,(int) par); // use the input number of images
		
	}
	else if(cs(key, "release")){
		buffer.release(); 
	}
	else if(cs(key, "queue")){
		buffer.queue(); 
	}
	else if(cs(key, "wait")){
		
		if(isnan(par)) buffer.wait();
		else buffer.wait((long int) par); 
		
		plhs[0]=buffer.getImageMatrix(); 
		
	}
	
}

BufferQueue::BufferQueue(){

}

BufferQueue::~BufferQueue(){

	release();

}

void BufferQueue::printout(){

	printf("hndl= %ld\n", (long int) hndl); 

}

void BufferQueue::allocate(AT_H camera_handle, int num_images){

	int rc=0; // return code

	hndl=camera_handle;

	release(); // first off release all the existing buffers! 
	
	rc=AT_InitialiseUtilityLibrary();
	if(rc) throw_error("Cannot initialize utilities library!", "allocate", rc); 
	
	num_arrays=num_images; // make sure to keep track of how many images were allocated! 
	
	rc=AT_GetInt(hndl, L"ImageSizeBytes", &im_size_bytes); // make sure to keep track of the size of each array
	if(rc) throw_error("Cannot get the image size from camera handle", "allocate", rc); 
	
	AT_BOOL bool_value=0;
	rc=AT_GetBool(hndl, L"CameraAcquiring", &bool_value); 
	if(rc) throw_error("Problem when checking if camera is running!", "allocate", rc); 

	if(bool_value){
		rc=AT_Command(hndl, L"AcquisitionStop"); 
		if(rc) throw_error("Problem when stopping camera!", "allocate", rc); 
	}
	
	rc=AT_Flush(hndl); 
	if(rc) throw_error("Problem when flushing buffers!", "allocate", rc); 
	
	arrays=new unsigned char*[num_images]; 

	index=0;
	for(int i=0;i<num_images;i++){
	
		arrays[i]=new unsigned char[im_size_bytes](); 
		rc=AT_QueueBuffer(hndl, arrays[i], im_size_bytes); 
		if(rc) throw_error("Problem when allocating buffers!", "allocate", rc); 
		index++; 
	}
	
	latest_image=0;
	
}

void BufferQueue::release(){

	int rc=0; // return code
		
	AT_BOOL bool_value=0;
	rc=AT_GetBool(hndl, L"CameraAcquiring", &bool_value); 
	if(rc) throw_error("Problem when checking if camera is running!", "release", rc); 

	if(bool_value){
		rc=AT_Command(hndl, L"AcquisitionStop"); 
		if(rc) throw_error("Problem when stopping camera!", "release", rc); 
	}
	
	rc=AT_FinaliseUtilityLibrary();
	if(rc) throw_error("Cannot release utilities library!", "allocate", rc); 
	
	rc=AT_Flush(hndl); 
	if(rc) throw_error("Problem when flushing buffers!", "release", rc); 
	
	if(arrays){
		for(int i=0;i<num_arrays;i++){
			
			delete[](arrays[i]); 
			arrays[i]=0; 
			
		}
		delete[](arrays); 
		arrays=0;
	}
	
	latest_image=0;

}

void BufferQueue::queue(){

	int rc=0; // return code
	
	latest_image=0;
	
	index++;
	if(index>=num_arrays) index=0;
	rc=AT_QueueBuffer(hndl, arrays[index], im_size_bytes); 
	if(rc) throw_error("Problem when queuing buffers!", "queue", rc); 
	
}

void BufferQueue::wait(long int timeout){

	int rc=0; // return code
	
	int buf_size=0;
	rc=AT_WaitBuffer(hndl, &latest_image, &buf_size, timeout); 
	if(rc) throw_error("Problem when waiting for buffers!", "wait", rc); 

	// printf("buf_size= %d\n", buf_size); 
	
	// do we want to verify that the same image we just got does not get queued again?

}

mxArray *BufferQueue::getImageMatrix(){
	
	int rc=0;
	
	AT_64 height=0;
	rc=AT_GetInt(hndl, L"AOIHeight", &height); 
	if(rc) throw_error("Cannot get the AOI height from camera!", "getImageMatrix", rc); 
	
	AT_64 width=0;
	rc=AT_GetInt(hndl, L"AOIWidth", &width); 
	if(rc) throw_error("Cannot get the AOI width from camera!", "getImageMatrix", rc); 
	
	AT_64 stride=0;
	rc=AT_GetInt(hndl, L"AOIStride", &stride); 
	if(rc) throw_error("Cannot get the AOI stride from camera!", "getImageMatrix", rc); 
	// printf("height= %d | width= %d | stride= %d | size= %d\n", height, width, stride, im_size_bytes); 
	
	// this will break if we are not using Mono16 encoding!
	AT_WC wide_str[STRLN]={0};
	rc=get_enum(hndl, L"PixelEncoding", wide_str);
	if(rc) throw_error("Cannot get the Pixel Encoding from camera!", "getImageMatrix", rc); 
	if(wcscmp(wide_str, L"Mono16")) throw_error("The camera is not in Mono16 encoding. Change the encoding or add options to convert the output buffers!"); 
	
	// height=height/2;
	// mxArray *matrix=mxCreateNumericMatrix(width, height, mxUINT16_CLASS, mxREAL); 
	// unsigned short *ptr=(unsigned short*) mxGetData(matrix);
	
	// mxArray *matrix=mxCreateNumericMatrix(im_size_bytes, 1, mxUINT8_CLASS, mxREAL); 
	// unsigned char *ptr=(unsigned char*) mxGetData(matrix);

	mxArray *matrix=mxCreateNumericMatrix(width, height, mxUINT16_CLASS, mxREAL); 
//	unsigned short *ptr=(unsigned short*) mxGetData(matrix);
	unsigned char *ptr=(unsigned char*) mxGetData(matrix);
	
	if(latest_image==0) throw_error("Array 'latest_image' has not been filled!", "getImageMatrix"); 
	
	// rc=AT_ConvertBuffer(latest_image, (AT_U8*) ptr, width, height, stride, L"Mono16", L"Mono16"); 
	
	// printf("sizeof(unsigned short)= %d | pixel values: %d %d %d %d %d\n", sizeof(unsigned short), latest_image[0], latest_image[1], latest_image[2], latest_image[3], latest_image[4]); 
	
	// printf("last pixel values in latest_image are: %d %d\n", latest_image[width*height*2-2], latest_image[width*height*2-1]); 
	
	// printf("last pixel values in matrix/ptr are: %d %d\n", latest_image[width*height*2-2], latest_image[width*height*2-1]); 
	
	// for(int i=0;i<height;i++) for(int j=0;j<width*2;j++) ptr[i*width*2+j]=latest_image[i*width*2+j]; 
	
	// memcpy(ptr, latest_image, width*height*2); 
	
	for(int i=0;i<height;i++){ // copy each row
		memcpy(ptr+(i*width*2), latest_image+i*stride, width*2);
	}
	
	return matrix; 

}

void throw_error(const char *description, const char *par_name, int error_code){

	char report[STRLN]={0};
	
	snprintf(report, STRLN, "%s", description); 
	if(par_name) snprintf(report, STRLN, "%s, parameter: %s", report, par_name);
	if(error_code) snprintf(report, STRLN, "%s, error code: %d", report, error_code); 
	
	mexErrMsgIdAndTxt("MATLAB:obs:cam:mex_new:buffer", report);

}

int get_enum(AT_H hndl, const AT_WC *par_name, AT_WC *str){

	int index=0;
	int rc=AT_GetEnumIndex(hndl, par_name, &index); 
	if(rc) throw_error("Cannot get enum index");  // need to thing of some safer way to transmit this error

	rc=AT_GetEnumStringByIndex(hndl, par_name, index, str, STRLN); 

	return rc;

}

bool cs(const char *keyword, const char *compare_str, int num_letters){ // compare two strings ignoring case and so on (for varargin parsing)
	
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
	
	return  success!=0;
	
}

bool cs(const char *keyword, const char *str1, const char *str2, int num_letters){ // compare two strings ignoring case and so on (for varargin parsing)

	return cs(keyword, str1, num_letters) || cs(keyword, str2, num_letters);

}

bool cs(const char *keyword, const char *str1, const char *str2, const char *str3, int num_letters){ // compare two strings ignoring case and so on (for varargin parsing)

	return cs(keyword, str1, num_letters) || cs(keyword, str2, num_letters) || cs(keyword, str3, num_letters);

}
