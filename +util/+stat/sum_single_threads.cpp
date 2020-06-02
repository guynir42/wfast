#include "mex.h"

#include <stdio.h>
#include <string.h>
#include <cctype>
#include <thread>
#include <chrono>
#include <vector>


#define STRLN 64 // maximum string length (for copying)

void sum_section_double(float *output, unsigned char *matrix, int start_idx, int finish_idx, int N, int N_frames, int N_cutouts, int bytes); 
void sum_section_single(float *output, unsigned char *matrix, int start_idx, int finish_idx, int N, int N_frames, int N_cutouts, int bytes); 
void sum_section_uint16(float *output, unsigned char *matrix, int start_idx, int finish_idx, int N, int N_frames, int N_cutouts, int bytes); 

// utility functions to compare strings
bool cs(const char *keyword, const char *compare_str, int num_letters=3);
bool cs(const char *keyword, const char *str1, const char *str2, int num_letters=3);
bool cs(const char *keyword, const char *str1, const char *str2, const char *str3, int num_letters=3);
bool parse_bool(mxArray *value);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){
	
	
	int num_threads=1; // the default is not to use multithreading
	
	// check inputs!
	if (nrhs==0){
		
		const char *string[1]={"util.stat.sum_single"};
		mxArray *array[1]={mxCreateCharMatrixFromStrings(1, string)};
		mexCallMATLAB(0,0,1,array,"help"); 
		return;
		
	}
	
	if(mxIsEmpty(prhs[0])){ // empty input, then just return with an empty output
		const mwSize dims[]={0,0};
		plhs[0]=mxCreateNumericArray(0,dims, mxSINGLE_CLASS, mxREAL); // return an empty array
		return;
	}


	for(int i=1;i<nrhs;i+=2){ // parse varargin (skip input 0 which is the data cube)
	
		char key[STRLN];
		if(mxIsChar(prhs[i])==0) mexErrMsgIdAndTxt("MATLAB:util:stat:sum_single:inputNotChar", "Input %d to sum_single is not a string...", i+1);
		mxGetString(prhs[i], key, STRLN); // copy the string data
		
		mxArray *val=0;
		if(i+1<nrhs) val=(mxArray*) prhs[i+1]; // if the varargin is odd numbered, leave val=0 as default
		
		if(cs(key, "threads")){
			if(val==0 || mxIsEmpty(val)) mexErrMsgIdAndTxt("MATLAB:util:stat:sum_single:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			else{
				if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:stat:sum_single:inputNotNumericScalar", "Input %d to sum_single is not a numeric scalar...", i+2);
				num_threads=(int) mxGetScalar(val);
			}
			
		}
	
	}// for i (varargin)
	
	mwSize *dims=(mwSize*) mxGetDimensions(prhs[0]); // dims[0] is the height while dims[1] is the width while N_frames is the number of frames
	int ndims=mxGetNumberOfDimensions(prhs[0]);	
	
	if(ndims>4) mexErrMsgIdAndTxt("MATLAB:util:stat:sum_singles:inputMoreThan4D", "Input 1 to sum_singles should have at most 4 dimensions...");
	
	mwSize N_frames=1;
	mwSize N_cutouts=1;
	if(ndims>=3) N_frames=dims[2];
	if(ndims>=4) N_cutouts=dims[3];
	
	mwSize dims_out[4]={dims[0], dims[1],1,N_cutouts};
	int N=dims_out[0]*dims_out[1]; // size of each image
	
	plhs[0]=mxCreateNumericArray(4, dims_out, mxSINGLE_CLASS, mxREAL); // create the output matrix (in matlab)
	float *output=(float*) mxGetData(plhs[0]); // get the C type array 
	
	unsigned char *matrix=(unsigned char*) mxGetData(prhs[0]);
		
	int bytes=0; 
	void (*sum_section)(float *output, unsigned char *matrix, int start_idx, int finish_idx, int N, int N_frames, int N_cutouts, int bytes);
	
	int num_pixels=N/num_threads;
	int start_idx=0;
	int finish_idx=0;
	
	if(mxIsClass(prhs[0], "double")){
		bytes=8;
		sum_section=&sum_section_double;
	}
	else if(mxIsClass(prhs[0], "single")){
		bytes=4;
		sum_section=&sum_section_single;
	}
	else if(mxIsClass(prhs[0], "uint16")){
		bytes=2;
		sum_section=&sum_section_uint16;
	}
	
	else mexErrMsgIdAndTxt("MATLAB:util:stat:sum_singles:inputTypeUnrecognized", "Input 1 to sum_singles is not a double, single or uint16 array...");
	
	std::vector<std::thread> threads;
	
	for(int t=0;t<num_threads;t++){
		
		if(t<num_threads-1){
			start_idx=t*num_pixels;
			finish_idx=(t+1)*num_pixels;
			threads.push_back(std::thread(sum_section, output, matrix, start_idx, finish_idx, N, N_frames, N_cutouts, bytes)); 
			// sum_section(output, matrix, start_idx, finish_idx, N, N_frames, N_cutouts, bytes);
		}
		else{ // last iteration goes on the main thread				
			start_idx=t*num_pixels;
			finish_idx=N; // last iteration gets the additional round-off pixels
			sum_section(output, matrix, start_idx, finish_idx, N, N_frames, N_cutouts, bytes);
		}
		
	}
	
	for(int t=0;t<num_threads-1;t++){
		threads[t].join();
	}

}

void sum_section_double(float *output, unsigned char *matrix, int start_idx, int finish_idx, int N, int N_frames, int N_cutouts, int bytes){
	
	double new_value=0;
		
	for(int k=0;k<N_cutouts;k++){
		
		for(int j=0;j<N_frames;j++){// go over all frames
			
			// printf("j= %d | new_value= %f\n", j, new_value);
			
			for(int i=start_idx;i<finish_idx;i++){ //  go over pixels in each image
				
				memcpy(&new_value, &matrix[k*N*N_frames*bytes+j*N*bytes+i*bytes], bytes);
				if(new_value==new_value) output[k*N+i]+=new_value; // condition is to remove NaN values
				
			} // for i (pixels in image)
			
		} // for j (frames)
		
	}// for k 
}

void sum_section_single(float *output, unsigned char *matrix, int start_idx, int finish_idx, int N, int N_frames, int N_cutouts, int bytes){
	
	float new_value=0;
		
	for(int k=0;k<N_cutouts;k++){
		
		for(int j=0;j<N_frames;j++){// go over all frames
			
			// printf("j= %d | new_value= %f\n", j, new_value);
			
			for(int i=start_idx;i<finish_idx;i++){ //  go over pixels in each image
				
				memcpy(&new_value, &matrix[k*N*N_frames*bytes+j*N*bytes+i*bytes], bytes);
				if(new_value==new_value) output[k*N+i]+=new_value; // condition is to remove NaN values
				
			} // for i (pixels in image)
			
		} // for j (frames)
		
	}// for k 
}

void sum_section_uint16(float *output, unsigned char *matrix, int start_idx, int finish_idx, int N, int N_frames, int N_cutouts, int bytes){
	
	unsigned short int new_value=0;
		
	for(int k=0;k<N_cutouts;k++){
		
		for(int j=0;j<N_frames;j++){// go over all frames
			
			// printf("j= %d | new_value= %f\n", j, new_value);
			
			for(int i=start_idx;i<finish_idx;i++){ //  go over pixels in each image
				
				memcpy(&new_value, &matrix[k*N*N_frames*bytes+j*N*bytes+i*bytes], bytes);
				if(new_value==new_value) output[k*N+i]+=new_value; // condition is to remove NaN values
				
			} // for i (pixels in image)
			
		} // for j (frames)
		
	}// for k 
}

bool cs(const char *keyword, const char *compare_str, int num_letters){ // compare two strings ignoring case and so one (for varargin parsing)
	
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

bool cs(const char *keyword, const char *str1, const char *str2, int num_letters){ // compare two strings ignoring case and so one (for varargin parsing)

	return cs(keyword, str1, num_letters) || cs(keyword, str2, num_letters);

}

bool cs(const char *keyword, const char *str1, const char *str2, const char *str3, int num_letters){ // compare two strings ignoring case and so one (for varargin parsing)

	return cs(keyword, str1, num_letters) || cs(keyword, str2, num_letters) || cs(keyword, str3, num_letters);

}

bool parse_bool(mxArray *value){ // if input is string of yes/no or on/off or number 1/0 output a boolean corresponding to the value
	
	if (mxIsEmpty(value)) return 0;
	if (mxIsScalar(value)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotScalar", "Input 1 to parse_bool must be a scalar!");
	if (mxIsLogical(value)) return mxGetScalar(value)!=0;
	if (mxIsNumeric(value)) return mxGetScalar(value)!=0;
	if (mxIsChar(value)) return cs(mxArrayToString(value), "yes", "on");
	return 0;
}


/*
		double *matrix=(double*) mxGetData(prhs[0]);
		
		for(int k=0;k<N_cutouts;k++){
			
			for(int j=0;j<N_frames;j++){// go over all frames
			
				for(int i=0;i<N;i++){ //  go over pixels in each image
					
					output[k*N+i]+=matrix[k*N*N_frames+j*N+i];
					
				} // for i (pixels in image)
				
			} // for j (frames)
			
		}// for k 
		
*/	