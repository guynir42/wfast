#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){
	
	
	// check inputs!
	if (nrhs==0){ 
		mexPrintf("Usage: S = sum_singles(matrix)\n");
		mexPrintf("Takes inputs of types double, single and uint16 and returns a sum along the 3rd dimension in single precision.\n");
		return;
	}
	
	if(mxIsEmpty(prhs[0])){ // no input, then just return with all empty outputs...
		const mwSize dims[]={0,0};
		plhs[0]=mxCreateNumericArray(0,dims, mxSINGLE_CLASS, mxREAL); // return an empty array
		return;
	}
	
	
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
	
	if(mxIsClass(prhs[0], "double")){
		
		double *matrix=(double*) mxGetData(prhs[0]);
		
		for(int k=0;k<N_cutouts;k++){
			
			for(int j=0;j<N_frames;j++){// go over all frames
			
				for(int i=0;i<N;i++){ //  go over pixels in each image
					
					output[k*N+i]+=matrix[k*N*N_frames+j*N+i];
					
				} // for i (pixels in image)
				
			} // for j (frames)
			
		}// for k 
		
	}
	else if(mxIsClass(prhs[0], "single")){
		
		float *matrix=(float*) mxGetData(prhs[0]);
		
		for(int k=0;k<N_cutouts;k++){
			
			for(int j=0;j<N_frames;j++){// go over all frames
			
				for(int i=0;i<N;i++){ //  go over pixels in each image
					
					output[k*N+i]+=matrix[k*N*N_frames+j*N+i];
					
				} // for i (pixels in image)
				
			} // for j (frames)
			
		}// for k 
		
	}
	else if(mxIsClass(prhs[0], "uint16")){
		
		short unsigned int *matrix=(short unsigned int*) mxGetData(prhs[0]);
		
		for(int k=0;k<N_cutouts;k++){
			
			for(int j=0;j<N_frames;j++){// go over all frames
			
				for(int i=0;i<N;i++){ //  go over pixels in each image
					
					output[k*N+i]+=matrix[k*N*N_frames+j*N+i];
					
				} // for i (pixels in image)
				
			} // for j (frames)
			
		}// for k 
		
	}
	else mexErrMsgIdAndTxt("MATLAB:util:stat:sum_singles:inputTypeUnrecognized", "Input 1 to sum_singles is not a double, single or uint16 array...");
	
}