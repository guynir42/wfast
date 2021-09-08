#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>

#include <thread>
#include <chrono>

void find_peaks(unsigned short **array, unsigned short *images, size_t rows, size_t cols, int start_idx, int end_idx, bool *mask, unsigned short threshold, int max_number);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){

	// check inputs!
	if (nrhs==0){ // no inputs, output the help
		
		const char *string[1]={"util.img.find_flares"};
		mxArray *array[1]={mxCreateCharMatrixFromStrings(1, string)};
		mexCallMATLAB(0,0,1,array,"help"); 
		return;
		
	}
				  
	// check argument 1 (input image)
	if(mxIsNumeric(prhs[0])==0 || mxIsClass(prhs[0], "uint16")==0) mexErrMsgIdAndTxt("util:img:find_flares:inputNotUINT16", "Input 1 to find_flares must be uint16 class!");
	if(mxIsEmpty(prhs[0])){ // empty input matrix just returns an empty matrix
		plhs[0]=mxDuplicateArray(prhs[0]);
		if(nlhs>1) plhs[1]=mxDuplicateArray(prhs[0]);
		return; 
	}
	if( mxGetN(prhs[0])<2 || mxGetM(prhs[0])<2 || mxGetNumberOfDimensions(prhs[0])>3 ) mexErrMsgIdAndTxt("util:img:find_flares:wrongInputDimensions", "Must input a (2D or 3D) matrix!");
	
	// get the dimensions of the input
	size_t rows=mxGetM(prhs[0]);
	size_t cols=mxGetN(prhs[0]);
	size_t pages = 1;
	if(mxGetNumberOfDimensions(prhs[0])>2){
		
		const mwSize *dims=mxGetDimensions(prhs[0]);
		rows=dims[0];
		cols=dims[1];
		pages=dims[2];
		
	}
	
	// check argument 2 (source/bad pixel mask)
	if(mxIsLogical(prhs[1])==0) mexErrMsgIdAndTxt("util:img:find_flares:inputNotLogical", "Input 2 to find_flares is not logical class!");
	if( mxGetN(prhs[1])!=cols || mxGetM(prhs[1])!=rows ) mexErrMsgIdAndTxt("util:img:find_flares:sizeMismatch", "The x and y size of input 1 and 2 must be the same!");
	
	// check argument 3 (threshold)
	int threshold=256;
	if(nrhs>2 && mxIsEmpty(prhs[2])==0){
		if(mxIsNumeric(prhs[2])==0 || mxIsScalar(prhs[2])==0) mexErrMsgIdAndTxt("util:img:find_flares:inputNotNumericScalar", "Input 3 to find_flares is not a numeric scalar!");
		threshold=(int) mxGetScalar(prhs[2]);
	}
	
	// check argument 4 (number of threads)
	int num_threads=1;
	if(nrhs>3 && mxIsEmpty(prhs[3])==0){
		if(mxIsNumeric(prhs[3])==0 || mxIsScalar(prhs[3])==0) mexErrMsgIdAndTxt("util:img:find_flares:inputNotNumericScalar", "Input 4 to find_flares is not a numeric scalar!");
		num_threads=(int) mxGetScalar(prhs[3]);
	}
	
	// check argument 5 (maximum number of peaks)
	int max_number=100;
	if(nrhs>4 && mxIsEmpty(prhs[4])==0){
		if(mxIsNumeric(prhs[4])==0 || mxIsScalar(prhs[4])==0) mexErrMsgIdAndTxt("util:img:find_flares:inputNotNumericScalar", "Input 5 to find_flares is not a numeric scalar!");
		max_number=(int) mxGetScalar(prhs[4]);
	}
	
	// check argument 6 (debug bit)
	int debug_bit=0;
	if(nrhs>5 && mxIsEmpty(prhs[5])==0){
		if(mxIsNumeric(prhs[5])==0 && mxIsLogical(prhs[5])==0 || mxIsScalar(prhs[5])==0) mexErrMsgIdAndTxt("util:img:find_flares:inputNotNumericScalar", "Input 6 to find_flares is not a numeric/logical scalar!");
		debug_bit=(int) mxGetScalar(prhs[5]);
	}
	
	// if(debug_bit) printf("cols= %d | rows= %d\n", cols, rows); 
	
	unsigned short ***array=new unsigned short**[num_threads]; // allocate space for the results from each thread
	for(int t=0;t<num_threads;t++){ 
		array[t]=new unsigned short*[max_number]; 
		for(int i=0;i<max_number;i++){ 
			array[t][i]=new unsigned short[4]; 
			for(int j=0;j<4;j++) array[t][i][j]=0; // allocate 4 places for x,y, frame number and peak value
		}
		
	}
	
	// array contains 4 values: x,y, frame num, peak value
	// array has 1st index for number of threads, 2nd for number of peaks, 3rd for those 4 values
	
	unsigned short *images=(unsigned short*) mxGetData(prhs[0]); // get the input image matrix
	bool *mask=mxGetLogicals (prhs[1]); // get the input mask
	
	// unsigned char new_thresh=(unsigned char) (threshold/256); 
	
	std::vector<std::thread> threads;
	
	int num_pages_per_thread=((int) pages)/num_threads;
	
	for(int t=0;t<num_threads;t++){
		
		int start_idx=t*num_pages_per_thread;
		
		if(t<num_threads-1){ // the first sections of the data go to asynchronuous threads
			int end_idx=start_idx+num_pages_per_thread;
			threads.push_back(std::thread(find_peaks, array[t], images, rows, cols, start_idx, end_idx, mask, threshold, max_number));  
		}
		else{ // last iteration goes on the main thread				
			int end_idx=(int) pages; // last iteration gets the additional round-off pages
			find_peaks(array[t], images, rows, cols, start_idx, end_idx, mask, threshold, max_number);
		}
		
	}
	
	for(int t=0;t<num_threads-1;t++){ // collect the results from all threads (this only completes after all threads are done)
		threads[t].join();
	}
	
	std::vector< std::vector<int> > results;
	int dist=3; // distance between existing peaks
	
	// consolidate the results
	for(int t=0;t<num_threads;t++){ // we can gain some speed improvement if we collect a separate vector for each thread, but this may be negligible
		
		for(int i=0;i<max_number;i++){
			
			unsigned short *new_result=array[t][i]; // get the row inside the large "array" so we can grab the results from that
			
			if(new_result[0]>0){
				
				// printf("i=%d | val= %05d | x= %04d | y= %04d | p= %02d\n", i, new_result[0], new_result[1], new_result[2], new_result[3]);
				
				int j=0; // need to query this to find if we broke the loop...
				for(j=0;j<results.size();j++){ // check the new result against all previously saved results
					if (new_result[3]==results[j][3] // check if two flares are on the same frame and close to each other
						&& abs((int) new_result[1] - (int) results[j][1])<dist 
						&& abs((int) new_result[2] - (int) results[j][2])<dist){ // I think we can do something more sophisticated with clustering but I don't have patience for this
						
						if(new_result[0]>results[j][0]) { // new results is brighter so it supercedes the old one
							for(int k=0;k<4;k++) results[j][k]=new_result[k]; // overwrite the old result
						} // new_result is better
						
						break; // skip the other results because we already see the new_result is too close to an existing result
						
					} // if close enough to an existing result
	
				} // for j (go over all previous results)
				
				if(j==results.size()) results.push_back(std::vector<int> (new_result, new_result+4) ); // new_result was not close to any existing result
				
			} // if new result is non-zero
			
		} // go over array i index
		
	} // go over array t index
	
	int N=results.size(); 
	
	if(debug_bit) for(int j=0;j<N;j++) printf("j= %02d | val= %d | x= %d | y= %d | p= %d \n", j, results[j][0],results[j][1],results[j][2],results[j][3]);
	
	plhs[0]=mxCreateDoubleMatrix(N, 4, mxREAL); // create the output matrix (in matlab)
	double *pos= mxGetPr(plhs[0]); // get the C type array 
	
	// plhs[1]=mxCreateDoubleMatrix(N, 1, mxREAL); // create the output matrix (in matlab)
	// double *peak= mxGetPr(plhs[1]); // get the C type array 
	
	for(int j=0;j<N;j++){
		for(int k=0;k<3;k++) pos[j+k*N]=results[j][k+1]; 
		pos[j+3*N]=results[j][0];
	} 
	
} // end mex function

void find_peaks(unsigned short **array, unsigned short *images, size_t rows, size_t cols, int start_idx, int end_idx, bool *mask, unsigned short threshold, int max_number){

	int counter=0;
	
	for (int p=start_idx; p<end_idx;p++){
		
		for(int i=0;i<rows*cols;i++){

			if(!mask[i]){
				
				//unsigned char value=images[p*rows*cols*2 + i*2 + 1]; // take only the significant byte
				
				unsigned short value=images[p*rows*cols + i]; // take one pixel value at a time
				
				if(value>=threshold){ 
				
					// save the information on this peak
					array[counter][0]=  value;
					array[counter][1]=  ( i/rows ) + 1; // x position
					array[counter][2]=  ( i%rows ) + 1; // y position
					array[counter][3]=  p + 1; // page or image number 
					
					counter++;
					if(counter>=max_number) return;
					
				}// if value>threshold
								
			} // if mask==0
			
		} //  for i (image position)
		
	}// for p (pages)
	
}