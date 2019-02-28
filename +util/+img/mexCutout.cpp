#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <math.h>

// usage M_out = mexCutout(M_in, positions, cut_size, [pad_value=0], [debug_bit=0])

int low_corner(int center, int size){
	
	int half_size = size/2;
	
	return center-half_size;
	
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){

    if(nrhs<3){ printf("USAGE: M_out = mexCutout(M_in, positions, cut_size, [pad_value=0], [debug_bit=0])\n"); return; }
	
	// check argument 0 
	if(mxIsNumeric(prhs[0])==0) mexErrMsgIdAndTxt("MATLAB:util:img:mexCutout:inputNotNumeric", "Input 1 to mexCutout is not numeric...");
	if( mxGetN(prhs[0])<2 || mxGetM(prhs[0])<2 || mxGetNumberOfDimensions(prhs[0])>3 ) mexErrMsgIdAndTxt("MATLAB:util:img:mexCutout:inputNotMatrix", "Must input a (2D or 3D) matrix to mexCutout!");

	// check argument 1
	if(mxIsNumeric(prhs[1])==0) mexErrMsgIdAndTxt("MATLAB:util:img:mexCutout:inputNotNumeric", "Input 2 to mexCutout is not numeric...");
	if( mxGetN(prhs[1])!=2) mexErrMsgIdAndTxt("MATLAB:util:img:mexCutout:inputWrongSize", "Second input must be a 2xN matrix of x and y center positions");
	
	// check argument 2
	if(mxIsNumeric(prhs[2])==0) mexErrMsgIdAndTxt("MATLAB:util:img:mexCutout:inputNotNumeric", "Input 3 to mexCutout is not numeric...");
	int cut_size=(int) mxGetScalar(prhs[2]);
		
	// check argument 3
	double pad_value=0;
	if(nrhs>3){
		if(mxIsNumeric(prhs[3])==0) mexErrMsgIdAndTxt("MATLAB:util:img:mexCutout:inputNotNumeric", "Input 4 to mexCutout is not numeric...");
		pad_value=mxGetScalar(prhs[3]);
	}
	
	int debug_bit=0;
	if(nrhs>4){
		if(mxIsNumeric(prhs[4])==0) mexErrMsgIdAndTxt("MATLAB:util:img:mexCutout:inputNotNumeric", "Input 5 to mexCutout is not numeric...");
		debug_bit=mxGetScalar(prhs[4]);
	}
	
	size_t rows=mxGetM(prhs[0]);
	size_t cols=mxGetN(prhs[0]);
	size_t pages = 1;
	if(mxGetNumberOfDimensions(prhs[0])>2){
		
		const mwSize *dims=mxGetDimensions(prhs[0]);
		rows=dims[0];
		cols=dims[1];
		pages=dims[2];
		
		
	}
	
	if(debug_bit) mexPrintf("rows= %d | cols= %d | pages= %d\n", rows, cols, pages);
	
	size_t num_cuts=mxGetM(prhs[1]);
	
	// get the cut positions into an array:
	double *pos_ptr=(double*) mxGetPr(prhs[1]);
	int *pos[2];
	
	for(int j=0; j<2; j++){
		
		pos[j]= (int*) mxCalloc(num_cuts, sizeof(int));
		
	}// for j
	
	for(int i=0;i<num_cuts; i++){
	
		for(int j=0; j<2; j++){
			
			pos[j][i]=(int) pos_ptr[num_cuts*j+i]-1; // minus one for conversion btw matlab indices and C indices
			
		}// for j
		
		if(debug_bit) mexPrintf("pos: %d %d\n", pos[0][i], pos[1][i]);
			
	}// for ii
	
	if(debug_bit) mexPrintf("input class: %s\n", mxGetClassName(prhs[0]));
	const size_t out_dims[4]={cut_size, cut_size, pages, num_cuts};
	
	if(mxIsClass(prhs[0], "double")){
		
		double *in_array=(double*) mxGetData(prhs[0]);
		double *out_array= (double*) mxCalloc(cut_size*cut_size*pages*num_cuts, sizeof(double)); // 
		for(int c=0; c<num_cuts; c++) for(int p=0; p<pages; p++) for(int j=0; j<cut_size; j++) for(int i=0; i<cut_size; i++) 
			out_array[c*cut_size*cut_size*pages+p*cut_size*cut_size+j*cut_size+i]=pad_value; // start by initializing everything to the pad value
	
		
		// copy the right part of the input matrix
		for(int c=0; c<num_cuts; c++){
			
			int x1= low_corner(pos[0][c], cut_size);
			int y1= low_corner(pos[1][c], cut_size);
			int cut_size_x = cut_size;
			int cut_size_y = cut_size;
			
			int push_x=0;
			int push_y=0;
			if(x1<0) push_x=-x1;
			if(y1<0) push_y=-y1;
			if(cut_size_x+x1>=cols) cut_size_x=cols-x1;
			if(cut_size_y+y1>=rows) cut_size_y=rows-y1;
					
			for(int p=0; p<pages; p++) for(int j=push_x; j<cut_size_x; j++) for(int i=push_y; i<cut_size_y; i++) 
				out_array[c*cut_size*cut_size*pages+p*cut_size*cut_size+j*cut_size+i]=in_array[rows*cols*p+x1*rows+y1+j*rows+i];
			
		}
		
		plhs[0]=mxCreateNumericArray(4, out_dims, mxDOUBLE_CLASS, mxREAL);
		mxSetData(plhs[0], out_array);		
				
		
	}
	else if(mxIsClass(prhs[0], "uint16")){
	
		if(isnan(pad_value) || pad_value<0) pad_value=0; // cannot have NaN or negative values in the uint16 matrix
		
		unsigned short int *in_array=(unsigned short int*) mxGetData(prhs[0]);
		unsigned short int *out_array= (unsigned short int*) mxCalloc(cut_size*cut_size*pages*num_cuts, sizeof(unsigned short int)); // 
		for(int c=0; c<num_cuts; c++) for(int p=0; p<pages; p++) for(int j=0; j<cut_size; j++) for(int i=0; i<cut_size; i++) 
			out_array[c*cut_size*cut_size*pages+p*cut_size*cut_size+j*cut_size+i]=pad_value; // start by initializing everything to the pad value
	
		// copy the right part of the input matrix
		for(int c=0; c<num_cuts; c++){
			
			int x1= low_corner(pos[0][c], cut_size);
			int y1= low_corner(pos[1][c], cut_size);
			int cut_size_x = cut_size;
			int cut_size_y = cut_size;
			
			int push_x=0;
			int push_y=0;
			if(x1<0) push_x=-x1;
			if(y1<0) push_y=-y1;
			if(cut_size_x+x1>=cols) cut_size_x=cols-x1;
			if(cut_size_y+y1>=rows) cut_size_y=rows-y1;
					
			for(int p=0; p<pages; p++) for(int j=push_x; j<cut_size_x; j++) for(int i=push_y; i<cut_size_y; i++) 
				out_array[c*cut_size*cut_size*pages+p*cut_size*cut_size+j*cut_size+i]=in_array[rows*cols*p+x1*rows+y1+j*rows+i];
			
		}

		plhs[0]=mxCreateNumericArray(4, out_dims, mxUINT16_CLASS, mxREAL);
		mxSetData(plhs[0], out_array);		
		
	}
	else{
		mexErrMsgIdAndTxt("MATLAB:util:img:mexCutout:inputTypeNotRecongnized", "Input must be double or uint16...");
    }
	
	for(int j=0; j<2; j++) mxFree(pos[j]);
	
	
				 
}