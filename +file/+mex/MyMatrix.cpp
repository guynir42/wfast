#include "MyMatrix.h"

MyMatrix::MyMatrix(){
	
}

MyMatrix::MyMatrix(const char *name, const mxArray *matrix, bool deflate){
	
	input(name, matrix, deflate);
	
}

void MyMatrix::clear(){
	
	//printf("clearing MyMatrix for %s\n", data_name); 
	
	matrix_uint16=0;
	matrix_double=0;
	matrix_float=0;
	
	ndims=0;
	dims=0;
	for(int i=0;i<4;i++) dims_c[i]=0;
	numel=0;	
	rows=0;
	cols=0;
	frames=1;
	cutouts=1;
	bits=0;
	use_deflate=0; 
	
	attributes.clear();
	

}

void MyMatrix::input(const char *name, const mxArray *matrix, bool deflate){
	
	clear(); // make sure that for any consecutive calls to "input", only the last is used... 
	
	if(mxIsEmpty(matrix)) return; // if we got an empty matrix, we can't really do anything except clear this object and keep going
	else if(mxIsClass(matrix, "uint16")){ matrix_uint16 = (unsigned short int*) mxGetData(matrix); bits=2; }
	else if(mxIsClass(matrix, "double")){ matrix_double = mxGetPr(matrix); bits=8; }
	else if(mxIsClass(matrix, "single")){ matrix_float = (float*) mxGetData(matrix); bits=4; }
	else mexErrMsgIdAndTxt( "MATLAB:file:mex:mexWrite:wrongDataType", "images must be uint16, float or double!");
	
	// data_name=std::string(name);
	snprintf(data_name, STRLN, "%s", name);
	
	ndims=(int) mxGetNumberOfDimensions(matrix);
	numel=mxGetNumberOfElements(matrix);
	dims=(size_t*) mxGetDimensions(matrix); // is this a memory leak? or worse: a pointer to an mxArray that is destroyed back in matlab???
	// dims_temp=(size_t*) mxGetDimensions(matrix);
	// for(int i=0;i<ndims;i++) dims[i]=dims_temp[i];
	rows=dims[0];
	cols=dims[1];
	if(ndims>=3) frames=dims[2];
	if(ndims>=4) cutouts=dims[3];
	
	for(int i=0;i<ndims;i++) dims_c[ndims-i-1]=dims[i];
	
	use_deflate=deflate;
	
}

bool MyMatrix::is_empty(){
	
	return matrix_uint16==0 && matrix_double==0 && matrix_float==0;
}

bool MyMatrix::is_double(){
	
	return matrix_double!=0;
}

bool MyMatrix::is_uint16(){
	
	return matrix_uint16!=0;
}

bool MyMatrix::is_float(){

	return matrix_float!=0;
}

void MyMatrix::printout(){
	
	if(is_empty()) return;
	printf("%s: %ld", data_name, dims[0]);
	for(int i=1; i<ndims; i++) printf("x%ld", dims[i]);
	if(is_double()) printf(" (double)\n");
	if(is_uint16()) printf(" (uint16)\n");
	if(is_float()) printf(" (single)\n");
	
	for(int i=0;i<attributes.size();i++) attributes[i].printout();
	
}