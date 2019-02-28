#include "MyAttribute.h"

MyAttribute::MyAttribute(){
	
	
}

MyAttribute::MyAttribute(const char *name){ // empty attribute constructor
	
	input(name);
	
}

void MyAttribute::input(const char *name){
	
	use_this=1;
	// att_name=std::string(name);	
	snprintf(att_name, STRLN, "%s", name);
		
}

MyAttribute::MyAttribute(const char *name, double value){ // scalar type attribute

	input(name, value);

}

void MyAttribute::input(const char *name, double value){

	input(name);
	is_scalar=1;
	scalar=value;
	
}

MyAttribute::MyAttribute(const char *name, double *data, int length){ // vector type attribute

	input(name, data, length);
	
}

void MyAttribute::input(const char *name, double *data, int length){

	input(name);
	is_vec=1;
	vec=std::vector<double>(data, data+length);

}

MyAttribute::MyAttribute(const char *name, const char *string){ // string type attribute

	input(name, string);

}

void MyAttribute::input(const char *name, const char *string){
	
	input(name);
	is_str=1;
	str=std::string(string);
	
}

MyAttribute::MyAttribute(const char *name, const mxArray *value){
	
	input(name, value);
	
}

void MyAttribute::input(const char *name, const mxArray *value){
	
	if(mxIsEmpty(value)) input(name);
	else if(mxIsChar(value)) input(name, mxArrayToString(value));	
	else if(mxIsScalar(value)) input(name, mxGetScalar(value)); 
	else if(mxIsNumeric(value)) input(name, mxGetPr(value), mxGetN(value));
	// other data types are not supported.
	else if(mxGetN(value)>1 && mxGetM(value)>1) mexErrMsgIdAndTxt( "MATLAB:file:mex:mexWrite:attributeMatrix", "Attribute cannot accept matices...");
	else if(mxIsObject(value)) mexErrMsgIdAndTxt( "MATLAB:file:mex:mexWrite:attributeObject", "Attribute cannot accept objects...");
	else if(mxIsStruct(value)) mexErrMsgIdAndTxt( "MATLAB:file:mex:mexWrite:attributeStruct", "Attribute cannot accept structs...");
	
}

void MyAttribute::printout(){
	
	if(is_empty()) return;
	printf("%s: ", att_name);
	if(is_scalar) printf("%f", scalar);
	if(is_vec) for(int i=0; i<vec.size(); i++) printf(" %f", vec[i]);
	if(is_str) printf("%s", str.c_str());
	printf("\n");
	
}

bool MyAttribute::is_empty(){
	
	return is_scalar==0 && is_vec==0 && is_str==0;
	
}