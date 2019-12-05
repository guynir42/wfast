#include "mex.h"
#include "matrix.h"

#include "atcore.h"

#define STRLN 256

// USAGE: command(hndl, par_name)

// this is used to report errors back to matlab (definition at the end)
void throw_error(const char *description, const char *par_name=0, int error_code=0); 

// utility functions to compare strings
bool cs(const char *keyword, const char *compare_str, int num_letters=3);
bool cs(const char *keyword, const char *str1, const char *str2, int num_letters=3);
bool cs(const char *keyword, const char *str1, const char *str2, const char *str3, int num_letters=3);
bool parse_bool(mxArray *value);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){

	// check inputs!
	if (nrhs==0){
		const char *string[1]={"obs.cam.mex_new.set"};
		mxArray *array[1]={mxCreateCharMatrixFromStrings(1, string)};
		mexCallMATLAB(0,0,1,array,"help"); 
		return;
	}
		
	if(mxIsScalar(prhs[0])==0 || mxIsNumeric(prhs[0])==0) throw_error("First input to command() must be numeric scalar (camera handle)"); 
	AT_H hndl = (AT_H) mxGetScalar(prhs[0]);
	
	int rc=0; // return value!
	
	char key[STRLN];
	if(mxIsChar(prhs[1])==0) throw_error("Second input to command() must be character array (name of feature/parameter)"); 
	mxGetString(prhs[1], key, STRLN); // copy the string data
	
	if(cs(key, "start acquisition")){
		
		rc=AT_Command(hndl, L"AcquisitionStart"); 

	}
	else if(cs(key, "stop acquisition")){
	
		rc=AT_Command(hndl, L"AcquisitionStop"); 
	
	}
	else if(cs(key, "capture")){
		
		rc=AT_Command(hndl, L"SoftwareTrigger"); 
		
	}
	if(rc) throw_error("Command failed", key, rc); 
	
	
}
void throw_error(const char *description, const char *par_name, int error_code){

	char report[STRLN]={0};
	
	snprintf(report, STRLN, "%s", description); 
	if(par_name) snprintf(report, STRLN, "%s, parameter: %s", report, par_name);
	if(error_code) snprintf(report, STRLN, "%s, error code: %d", report, error_code); 
	
	mexErrMsgIdAndTxt("MATLAB:obs:cam:mex_new:command", report);

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