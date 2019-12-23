#include "mex.h"
#include "matrix.h"
#include <math.h>
#include "atcore.h"

#define STRLN 256

// USAGE: set(hndl, par_name, val) 

// this is used to report errors back to matlab (definition at the end)
void throw_error(const char *description, const char *par_name=0, int error_code=0); 



// utility functions to compare strings
bool cs(const char *keyword, const char *compare_str, int num_letters=3);
bool cs(const char *keyword, const char *str1, const char *str2, int num_letters=3);
bool cs(const char *keyword, const char *str1, const char *str2, const char *str3, int num_letters=3);
bool parse_bool(mxArray *value);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){

	// see if there are no inputs just print the help
	if (nrhs==0){
		const char *string[1]={"obs.cam.mex_new.set"};
		mxArray *array[1]={mxCreateCharMatrixFromStrings(1, string)};
		mexCallMATLAB(0,0,1,array,"help"); 
		return;
	}
	
	// verify first input is a camera handle
	if(mxIsScalar(prhs[0])==0 || mxIsNumeric(prhs[0])==0) throw_error("First input to get() must be numeric scalar (camera handle)"); 
	AT_H hndl = (AT_H) mxGetScalar(prhs[0]);
	
	// verify the second input is a string
	char key[STRLN]={0};
	if(mxIsChar(prhs[1])==0) throw_error("Second input to get() must be character array (name of feature/parameter)"); 
	mxGetString(prhs[1], key, STRLN); // copy the string data
		
	// get the value from the 3rd input or use the default NaN
	int rc=0; // return value!
	double value=0; // input value
	char str[STRLN]={0};
	AT_WC wide_str[STRLN]={0};
	
	if(nrhs<3 || mxIsEmpty(prhs[2])){ // default value
		value=NAN;
	}
	else if(mxIsScalar(prhs[2]) && mxIsNumeric(prhs[0])){ // take the input scalar
		value=mxGetScalar(prhs[2]); 
	}
	else if(mxIsChar(prhs[2])){ // take the input string
		mxGetString(prhs[2], str, STRLN);
		mbstowcs(wide_str, str, STRLN); 
	}
	else throw_error("Input 3 must be a numeric-scalar or char vector", key); // cannot read this parameter! 
	
	
	// choose which parameter to change! 
	if(cs(key, "temperature", "target temperature")){
	
		if(isnan(value)) throw_error("Input 3 must be a non-empty numeric-scalar!"); 
		rc=AT_SetFloat(hndl, L"TargetSensorTemperature", value); 
	
	}
	else if(cs(key, "exposure time", "exp time")){

		if(isnan(value)) value=0.03; // the W-FAST default
		rc=AT_SetFloat(hndl, L"ExposureTime", value); 
			
	}
	else if(cs(key, "frame rate")){
	
		if(isnan(value)) value=25; // the W-FAST default
		rc=AT_SetFloat(hndl, L"FrameRate", value); 
		
	}
	else if(cs(key, "binning")){
	
		if(isnan(value)) throw_error("Input 3 must be a non-empty character vector!\n Use '1x1' or '2x2' or '3x3' or '4x4' or '8x8'."); 
		rc=AT_SetEnumString(hndl, L"AOIBinning", wide_str);
	
	}
	else if(cs(key, "height", "AOI height", 4)){
		
		AT_64 int_value=0;
		
		if(isnan(value)){ // no input means get the default maximum height
			rc=AT_GetIntMax(hndl, L"AOIHeight", &int_value); // if no input given, set the AOI to default value
			if(rc) throw_error("Problem getting max height", key, rc); 
			value=(double) int_value;
		}
		
		rc=AT_SetInt(hndl, L"AOIHeight", (AT_64) value); 
		
	}
	else if(cs(key, "width", "AOI width", 4)){
		
		AT_64 int_value=0;
		if(isnan(value)){ // no input means get the default maximum width
			rc=AT_GetIntMax(hndl, L"AOIWidth", &int_value); // if no input given, set the AOI to default value
			if(rc) throw_error("Problem getting max width", key, rc); 
			value=(double) int_value;
		}
		
		rc=AT_SetInt(hndl, L"AOIWidth", (AT_64) value); 
	
	}
	else if(cs(key, "top", "AOI top", 4)){
		
		if(isnan(value)) value=1;
		
		rc=AT_SetInt(hndl, L"AOITop", (AT_64) value); 
		
	}
	else if(cs(key, "left", "AOI left", 4)){
		
		if(isnan(value)) value=1;
		
		rc=AT_SetInt(hndl, L"AOILeft", (AT_64) value); 
		
	}else if(cs(key, "cycle mode")){

		if(isnan(value)) throw_error("Input 3 must be a non-empty character vector!"); 
		rc=AT_SetEnumString(hndl, L"CycleMode", wide_str);
		
	}
	else if(cs(key, "shutter mode")){
		
		if(isnan(value)) throw_error("Input 3 must be a non-empty character vector!\n Use 'Global - 100% Duty Cycle' or 'Rolling - 100% Duty Cycle'."); 
		rc=AT_SetEnumString(hndl, L"ElectronicShutteringMode", wide_str);
		
	}
	else if(cs(key, "trigger mode")){
	
		if(isnan(value)) throw_error("Input 3 must be a non-empty character vector!\n Use 'Software' or 'Internal'."); 
		rc=AT_SetEnumString(hndl, L"TriggerMode", wide_str); 
		
	}
	else if(cs(key, "cycle mode")){
	
		if(isnan(value)) throw_error("Input 3 must be a non-empty character vector!\n Use 'Fixed' or 'Continuous'."); 
		rc=AT_SetEnumString(hndl, L"CycleMode", wide_str); 
		
	}
	else if(cs(key, "count")){
		rc=AT_SetInt(hndl, L"FrameCount", (AT_64) value); 
	}	
	else if(cs(key, "encoding")){
		
		if(isnan(value)) throw_error("Input 3 must be a non-empty character vector!\n Use 'Mono12', 'Mono12Packed', 'Mono16', 'Mono32, 'Mono12CXP', 'Mono18'."); 
		rc=AT_SetEnumString(hndl, L"PixelEncoding", wide_str);
		
	}
	else if(cs(key, "noise filter")){
		
		if(isnan(value)) value=0; // the default is not to use this feature!
		rc=AT_SetBool(hndl, L"SpuriousNoiseFilter", (bool) value); 
		
	}
	else if(cs(key, "blemish correction")){
		
		if(isnan(value)) value=0; // the default is not to use this feature!
		rc=AT_SetBool(hndl, L"StaticBlemishCorrection", (bool) value); 
		
	}
	else if(cs(key, "simulator mode")){
		
		// if(isnan(value)) swprintf(wide_str, STRLN, L"Off"); 
		if(isnan(value)) throw_error("Input 3 must be a non-empty character vector!\n Use 'Off', 'ColumnCount', 'RowCount', 'FrameCount' or 'FixedPixel'."); 
		rc=AT_SetEnumString(hndl, L"FrameGenMode", wide_str);
		
	}
	else if(cs(key, "cooling")){
		
		if(isnan(value)) value=1; // the default is not use this feature!
		rc=AT_SetBool(hndl, L"SensorCooling", (bool) value); 
				
	}
	else if(cs(key, "fan mode")){
		
		if(isnan(value)) throw_error("Input 3 must be a non-empty character vector!\n Use 'On', 'Off', or 'Low'. "); 
		rc=AT_SetEnumString(hndl, L"FanSpeed", wide_str);
				
	}
	else if(cs(key, "metadata")){

		if(isnan(value)) value=1; // the default is to use this feature!
		rc=AT_SetBool(hndl, L"MetadataEnable", (bool) value); 
				
	}
	else throw_error("Unknown parameter", key); 
	
	if(rc) throw_error("Set failed", key, rc); 
	
	
	
}

void throw_error(const char *description, const char *par_name, int error_code){

	char report[STRLN]={0};
	
	snprintf(report, STRLN, "%s", description); 
	if(par_name) snprintf(report, STRLN, "%s, parameter: %s", report, par_name);
	if(error_code) snprintf(report, STRLN, "%s, error code: %d", report, error_code); 
	
	mexErrMsgIdAndTxt("MATLAB:obs:cam:mex_new:set", report);

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

bool parse_bool(mxArray *value){ // if input is string of yes/no or on/off or number 1/0 output a boolean corresponding to the value
	
	if (mxIsEmpty(value)) return 0;
	if (mxIsScalar(value)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotScalar", "Input 1 to parse_bool must be a scalar!");
	if (mxIsLogical(value)) return mxGetScalar(value)!=0;
	if (mxIsNumeric(value)) return mxGetScalar(value)!=0;
	if (mxIsChar(value)) return cs(mxArrayToString(value), "yes", "on");
	return 0;
}