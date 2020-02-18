#include "mex.h"
#include "matrix.h"
#include <math.h>
#include "atcore.h"

#define STRLN 256

// USAGE: val = get(hndl, par_name, qualifier) 

// this is used to report errors back to matlab (definition at the end)
void throw_error(const char *description, const char *par_name=0, int error_code=0); 

// shortcut to just get an enum string
int get_string(AT_H hndl, const AT_WC *par_name, AT_WC *str, int length=STRLN);
int get_enum(AT_H hndl, const AT_WC *par_name, AT_WC *str);

// these are used to turn various data types into mxArray outputs
mxArray *output_scalar(double value); 
mxArray *output_scalar(AT_64 value); 
mxArray *output_scalar(AT_BOOL value); 
mxArray *output_vector(double *values, int length);
mxArray *output_string(const char *str); 
mxArray *output_string(const AT_WC *str);

// utility functions to compare strings
bool cs(const char *keyword, const char *compare_str, int num_letters=3);
bool cs(const char *keyword, const char *str1, const char *str2, int num_letters=3);
bool cs(const char *keyword, const char *str1, const char *str2, const char *str3, int num_letters=3);
bool parse_bool(mxArray *value);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){

	// if no inputs are given, just show the help file
	if (nrhs==0){
		const char *string[1]={"obs.cam.mex_new.get"};
		mxArray *array[1]={mxCreateCharMatrixFromStrings(1, string)};
		mexCallMATLAB(0,0,1,array,"help"); 
		return;
	}
	
	// check the 1st input is a camera handle
	if(mxIsScalar(prhs[0])==0 || mxIsNumeric(prhs[0])==0) throw_error("First input to get() must be numeric scalar (camera handle)"); 
	AT_H hndl = (AT_H) mxGetScalar(prhs[0]);
	
	// check 2nd input is a parameter name
	char key[STRLN]={0};
	if(mxIsChar(prhs[1])==0) throw_error("Second input to get() must be character array (name of feature/parameter)"); 
	mxGetString(prhs[1], key, STRLN); // copy the string data
	
	// check if there is a 3rd argument string
	char qualifier[STRLN]={0};
	if(nrhs>2){
		if(mxIsChar(prhs[2])==0) throw_error("Third input to get() must be character array (qualifier like 'min' or 'max')", key); 
		mxGetString(prhs[2], qualifier, STRLN); // copy the string data
	}
	
	int rc=0; // return value!
	
	// AT_WC *wide_str=new AT_WC[STRLN]();
	AT_WC wide_str[STRLN]={0};
	char str[STRLN]={0};
	double value=NAN; 
	AT_BOOL bool_value=0;	
	AT_64 int_value=0;
	
	if(cs(key, "running")){
		
		rc=AT_GetBool(hndl, L"CameraAcquiring", &bool_value); 
		plhs[0]=output_scalar(bool_value);
		
	}
	else if(cs(key, "temperature")){
	
		rc=AT_GetFloat(hndl, L"SensorTemperature", &value); 
		plhs[0]=output_scalar(value); 
	
	}
	else if(cs(key, "target temperature")){
		
		if (nrhs<3) rc=AT_GetFloat(hndl, L"TargetSensorTemperature", &value); 
		else if(cs(qualifier, "max")) rc=AT_GetFloatMax(hndl, L"TargetSensorTemperature", &value); 
		else if(cs(qualifier, "min")) rc=AT_GetFloatMin(hndl, L"TargetSensorTemperature", &value); 
		
		plhs[0]=output_scalar(value); 
	
	}
	else if(cs(key, "exposure time", "exp time")){
	
		if(nrhs<3) rc=AT_GetFloat(hndl, L"ExposureTime", &value); 
		else if(cs(qualifier, "max")) rc=AT_GetFloatMax(hndl, L"ExposureTime", &value); 
		else if(cs(qualifier, "min")) rc=AT_GetFloatMin(hndl, L"ExposureTime", &value); 
		
		plhs[0]=output_scalar(value); 
	
	}
	else if(cs(key, "frame rate")){
	
		if(nrhs<3) rc=AT_GetFloat(hndl, L"FrameRate", &value); 
		else if(cs(qualifier, "max")) rc=AT_GetFloatMax(hndl, L"FrameRate", &value); 
		else if(cs(qualifier, "min")) rc=AT_GetFloatMin(hndl, L"FrameRate", &value); 
		plhs[0]=output_scalar(value); 
	
	}
	else if(cs(key, "height", "AOI height", 4)){
	
		if(nrhs<3) rc=AT_GetInt(hndl, L"AOIHeight", &int_value);
		else if(cs(qualifier, "max")) rc=AT_GetIntMax(hndl, L"AOIHeight", &int_value);
		else if(cs(qualifier, "min")) rc=AT_GetIntMin(hndl, L"AOIHeight", &int_value);
		plhs[0]=output_scalar(int_value); 
	
	}
	else if(cs(key, "width", "AOI width", 4)){
	
		if(nrhs<3) rc=AT_GetInt(hndl, L"AOIWidth", &int_value);
		else if(cs(qualifier, "max")) rc=AT_GetIntMax(hndl, L"AOIWidth", &int_value);
		else if(cs(qualifier, "min")) rc=AT_GetIntMin(hndl, L"AOIWidth", &int_value);
		plhs[0]=output_scalar(int_value); 
	
	}
	else if(cs(key, "top", "AOI top", 4)){
	
		if(nrhs<3) rc=AT_GetInt(hndl, L"AOITop", &int_value);
		else if(cs(qualifier, "max")) rc=AT_GetIntMax(hndl, L"AOITop", &int_value);
		else if(cs(qualifier, "min")) rc=AT_GetIntMin(hndl, L"AOITop", &int_value);
		plhs[0]=output_scalar(int_value); 
	
	}
	else if(cs(key, "left", "AOI left", 4)){
	
		if(nrhs<3) rc=AT_GetInt(hndl, L"AOILeft", &int_value);
		else if(cs(qualifier, "max")) rc=AT_GetIntMax(hndl, L"AOILeft", &int_value);
		else if(cs(qualifier, "min")) rc=AT_GetIntMin(hndl, L"AOILeft", &int_value);
		plhs[0]=output_scalar(int_value); 
	
	}
	else if(cs(key, "stride", "AOI stride", 4)){
	
		rc=AT_GetInt(hndl, L"AOIStride", &int_value);
		plhs[0]=output_scalar(int_value); 
	
	}
	else if(cs(key, "size", "image size bytes", 4)){
	
		rc=AT_GetInt(hndl, L"ImageSizeBytes", &int_value);
		plhs[0]=output_scalar(int_value); 
	
	}
	else if(cs(key, "binning")){
		
		rc=get_enum(hndl, L"AOIBinning", wide_str);
		plhs[0]=output_string(wide_str); 
		
	}
	else if(cs(key, "vertical binning")){
	
		rc=AT_GetInt(hndl, L"AOIVBin", &int_value);
		plhs[0]=output_scalar(int_value); 
	
	}
	else if(cs(key, "horizontal binning")){
	
		rc=AT_GetInt(hndl, L"AOIHBin", &int_value);
		plhs[0]=output_scalar(int_value); 
	
	}
	else if(cs(key, "bytes per pixel")){
		
		rc=AT_GetFloat(hndl, L"BytesPerPixel", &value);
		plhs[0]=output_scalar(value); 
		
	}
	else if(cs(key, "baseline")){
	
		rc=AT_GetInt(hndl, L"Baseline", &int_value); 
		plhs[0]=output_scalar(int_value); 
	
	}
	else if(cs(key, "status")){
				
		rc=get_enum(hndl, L"CameraStatus", wide_str);
		
		plhs[0]=output_string(wide_str); 
		
	}
	else if(cs(key, "cycle mode")){
		
		if(nrhs<3){
			rc=get_enum(hndl, L"CycleMode", wide_str);
			plhs[0]=output_string(wide_str); 
		}
		else if(cs(qualifier, "options", "list")){ 
			snprintf(str, STRLN, "Fixed | Continuous"); 
			plhs[0]=output_string(str);
		}
		
	}
	else if(cs(key, "count")){
	
		rc=AT_GetInt(hndl, L"FrameCount", &int_value);
		plhs[0]=output_scalar(int_value); 
		
	}
	else if(cs(key, "shutter mode")){
		
		if(nrhs<3){ 
			rc=get_enum(hndl, L"ElectronicShutteringMode", wide_str);
			plhs[0]=output_string(wide_str);
		}
		else if(cs(qualifier, "options", "list")){ 
			snprintf(str, STRLN, "Rolling | Rolling - 100% Duty Cycle | Global | Global - 100% Duty Cycle"); 
			plhs[0]=output_string(str);
		}
		 
		
	}
	else if(cs(key, "trigger mode")){
	
		if(nrhs<3){ 
			rc=get_enum(hndl, L"TriggerMode", wide_str); 
			plhs[0]=output_string(wide_str); 
		}
		else if(cs(qualifier, "options", "list")){ 
			snprintf(str, STRLN, "Internal | Software"); 
			plhs[0]=output_string(str); 
		}
		
	}
	else if(cs(key, "encoding")){
		
		if(nrhs<3){ 
			rc=get_enum(hndl, L"PixelEncoding", wide_str);
			plhs[0]=output_string(wide_str); 
		}
		else if(cs(qualifier, "options", "list")){ 
			snprintf(str, STRLN, "Mono12 | Mono12Packed | Mono16 | Mono32 | Mono12CXP | Mono18"); 
			plhs[0]=output_string(str); 
		}
		
	}
	else if(cs(key, "noise filter")){
		
		rc=AT_GetBool(hndl, L"SpuriousNoiseFilter", &bool_value); 
		plhs[0]=output_scalar(bool_value);
		
	}
	else if(cs(key, "blemish correction")){
		
		rc=AT_GetBool(hndl, L"StaticBlemishCorrection", &bool_value); 
		plhs[0]=output_scalar(bool_value);
		
	}
	else if(cs(key, "pixel width", 7)){
		
		rc=AT_GetFloat(hndl, L"PixelWidth", &value); 
		plhs[0]=output_scalar(value); 
	
	}
	else if(cs(key, "pixel height", 7)){
		
		rc=AT_GetFloat(hndl, L"PixelHeight", &value); 
		plhs[0]=output_scalar(value); 
	
	}
	else if(cs(key, "readout time")){
		
		rc=AT_GetFloat(hndl, L"ReadoutTime", &value); 
		plhs[0]=output_scalar(value); 
	
	}
	else if(cs(key, "simulator mode")){
		
		if(nrhs<3){
			rc=get_enum(hndl, L"FrameGenMode", wide_str);
			plhs[0]=output_string(wide_str); 
		}
		else if(cs(qualifier, "options", "list")){ 
			snprintf(str, STRLN, "Off | ColumnCount | RowCount | FrameCount | FixedPixel"); 
			plhs[0]=output_string(str); 
		}
		
	}
	else if(cs(key, "cooling")){
		
		rc=AT_GetBool(hndl, L"SensorCooling", &bool_value); 
		plhs[0]=output_scalar(bool_value);
		
	}
	else if(cs(key, "fan mode")){
		
		if(nrhs<3){ 
			rc=get_enum(hndl, L"FanSpeed", wide_str);
			plhs[0]=output_string(wide_str); 
		}
		else if(cs(qualifier, "options", "list")){ 
			snprintf(str, STRLN, "Off | On | Low"); 
			plhs[0]=output_string(str); 
		}
		
	}
	else if(cs(key, "clock frequency", "frequency")){
		
		rc=AT_GetInt(hndl, L"TimestampClockFrequency", &int_value);
		plhs[0]=output_scalar(int_value);
		
	}
	else if(cs(key, "timestamp")){
	
		rc=AT_GetInt(hndl, L"TimestampClock", &int_value);
		plhs[0]=output_scalar(int_value);
	
	}
	else if(cs(key, "metadata")){
		
		rc=AT_GetBool(hndl, L"MetadataEnable", &bool_value); 
		plhs[0]=output_scalar(bool_value);
		
	}
	else if(cs(key, "layout")){
	
		rc=get_enum(hndl, L"AOILayout", wide_str);
		plhs[0]=output_string(wide_str); 
	
	}
	else if(cs(key, "name")){
	
		rc=get_string(hndl, L"CameraName", wide_str); 
				
		plhs[0]=output_string(wide_str); 
		
	}
	else if(cs(key, "information", "camera information")){
	
		// AT_WC *str1=new AT_WC[STRLN*10];
		
		AT_WC wide_long_str[STRLN*10]; 
		rc=get_string(hndl, L"CameraInformation", wide_long_str, STRLN*10); 
		
		plhs[0]=output_string(wide_long_str); 
		
		// delete[](str1);
	
	}
	else throw_error("Unknown parameter", key); 

	// delete[](str); 
	
	if(rc) throw_error("Get failed", key, rc); 
	
}

void throw_error(const char *description, const char *par_name, int error_code){

	char report[STRLN]={0};
	
	snprintf(report, STRLN, "%s", description); 
	if(par_name) snprintf(report, STRLN, "%s, parameter: %s", report, par_name);
	if(error_code) snprintf(report, STRLN, "%s, error code: %d", report, error_code); 
	
	mexErrMsgIdAndTxt("MATLAB:obs:cam:mex_new:get", report);

}

int get_string(AT_H hndl, const AT_WC *par_name, AT_WC *str, int length){

	// int length=0;
	// int rc=AT_GetStringMaxLength(hndl, par_name, &length);
	// printf("length= %d\n", length); 
	// if(rc) throw_error("Cannot get max string length"); // need to thing of some safer way to transmit this error
	
	int rc=0;
	
	rc=AT_GetString(hndl, par_name, str, length); 

	return rc;

}

int get_enum(AT_H hndl, const AT_WC *par_name, AT_WC *str){

	int index=0;
	int rc=AT_GetEnumIndex(hndl, par_name, &index); 
	if(rc) throw_error("Cannot get enum index");  // need to thing of some safer way to transmit this error

	rc=AT_GetEnumStringByIndex(hndl, par_name, index, str, STRLN); 

	return rc;

}

mxArray *output_scalar(double value){

	return mxCreateDoubleScalar(value);

}

mxArray *output_scalar(AT_64 value){

	return mxCreateDoubleScalar((double)value);

}

mxArray *output_scalar(AT_BOOL value){

	return mxCreateDoubleScalar((double)value);

}

mxArray *output_vector(double *values, int length){

	mxArray *array=mxCreateNumericMatrix(1,length,mxDOUBLE_CLASS, mxREAL);
	
	double *dbl_ptr=mxGetPr(array);
	
	for(int i=0;i<length;i++) dbl_ptr[i]=values[i];
	
	return array;
	
}

mxArray *output_string(const char *str){

	return mxCreateCharMatrixFromStrings(1, &str);

}

mxArray *output_string(const AT_WC *str){
	
	int N=wcslen(str); 
	
	char *str2=new char[N+1];
	
	wcstombs(str2, str, N+1);	

	mxArray *array=mxCreateCharMatrixFromStrings(1, (const char**) &str2);
	
	delete[](str2);
	
	return array;

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