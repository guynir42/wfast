#include "SaveData.h"

// default values...
double *SaveData::mex_flag=0;
int SaveData::debug_bit=0;
int SaveData::deflate=0;
size_t SaveData::chunk_size=64;
bool SaveData::async_write=0;


SaveData::SaveData(const char *filename){
	
	strncpy(this->filename, filename, STRLN);
	
}

SaveData::~SaveData(){
	
	
}

void SaveData::write(){
	
	mex_flag[0]=1; // lock 
	// mex_flag[1]=0;
	writeData();
	writePars();
	mex_flag[0]=0; // release
	delete this; // must do a cleanup in the asynchronous function! 
	
}
void SaveData::parseVararginPairs(int N, const mxArray *vars[]){
		
	// if(nrhs%2==1) mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:vararginNotPairs", "varargin must be given as pairs!");
	for(int i=2; i<N;i+=2){
		
		if(mxIsChar(vars[i])==0) mexErrMsgIdAndTxt( "MATLAB:file:mex:mexWrite:inputNotString", "keyword must be a string.");
		const char *keyword=mxArrayToString(vars[i]);
		
		
		// if we decide this pair is worth parsing, use keyword and value variables:
		const mxArray *value=mxCreateDoubleScalar(1); // the positive approach
		if(i+1<N) value=vars[i+1];	
		if(mxIsEmpty(value)) continue; // just skip empty inputs
		
		// the third parameters tells MyMatrix if we want to deflate it (ignored if we are not using deflate at all)
		if(cs(keyword, "images")) images.input("images", value, 1); 
		else if(cs(keyword, "images_raw")) images_raw.input("images_raw", value, 1);
		else if(cs(keyword, "images_cal")) images_cal.input("images_cal", value, 1);
		else if(cs(keyword, "cutouts_raw")) cutouts_raw.input("cutouts_raw", value, 1);
		else if(cs(keyword, "cutouts_cal")) cutouts_cal.input("cutouts_cal", value, 1);
		else if(cs(keyword, "positions")) positions.input("positions", value, 0);
		else if(cs(keyword, "full_sum")) full_sum.input("full_sum", value, 1);
		else if(cs(keyword, "num_sum")) full_sum.attributes.push_back(MyAttribute("num_sum", value));
		else if(cs(keyword, "timestamps")) timestamps.input("timestamps", value, 0);
		else if(cs(keyword, "t_start", "file_start_datestring", 16)) timestamps.attributes.push_back(MyAttribute("t_start", value));
		else if(cs(keyword, "t_end", "file_write_datestring", 16)) timestamps.attributes.push_back(MyAttribute("t_end", value));
		else if(cs(keyword, "t_end_stamp", "file_write_timestamp", 17)) timestamps.attributes.push_back(MyAttribute("t_end_stamp", value));
		else if(cs(keyword, "psfs", 4)) psfs.input("psfs", value, 0);		
		else if(cs(keyword, "psf_sampling",4)) timestamps.attributes.push_back(MyAttribute("psf_sampling", value));
		else if(cs(keyword, "lightcurves")) lightcurves.input("lightcurves", value, 0);
		else if(cs(keyword, "pars", "parameters", "cell")) readParsCellArray(value);
		else if(cs(keyword, "buffers", "buffer wheel")) parseBufferWheelObject(value); // this can replace some of the definitions (not yet implemented...)
		else if(cs(keyword, "debug_bit")) debug_bit=(int) mxGetScalar(value);
		else if(cs(keyword, "deflate")) deflate=(int) mxGetScalar(value);
		else if(cs(keyword, "chunk")) chunk_size=(int) mxGetScalar(value);
		else if(cs(keyword, "async_write")) async_write=mxGetScalar(value);
		
	}

}

void SaveData::dataChecks(){
	
	if(positions.cols!=2) mexErrMsgIdAndTxt( "MATLAB:file:mex:mexWrite:positionsWrongDimension", "positions must be a two-column matrix!");
	// add all other data checks...
}

void SaveData::parseBufferWheelObject(const mxArray *buf){
	
}

void SaveData::readParsCellArray(const mxArray *cell){
	
	size_t N=0;
	if(mxIsEmpty(cell)) return; 
	if(mxGetN(cell)>1 && mxGetM(cell)>1) mexErrMsgIdAndTxt( "MATLAB:file:mex:SaveData:readParsCellArray", "Cell array must be 1D...");
	if(mxGetN(cell)>1) N=mxGetN(cell);
	if(mxGetM(cell)>1) N=mxGetM(cell);
	
	for(int i=0;i<N;i++){
		
		const mxArray *s=mxGetCell(cell, i);
		
		if(mxIsEmpty(s) || mxIsScalar(s)==0 || mxIsStruct(s)==0 || mxGetNumberOfFields(s)<1) 
			mexErrMsgIdAndTxt( "MATLAB:file:mex:SaveData:readParsCellArray", "Cell content must be a single struct with one field at least..");
		
		parameter_addresses_vector.push_back(std::string(mxArrayToString(mxGetFieldByNumber(s,0,0)))); // the address of the parameter object (should be first field)
		parameter_attributes_2D_vector.push_back(std::vector<MyAttribute>()); // an empty vector to be populated from the struct fields
		
		for(int j=1;j<mxGetNumberOfFields(s);j++){
			
			const mxArray *value=mxGetFieldByNumber(s,0,j);
			const char *name=mxGetFieldNameByNumber(s,j);
			parameter_attributes_2D_vector[i].push_back(MyAttribute(name, value));
			
		}
		
	}
	
	
}

void SaveData::setFilename(const char *name){
	
	snprintf(filename, STRLN, "%s", name);
	
}

void SaveData::print_help(){
	
	mexPrintf("Usage: writeHDF5(filename, mex_flag, varargin)\n");
	mexPrintf("mex_flag should be a 2 element vector, flag[0] means locked for writing, flag[1] is a debug counter.\n");
	
	mexPrintf("\n");
	mexPrintf("OPTIONAL ARGUMENTS \n");
	mexPrintf("--------------------------------\n");
	mexPrintf("Data inputs (must be given explicitely): \n");
	mexPrintf("-images: old input, try not to use this anymore. Can be uint16 or double\n");
	mexPrintf("-images_raw: new variable for uncalibrated files from camera. Usually uint16. \n");
	mexPrintf("-images_cal: full set of calibrated images. Usually double. \n");
	mexPrintf("cutouts_raw: when saving only cutouts around stars. Prefer to save raw cutouts, usually uint16. \n");	
	mexPrintf("cutouts_cal: when saving only cutouts around stars. Already calibrated, usually double. \n");
	mexPrintf("-positions: positions of centers of cutouts (x,y pairs).\n");
	mexPrintf("-timestamps: time of each frame. \n");
	mexPrintf("-psfs: if PSF data is given from WFS or from simulation.\n");
	mexPrintf("-lightcurves: photometry data for each cutout. \n");
	mexPrintf("-buffer: BufferWheel object used to extract add-on metadata.\n");
	mexPrintf("-pars_cell: cell array of address-struct pairs with metadata for this exposure. \n");
	mexPrintf("            (Use util.oop.save with output type 'struct').\n");
	
	mexPrintf("\n");
	mexPrintf("Metadata addons for the data-sets (or just give the BufferWheel object):");
	mexPrintf("-file_write_timestamp: when the file was sent to be written, relative to timestamps.\n");
	mexPrintf("-file_write_datestr: when the file was sent to be written, relative to system clock. \n");
	mexPrintf("-file_start_datestr: when the first exposure on file was started, relative to system clock. \n");
	mexPrintf("-psf_sampling: in units of pixels per lambda/D.\n");	
	
	mexPrintf("\n");
	mexPrintf("Parameters and controls: \n");
	mexPrintf("-debug_bit: controls level of output verbosity\n");
	mexPrintf("-async_write: uses separate thread to write data to file.\n");
	mexPrintf("-defalte: controls level of compression.\n");
	mexPrintf("-chunk: sets the size of tiles used for compression\n");
	mexPrintf("-file_type: can choose hdf5, fits, mat (override filename extension). \n");
	mexPrintf("            If no file type is given, extract from filename. \n");
	mexPrintf("            If no extension is given, the default is HDF5.\n");
	
}

void SaveData::handle_errors(const char *short_desc, const char *desc){
	
	char error_code[256];
	strncat("MATLAB:file:mex:mexWrite:", short_desc, 256);
	mexErrMsgIdAndTxt(error_code, desc);
	
}

bool SaveData::cs(const char *keyword, const char *compare_str, int num_letters){
	
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
	bool success=1;
	
	N=num_letters;
	if(strlen(str1)>N) N=strlen(str1); // number of letters to compare (minimum 3, or length of keyword).
	
	for(int i=0;i<N;i++){
		
		if(str1[i]!=str2[i]){
			success=0;
			break;
		}
		
	}
	
	return success;
	
}

bool SaveData::cs(const char *keyword, const char *str1, const char *str2, int num_letters){

	return cs(keyword, str1, num_letters) || cs(keyword, str2, num_letters);

}

bool SaveData::cs(const char *keyword, const char *str1, const char *str2, const char *str3, int num_letters){

	return cs(keyword, str1, num_letters) || cs(keyword, str2, num_letters) || cs(keyword, str3, num_letters);

}
