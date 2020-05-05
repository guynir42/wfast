#include "SaveData.h"

// default values...
double *SaveData::mex_flag=0;
int SaveData::debug_bit=0;
int SaveData::deflate=0;
size_t SaveData::chunk_size=64;
int SaveData::async_write=0;


SaveData::SaveData(const char *filename){
	
	strncpy(this->filename, filename, STRLN);
	
}

SaveData::~SaveData(){
	
	
}

void SaveData::write(){
	
	mex_flag[0]=1; // lock 
	// mex_flag[1]=0;
	writeData();
	writeHeader();
	mex_flag[0]=0; // release
	delete this; // must do a cleanup in the asynchronous function! 
	
}

void SaveData::parseVararginPairs(int N, const mxArray *vars[]){
		
	// if(nrhs%2==1) mexErrMsgIdAndTxt( "MATLAB:obs:mexWrite:vararginNotPairs", "varargin must be given as pairs!");
	for(int i=2; i<N;i+=2){
		
		if(mxIsChar(vars[i])==0 && mxIsStruct(vars[i])==0) mexErrMsgIdAndTxt( "MATLAB:file:mex:mexWrite:inputNotStringOrStruct", "Keyword must be a string or a struct.");
		if(mxIsStruct(vars[i])){ readStruct(vars[i]); i--; continue; } // skip this element and see if there are any more inputs to parse... 
		
		const char *keyword=mxArrayToString(vars[i]);
		
		// if we decide this pair is worth parsing, use keyword and value variables:
		const mxArray *value=mxCreateDoubleScalar(1); // the positive approach
		if(i+1<N) value=vars[i+1];	
		// if(mxIsEmpty(value)) continue; // just skip empty inputs
		
		// the third parameters tells MyMatrix if we want to deflate it (ignored if we are not using deflate at all)
		if(cs(keyword, "images")) images.input("images", value, 1); 
		
		else if(cs(keyword, "timestamps")) timestamps.input("timestamps", value, 0);
		else if(cs(keyword, "t_start", "file_start_datestring", 16)) timestamps.attributes.push_back(MyAttribute("t_start", value));
		else if(cs(keyword, "t_end", "file_write_datestring", 16)) timestamps.attributes.push_back(MyAttribute("t_end", value));
		else if(cs(keyword, "t_end_stamp", "file_write_timestamp", 17)) timestamps.attributes.push_back(MyAttribute("t_end_stamp", value));
		
		else if(cs(keyword, "cutouts")) cutouts.input("cutouts", value, 1);
		
		else if(cs(keyword, "positions")) positions.input("positions", value, 0);
		else if(cs(keyword, "obj_idx", "object_idx")) positions.attributes.push_back(MyAttribute("obj_idx", value));
		else if(cs(keyword, "coordinates")) coordinates.input("coordinates", value, 0);
		else if(cs(keyword, "magnitudes")) magnitudes.input("magnitudes", value, 0);
		else if(cs(keyword, "temperatures")) temperatures.input("temperature", value, 0);
		
		else if(cs(keyword, "cutouts_bg")) cutouts_bg.input("cutouts_bg", value, 1);
		else if(cs(keyword, "positions_bg")) positions_bg.input("positions_bg", value, 0);
		
		else if(cs(keyword, "stack")) stack.input("stack", value, 1);
		else if(cs(keyword, "num_sum")) stack.attributes.push_back(MyAttribute("num_sum", value));
		else if(cs(keyword, "psfs", 4)) psfs.input("psfs", value, 0);		
		else if(cs(keyword, "sampling_psf",4)) psfs.attributes.push_back(MyAttribute("sampling_psf", value));
		
		else if(cs(keyword, "fluxes")) fluxes.input("fluxes", value, 0); 
		else if(cs(keyword, "errors")) errors.input("errors", value, 0);
		else if(cs(keyword, "areas")) areas.input("areas", value, 0);
		else if(cs(keyword, "backgrounds")) backgrounds.input("backgrounds", value, 0);
		else if(cs(keyword, "variances")) variances.input("variances", value, 0);
		else if(cs(keyword, "offsets_x")) offsets_x.input("offsets_x", value, 0);
		else if(cs(keyword, "offsets_y")) offsets_y.input("offsets_y", value, 0);
		else if(cs(keyword, "centroids_x")) centroids_x.input("centroids_x", value, 0);
		else if(cs(keyword, "centroids_y")) centroids_y.input("centroids_y", value, 0);
		else if(cs(keyword, "widths")) widths.input("widths", value, 0);
		else if(cs(keyword, "bad_pixels")) bad_pixels.input("bad_pixels", value, 0);
		else if(cs(keyword, "flags")) flags.input("flags", value, 0);
		
		// add more names to the list if new data fields are added
		else if(cs(keyword, "header", "pars", "parameters", "cell")) readHeaderCellArray(value);
		else if(cs(keyword, "debug_bit")) debug_bit=(int) mxGetScalar(value);
		else if(cs(keyword, "deflate")) deflate=(int) mxGetScalar(value);
		else if(cs(keyword, "chunk")) chunk_size=(int) mxGetScalar(value);
		else if(cs(keyword, "async_write")) async_write=(int) mxGetScalar(value);
		
	}

}

void SaveData::readStruct(const mxArray *buf){
	
	// mexPrintf("reading data from struct...\n");
	
	std::vector<std::string> names;
	names.push_back("images");
	
	names.push_back("timestamps");
	names.push_back("t_start");
	names.push_back("t_end");
	names.push_back("t_end_stamp");
	
	names.push_back("cutouts");
	
	names.push_back("positions");
	names.push_back("obj_idx");
	names.push_back("coordinates");
	names.push_back("magnitudes");
	names.push_back("temperatures");
	
	names.push_back("cutouts_bg");
	names.push_back("positions_bg");
	
	names.push_back("stack");
	names.push_back("num_sum");
	names.push_back("psfs");
	names.push_back("sampling_psf");
	
	names.push_back("fluxes");
	names.push_back("errors"); 
	names.push_back("areas"); 
	names.push_back("backgrounds");
	names.push_back("variances"); 
	names.push_back("offsets_x"); 
	names.push_back("offsets_y"); 
	names.push_back("centroids_x"); 
	names.push_back("centroids_y");
	names.push_back("widths");
	names.push_back("bad_pixels"); 
	names.push_back("flags"); 
	
	// add more names to the list if new data fields are added
	
	mxArray *value=0;
	
	for(int i=0; i<names.size(); i++){
		
		// mexPrintf("i= %02d | name: %s\n", i, names[i].c_str());
		
		value=mxGetField(buf, 0, names[i].c_str());
		if(value){
			
			const char *keyword=names[i].c_str();
			
			if(cs(keyword, "images")) images.input("images", value, 1); 
			
			else if(cs(keyword, "timestamps")) timestamps.input("timestamps", value, 0);
			else if(cs(keyword, "t_start", "file_start_datestring", 16)) timestamps.attributes.push_back(MyAttribute("t_start", value));
			else if(cs(keyword, "t_end", "file_write_datestring", 16)) timestamps.attributes.push_back(MyAttribute("t_end", value));
			else if(cs(keyword, "t_end_stamp", "file_write_timestamp", 17)) timestamps.attributes.push_back(MyAttribute("t_end_stamp", value));
			
			else if(cs(keyword, "cutouts", 8)) cutouts.input("cutouts", value, 1);
			
			else if(cs(keyword, "positions", 9)) positions.input("positions", value, 0);
			else if(cs(keyword, "obj_idx", "object_idx")) positions.attributes.push_back(MyAttribute("obj_idx", value)); 
			else if(cs(keyword, "coordinates")) coordinates.input("coordinates", value, 0);
			else if(cs(keyword, "magnitudes")) magnitudes.input("magnitudes", value, 0);
			else if(cs(keyword, "temperatures")) temperatures.input("temperature", value, 0);
			
			else if(cs(keyword, "cutouts_bg", 8)) cutouts_bg.input("cutouts_bg", value, 1);
			else if(cs(keyword, "positions_bg", 9)) positions_bg.input("positions_bg", value, 0);
			
			else if(cs(keyword, "stack")) stack.input("stack", value, 1);
			else if(cs(keyword, "num_sum")) stack.attributes.push_back(MyAttribute("num_sum", value));
			else if(cs(keyword, "psfs", 4)) psfs.input("psfs", value, 0);		
			else if(cs(keyword, "sampling_psf",4)) psfs.attributes.push_back(MyAttribute("psf_sampling", value));
			
			else if(cs(keyword, "fluxes")) fluxes.input("fluxes", value, 0);
			else if(cs(keyword, "errors")) errors.input("errors", value, 0);
			else if(cs(keyword, "areas")) areas.input("areas", value, 0);
			else if(cs(keyword, "backgrounds")) backgrounds.input("backgrounds", value, 0);
			else if(cs(keyword, "variances")) variances.input("variances", value, 0);
			else if(cs(keyword, "offsets_x")) offsets_x.input("offsets_x", value, 0);
			else if(cs(keyword, "offsets_y")) offsets_y.input("offsets_y", value, 0);
			else if(cs(keyword, "centroids_x")) centroids_x.input("centroids_x", value, 0);
			else if(cs(keyword, "centroids_y")) centroids_y.input("centroids_y", value, 0);
			else if(cs(keyword, "widths")) widths.input("widths", value, 0);
			else if(cs(keyword, "bad_pixels")) bad_pixels.input("bad_pixels", value, 0);
			else if(cs(keyword, "flags")) flags.input("flags", value, 0);
			
			
			// add more names to the list if new data fields are added
		}
		
	}
}

void SaveData::dataChecks(){
	
	if(positions.cols!=2) mexErrMsgIdAndTxt( "MATLAB:file:mex:mexWrite:positionsWrongDimension", "positions must be a two-column matrix!");
	// add all other data checks...
}

void SaveData::readHeaderCellArray(const mxArray *cell){
	
	size_t N=0;
	if(mxIsEmpty(cell)) return; 
	if(mxGetN(cell)>1 && mxGetM(cell)>1) mexErrMsgIdAndTxt( "MATLAB:file:mex:SaveData:readHeaderCellArray", "Cell array must be 1D...");
	if(mxGetN(cell)>1) N=mxGetN(cell);
	if(mxGetM(cell)>1) N=mxGetM(cell);
	
	for(int i=0;i<N;i++){
		
		const mxArray *s=mxGetCell(cell, i);
		
		if(mxIsEmpty(s) || mxIsScalar(s)==0 || mxIsStruct(s)==0 || mxGetNumberOfFields(s)<1) 
			mexErrMsgIdAndTxt( "MATLAB:file:mex:SaveData:readHeaderCellArray", "Cell content must be a single struct with one field at least..");
		
		parameter_addresses_vector.push_back(std::string(mxArrayToString(mxGetFieldByNumber(s,0,0)))); // the address of the parameter object (should be first field)
		parameter_attributes_2D_vector.push_back(std::vector<MyAttribute>()); // an empty vector to be populated from the struct fields
		
		for(int j=1;j<mxGetNumberOfFields(s);j++){
			
			const mxArray *value=mxGetFieldByNumber(s,0,j);
			const char *name=mxGetFieldNameByNumber(s,j);
			
			if(mxIsCell(value)){
				
				for(int k=0;k<mxGetNumberOfElements(value); k++){
			
					char name_ext[256]; 
					snprintf(name_ext, 256, "%s{%d}", name, k+1);
					parameter_attributes_2D_vector[i].push_back(MyAttribute(name_ext, mxGetCell(value, k))); // this does not support imbedded structs or cells!
				}
				
			}
			else if(mxIsStruct(value)){
				
				for(int k=0;k<mxGetNumberOfFields(value); k++){
					
					const mxArray *sub_value=mxGetFieldByNumber(value,0,k);
					const char *sub_name=mxGetFieldNameByNumber(value,k);
					char name_ext[256]; 
					snprintf(name_ext, 256, "%s.%s", name, sub_name, k);
					parameter_attributes_2D_vector[i].push_back(MyAttribute(name_ext, sub_value)); // this does not support imbedded structs or cells!
				}
				
			}
			else
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
	mexPrintf("-images: full frame images.\n");
	mexPrintf("-timestamps: time of each frame. \n");
	mexPrintf("cutouts: when saving only cutouts around stars. Not calibrated, uint16. \n");
	mexPrintf("-positions: positions of centers of cutouts (x,y pairs).\n");
	mexPrintf("-coordinates: positions of centers of cutouts (RA,DE pairs).\n");
	mexPrintf("-magnitudes: of each star in each cutout, based on some catalog.\n");
	mexPrintf("-temperatures: of each star in each cutout, based on some catalog (degrees Kelvin).\n");
	mexPrintf("-stack: sum of  full franme images (calibrated)\n");
	mexPrintf("-psfs: if PSF data is given from WFS or from simulation.\n");
	mexPrintf("-fluxes: photometry data for each cutout. \n");
	// mexPrintf("-buffer: BufferWheel object used to extract add-on metadata.\n");
	mexPrintf("-header_cell: cell array of address-struct pairs with metadata for this exposure. \n");
	mexPrintf("            (Use util.oop.save with output type 'struct').\n");
	
	// mexPrintf("\n");
	// mexPrintf("Metadata addons for the data-sets (or just give the BufferWheel object):");
	// mexPrintf("-file_write_timestamp: when the file was sent to be written, relative to timestamps.\n");
	// mexPrintf("-file_write_datestr: when the file was sent to be written, relative to system clock. \n");
	// mexPrintf("-file_start_datestr: when the first exposure on file was started, relative to system clock. \n");
	// mexPrintf("-psf_sampling: in units of pixels per lambda/D.\n");	
	
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

bool SaveData::cs(const char *keyword, const char *str1, const char *str2, const char *str3, const char *str4, int num_letters){

	return cs(keyword, str1, num_letters) || cs(keyword, str2, num_letters) || cs(keyword, str3, num_letters) || cs(keyword, str4, num_letters);

}
