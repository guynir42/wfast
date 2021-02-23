#include "mex.h"
#include "matrix.h"
#include "photometry2.h"

#define NUM_GLOBAL_OBJECTS 10

Photometry photometry[NUM_GLOBAL_OBJECTS]; 

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){ // entry point from matlab

	/**********  check inputs! ***********/
	
	if (nrhs==0){ // show help if nargin==0
		
		const char *string[1]={"util.img.photometry2"};
		mxArray *array[1]={mxCreateCharMatrixFromStrings(1, string)};
		mexCallMATLAB(0,0,1,array,"help"); 
		return;
		
	}
	
	// read the input data and parameters
	if(mxIsEmpty(prhs[0])) return; // no input, then just return with all empty outputs...
	
	int idx=parseIndex(nrhs, prhs); // figure out which instance (index) of these objects to use
	
	photometry[idx].parseInputs(nrhs, prhs); // get varargin pairs etc
	photometry[idx].makeArrays(); // make the masks and indices vectors (unless they're already made)
	
	photometry[idx].run(); // do the actual calculations (possibly in multithread)

	char *field_names[7];
	for(int i=0;i<7;i++) field_names[i]=(char*) mxMalloc(STRLN);
	
	int counter=0;
	
	if(photometry[idx].pars.use_raw) snprintf(field_names[counter++], STRLN, "raw_photometry"); 
	if(photometry[idx].pars.use_forced) snprintf(field_names[counter++], STRLN, "forced_photometry"); 
	if(photometry[idx].pars.use_apertures) snprintf(field_names[counter++], STRLN, "apertures_photometry"); 
	if(photometry[idx].pars.use_gaussian) snprintf(field_names[counter++], STRLN, "gaussian_photometry"); 
	if(photometry[idx].pars.use_gaussian && photometry[idx].pars.use_forced) snprintf(field_names[counter++], STRLN, "forced_gaussian_photometry"); 
	snprintf(field_names[counter++], STRLN, "averages"); 
	snprintf(field_names[counter++], STRLN, "parameters"); 
	
	plhs[0]=mxCreateStructMatrix(1,1,counter, (const char**) field_names); 
	
	counter=0;
	if(photometry[idx].pars.use_raw) mxSetFieldByNumber(plhs[0], 0, counter++, photometry[idx].outputStruct(photometry[idx].output_raw)); 
	if(photometry[idx].pars.use_forced) mxSetFieldByNumber(plhs[0], 0, counter++, photometry[idx].outputStruct(photometry[idx].output_forced, photometry[idx].pars.num_radii)); 
	if(photometry[idx].pars.use_apertures) mxSetFieldByNumber(plhs[0], 0, counter++, photometry[idx].outputStruct(photometry[idx].output_apertures, photometry[idx].pars.num_radii)); 
	if(photometry[idx].pars.use_gaussian) mxSetFieldByNumber(plhs[0], 0, counter++, photometry[idx].outputStruct(photometry[idx].output_gaussian)); 
	if(photometry[idx].pars.use_gaussian && photometry[idx].pars.use_forced) mxSetFieldByNumber(plhs[0], 0, counter++, photometry[idx].outputStruct(photometry[idx].output_forced_gaussian)); 
	
	mxSetFieldByNumber(plhs[0], 0, counter++, photometry[idx].outputAverages()); 
	mxSetFieldByNumber(plhs[0], 0, counter++, photometry[idx].outputMetadataStruct());

	if(nlhs>1) plhs[1]=photometry[idx].outputArraysStruct();
	
}

int parseIndex(int nrhs, const mxArray *prhs[]){ // choose index for persistent photometry object (to reuse masks and output arrays)
	
	int index=1; // the default is to use the first (zero index after conversion to C++)
	mxArray *val=0;
		
	for(int i=1;i<nrhs;i+=2){ // parse varargin
		
		char key[STRLN];
		if(mxIsChar(prhs[i])==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotChar", "Input %d to photometry is not a string...", i+1);
		mxGetString(prhs[i], key, STRLN); // copy the string data
		
		if(i+1<nrhs) val=(mxArray*) prhs[i+1]; // if the varargin is odd numbered, leave val=0 as default
		
		if(cs(key, "index")){
			
			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for '%s' at input", key, i+2);
			if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not a numeric scalar...", i+2);
			index=(int) mxGetScalar(val);
			
		}
		
	}
	
	index--; // convert matlab one based index to C++ indexing
	
	if(index<0 || index>=NUM_GLOBAL_OBJECTS) 
		mexErrMsgIdAndTxt("MATLAB:util:img:photometry:indexOutOfBounds", "Index given (%d) exceeds the range of the global object array (%d)", index, NUM_GLOBAL_OBJECTS);
	
	return index; 
	
}

Photometry::Photometry(){ // class constructor

	// cannot instantiate an array in the header file!
	// consider changing these defaults at some point... 
	pars.ap_radii[0]=3;
	pars.ap_radii[1]=5;
	pars.ap_radii[2]=7; 
	
}

Photometry::~Photometry(){ // destructor cleans up intermidiate arrays
	
	if(pars.debug_bit) printf("now deleting photometry object!\n"); 
	
	if(pars.ap_radii) delete(pars.ap_radii); // make sure to free this little array
	
	deleteArrays();
	
}

mxArray *Photometry::outputStruct(float **output, int num_fluxes){ // wrap up the output matrices as a nice matlab style array

	// copy the field names into something persistent! 
	char *field_names[NUM_DATA_TYPES];
	for(int i=0;i<NUM_DATA_TYPES;i++){ 
		field_names[i]=(char*) mxMalloc(STRLN);
		snprintf(field_names[i], STRLN, data_types[i]);
	}
	
	mxArray *struct_array=mxCreateStructMatrix(1,1, NUM_DATA_TYPES, (const char**) field_names); 

	mwSize out_size[3]={(mwSize) output_size[0], (mwSize) output_size[1], (mwSize) num_fluxes}; 
	
	mxArray *matrix=0;
	for(int i=0;i<NUM_DATA_TYPES; i++){
		matrix=mxCreateNumericArray(3, out_size, mxSINGLE_CLASS, mxREAL);
		memcpy((float*)mxGetData(matrix), output[i], output_size[0]*output_size[1]*num_fluxes*sizeof(float));
		mxSetFieldByNumber(struct_array, 0, i, matrix);
	}
	
	return struct_array; 
	
}

mxArray *Photometry::outputAverages(){ // add a struct with the average offsets and widths

	char *field_names[3];
	for(int i=0;i<3;i++) field_names[i]=(char*)mxMalloc(STRLN);
	
	snprintf(field_names[0], STRLN, "offset_x"); 
	snprintf(field_names[1], STRLN, "offset_y"); 
	snprintf(field_names[2], STRLN, "width"); 
	
	mxArray *struct_array=mxCreateStructMatrix(1,1,3, (const char**) field_names);
	
	mxArray *matrix_offset_x=mxCreateNumericMatrix(dims[2], dims[3], mxSINGLE_CLASS, mxREAL); 
	memcpy((float*)mxGetData(matrix_offset_x), average_offset_x, num_cutouts*sizeof(float));
	mxSetFieldByNumber(struct_array, 0, 0, matrix_offset_x);
	
	mxArray *matrix_offset_y=mxCreateNumericMatrix(dims[2], dims[3], mxSINGLE_CLASS, mxREAL); 
	memcpy((float*)mxGetData(matrix_offset_y), average_offset_y, num_cutouts*sizeof(float));
	mxSetFieldByNumber(struct_array, 0, 1, matrix_offset_y);
	
	mxArray *matrix_width=mxCreateNumericMatrix(dims[2], dims[3], mxSINGLE_CLASS, mxREAL); 
	memcpy((float*)mxGetData(matrix_width), average_width, num_cutouts*sizeof(float));
	mxSetFieldByNumber(struct_array, 0, 2, matrix_width);
	
	return struct_array; 
	
}

mxArray *Photometry::outputMetadataStruct(){ // add a struct with some of the parameters and the different aperture masks used

	std::vector<std::string> field_names_vector;
	field_names_vector.push_back("cutout_size");
	field_names_vector.push_back("cutout_type"); 
	field_names_vector.push_back("aperture_radius"); 
	field_names_vector.push_back("forced_radius");
	field_names_vector.push_back("gauss_sigma");
	field_names_vector.push_back("annulus_radii");
	field_names_vector.push_back("iterations"); 
	field_names_vector.push_back("gain");
	field_names_vector.push_back("shift_resolution");
	
	
	const char **field_names=new const char*[field_names_vector.size()];
	for(int i=0;i<field_names_vector.size(); i++) field_names[i]=field_names_vector[i].c_str();
	mxArray *struct_array=mxCreateStructMatrix(1,1,field_names_vector.size(), field_names); 
	
	delete(field_names);
	
	int num=0;
	mxArray *array=0;
	
	// the size of the cutouts matrix
	array=mxCreateNumericMatrix(1,4,mxDOUBLE_CLASS, mxREAL);
	double *dbl_ptr=mxGetPr(array);
	for(int i=0;i<4;i++) dbl_ptr[i]=(double) dims[i];
	mxSetFieldByNumber(struct_array, 0, num++, array);
	
	// this contains the class of the cutouts (by default this is single)
	const char *cutouts_class[1]={"single"};
	array=mxCreateCharMatrixFromStrings(1, cutouts_class);
	mxSetFieldByNumber(struct_array, 0, num++, array);
	
	// this contains the aperture radius / radii used
	array=mxCreateNumericMatrix(1, pars.num_radii, mxDOUBLE_CLASS, mxREAL);
	dbl_ptr=mxGetPr(array);
	for(int i=0;i<pars.num_radii;i++) dbl_ptr[i]=pars.ap_radii[i];
	mxSetFieldByNumber(struct_array, 0, num++, array);
			
	// this contains the forced radius / radii used (this is the same as aperture)
	array=mxCreateNumericMatrix(1, pars.num_radii, mxDOUBLE_CLASS, mxREAL);
	dbl_ptr=mxGetPr(array);
	for(int i=0;i<pars.num_radii;i++) dbl_ptr[i]=pars.ap_radii[i];
	mxSetFieldByNumber(struct_array, 0, num++, array);
	
	// this contains the gauss sigma scalar
	array=mxCreateDoubleScalar(pars.gauss_sigma);
	mxSetFieldByNumber(struct_array, 0, num++, array);
	
	// this contains the inner and outer radii of the annulus
	array=mxCreateNumericMatrix(1, 2, mxDOUBLE_CLASS, mxREAL);
	dbl_ptr=mxGetPr(array);
	dbl_ptr[0]=pars.inner_radius;
	dbl_ptr[1]=pars.outer_radius;
	mxSetFieldByNumber(struct_array, 0, num++, array);
	
	// the number of iterations
	array=mxCreateDoubleScalar(pars.num_iterations);
	mxSetFieldByNumber(struct_array, 0, num++, array);
	
	// the gain used in estimating the errors
	array=mxCreateDoubleScalar(pars.gain);
	mxSetFieldByNumber(struct_array, 0, num++, array);
	
	// the shift resolution
	array=mxCreateDoubleScalar(pars.resolution);
	mxSetFieldByNumber(struct_array, 0, num++, array);
	
	return struct_array; 

}

mxArray *Photometry::outputArraysStruct(){ // add a struct with the actual masks and grid arrays (for debugging)

	std::vector<std::string> field_names_vector;
	field_names_vector.push_back("aperture_masks");
	field_names_vector.push_back("annulus_masks");
	field_names_vector.push_back("gaussian_masks");
	field_names_vector.push_back("aperture_indices");
	field_names_vector.push_back("annulus_indices");
	field_names_vector.push_back("X");
	field_names_vector.push_back("Y");
	field_names_vector.push_back("dx");
	field_names_vector.push_back("dy");
	
	const char **field_names=new const char*[field_names_vector.size()];
	for(int i=0;i<field_names_vector.size(); i++) field_names[i]=field_names_vector[i].c_str();
	
	mxArray *struct_array=mxCreateStructMatrix(1,1,field_names_vector.size(), field_names); 
	
	delete(field_names); field_names=0;
	
	// copy over the mask arrays
	mxArray *array=0;
	mwSize mask_dims[4]={dims[0],dims[1],num_shifts,pars.num_radii};
	int counter=0; 
	
	// save the aperture masks
	array=mxCreateNumericArray(4, mask_dims, mxSINGLE_CLASS, mxREAL);
	float *flt_ptr=(float*) mxGetData(array);
	memcpy(flt_ptr, apertures, N*num_shifts*pars.num_radii*sizeof(float));
	mxSetFieldByNumber(struct_array, 0, counter++, array);
	
	// save the annulus masks
	mask_dims[3]=1;
	array=mxCreateNumericArray(4, mask_dims, mxSINGLE_CLASS, mxREAL);
	flt_ptr=(float*) mxGetData(array);
	memcpy(flt_ptr, annulii, N*num_shifts*sizeof(float));
	mxSetFieldByNumber(struct_array, 0, counter++, array);
	
	// save the gaussian masks
	array=mxCreateNumericArray(4, mask_dims, mxSINGLE_CLASS, mxREAL);
	flt_ptr=(float*) mxGetData(array);
	memcpy(flt_ptr, gaussians, N*num_shifts*sizeof(float));
	mxSetFieldByNumber(struct_array, 0, counter++, array);
	
	// copy over the indices vectors
	array=outputIndicesVectors(aperture_indices, pars.num_radii); 
	mxSetFieldByNumber(struct_array, 0, counter++, array);
	
	array=outputIndicesVectors(annulus_indices); 
	mxSetFieldByNumber(struct_array, 0, counter++, array);

	// copy out the X/Y and dx/dy arrays
	array=mxCreateNumericMatrix(dims[0],dims[1], mxSINGLE_CLASS, mxREAL); 
	flt_ptr=(float*) mxGetData(array); 
	memcpy(flt_ptr, X, N*sizeof(float));
	mxSetFieldByNumber(struct_array, 0, counter++, array);
	
	array=mxCreateNumericMatrix(dims[0],dims[1], mxSINGLE_CLASS, mxREAL); 
	flt_ptr=(float*) mxGetData(array); 
	memcpy(flt_ptr, Y, N*sizeof(float));
	mxSetFieldByNumber(struct_array, 0, counter++, array);
	
	array=mxCreateNumericMatrix(shift_dims[0],shift_dims[1], mxSINGLE_CLASS, mxREAL); 
	flt_ptr=(float*) mxGetData(array); 
	memcpy(flt_ptr, dx, num_shifts*sizeof(float));
	mxSetFieldByNumber(struct_array, 0, counter++, array);	
	
	array=mxCreateNumericMatrix(shift_dims[0],shift_dims[1], mxSINGLE_CLASS, mxREAL); 
	flt_ptr=(float*) mxGetData(array); 
	memcpy(flt_ptr, dy, num_shifts*sizeof(float));
	mxSetFieldByNumber(struct_array, 0, counter++, array);
	
	return struct_array; 

}

mxArray *Photometry::outputIndicesVectors(std::vector<int> *vectors, int num_radii){ // produce a cell array with size num_shifts*num_radii (default=1), each with the list of indices for that mask
	
	mxArray *cell=mxCreateCellMatrix (num_shifts, num_radii);
	
	for(int k=0;k<num_radii;k++){ // go over aperture radii (or nothing if the default num_radii=1)
		
		for(int j=0;j<num_shifts; j++){ // go over the different shifts of this mask
			
			// get the data from the vector into a matlab vector
			mxArray *value=mxCreateNumericMatrix(vectors[j+k*num_shifts].size(), 1, mxDOUBLE_CLASS, mxREAL);
			double *data=mxGetPr(value); 
			for(int i=0; i<vectors[j+k*num_shifts].size(); i++) data[i]=vectors[j+k*num_shifts][i]+1; // copy the data and add one for the conversion to matlab indices
			
			mxSetCell(cell, j+k*num_shifts, value);
			
		} // for j
		
	} // for k
	
	return cell;
	
}

void Photometry::parseInputs(int nrhs, const mxArray *prhs[]){ // take the cutouts input and the matlab style varargin and parse into c++

	if(mxIsClass(prhs[0], "single")==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotFloat", "Input 1 to photometry is not a single/float array...");
	cutouts=(float*) mxGetData(prhs[0]);
	mwSize *new_dims=(mwSize*)mxGetDimensions(prhs[0]);
	ndims=(int) mxGetNumberOfDimensions(prhs[0]);	
	
	if(new_dims[0]!=dims[0] || new_dims[1]!=dims[1]){ 
		deleteArrays(); 
		for(int j=0;j<2; j++) dims[j]=new_dims[j];
		N=(int) (dims[0]*dims[1]); // dims[0] is the height while dims[1] is the width
	}
	
	if(ndims>2 && new_dims[2]!=dims[2]){
		deleteOutputs();
		dims[2]=new_dims[2];		
	}
	else if(ndims<=2 && dims[2]!=1){
		deleteOutputs();
		dims[2]=1;
	}
	
	if(ndims>3 && new_dims[3]!=dims[3]){
		deleteOutputs();
		dims[3]=new_dims[3];		
	}
	else if(ndims<=3 && dims[3]!=1){
		deleteOutputs();
		dims[3]=1;
	}
	
	output_size[0]=dims[2];
	output_size[1]=dims[3];
	
	num_cutouts=(int) (output_size[0]*output_size[1]); // how many values in each of the above arrays... 
	
	// before parsing varargin, must setup defaults
	Photometry::Parameters new_pars; 
	
	for(int i=1;i<nrhs;i+=2){ // parse varargin
		
		char key[STRLN];
		if(mxIsChar(prhs[i])==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotChar", "Input %d to photometry is not a string...", i+1);
		mxGetString(prhs[i], key, STRLN); // copy the string data
		
		mxArray *val=0;
		if(i+1<nrhs) val=(mxArray*) prhs[i+1]; // if the varargin is odd numbered, leave val=0 as default
		
		if(cs(key, "apertures", "radius", "radii")){

			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsNumeric(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not numeric...", i+2);
			
			// temporary values
			new_pars.ap_radii=(double*) mxGetData(val);
			new_pars.num_radii=(int) mxGetNumberOfElements(val);
			std::sort(new_pars.ap_radii, new_pars.ap_radii+new_pars.num_radii); 
			
		}
		else if(cs(key, "gaussian", "sigma", "gauss_sigma")){
			
			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not a numeric scalar...", i+2);
			
			new_pars.gauss_sigma=mxGetScalar(val); 
			
		}
		else if(cs(key, "annulus", "annuli")){

			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsNumeric(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not numeric...", i+2);
			
			int num=(int) mxGetNumberOfElements(val);
			
			if(mxGetNumberOfElements(val)==2){ // we got two inputs for annulus
				
				double *new_annuli=mxGetPr(val);
				new_pars.inner_radius=new_annuli[0];
				new_pars.outer_radius=new_annuli[1]; 
				
			}
			else if(mxIsScalar(val)){ // if no second input is given, assume it is 0 (=infinity)
				
				double new_annulus=mxGetScalar(val);
				new_pars.inner_radius=new_annulus; 
				new_pars.outer_radius=0; // the default is 0 (=infinity)
				
			}			
			else mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputWrongLength", "Input %d to photometry must be a scalar or 2-element vector...", i+1);
			
		}
		else if(cs("key", "gain")){
			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsEmpty(val)) continue;
			if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericSingle", "Input %d to photometry is not a numeric scalar!", i+2);
			pars.gain=mxGetScalar(val);
			
		}
		else if(cs("key", "scintillation_fraction")){
			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsEmpty(val)) continue;
			if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericSingle", "Input %d to photometry is not a numeric scalar!", i+2);
			pars.scintillation_fraction=mxGetScalar(val);
			
		}
		else if(cs(key, "resolution", "shift resolution")){
			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsEmpty(val)) continue;
			if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericSingle", "Input %d to photometry is not a numeric scalar!", i+2);
			
			new_pars.resolution=(int) mxGetScalar(val); 
			
		}
		else if(cs(key, "threads")){
			if(val==0 || mxIsEmpty(val)) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			else{
				if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry is not a numeric scalar...", i+2);
				pars.num_threads=(int) mxGetScalar(val);
			}
			
		}
		else if(cs(key, "iterations", "num_iterations", "number_iterations")){
			if(val==0 || mxIsEmpty(val)) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			else{
				if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry is not a numeric scalar...", i+2);
				pars.num_iterations=(int) mxGetScalar(val);
			}
			
		}
		else if(cs(key, "use_raw")){
			if(val==0 || mxIsEmpty(val)) pars.use_raw=1; // if no input, assume positive
			else{
				if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry is not a numeric scalar...", i+2);
				pars.use_raw=parse_bool(val);
			}
			
		}
		else if(cs(key, "use_centering_aperture")){
			if(val==0 || mxIsEmpty(val)) pars.use_centering_aperture=1; // if no input, assume positive
			else{
				if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry is not a numeric scalar...", i+2);
				pars.use_centering_aperture=parse_bool(val);
			}
			
		}
		else if(cs(key, "use_gaussian")){
			if(val==0 || mxIsEmpty(val)) pars.use_gaussian=1; // if no input, assume positive
			else{

				bool value=0;
				if(mxIsScalar(val) && (mxIsNumeric(val) || mxIsLogical(val))) value=(bool) mxGetScalar(val); 
				else if(mxIsChar(val)) value=parse_bool(val); 
				else mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry must be scalar or string", i+2);
				pars.use_gaussian=value;
				
			}
			
		}
		else if(cs(key, "use_apertures")){
			if(val==0 || mxIsEmpty(val)) pars.use_apertures=1; // if no input, assume positive
			else{
				
				bool value=0;
				if(mxIsScalar(val) && (mxIsNumeric(val) || mxIsLogical(val))) value=(bool) mxGetScalar(val); 
				else if(mxIsChar(val)) value=parse_bool(val); 
				else mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry must be scalar or string", i+2);
				pars.use_apertures=value;
				
			}
		}
		else if(cs(key, "use_forced")){
			if(val==0 || mxIsEmpty(val)) pars.use_forced=1; // if no input, assume positive
			else{
				
				bool value=0;
				if(mxIsScalar(val) && (mxIsNumeric(val) || mxIsLogical(val))) value=(bool) mxGetScalar(val); 
				else if(mxIsChar(val)) value=parse_bool(val); 
				else mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry must be scalar or string", i+2);
				pars.use_forced=value;
				
			}
			
		}
		else if(cs(key, "use_median", "median")){
			if(val==0 || mxIsEmpty(val)) pars.use_median=1; // if no second input, assume positive
			else{
				
				bool value=0;
				if(mxIsScalar(val) && (mxIsNumeric(val) || mxIsLogical(val))) value=(bool) mxGetScalar(val); 
				else if(mxIsChar(val)) value=parse_bool(val); 
				else mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry must be scalar or string", i+2);
				pars.use_median=value;
				
			}
			
		}
		else if(cs(key, "use_positives", "positives")){
			if(val==0 || mxIsEmpty(val)) pars.use_positives=1; // if no input, assume positive
			else{
				
				bool value=0;
				if(mxIsScalar(val) && (mxIsNumeric(val) || mxIsLogical(val))) value=(bool) mxGetScalar(val); 
				else if(mxIsChar(val)) value=parse_bool(val); 
				else mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry must be scalar or string", i+2);
				pars.use_positives=value;
				
			}
			
		}
		else if(cs(key, "debug_bit")){
			if(val==0 || mxIsEmpty(val)) pars.debug_bit=1; // if no input, assume positive
			else{
				if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry is not a numeric scalar...", i+2);
				pars.debug_bit=(int) mxGetScalar(val);
			}
			
		}
	}
	
	ingestParameters(new_pars); // take these parameters and figure out what can stay and what needs to change
	
	if(pars.debug_bit>1){ // check that all inputs have been received! 
		mexPrintf("cutouts: [");
		for(int i=0;i<ndims;i++){if(i>0) mexPrintf("x"); mexPrintf("%d", dims[i]); }
		mexPrintf("] (%s) \n", mxGetClassName(prhs[0]));
		
		mexPrintf("use: raw= %d | circle= %d | apertrues= %d | forced= %d | gaussian= %d \n", 
			pars.use_raw, pars.use_centering_aperture, pars.use_apertures, pars.use_forced, pars.use_gaussian); 
			
		mexPrintf("apertures: ");
		for(int i=0;i<pars.num_radii;i++) mexPrintf(" %f", pars.ap_radii[i]); 
				
		mexPrintf(" | annuli: %f %f | sigma= %f \n", pars.inner_radius, pars.outer_radius, pars.gauss_sigma);
		
		mexPrintf("iterations: %d | resolution= %d | num_threads= %d \n", pars.num_iterations, pars.resolution, pars.num_threads);
		
		mexPrintf("use_median= %d | use_positives= %d | debug_bit= %d\n", pars.use_median, pars.use_positives, pars.debug_bit); 
		
	}
	
}

void Photometry::ingestParameters(Photometry::Parameters new_pars){ // for each parameter that has changed, need to re-allocate mask/output array

	// update the aperture/forced radii
	if(new_pars.num_radii!=pars.num_radii){ // different number of radii, must recreate the radii array
		deleteOutputArray(output_apertures); 
		deleteOutputArray(output_forced); 
		deleteApertureMasks();
		pars.num_radii=new_pars.num_radii;
		if(pars.ap_radii) delete(pars.ap_radii); 
		pars.ap_radii=new double[pars.num_radii]; 
		for(int j=0;j<pars.num_radii; j++) pars.ap_radii[j]=new_pars.ap_radii[j]; 
	}
	else{ // can re-use the same array
		for(int j=0;j<pars.num_radii; j++){
			if(new_pars.ap_radii[j]!=pars.ap_radii[j]) deleteApertureMasks(); // at least one of the radii is different... need to re-allocate aperture masks!
			pars.ap_radii[j]=new_pars.ap_radii[j];
		}
	}

	// update the gaussian mask
	if(new_pars.gauss_sigma!=pars.gauss_sigma){
		
		deleteGaussianMasks();
		pars.gauss_sigma=new_pars.gauss_sigma;
		
	}
	
	// update the annulus inner and outer radius
	if(pars.inner_radius!=new_pars.inner_radius){
		deleteAnnulusMasks();
		pars.inner_radius=new_pars.inner_radius; 
	}
	
	if(new_pars.outer_radius>1e4 || new_pars.outer_radius!=new_pars.outer_radius) // if infinite or NaN
		new_pars.outer_radius=0; // the convention is to set the outer radius to zero when representing infinity 
	
	if(pars.outer_radius!=new_pars.outer_radius){
		deleteAnnulusMasks();
		pars.outer_radius=new_pars.outer_radius;
	} 
	
	// update the shift resolution
	if(new_pars.resolution!=pars.resolution){
		deleteGrids(); 
		deleteMasks(); 
		pars.resolution=new_pars.resolution;
	}
	
}

// INITIALIZATIONS

void Photometry::makeArrays(){ // generate all the one-time arrays in memory
	
	makeAllOutputs();
	makeGrids();
	makeMasks(); 
	
}

void Photometry::makeGrids(){ // make the X/Y coordinate grid and the dx/dy possible aperture shift matrices
	
	if(X==0 || Y==0 || dx==0 || dy==0) deleteGrids(); 
	
	X=new float[N];
	Y=new float[N];

	meshgrid(X,Y, dims); // make the static grid points
	
	for(int j=0;j<2;j++) shift_dims[j]=(dims[j]-1)*pars.resolution+1; // add "resolution"-1 number of points between each existing point
	num_shifts=(int) (shift_dims[0]*shift_dims[1]); 
	
	dx=new float[num_shifts];
	dy=new float[num_shifts];
	
	meshgrid(dx, dy, shift_dims, pars.resolution);
	
}

void Photometry::meshgrid(float *x, float *y, mwSize *dims, int resolution){ // fill the two matrices with dimensions "dims" with centered grid points (use resolution>1 to divide each pixel into multiple shifts)
	
	int NN=(int) (dims[0]*dims[1]); // the dimensions must be larger if you are using resolution>1
	
	for(int i=0;i<NN;i++){ // meshgrid
		
		x[i]=(float)(i/dims[0])-dims[1]/2;
		x[i]/=(float) resolution;
		y[i]=(float)(i%dims[0])-dims[0]/2;
		y[i]/=(float) resolution;
		
	}
	
}

void Photometry::makeAllOutputs(){ // generate all the output matrices needed
	 
	if(output_raw==0) allocateOutputArray(output_raw);
	initializeOutputArray(output_raw);
	
	if(output_apertures==0) allocateOutputArray(output_apertures, pars.num_radii);
	initializeOutputArray(output_apertures, pars.num_radii);
	
	if(output_forced==0) allocateOutputArray(output_forced, pars.num_radii);	
	initializeOutputArray(output_forced, pars.num_radii);
	
	if(output_gaussian==0) allocateOutputArray(output_gaussian); 
	initializeOutputArray(output_gaussian);
	
	if(output_forced_gaussian==0) allocateOutputArray(output_forced_gaussian); 
	initializeOutputArray(output_forced_gaussian);
	
	if(best_offset_x==0) best_offset_x=new float[num_cutouts]; 
	for(int j=0; j<num_cutouts; j++) best_offset_x[j]=NAN;
	
	if(best_offset_y==0) best_offset_y=new float[num_cutouts]; 
	for(int j=0; j<num_cutouts; j++) best_offset_y[j]=NAN;
	
	if(average_offset_x==0) average_offset_x=new float[num_cutouts];
	for(int j=0; j<num_cutouts; j++) average_offset_x[j]=NAN;
	
	if(average_offset_y==0) average_offset_y=new float[num_cutouts];
	for(int j=0; j<num_cutouts; j++) average_offset_y[j]=NAN;
	
	if(average_width==0) average_width=new float[num_cutouts];
	for(int j=0; j<num_cutouts; j++) average_width[j]=NAN;
	
}

void Photometry::allocateOutputArray(float **&output, int num_fluxes){ // create new memory for flux, background etc...

	output=new float*[NUM_DATA_TYPES];
	for(int i=0;i<NUM_DATA_TYPES;i++) output[i]=new float[output_size[0]*output_size[1]*num_fluxes]; 
	
}

void Photometry::initializeOutputArray(float **output, int num_fluxes){ // clear the content of existing memory for flux, background etc...

	for(int i=0;i<NUM_DATA_TYPES; i++){

		for(int k=0;k<num_fluxes; k++){
			
			for(int j=0;j<num_cutouts; j++) output[i][k*num_cutouts+j]=0;
			
		}
		
	}
	
}

void Photometry::makeMasks(){ // make a bank of masks to be used for different shifts and types of photometry

	if(apertures==0) makeApertureMasks(); 
	if(annulii==0) makeAnnulusMasks();
	if(gaussians==0) makeGaussianMasks(); 

}

void Photometry::makeApertureMasks(){ // make the aperture masks for e.g., forced photometry, with all possible shifts (dx/dy)
	
	apertures=new float[N*num_shifts*pars.num_radii]; 
	aperture_indices=new std::vector<int>[num_shifts*pars.num_radii]; // keep a running list of all the aperture indices for quick summing
	
	for(int k=0;k<pars.num_radii;k++){ // go over all aperture radii
	
		double radius=pars.ap_radii[k]; 
	
		for(int j=0; j<num_shifts;j++){ // go over the different shifts
		
			for(int i=0;i<N;i++){ // go over the pixels in each mask
				
				float r=(float) sqrt(pow(X[i]-dx[j],2)+pow(Y[i]-dy[j],2));
				
				float value=(float) (radius+0.5-r);
				if(value<1) value=0;
				else if(value>1) value=1;
				apertures[k*num_shifts*N + j*N + i] = value;
				
				if(k>0){ // remove values that are included in the previous aperture
					for(int m=1;m<=k; m++) if(apertures[(k-m)*num_shifts*N + j*N + i]>0) value=0; // don't include values that were in the previous apertures
				}
				
				if(value>0) aperture_indices[k*num_shifts + j].push_back(i); // keep a running list of all the aperture indices for quick summing
				
			} // for i
			
		} // for j 
	
	} // for k 
	
}

void Photometry::makeAnnulusMasks(){ // make a mask for the background annulus for each shift
	
	annulii=new float[N*num_shifts]; 
	annulus_indices=new std::vector<int>[num_shifts]; // keep a running list of all the annulus indices for quick summing
	
	float r1=pars.inner_radius;
	float r2=pars.outer_radius;
	if(r2<r1) r2=1e10; // replace for infinity
	
	for(int j=0; j<num_shifts;j++){ // go over the different shifts
	
		for(int i=0;i<N;i++){ // go over the pixels in each mask
			
			float r=(float) sqrt(pow(X[i]-dx[j],2)+pow(Y[i]-dy[j],2));
			float value=0;
			if(r>=r1 && r<r2) value=1;
			annulii[j*N + i] = value;
			
			if(value>0) annulus_indices[j].push_back(i); // keep a running list of all the annulus indices for quick summing			
			
		} // for i
		
	} // for j 
	
}

void Photometry::makeGaussianMasks(){ // make a gaussian weighted mask for each shift dx/dy
	
	gaussians=new float[N*num_shifts];
	
	for(int j=0; j<num_shifts;j++){ // go over the different shifts
	
		for(int i=0;i<N;i++){ // go over the pixels in each mask
			
			float r=(float) sqrt(pow(X[i]-dx[j],2)+pow(Y[i]-dy[j],2));
			
			gaussians[j*N + i] = (float) exp(-0.5*pow(r/pars.gauss_sigma,2)) / (2*PI*pars.gauss_sigma*pars.gauss_sigma); // normalized to unity
			
		} // for i
		
	} // for j 
	
}

// CLEANUP and DELETION

void Photometry::deleteArrays(){ // get rid of all the masks and output arrays (when destroying this object)

	// printf("deleting arrays...\n");

	deleteOutputs();
	deleteMasks();
	deleteGrids();
	
}

void Photometry::deleteGrids(){ // get rid of the meshgrid arrays
	
	if(X){ delete[](X); X=0; }
	if(Y){ delete[](Y); Y=0; }
	if(dx){ delete[](dx); dx=0; }
	if(dy){ delete[](dy); dy=0; }
	
}

void Photometry::deleteOutputs(){ // get rid of output arrays only

	// printf("deleting outputs...\n");

	deleteOutputArray(output_raw);
	deleteOutputArray(output_forced);
	deleteOutputArray(output_apertures);
	deleteOutputArray(output_gaussian);
	deleteOutputArray(output_forced_gaussian);
	
	if(best_offset_x){ delete[](best_offset_x); best_offset_x=0;}
	if(best_offset_y){ delete[](best_offset_y); best_offset_y=0;}
	
	if(average_offset_x){ delete[](average_offset_x); average_offset_x=0;}
	if(average_offset_y){ delete[](average_offset_y); average_offset_y=0;}
	if(average_width){ delete[](average_width); average_width=0;}

}

void Photometry::deleteOutputArray(float **&output){ // go over and deallocate the memory for one of the outputs
	
	if(output){
		
		for(int i=0;i<NUM_DATA_TYPES;i++){
			
			delete[](output[i]);
					
		}
		
		delete[](output); 
		output=0;
		
	}

}

void Photometry::deleteMasks(){ // go over and deallocate the memory for all the masks

	deleteApertureMasks();
	deleteAnnulusMasks();
	deleteGaussianMasks();

}

void Photometry::deleteApertureMasks(){ // get rid of the arrays for apertures
	
	if(apertures){ delete[](apertures); apertures=0; }
	if(aperture_indices){ delete[](aperture_indices); aperture_indices=0; }
}

void Photometry::deleteAnnulusMasks(){
	
	if(annulii){ delete[](annulii); annulii=0; }
	if(annulus_indices){ delete[](annulus_indices); annulus_indices=0; }

}

void Photometry::deleteGaussianMasks(){
		
	if(gaussians){ delete[](gaussians); gaussians=0; }
	
}

// ACTUAL CALCULATIONS!

void Photometry::run(){ // run photometry on all cutouts! the results are put into various sub-arrays of "output"

	if(pars.num_threads<=1){ // single threaded case
		for(int j=0;j<num_cutouts;j++){ // number of cutouts
			
			calculate(j);
			
		} // for j
	}
	else{ // multi-threaded case
		
		int step=num_cutouts/pars.num_threads;
		int current_idx=0;
		
		std::vector<std::thread> t;
				
		for(int i=0;i<pars.num_threads-1;i++){
			
			if(pars.debug_bit>2) mexPrintf("Sending a thread for indices %d to %d\n", current_idx, current_idx+step);
			t.push_back(std::thread(&Photometry::run_idx, this, current_idx, current_idx+step));
			current_idx+=step; 
			
		}
		
		if(pars.debug_bit>2) mexPrintf("Running on main thread indices %d to %d\n", current_idx, num_cutouts);
		run_idx(current_idx,num_cutouts); // run the remaining cutouts on the main thread! 
		
		for(int i=0;i<pars.num_threads-1;i++){
			t[i].join();
		}
		
	}
	
	if(pars.use_forced) runForced(); 
	
}

void Photometry::run_idx(int start_idx, int end_idx){ // run on subsets of cutouts, possibly in separate threads...

	for(int j=start_idx;j<end_idx;j++){ // partial list of cutouts
			
			calculate(j);
			
	} // for j

}

void Photometry::calculate(int j){ // do the actual calculations on a single cutout
	
	float *image=&cutouts[j*N]; // a pointer to the raw data
	
	float **output=output_raw; // make sure to get pointers to the raw output data	
	float *flux=output[IDX_FLUX];
	float *area=output[IDX_AREA];
	float *error=output[IDX_ERR];
	float *background=output[IDX_BG];
	float *variance=output[IDX_VAR];
	float *offset_x=output[IDX_DX];
	float *offset_y=output[IDX_DY];
	float *width=output[IDX_WD];
	float *bad_pixels=output[IDX_BAD];
	float *flag=output[IDX_FLAG];
	
	bad_pixels[j]=(float) countNaNs(image);
	
	int idx=getShiftIndex(0,0); // for the first order, raw photometry estimate of the background, use centered annulus
	int annulus_pixels=0; 
	float norm=0;
	float m1x=0;
	float m1y=0;
	float norm2=0;
	float m2x=0;
	float m2y=0;
	float mxy=0;
	
	if(bad_pixels[j]==N){ // image is all NaN!
		
		flux[j]=NAN;
		area[j]=NAN;
		error[j]=NAN;
		background[j]=NAN;
		variance[j]=NAN;
		offset_x[j]=NAN;
		offset_y[j]=NAN;
		width[j]=NAN;
		flag[j]=1;
		
		// make sure to set the other data types to NaN and flag them
		for(int i=0;i<NUM_DATA_TYPES;i++) output_gaussian[i][j]=NAN;
		output_gaussian[IDX_FLAG][j]=1;
		
		for(int i=0;i<NUM_DATA_TYPES;i++) output_forced[i][j]=NAN;
		output_gaussian[IDX_FLAG][j]=1;
		
		for(int k=0;k<pars.num_radii;k++){ 
			for(int i=0;i<NUM_DATA_TYPES;i++) output_apertures[i][k*num_cutouts + j]=NAN;
			output_apertures[IDX_FLAG][k*num_cutouts+j]=1;
		}
	
		return; // don't bother to do anything more with this image
		
	}
	
	best_offset_x[j]=0;
	best_offset_y[j]=0;
	
	if(pars.use_raw){

		flux[j]=sumArrays(image); // simple raw photometry! 
		
		area[j]=N-bad_pixels[j]; 
		
		int idx=getShiftIndex(0,0); // for the first order, raw photometry estimate of the background, use centered annulus
		annulus_pixels=countNonNaNsIndices(image, annulus_indices, idx); // how many non-NaN do we have in this annulus?
		
		// what happens if annulus_pixels is zero?
		if(pars.use_median) background[j]=medianIndices(image, annulus_indices, idx); 
		else background[j]=sumIndices(image, annulus_indices, idx)/annulus_pixels;
		
		variance[j]=sumIndices(image, background[j], image, background[j], annulus_indices, idx)/annulus_pixels; // subtract the average b/g then take the square, and sum all the values (then divide by number of pixels)
		
		error[j]=getError(flux[j]-area[j]*background[j], area[j]*variance[j], annulus_pixels*variance[j]); 
		
		// first moments:
		norm=flux[j]-area[j]*background[j]; // what if norm is zero or very close to zero?
		m1x=sumArrays(image, background[j], X)/norm;
		m1y=sumArrays(image, background[j], Y)/norm;
		
		offset_x[j]=m1x;
		offset_y[j]=m1y;

		// second moments
		norm2=sumArrays(image, background[j]); // the total flux, subtracting background, possibly ignoring negative pixels
		m2x=sumArrays(image, background[j], X, m1x, X, m1x)/norm; // I*(X-m1x)^2 / norm
		m2y=sumArrays(image, background[j], Y, m1y, Y, m1y)/norm; // I*(Y-m1y)^2 / norm
		mxy=sumArrays(image, background[j], X, m1x, Y, m1y)/norm; // I*(X-m1x)*(Y-m1y) / norm
		
		width[j]=getWidthFromMoments(m2x, m2y, mxy); 
		
		flag[j]=checkMoments(offset_x[j], offset_y[j], width[j]); 
		
		// keep track of the best estimate offsets
		best_offset_x[j]=offset_x[j];
		best_offset_y[j]=offset_y[j];
		
	}
	
	if(pars.use_centering_aperture){ // do one iteration with large aperture to improve centroid positions
		
		idx=getShiftIndex(best_offset_x[j],best_offset_y[j]); // the updated centroid

		annulus_pixels=countNonNaNsIndices(image, annulus_indices, idx); // how many non-NaN do we have in this annulus?

		float B=0; // temporary background value
		if(pars.use_median) B=medianIndices(image, annulus_indices, idx); 
		else B=sumIndices(image, annulus_indices, idx)/annulus_pixels; 
		
		float F=sumIndices(image, aperture_indices, idx); 
		int   A=countNonNaNsIndices(image, aperture_indices, idx);

		m1x=sumIndices(image, B, X, aperture_indices, idx)/(F-B*A); // the weighted centroid only on the selected pixels
		m1y=sumIndices(image, B, Y, aperture_indices, idx)/(F-B*A); // the weighted centroid only on the selected pixels
		
		// keep track of the best estimate offsets
		best_offset_x[j]=m1x;
		best_offset_y[j]=m1y;
		
	} // finished centering aperture
	
	if(pars.use_gaussian){// do gaussian photometry! 
		
		output=output_gaussian; // make sure to get pointers to the gaussian output data	
		flux=output[IDX_FLUX];
		area=output[IDX_AREA];
		error=output[IDX_ERR];
		background=output[IDX_BG];
		variance=output[IDX_VAR];
		offset_x=output[IDX_DX];
		offset_y=output[IDX_DY];
		width=output[IDX_WD];
		bad_pixels=output[IDX_BAD];
		flag=output[IDX_FLAG];
		
		for(int k=0; k<pars.num_iterations; k++){
		
			idx=getShiftIndex(best_offset_x[j],best_offset_y[j]); // the updated centroid
			
			annulus_pixels=countNonNaNsIndices(image, annulus_indices, idx); // how many non-NaN do we have in this annulus?
			
			// area[j]=1.0/sumArrays(&gaussians[idx*N],&gaussians[idx*N]); 
			area[j]=sumArrays(&gaussians[idx*N])/sumArrays(&gaussians[idx*N],&gaussians[idx*N]); 
			
			flux[j]=sumArrays(image, &gaussians[idx*N])*area[j]; // we must multiply this and all moments by area to compensate for the gaussian tapering
			
			if(pars.use_median) background[j]=medianIndices(image, annulus_indices, idx); 
			else background[j]=sumIndices(image, annulus_indices, idx)/annulus_pixels; 	
			
			norm=flux[j]-area[j]*background[j];
			
			m1x=sumArrays(image, background[j], X, &gaussians[idx*N])*area[j]/norm;
			m1y=sumArrays(image, background[j], Y, &gaussians[idx*N])*area[j]/norm;
			
			offset_x[j]=m1x;
			offset_y[j]=m1y;
			
			best_offset_x[j]=m1x;
			best_offset_y[j]=m1y;
			
		}// for k (iterations)
		
		// only get the variance/error/width after choosing the best spot
		variance[j]=sumIndices(image, background[j], image, background[j], annulus_indices, idx)/annulus_pixels; // subtract the average b/g then take the square, and sum all the values (then divide by number of pixels)
		
		error[j]=getError(flux[j]-area[j]*background[j], area[j]*variance[j], annulus_pixels*variance[j]);
		
		norm2=sumArrays(image, background[j], &gaussians[idx*N])*area[j]; 
		
		m2x=sumArrays(image, background[j], X, m1x, X, m1x, &gaussians[idx*N])*area[j]/norm2;
		m2y=sumArrays(image, background[j], Y, m1y, Y, m1y, &gaussians[idx*N])*area[j]/norm2;
		mxy=sumArrays(image, background[j], X, m1x, Y, m1y, &gaussians[idx*N])*area[j]/norm2;
		
		width[j]=getWidthFromMoments(m2x, m2y, mxy); 
				
		bad_pixels[j]=countNaNs(image, &gaussians[idx*N]); 
		
		flag[j]=checkMoments(offset_x[j], offset_y[j], width[j]); 

	} // finished doing gaussian photometry
	
	if(pars.use_apertures){ // do multi-aperture photometry (wedding cake)
		
		output=output_apertures; // make sure to get pointers to the apertures output data	
		flux=output[IDX_FLUX];
		area=output[IDX_AREA];
		error=output[IDX_ERR];
		background=output[IDX_BG];
		variance=output[IDX_VAR];
		offset_x=output[IDX_DX];
		offset_y=output[IDX_DY];
		width=output[IDX_WD];
		bad_pixels=output[IDX_BAD];
		flag=output[IDX_FLAG];
		
		idx=getShiftIndex(best_offset_x[j],best_offset_y[j]); // the updated centroid
		
		// get the background reading for all concentric apertures
		annulus_pixels=countNonNaNsIndices(image, annulus_indices, idx); // how many non-NaN do we have in this annulus?
		
		if(pars.use_median) background[j]=medianIndices(image, annulus_indices, idx); 
		else background[j]=sumIndices(image, annulus_indices, idx)/annulus_pixels; 	
		
		variance[j]=sumIndices(image, background[j], image, background[j], annulus_indices, idx)/annulus_pixels; // subtract the average b/g then take the square, and sum all the values (then divide by number of pixels)
		
		float *partial_m1x=new float[pars.num_radii];
		float *partial_m1y=new float[pars.num_radii];
		float *partial_m2x=new float[pars.num_radii];
		float *partial_m2y=new float[pars.num_radii];
		float *partial_mxy=new float[pars.num_radii];
		
		for(int k=0;k<pars.num_radii;k++){
		
			// the background/variance is the same for each aperture size, but 
			// it is easier to just copy these out for each aperture than to clip the data or leave zeros... 
			background[j+num_cutouts*k]=background[j];
			variance[j+num_cutouts*k]=variance[j];
		
			flux[j+num_cutouts*k]=sumIndices(image, aperture_indices, idx+num_shifts*k); 
			area[j+num_cutouts*k]=countNonNaNsIndices(image, aperture_indices, idx+num_shifts*k); 
		
			if(k>0){// sum the previous aperture into this one... 
			
				flux[j+num_cutouts*k]+=flux[j+num_cutouts*(k-1)];
				area[j+num_cutouts*k]+=area[j+num_cutouts*(k-1)];
			
			}
		
			partial_m1x[k]=sumIndices(image, background[j], X, aperture_indices, idx+num_shifts*k);
			partial_m1y[k]=sumIndices(image, background[j], Y, aperture_indices, idx+num_shifts*k);
			
			m1x=0;
			m1y=0;
			
			for(int m=0; m<k+1; m++){ // sum all the partial moments
				m1x+=partial_m1x[k-m];
				m1y+=partial_m1y[k-m];
			}
			
			norm=flux[j+num_cutouts*k]-area[j+num_cutouts*k]*background[j];
			
			m1x/=norm;
			m1y/=norm;
			
			offset_x[j+num_cutouts*k]=m1x;
			offset_y[j+num_cutouts*k]=m1y;
			
			m2x=0;
			m2y=0;
			mxy=0;
			
			partial_m2x[k]=sumIndices(image, background[j], X, m1x, X, m1x, aperture_indices, idx+num_shifts*k);
			partial_m2y[k]=sumIndices(image, background[j], Y, m1y, Y, m1y, aperture_indices, idx+num_shifts*k);
			partial_mxy[k]=sumIndices(image, background[j], X, m1x, Y, m1y, aperture_indices, idx+num_shifts*k);
			
			for(int m=0; m<k+1; m++){ // sum all the partial moments
				m2x+=partial_m2x[k-m];
				m2y+=partial_m2y[k-m];
				mxy+=partial_mxy[k-m];
			}
			
			float norm2=sumIndices(image, background[j], aperture_indices, idx+num_shifts*k);

			m2x/=norm2;
			m2y/=norm2;
			mxy/=norm2;
			
			error[j+num_cutouts*k]=getError(flux[j+num_cutouts*k]-area[j+num_cutouts*k]*background[j], area[j+num_cutouts*k]*variance[j], annulus_pixels*variance[j]);
			width[j+num_cutouts*k]=getWidthFromMoments(m2x, m2y, mxy); 

			if(k==0) bad_pixels[j]=aperture_indices[idx].size()-countNonNaNsIndices(image, aperture_indices, idx);
			else bad_pixels[j+num_cutouts*k]+=aperture_indices[idx+num_shifts*k].size()
			                                   -countNonNaNsIndices(image, aperture_indices, idx+num_shifts*k);
			
			flag[j+num_cutouts*k]=checkMoments(offset_x[j], offset_y[j], width[j]); 

		}// for k

		delete [] partial_m1x;
		delete [] partial_m1y;
		delete [] partial_m2x;
		delete [] partial_m2y;
		delete [] partial_mxy;
		
	}
	
	
	// forced photometry must be done after all of these are done! 
	
}

void Photometry::runForced(){
	
	if(pars.use_gaussian) calcFrameAverages(output_gaussian);
	else if(pars.use_apertures) calcFrameAverages(output_apertures, pars.num_radii);
	else if(pars.use_raw) calcFrameAverages(output_raw);
	else resetFrameAverages(); // no other photometry was used, just center the forced apertures in the middle of the cutouts (offsets=0)
	
	if(pars.num_threads<=1){
		for(int j=0;j<num_cutouts;j++){ // number of cutouts
			
			calculateForced(j);
			
		} // for j
	}
	else{
		
		int step=num_cutouts/pars.num_threads;
		int current_idx=0;
		
		std::vector<std::thread> t;
				
		for(int i=0;i<pars.num_threads-1;i++){
			
			if(pars.debug_bit>2) mexPrintf("Sending a thread for indices %d to %d\n", current_idx, current_idx+step);
			t.push_back(std::thread(&Photometry::runForced_idx, this, current_idx, current_idx+step));
			current_idx+=step; 
			
		}
		
		if(pars.debug_bit>2) mexPrintf("Running on main thread indices %d to %d\n", current_idx, num_cutouts);
		runForced_idx(current_idx,num_cutouts); // run the remaining cutouts on the main thread! 
		
		for(int i=0;i<pars.num_threads-1;i++){
			t[i].join();
		}
		
	}
	
}

void Photometry::runForced_idx(int start_idx, int end_idx){
	
	for(int j=start_idx;j<end_idx;j++){ // partial list of cutouts
			
			calculateForced(j);
			
	} // for j
	
}

void Photometry::calculateForced(int j){

	float *image=&cutouts[j*N]; // a pointer to the raw data

	// setup the correct outputs
	float **output=output_forced;
	float *flux=output[IDX_FLUX];
	float *area=output[IDX_AREA];
	float *error=output[IDX_ERR];
	float *background=output[IDX_BG];
	float *variance=output[IDX_VAR];
	float *offset_x=output[IDX_DX];
	float *offset_y=output[IDX_DY];
	float *width=output[IDX_WD];
	float *bad_pixels=output[IDX_BAD];
	float *flag=output[IDX_FLAG];
		
	int idx=getShiftIndex(average_offset_x[j], average_offset_y[j]); // the average centroid
	// get the background reading for all concentric apertures
		
	float annulus_pixels=countNonNaNsIndices(image, annulus_indices, idx); // how many non-NaN do we have in this annulus?
	
	if(pars.use_median) background[j]=medianIndices(image, annulus_indices, idx); 
	else background[j]=sumIndices(image, annulus_indices, idx)/annulus_pixels; 	
	
	variance[j]=sumIndices(image, background[j], image, background[j], annulus_indices, idx)/annulus_pixels; // subtract the average b/g then take the square, and sum all the values (then divide by number of pixels)
	
	float *partial_m1x=new float[pars.num_radii];
	float *partial_m1y=new float[pars.num_radii];
	float *partial_m2x=new float[pars.num_radii];
	float *partial_m2y=new float[pars.num_radii];
	float *partial_mxy=new float[pars.num_radii];
	
	for(int k=0;k<pars.num_radii;k++){
	
		// the background/variance is the same for each aperture size, but 
		// it is easier to just copy these out for each aperture than to clip the data or leave zeros... 
		background[j+num_cutouts*k]=background[j];
		variance[j+num_cutouts*k]=variance[j];
	
		flux[j+num_cutouts*k]=sumIndices(image, aperture_indices, idx+num_shifts*k); 
		area[j+num_cutouts*k]=countNonNaNsIndices(image, aperture_indices, idx+num_shifts*k); 
	
		if(k>0){
		
			flux[j+num_cutouts*k]+=flux[j+num_cutouts*(k-1)];
			area[j+num_cutouts*k]+=area[j+num_cutouts*(k-1)];
		
		}
	
		partial_m1x[k]=sumIndices(image, background[j], X, aperture_indices, idx+num_shifts*k);
		partial_m1y[k]=sumIndices(image, background[j], Y, aperture_indices, idx+num_shifts*k);
		
		float m1x=0;
		float m1y=0;
		
		for(int m=0; m<k+1; m++){ // sum all the partial moments
			m1x+=partial_m1x[k-m];
			m1y+=partial_m1y[k-m];
		}
		
		float norm=flux[j+num_cutouts*k]-area[j+num_cutouts*k]*background[j];
		
		m1x/=norm;
		m1y/=norm;
		
		offset_x[j+num_cutouts*k]=m1x;
		offset_y[j+num_cutouts*k]=m1y;
		
		float m2x=0;
		float m2y=0;
		float mxy=0;
		
		partial_m2x[k]=sumIndices(image, background[j], X, m1x, X, m1x, aperture_indices, idx+num_shifts*k);
		partial_m2y[k]=sumIndices(image, background[j], Y, m1y, Y, m1y, aperture_indices, idx+num_shifts*k);
		partial_mxy[k]=sumIndices(image, background[j], X, m1x, Y, m1y, aperture_indices, idx+num_shifts*k);
		
		for(int m=0; m<k+1; m++){ // sum all the partial moments
			m2x+=partial_m2x[k-m];
			m2y+=partial_m2y[k-m];
			mxy+=partial_mxy[k-m];
		}
		
		float norm2=sumIndices(image, background[j], aperture_indices, idx+num_shifts*k);

		m2x/=norm2;
		m2y/=norm2;
		mxy/=norm2;
		
		error[j+num_cutouts*k]=getError(flux[j+num_cutouts*k]-area[j+num_cutouts*k]*background[j], area[j+num_cutouts*k]*variance[j], annulus_pixels*variance[j]);
		width[j+num_cutouts*k]=getWidthFromMoments(m2x, m2y, mxy); 

		if(k==0) bad_pixels[j]=aperture_indices[idx].size()-countNonNaNsIndices(image, aperture_indices, idx);
		else bad_pixels[j+num_cutouts*k]+=aperture_indices[idx+num_shifts*k].size()
										   -countNonNaNsIndices(image, aperture_indices, idx+num_shifts*k);
		
		flag[j+num_cutouts*k]=checkMoments(offset_x[j], offset_y[j], width[j]); 

	}// for k

	delete [] partial_m1x;
	delete [] partial_m1y;
	delete [] partial_m2x;
	delete [] partial_m2y;
	delete [] partial_mxy;
	
	if(pars.use_gaussian){// do gaussian photometry! 
		
		output=output_forced_gaussian; // make sure to get pointers to the FORCED-gaussian output data	
		flux=output[IDX_FLUX];
		area=output[IDX_AREA];
		error=output[IDX_ERR];
		background=output[IDX_BG];
		variance=output[IDX_VAR];
		offset_x=output[IDX_DX];
		offset_y=output[IDX_DY];
		width=output[IDX_WD];
		bad_pixels=output[IDX_BAD];
		flag=output[IDX_FLAG];
		
		annulus_pixels=countNonNaNsIndices(image, annulus_indices, idx); // how many non-NaN do we have in this annulus?
		
		// area[j]=1.0/sumArrays(&gaussians[idx*N],&gaussians[idx*N]); 
		area[j]=sumArrays(&gaussians[idx*N])/sumArrays(&gaussians[idx*N],&gaussians[idx*N]); 
		
		flux[j]=sumArrays(image, &gaussians[idx*N])*area[j]; // we must multiply this and all moments by area to compensate for the gaussian tapering
		
		if(pars.use_median) background[j]=medianIndices(image, annulus_indices, idx); 
		else background[j]=sumIndices(image, annulus_indices, idx)/annulus_pixels; 	
		
		float norm=flux[j]-area[j]*background[j];
		
		float m1x=sumArrays(image, background[j], X, &gaussians[idx*N])*area[j]/norm;
		float m1y=sumArrays(image, background[j], Y, &gaussians[idx*N])*area[j]/norm;
		
		offset_x[j]=m1x;
		offset_y[j]=m1y;
		
		// best_offset_x[j]=m1x;
		// best_offset_y[j]=m1y;
			
		variance[j]=sumIndices(image, background[j], image, background[j], annulus_indices, idx)/annulus_pixels; // subtract the average b/g then take the square, and sum all the values (then divide by number of pixels)
		
		error[j]=getError(flux[j]-area[j]*background[j], area[j]*variance[j], annulus_pixels*variance[j]);
		
		float norm2=sumArrays(image, background[j], &gaussians[idx*N])*area[j]; 
		
		float m2x=sumArrays(image, background[j], X, m1x, X, m1x, &gaussians[idx*N])*area[j]/norm2;
		float m2y=sumArrays(image, background[j], Y, m1y, Y, m1y, &gaussians[idx*N])*area[j]/norm2;
		float mxy=sumArrays(image, background[j], X, m1x, Y, m1y, &gaussians[idx*N])*area[j]/norm2;
		
		width[j]=getWidthFromMoments(m2x, m2y, mxy); 
				
		bad_pixels[j]=countNaNs(image, &gaussians[idx*N]); 
		
		flag[j]=checkMoments(offset_x[j], offset_y[j], width[j]); 

	} // finished doing gaussian photometry
	
}

int Photometry::getShiftIndex(float x, float y){ // find the index closest to the specific shift value x and y in the shift matrices (dx and dy)
	
	// first make sure there are not shifts that are too big or too small
	if(x<dx[0]) x=dx[0];
	if(y<dy[0]) y=dy[0];
	if(x>dx[num_shifts-1]) x=dx[num_shifts-1];
	if(y>dy[num_shifts-1]) y=dy[num_shifts-1];
	
	// replace NaN values with zero (I don't have a better idea but this should be rare)
	if(x!=x) x=0;
	if(y!=y) y=0; 
	
	int idx_x=(int) round(pars.resolution*x+shift_dims[1]/2);
	int idx_y=(int) round(pars.resolution*y+shift_dims[0]/2);
	
	int idx=(int) (idx_x*shift_dims[0]+idx_y);
	
	return idx; 
	
}

float Photometry::getError(float reduced_flux, float aperture_variance, float background_variance){ // calculate the best estimate for the noise, including background noise, source noise, and scintillation
	
	float V=pow(reduced_flux,2)*pars.scintillation_fraction + reduced_flux*pars.gain + aperture_variance + background_variance;
	
	if(V<=0) return NAN;
	else return sqrt(V); 
	
}

bool Photometry::checkMoments(float offset_x, float offset_y, float width){ // returns 1 if there is a problem with the offsets or width
	
	if(offset_x>dims[1]/2) return 1;
	if(offset_y>dims[0]/2) return 1;
	if(width>dims[0]/2 || width>dims[1]/2) return 1;
	if(width<=0) return 1;
	if(offset_x!=offset_x || offset_y!=offset_y || width!=width) return 1; // if one of them are NaN
	
	return 0; // if all checks are not triggered, return with the ok flag
	
	
}

// GET WIDTHS AND OFFSETS

float Photometry::getWidthFromMoments(float m2x, float m2y, float mxy){ // calculate the eigenvalues of the 2nd moments and from that find the average width
// got this little nugget from: https://yutsumura.com/express-the-eigenvalues-of-a-2-by-2-matrix-in-terms-of-the-trace-and-determinant/

	float tr=m2x+m2y; // trace of the matrix
	float det=m2x*m2y - mxy*mxy; // determinant of the matrix
	
	if( (tr*tr-4*det) > 0){ // matrix has two different eigenvalues 

		float eig1=(tr-sqrt(tr*tr-4*det))/2; // the eigenvalues of this matrix
		float eig2=(tr+sqrt(tr*tr-4*det))/2; // help define the ellipse major/minor semi-axes

		 // printf("m2x= %f | m2y= %f | tr= %f | det= %f | eig1= %f | eig2=%f\n", m2x, m2y, tr, det, eig1, eig2);
	
		if (eig1<0 || eig2<0) return NAN; // std::nanf("") if one of the eigenvalues is negative we are in trouble
		// if (eig1<0 && eig2<0) return std::nanf(""); // if both of the eigenvalues is negative we are in trouble
		// else if (eig1<0) return sqrt(eig2); // use only one eigenvalue
		// else if (eig2<0) return sqrt(eig1); // use only one eigenvalue
		else return (sqrt(eig1)+sqrt(eig2))/2; // the sqrt of the eigenvalues are the minor/major semi-axes (return their average)
		
	} 
	else{ // matrix has two very similar eigenvalues, we should use only the trace
		
		float eig=tr/2; 
		
		 // printf("m2x= %f | m2y= %f | tr= %f | det= %f | eig= %f\n", m2x, m2y, tr, det);
	
		if (eig<0) return NAN; // std::nanf(""); 
		else return sqrt(eig); 
		
	}
	
}

void Photometry::calcFrameAverages(float **output, int num_radii){ // save the average offset_x, offset_y and width for each frame

	int num_frames=dims[2]; // assume dim 3 is the frame counter! 
	int num_stars=dims[3]; // which also means each star's light is continuous in memory
	
	for(int i=0;i<num_frames;i++){ // this is an inefficient way to traverse the memory but I don't think it is very important
	
		// these will be used to find the average offsets/width for this frame
		float sum_x=0;
		float sum_y=0;
		float sum_w=0;
		float sum_f=0; // total flux in each frame
		
		for(int j=0;j<num_stars;j++){
			
			if(output[IDX_FLAG][i+j*num_frames+(num_radii-1)*num_cutouts]==0){ // don't take any flagged frames

				// notice that the stars are picked from each frame
				float flux=output[IDX_FLUX][i+j*num_frames+(num_radii-1)*num_cutouts]; 
				float dx=output[IDX_DX][i+j*num_frames+(num_radii-1)*num_cutouts];
				float dy=output[IDX_DY][i+j*num_frames+(num_radii-1)*num_cutouts];
				float wd=output[IDX_WD][i+j*num_frames+(num_radii-1)*num_cutouts];
				
				if(flux==flux && dx==dx && dy==dy && wd==wd){ // do not add NaNs! 
					sum_f+=flux;
					sum_x+=dx*flux; 
					sum_y+=dy*flux; 
					sum_w+=wd*flux;
				}
				
			}
			
		}// for j (stars)
	
		if(sum_f>0){ 
			sum_x/=sum_f;
			sum_y/=sum_f;
			sum_w/=sum_f;
		}
		else{
			sum_x=0;
			sum_y=0;
			sum_w=0;
		}
	
		for(int j=0;j<num_stars;j++){
		
			average_offset_x[i+j*num_frames]=sum_x;
			average_offset_y[i+j*num_frames]=sum_y;
			average_width[i+j*num_frames]=sum_w;
		
		}// for j (stars)
	
	}// for i (frames)

}

void Photometry::resetFrameAverages(){
	
	int num_frames=dims[2]; // assume dim 3 is the frame counter! 
	int num_stars=dims[3]; // which also means each star's light is continuous in memory
	
	for(int i=0;i<num_frames;i++){ // this is an inefficient way to traverse the memory but I don't think it is very important
	
		for(int j=0;j<num_stars;j++){
		
			average_offset_x[i+j*num_frames]=0;
			average_offset_y[i+j*num_frames]=0;
			average_width[i+j*num_frames]=0;
		
		}// for j (stars)
	
	}// for i (frames)
	
}

// ARRAY ARITHMETIC

int Photometry::countNaNs(const float *array){ // count the number of NaNs in the array
	
	int num_nans=0;
	
	for(int i=0;i<N; i++) if(array[i]!=array[i]) { num_nans++; } // NOTE: the array[i]!=array[i] test is much MUCH faster than testing isnan(array[i]) or even _isnanf(array[i]). 
		
	return num_nans;
	
}

float Photometry::countNaNs(const float *array, const float *array2){ // count the number of NaNs weighted by a mask in array2
	
	float num_nans=0;
	
	for(int i=0;i<N; i++) if(array[i]!=array[i]) { num_nans+=array2[i]; } // NOTE: the array[i]!=array[i] test is much MUCH faster than testing isnan(array[i]) or even _isnanf(array[i]). 
		
	return num_nans;
	
}

float Photometry::sumArrays(const float *array1){ // just the sum of all non-NaN pixels in this array
	
	float S=0;
	
	for(int i=0;i<N;i++) if(array1[i]==array1[i]) S+=array1[i];
	
	return S;
	
}

float Photometry::sumArrays(const float *array1, float offset1){
	
	float S=0;
	
	if(pars.use_positives)
		for(int i=0;i<N;i++) if(array1[i]>offset1) S+=array1[i]-offset1;
	else
		for(int i=0;i<N;i++) if(array1[i]==array1[i]) S+=array1[i]-offset1;
	
	return S;
	
}

float Photometry::sumArrays(const float *array1, const float *array2){ // sum of all non-NaN values of array1*array2
	
	float S=0;
	
	for(int i=0;i<N;i++) 
		if(array1[i]==array1[i])
			S+=array1[i]*array2[i];
	
	return S;
	
}

float Photometry::sumArrays(const float *array1, float offset1, const float *array2){ // sum of all non-NaN values of array1*array2
	
	float S=0;
	
	for(int i=0;i<N;i++) if(array1[i]==array1[i]) S+=(array1[i]-offset1)*array2[i];
	
	return S;
	
}

float Photometry::sumArraysPos(const float *array1, float offset1, const float *array2){ // sum of all non-NaN values of array1*array2
	
	float S=0;
	
	for(int i=0;i<N;i++){ 
	
		float value = array1[i]-offset1;
		if(pars.use_positives){
			if(array1[i]>0) S+=value*array2[i];
		}
		else{
			if(array1[i]==array1[i]) S+=value*array2[i];
		}
	}
	
	return S;
	
}

float Photometry::sumArrays(const float *array1, float offset1, const float *array2, const float *array3){ // sum of all non-NaN values of array1*array2*array3
	
	float S=0;
	
	for(int i=0;i<N;i++) 
		if(array1[i]==array1[i])
			S+=(array1[i]-offset1)*array2[i]*array3[i];
	
	return S;
	
}

float Photometry::sumArrays(const float *array1, float offset1, const float *array2, float offset2, const float *array3, float offset3){ // sum of all non-NaN values of array1*array2*array3
	
	// use this for calculating the 2nd moments with simple (raw) photometry: 
	// array1 is the image, offset1 is the background,
	// array2 and array3 are X or Y, their offsets are the 1st moments. 
	
	float S=0;
	
	if(pars.use_positives){
		for(int i=0;i<N;i++) 
			if(array1[i]>offset1) // when true, array1-offset1 is positive (and non-NaN)
				S+=(array1[i]-offset1)*(array2[i]-offset2)*(array3[i]-offset3);
	}
	else{
			
		for(int i=0;i<N;i++) 
			if(array1[i]==array1[i]) // only sum non-NaN numbers
				S+=(array1[i]-offset1)*(array2[i]-offset2)*(array3[i]-offset3);
	}
	
	return S;
	
}

float Photometry::sumArrays(const float *array1, float offset1, const float *array2, float offset2, const float *array3, float offset3, const float *array4){ // sum of all non-NaN values of array1*array2*array3*array4
	
	// use this for calculating the 2nd moments with gaussian weights: 
	// array1 is the image, offset1 is the background,
	// array2 and array3 are X or Y, their offsets are the 1st moments, 
	// array4 is the gaussian (so it doesn't need an offset). 
	
	float S=0;
	
	// this is used to calculate 2nd moment, so we can have two options:
	if(pars.use_positives){
		for(int i=0;i<N;i++) 
			if(array1[i]>offset1) // when true, array1-offset1 is positive
				S+=(array1[i]-offset1)*(array2[i]-offset2)*(array3[i]-offset3)*array4[i];
	
	}
	else{ 
		for(int i=0;i<N;i++) 

			if(array1[i]==array1[i])
				S+=(array1[i]-offset1)*(array2[i]-offset2)*(array3[i]-offset3)*array4[i];
	
	}
	
	return S;
	
}

int Photometry::countNonNaNsIndices(const float *array1, const std::vector<int> *vector, int idx){ // sum the number of non-NaN values in array1 on the indices in vector[idx]

	int num=0;
	
	for(int i=0; i<vector[idx].size(); i++){
		
		float value=array1[vector[idx][i]]; 
		
		if(value==value) num++;
		
	}

	return num;

}

float Photometry::sumIndices(const float *array1, const std::vector<int> *vector, int idx){ // sum the values of array1 on the indices in vector[idx]
	
	float sum=0;

	for(int i=0; i<vector[idx].size(); i++){
		
		float value=array1[vector[idx][i]]; 
		
		if(value==value) sum+=value;
		
	}

	return sum;
	
}

float Photometry::sumIndices(const float *array1, float offset1, const std::vector<int> *vector, int idx){
	
	
	float sum=0;
	
	for(int i=0; i<vector[idx].size(); i++){
		
		float value1=array1[vector[idx][i]]-offset1; 
		
		if(pars.use_positives){
			if(value1>0) sum+=value1;
		}
		else{
			if(value1==value1) sum+=value1;
		}
		
	}

	return sum;
	
}

float Photometry::sumIndices(const float *array1, const float *array2, const std::vector<int> *vector, int idx){
	
	
	float sum=0;
	
	for(int i=0; i<vector[idx].size(); i++){
		
		float value1=array1[vector[idx][i]]; 
		float value2=array2[vector[idx][i]];

		if(value1==value1) sum+=value1*value2;
		
	}

	return sum;
	
}

float Photometry::sumIndices(const float *array1, float offset1, const float *array2, const std::vector<int> *vector, int idx){
	
	
	float sum=0;
	
	for(int i=0; i<vector[idx].size(); i++){
		
		float value1=array1[vector[idx][i]]-offset1; 
		float value2=array2[vector[idx][i]];
		if(value1==value1) sum+=value1*value2;
		
	}

	return sum;
	
}

float Photometry::sumIndices(const float *array1, float offset1, const float *array2, float offset2, const std::vector<int> *vector, int idx){
	
	
	float sum=0;
	
	for(int i=0; i<vector[idx].size(); i++){
		
		float value1=array1[vector[idx][i]]-offset1; 
		float value2=array2[vector[idx][i]]-offset2;
		if(value1==value1) sum+=value1*value2;
		
	}

	return sum;
	
}

float Photometry::sumIndices(const float *array1, float offset1, const float *array2, float offset2, const float *array3, float offset3, const std::vector<int> *vector, int idx){
	
	// use this for calculating the 2nd moments with different apertures: 
	// array1 is the image, offset1 is the background,
	// array2 and array3 are X or Y, their offsets are the 1st moments, 
	// then it is just an indices vector
	
	float sum=0;

	for(int i=0; i<vector[idx].size(); i++){
		
		float value1=array1[vector[idx][i]]-offset1; 
		float value2=array2[vector[idx][i]]-offset2; 
		float value3=array3[vector[idx][i]]-offset3; 
		
		if(pars.use_positives){
			if(value1>0) sum+=value1*value2*value3; // when use_positive is true, only sum positive pixels for 2nd moments
		}
		else{ 
			if(value1==value1) sum+=value1*value2*value3; // otherwise sum any non-NaN pixels
		}
		
	}

	return sum;
	
}
	
float Photometry::medianIndices(const float *array, const std::vector<int> *vector, int idx){ // find the median of the array points indicated by vector[idx]
	
	std::vector<float> values=std::vector<float>(); // an empty vector 
	
	for(int i=0;i<vector[idx].size(); i++) if(array[vector[idx][i]]==array[vector[idx][i]]) values.push_back(array[vector[idx][i]]); 
	
	if(values.size()==0) return 0; 
	
	std::sort(values.begin(), values.end()); 
	
	if(values.size()%2==1){ // odd number (pick middle value)
		
		return values[values.size()/2];
		
	}
	else{ // even number (take average of two middle values)
		
		return (values[values.size()/2-1] + values[values.size()/2])/2;
		
	}
}

// UTILITIES

void Photometry::printMatrix(const int *array, const char *name){ // debugging output of a cutout matrix
	
	mexPrintf("%s= \n", name);
	for(int j=0;j<dims[0];j++){ 
		for(int i=0;i<dims[1];i++) mexPrintf("%d ", array[i*dims[0]+j]);
		mexPrintf("\n");
	}
	
	mexPrintf("\n");
	
	
}

void Photometry::printMatrix(const float *array, const char *name){ // debugging output of a cutout matrix
	
	mexPrintf("%s= \n", name);
	for(int j=0;j<dims[0];j++){ 
		for(int i=0;i<dims[1];i++) mexPrintf("%5.3f ", array[i*dims[0]+j]);
		mexPrintf("\n");
	}
	
	mexPrintf("\n");
	
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