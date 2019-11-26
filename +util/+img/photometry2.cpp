#include "mex.h"
#include "matrix.h"
#include "photometry2.h"

Photometry photometry; 

const char Photometry::data_types[NUM_DATA_TYPES][STRLN]={"flux", "area", "error", "background", "variance", "offset_x", "offset_y", "width", "bad_pixels", "flag"};

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){

	// check inputs!
	if (nrhs==0){
		mexPrintf("Usage: phot_struct = photometry(cutouts, varargin)\n"); 
		mexPrintf("OPTIONAL ARGUMENTS:\n");
		mexPrintf("-------------------\n");
		mexPrintf("...\n");
		return;
	}
	
	// read the input data and parameters
	if(mxIsEmpty(prhs[0])){ // no input, then just return with all empty outputs...
		// const mwSize dims[]={0,0};
		// for(int i=0;i<9;i++) plhs[i]=mxCreateNumericArray(0,dims, mxSINGLE_CLASS, mxREAL); // return all empty arrays...
		return;
	}

	photometry.parseInputs(nrhs, prhs);
	photometry.makeArrays(); 
	
	photometry.run();

	char *field_names[6];
	for(int i=0;i<6;i++) field_names[i]=(char*) mxMalloc(STRLN);
	
	snprintf(field_names[0], STRLN, "raw_photometry"); 
	snprintf(field_names[1], STRLN, "forced_photometry"); 
	snprintf(field_names[2], STRLN, "apertures_photometry"); 
	snprintf(field_names[3], STRLN, "gaussian_photometry"); 
	snprintf(field_names[4], STRLN, "parameters"); 
	snprintf(field_names[5], STRLN, "arrays"); 
	
	plhs[0]=mxCreateStructMatrix(1,1,5, (const char**) field_names); 
	
	mxSetFieldByNumber(plhs[0], 0, 0, photometry.outputStruct(photometry.output_raw)); 
	mxSetFieldByNumber(plhs[0], 0, 1, photometry.outputStruct(photometry.output_forced)); 
	mxSetFieldByNumber(plhs[0], 0, 2, photometry.outputStruct(photometry.output_apertures, photometry.num_radii)); 
	mxSetFieldByNumber(plhs[0], 0, 3, photometry.outputStruct(photometry.output_gaussian)); 
	mxSetFieldByNumber(plhs[0], 0, 4, photometry.outputMetadataStruct());
	
	if(nlhs>1) plhs[1]=photometry.outputArraysStruct();
	
}

Photometry::Photometry(){ // class constructor

	// start by initializing any default arrays we may have
	ap_radii=new double[3];
	ap_radii[0]=3;
	ap_radii[1]=4;
	ap_radii[2]=5; 
	num_radii=3; 

}

Photometry::~Photometry(){ // destructor cleans up intermidiate arrays
	
	if(debug_bit) printf("now deleting photometry object!\n"); 
	
	if(ap_radii) delete(ap_radii); // make sure to free this little array
	
	deleteArrays();
	
}

mxArray *Photometry::outputStruct(float **output, int num_fluxes){ // wrap up the output matrices as a nice matlab style array

	// copy the field names into something persistent! 
	char *field_names[NUM_DATA_TYPES];
	for(int i=0;i<NUM_DATA_TYPES;i++){ 
		field_names[i]=(char*) mxMalloc(STRLN);
		snprintf(field_names[i], STRLN, data_types[i]);
	}
		
	mxArray *struct_array=mxCreateStructMatrix(1,1,NUM_DATA_TYPES, (const char**) field_names); 

	mwSize out_size[3]={(mwSize) output_size[0], (mwSize) output_size[1], (mwSize) num_fluxes}; 
	
	mxArray *matrix=0;
	for(int i=0;i<NUM_DATA_TYPES; i++){
		matrix=mxCreateNumericArray(3, out_size, mxSINGLE_CLASS, mxREAL);
		memcpy((float*)mxGetData(matrix), output[i], output_size[0]*output_size[1]*num_fluxes*sizeof(float));
		mxSetFieldByNumber(struct_array, 0, i, matrix);
	}
	
	return struct_array; 
	
}

mxArray *Photometry::outputMetadataStruct(){ // add a struct with some of the parameters and the different aperture masks used

	std::vector<std::string> field_names_vector;
	field_names_vector.push_back("cutout_size");
	field_names_vector.push_back("cutout_type"); 
	field_names_vector.push_back("aperture_radii"); 
	field_names_vector.push_back("annulus_radii");
	field_names_vector.push_back("forced_radius");
	field_names_vector.push_back("gauss_sigma");
	field_names_vector.push_back("gain");
	field_names_vector.push_back("shift_resolution");
	
	
	const char **field_names=new const char*[field_names_vector.size()];
	for(int i=0;i<field_names_vector.size(); i++) field_names[i]=field_names_vector[i].c_str();
	mxArray *struct_array=mxCreateStructMatrix(1,1,field_names_vector.size(), field_names); 
	
	delete(field_names);
	
	mxArray *array=0;
	array=mxCreateNumericMatrix(1,4,mxDOUBLE_CLASS, mxREAL);
	double *dbl_ptr=mxGetPr(array);
	for(int i=0;i<4;i++) dbl_ptr[i]=(double) dims[i];
	mxSetFieldByNumber(struct_array, 0, 0, array);
	
	const char *cutouts_class[1]={"single"};
	array=mxCreateCharMatrixFromStrings(1, cutouts_class);
	mxSetFieldByNumber(struct_array, 0, 1, array);
	
	array=mxCreateNumericMatrix(1, num_radii, mxDOUBLE_CLASS, mxREAL);
	dbl_ptr=mxGetPr(array);
	for(int i=0;i<num_radii;i++) dbl_ptr[i]=ap_radii[i];
	mxSetFieldByNumber(struct_array, 0, 2, array);
	
	array=mxCreateNumericMatrix(1, 2, mxDOUBLE_CLASS, mxREAL);
	dbl_ptr=mxGetPr(array);
	dbl_ptr[0]=inner_radius;
	dbl_ptr[1]=outer_radius;
	mxSetFieldByNumber(struct_array, 0, 3, array);
		
	array=mxCreateDoubleScalar(forced_radius);
	mxSetFieldByNumber(struct_array, 0, 4, array);
	
	array=mxCreateDoubleScalar(gauss_sigma);
	mxSetFieldByNumber(struct_array, 0, 5, array);
	
	array=mxCreateDoubleScalar(gain);
	mxSetFieldByNumber(struct_array, 0, 6, array);
	
	array=mxCreateDoubleScalar(resolution);
	mxSetFieldByNumber(struct_array, 0, 7, array);
	
	
	return struct_array; 

}

mxArray *Photometry::outputArraysStruct(){ // add a struct with the actual masks and grid arrays

	std::vector<std::string> field_names_vector;
	field_names_vector.push_back("aperture_masks");
	field_names_vector.push_back("annulus_masks");
	field_names_vector.push_back("gaussian_masks");
	field_names_vector.push_back("aperture_indices");
	field_names_vector.push_back("forced_indices");
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
	mwSize mask_dims[4]={dims[0],dims[1],num_shifts,num_radii};
	
	// save the aperture masks
	array=mxCreateNumericArray(4, mask_dims, mxSINGLE_CLASS, mxREAL);
	float *flt_ptr=(float*) mxGetData(array);
	memcpy(flt_ptr, apertures, N*num_shifts*num_radii*sizeof(float));
	mxSetFieldByNumber(struct_array, 0, 0, array);
	
	// save the annulus masks
	mask_dims[3]=1;
	array=mxCreateNumericArray(4, mask_dims, mxSINGLE_CLASS, mxREAL);
	flt_ptr=(float*) mxGetData(array);
	memcpy(flt_ptr, annulii, N*num_shifts*sizeof(float));
	mxSetFieldByNumber(struct_array, 0, 1, array);
	
	// save the gaussian masks
	array=mxCreateNumericArray(4, mask_dims, mxSINGLE_CLASS, mxREAL);
	flt_ptr=(float*) mxGetData(array);
	memcpy(flt_ptr, gaussians, N*num_shifts*sizeof(float));
	mxSetFieldByNumber(struct_array, 0, 2, array);
	
	// copy over the indices vectors
	array=outputIndicesVectors(aperture_indices, num_radii); 
	mxSetFieldByNumber(struct_array, 0, 3, array);
	
	array=outputIndicesVectors(forced_indices); 
	mxSetFieldByNumber(struct_array, 0, 4, array);
	
	array=outputIndicesVectors(annulus_indices); 
	mxSetFieldByNumber(struct_array, 0, 5, array);

	// copy out the X/Y and dx/dy arrays
	array=mxCreateNumericMatrix(dims[0],dims[1], mxSINGLE_CLASS, mxREAL); 
	flt_ptr=(float*) mxGetData(array); 
	memcpy(flt_ptr, X, N*sizeof(float));
	mxSetFieldByNumber(struct_array, 0, 6, array);
	
	array=mxCreateNumericMatrix(dims[0],dims[1], mxSINGLE_CLASS, mxREAL); 
	flt_ptr=(float*) mxGetData(array); 
	memcpy(flt_ptr, Y, N*sizeof(float));
	mxSetFieldByNumber(struct_array, 0, 7, array);
	
	array=mxCreateNumericMatrix(shift_dims[0],shift_dims[1], mxSINGLE_CLASS, mxREAL); 
	flt_ptr=(float*) mxGetData(array); 
	memcpy(flt_ptr, dx, num_shifts*sizeof(float));
	mxSetFieldByNumber(struct_array, 0, 8, array);	
	
	array=mxCreateNumericMatrix(shift_dims[0],shift_dims[1], mxSINGLE_CLASS, mxREAL); 
	flt_ptr=(float*) mxGetData(array); 
	memcpy(flt_ptr, dy, num_shifts*sizeof(float));
	mxSetFieldByNumber(struct_array, 0, 9, array);
	
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
		for(int j=0;j<ndims; j++) dims[j]=new_dims[j];
		N=(int) (dims[0]*dims[1]); // dims[0] is the height while dims[1] is the width
	}
	
	mwSize new_output_size[2]={0};
	// size of non-empty output is [size(cutouts,3), size(cutouts,4)]
	if(ndims>2) new_output_size[0]=dims[2];
	else new_output_size[0]=1;
	if(ndims>3) new_output_size[1]=dims[3];
	else new_output_size[1]=1;
	
	if(new_output_size[0]!=output_size[0] || new_output_size[1]!=output_size[1]){
		
		deleteArrays(); 
		output_size[0]=new_output_size[0];
		output_size[1]=new_output_size[1];
		
	}
	
	num_cutouts=(int) (output_size[0]*output_size[1]); // how many values in each of the above arrays... 
	
	for(int i=1;i<nrhs;i+=2){ // parse varargin
		
		char key[STRLN];
		if(mxIsChar(prhs[i])==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotChar", "Input %d to photometry is not a string...", i+1);
		mxGetString(prhs[i], key, STRLN); // copy the string data
		
		mxArray *val=0;
		if(i+1<nrhs) val=(mxArray*) prhs[i+1]; // if the varargin is odd numbered, leave val=0 as default
		
		if(cs(key, "aperture", "radius")){

			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsNumeric(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not numeric...", i+2);
			
			// temporary values
			double *new_ap_radii=(double*) mxGetData(val);
			int new_num_radii=(int) mxGetNumberOfElements(val);
			
			if(new_num_radii!=num_radii){
				deleteOutputArray(output_apertures); 
				deleteApertureMasks();
				num_radii=new_num_radii;
				if(ap_radii) delete(ap_radii); 
				ap_radii=new double[num_radii]; 
				for(int j=0;j<num_radii; j++) ap_radii[j]=new_ap_radii[j]; 
			}
			else{
				for(int j=0;j<num_radii; j++){
					if(new_ap_radii[j]!=ap_radii[j]) deleteApertureMasks(); 
					ap_radii[j]=new_ap_radii[j];
				}
			}
			
			forced_radius=ap_radii[num_radii-1]; // by default we use the biggest wedding cake radius as the forced photometry radius!

		}
		else if(cs(key, "gaussian", "sigma", "gauss_sigma")){
			
			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not a numeric scalar...", i+2);
			double new_gauss_sigma=mxGetScalar(val);
			
			if(new_gauss_sigma!=gauss_sigma){
				
				deleteGaussianMasks();
				gauss_sigma=new_gauss_sigma;
				
			}
			
		}
		else if(cs(key, "annulus")){

			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsNumeric(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not numeric...", i+2);
			
			int num=(int) mxGetNumberOfElements(val);
			double *new_annuli=(double*) mxGetData(val);
			
			if(num>0 && inner_radius!=new_annuli[0]){
				deleteAnnulusMasks();
				inner_radius=new_annuli[0];
			}
			
			if(num>1 && outer_radius!=new_annuli[1]){
				
				outer_radius=new_annuli[1];
			} 
			else if(outer_radius!=0){ 
				deleteAnnulusMasks();
				outer_radius=0;
			}
			else outer_radius=0; 
			
		}
		else if(cs("key", "gain")){
			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsEmpty(val)) continue;
			if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericSingle", "Input %d to photometry is not a numeric scalar!", i+2);
			gain=mxGetScalar(val);
			
		}
		else if(cs(key, "resolution")){
			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsEmpty(val)) continue;
			if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericSingle", "Input %d to photometry is not a numeric scalar!", i+2);
			int new_resolution=(int) mxGetScalar(val);
			
			if(new_resolution!=resolution){
				deleteGrids(); 
				deleteMasks(); 
				resolution=new_resolution;
			}
			
		}
		else if(cs(key, "debug_bit")){
			if(val==0 || mxIsEmpty(val)) debug_bit=1; // if no input, assume positive
			else{
				if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry is not a numeric scalar...", i+2);
				debug_bit=(int) mxGetScalar(val);
			}
			
		}
		else if(cs(key, "threads")){
			if(val==0 || mxIsEmpty(val)) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			else{
				if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry is not a numeric scalar...", i+2);
				num_threads=(int) mxGetScalar(val);
			}
			
		}
		
	}
	
	
	if(debug_bit>1){ // check that all inputs have been received! 
		mexPrintf("cutouts: [");
		for(int i=0;i<ndims;i++){if(i>0) mexPrintf("x"); mexPrintf("%d", dims[i]); }
		mexPrintf("] (%s) \n", mxGetClassName(prhs[0]));
		// additional printbacks come here... 
	}
	
}

// INITIALIZATIONS

void Photometry::clear(){ // make sure all the output arrays are NaN before filling them

	// to be added later...

}

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
	
	for(int j=0;j<2;j++) shift_dims[j]=dims[j]*resolution;
	num_shifts=(int) (shift_dims[0]*shift_dims[1]); 
	
	dx=new float[num_shifts];
	dy=new float[num_shifts];
	
	meshgrid(dx, dy, shift_dims, resolution);
	
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
	
	if(output_raw) initializeOutputArray(output_raw);
	else allocateOutputArray(output_raw);
	
	if(output_forced) initializeOutputArray(output_forced);
	else allocateOutputArray(output_forced);
	
	if(output_apertures) initializeOutputArray(output_apertures, num_radii);
	else allocateOutputArray(output_apertures, num_radii);
	
	if(output_gaussian) initializeOutputArray(output_gaussian);
	else allocateOutputArray(output_gaussian); 
		
}

void Photometry::allocateOutputArray(float **&output, int num_fluxes){ // create new memory for flux, background etc...

	output=new float*[NUM_DATA_TYPES];
	for(int i=0;i<NUM_DATA_TYPES;i++) output[i]=new float[output_size[0]*output_size[1]*num_fluxes]; 
	
}

void Photometry::initializeOutputArray(float **output, int num_fluxes){ // clear the content of existing memory for flux, background etc...

	

}

void Photometry::makeMasks(){ // make a bank of masks to be used for different shifts and types of photometry

	if(apertures==0) makeApertureMasks(); 
	if(annulii==0) makeAnnulusMasks();
	if(gaussians==0) makeGaussianMasks(); 

}

void Photometry::makeApertureMasks(){ // make the aperture masks for e.g., forced photometry, with all possible shifts (dx/dy)
	
	apertures=new float[N*num_shifts*num_radii]; 
	aperture_indices=new std::vector<int>[num_shifts*num_radii]; // keep a running list of all the aperture indices for quick summing
	forced_indices=new std::vector<int>[num_shifts]; // same only for the last radii we keep the entire circle (for forced photometry)
	
	for(int k=0;k<num_radii;k++){ // go over all aperture radii
	
		double radius=ap_radii[k]; 
	
		for(int j=0; j<num_shifts;j++){ // go over the different shifts
		
			for(int i=0;i<N;i++){ // go over the pixels in each mask
				
				float r=(float) sqrt(pow(X[i]-dx[j],2)+pow(Y[i]-dy[j],2));
				
				float value=(float) (radius+0.5-r);
				if(value<1) value=0;
				else if(value>1) value=1;
				apertures[k*num_shifts*N + j*N + i] = value;
				
				if(k==num_radii-1){ // do the forced aperture only on the last k*num_shifts
					if(value>0) forced_indices[j].push_back(i);
				}
				
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
	
	for(int j=0; j<num_shifts;j++){ // go over the different shifts
	
		for(int i=0;i<N;i++){ // go over the pixels in each mask
			
			float r=(float) sqrt(pow(X[i]-dx[j],2)+pow(Y[i]-dy[j],2));
			
			float value=0;
			if(r>=inner_radius && r<outer_radius) value=1;
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
			
			gaussians[j*N + i] = (float) exp(-0.5*pow(r/gauss_sigma,2));
			
		} // for i
		
	} // for j 
	
}

// CLEANUP and DELETION

void Photometry::deleteArrays(){ // get rid of all the masks and output arrays (when destroying this object)

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

	deleteOutputArray(output_raw);
	deleteOutputArray(output_forced);
	deleteOutputArray(output_apertures);
	deleteOutputArray(output_gaussian);

}

void Photometry::deleteOutputArray(float **output){ // go over and deallocate the memory for one of the outputs

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

	if(num_threads<=1){
		for(int j=0;j<num_cutouts;j++){ // number of cutouts
			
			calculate(j);
			
		} // for j
	}
	else{
		
		int step=num_cutouts/num_threads;
		int current_idx=0;
		
		std::vector<std::thread> t;
				
		for(int i=0;i<num_threads-1;i++){
			
			if(debug_bit>2) mexPrintf("Sending a thread for indices %d to %d\n", current_idx, current_idx+step);
			t.push_back(std::thread(&Photometry::run_idx, this, current_idx, current_idx+step));
			current_idx+=step; 
			
		}
		
		if(debug_bit>2) mexPrintf("Running on main thread indices %d to %d\n", current_idx, num_cutouts);
		run_idx(current_idx,num_cutouts); // run the remaining cutouts on the main thread! 
		
		for(int i=0;i<num_threads-1;i++){
			t[i].join();
		}
		
	}
	
	if(debug_bit>2) printMatrix(cutouts, "cutouts");
	
}

void Photometry::run_idx(int start_idx, int end_idx){ // run on subsets of cutouts, possibly in separate threads...

	for(int j=start_idx;j<end_idx;j++){ // partial list of cutouts
			
			calculate(j);
			
	} // for j

}

void Photometry::calculate(int j){ // do the actual calculations on a single cutout
	
	float *image=&cutouts[j*N]; // a pointer to the raw data
	
	float **output=output_raw; // make sure to get pointers to the raw output data	
	float *flux=output[0];
	float *error=output[1];
	float *area=output[2];
	float *background=output[3];
	float *variance=output[4];
	float *offset_x=output[5];
	float *offset_y=output[6];
	float *width=output[7];
	float *bad_pixels=output[8];
	float *flag=output[9];
	
	bad_pixels[j]=countNaNs(image);
	
	if(bad_pixels[j]==N){ // image is all NaN!
		
		// mexPrintf("skipping image, all nans\n");
		flux[j]=NAN;
		area[j]=NAN;
		error[j]=NAN;
		background[j]=NAN;
		variance[j]=NAN;
		offset_x[j]=NAN;
		offset_y[j]=NAN;
		width[j]=NAN;
		flag[j]=1;
		
		for(int i=0;i<9;i++) output_gaussian[i][j]=NAN;
		output_gaussian[9][j]=1;
		
		for(int i=0;i<9;i++) output_forced[i][j]=NAN;
		output_gaussian[9][j]=1;
		
		for(int k=0;k<num_radii;k++) for(int i=0;i<9;i++) output_apertures[i][k*num_cutouts + j]=NAN;
		output_gaussian[9][j]=1;
		
		return;
		
	}
	
	flux[j]=sumArrays(image); // simple raw photometry! 
	area[j]=N-bad_pixels[j]; 
	
	// error[j]
	// background[j]
	// variance[j]
	
	// first moments:
	float m1x=sumArrays(image, X);
	float m1y=sumArrays(image, Y);

	offset_x[j]=m1x;
	offset_y[j]=m1y;
	
	// second moments
	float m2x=sumArrays(image, X, m1x, X, m1x); // I*(X-m1x)^2
	float m2y=sumArrays(image, Y, m1y, Y, m1y); // I*(Y-m1y)^2
	float mxy=sumArrays(image, X, m1x, Y, m1y); // I*(X-m1x)*(Y-m1y)
	
	width[j]=getWidthFromMoments(m2x, m2y, mxy); 
	
	flag[j]=checkMoments(offset_x[j], offset_y[j], width[j]); 
	
	// start doing the gaussian photometry! 
	
	
}

void Photometry::runForced(){
	
	
}

void Photometry::runForced_idx(int start_idx, int end_idx){
	
	
}

void Photometry::calculateForced(int j){
	
	
}

bool Photometry::checkMoments(float offset_x, float offset_y, float width){ // returns 1 if there is a problem with the offsets or width
	
	if(offset_x>dims[1]/2) return 1;
	if(offset_y>dims[0]/2) return 1;
	if(width>dims[0]/2 || width>dims[1]/2) return 1;
	if(width<0) return 1;
	if(isnan(offset_x) || isnan(offset_y) || isnan(width)) return 1;
	
	return 0; // if all checks are not triggered, return with the ok flag
	
	
}

// GET WIDTHS AND OFFSETS

float Photometry::getWidthFromMoments(float m2x, float m2y, float mxy){ // calculate the eigenvalues of the 2nd moments and from that find the average width
// got this little nugget from: https://yutsumura.com/express-the-eigenvalues-of-a-2-by-2-matrix-in-terms-of-the-trace-and-determinant/

	float tr=m2x+m2y;
	float det=m2x*m2y - mxy*mxy;
	
	float r1=(tr-sqrt(tr*tr-4*det))/2;
	float r2=(tr+sqrt(tr*tr-4*det))/2;
	
	return (sqrt(r1)+sqrt(r2))/2;

}

float Photometry::getAverageWidth(float **output){ // calculate the average width from all the good measurements in this output type

	// need to add a check for flag too! 

	float *width=output[7];
	float *flag=output[9];

	float S=0;
	int counter=0;
	for(int j=0;j<num_cutouts;j++) if(flag[j]==0) { S+=width[j]; counter++; }
	return S/counter; // mean of the widths, excluding NaNs

}

float Photometry::getAverageOffsetX(float **output){ // to be implemented!

	return 0;

}

float Photometry::getAverageOffsetY(float **output){ // to be implemented!

	return 0;

}

// ARRAY ARITHMETIC

bool Photometry::countNaNs(const float *array){ // check if an array has all NaN values
	
	int num_nans=0;
	
	for(int i=0;i<N; i++) if(array[i]!=array[i]) { num_nans++; }
		
	return num_nans;
	
}

float Photometry::sumArrays(const float *array1){ // just the sum of all non-NaN pixels in this array
	
	float S=0;
	
	for(int i=0;i<N;i++) if(array1[i]!=array1[i]) S+=array1[i];
	
	return S;
	
}

float Photometry::sumArrays(const float *array1, const float *array2){ // sum of all non-NaN values of array1*array2
	
	float S=0;
	
	for(int i=0;i<N;i++) 
		// if(_isnanf(array1[i])==0 && _isnanf(array2[i])==0) 
		//if(_isnanf(array1[i])==0) 
		if(array1[i]!=array1[i])
			S+=array1[i]*array2[i];
	
	return S;
	
}

float Photometry::sumArrays(const float *array1, float offset1, const float *array2){ // sum of all non-NaN values of array1*array2
	
	float S=0;
	
	for(int i=0;i<N;i++) if(array1[i]!=array1[i]) S+=(array1[i]-offset1)*array2[i];
	
	return S;
	
}

float Photometry::sumArrays(const float *array1, const float *array2, const float *array3){ // sum of all non-NaN values of array1*array2*array3
	
	float S=0;
	
	for(int i=0;i<N;i++) 
		if(array1[i]!=array1[i])
			S+=array1[i]*array2[i]*array3[i];
	
	return S;
	
}

float Photometry::sumArrays(const float *array1, const float *array2, float offset2, const float *array3, float offset3){ // sum of all non-NaN values of array1*array2*array3
	
	float S=0;
	
	for(int i=0;i<N;i++) 
		//if(_isnanf(array1[i])==0 && _isnanf(array2[i])==0 && _isnanf(array3[i])==0) 
		//if(_isnanf(array1[i]))
		if(array1[i]!=array1[i])
			S+=array1[i]*(array2[i]-offset2)*(array3[i]-offset3);
	
	return S;
	
}

float Photometry::sumArrays(const float *array1, const float *array2, const float *array3, const float *array4){ // sum of all non-NaN values of array1*array2*array3*array4
	
	float S=0;
	
	for(int i=0;i<N;i++) 
		// if(_isnanf(array1[i])==0 && _isnanf(array2[i])==0 && _isnanf(array3[i])==0 && _isnanf(array4[i])==0) 
		if(array1[i]!=array1[i])
			S+=array1[i]*array2[i]*array3[i]*array4[i];
	
	return S;
	
}

float Photometry::sumArrays(const float *array1, float offset1, const float *array2, float offset2, const float *array3, float offset3, const float *array4){ // sum of all non-NaN values of array1*array2*array3*array4
	
	float S=0;
	
	for(int i=0;i<N;i++) 
		// if(_isnanf(array1[i])==0 && _isnanf(array2[i])==0 && _isnanf(array3[i])==0 && _isnanf(array4[i])==0) 
		if(array1[i]!=array1[i])
			S+=(array1[i]-offset1)*(array2[i]-offset2)*(array3[i]-offset3)*array4[i];
	
	return S;
	
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