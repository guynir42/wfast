#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <cctype>
#include <thread>
#include <chrono>

#define STRLN 64 // maximum string length (for copying)
#define NUM_DATA_TYPES 10 // flux, error, area, background, variance, offset_x, offset_y, width, bad_pixels, flag

// utility functions to compare strings
bool cs(const char *keyword, const char *compare_str, int num_letters=3);
bool cs(const char *keyword, const char *str1, const char *str2, int num_letters=3);
bool cs(const char *keyword, const char *str1, const char *str2, const char *str3, int num_letters=3);
bool parse_bool(mxArray *value);

// Usage: [flux, error, area, background, variance, offset_x, offset_y, width] = photometry(cutouts, varargin)

class Photometry{
	
	public:
	
	float *cutouts=0; // the input matrix on which we calculate everything (dim 1&2 are y&x, dim 3 is frame number, dim 4 is star number)
	mwSize *dims=0; // array with the sizes of cutouts
	mwSize ndims=0; // number of dimensions of cutouts
	int num_cutouts=0;  // number of cutouts (dim 3 times dim 4)
	int N=0; // number of pixels in each cutout (dim 1 times dim 2)
	
	double gain_scalar=1; // calculate the source noise using the gain
	int num_threads=1; // for future use with built-in multithreading
	
	int debug_bit=0;
	
	// these can be shared with all calculations (i.e., they are kept in memory between calls to photometry2!)
	float *X=0; // grid points not including shifts. 
	float *Y=0; // grid points not including shifts. 
	
	int num_shifts=0; // number of shifts we need to cover all possible star positions (i.e., the length of dx and dy)
	float *dx=0; // list of offsets in x for the center of the mask
	float *dy=0; // list of offsets in y for the center of the mask
		
	float *apertures=0; // 3D matrix of apertures used for forced photometry, one for each dx/dy shift
	double radius=5; // the default radius used for forced photometry (pixels)

	float *wedding_cakes=0; // 4D matrix of aperture+annuli for the wedding cake photometry
	int num_radii=0; // how many different aperture arrays do we have for the wedding cake
	double *ap_radii=0; // radii of different apertures, given in pixel units
	
	float *annulus=0; // 3D matrix of the annulus used for all kinds of photometry, one for each dx/dy shift
	double inner_radius=0; // of the annulus (pixels)
	double outer_radius=0; // of the annulus (pixels)

	float *gaussians=0; // 3D matrix of gaussian weighted apertures for "PSF" photometry, one for each dx/dy shift
	double gauss_sigma=2; // the width of the gaussian (in pixels)
	
	// output arrays are defined here (in C++)
	mwSize output_size[2]={0,0}; // the same as dim 3 and dim 4 of cutouts
	float **output_raw=0;
	float **output_forced=0;
	float **output_wedding_cake=0;
	float **output_gaussian=0;
	const static char data_types[NUM_DATA_TYPES][STRLN]; // this holds strings containing: flux, error, area, background, variance, offset_x, offset_y, width, bad_pixels, flag
	
	// function prototypes (implementation at the end)
	Photometry();
	~Photometry();
	void clear(); // make sure all the output arrays are NaN before filling them
	void parseInputs(int nrhs, const mxArray *prhs[]);
	
	void makeArrays(); // create all the required memory if it isn't already allocated, make all the required masks
	void allocateOutputArray(float **output, int num_fluxes); // allocate memory for a set of outputs: flux, error, area, background, variance, offset_x/y, width, bad_pixels, flag
	void initializeOutputArray(float **output, int num_fluxes); // re-use memory but make sure it is intialized to NaN first
	void makeMasks(); // make a bank of masks to be used for different shifts and types of photometry
	void deleteArrays(); // get rid of all the masks and output arrays (when destroying this object)
	void deleteOutputArray(float **output); // go over and deallocate the memory for all the outputs
	void deleteMasks(); // go over and deallocate the memory for all the masks
	void circle_mask(float *array, float dx, float dy, double radius);
	void annulus_mask(float *array, float dx, float dy, double radius1, double radius2);
	void gaussian_mask(float *array, float dx, float dy);
	
	void run(float **output);
	void run_idx(int start_idx, int end_idx, float **output); // use this to run only a subset of the cutouts (for multithreading)
	void calculate(int j, float **output); // do the actual work on some cutout array
	float getWidthFromMoments(float m2x, float m2y, float mxy); // from the eigenvalues of the 2nd moments
	
	// make averages over the results
	float getAverageWidth(float **output);
	float getAverageOffsetX(float **output);	
	float getAverageOffsetY(float **output);
	
	// the sum of the product of array1...
	float sumArrays(const float *array1);
	float sumArrays(const float *array1, const float *array2);
	float sumArrays(const float *array1, const float *array2, const float *array3);
	float sumArrays(const float *array1, const float *array2, const float *array3, const float *array4);
	
	// just to catch some bugs:
	bool isAllNaNs(const float *array); 
	
	// print on screen
	void printMatrix(const int *array, const char *name);
	void printMatrix(const float *array, const char *name);

} photometry; 

const char Photometry::data_types[NUM_DATA_TYPES][STRLN]={"flux", "error", "area", "background", "variance", "offset_x", "offset_y", "width", "bad_pixels", "flag"};

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
	// printf("flux(1,1)= %f\n", photometry.flux[0]); 
	// photometry.flux[0]++;
		
	// Photometry phot(nrhs, prhs);
	
	// phot.run();

	
}

Photometry::Photometry(){ // class constructor

	// allocate intermediate arrays
	// X=(float *) mxCalloc(N, sizeof(float));
	// Y=(float *) mxCalloc(N, sizeof(float)); 
	
	// for(int i=0;i<N;i++){ // meshgrid
		
		// X[i]=(float)(i/dims[0])-dims[1]/2;
		// Y[i]=(float)(i%dims[0])-dims[0]/2;
		
	// }
	
	// if(debug_bit>2){ // check that meshgrid returned what we expect
		// printMatrix(X, "X");
		// printMatrix(Y, "Y");
	// }
	
}

Photometry::~Photometry(){ // destructor cleans up intermidiate arrays
	
	// mxFree(X);
	// mxFree(Y);
	
	
}

void Photometry::parseInputs(int nrhs, const mxArray *prhs[]){ // take the cutouts input and the matlab style varargin and parse into c++

	if(mxIsClass(prhs[0], "single")==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotFloat", "Input 1 to photometry is not a single/float array...");
	cutouts=(float*) mxGetData(prhs[0]);
	dims=(mwSize*)mxGetDimensions(prhs[0]);
	ndims=mxGetNumberOfDimensions(prhs[0]);	
	N=dims[0]*dims[1]; // dims[0] is the height while dims[1] is the width
	
	// size of non-empty output is [size(cutouts,3), size(cutouts,4)]
	if(ndims>2) output_size[0]=dims[2];
	else output_size[0]=1;
	if(ndims>3) output_size[1]=dims[3];
	else output_size[1]=1;
	num_cutouts=output_size[0]*output_size[1]; // how many values in each of the above arrays... 
	
	for(int i=1;i<nrhs;i+=2){ // parse varargin
		
		char key[STRLN];
		if(mxIsChar(prhs[i])==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotChar", "Input %d to photometry is not a string...", i+1);
		mxGetString(prhs[i], key, STRLN); // copy the string data
		
		mxArray *val=0;
		if(i+1<nrhs) val=(mxArray*) prhs[i+1]; // if the varargin is odd numbered, leave val=0 as default
		
		if(cs(key, "aperture", "radius")){

			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsNumeric(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not numeric...", i+2);
			ap_radii=(double*) mxGetData(val);
			num_radii=mxGetNumberOfElements(val);
			radius=ap_radii[num_radii-1]; // by default we use the biggest wedding cake radius as the forced photometry radius!
			
		}
		else if(cs(key, "gaussian", "sigma", "gauss_sigma")){
			
			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not a numeric scalar...", i+2);
			gauss_sigma=mxGetScalar(val);
			
		}
		else if(cs(key, "annulus")){

			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsNumeric(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not numeric...", i+2);
			
			int num=mxGetNumberOfElements(val);
			double *annuli=(double*) mxGetData(val);
			if(num>0) inner_radius=annuli[0];
			if(num>1) outer_radius=annuli[1];
			else outer_radius=0;
			
		}
		else if(cs("key", "gain_scalar")){
			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsEmpty(val)) continue;
			if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericSingle", "Input %d to photometry is not a numeric scalar!", i+2);
			gain_scalar=mxGetScalar(val);
			
		}
		else if(cs(key, "debug_bit")){
			if(val==0 || mxIsEmpty(val)) debug_bit=1; // if no input, assume positive
			else{
				if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry is not a numeric scalar...", i+2);
				debug_bit=mxGetScalar(val);
			}
			
		}
		else if(cs(key, "threads")){
			if(val==0 || mxIsEmpty(val)) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			else{
				if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry is not a numeric scalar...", i+2);
				num_threads=mxGetScalar(val);
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

void Photometry::makeArrays(){ // generate all the one-time arrays in memory
	
	if(output_raw) initializeOutputArray(output_raw);
	else allocateOutputArray(output_raw);
	
	if(output_forced) initializeOutputArray(output_forced);
	else allocateOutputArray(output_forced);
	
	if(output_wedding_cake) initializeOutputArray(output_wedding_cake);
	else allocateOutputArray(output_wedding_cake);
	
	if(output_gaussian) initializeOutputArray(output_gaussian);
	else allocateOutputArray(output_gaussian); 
		
}

void Photometry::initializeOutputArray(float **output){ // clear the content of existing memory for flux, background etc...

	

}

void Photometry::allocateOutputArray(float **output){ // create new memory for flux, background etc...

}

void Photometry::run(float **output){ // run photometry on all cutouts! the results are put into various sub-arrays of "output"

	// mexPrintf("num_cutouts= %d\n", num_cutouts);

	if(num_threads<=1){
		for(int j=0;j<num_cutouts;j++){ // number of cutouts
			
			calculate(j, output);
			
		} // for j
	}
	else{
		
		int step=num_cutouts/num_threads;
		int current_idx=0;
		
		std::vector<std::thread> t;
				
		for(int i=0;i<num_threads-1;i++){
			
			if(debug_bit>2) mexPrintf("Sending a thread for indices %d to %d\n", current_idx, current_idx+step);
			t.push_back(std::thread(&Photometry::run_idx, this, current_idx, current_idx+step, output));
			current_idx+=step; 
			
		}
		
		if(debug_bit>2) mexPrintf("Running on main thread indices %d to %d\n", current_idx, num_cutouts);
		run_idx(current_idx,num_cutouts, output); // run the remaining cutouts on the main thread! 
		
		for(int i=0;i<num_threads-1;i++){
			t[i].join();
		}
		
	}
	
	if(debug_bit>2) printMatrix(cutouts, "cutouts");
	
}

void Photometry::run_idx(int start_idx, int end_idx, float **output){ // run on subsets of cutouts, possibly in separate threads...

	for(int j=start_idx;j<end_idx;j++){ // partial list of cutouts
			
			calculate(j, output);
			
	} // for j

}

void Photometry::calculate(int j, float **output){ // do the actual calculations on a single cutout
	
	// float *x=(float *)mxCalloc(N, sizeof(float)); // grid with shift added using offset_x
	// float *y=(float *)mxCalloc(N, sizeof(float)); // grid with shift added using offset_y
	// float *x2=(float *)mxCalloc(N, sizeof(float)); // grid with shift added using found moment
	// float *y2=(float *)mxCalloc(N, sizeof(float)); // grid with shift added using found moment
	// float *ap_array=(float *)mxCalloc(N, sizeof(float)); // grid with the aperture shape
	// float *bg_array=(float *)mxCalloc(N, sizeof(float)); // grid with the background shape
	// float *weight=(float *)mxCalloc(N, sizeof(float)); // weight includes PSF, aperture, variance map
	// float *error_array=(float *)mxCalloc(N, sizeof(float)); // the denominator for the expression of weight
	// float *bad_array=(float *)mxCalloc(N, sizeof(float)); // grid with the bad pixels marked
	// float *image_sub=(float *)mxCalloc(N, sizeof(float));
	// float *image=(float *)mxCalloc(N, sizeof(float)); // can be a copy of image_raw or image_sub depending on value of "subtract" option
	float *image_raw=&cutouts[j*N]; // a pointer to the raw data
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
	
	if(isAllNaNs(image_raw)){
		// mexPrintf("skipping image, all nans\n");
		flux[j]=NAN;
		error[j]=NAN;
		area[j]=0;
		background[j]=NAN;
		variance[j]=NAN;
		offset_x[j]=NAN;
		offset_y[j]=NAN;
		width[j]=NAN;
		bad_pixels[j]=N;
		flag[j]=1;
		return;
	}
	
	flux[j]=sumArrays(image_raw); 
	
	// mxFree(x);
	// mxFree(y);
	// mxFree(x2);
	// mxFree(y2);
	// mxFree(ap_array);
	// mxFree(bg_array);
	// mxFree(weight);
	// mxFree(error_array);
	// mxFree(bad_array);
	// mxFree(image);
	// mxFree(image_sub);
	
}

float Photometry::getWidthFromMoments(float m2x, float m2y, float mxy){ // calculate the eigenvalues of the 2nd moments and from that find the average width
// got this little nugget from: https://yutsumura.com/express-the-eigenvalues-of-a-2-by-2-matrix-in-terms-of-the-trace-and-determinant/

	float tr=m2x+m2y;
	float det=m2x*m2y - mxy*mxy;
	
	float r1=(tr-sqrt(tr*tr-4*det))/2;
	float r2=(tr+sqrt(tr*tr-4*det))/2;
	
	return (sqrt(r1)+sqrt(r2))/2;

}

void Photometry::circle_mask(float *array, float dx, float dy, double radius){ // make the aperture mask for e.g., forced photometry, with shift dx/dy
	
	for(int i=0;i<N;i++){
		
		float r=(float) sqrt(pow(X[i]-dx,2)+pow(Y[i]-dy,2));
		
		array[i]=radius+0.5-r;
		if(array[i]<1) array[i]=0;
		else if(array[i]>1) array[i]=1;
		
	}
}

void Photometry::annulus_mask(float *array, float dx, float dy, double radius1, double radius2){ // make a mask for the background annulus or wedding cake photometry
	
	if(radius2<=0) radius2=1e10; // replace for Inf
	
	for(int i=0;i<N;i++){
		
		float r=(float) sqrt(pow(X[i]-dx,2)+pow(Y[i]-dy,2));
		
		if (r>radius1 && r<radius2)	array[i]=1;
		else array[i]=0; 
		
	}
	
}

void Photometry::gaussian_mask(float *array, float dx, float dy){ // make a gaussian weighted mask with shift dx/dy
	
	for(int i=0;i<N;i++){
		
		float r=(float) sqrt(pow(X[i]-dx,2)+pow(Y[i]-dy,2));
		
		array[i]=exp(-0.5*pow(r/gauss_sigma,2));
		
	}
	
	
}

float Photometry::getAverageWidth(float **output){ // calculate the average width from all the good measurements in this output type

	// need to add a check for flag too! 

	float *width=output[7];

	float S=0;
	int counter=0;
	for(int j=0;j<num_cutouts;j++) if(isnan(width[j])==0) { S+=width[j]; counter++; }
	return S/counter; // mean of the widths, excluding NaNs

}

float Photometry::getAverageOffsetX(float **output){ // to be implemented!

	return 0;

}

float Photometry::getAverageOffsetY(float **output){ // to be implemented!

	return 0;

}

float Photometry::sumArrays(const float *array1){ // just the sum of all non-NaN pixels in this array
	
	float S=0;
	
	for(int i=0;i<N;i++) if(isnan(array1[i])==0) S+=array1[i];
	
	return S;
	
}

float Photometry::sumArrays(const float *array1, const float *array2){ // sum of all non-NaN values of array1*array2
	
	float S=0;
	
	for(int i=0;i<N;i++) if(isnan(array1[i])==0 && isnan(array2[i])==0) S+=array1[i]*array2[i];
	
	return S;
	
}

float Photometry::sumArrays(const float *array1, const float *array2, const float *array3){ // sum of all non-NaN values of array1*array2*array3
	
	float S=0;
	
	for(int i=0;i<N;i++) 
		if(isnan(array1[i])==0 && isnan(array2[i])==0 && isnan(array3[i])==0) 
			S+=array1[i]*array2[i]*array3[i];
	
	return S;
	
}

float Photometry::sumArrays(const float *array1, const float *array2, const float *array3, const float *array4){ // sum of all non-NaN values of array1*array2*array3*array4
	
	float S=0;
	
	for(int i=0;i<N;i++) 
		if(isnan(array1[i])==0 && isnan(array2[i])==0 && isnan(array3[i])==0 && isnan(array4[i])==0) 
			S+=array1[i]*array2[i]*array3[i]*array4[i];
	
	return S;
	
}

bool Photometry::isAllNaNs(const float *array){ // check if an array has all NaN values
	
	int num_vals=0;
	int num_nans=0;
	
	for(int i=0;i<N; i++){
		
		num_vals++;
		if (isnan(array[i])) num_nans++;
		
	}
	
	return num_vals==num_nans;
	
}

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
	
	return success;
	
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