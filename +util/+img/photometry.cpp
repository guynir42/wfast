#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <cctype>

#define STRLN 64 // maximum string length (for copying)

// utility functions to compare strings
bool cs(const char *keyword, const char *compare_str, int num_letters=3);
bool cs(const char *keyword, const char *str1, const char *str2, int num_letters=3);
bool cs(const char *keyword, const char *str1, const char *str2, const char *str3, int num_letters=3);
bool parse_bool(mxArray *value);

// Usage: [flux, weight, background, variance, offset_x, offset_y, width] = photometry(cutouts, varargin)

class Photometry{
	
	public:
	
	float *cutouts=0;
	
	float epsilon=0.1f; // minimal value for both dx, dy to change. If change is smaller than epsilon, further iterations are skipped... 
	mwSize *dims=0;
	mwSize ndims=0;
	int N=0; 
	int num_iter=1;
	int subtract=1;
	int debug_bit=0;
	mwSize out_dims[2]={0,0};
	int num_cutouts=0;
	float multiplier=3.5; // default radius of aperture is previous width time this multiplier
	float seeing=2; // default gauss width is seeing*2.355

	// definitions of aperture/annulus types
	enum ap_enum {SQUARE=0, CIRCLE, GAUSSIAN};
	enum bg_enum {CORNERS=0, ANNULUS};

	char aperture_string[STRLN] = "square";
	ap_enum ap_type=SQUARE;
	double *ap_pars=0;
	int num_ap_pars=0;
	
	char background_string[STRLN] = "corners";
	bg_enum bg_type=CORNERS;
	double *bg_pars=0;
	int num_bg_pars=0;

	// function pointers
	// void (Photometry::*ap_func)(float *array, float *x, float *y);
	// void (Photometry::*bg_func)(float *array, float *x, float *y);

	// these can be shared with all calculations
	float *X=0; // grid points not including shifts. 
	float *Y=0; // grid points not including shifts. 

	// output arrays are defined here (in C++)
	float *flux=0;
	float *weight=0;
	float *background=0;
	float *variance=0;
	float *offset_x=0;
	float *offset_y=0;
	float *width=0;
	float *bad_pixel=0;

	// output arrays are defined here (in matlab pointers)
	mxArray *flux_ptr=0;
	mxArray *weight_ptr=0;
	mxArray *background_ptr=0;
	mxArray *variance_ptr=0;
	mxArray *offset_x_ptr=0;
	mxArray *offset_y_ptr=0;
	mxArray *width_ptr=0;
	mxArray *bad_pixel_ptr=0;
	
	// function prototypes (implementation at the end)
	Photometry(int nrhs, const mxArray *prhs[]);
	~Photometry();
	void parseInputs(int nrhs, const mxArray *prhs[]);
	void run();
	void calculate(int j);
	float getWidthFromMoments(float m2x, float m2y, float mxy);
	
	// these are various shapes that can be used to intersect the image/cutout
	void main_mask(float *array, float *x, float *y);
	void square_mask(float *array, float *x, float *y);
	void circle_mask(float *array, float *x, float *y);
	void gaussian_mask(float *array, float *x, float *y);
	void secondary_mask(float *array, float *x, float *y);
	void corners_mask(float *array, float *x, float *y);
	void annulus_mask(float *array, float *x, float *y);

	// utility to get number of pixesl from fractions / pixels
	float pixels(double input);
	float getAverageWidth();
	
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

};

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){

	// check inputs!
	if (nrhs==0){
		mexPrintf("Usage: [flux, weight, background, variance, offset_x, offset_y, width] = photometry(cutouts, varargin)\n"); 
		mexPrintf("OPTIONAL ARGUMENTS:\n");
		mexPrintf("-------------------\n");
		mexPrintf("...\n");
		return;
	}
	
	// read the input data and parameters
	if(mxIsEmpty(prhs[0])){ // no input, then just return with all empty outputs...
		const mwSize dims[]={0,0};
		for(int i=0;i<7;i++) plhs[i]=mxCreateNumericArray(0,dims, mxSINGLE_CLASS, mxREAL); // return all empty arrays...
		return;
	}

	Photometry phot(nrhs, prhs);
	
	phot.run();
	plhs[0]=phot.flux_ptr; 
	plhs[1]=phot.weight_ptr; 
	plhs[2]=phot.background_ptr; 
	plhs[3]=phot.variance_ptr;
	plhs[4]=phot.offset_x_ptr; 
	plhs[5]=phot.offset_y_ptr; 
	plhs[6]=phot.width_ptr; 
	plhs[7]=phot.bad_pixel_ptr;
	
}

Photometry::Photometry(int nrhs, const mxArray *prhs[]){ // class constructor

	parseInputs(nrhs, prhs);
	
	// allocate intermediate arrays
	X=(float *) mxCalloc(N, sizeof(float));
	Y=(float *) mxCalloc(N, sizeof(float)); 
	
	for(int i=0;i<N;i++){ // meshgrid
		
		X[i]=(float)(i/dims[0])-dims[1]/2;
		Y[i]=(float)(i%dims[0])-dims[0]/2;
		
	}
	
	if(debug_bit>2){ // check that meshgrid returned what we expect
		printMatrix(X, "X");
		printMatrix(Y, "Y");
	}
	
}

Photometry::~Photometry(){ // destructor cleans up intermidiate arrays
	
	mxFree(X);
	mxFree(Y);
	
}

void Photometry::parseInputs(int nrhs, const mxArray *prhs[]){

	if(mxIsClass(prhs[0], "single")==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotFloat", "Input 1 to photometry is not a single/float array...");
	cutouts=(float*) mxGetData(prhs[0]);
	dims=(mwSize*)mxGetDimensions(prhs[0]);
	ndims=mxGetNumberOfDimensions(prhs[0]);	
	N=dims[0]*dims[1]; // dims[0] is the height while dims[1] is the width
	
	// size of non-empty output is [size(cutouts,3), size(cutouts,4)]
	if(ndims>2) out_dims[0]=dims[2];
	else out_dims[0]=1;
	if(ndims>3) out_dims[1]=dims[3];
	else out_dims[1]=1;
	num_cutouts=out_dims[0]*out_dims[1]; // how many values in each of the above arrays... 
	
	for(int i=1;i<nrhs;i+=2){ // parse varargin
		
		char key[STRLN];
		if(mxIsChar(prhs[i])==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotChar", "Input %d to photometry is not a string...", i+1);
		mxGetString(prhs[i], key, STRLN); // copy the string data
		
		mxArray *val=0;
		if(i+1<nrhs) val=(mxArray*) prhs[i+1]; // if the varargin is odd numbered, leave val=0 as default
		
		if(cs(key, "square")){
			
			// ap_func=&Photometry::square_mask; 
			ap_type=SQUARE;
			snprintf(aperture_string, STRLN, "SQUARE");
			
			if(val && mxIsEmpty(val)==0){ // check if there are any numerical parameters passed to "square" function
				if(mxIsNumeric(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not numeric...", i+2);
				ap_pars=(double*) mxGetData(val);
				num_ap_pars=mxGetNumberOfElements(val);
			} 
			
		}
		else if(cs(key, "circle", "aperture")){
			
			// ap_func=&Photometry::circle_mask; 
			ap_type=CIRCLE;
			snprintf(aperture_string, STRLN, "CIRCLE");
			
			if(val && mxIsEmpty(val)==0){ // check if there are any numerical parameters passed to "circle" function
				if(mxIsNumeric(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not numeric...", i+2);
				ap_pars=(double*) mxGetData(val);
				num_ap_pars=mxGetNumberOfElements(val);
			} 
			
			
		}		
		else if(cs(key, "gaussian", "psf")){
			
			// ap_func=&Photometry::gaussian_mask; 
			ap_type=GAUSSIAN;
			snprintf(aperture_string, STRLN, "GAUSSIAN");
			
			if(val && mxIsEmpty(val)==0){ // check if there are any numerical parameters passed to "gaussian" function
				if(mxIsNumeric(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not numeric...", i+2);
				ap_pars=(double*) mxGetData(val);
				num_ap_pars=mxGetNumberOfElements(val);
			} 
			
			
		}
		else if(cs(key, "corners")){
			
			// g_func=&Photometry::corners_mask; 
			bg_type=CORNERS;
			snprintf(background_string, STRLN, "CORNERS");
			
			if(val && mxIsEmpty(val)==0){ // check if there are any numerical parameters passed to "corners" function
				if(mxIsNumeric(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not numeric...", i+2);
				bg_pars=(double*) mxGetData(val);
				num_bg_pars=mxGetNumberOfElements(val);
			} 
			
			
		}
		else if(cs(key, "annulus")){
			
			// bg_func=&Photometry::annulus_mask; 
			bg_type=ANNULUS;
			snprintf(background_string, STRLN, "ANNULUS");
			
			if(val && mxIsEmpty(val)==0){ // check if there are any numerical parameters passed to "annulus" function
				if(mxIsNumeric(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not numeric...", i+2);
				bg_pars=(double*) mxGetData(val);
				num_bg_pars=mxGetNumberOfElements(val);
			} 
			
		}
		else if(cs(key, "fluxes")){
			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsEmpty(val)) continue;
			if(mxIsNumeric(val)==0 || mxIsSingle(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericSingle", "Input %d to photometry is not a numeric single precision float...", i+2);
			
			flux_ptr=mxDuplicateArray(val);
			
		}
		else if(cs(key, "weights")){
			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsEmpty(val)) continue;
			if(mxIsNumeric(val)==0 || mxIsSingle(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not a numeric single precision float...", i+2);
			
			weight_ptr=mxDuplicateArray(val);
			
		}
		else if(cs(key, "backgrounds")){
			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsEmpty(val)) continue;
			if(mxIsNumeric(val)==0 || mxIsSingle(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not a numeric single precision float...", i+2);
			
			background_ptr=mxDuplicateArray(val);
			
		}
		else if(cs(key, "variances")){
			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsEmpty(val)) continue;
			if(mxIsNumeric(val)==0 || mxIsSingle(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not a numeric single precision float...", i+2);
			
			variance_ptr=mxDuplicateArray(val);
			
		}
		else if(cs(key, "dx") || cs(key, "offset_x", 8) || cs(key, "offsets_x", 9)){
			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsEmpty(val)) continue;
			if(mxIsNumeric(val)==0 || mxIsSingle(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not a numeric single precision float...", i+2);
			
			offset_x_ptr=mxDuplicateArray(val);
			
		}
		else if(cs(key, "dy") || cs(key, "offset_y", 8) || cs(key, "offsets_y", 9)){
			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsEmpty(val)) continue;
			if(mxIsNumeric(val)==0 || mxIsSingle(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not a numeric single precision float...", i+2);
			
			offset_y_ptr=mxDuplicateArray(val);
			
		}
		else if(cs(key, "widths")){
			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsEmpty(val)) continue;
			if(mxIsNumeric(val)==0 || mxIsSingle(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not a numeric single precision float...", i+2);
			
			width_ptr=mxDuplicateArray(val);
			
		}
		else if(cs(key, "bad_pixels")){
			if(val==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			if(mxIsEmpty(val)) continue;
			if(mxIsNumeric(val)==0 || mxIsSingle(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumeric", "Input %d to photometry is not a numeric single precision float...", i+2);
			
			bad_pixel_ptr=mxDuplicateArray(val);
			
		}
		else if(cs(key, "subtract")){
			if(val==0 || mxIsEmpty(val)) subtract=1; // if no input, assume positive
			else{
				subtract=parse_bool(val); 
			}
			
		}
		else if(cs(key, "iterations", "num_iterations", "number_iterations")){
			if(val==0 || mxIsEmpty(val)) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			else{
				if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry is not a numeric scalar...", i+2);
				num_iter=mxGetScalar(val);
			}
			
		}
		else if(cs(key, "epsilon")){
			if(val==0 || mxIsEmpty(val)) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			else{
				if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry is not a numeric scalar...", i+2);
				epsilon=mxGetScalar(val);
			}
				
		}
		else if(cs(key, "debug_bit")){
			if(val==0 || mxIsEmpty(val)) debug_bit=1; // if no input, assume positive
			else{
				if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry is not a numeric scalar...", i+2);
				debug_bit=mxGetScalar(val);
			}
			
		}
		else if(cs(key, "multiplier")){
			if(val==0 || mxIsEmpty(val)) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			else{
				if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry is not a numeric scalar...", i+2);
				multiplier=mxGetScalar(val);
			}
			
		}
		else if(cs(key, "seeing")){
			if(val==0 || mxIsEmpty(val)) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:notEnoughInputs", "Expected varargin pair for %s at input", key, i+2);
			else{
				if(mxIsNumeric(val)==0 || mxIsScalar(val)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotNumericScalar", "Input %d to photometry is not a numeric scalar...", i+2);
				seeing=mxGetScalar(val);
			}
			
		}
		
	}
	
	if(flux_ptr==0) flux_ptr=mxCreateNumericArray(2, (const mwSize*) out_dims, mxSINGLE_CLASS, mxREAL);
	flux=(float*) mxGetData(flux_ptr);
	if(weight_ptr==0) weight_ptr=mxCreateNumericArray(2, (const mwSize*) out_dims, mxSINGLE_CLASS, mxREAL);
	weight=(float*) mxGetData(weight_ptr);
	if(background_ptr==0) background_ptr=mxCreateNumericArray(2, (const mwSize*) out_dims, mxSINGLE_CLASS, mxREAL);
	background=(float*) mxGetData(background_ptr);
	if(variance_ptr==0) variance_ptr=mxCreateNumericArray(2, (const mwSize*) out_dims, mxSINGLE_CLASS, mxREAL);
	variance=(float*) mxGetData(variance_ptr);
	if(offset_x_ptr==0) offset_x_ptr=mxCreateNumericArray(2, (const mwSize*) out_dims, mxSINGLE_CLASS, mxREAL);
	offset_x=(float*) mxGetData(offset_x_ptr);
	if(offset_y_ptr==0) offset_y_ptr=mxCreateNumericArray(2, (const mwSize*) out_dims, mxSINGLE_CLASS, mxREAL);
	offset_y=(float*) mxGetData(offset_y_ptr);
	if(width_ptr==0) width_ptr=mxCreateNumericArray(2, (const mwSize*) out_dims, mxSINGLE_CLASS, mxREAL);
	width=(float*) mxGetData(width_ptr);
	if(bad_pixel_ptr==0) bad_pixel_ptr=mxCreateNumericArray(2, (const mwSize*) out_dims, mxSINGLE_CLASS, mxREAL);
	bad_pixel=(float*) mxGetData(bad_pixel_ptr);
	
	if(debug_bit){ // check that all inputs have been received! 
		mexPrintf("cutouts: [");
		for(int i=0;i<ndims;i++){if(i>0) mexPrintf("x"); mexPrintf("%d", dims[i]); }
		mexPrintf("] (%s) ", mxGetClassName(prhs[0]));
		mexPrintf("| aperture: %s (%d) | pars= ", aperture_string, ap_type);
		for(int i=0;i<num_ap_pars;i++) mexPrintf("%4.2f ", ap_pars[i]);
		mexPrintf("| background: %s (%d) | pars= ", background_string, bg_type);
		for(int i=0;i<num_bg_pars;i++) mexPrintf("%4.2f ", bg_pars[i]);
		mexPrintf("| iter= %d | debug_bit= %d | subtract= %d | N= %d\n", num_iter, debug_bit, subtract, N);
	}
	
}

void Photometry::run(){

	// mexPrintf("num_cutouts= %d\n", num_cutouts);

	for(int j=0;j<num_cutouts;j++){ // number of cutouts
		
		calculate(j);
		
	} // for j
	
	if(debug_bit>2) printMatrix(cutouts, "cutouts");
	
}

void Photometry::calculate(int j){
	
	// mexPrintf("N= %d | iter= %d\n", N, num_iter);
	
	float *x=(float *)mxCalloc(N, sizeof(float)); // grid with shift added using offset_x
	float *y=(float *)mxCalloc(N, sizeof(float)); // grid with shift added using offset_y
	float *x2=(float *)mxCalloc(N, sizeof(float)); // grid with shift added using found moment
	float *y2=(float *)mxCalloc(N, sizeof(float)); // grid with shift added using found moment
	float *ap_array=(float *)mxCalloc(N, sizeof(float)); // grid with the aperture shape
	float *bg_array=(float *)mxCalloc(N, sizeof(float)); // grid with the background shape
	float *bad_array=(float *)mxCalloc(N, sizeof(float)); // grid with the bad pixels marked
	float *image_sub=(float *)mxCalloc(N, sizeof(float));
	float *image=(float *)mxCalloc(N, sizeof(float)); // can be a copy of image_raw or image_sub depending on value of "subtract" option
	float *image_raw=&cutouts[j*N]; // a pointer to the raw data
	
	if(ap_type==0 && bg_type==0) num_iter=1; // if using the simplest aperture/background we have nothing to gain from running iterations...
	
	for(int k=0;k<num_iter;k++){ // number of iterations 
		
		// if offsets exists from last iteration or last batch, use them
		float dx=offset_x[j];
		if(isnan(dx)) dx=0; // NaN values must be replaced with null position
		float dy=offset_y[j];
		if(isnan(dy)) dy=0; // NaN values must be replaced with null position
	
		for(int i=0;i<N;i++){ // x,y grid with offsets
			x[i]=X[i]-dx;
			y[i]=Y[i]-dy;
		}
		
		// this is the ugliest syntax I've ever seen
		//(this->*bg_func)(bg_array, x,y); // make an offset background mask
		//(this->*ap_func)(ap_array, x,y); // make an offset aperture mask
		
		main_mask(ap_array, x, y);
		secondary_mask(bg_array, x, y);
		
		// printMatrix(ap_array, "ap_array");
		// printMatrix(bg_array, "bg_array");
	
		// first calculate the b/g so we can subtract it! 
		float bg_weight=sumArrays(bg_array);

		if(bg_weight>0)	background[j]=sumArrays(image_raw, bg_array)/bg_weight; // average background value per pixels
		else background[j]=0;
		
		// now make a background subtracted image
		for(int i=0;i<N;i++) image_sub[i]=image_raw[i]-background[j];
		if(subtract) memcpy(image, image_sub, N*sizeof(float));
		else memcpy(image, image_raw, N*sizeof(float));
		variance[j]=sumArrays(image_sub, image_sub, bg_array)/bg_weight; // variance of the background area
		
		// use the variance to find negative outlier pixels and mark them as bad pixels
		float std=sqrt(variance[j]);
		for(int i=0;i<N;i++) if(image[i]<-3*std) image[i]=NAN; 
		
		// find the number of bad pixels in the aperture mask... 				
		bad_pixel[j]=0;
		for(int i=0;i<N;i++){
		
			bad_array[i]=0; // initialization 
			if(isnan(image[i]) && ap_array[i]>0){
				bad_array[i]=1; // find all the bad pixels in the image
				bad_pixel[j]++;
			}
		}
		
		for(int i=0;i<N;i++) if(bad_array[i]){ ap_array[i]=NAN; bg_array[i]=NAN; } // make sure the ap/bg arrays have NaNs everywhere the original image has NaNs
		
		weight[j]=sumArrays(ap_array); // number of pixels in this aperture
		
		for(int i=0;i<N;i++) ap_array[i]/=weight[j]; // normalize aperture
		float sum_ap_square=sumArrays(ap_array, ap_array); // normaliztion by sum(ap^2)
		
		for(int i=0;i<N;i++) image[i]*=ap_array[i]/sum_ap_square; // weigh by the normalized aperture
		
		float m0=sumArrays(image); // flux after going through aperture, normalized by sum(ap^2)

		 // debug output! 
		// mexPrintf("j= %d | ap_array[312]= %f | weight[j]= %f | sum_ap_square= %f | m0= %f | offset_x= %f | offset_y= %f\n", j, ap_array[312], weight[j], sum_ap_square, m0, offset_x[j], offset_y[j]);
		// if (isAllNaNs(image)) mexPrintf("Found image j= %d with all nans!\n", j);
		
		// calculate first moments
		float m1x=sumArrays(image, X)/m0;
		float m1y=sumArrays(image, Y)/m0;
		
		for(int i=0;i<N;i++){ // x,y grid with offsets using new moments
			if(isnan(m1x)==0) x2[i]=X[i]-m1x;
			if(isnan(m1y)==0) y2[i]=Y[i]-m1y;
		}
		
		// second moments (using the new grids)
		float m2x=sumArrays(image, x2, x2)/m0;
		float m2y=sumArrays(image, y2, y2)/m0;
		float mxy=sumArrays(image, x2, y2)/m0;
		
		flux[j]=m0;
		
		if(m0>0 && fabs(m1x)<dims[1]/2 && fabs(m1y)<dims[0]/2){ // values are reliable enough to fill the other parameters
			offset_x[j]=m1x;
			offset_y[j]=m1y;
			// width[j]=sqrt((m2x+m2y)/2);
			width[j]=getWidthFromMoments(m2x, m2y, mxy);
		}
		else{
			offset_x[j]=NAN;
			offset_y[j]=NAN;
			width[j]=NAN;
		}
		
		if(fabs(dx-m1x)<epsilon && fabs(dy-m1y)<epsilon) break; // if the iterations don't move the aperture by much, just stop iterating. 
		
	} // for k
	
	mxFree(x);
	mxFree(y);
	mxFree(x2);
	mxFree(y2);
	mxFree(ap_array);
	mxFree(bg_array);
	mxFree(bad_array);
	mxFree(image);
	mxFree(image_sub);
	
}

float Photometry::getWidthFromMoments(float m2x, float m2y, float mxy){
// got this little nugget from: https://yutsumura.com/express-the-eigenvalues-of-a-2-by-2-matrix-in-terms-of-the-trace-and-determinant/

	float tr=m2x+m2y;
	float det=m2x*m2y - mxy*mxy;
	
	float r1=(tr-sqrt(tr*tr-4*det))/2;
	float r2=(tr+sqrt(tr*tr-4*det))/2;
	
	return (sqrt(r1)+sqrt(r2))/2;

}

void Photometry::main_mask(float *array, float *x, float *y){

	if(ap_type==SQUARE) square_mask(array, x, y);
	else if(ap_type==CIRCLE) circle_mask(array, x, y);
	else if(ap_type==GAUSSIAN) gaussian_mask(array, x, y);

	// printMatrix(array, "main mask");
	
}

void Photometry::square_mask(float *array, float *x, float *y){
	
	for(int i=0;i<N;i++) array[i]=1;
	
}

void Photometry::circle_mask(float *array, float *x, float *y){
	
	float radius=sqrt(N)/2; // default radius of aperture
	
	if(num_ap_pars>0){ // got an override to default radius
		if(isnan(ap_pars[0])){ // NaN override means get the width from last time
			float w=getAverageWidth();
			if(w>0) radius=w*multiplier;
		}
		else radius=(float) pixels(ap_pars[0]);
	}
	
	
	
	// mexPrintf("radius= %f\n", radius);
	for(int i=0;i<N;i++){
		
		float r=(float) sqrt(x[i]*x[i]+y[i]*y[i]);
		
		array[i]=radius+0.5-r;
		//if(array[i]<0) array[i]=0;
		if(array[i]<1) array[i]=0;
		if(array[i]>1) array[i]=1;
		
	}
}

void Photometry::gaussian_mask(float *array, float *x, float *y){
	
	
	float sigma=sqrt(N)/2; // default gaussian width parameter
	float threshold=1e-6; // under this value is replaced with NaN (do we need this?)

	if(num_ap_pars>0){ // got an override to default sigma
		if(isnan(ap_pars[0])){ // NaN override means get the width from last time
			float w=getAverageWidth();
			if(w>0) sigma=w*multiplier;
		}
		else sigma=(float) pixels(ap_pars[0]);
	}
	
	for(int i=0;i<N;i++){
		
		float r=(float) sqrt(x[i]*x[i]+y[i]*y[i]);
		
		array[i]=exp(-0.5*r*r/sigma/sigma);
		
		// if(array[i]<threshold) array[i]=NAN;
		
	}
	
	
}

void Photometry::secondary_mask(float *array, float *x, float *y){

	if(bg_type==CORNERS) corners_mask(array, x, y);
	else if(bg_type==ANNULUS) annulus_mask(array, x, y);
	// printMatrix(array, "secondary mask");
	
}

void Photometry::corners_mask(float *array, float *x, float *y){
	
	float corner_size=0.15f;
	
	if(num_bg_pars>0) corner_size=bg_pars[0];
	
	corner_size=round(pixels(corner_size));
	
	for(int i=0;i<N;i++){
		if(i/dims[0]<corner_size && i%dims[0]<corner_size) array[i]=1;
		else if(i/dims[0]<corner_size && i%dims[0]>=dims[1]-corner_size) array[i]=1;
		else if(i/dims[0]>=dims[0]-corner_size && i%dims[0]<corner_size) array[i]=1;
		else if(i/dims[0]>=dims[0]-corner_size && i%dims[0]>=dims[1]-corner_size) array[i]=1;
		else array[i]=0;
	}
	
}

void Photometry::annulus_mask(float *array, float *x, float *y){
	
	
	float radius1=sqrt(N)/2-1;
	float radius2=1e10; // inf

	if(num_bg_pars>0) radius1=pixels(bg_pars[0]);
	if(num_bg_pars>1) radius2=pixels(bg_pars[1]);
	
	for(int i=0;i<N;i++){
		
		float r=(float) sqrt(x[i]*x[i]+y[i]*y[i]);
		
		if (r>radius1 && r<radius2)	array[i]=1;
		else array[i]=0; 
		
	}
	
}

float Photometry::pixels(double input){
	
	if(isnan(input)) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNaN", "Input to pixels is NaN!", input);
	if(input<=0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotPositive", "Input to pixels is %f<0", input);
	if(input>1)	return input;
	
	// assume 0<input<1
	// get minimum of dims[0] and dims[1]
	int d=dims[0];
	if(dims[1]<d) d=dims[1]; 
	return d*input; // get fraction of smaller dimension of cutouts
	
}

float Photometry::getAverageWidth(){

	if(width_ptr){
		float S=0;
		int counter=0;
		for(int j=0;j<num_cutouts;j++) if(isnan(width[j])==0) { S+=width[j]; counter++; }
		return S/counter; // mean of the widths, excluding NaNs
	}
	else return 0;

}

float Photometry::sumArrays(const float *array1){
	
	float S=0;
	
	for(int i=0;i<N;i++) if(isnan(array1[i])==0) S+=array1[i];
	
	return S;
	
}

float Photometry::sumArrays(const float *array1, const float *array2){
	
	float S=0;
	
	for(int i=0;i<N;i++) if(isnan(array1[i])==0 && isnan(array2[i])==0) S+=array1[i]*array2[i];
	
	return S;
	
}

float Photometry::sumArrays(const float *array1, const float *array2, const float *array3){
	
	float S=0;
	
	for(int i=0;i<N;i++) 
		if(isnan(array1[i])==0 && isnan(array2[i])==0 && isnan(array3[i])==0) 
			S+=array1[i]*array2[i]*array3[i];
	
	return S;
	
}

float Photometry::sumArrays(const float *array1, const float *array2, const float *array3, const float *array4){
	
	float S=0;
	
	for(int i=0;i<N;i++) 
		if(isnan(array1[i])==0 && isnan(array2[i])==0 && isnan(array3[i])==0 && isnan(array4[i])==0) 
			S+=array1[i]*array2[i]*array3[i]*array4[i];
	
	return S;
	
}

bool Photometry::isAllNaNs(const float *array){
	
	int num_vals=0;
	int num_nans=0;
	
	for(int i=0;i<N; i++){
		
		num_vals++;
		if (isnan(array[i])) num_nans++;
		
	}
	
	return num_vals==num_nans;
	
}

void Photometry::printMatrix(const int *array, const char *name){
	
	mexPrintf("%s= \n", name);
	for(int j=0;j<dims[0];j++){ 
		for(int i=0;i<dims[1];i++) mexPrintf("%d ", array[i*dims[0]+j]);
		mexPrintf("\n");
	}
	
	mexPrintf("\n");
	
	
}

void Photometry::printMatrix(const float *array, const char *name){
	
	mexPrintf("%s= \n", name);
	for(int j=0;j<dims[0];j++){ 
		for(int i=0;i<dims[1];i++) mexPrintf("%5.3f ", array[i*dims[0]+j]);
		mexPrintf("\n");
	}
	
	mexPrintf("\n");
	
}

bool cs(const char *keyword, const char *compare_str, int num_letters){
	
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

bool cs(const char *keyword, const char *str1, const char *str2, int num_letters){

	return cs(keyword, str1, num_letters) || cs(keyword, str2, num_letters);

}

bool cs(const char *keyword, const char *str1, const char *str2, const char *str3, int num_letters){

	return cs(keyword, str1, num_letters) || cs(keyword, str2, num_letters) || cs(keyword, str3, num_letters);

}

bool parse_bool(mxArray *value){
	
	if (mxIsEmpty(value)) return 0;
	if (mxIsScalar(value)==0) mexErrMsgIdAndTxt("MATLAB:util:img:photometry:inputNotScalar", "Input 1 to parse_bool must be a scalar!");
	if (mxIsLogical(value)) return mxGetScalar(value)!=0;
	if (mxIsNumeric(value)) return mxGetScalar(value)!=0;
	if (mxIsChar(value)) return cs(mxArrayToString(value), "yes", "on");
	
}