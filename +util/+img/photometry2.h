#ifndef PHOTOMETRY2_H
#define PHOTOMETRY2_H

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <cctype>
#include <thread>
#include <chrono>
#include <algorithm> 

#define STRLN 64 // maximum string length (for copying)
#define NUM_DATA_TYPES 10 // flux, area, error, background, variance, offset_x, offset_y, width, bad_pixels, flag

// utility functions to compare strings
bool cs(const char *keyword, const char *compare_str, int num_letters=3);
bool cs(const char *keyword, const char *str1, const char *str2, int num_letters=3);
bool cs(const char *keyword, const char *str1, const char *str2, const char *str3, int num_letters=3);
bool parse_bool(mxArray *value);

// Usage: [flux, error, area, background, variance, offset_x, offset_y, width] = photometry(cutouts, varargin)

class Photometry{
	
	public:
		
	float *cutouts=0; // the input matrix on which we calculate everything (dim 1&2 are y&x, dim 3 is frame number, dim 4 is star number)
	mwSize dims[4]={0}; // array with the sizes of cutouts
	mwSize ndims=0; // number of dimensions of cutouts
	int num_cutouts=0;  // number of cutouts (dim 3 times dim 4)
	int N=0; // number of pixels in each cutout (dim 1 times dim 2)
	
	double gain=1; // calculate the source noise using the gain
	double scintillation_fraction=0; // use this to add the estimated scintillation noise (in fractions of the reduced flux)
	int num_threads=1; // for future use with built-in multithreading
	int num_iterations=2; // how many iterations of repositioning should we do
	bool use_centering_aperture=1; // run one level of aperture photometry (centroids only) before the first gaussian iterations
	bool use_gaussian=1; // decide if you want to use gaussians photometry at all
	bool use_apertures=1; // decide if you want to use aperture photometry (wedding cake)
	bool use_forced=1; // decide if you want to use forced photometry after finding the best offsets
	
	int debug_bit=0;
	
	// these can be shared with all calculations (i.e., they are kept in memory between calls to photometry2!)
	float *X=0; // grid points not including shifts. 
	float *Y=0; // grid points not including shifts. 
	
	int resolution=1; // how many different aperture/gaussian shifts we want in each pixel
	mwSize shift_dims[2]={0}; // dimensions of shift matrices dx/dy	
	int num_shifts=0; // number of shifts we need to cover all possible star positions (i.e., the length of dx and dy, shift_dims[0]*shift_dims[1])
	float *dx=0; // list of offsets in x for the center of the mask
	float *dy=0; // list of offsets in y for the center of the mask
	
	float *forced=0; // 3D matrix of apertures used for forced photometry, one for each dx/dy shift
	std::vector<int> *forced_indices=0; // array of length "num_shifts" of index vectors telling what part of each matrix to sum in forced photometry
	double forced_radius=5; // the default radius used for forced photometry (pixels)

	float *apertures=0; // 4D matrix of aperture+annuli for the wedding cake photometry
	std::vector<int> *aperture_indices=0; // array of length "num_shifts" of index vectors telling what part of each matrix to sum in wedding cake photometry
	int num_radii=0; // how many different aperture arrays do we have for the wedding cake
	double *ap_radii=0; // radii of different apertures, given in pixel units
	
	float *annulii=0; // 3D matrix of the annulus used for all kinds of photometry, one for each dx/dy shift
	std::vector<int> *annulus_indices=0; // array of length "num_shifts" of index vectors telling what part of each matrix to sum in annulus calculation
	double inner_radius=10; // of the annulus (pixels)
	double outer_radius=0; // of the annulus (pixels)

	float *gaussians=0; // 3D matrix of gaussian weighted apertures for "PSF" photometry, one for each dx/dy shift
	double gauss_sigma=2; // the width of the gaussian (in pixels)
	
	// output arrays are defined here (in C++)
	mwSize output_size[2]={0}; // the same as dim 3 and dim 4 of cutouts
	float **output_raw=0;
	float **output_forced=0;
	float **output_apertures=0;
	float **output_gaussian=0;
	float *best_offset_x=0;
	float *best_offset_y=0;
	const static char data_types[NUM_DATA_TYPES][STRLN]; // this holds strings containing: flux, error, area, background, variance, offset_x, offset_y, width, bad_pixels, flag
	
	// function prototypes (implementation at the end)
	Photometry();
	~Photometry();
	void parseInputs(int nrhs, const mxArray *prhs[]);
	mxArray *outputStruct(float **output, int num_fluxes=1); // wrap up the output matrices as a nice matlab style array
	mxArray *outputMetadataStruct(); // add a struct with some of the parameters and the different aperture masks used
	mxArray *outputArraysStruct(); // add a struct with the actual masks and grid arrays
	mxArray *outputIndicesVectors(std::vector<int> *vectors, int num_radii=1); // produce a cell array with size num_shifts*num_radii (default=1), each with the list of indices for that mask
	
	void clear(); // make sure all the output arrays are NaN before filling them
	
	void makeArrays(); // create all the required memory if it isn't already allocated, make all the required masks
	void makeAllOutputs(); // generate all the output matrices needed
	void allocateOutputArray(float **&output, int num_fluxes=1); // allocate memory for a set of outputs: flux, area, error, background, variance, offset_x/y, width, bad_pixels, flag
	void initializeOutputArray(float **output, int num_fluxes=1); // re-use memory but make sure it is intialized to NaN first
	
	void makeGrids(); // make the X/Y coordinate grid and the dx/dy possible aperture shift matrices
	void makeMasks(); // make a bank of masks to be used for different shifts and types of photometry
	void meshgrid(float *x, float *y, mwSize *dims, int resolution=1);  // fill the two matrices with dimensions "dims" with centered grid points (use resolution>1 to divide each pixel into multiple shifts)
	void makeApertureMasks();
	void makeAnnulusMasks();
	void makeGaussianMasks();
	
	void deleteArrays(); // get rid of all the masks and output arrays (when destroying this object)
	void deleteGrids(); // get rid of the grid matrices
	void deleteOutputs(); // get rid of output arrays only
	void deleteOutputArray(float **&output); // go over and deallocate the memory for one of the outputs
	void deleteMasks(); // go over and deallocate the memory for all the masks
	void deleteApertureMasks();
	void deleteAnnulusMasks();
	void deleteGaussianMasks(); 
	
	void run();
	void run_idx(int start_idx, int end_idx); // use this to run only a subset of the cutouts (for multithreading)
	void calculate(int j); // do the actual work on some cutout array
	int getShiftIndex(float x, float y); // find the index closest to the specific shift value x and y in the shift matrices (dx and dy)
	float getError(float variance, float reduced_flux); // calculate the best estimate for the noise, including background noise, source noise, and scintillation
	void runForced(); 
	void runForced_idx(int start_idx, int end_idx);
	void calculateForced(int j); 
	bool checkMoments(float offset_x, float offset_y, float width); // returns 1 if there is a problem with the offsets or width

	// make averages over the results
	float getWidthFromMoments(float m2x, float m2y, float mxy); // from the eigenvalues of the 2nd moments
	float getAverageWidth(float **output);
	float getAverageOffsetX(float **output);	
	float getAverageOffsetY(float **output);
	
	// the sum of the product of array1...
	int countNaNs(const float *array); 
	float sumArrays(const float *array1);
	float sumArrays(const float *array1, const float *array2);
	float sumArrays(const float *array1, float offset1, const float *array2);
	float sumArrays(const float *array1, const float *array2, const float *array3);
	float sumArrays(const float *array1, float offset1, const float *array2, const float *array3);
	float sumArrays(const float *array1, const float *array2, float offset2, const float *array3, float offset3);	
	float sumArrays(const float *array1, float offset1, const float *array2, float offset2, const float *array3, float offset3);
	float sumArrays(const float *array1, const float *array2, const float *array3, const float *array4);
	float sumArrays(const float *array1, float offset1, const float *array2, float offset2, const float *array3, float offset3, const float *array4);

	int countNonNaNsIndices(const float *array1, const std::vector<int> *vector, int idx); // sum the number of non-NaN values in array1 on the indices in vector[idx]
	float sumIndices(const float *array1, const std::vector<int> *vector, int idx); // sum the values of array1 on the indices in vector[idx]
	float sumIndices(const float *array1, float offset1, const std::vector<int> *vector, int idx); 
	float sumIndices(const float *array1, const float *array2, const std::vector<int> *vector, int idx);
	float sumIndices(const float *array1, float offset1, const float *array2, const std::vector<int> *vector, int idx);
	float sumIndices(const float *array1, float offset1, const float *array2, float offset2, const std::vector<int> *vector, int idx);
	float sumIndices(const float *array1, float offset1, const float *array2, float offset2, const float *array3, float offset3, const std::vector<int> *vector, int idx);
	float medianIndices(const float *array, const std::vector<int> *vector, int idx); // find the median of the array points indicated by vector[idx]
	
	// print on screen
	void printMatrix(const int *array, const char *name);
	void printMatrix(const float *array, const char *name);

};

#endif