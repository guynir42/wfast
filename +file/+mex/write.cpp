#include "mex.h"
#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <thread>
#include <chrono>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "SaveData.h"
#include "SaveDataHDF5.h"
#include "MyMatrix.h"
#include "MyAttribute.h"


const char *parse_extension(const char *filename);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){
					  
	if(nrhs<1){ // usage is detailed when function is called without arguments!
		SaveData::print_help();
		return;
	}
	
	// check minimal inputs
	if(nrhs<2) mexErrMsgIdAndTxt("MATLAB:file:mex:mexWrite:invalidNumInputs", "At least 2 arguments are required!");
	if(mxIsChar(prhs[0])!=1) mexErrMsgIdAndTxt( "MATLAB:file:mex:mexWrite:inputNotString", "filename must be a string.");
	if(mxGetM(prhs[1])!=1) mexErrMsgIdAndTxt( "MATLAB:file:mex:mexWrite:inputNotVector", "mex_flag must be a row vector.");
	if(mxGetN(prhs[1])!=2) mexErrMsgIdAndTxt( "MATLAB:file:mex:mexWrite:inputSizeMismatch", "mex_flag must have exactly 2 elements!");
	
	SaveData *s=0;
	const char *filename=mxArrayToString(prhs[0]);
	const char *extension=parse_extension(filename);
	
	if (strlen(extension)==0 || SaveData::cs(extension, "h5", "h5z", "hdf5")) s=new SaveDataHDF5(filename);
	else mexErrMsgIdAndTxt( "MATLAB:file:mex:mexWrite:unknownExtension", "filename extension is unknown. Try 'h5' or 'hdf5'...");
	
	s->mex_flag=mxGetPr(prhs[1]);
	s->parseVararginPairs(nrhs, prhs);
	
	if(s->debug_bit>2){// get back some feedback on inputs...

		mexPrintf("\nfilename: '%s'\n", s->filename);
		mexPrintf("mex_flag_write= %g %g\n", s->mex_flag[0], s->mex_flag[1]);
		
		s->images.printout();
		s->timestamps.printout();
		s->juldates.printout(); 
		s->cutouts.printout();
		s->positions.printout();
		s->coordinates.printout();
		s->magnitudes.printout();
		s->temperatures.printout();
		s->fluxes.printout();
		s->cutouts_bg.printout();
		s->positions_bg.printout();
		s->backgrounds.printout();
		s->stack.printout();
		s->psfs.printout();
		
		mexPrintf("debug= %d | deflate= %d | chunk= %d | async= %d \n\n", s->debug_bit, s->deflate, s->chunk_size, s->async_write);
		
	}
	
	if(s->debug_bit>5){ // print all the parameters in the cell array too
		
		for(int i=0;i<s->parameter_addresses_vector.size();i++){
			
			printf("METADATA OBJECT IN: %s\n", s->parameter_addresses_vector[i].c_str());
			
			for(int j=0;j<s->parameter_attributes_2D_vector[i].size(); j++){
				
				s->parameter_attributes_2D_vector[i][j].printout();
				
			}
			
		}
		
	}
	
	if(s->async_write==0) s->write();
	else{
		std::thread mythread(&SaveData::write, s);
		mythread.detach();
	}
	
	
}

const char *parse_extension(const char *filename){
	
	const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
	
}