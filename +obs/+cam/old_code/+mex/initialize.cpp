#include "mex.h"
#include "matrix.h"

#include <stdio.h>
#include "AndorCamera.h"

AndorCamera *ac=0; 

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){

	printf("Trying to connect to Andor camera!\n"); 
	
	if(ac==0) ac=new AndorCamera(); 
	
	printf("ImageSizeBytes= %f\n", ac->getImageSize()); 
	
}