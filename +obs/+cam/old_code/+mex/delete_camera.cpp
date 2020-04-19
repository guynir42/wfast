
#include <stdio.h>
#include "mex.h"
#include "matrix.h"

#include "AndorCamera.h"

extern AndorCamera *ac; 

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){

	printf("Deleting Andor camera!\n"); 
	
	if(ac==0){ delete ac; ac=0; }
	
	printf("Deleting complete!\n"); 
	
}