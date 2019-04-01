#ifndef MYMATRIX_H
#define MYMATRIX_H

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

#include "MyAttribute.h"

class MyMatrix{ // used to track the size of an image dataset
	
public:
	
	MyMatrix();
	MyMatrix(const char *name, const mxArray *image, bool deflate=0);
	void input(const char *name, const mxArray *image, bool deflate=0);
	void printout();
	bool is_empty();
	bool is_double();
	bool is_uint16();
	bool is_float();
	
	unsigned short int *matrix_uint16=0;
	double *matrix_double=0;
	float *matrix_float=0;
	
	
	char data_name[STRLN];
	int ndims=0;
	size_t *dims;
	size_t dims_c[4];
	size_t numel=0;	
	size_t rows=0;
	size_t cols=0;
	size_t frames=1;
	size_t cutouts=1;
	size_t bits=0;
	bool use_deflate=0; 
	
	std::vector<MyAttribute> attributes;
	
};

#endif