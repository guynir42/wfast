#ifndef MYATTRIBUTE_H
#define MYATTRIBUTE_H

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

#define STRLN 256

class MyAttribute {

public:

	MyAttribute(); // default constructor makes an empty attribute with no name (it will not get written)
	MyAttribute(const char *name); // empty attribute constructor	
	void input(const char *name);
	MyAttribute(const char *name, double value); // scalar type attribute
	void input(const char *name, double value);
	MyAttribute(const char *name, double *data, int length); // vector type attribute
	void input(const char *name, double *data, int length);
	MyAttribute(const char *name, const char *string); // string type attribute
	void input(const char *name, const char *string);
	MyAttribute(const char *name, const mxArray *value);
	void input(const char *name, const mxArray *value);
	
	bool is_empty();
	void printout();
	
	char att_name[STRLN];
	int is_scalar=0;
	int is_vec=0;
	int is_str=0;
	
	bool use_this=0; // if this is zero we will not write this attribute to disk... 
	
	std::string str;
	double scalar;
	std::vector<double> vec;
	

};

#endif