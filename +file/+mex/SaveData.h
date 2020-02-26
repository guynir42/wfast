#ifndef SAVEDATA_H
#define SAVEDATA_H

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

#include "MyMatrix.h"
#include "MyAttribute.h"

class SaveData {
	
	public:
	
	SaveData(const char *filename);
	virtual ~SaveData();
	
	virtual void write();
	virtual void writeData()=0;
	virtual void writeHeader()=0;
    	
	void parseVararginPairs(int N, const mxArray *vars[]); // get data, metadata and switches from the keyword-value pairs
	void readStruct(const mxArray *buf); // parse the metadata add-on varaibles kept inside the BufferWheel
	void readHeaderCellArray(const mxArray *cell); // parse the output of util.oop.save when generating cell of address-struct pairs	
	void dataChecks();
	void setFilename(const char *name);
	static void print_help();
	void handle_errors(const char *short_desc, const char *desc);
	
	char filename[STRLN];	
	static double *mex_flag; // flag[0]: started writing, flag[1]: finished writing, flag[2]: error
	
	MyMatrix images;
	MyMatrix timestamps; // 1D vector of times	
	MyMatrix cutouts;
	MyMatrix positions;
	MyMatrix coordinates;
	MyMatrix magnitudes;
	MyMatrix temperatures;	
	MyMatrix fluxes;
	MyMatrix cutouts_bg;
	MyMatrix positions_bg;
	MyMatrix backgrounds;
	MyMatrix stack;
	MyMatrix psfs;
	
	// for writing parameters objects and other attribute/metadata classes:
	std::vector< std::string > parameter_addresses_vector; // names of each struct 	
	std::vector< std::vector<MyAttribute> > parameter_attributes_2D_vector; // a list of parameter objects (converted to structures) 
	
	static int debug_bit;
	static int deflate;
	static size_t chunk_size;	
	static int async_write;
		
	static bool cs(const char *keyword, const char *compare_str, int num_letters=3);
	static bool cs(const char *keyword, const char *str1, const char *str2, int num_letters=3);
	static bool cs(const char *keyword, const char *str1, const char *str2, const char *str3, int num_letters=3);
	static bool cs(const char *keyword, const char *str1, const char *str2, const char *str3, const char *str4, int num_letters=3);
	
};

#endif
