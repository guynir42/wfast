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
	
	mxArray *buf_struct=0;
	
	MyMatrix images;
	MyMatrix timestamps; // 1D vector of times	
	MyMatrix cutouts;
	MyMatrix positions;
	MyMatrix coordinates;
	MyMatrix magnitudes;
	MyMatrix temperatures;	
	MyMatrix cutouts_bg;
	MyMatrix positions_bg;
	MyMatrix stack;
	MyMatrix psfs;
	
	MyMatrix fluxes;
	MyMatrix errors;
	MyMatrix areas;
	MyMatrix backgrounds;
	MyMatrix variances;
	MyMatrix offsets_x;
	MyMatrix offsets_y;
	MyMatrix centroids_x;
	MyMatrix centroids_y;
	MyMatrix widths;
	MyMatrix bad_pixels;
	MyMatrix flags;
	
	// for writing parameters objects and other attribute/metadata classes:
	std::vector< std::string > parameter_addresses_vector; // names of each struct 	
	std::vector< std::vector<MyAttribute> > parameter_attributes_2D_vector; // a list of parameter objects (converted to structures) 
	
	// these are all static so that objects like MyDataset (defined in SaveDataHDF5.h) can access them
	static int debug_bit; // default=0. Control level of verbosity
	static int deflate; // default=0. Choose level of deflation/compression (usually 1 is enough)
	static size_t chunk_size; // default=64. Size of chunk will automatically shrink if it is bigger than image size (assume square chunks in image plane)
	static int async_write; // default=0. Use threads to write data while continuing work
	static int photometric_write; // default=1. Use this to save fluxes, backgrounds, widths, etc... 
	static int full_save; // default=1. Use this to save full-frame data cubes
	static int prefer_images_to_stack; // default=1. If true, and if stacks and images have only 2 dimensions, then will only save images, not stack. 
	// since these are static and non-const, we initialize them in the SaveData.cpp source file! 
	
	static bool cs(const char *keyword, const char *compare_str, int num_letters=3);
	static bool cs(const char *keyword, const char *str1, const char *str2, int num_letters=3);
	static bool cs(const char *keyword, const char *str1, const char *str2, const char *str3, int num_letters=3);
	static bool cs(const char *keyword, const char *str1, const char *str2, const char *str3, const char *str4, int num_letters=3);
	
};

#endif
