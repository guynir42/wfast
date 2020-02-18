#ifndef ANDORCAMERA_H 
#define ANDORCAMERA_H

#include "atcore.h"
#include "atutility.h"

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
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <random>

class AndorCamera {
	
	public:
	
	AndorCamera();
	virtual ~AndorCamera();
	
	// SDK specific calls
	void startup();
	void finishup();
	void batch(int idx);
	void restart();
	
	unsigned long long int getTimestamps(unsigned char *buf, int ImageSizeBytes);
	unsigned long long int extractIntValue(unsigned char *buf, unsigned long int offset, int numByes);
	
	// general purpose functions:	
	void loadFromCamera(const mxArray *camera);
	void loadFromBuffers(mxArray *buffers);
	void loadFromOptionsStruct(const mxArray *options);
	void loop();
	
	// int waitForRecording(int idx);
	int waitForReading(int idx);
	int waitForWriting(int idx);
	
	// utilities
	double getPosixTime();	
	unsigned int getTimeout();
	
	void printout();
	void report_error(const char *string, int error_value, double *mex_flag_cam);
	void save_error(const char *string, int error_value, int batch_number);
	
	// parameters:
	unsigned int num_batches=1;
	unsigned int batch_size;
	
	AT_64 height;
	AT_64 width;
	AT_64 stride=0;
	AT_64 clockFreq=0;
	
	long int hndl;
	
	static const int NBUF=256; // maximum number of buffers ever possible... 
	
	int num_buffers;
	unsigned short int *images_ptrs[NBUF];
	double *timestamps_ptrs[NBUF];
	double *t_vec_ptrs[NBUF];
	double *mex_flag_cam; // just one vector
	double *mex_flag_record_ptrs[NBUF];
	double *mex_flag_read_ptrs[NBUF];
	double *mex_flag_write_ptrs[NBUF];
	
	double *index_rec_vector;
	
	int timeout=2000; // single frame timeout in milliseconds
	
	int use_async;
	int debug_bit;
	
	bool restart_clock=1;
	double expT=0;
	double *error_log=0;
	unsigned int current_batch=0;
	
	// internal data used by the SDK
	static const int NTEMPBUF=30; // how many buffers to give the zyla...
	
	alignas(8) unsigned char *temp_buffers[NTEMPBUF]={NULL}; // number of empty pointers for Zyla
	int image_size_bytes;
	
	
	
};

#endif