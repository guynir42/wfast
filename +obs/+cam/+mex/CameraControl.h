#ifndef CAMERACONTROL_H 
#define CAMERACONTROL_H

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

class CameraControl {
	
	public:
	
	CameraControl();
	virtual ~CameraControl();
	
	
	// these are overriden in each sub-class:
	virtual void startup()=0;
	virtual void finishup()=0;
	virtual void record(int idx)=0;
	virtual std::string cam_name()=0;
	
	// general purpose functions:	
	void loadFromCamera(const mxArray *camera);
	void loadFromBuffers(mxArray *buffers);
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
	unsigned int height;
	unsigned int width;
	unsigned int batch_size;
	
	long int hndl;
	
	static const int NBUF=30; // maximum number of buffers ever possible... 
	
	int num_buffers;
	unsigned short int *images_ptrs[NBUF];
	double *timestamps_ptrs[NBUF];
	double *t_vec_ptrs[NBUF];
	double *mex_flag_cam; // just one vector
	double *mex_flag_record_ptrs[NBUF];
	double *mex_flag_read_ptrs[NBUF];
	double *mex_flag_write_ptrs[NBUF];
	
	double *index_rec_vector;
	
	int timeout=1000; // single frame timeout in milliseconds
	
	int use_async;
	int debug_bit;
	
	double expT=0;
	double *error_log=0;
	int current_batch=0;
	// int frame_delay_ms=0;
		
};

#endif