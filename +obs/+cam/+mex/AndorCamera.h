#ifndef ANDORCAMERA_H
#define ANDORCAMERA_H

#include <stdio.h>
//#include <math.h>

#include "atcore.h"

#define STRLN 256

class __declspec(dllexport) AndorCamera {

public:
	
	AndorCamera(); // constructor also calls the initialization of library and camera
	~AndorCamera(); // destructor deallocates the camera handle and shuts down the library
	void connect(); // try this if you can't connect to camera on the first try in the constructor (doesn't initialize library!)
	void disconnect(); // deallocates the camera handle (not the library!). 

	// getters
	int getLatestErrorValue(); 
	char *getLatestErrorString(); 
	
	double isRunning(); 
	double getImageSize(); 
	
	
	// setters
	
	// commands
	void start(); 
	void stop(); 
	
	// private
	AT_H _hndl=0; // camera handle
	int _err_code=0; // from the latest call to an SDK function
	char *_err_str=0; // allocated in the constructor
	double _mex_flag=0; // pointer to array of numbers to communicate between matlab and C++ threads
	
	
};

#endif