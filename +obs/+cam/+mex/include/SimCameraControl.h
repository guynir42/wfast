#ifndef SIMCAMERACONTROL_H
#define SIMCAMERACONTROL_H

#include "CameraControl.h"
#include "mex.h"
#include "matrix.h"

#include <random>

class SimCameraControl : public CameraControl{
	
	public:
	SimCameraControl();
	virtual ~SimCameraControl();
	
	virtual void startup();
	virtual void finishup();
	virtual void record(int idx);
	virtual std::string cam_name();
	
	std::default_random_engine generator;
	
	long int batch_counter=1;
	
};

#endif