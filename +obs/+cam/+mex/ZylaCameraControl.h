#ifndef ZYLACAMERACONTROL_H
#define ZYLACAMERACONTROL_H

#include "CameraControl.h"
#include "atcore.h"
#include "atutility.h"



class ZylaCameraControl : public CameraControl{
	
	public: 
	
	ZylaCameraControl();
	virtual ~ZylaCameraControl();
	
	void startup();
	void finishup();
	void record(int idx);
	void restart();
	unsigned long long int getTimestamps(unsigned char *buf, int ImageSizeBytes);
	unsigned long long int extractIntValue(unsigned char *buf, unsigned long int offset, int numByes);
	std::string cam_name(); // returns "Zyla"
		
	static const int NTEMPBUF=30; // how many buffers to give the zyla...
	
	alignas(8) unsigned char *temp_buffers[NTEMPBUF]={NULL}; // number of empty pointers for Zyla
	int image_size_bytes;
	
	private:
	AT_64 this_stride=0;
	AT_64 this_height=0;
	AT_64 this_width=0;
	AT_64 clockFreq=0;
	
};
	
#endif