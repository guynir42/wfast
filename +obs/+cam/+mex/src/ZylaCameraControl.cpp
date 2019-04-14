#include "ZylaCameraControl.h"

ZylaCameraControl::ZylaCameraControl() : CameraControl(){
	
	
}

ZylaCameraControl::~ZylaCameraControl(){
	
	
}

void ZylaCameraControl::startup(){
	
	if(debug_bit) printf("starting Zyla...\n");
	
	AT_64 ImageSizeBytes; AT_GetInt((AT_H) hndl, L"ImageSizeBytes", &ImageSizeBytes); 
	//cast so that the value can be used in the AT_QueueBuffer function 
	image_size_bytes = (int) ImageSizeBytes;	
	int ret=0;
	
	mex_flag_cam[0]=1;
	
	// from the SDK manual page 47:
	// AT_64 this_stride=0;
	ret=AT_GetInt((AT_H) hndl, L"AOIStride", &this_stride);
	if(ret!=AT_SUCCESS){ report_error("startup>get stride", ret, mex_flag_cam); return; }
	
	if(ret!=AT_SUCCESS){ report_error("startup>get height", ret, mex_flag_cam); return; }
	
	// AT_64 this_width=0;
	// AT_64 this_height=0;	
	ret=AT_GetInt((AT_H) hndl, L"AOIHeight", &this_height);
	ret=AT_GetInt((AT_H) hndl, L"AOIWidth", &this_width);
	if(ret!=AT_SUCCESS){ report_error("startup>get width", ret, mex_flag_cam); return; }
	
	ret=AT_GetFloat((AT_H) hndl, L"ExposureTime", &expT);
	if(ret!=AT_SUCCESS){ report_error("startup>get expT", ret, mex_flag_cam); return; }

	// AT_64 clockFreq=0;
	ret=AT_GetInt((AT_H) hndl, L"TimestampClockFrequency", &clockFreq);
	if(ret!=AT_SUCCESS){ report_error("startup>get clockFreq", ret, mex_flag_cam); return; }
	
	ret=AT_Flush((AT_H) hndl);
	if(ret!=AT_SUCCESS){ report_error("startup>flush buffers", ret, mex_flag_cam); return; }
	
	for(int i=0;i<NTEMPBUF;i++){
		
		temp_buffers[i]=new unsigned char[image_size_bytes+8]; // add 8 for byte alignment
		ret=AT_QueueBuffer((AT_H) hndl, temp_buffers[i], image_size_bytes);
		if(ret!=AT_SUCCESS){ report_error("startup>queue buffers", ret, mex_flag_cam); return; }
		
	}
	
	ret = AT_Command((AT_H) hndl, L"TimestampClockReset");
	if(ret!=AT_SUCCESS){ report_error("startup>clock reset", ret, mex_flag_cam); return; }
	
	// unsigned long long milliseconds_since_epoch =
    // std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
	// printf("Current POSIX time is: %lld.%lld seconds\n", milliseconds_since_epoch/1000, milliseconds_since_epoch%1000);
	
	ret=AT_Command((AT_H) hndl, L"ACquisitionStart");
	if(ret!=AT_SUCCESS){ report_error("startup>acquisition start", ret, mex_flag_cam); return; }
		
}

void ZylaCameraControl::finishup(){

	int ret = 0;	
	
	if (debug_bit) printf("stopping Zyla...\n");
	
	ret=AT_Command((AT_H) hndl, L"ACquisitionStop");	
	if(ret!=AT_SUCCESS){ report_error("finishup>stop acquisition", ret, mex_flag_cam); return; }
	
	ret=AT_Flush((AT_H) hndl);
	if(ret!=AT_SUCCESS){ report_error("finishup>flush buffers", ret, mex_flag_cam); return; }
	
	for(int i=0;i<NTEMPBUF;i++){ 
		if(temp_buffers[i]) delete temp_buffers[i];
		temp_buffers[i]=0;
	}
	
	mex_flag_cam[0]=0; // mark that camera is done recording... 
	
}

void ZylaCameraControl::record(int idx){
	
	// ANDOR SDK CODE HERE 
	//assumes library is initialized when CameraControl is constructed (and handle is openned)
	// and buffers are allocated in startup()		
	AT_64 ImageSizeBytes = 0;
	AT_GetInt((AT_H) hndl, L"ImageSizeBytes", &ImageSizeBytes);
	int ret=0;
	
	unsigned char *buf=0;
	int return_buf_size=0;
	
	ret = AT_Command((AT_H) hndl, L"SoftwareTrigger");
	if(ret!=AT_SUCCESS){ report_error("record>software trigger", ret, mex_flag_cam); return; }
	
	AT_64 clock=0;
		
	for(unsigned int i=0; i<batch_size; i++){
		
		int error_counter=0; // to make sure we are not in an endless loop... 
		
		ret=AT_WaitBuffer((AT_H) hndl, &buf, &return_buf_size, getTimeout());
		if(ret!=AT_SUCCESS){ // This is a special error case. We report it but also restart the acquisition... 
			// report_error("record>wait buffer", ret, mex_flag_cam); 
			// printf("WARNING: record>wait buffer: return value %d... restarting acquisition!\n", ret); 
			save_error("record>wait buffer", ret, current_batch); 
			restart();
			i--;
			error_counter++;
			if(error_counter>10) report_error("record>wait buffer", ret, mex_flag_cam);
			continue;
		}
		
		
		ret = AT_Command((AT_H) hndl, L"SoftwareTrigger");
		if(ret!=AT_SUCCESS){ report_error("record>software trigger", ret, mex_flag_cam); return; }
		
		// if(frame_delay_ms) std::this_thread::sleep_for(std::chrono::milliseconds(frame_delay_ms));
		
		int pos=i*height*width; // where inside the buffer to start writing now... 		
				
		ret=AT_ConvertBuffer(buf, (AT_U8*)(images_ptrs[idx]+pos), this_width, this_height, this_stride, L"Mono16", L"Mono16");
		if(ret!=AT_SUCCESS){ report_error("record>convert buffer", ret, mex_flag_cam); return; }
		
		clock=getTimestamps(buf, (int) ImageSizeBytes);
		timestamps_ptrs[idx][i]=((double)clock)/clockFreq;
		
		// unsigned short int *temp=(images_ptrs[idx]+pos); // debugging only!!!
		// temp=(unsigned short int*) buf; // debugging only!!!
		
		// printf("i= %d | queue buffer...\n", i);
		ret=AT_QueueBuffer((AT_H) hndl, buf, static_cast<int> (ImageSizeBytes));
		if(ret!=AT_SUCCESS){ report_error("record>queue buffer", ret, mex_flag_cam); return; }
		
	}
	
	// get the system clock value (in POSIX time)
	t_vec_ptrs[idx][1]=getPosixTime();
	
	// add a final timestamp that corresponds to this system time...
	// clock=getTimestamps(buf, (int) ImageSizeBytes);	// this gets the start of last frame!
	
	AT_GetInt((AT_H) hndl, L"TimestampClock", &clock); // current time! 
	t_vec_ptrs[idx][2]=((double)clock)/clockFreq;
	
}

void ZylaCameraControl::restart(){
	
	if (debug_bit) printf("restarting Zyla...\n");
	
	int ret=0;
	
	ret=AT_Command((AT_H) hndl, L"ACquisitionStop");	
	if(ret!=AT_SUCCESS){ report_error("restart>stop acquisition", ret, mex_flag_cam); return; }
	
	ret=AT_Flush((AT_H) hndl);
	if(ret!=AT_SUCCESS){ report_error("restart>flush buffers", ret, mex_flag_cam); return; }

	for(int i=0;i<NTEMPBUF;i++){
		
		ret=AT_QueueBuffer((AT_H) hndl, temp_buffers[i], image_size_bytes);
		if(ret!=AT_SUCCESS){ report_error("restart>queue buffers", ret, mex_flag_cam); return; }
				
	}
			
	ret=AT_Command((AT_H) hndl, L"ACquisitionStart");
	if(ret!=AT_SUCCESS){ report_error("restart>acquisition start", ret, mex_flag_cam); return; }
	
	ret = AT_Command((AT_H) hndl, L"SoftwareTrigger");
	if(ret!=AT_SUCCESS){ report_error("record>software trigger", ret, mex_flag_cam); return; }
	
}

unsigned long long int ZylaCameraControl::getTimestamps(unsigned char *buf, int ImageSizeBytes){
// this is copied pretty much as-is from AT_GetTimeStamp.m from the Andor matlab SDK.
// I imagine their code is copied directly from the Solis source code (it looks a lot like C code to me)

	unsigned long long int clock=0;

	int LENGTH_FIELD_SIZE = 4;
	int CID_FIELD_SIZE = 4;
	int CID_FPGA_TICKS = 1;
	// int AT_SUCCESS = 0;
	// int AT_ERR_NODATA = 11;

	// Get length of timestamp
	unsigned long int offset = (int) ImageSizeBytes - LENGTH_FIELD_SIZE;
	int timeStampLength = (int) extractIntValue(buf, offset,LENGTH_FIELD_SIZE)- CID_FIELD_SIZE;

	// Get Metadata identifier
	offset = offset - CID_FIELD_SIZE;
	int cid = (int) extractIntValue(buf, offset,CID_FIELD_SIZE);

	// Extract Metadata if present
	int rc=0;
	if (cid == CID_FPGA_TICKS){
		offset = offset - timeStampLength;
		clock =  extractIntValue(buf, offset,timeStampLength);		
	}
	else{
		clock=0;
		report_error("record>get timestamp", AT_ERR_NODATA, mex_flag_cam);
	}

	return clock;
	
}

unsigned long long int ZylaCameraControl::extractIntValue(unsigned char *buf, unsigned long int offset, int numBytes){
	
	unsigned long long int value=0;
	
	for(int i=0; i<numBytes; i++){ // value = value + uint64(array(offset+i)) * uint64(power(2, 8*(i-1)));
	
		unsigned long long int number=buf[offset+i];
		value+=number*(unsigned long long int)pow(2,(8*i));
	
	}

	return value;
	
}

std::string ZylaCameraControl::cam_name(){
	
	return std::string("Zyla");
	
}