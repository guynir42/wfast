#include "AndorCamera.h"

AndorCamera::AndorCamera(){
	
	// printf("this is AndorCamera default constructor!\n");

}

AndorCamera::~AndorCamera(){
	
	// printf("this is AndorCamera destructor...\n");
	
}

void AndorCamera::loadFromCamera(const mxArray *camera){ // load control switches from camera object
	
	mxArray *prop=0;
	
	prop=(mxArray *)mxGetProperty(camera, 0, "debug_bit");
	if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:AndorCamera:cannotFindProperty", "Cannot find property 'debug_bit' inside AndorCamera!");
	debug_bit=(int) mxGetScalar(prop);
	
	prop=(mxArray *)mxGetProperty(camera, 0, "use_async");
	if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:AndorCamera:cannotFindProperty", "Cannot find property 'use_async' inside AndorCamera!");
	use_async=(int) mxGetScalar(prop);
	
	// prop=(mxArray *)mxGetProperty(camera, 0, "im_size");
	// if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:AndorCamera:cannotFindProperty", "Cannot find property 'im_size' inside AndorCamera!");
	// if(mxIsEmpty(prop)) mexErrMsgIdAndTxt( "MATLAB:obs:cam:AndorCamera:imageSizeEmpty", "Must define an image size in AndorCamera!");
	// double *im_size=mxGetPr(prop);
	// height = (unsigned int) im_size[0];
	// width = (unsigned int) im_size[1];
	
	prop=(mxArray *)mxGetProperty(camera, 0, "batch_size");
	if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:AndorCamera:cannotFindProperty", "Cannot find property 'batch_size' inside AndorCamera!");
	batch_size=(unsigned int) mxGetScalar(prop);
			
	prop=(mxArray *)mxGetProperty(camera, 0, "hndl");
	if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:AndorCamera:cannotFindProperty", "Cannot find property 'hndl' inside AndorCamera!");
	hndl=(long int) mxGetScalar(prop);
	
	int ret=0;
	ret=AT_GetInt((AT_H) hndl, L"AOIHeight", &height); 
	if(ret) mexErrMsgIdAndTxt( "MATLAB:obs:cam:AndorCamera:getHeight", "Cannot get heigh from device.");
	ret=AT_GetInt((AT_H) hndl, L"AOIWidth", &width);
	if(ret) mexErrMsgIdAndTxt( "MATLAB:obs:cam:AndorCamera:getWidth", "Cannot get width from device.");
	ret=AT_GetInt((AT_H) hndl, L"AOIStride", &stride);
	if(ret) mexErrMsgIdAndTxt( "MATLAB:obs:cam:AndorCamera:getStride", "Cannot get stride from device.");

	int trigger_index=0;
	ret=AT_GetEnumIndex((AT_H) hndl, L"TriggerMode", &trigger_index); 
	if(ret) mexErrMsgIdAndTxt( "MATLAB:obs:cam:AndorCamera:getTriggerIndex", "Cannot get trigger mode index from device.");
	wchar_t trigger_mode[256];
	ret=AT_GetEnumStringByIndex((AT_H) hndl, L"TriggerMode", trigger_index, trigger_mode, 256);
	if(ret) mexErrMsgIdAndTxt( "MATLAB:obs:cam:AndorCamera:getTriggerMode", "Cannot get trigger mode from device.");
	// printf("trigger mode= %S\n", trigger_mode); 
	if(wcscmp(L"Software", trigger_mode)==0) use_software_trigger=1;
	else use_software_trigger=0;
	
}

void AndorCamera::loadFromBuffers(mxArray *buffers){ // load mex flags and data arrays from buffer struct
	
	num_buffers=(int) mxGetNumberOfElements(buffers);
	mxArray *prop=0;
	
	// get the info we need from buffers into pure C++ (and preallocate any buffers/images)
	for(int b=0;b<num_buffers;b++){
		
		if(mex_flag_cam[1]) break; // get the signal to stop...
		
		prop=(mxArray *) mxGetField(buffers, b, "mex_flag_record");
		if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:AndorCamera:cannotFindProperty", "Cannot find property 'mex_flag_record' inside buffers struct!");
		mex_flag_record_ptrs[b]=mxGetPr(prop);
		
		prop=(mxArray *) mxGetField(buffers, b, "mex_flag_read");
		if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:AndorCamera:cannotFindProperty", "Cannot find property 'mex_flag_read' inside buffers struct!");
		mex_flag_read_ptrs[b]=mxGetPr(prop);
		
		prop=(mxArray *) mxGetField(buffers, b, "mex_flag_write");
		if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:AndorCamera:cannotFindProperty", "Cannot find property 'mex_flag_write' inside buffers struct!");
		mex_flag_write_ptrs[b]=mxGetPr(prop);
		
		prop=(mxArray *) mxGetField(buffers, b, "images");
		if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:AndorCamera:cannotFindProperty", "Cannot find property 'images' inside buffers struct!");
		const size_t *dims=mxGetDimensions(prop);
		size_t batch_size_temp=1;
		if(mxGetNumberOfDimensions(prop)>=3) batch_size_temp=dims[2];
		
		if(mxIsEmpty(prop) || height!=dims[1] || width!=dims[0] || batch_size!=batch_size_temp || mxIsClass(prop, "uint16")==0){ // choose uint16 or other data type...
			
			if(waitForReading(b)) return; // just making sure we are not allocating over data that is getting read			
			if(waitForWriting(b)) return; // just making sure we are not allocating over data that is getting written
			
			// assuming MATLAB uses memory aligned to 8-byte boundary...
			if(debug_bit>2) printf("Allocating a new matrix...\n");

			//images_ptrs[b]=(unsigned short int*) mxCalloc((int) (height*width*batch_size), 2); // use 2 to initialize a uint16
			mwSize dims[3] = {(mwSize) width, (mwSize) height, (mwSize) batch_size}; // width <-> height exchanged for C <-> matlab conversion
			mxArray *array=mxCreateNumericArray(3, dims, mxUINT16_CLASS, mxREAL);
			//mxSetData(array, images_ptrs[b]);
			images_ptrs[b]=(unsigned short int *) mxGetData(array);
			mxSetField(buffers, b, "images", array);
		}
		else{ 
			images_ptrs[b]=(unsigned short int*) mxGetData(prop);
		}
		
		prop=(mxArray *) mxGetField(buffers, b, "timestamps");
		if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:AndorCamera:cannotFindProperty", "Cannot find property 'timestamps' inside buffers struct!");
		dims=mxGetDimensions(prop);
		
		if(mxIsEmpty(prop) || batch_size!=dims[0] || mxIsClass(prop, "double")==0){ // timestamps must be a double vector the same length as "batch_size"
			
			if(waitForReading(b)) return; // just making sure we are not allocating over data that is getting read
			if(waitForWriting(b)) return; // just making sure we are not allocating over data that is getting written
			
			// timestamps_ptrs[b]=(double*) mxCalloc((int) batch_size, sizeof(double));
			mwSize dims[2] = {(mwSize) batch_size, (mwSize) 1};
			mxArray *array=mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
			// mxSetData(array, timestamps_ptrs[b]);
			timestamps_ptrs[b]=(double*) mxGetData(array);
			mxSetField(buffers, b, "timestamps", array);
		}
		else{ 
			timestamps_ptrs[b]=(double*) mxGetData(prop);
		}
		
		prop=(mxArray *) mxGetField(buffers, b, "t_vec");
		if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:AndorCamera:cannotFindProperty", "Cannot find property 't_vec' inside buffers struct!");
		dims=mxGetDimensions(prop);
		
		if(mxIsEmpty(prop) || dims[0]!=1 || dims[1]!=3 || mxIsClass(prop, "double")==0){ // t_vec must be a 3 element double vector
			
			if(waitForReading(b)) return; // just making sure we are not allocating over data that is getting read
			if(waitForWriting(b)) return; // just making sure we are not allocating over data that is getting written
			
			// t_vec_ptrs[b]=(double*) mxCalloc(3, sizeof(double));
			mwSize dims[2] = {(mwSize) 1, (mwSize) 3};
			mxArray *array=mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
			// SetData(array, t_vec_ptrs[b]);
			t_vec_ptrs[b]=(double*) mxGetData(array);
			mxSetField(buffers, b, "t_vec", array);
		}
		else{ 
			t_vec_ptrs[b]=(double*) mxGetData(prop);
		}
		
	}// loop over different buffers...
	
}

void AndorCamera::loadFromOptionsStruct(const mxArray *options){ // load additional options to override the Camera object (e.g., expT, ROI, batch_size)

	// need to implement this! 

}

void AndorCamera::startup(){
	
	if(debug_bit>1) printf("starting Andor camera...\n");
	
	int ret=0;
	
	ret=AT_SetBool(hndl, L"MetadataEnable", 1); 
	if(ret!=AT_SUCCESS){ report_error("startup>set metadata", ret, mex_flag_cam); return; }

	ret=AT_SetBool(hndl, L"MetadataTimestamp", 1); 
	if(ret!=AT_SUCCESS){ report_error("startup>set metadata timestamps", ret, mex_flag_cam); return; }
	
	AT_64 ImageSizeBytes; AT_GetInt((AT_H) hndl, L"ImageSizeBytes", &ImageSizeBytes); 
	//cast so that the value can be used in the AT_QueueBuffer function 
	image_size_bytes = (int) ImageSizeBytes;	
	
	mex_flag_cam[0]=1;
	
	// ret=AT_GetInt((AT_H) hndl, L"AOIStride", &this_stride);
	// if(ret!=AT_SUCCESS){ report_error("startup>get stride", ret, mex_flag_cam); return; }
	
	
	//ret=AT_GetInt((AT_H) hndl, L"AOIHeight", &this_height);
	// if(ret!=AT_SUCCESS){ report_error("startup>get height", ret, mex_flag_cam); return; }
	
	//ret=AT_GetInt((AT_H) hndl, L"AOIWidth", &this_width);
	//if(ret!=AT_SUCCESS){ report_error("startup>get width", ret, mex_flag_cam); return; }
	
	ret=AT_GetFloat((AT_H) hndl, L"ExposureTime", &expT);
	if(ret!=AT_SUCCESS){ report_error("startup>get expT", ret, mex_flag_cam); return; }

	// AT_64 clockFreq=0;
	ret=AT_GetInt((AT_H) hndl, L"TimestampClockFrequency", &clockFreq);
	if(ret!=AT_SUCCESS){ report_error("startup>get clockFreq", ret, mex_flag_cam); return; }
	
	ret=AT_Flush((AT_H) hndl);
	if(ret!=AT_SUCCESS){ report_error("startup>flush buffers", ret, mex_flag_cam); return; }
	
	for(int i=0;i<NTEMPBUF;i++){
		
		// temp_buffers[i]=new unsigned char[image_size_bytes+8]; // add 8 for byte alignment
		temp_buffers[i]=(unsigned char*) mxCalloc(image_size_bytes,1); 
		// aligned_buffers[i]=reinterpret_cast<unsigned char*>((reinterpret_cast<unsigned long long>(temp_buffers[i]) + 7) & ~0x7);
		ret=AT_QueueBuffer((AT_H) hndl, temp_buffers[i], image_size_bytes);
		if(ret!=AT_SUCCESS){ report_error("startup>queue buffers", ret, mex_flag_cam); return; }
		
	}
	
	if(restart_clock){
		ret = AT_Command((AT_H) hndl, L"TimestampClockReset");
		if(ret!=AT_SUCCESS){ report_error("startup>clock reset", ret, mex_flag_cam); return; }
	}
	
	// unsigned long long milliseconds_since_epoch =
    // std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
	// printf("Current POSIX time is: %lld.%lld seconds\n", milliseconds_since_epoch/1000, milliseconds_since_epoch%1000);
	
	ret=AT_Command((AT_H) hndl, L"AcquisitionStart");
	std::this_thread::sleep_for(std::chrono::milliseconds(100));
	
	if(use_software_trigger){
		ret = AT_Command((AT_H) hndl, L"SoftwareTrigger");
		if(ret!=AT_SUCCESS){ report_error("batch>software trigger", ret, mex_flag_cam); return; }
	}
	
	if(ret!=AT_SUCCESS){ report_error("startup>acquisition start", ret, mex_flag_cam); return; }
		
}

void AndorCamera::finishup(){

	int ret = 0;	
	
	if (debug_bit>1) printf("stopping Andor camera...\n");
	
	ret=AT_Command((AT_H) hndl, L"ACquisitionStop");	
	if(ret!=AT_SUCCESS){ report_error("finishup>stop acquisition", ret, mex_flag_cam); return; }
	
	ret=AT_Flush((AT_H) hndl);
	if(ret!=AT_SUCCESS){ report_error("finishup>flush buffers", ret, mex_flag_cam); return; }
	
	for(int i=0;i<NTEMPBUF;i++){ 
		// if(temp_buffers[i]) delete temp_buffers[i];
		if(temp_buffers[i]) mxFree(temp_buffers[i]);
		temp_buffers[i]=0;
		//aligned_buffers[i]=0; 
	}
	
	mex_flag_cam[0]=0; // mark that camera is done recording... 
	
}

void AndorCamera::loop(){ // start producing images in a loop over batches
	
	if(mex_flag_cam[0]){
		if(debug_bit>1) printf("Camera is already recording... \n"); 
		return; 
	}
	
	mex_flag_cam[0]=1; // this tells us the camera is now recording...
	startup();
	
	// now go into the long loop of how many batches you want to shoot
	for(current_batch=0; current_batch<num_batches; current_batch++){
		
		int idx = (int) index_rec_vector[0];
		
		// printf("IDX= 2 | flag: %d %d %d %d\n",(int) mex_flag_read_ptrs[1][0], (int) mex_flag_read_ptrs[1][1], (int) mex_flag_read_ptrs[1][2], (int) mex_flag_read_ptrs[1][3]);
		
		if(mex_flag_cam[1]){ // get the signal to stop...
			mex_flag_cam[0]=0; // unlock for later
			mex_flag_cam[1]=0; // unlock for later
			if(debug_bit>1) printf("Camera: Received stop flag after %d batches... \n", current_batch);
			break;
		} 
		
		// make sure buffer is available for recording...
		if(waitForReading(idx)) break; // if buffer times out on reading, break the loop
		if(waitForWriting(idx)) break; // if buffer times out on writing, break the loop
		
		t_vec_ptrs[idx][0]=getPosixTime();
		
		// if(mex_flag_read_ptrs[idx][0]==1) printf("THIS SHOULDNT HAPPEN!!!\n\n"); // we can only reach this point after buffer is released from reading
		
		// mex_flag_record_ptrs[idx][1]=; // not finished recording
		mex_flag_record_ptrs[idx][0]=1; // recording
		
		if(debug_bit>5) printf("Camera: Recording batch % 4d on buffer %d | record flag: %d %d | read flag: %d %d\n", current_batch+1, idx+1, 
			(int) mex_flag_record_ptrs[idx][0], (int) mex_flag_record_ptrs[idx][1], 
			(int) mex_flag_read_ptrs[idx][0], (int) mex_flag_read_ptrs[idx][1]);
		
		batch(idx);
		
		// make sure to flag this buffer as UNREAD
		//mex_flag_read_ptrs[idx][0]=0; // not started reading
		//mex_flag_read_ptrs[idx][2]=1; // has data for reading
		//mex_flag_read_ptrs[idx][1]=0; // not done reading
		
		// mex_flag_record_ptrs[idx][1]=1; // finished recording
		mex_flag_record_ptrs[idx][0]=0; // not recording
		mex_flag_read_ptrs[idx][0]=1; // lock the buffer until it is read out...
		
		if(debug_bit>5) printf("Camera: finished recording to buffer   %d | record_flag= %d %d | read_flag= %d %d\n", idx+1, 
						(int) mex_flag_record_ptrs[idx][0], (int) mex_flag_record_ptrs[idx][1], 
						(int) mex_flag_read_ptrs[idx][0], (int) mex_flag_read_ptrs[idx][1]);

		if(mex_flag_cam[2]) break;
		
		index_rec_vector[0]++; // this is so the index_rec in the matlab class outside is also updated...
		if(index_rec_vector[0]>=num_buffers) index_rec_vector[0]=0;
		
	}
	
	finishup();
	mex_flag_cam[0]=0;// this tells use the camera is now ready to start recording again...
	if(debug_bit>1) printf("Camera: Done looping on %d / %d batches...\n", current_batch, num_batches);
	delete this; // must do a cleanup in the asynchronous function! 
	
}

void AndorCamera::batch(int idx){
	
	// ANDOR SDK CODE HERE 
	// assumes library is initialized when CameraControl is constructed (and handle is openned)
	// and buffers are allocated in startup()		
	AT_64 ImageSizeBytes = 0;
	AT_GetInt((AT_H) hndl, L"ImageSizeBytes", &ImageSizeBytes);
	int ret=0;
	
	int return_buf_size=0; // the WaitBuffer command returns a buffer size that should be the same as ImageSizeBytes
	
	AT_64 clock=0;

	//ret = AT_Command((AT_H) hndl, L"SoftwareTrigger");
	//if(ret!=AT_SUCCESS){ report_error("batch>software trigger", ret, mex_flag_cam); return; }

	for(unsigned int i=0; i<batch_size; i++){ // loop over a 100 images in a batch
		
		unsigned char *buf=0;
		int error_counter=0; // to make sure we are not in an endless loop... 

		ret=AT_WaitBuffer((AT_H) hndl, &buf, &return_buf_size, getTimeout());
		
		if(ret!=AT_SUCCESS){ // This is a special error case. We report it but also restart the acquisition... 
			// report_error("batch>wait buffer", ret, mex_flag_cam); 
			if(debug_bit>1) printf("WARNING: batch>wait buffer: return value %d... restarting acquisition!\n", ret); 
			save_error("batch>wait buffer", ret, current_batch); 
			restart();
			i--;
			error_counter++;
			if(error_counter>12) report_error("batch>wait buffer", ret, mex_flag_cam);
			continue;
		}
		
		unsigned long long int pos=i*height*width; // where inside the buffer to start writing now... 		
		
		// ret=AT_ConvertBuffer(buf, (AT_U8*)(images_ptrs[idx]+pos), width, height, stride, L"Mono16", L"Mono16");
		//if(ret!=AT_SUCCESS){ report_error("batch>convert buffer", ret, mex_flag_cam); return; }
		
		if(use_software_trigger){ // only true in software trigger mode
			ret = AT_Command((AT_H) hndl, L"SoftwareTrigger");
			if(ret!=AT_SUCCESS){ report_error("batch>software trigger", ret, mex_flag_cam); return; }
		}
		
		memset(images_ptrs[idx]+pos, 0, height*width*2); // initialize the output array image_ptrs that is used by matlab
		
		// memcpy(images_ptrs[idx]+pos, buf, height*width*2);
		for(int j=0;j<height;j++){ // copy each row
			// memcpy(images_ptrs[idx]+pos+j*width, buf+j*stride, width*2);
			memmove(images_ptrs[idx]+pos+width*j, buf+stride*j, width*2);
		}
		
		// for(int j=0;j<height;j++){ // add a signal to each row to verify which image belongs to which frame
			// *(images_ptrs[idx]+pos+j*width+i)=128; // for debugging only
		// }
		
		clock=getTimestamps(buf, (int) ImageSizeBytes); // timestamp in FPGA clock ticks
		timestamps_ptrs[idx][i]=((double)clock)/clockFreq; // timestamps in seconds (passed out to matlab)
		
		// re-initialize the values in the buf to zero
		memset(buf, 0, ImageSizeBytes);

		ret=AT_QueueBuffer((AT_H) hndl, buf, static_cast<int> (ImageSizeBytes));
		if(ret!=AT_SUCCESS){ report_error("batch>queue buffer", ret, mex_flag_cam); return; }

		//memset(images_ptrs[idx]+pos, i, height*width*2); // this is for debugging only! 
	
	} // for i (each frame in the batch)
	
	// to pass in the first pixel of each image the number of repeated pixels (debugging only!)
	// for(int i=0;i<80;i++){
	
		// long long int S=0;
		// for(int j=0;j<width*height;j++){
			// if(*(images_ptrs[idx]+(i+10)*width*height+j)==*(images_ptrs[idx]+i*width*height+j)) S++; 
		// }
		// *(images_ptrs[idx]+i*width*height)=S/1000; 
	// }
	
	// get the system clock value (in POSIX time)
	t_vec_ptrs[idx][1]=getPosixTime();
	
	// add a final timestamp that corresponds to this system time...
	// clock=getTimestamps(buf, (int) ImageSizeBytes);	// this gets the start of last frame!
	
	AT_GetInt((AT_H) hndl, L"TimestampClock", &clock); // current time! 
	t_vec_ptrs[idx][2]=((double)clock)/clockFreq;
	
}

void AndorCamera::restart(){
	
	if (debug_bit>1) printf("restarting Andor camera...\n");
	
	int ret=0;
	
	ret=AT_Command((AT_H) hndl, L"ACquisitionStop");	
	if(ret!=AT_SUCCESS){ report_error("restart>stop acquisition", ret, mex_flag_cam); return; }
	
	ret=AT_Flush((AT_H) hndl);
	if(ret!=AT_SUCCESS){ report_error("restart>flush buffers", ret, mex_flag_cam); return; }

	for(int i=0;i<NTEMPBUF;i++){
		
		ret=AT_QueueBuffer((AT_H) hndl, temp_buffers[i], image_size_bytes);
		if(ret!=AT_SUCCESS){ report_error("restart>queue buffers", ret, mex_flag_cam); return; }
				
	}
			
	ret=AT_Command((AT_H) hndl, L"AcquisitionStart");
	if(ret!=AT_SUCCESS){ report_error("restart>acquisition start", ret, mex_flag_cam); return; }
	if(use_software_trigger){
		ret = AT_Command((AT_H) hndl, L"SoftwareTrigger");
		if(ret!=AT_SUCCESS){ report_error("batch>software trigger", ret, mex_flag_cam); return; }
	}
	
}

unsigned long long int AndorCamera::getTimestamps(unsigned char *buf, int ImageSizeBytes){ // recover the timestamps from the metadata
// this is copied pretty much as-is from AT_GetTimeStamp.m from the Andor matlab SDK.
// I imagine their code is copied directly from the Solis source code (it looks a lot like C code to me)
// Improved this for the Balor using the manual chapter 4.6 about metadata

	unsigned long long int clock=0;

	int LENGTH_FIELD_SIZE=4;
	int CID_FIELD_SIZE=4;
	int CID_FRAME_INFO=7; 
	int CID_FPGA_TICKS=1;
	// int AT_SUCCESS = 0;
	// int AT_ERR_NODATA = 11;
	
	unsigned long int offset = (int) ImageSizeBytes; // the initial position is at the end of the buffer! 
	
	for(int i=0;i<10;i++){ // how many chunks can we expect to read here...?
	
		// Get length of this metadata chunk
		offset-=LENGTH_FIELD_SIZE; // positition (in bytes) inside image buffer where we want to read the "length" block
		int chunk_length = (int) extractIntValue(buf, offset, LENGTH_FIELD_SIZE); // read how many bytes are in the rest of the block, includng metadata+CID
		if(offset<0) return 0;
		
		// Get Metadata identifier
		offset-= CID_FIELD_SIZE;
		int cid = (int) extractIntValue(buf, offset,CID_FIELD_SIZE);
		if(offset<0) return 0;
		
		chunk_length-=CID_FIELD_SIZE; // now we want to know how much to go back, after already having moved back CID_FIELD_SIZE
		offset-=chunk_length; // now we are at the position to read that metadata! 
		if(offset<0) return 0;
		
		// Extract Metadata if present
		int rc=0;
		if(cid == CID_FRAME_INFO){ // includes frame height/width/stride and other useless stuff
			continue;
		}
		else if (cid == CID_FPGA_TICKS){// this is what we actually want, the timestamp data! 
			
			clock =  extractIntValue(buf, offset, chunk_length); // chunk length should be 8 bytes, just for the FPGA clock counter
			return clock; 
		}
		else if(cid == 0){ // we reached the image data block! 
			clock=0;
			report_error("batch>get timestamp", AT_ERR_NODATA, mex_flag_cam);
			return 0;
		}

		return 0; // this shouldn't happen!
		
	}
		
}

unsigned long long int AndorCamera::extractIntValue(unsigned char *buf, unsigned long int offset, int numBytes){ // get values from bytes 
	
	unsigned long long int value=0;
	
	for(int i=0; i<numBytes; i++){ // value = value + uint64(array(offset+i)) * uint64(power(2, 8*(i-1)));
	
		unsigned long long int number=buf[offset+i];
		value+=number*(unsigned long long int)pow(2,(8*i));
	
	}

	return value;
	
}

int AndorCamera::waitForReading(int idx){ // make sure buffer is not locked while running analysis on data
	
	// printf("getTimeout()= %u | batch_size= %u\n", getTimeout(), batch_size);
	unsigned int N=((unsigned int) getTimeout()*batch_size); // in miliseconds
	unsigned int mil=10; // time resolution in milliseconds
	
	for(unsigned int i=0;i<N/mil; i++){
		 // if we haven't started reading or already finished reading this buffer, or there is no data in there, just get out with success
		// if((mex_flag_read_ptrs[idx][0]==0 || mex_flag_read_ptrs[idx][1]==1) && mex_flag_read_ptrs[idx][2]==0) return 0;
		if(mex_flag_read_ptrs[idx][0]==0) return 0; // if released for reading (or never started reading)
		
		if(mex_flag_cam[1]==1) return 1; // camera stop command 

		mex_flag_read_ptrs[idx][1]++; // increase counter (how many times we got delayed over reading)

		if(i==0){ 
		
			if(debug_bit>2) printf("Camera: waiting for buffer %d to finish reading... read flag: %d %d | record flag: %d %d\n", 
			idx+1, (int) mex_flag_read_ptrs[idx][0], (int) mex_flag_read_ptrs[idx][1],
			(int) mex_flag_record_ptrs[idx][0], (int) mex_flag_record_ptrs[idx][0]);
			
		}
		
		std::this_thread::sleep_for(std::chrono::milliseconds(mil));
			
	}
		
	if (debug_bit>5) printf("Camera: waitForReading timeout after %d milliseconds...\n", N);
	return 1;
	
}

int AndorCamera::waitForWriting(int idx){ // make sure buffer is not locked while dumping data to disk
	
	int N=((int)getTimeout()*batch_size); // in miliseconds
	int mil=10; // time resolution in milliseconds
	
	for(int i=0;i<N/mil; i++){
		
		//if(mex_flag_write_ptrs[idx][0]==0 || mex_flag_write_ptrs[idx][1]==1) return 0; // if we haven't started writing or already finished writing this buffer, just get out with success
		if(mex_flag_write_ptrs[idx][0]==0) return 0; // not locked for writing		
		
		if(mex_flag_cam[1]==1) return 1; // camera stop command 

		mex_flag_write_ptrs[idx][1]++; // increase counter (how many times we got delayed over writing)

		if(i==0){ 
			if(debug_bit>2) printf("Camera: waiting for buffer %d to finish writing... write flag: %d %d\n", idx+1, (int) mex_flag_write_ptrs[idx][0], (int) mex_flag_write_ptrs[idx][1]);
			
		}

		std::this_thread::sleep_for(std::chrono::milliseconds(mil));
			
	}
		
	if(debug_bit>1) printf("Camera: waitForWriting timeout after %d milliseconds...\n", N);
	return 1;
	
}

double AndorCamera::getPosixTime(){ // get time from computer's clock 
	
	unsigned long long milliseconds_since_epoch = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
	return ((double) milliseconds_since_epoch)/1000;
	
}

unsigned int AndorCamera::getTimeout(){ // in milliseconds
	
	if (timeout>expT*10000) return timeout;
	else return (unsigned int) (expT*10000);

}

void AndorCamera::printout(){ // print to screen some info on the camera switches
	
	printf("Camera: debug_bit= %d | use_async= %d | im_size= %dx%d | batch_size= %d | hndl= %d | mex_flag_cam= %d %d %d %d\n", 
		debug_bit, use_async, height, width, batch_size, hndl, (int) mex_flag_cam[0], (int) mex_flag_cam[1], (int) mex_flag_cam[2], (int) mex_flag_cam[3]);
	
}

void AndorCamera::report_error(const char *string, int error_value, double *mex_flag_cam){ // print error on screen and set error flag
	
	if(debug_bit>1) printf("ERROR in %s: return value %d... stopping camera! \n", string, error_value);
	mex_flag_cam[1]=1;
	mex_flag_cam[2]=error_value;
	
}

void AndorCamera::save_error(const char *string, int error_value, int batch_number){ // save error flag in vector of errors
	
	return; // debug only
	
	printf("ERROR in %s: return value %d... ", string, error_value);
	
	if(error_log){
		printf("logging error for batch %d! \n", batch_number);
		error_log[batch_number-1]=error_value;
	}
	else printf("\n");
	
}

