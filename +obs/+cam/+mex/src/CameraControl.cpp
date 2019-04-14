#include "CameraControl.h"

CameraControl::CameraControl(){
	
	// mexPrintf("this is CameraControl default constructor!\n");

}

CameraControl::~CameraControl(){
	
	// printf("this is CameraControl destructor...\n");
	
}

void CameraControl::loadFromCamera(const mxArray *camera){
	
	mxArray *prop=0;
	
	prop=(mxArray *)mxGetProperty(camera, 0, "debug_bit");
	if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:CameraControl:cannotFindProperty", "Cannot find property 'debug_bit' inside CameraControl!");
	debug_bit=(int) mxGetScalar(prop);
	
	prop=(mxArray *)mxGetProperty(camera, 0, "use_async");
	if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:CameraControl:cannotFindProperty", "Cannot find property 'use_async' inside CameraControl!");
	use_async=(int) mxGetScalar(prop);
	
	prop=(mxArray *)mxGetProperty(camera, 0, "im_size");
	if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:CameraControl:cannotFindProperty", "Cannot find property 'use_async' inside CameraControl!");
	if(mxIsEmpty(prop)) mexErrMsgIdAndTxt( "MATLAB:obs:cam:CameraControl:imageSizeEmpty", "Must define an image size in CameraControl!");
	double *im_size=mxGetPr(prop);
	height = (unsigned int) im_size[0];
	width = (unsigned int) im_size[1];
	
	prop=(mxArray *)mxGetProperty(camera, 0, "batch_size");
	if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:CameraControl:cannotFindProperty", "Cannot find property 'batch_size' inside CameraControl!");
	batch_size=(unsigned int) mxGetScalar(prop);
			
	prop=(mxArray *)mxGetProperty(camera, 0, "hndl");
	if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:CameraControl:cannotFindProperty", "Cannot find property 'hndl' inside CameraControl!");
	hndl=(long int) mxGetScalar(prop);
	
	// prop=(mxArray *)mxGetProperty(camera, 0, "frame_delay_ms");
	// if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:CameraControl:cannotFindProperty", "Cannot find property 'frame_delay_ms' inside CameraControl!");
	// frame_delay_ms=(int) mxGetScalar(prop);
		
}

void CameraControl::loadFromBuffers(mxArray *buffers){
	
	num_buffers=(int) mxGetNumberOfElements(buffers);
	mxArray *prop=0;
	
	// get the info we need from buffers into pure C++ (and preallocate any buffers/images)
	for(int b=0;b<num_buffers;b++){
		
		if(mex_flag_cam[1]) break; // get the signal to stop...
		
		prop=(mxArray *) mxGetField(buffers, b, "mex_flag_record");
		if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:CameraControl:cannotFindProperty", "Cannot find property 'mex_flag_record' inside buffers struct!");
		mex_flag_record_ptrs[b]=mxGetPr(prop);
		
		prop=(mxArray *) mxGetField(buffers, b, "mex_flag_read");
		if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:CameraControl:cannotFindProperty", "Cannot find property 'mex_flag_read' inside buffers struct!");
		mex_flag_read_ptrs[b]=mxGetPr(prop);
		
		prop=(mxArray *) mxGetField(buffers, b, "mex_flag_write");
		if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:CameraControl:cannotFindProperty", "Cannot find property 'mex_flag_write' inside buffers struct!");
		mex_flag_write_ptrs[b]=mxGetPr(prop);
		
		prop=(mxArray *) mxGetField(buffers, b, "images");
		if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:CameraControl:cannotFindProperty", "Cannot find property 'images' inside buffers struct!");
		const size_t *dims=mxGetDimensions(prop);
		size_t batch_size_temp=1;
		if(mxGetNumberOfDimensions(prop)>=3) batch_size_temp=dims[2];
		
		if(mxIsEmpty(prop) || height!=dims[0] || width!=dims[1] || batch_size!=batch_size_temp || mxIsClass(prop, "uint16")==0){ // choose uint16 or other data type...
			
			if(waitForReading(b)) return; // just making sure we are not allocating over data that is getting read			
			if(waitForWriting(b)) return; // just making sure we are not allocating over data that is getting written
			
			// assuming MATLAB uses memory aligned to 8-byte boundary...
			if(debug_bit>2) mexPrintf("Allocating a new matrix...\n");

			images_ptrs[b]=(unsigned short int*) mxCalloc((int) (height*width*batch_size), 2); // use 2 to initialize a uint16
			mwSize dims[3] = {(mwSize) height, (mwSize) width, (mwSize) batch_size};
			mxArray *array=mxCreateNumericArray(3, dims, mxUINT16_CLASS, mxREAL);
			mxSetData(array, images_ptrs[b]);
			mxSetField(buffers, b, "images", array);
		}
		else{ 
			images_ptrs[b]=(unsigned short int*) mxGetData(prop);
		}
		
		prop=(mxArray *) mxGetField(buffers, b, "timestamps");
		if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:CameraControl:cannotFindProperty", "Cannot find property 'timestamps' inside buffers struct!");
		dims=mxGetDimensions(prop);
		
		if(mxIsEmpty(prop) || batch_size!=dims[0] || mxIsClass(prop, "double")==0){ // timestamps must be a double vector the same length as "batch_size"
			
			if(waitForReading(b)) return; // just making sure we are not allocating over data that is getting read
			if(waitForWriting(b)) return; // just making sure we are not allocating over data that is getting written
			
			timestamps_ptrs[b]=(double*) mxCalloc((int) batch_size, sizeof(double));
			mwSize dims[2] = {(mwSize) batch_size, (mwSize) 1};
			mxArray *array=mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
			mxSetData(array, timestamps_ptrs[b]);
			mxSetField(buffers, b, "timestamps", array);
		}
		else{ 
			timestamps_ptrs[b]=(double*) mxGetData(prop);
		}
		
		prop=(mxArray *) mxGetField(buffers, b, "t_vec");
		if(prop==NULL) mexErrMsgIdAndTxt( "MATLAB:obs:cam:CameraControl:cannotFindProperty", "Cannot find property 't_vec' inside buffers struct!");
		dims=mxGetDimensions(prop);
		
		if(mxIsEmpty(prop) || dims[0]!=1 || dims[1]!=3 || mxIsClass(prop, "double")==0){ // t_vec must be a 3 element double vector
			
			if(waitForReading(b)) return; // just making sure we are not allocating over data that is getting read
			if(waitForWriting(b)) return; // just making sure we are not allocating over data that is getting written
			
			t_vec_ptrs[b]=(double*) mxCalloc(3, sizeof(double));
			mwSize dims[2] = {(mwSize) 1, (mwSize) 3};
			mxArray *array=mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
			mxSetData(array, t_vec_ptrs[b]);
			mxSetField(buffers, b, "t_vec", array);
		}
		else{ 
			t_vec_ptrs[b]=(double*) mxGetData(prop);
		}
		
	}// loop over different buffers...
	
}

void CameraControl::loop(){
	
	if(mex_flag_cam[0]){
		if(debug_bit) printf("Camera is already recording... \n"); 
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
			if(debug_bit) printf("Camera: Received stop flag after %d batches... \n", current_batch);
			break;
		} 
		
		// make sure buffer is available for recording...
		if(waitForReading(idx)) break; // if buffer times out on reading, break the loop
		if(waitForWriting(idx)) break; // if buffer times out on writing, break the loop
		
		t_vec_ptrs[idx][0]=getPosixTime();
		
		if(mex_flag_read_ptrs[idx][0]==1) printf("THIS SHOULDNT HAPPEN!!!\n\n"); // we can only reach this point after buffer is released from reading
		
		// mex_flag_record_ptrs[idx][1]=; // not finished recording
		mex_flag_record_ptrs[idx][0]=1; // recording
		
		if(debug_bit>5) printf("Camera: Recording batch %d on buffer %d | record flag: %d %d | read flag: %d %d\n", current_batch+1, idx+1, 
			(int) mex_flag_record_ptrs[idx][0], (int) mex_flag_record_ptrs[idx][1], 
			(int) mex_flag_read_ptrs[idx][0], (int) mex_flag_read_ptrs[idx][1]);
		
		record(idx);
		
		// make sure to flag this buffer as UNREAD
		//mex_flag_read_ptrs[idx][0]=0; // not started reading
		//mex_flag_read_ptrs[idx][2]=1; // has data for reading
		//mex_flag_read_ptrs[idx][1]=0; // not done reading
		
		// mex_flag_record_ptrs[idx][1]=1; // finished recording
		mex_flag_record_ptrs[idx][0]=0; // not recording
		mex_flag_read_ptrs[idx][0]=1; // lock the buffer until it is read out...
		
		if(debug_bit>5) printf("Camera: finished recording to buffer %d... record_flag= %d %d | read_flag= %d %d\n", idx+1, 
						(int) mex_flag_record_ptrs[idx][0], (int) mex_flag_record_ptrs[idx][1], 
						(int) mex_flag_read_ptrs[idx][0], (int) mex_flag_read_ptrs[idx][1]);

		if(mex_flag_cam[2]) break;
		
		index_rec_vector[0]++; // this is so the index_rec in the matlab class outside is also updated...
		if(index_rec_vector[0]>=num_buffers) index_rec_vector[0]=0;
		
	}
	
	finishup();
	mex_flag_cam[0]=0;// this tells use the camera is now ready to start recording again...
	if(debug_bit) printf("Camera: Done looping on %d / %d batches...\n", current_batch, num_batches);
	delete this; // must do a cleanup in the asynchronous function! 
	
}

int CameraControl::waitForReading(int idx){
	
	// printf("getTimeout()= %u | batch_size= %u\n", getTimeout(), batch_size);
	unsigned int N=((unsigned int) getTimeout()*batch_size); // in miliseconds
	int mil=10; // time resolution in milliseconds
	
	for(int i=0;i<N/mil; i++){
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

int CameraControl::waitForWriting(int idx){
	
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
		
	printf("Camera: waitForWriting timeout after %d milliseconds...\n", N);
	return 1;
	
}

double CameraControl::getPosixTime(){
	
	unsigned long long milliseconds_since_epoch = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
	return ((double) milliseconds_since_epoch)/1000;
	
}

unsigned int CameraControl::getTimeout(){ // in milliseconds
	// printf("timeout= %d | expT= %f\n", timeout, expT);
	if (timeout>expT*5000) return timeout;
	else return (unsigned int) expT*5000;

}

void CameraControl::printout(){
	
	printf("Camera: debug_bit= %d | use_async= %d | im_size= %dx%d | batch_size= %d | hndl= %d | cam_name= %s | mex_flag_cam= %d %d %d %d\n", 
		debug_bit, use_async, height, width, batch_size, hndl, cam_name().c_str(), (int) mex_flag_cam[0], (int) mex_flag_cam[1], (int) mex_flag_cam[2], (int) mex_flag_cam[3]);
		
}

void CameraControl::report_error(const char *string, int error_value, double *mex_flag_cam){
	
	printf("ERROR in %s: return value %d... stopping camera! \n", string, error_value);
	mex_flag_cam[1]=1;
	mex_flag_cam[2]=error_value;
	
}

void CameraControl::save_error(const char *string, int error_value, int batch_number){
	
	printf("ERROR in %s: return value %d... ", string, error_value);
	
	if(error_log){
		printf("logging error for batch %d! \n", batch_number);
		error_log[batch_number-1]=error_value;
	}
	else printf("\n");
	
}

