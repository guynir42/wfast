#include "SimCameraControl.h"

SimCameraControl::SimCameraControl() : CameraControl(){

	

}

SimCameraControl::~SimCameraControl(){

}

void SimCameraControl::startup(){

	// mex_flag_cam[0]=1; // this tells us the camera is now recording...
	generator.seed( (unsigned int) std::chrono::system_clock::now().time_since_epoch().count());
	
}

void SimCameraControl::finishup(){

	if(debug_bit>2) printf("stopping camera now!\n");
	// mex_flag_cam[0]=0;// this tells use the camera is now ready to start recording again...
	
}

void SimCameraControl::record(int idx){
	
	std::normal_distribution<double> distribution(104.0,2.5);
	// std::normal_distribution<double> distribution(100.0,1);
	
	for (int i=0; i<height*width*batch_size; i++) {
		short int number= (short int) round(distribution(generator));
		// short int number=batch_counter; // debugging only...
		images_ptrs[idx][i]=number;
	}
	
	images_ptrs[idx][0]=200+batch_counter;
	
	batch_counter++;
	
}

std::string SimCameraControl::cam_name(){
	
	return std::string("SimCam");
	
}	
	