function ret=AT_RegisterFeatureCallback(hndl,featurename,callbackfunction)

%unsigned int AT_RegisterFeatureCallback(AT_H Hndl, const AT_WC* Feature, FeatureCallback EvCallback, void* Context)
%
%Description :	Registers a callback function to be called when a given feature changes
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature with which to register the callback
%                       callbackfunction - The function to be called when the feature is updated 
%
%Return		 :  ret - Check the help for return code meanings
%                       

ret = obs.cam.sdk.andorsdk3functions('AT_RegisterFeatureCallback',hndl,featurename,callbackfunction);