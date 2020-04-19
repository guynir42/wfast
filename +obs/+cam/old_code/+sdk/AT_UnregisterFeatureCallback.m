function ret =AT_UnegisterFeatureCallback(hndl,featurename,callbackfunction)

%unsigned int AT_UnegisterFeatureCallback(AT_H Hndl, const AT_WC* Feature, FeatureCallback EvCallback, void* Context)
%
%Description :	Unregisters a callback function for a given feature changes
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature with which to unregister the callback
%                       callbackfunction - The function you wish to stop recieving callbacks 
%
%Return		 :  ret - Check the help for return code meanings
%                       

ret = obs.cam.sdk.andorsdk3functions('AT_UnregisterFeatureCallback',hndl,featurename,callbackfunction);