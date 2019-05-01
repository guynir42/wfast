function ret=AT_SetBool(hndl,featurename,value)

%unsigned int AT_SetBool(AT_H Hndl, const AT_WC* Feature, AT_BOOL Value)
%
%Description :	Sets the value of a given Boolean Feature
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature you wish to check
%                       value - The value that the feature is to be set to
%
%Return		 :  ret - Check the help for return code meanings
%                      

ret = obs.cam.sdk.andorsdk3functions('AT_SetBool',hndl,featurename,value);