function ret=AT_SetInt(hndl,featurename,value)

%unsigned int AT_SetInt(AT_H Hndl, const AT_WC* Feature, AT_64 Value)
%
%Description :	Sets the value of a given Integer Feature
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature you wish to check
%                       value - The value to set the feature to
%
%Return		 :  ret - Check the help for return code meanings
%                       

ret = obs.cam.sdk.andorsdk3functions('AT_SetInt',hndl,featurename,value);