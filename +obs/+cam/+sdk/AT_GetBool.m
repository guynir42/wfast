function [ret,value]=AT_GetBool(hndl,featurename)

%unsigned int AT_GetBool(AT_H Hndl, const AT_WC* Feature, AT_BOOL* Value)
%
%Description :	Gets the value of a given Boolean Feature
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature you wish to check
%
%Return		 :  ret - Check the help for return code meanings
%                       value - The current value that the feature is set to

[ret,value] = obs.cam.sdk.andorsdk3functions('AT_GetBool',hndl,featurename);