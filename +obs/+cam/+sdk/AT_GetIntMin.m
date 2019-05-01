function [ret,value]=AT_GetIntMin(hndl,featurename)

%unsigned int AT_GetIntMin(AT_H Hndl, const AT_WC* Feature, AT_64* Value)
%
%Description :	Gets the minimum value that a given Integer Feature can be set to
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature you wish to check
%
%Return		 :  ret - Check the help for return code meanings
%                       value - The minimum value that the feature can be set to

[ret,value] = obs.cam.sdk.andorsdk3functions('AT_GetIntMin',hndl,featurename);