function [ret,value]=AT_GetIntMax(hndl,featurename)

%unsigned int AT_GetIntMax(AT_H Hndl, const AT_WC* Feature, AT_64* Value)
%
%Description :	Gets the maximum value that a given Integer Feature can be set to
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature you wish to check
%
%Return		 :  ret - Check the help for return code meanings
%                       value - The maximum value that the feature can be set to

[ret,value] = obs.cam.sdk.andorsdk3functions('AT_GetIntMax',hndl,featurename);