function [ret,maxLength]=AT_GetStringMaxLength(hndl,featurename)

%unsigned int AT_GetStringMaxLength(AT_H Hndl, const AT_WC* Feature, int* MaxStringLength)
%
%Description :	Gets the maximum length of a given String Feature
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature you wish to check
%
%Return		 :  ret - Check the help for return code meanings
%                       maxLength - The current value that the feature is set to

[ret,maxLength] = obs.cam.sdk.andorsdk3functions('AT_GetStringMaxLength',hndl,featurename);