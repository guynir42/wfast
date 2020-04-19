function [ret,count]=AT_GetEnumCount(hndl,featurename)

%unsigned int AT_GetEnumCount(AT_H Hndl,const  AT_WC* Feature, int* Count);
%
%Description :	Returns the count of a given Enumerated Feature
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature you wish to check
%
%Return		 :  ret - Check the help for return code meanings
%                       count - The current value that the feature is set to

[ret,count] = obs.cam.sdk.andorsdk3functions('AT_GetEnumCount',hndl,featurename);