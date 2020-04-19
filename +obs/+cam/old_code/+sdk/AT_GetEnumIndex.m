function [ret,index]=AT_GetEnumIndex(hndl,featurename)

%unsigned int AT_GetEnumIndex(AT_H Hndl, const AT_WC* Feature, int* Value);
%
%Description :	Gets the curent selected index of a given Enumerated Feature
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature you wish to check
%
%Return		 :  ret - Check the help for return code meanings
%                       index - The current index that the feature is set to

[ret,index] = obs.cam.sdk.andorsdk3functions('AT_GetEnumIndex',hndl,featurename);