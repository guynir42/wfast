function [ret,isWritable]=AT_IsWritable(hndl,featurename)

%unsigned int AT_IsWritable(AT_H Hndl, const AT_WC* Feature, AT_BOOL* Writable)
%
%Description :	Checks if a given feature is writable
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature you wish to check
%
%Return		 :  ret - Check the help for return code meanings
%                       isWritable - Boolean value indicating if the feature is writable

[ret,isWritable] = obs.cam.sdk.andorsdk3functions('AT_IsWritable',hndl,featurename);