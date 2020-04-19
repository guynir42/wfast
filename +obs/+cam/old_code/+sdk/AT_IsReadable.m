function [ret,isReadable]=AT_IsReadable(hndl,featurename)

%unsigned int AT_IsReadable(AT_H Hndl, const AT_WC* Feature, AT_BOOL* Readable)
%
%Description :	Checks if a given feature is readable
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature you wish to check
%
%Return		 :  ret - Check the help for return code meanings
%                       isReable - Boolean value indicating if the feature is readable

[ret,isReadable] = obs.cam.sdk.andorsdk3functions('AT_IsReadable',hndl,featurename);