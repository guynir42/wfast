function [ret,isImplemented]=AT_IsImplemented(hndl,featurename)

%unsigned int AT_IsImplemented(AT_H Hndl, const AT_WC* Feature, AT_BOOL* Implemented)
%
%Description :	Checks if a given feature has been implemented
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature you wish to check
%
%Return		 :  ret - Check the help for return code meanings
%                       isImplemented - Boolean value indicating if the feature has been implemented

[ret,isImplemented] = obs.cam.sdk.andorsdk3functions('AT_IsImplemented',hndl,featurename);