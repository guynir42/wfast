function [ret,isReadOnly]=AT_IsReadOnly(hndl,featurename)

%unsigned int AT_IsReadOnly(AT_H Hndl, const AT_WC* Feature, AT_BOOL* ReadOnly)
%
%Description :	Checks if a given feature is read only
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature you wish to check
%
%Return		 :  ret - Check the help for return code meanings
%                       isReadOnly - Boolean value indicating if the feature is read only

[ret,isReadOnly] = obs.cam.sdk.andorsdk3functions('AT_IsReadOnly',hndl,featurename);