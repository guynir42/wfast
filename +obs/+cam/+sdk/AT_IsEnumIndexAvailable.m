function [ret,isAvailable]=AT_IsEnumIndexAvailable(hndl,featurename,index)

%unsigned int AT_IsEnumIndexAvailable(AT_H Hndl, const AT_WC* Feature, int Index, AT_BOOL* Available)
%
%Description :	Checks if a given enumeration of this feature is available
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature you wish to check
%                       index - The index you wish to check
%
%Return		 :  ret - Check the help for return code meanings
%                       isAvailable - Boolean value indicating if the enumeration is available

[ret,isAvailable] = obs.cam.sdk.andorsdk3functions('AT_IsEnumIndexAvailable',hndl,featurename,index);