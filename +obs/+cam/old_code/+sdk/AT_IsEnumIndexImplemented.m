function [ret,isImplemented]=AT_IsEnumIndexImplemented(hndl,featurename,index)

%unsigned int AT_IsEnumIndexImplemented(AT_H Hndl, const AT_WC* Feature, int Index, AT_BOOL* Implemented)
%
%Description :	Checks if a given enumeration of this feature has been implemented
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature you wish to check
%                       index - The index you wish to check
%
%Return		 :  ret - Check the help for return code meanings
%                       isImplemented - Boolean value indicating if the enumeration has been implemented

[ret,isImplemented] = obs.cam.sdk.andorsdk3functions('AT_IsEnumIndexImplemented',hndl,featurename,index);