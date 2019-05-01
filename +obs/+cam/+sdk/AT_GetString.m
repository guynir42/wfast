function [ret,string]=AT_GetString(hndl,featurename,strLength)

%unsigned int AT_GetString(AT_H Hndl, const AT_WC* Feature, AT_WC* String, int StringLength)
%
%Description :	Gets the value of a given String Feature
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature you wish to check
%                       strLength - The maximum length of the string
%
%Return		 :  ret - Check the help for return code meanings
%                       string - The current value that the feature is set to

[ret,string] = obs.cam.sdk.andorsdk3functions('AT_GetString',hndl,featurename,strLength);