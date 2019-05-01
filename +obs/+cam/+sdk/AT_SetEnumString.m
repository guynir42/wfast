function ret=AT_SetEnumString(hndl,featurename,string)

%unsigned int AT_SetEnumString(AT_H Hndl, const AT_WC* Feature, const AT_WC* String);
%
%Description :	Sets the value of a given Enumerated Feature
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature you wish to check
%                       string - The value to set the feature to
%
%Return		 :  ret - Check the help for return code meanings
%                       

ret = obs.cam.sdk.andorsdk3functions('AT_SetEnumString',hndl,featurename,string);