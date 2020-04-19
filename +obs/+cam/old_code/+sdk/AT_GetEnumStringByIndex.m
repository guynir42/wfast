function [ret,string]=AT_GetEnumStringByIndex(hndl,featurename,index,strLength)

%unsigned int AT_GetEnumStringByIndex(AT_H Hndl, const AT_WC* Feature, int Index, AT_WC* String, int StringLength)
%
%Description :	Gets the value of a given Integer Feature
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature you wish to check
%                       index - The enumeration index which you wish to query
%                       strLength - The maximum string length
%
%Return		 :  ret - Check the help for return code meanings
%                       string - The current string value of the enumeration for the given index

[ret,string] = obs.cam.sdk.andorsdk3functions('AT_GetEnumStringByIndex',hndl,featurename,index,strLength);