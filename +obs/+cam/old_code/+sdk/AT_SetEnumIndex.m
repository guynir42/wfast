function ret=AT_SetEnumIndex(hndl,featurename,index)

%unsigned int AT_SetEnumIndex(AT_H Hndl, const AT_WC* Feature, AT_64 Value)
%
%Description :	Sets the value of a given Enumerated Feature to the selected index
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the feature you wish to check
%                       index - The index to set the feature to
%
%Return		 :  ret - Check the help for return code meanings
%                       

ret = obs.cam.sdk.andorsdk3functions('AT_SetEnumIndex',hndl,featurename,index);