function ret=AT_Command(hndl,featurename)

%unsigned int AT_Command(AT_H Hndl, const AT_WC* Feature)
%
%Description :	Sends a command to the SDK represented by the given feature
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       featurename - The name of the command feature 
%
%Return		 :  ret - Check the help for return code meanings
%                       

ret = obs.cam.sdk.andorsdk3functions('AT_Command',hndl,featurename);