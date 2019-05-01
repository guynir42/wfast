function ret=AT_Close(hndl)

%unsigned int AT_Close(AT_H *Hndl) 
%
%Description :	Closes the connection with a given camera
%
%Arguments	 :  hndl - The camera handle you wish to close
%
%Return		 :  Check the help for return code meanings
%

ret = obs.cam.sdk.andorsdk3functions('AT_Close', hndl);