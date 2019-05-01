function ret=AT_FinaliseLibrary()

%unsigned int AT_FinaliseLibrary()
%
%Description :	Finalises the Andor SDK3 library
%
%Arguments	 :  NONE
%
%Return		 :  Check the help for return code meanings
%

ret = obs.cam.sdk.andorsdk3functions('AT_FinaliseLibrary');