function ret=AT_InitialiseLibrary()

%unsigned int AT_InitialiseLibrary()
%
%Description :	Initialises the Andor SDK3 library
%
%Arguments	 :  NONE
%
%Return		 :  Check the help for return code meanings
%

ret = obs.cam.sdk.andorsdk3functions('AT_InitialiseLibrary');