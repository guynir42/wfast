function [ret,hndl]=AT_Open(index)

%unsigned int AT_Open(int CameraIndex, AT_H *Hndl)
%
%Description :	Opens the connection to an Andor SDK3 camera
%
%Arguments	 :  index - The camera index you wish to open
%
%Return		 :  ret - Check the help for return code meanings
%                       hndl - The camera handle that has been opened

[ret, hndl] = obs.cam.sdk.andorsdk3functions('AT_Open',index);