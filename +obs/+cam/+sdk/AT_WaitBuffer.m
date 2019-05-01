function [ret,buf]=AT_WaitBuffer(hndl,timeout)

%unsigned int AT_WaitBuffer(AT_H Hndl, AT_U8** Ptr, int* PtrSize, unsigned int Timeout)
%
%Description :	Waits for an acquisition for a given time and returns the acquired data
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       timeout - The time in ms to wait for an acquisition before returning a timeout
%
%Return		 :  ret - Check the help for return code meanings
%                       buf - The acquired image data as a byte array
%

[ret,buf] = obs.cam.sdk.andorsdk3functions('AT_WaitBuffer',hndl,timeout);