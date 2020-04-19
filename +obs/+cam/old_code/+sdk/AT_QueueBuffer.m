function ret=AT_QueueBuffer(hndl,imageSize)

%unsigned int AT_QueueBuffer(AT_H Hndl, AT_U8* Ptr, int PtrSize)
%
%Description :	Queues an image buffer of a given size in preperation for an acquisition
%
%Arguments	 :  hndl - The handle for the selected camera 
%                       imageSize - The size in bytes of the image buffer you want to be queued
%
%Return		 :  ret - Check the help for return code meanings
%                       

ret = obs.cam.sdk.andorsdk3functions('AT_QueueBuffer',hndl,imageSize);