function ret=AT_Flush(hndl)

%unsigned int AT_Flush(AT_H Hndl)
%
%Description :	Clears any buffers in both the input and output queues
%
%Arguments	 :  hndl - The handle for the selected camera 
%
%Return		 :  ret - Check the help for return code meanings
%                       

ret = obs.cam.sdk.andorsdk3functions('AT_Flush',hndl);