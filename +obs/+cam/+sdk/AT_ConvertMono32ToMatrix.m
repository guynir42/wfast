function [ret,matrix]=AT_ConvertMono32ToMatrix(buf,height,width,stride)

%Description :	Converts a byte array of Mono32 data to a Matrix of pixel values
%
%Arguments	 :  buf - The image data in bytes to be converted 
%                       height - The height of the image data in pixels
%                       width - The width of the image data in pixels
%                       stride - The width of the image data in bytes
%
%Return		 :  ret - Check the help for return code meanings
%                       matrix - The converted matrix of pixel values
%
if nargin < 4
    stride = 0;
end

[ret,matrix] = obs.cam.sdk.andorsdk3functions('AT_ConvertMono32ToMatrix',buf,height,width,stride);