% Usage: pos = find_cosmic_rays(images, mask, threshold=256, threads=1, max_number=100, debug_bit=0)
% Finds hot pixels (above the threshold) in 2D or 3D uint16 raw images. 
% Especially useful for detecting cosmic rays and other short lived flares
% in the data. Ideally the positions should be used to make cutouts around
% each peak, for saving the full frame rate data around those pixels. 
%
% The second argument is a logical mask of the same size as a single image. 
% This should include true wherever there is a bad pixel and also if there
% are any known sources like stars. 
%
% The minimal threshold is 256, which lets us search only the second byte
% in each pixel value. The max is 65535 which is full saturation. 
% 
% The "threads" argument can be used to split the images into multiple 
% CPU threads. The number of threads must be lower than the number of 
% physical cores and also should be lower than the number of images in the
% data cube. 
%
% The last input gives the maximum number of peaks PER THREAD. This helps 
% prevent exceeding the allowed runtime when an image has exceptionally many 
% peaks. 
%
% The output is a matrix with four columns, containing the x,y,z position 
% and value of each peak in the data cube (one peak in each row).
% The position is for the brightest pixel in the nearby region (3 pixels in 
% each direction). 
% 
% If you want to rearrange the peaks in order of brightest to faintest,
% just do sortrows(pos, 4, 'descend'). 
%