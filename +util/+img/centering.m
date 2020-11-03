function I_out = centering(I, varargin)
% Usage: I_out = centering(I, varargin)
% Move each star in the given images/cutouts "I" to be centered at the
% middle of the frame. 
% 
% OPTIONAL PARAMETERS:
%   -peak: how to find the peak. Default is "max", find the brightest
%          pixel. Alternative is "centroid". 
%   -photometery_parameters: pass extra parameters, as varargin cell array
%           with key-value pairs. These overwrite the defaults used for the
%           photometry function. This is only used if choosing
%           peak=centroids. 
%   -shift: how to move the data. Default is circshift, with integer pixel
%           moves. Alternative is "fft" for fractional moves with the
%           fft-shift method (util.img.fftshift). 

    import util.text.cs;
    
    if nargin==0, help('util.img.centering'); return; end 

    input = util.text.InputVars;
    input.input_var('peak', 'max'); % can also choose "centroid"
    input.input_var('photometry_parameters', {}); % a cell array with additional (overwriting) values for photometery2 (in case we use "centroids" method)
    input.input_var('shift', 'circ'); % can also choose "fft" for fractional moves
    input.scan_vars(varargin{:}); 

    if ndims(I)>4
        error('Cannot handle images with more than 4 dimensions... size(I)= %s.', util.text.print_vec(size(I), 'x')); 
    end

    if cs(input.peak, 'maximum')
        [~, idx] = util.stat.max2(I); 
        
        cen = floor(size(I)/2) + 1; % image center
        cen = cen(1:2); % keep only 2 first dimensions
        
        offset = cen-idx; 
        
    elseif cs(input.peak, 'centroids')
        s = util.img.photometry2(I, 'use_aperture', 0, 'use_centering', 1, ...
            'aperture', 5, 'use_gaussian', 1, 'gauss', 3, 'use_forced', 0, ...
            'index', 4); 
        offset = -[permute(s.gaussian_photometry.offset_x, [3,4,1,2]), ...
                permute(s.gaussian_photometry.offset_y, [3,4,1,2])];
    else
        error('Unknown "peak" option "%s". Use "max" or "centroid". ', input.peak); 
    end
    
    I_out = zeros(size(I), 'like', I); 
    
    if cs(input.shift, 'circshift')
    
        for ii = 1:size(I,3)
            for jj = 1:size(I,4)
                I_out(:,:,ii,jj) = circshift(I(:,:,ii,jj), round(offset(:,:,ii,jj))); 
            end
        end
    
    elseif cs(input.shift, 'fftshift')
        for ii = 1:size(I,3)
            for jj = 1:size(I,4)
                I_out(:,:,ii,jj) = util.img.FourierShift2D(I(:,:,ii,jj), offset(:,:,ii,jj)); 
            end
        end
    else
        error('Unknown "shift" option "%s". Use "circ" or "fft"', input.shift); 
    end
    
end




