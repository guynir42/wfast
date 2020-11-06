function val = fwhm(I, varargin)
% Usage: val = fwhm(I, varargin)
% Calculate the Full Width at Half-Maximum (FWHM) for an image.
%
% Method "radial":
% Uses the radialStats function to calculate the mean of each annulus, 
% then finds the point where the value reaches half. 
% Assumes image is a centered, radial distribution!
% 
% Method "moments":
% Call the photometry2 function, using an aperture that fits into the image
% and calculate the second moment, multiplying by 2.355. 
%
% Method "filters":
% Generate a filter bank with gaussians of increasing width, find the one
% that gives the highest S/N. 
%
% OPTIONAL ARGUMENTS:
%   -oversample: increase the resolution using sinc interpolation in Fourier
%                space. Default is [] which means no oversampling. 
%   -method: choose the different methods to calculate FWHM. 
%            Options: "radial", "moments", "filters". 
%   -gaussian: what power of the generalized gaussian to use (default=2 for 
%              normal gaussian). Use bigger values to get flatter top, like
%              when the telescope is de-focused. This is ONLY used when
%              setting method=filters! 
%   -number_interp: minimal number of points along slope. If the number of
%                   pixels in the image is larger than this, we don't need
%                   to interpolate. Otherwise it will interpolate to this 
%                   number of points. Default is 100. Only for "radial". 
%   -aperture: set the aperture radius for photometry2. Default is to use 
%              half the smaller image size as radius. Only for "moments". 
%   -step_size: the step size to use when using different filters. The 
%               default is to use 0.05 pixels. Only for "filters".
%   -max_size: the maximal filter size, in pixels. Default is half of the 
%              smaller image size, in pixels. Only for "filters". 
%   -min_size: the minimal filter size, in pixels. Default is 0.5 pixels. 
%              Only for "filters". 
%   -fft: use FFT convolution instead of filter2. Default false. 
%   NOTE: the parameters used for "filters" are the gaussian width, not the 
%         FWHM, so the end result (the best filter) is multiplied by 2.355. 
%

    import util.text.cs;

    if nargin==0, help('util.img.fwhm'); return; end

    if ndims(I)>4
        error('This function does not support arrays with more than 4 dimensions!'); 
    end
    
    input = util.text.InputVars;
    input.input_var('oversample', []); 
    input.input_var('method', 'radial'); 
    input.input_var('number_interp', 100);
    input.input_var('aperture', []); 
    input.input_var('step_size', 0.05);
    input.input_var('max_size', []); 
    input.input_var('min_size', 1);
    input.input_var('fft', false, 'use_fft'); 
    input.input_var('gaussian', 2); % what power of the generalized gaussian to use (default=2 for normal gaussian)
    input.scan_vars(varargin{:});
    
    % assume image is centered and padded to odd number of pixels
    
    S = size(I);
    S = S(1:2);
    
    if cs(input.method, 'radial')

        r_max = max(S);

        for ii = 1:size(I,3)
            for jj = 1:size(I,4)

                [values, radii] = util.img.radialStats(I, 'mean', 'r_min', 0, 'r_max', r_max, 'r_delta', 1); 
                if length(radii)<input.number_interp
                    radii_interp = linspace(0, r_max, input.number_interp);
                    values_interp = interp1(radii, values, radii_interp);
                else
                    radii_interp = radii;
                    values_interp = values;
                end

                idx = find(values_interp<max(values_interp)/2, 1, 'first');

                if isempty(idx)
                    val(ii,jj) = NaN;
                else
                    val(ii,jj) = radii_interp(idx);
                end

    %             plot(radii, values, '-', radii_interp, values_interp, ':');

            end
        end

        val = val*2; % convert half-width at half-maximum to FWHM

    elseif cs(input.method, 'moments')
        
        if isempty(input.aperture)
            input.aperture = ceil(min(S)/2); 
        end
        
        % to be continued... 
        
    elseif cs(input.method, 'filters')
        
        if isempty(input.max_size)
            input.max_size = ceil(min(S)/1.5); 
        end
        
%         I2 = I-util.stat.mean2(I); 
        
        I2 = reshape(I, [size(I,1), size(I,2).*size(I,3).*size(I,4)]); % flatten the array into 2D
        I2 = reshape(regionfill(I2, isnan(I2)), size(I)); % reshape it back after removing NaNs
        
        sig = input.min_size:input.step_size:input.max_size;
        mx = nan(size(I,3), size(I,4), size(sig,2)); 
        
        for ii = 1:length(sig)
            
%             g = util.img.gaussian2(sig(ii)./2.355, 'size', S, 'norm', 2); % gaussian normalized for matched filtering 
            g = util.img.generalized_gaussian('sigma_x', sig(ii)./2.355, 'size', S, 'norm', 2, 'power', input.gaussian);  
            
            if input.fft
                If = util.img.conv_f(g, I2); 
            else
                If = reshape(filter2(g, reshape(I2, [size(I,1), size(I,2).*size(I,3).*size(I,4)])), size(I)); % reshape back to original size
            end
            
            mx(:,:,ii) = permute(util.stat.max2(If), [3,4,1,2]); % the maximum in each cutout for each sigma
            
        end
        
        [~, idx] = nanmax(mx,[], 3); % find the maximum sigma for each ii and jj
        
        val = sig(idx); 
        
    else
        error('Unknown method "%s". Choose "radial" or "moments" or "filters". ', input.method); 
    end
    
    
    
end





