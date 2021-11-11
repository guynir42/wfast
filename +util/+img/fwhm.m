function val = fwhm(I, varargin)
% Usage: val = fwhm(I, varargin)
% Calculate the Full Width at Half-Maximum (FWHM) for a set of images. 
%
% INPUT: images in a 2D, 3D or 4D matrix. The function can run on multiple 
%        images, returning a value in a 2D matrix, where the 1st dimension 
%        is for the original 3th, and 2nd for the original 4th. 
% 
% NOTE: The images must be background subtracted for this to work! 
%
% OUTPUT: an estimate of the width. This is as close as we can get to the 
%         Full Width Half Maximum (FWHM) for each image. 
% 
% 
% Method "gaussian":
% Generate a filter bank with gaussians of increasing width, find the one
% that gives the highest S/N. 
% 
% Method "generalized gaussian":
% Same as gaussians, only using the generalized gaussian. See also the 
%
% OPTIONAL ARGUMENTS:
%   -method: which general method is used (default is "filters"):
%       *Method "radial": Uses the radialStats function to calculate the 
%        mean of each annulus, then finds the point where the value reaches 
%        half. Assumes image is a centered, radial distribution!
% 
%       *Method "moments": Call the photometry2 function, using an aperture 
%        that fits into the image and calculate the second moment, 
%        multiplying by 2.355. 
%       
%       *Method "filters": run a filter/template bank over the images, and 
%        choose the filter that gives the highest S/N. That filter's width
%        is used as the FWHM. The default is to use only gaussian kernels, 
%        but you can add also generalized gaussians by specifying the 
%        "generalized" parameter and a defocus annulus by specifying the 
%        "defocus" parameter. You can even add multiple values in each of 
%        those parameters, which makes even more templates. In either case, 
%        the best S/N for each image is chosen from all templates. 
%
%   -oversample: increase the resolution using sinc interpolation in Fourier
%                space. Default is [] which means no oversampling. 
%   -generalized: what power of the generalized gaussian to use. 
%                 Use bigger values to get flatter top, and sharper edges,
%                 This is ONLY used when setting method=filters. 
%                 Add multiple values to use more templates, added on top 
%                 of the gaussian templates. 
%   -defocus: when non empty, adds defocus annulus templates. 
%             The parameter value gives the sigma of the defocus_annulus. 
%             Only used when method=filters. 
%             Add multiple values to use more templates, added on top 
%             of the gaussian templates. 
%   -number_interp: minimal number of points along slope. If the number of
%                   pixels in the image is larger than this, we don't need
%                   to interpolate. Otherwise it will interpolate to this 
%                   number of points. Default is 100. Only for "radial". 
%   -aperture: set the aperture radius for photometry2. Default is to use 
%              half the smaller image size as radius. 
%              Only used when method=moments. 
%   -step_size: the step size to use when using different filters. The 
%               default is to use 0.05 pixels. 
%   -max_size: the maximal filter size, in pixels. Default is half of the 
%              smaller image size, in pixels. 
%   -min_size: the minimal filter size, in pixels. Default is 0.5 pixels. 
%   -fft: use FFT convolution instead of filter2. Default false. 
%
%

    import util.text.cs;

    if nargin==0, help('util.img.fwhm'); return; end

    if ndims(I)>4
        error('This function does not support arrays with more than 4 dimensions!'); 
    end
    
    input = util.text.InputVars;
    input.input_var('oversample', []); % before running the calculation, oversample using sinc interpolation
    input.input_var('method', 'filters'); % which method to choose: "gaussian", "generalized", "defocus", "radial", "moments"
    input.input_var('number_interp', 100); 
    input.input_var('aperture', []); % aperture to use for moments (default is what fits in the image size
    input.input_var('step_size', 0.05); % size steps for scanning multiple FWHM values
    input.input_var('max_size', []); % minimal FWHM value to scan (for filters only). Default is smallest image axis / 1.5
    input.input_var('min_size', 1); % maximal FWHM value to scan (for filters only)
    input.input_var('fft', false, 'use_fft'); % use FFT convolution instead of brute-force convolution (for filters only)
    input.input_var('generalized', [], 'power', 'generalized_power', 'generalized_gaussian_power'); % what power of the generalized gaussian to use
    input.input_var('defocus', [], 'defocus_sigma'); % The parameter gives the sigma of the defocus_annulus.  
    input.input_var('plot', false);
    input.input_var('ax', [], 'axis', 'axes');
    input.input_var('font_size', 18); 
    input.input_var('curves', 10); % how many curves to show when plotting
    input.scan_vars(varargin{:});
    
    % assume image is centered and padded to odd number of pixels
    
    S = size(I);
    S = S(1:2);
    
    I2 = reshape(I, [size(I,1), size(I,2).*size(I,3).*size(I,4)]); % flatten the array into 2D
    I2 = reshape(regionfill(I2, isnan(I2)), [S, size(I,3).*size(I,4)]); % reshape it back after removing NaNs, leave dim 3 and 4 together
    
    if isempty(input.max_size)
        input.max_size = ceil(min(S)/1.5); 
    end

    if input.oversample
        
        I2 = real(util.img.oversample(I2, input.oversample));
        
        S = size(I);
        S = S(1:2);
    
    else
        input.oversample = 1; % make calculations easier by setting this to a non-empty value
    end
    
    if cs(input.method, 'filters', 'templates') % all sorts of filter bank based calculations
        
        % g is a cell array of function handles, for each kind of template
        conversion = 2.55; % =2*sqrt(2*log(2));
        g{1} = @(w) util.shapes.gaussian(w./conversion, 'size', S, 'norm', 2); % gaussian normalized for matched filtering 
        l{1} = 'gaussians'; 
        
        for kk = 1:length(input.generalized)
            p = input.generalized(kk);
            conversion = 2*(p*log(2)).^(1/p); % replace 2.355 with the correct factor
            g{end+1} = @(w) util.shapes.generalized_gaussian('sigma_x', w./2.355, 'size', S, 'norm', 2, 'power', p);  
            l{end+1} = sprintf('generalized= %4.2f', input.generalized(kk)); 
        end
        
        for kk = 1:length(input.defocus)
            g{end+1} = @(w) defocus_func(w, input.defocus(kk), S); % use helper function to cover both annulus and small gaussian cases
            l{end+1} = sprintf('defocus sig= %4.2f', input.defocus(kk)); 
        end
        
        width = input.min_size:input.step_size:input.max_size; % this is in units of FWHM
        width = width.*input.oversample; % if the image has been enlarged, the scan should be enlarged as well
        snr = nan(size(width,2), length(g), size(I2,3)); % the results: maximal S/N for each width, each function, and each cutout 
        
        for kk = 1:length(g)

            for ii = 1:length(width)

                if input.fft
                    If = util.img.conv_f(g{kk}(width(ii)), I2); 
                else
                    I3 = reshape(I2, [size(I2,1), size(I2,2).*size(I2,3).*size(I2,4)]); % shape the 4D array into a long 2D matrix
                    If2 = filter2(g{kk}(width(ii)), I3);
                    If = reshape(If2, size(I2)); % reshape back to original size after filtering
                end

                snr(ii,kk,:) = util.stat.max2(If); % the maximum in each cutout (for each width and each kernel)

            end % for ii (widths)

        end % for kk (kernels)
        
        [~, idx] = util.stat.max2(snr); % find the maximum sigma for each ii
        
        val = util.vec.tocolumn(width(idx(1,1,:))); 
        
    elseif cs(input.method, 'radial')

        r_max = max(S);

        for ii = 1:size(I2,3)
            
            [values, radii] = util.img.radialStats(I2(:,:,ii), 'mean', 'r_min', 0, 'r_max', r_max, 'r_delta', 1); 
            
            N = max(length(radii), input.number_interp); 
            
            if ii==1
                radii_interp = zeros(size(I2,3),N);
                values_interp = zeros(size(I2,3),N);
            end
            
            if length(radii)<input.number_interp
                radii_interp(ii,:) = linspace(0, r_max, input.number_interp);
                values_interp(ii,:) = interp1(radii, values, radii_interp(ii,:));
            else
                radii_interp(ii,:) = radii;
                values_interp(ii,:) = values;
            end

            idx(ii) = find(values_interp(ii,:)<max(values_interp(ii,:))/2, 1, 'first');

            if isempty(idx(ii))
                val(ii) = NaN;
            else
                val(ii) = radii_interp(ii,idx(ii));
            end

        end

        val = val*2; % convert half-width at half-maximum to FWHM

    elseif cs(input.method, 'moments')
        
        if isempty(input.aperture)
            input.aperture = ceil(min(S)/2); 
        end
        
        % to be continued... 
        
    else
        error('Unknown "method" option "%s". Choose "gaussian", "generalized", "defocus", "radial" or "moments".', input.method); 
    end
    
    val = reshape(val, [size(I,3),size(I,4)]); 
        
    val = val./input.oversample;
    
    if input.plot
        
        if isempty(input.ax)
            input.ax = gca;
        end
        
        N = min(size(snr,3), input.curves); % don't draw more lines than given by "curves"
        
        if cs(input.method, 'filters', 'templates')
            
            cla(input.ax);
            hold(input.ax, 'on'); 
            highest_point = nanmax(snr(:));
            
            r = width./input.oversample;
            
            line_types = {'-', '--', ':'}; 
            marker_types = {'x', 's', 'v', 'p', 'v'}; 
            
            for ii = 1:N
                
                for kk = 1:size(snr,2)
            
                    mark = marker_types{mod(kk-1,length(marker_types))+1};
                    line = line_types{ceil(kk./length(marker_types))};
                    color = input.ax.ColorOrder(ii,:); 
                    
                    s = snr(:,kk,ii);
                    h = plot(input.ax, r, s, [mark line], 'Color', color); 
                    if ii>1
                        h.HandleVisibility = 'off'; 
                    end
                    
                    [mx,idx] = nanmax(s); 
                    
                    if mx==util.stat.max2(snr(:,:,ii))
                        plot(input.ax, r(idx), s(idx), 'o', 'Color', color, 'MarkerSize', 10, 'HandleVisibility', 'off');
                        text(input.ax, r(idx), s(idx)+0.05*highest_point, sprintf('%4.2f', r(idx)), ...
                            'Color', color, 'HorizontalAlignment', 'center', 'FontSize', input.font_size-2); 
                    end
                    
                end % for kk (kernels)
                
            end % for ii (stars)
            
            xlabel(input.ax, 'FWHM [pixels]');
            ylabel(input.ax, 'S/N'); 
            
            legend(input.ax, l, 'Location', 'Best'); 
            
            hold(input.ax, 'off');
            
        elseif cs(input.method, 'radial')
            
            hold(input.ax, 'off');
            
            h = plot(input.ax, radii_interp(1:N,:)'./input.oversample, values_interp(1:N,:)'); 
            
            hold(input.ax, 'on'); 
            
            for ii = 1:length(h)
                r = radii_interp(ii,idx(ii))./input.oversample;
                v = values_interp(ii,idx(ii));
                plot(input.ax, r, v, 'o', 'Color', h(ii).Color);
                plot(input.ax, [1,1].*r, [0, v], '--', 'Color', h(ii).Color);
                text(input.ax, r+0.2, v, sprintf('%4.2f', r), ...
                    'Color', h(ii).Color, 'HorizontalAlignment', 'left', 'FontSize', input.font_size-2); 
            end
            
            xlabel(input.ax, 'Radius [pixels]');
            ylabel(input.ax, 'Average light [counts]');
            
            hold(input.ax, 'off');
            
        end
        
        input.ax.FontSize = input.font_size; 
        box(input.ax, 'on'); 
        
    end
    
end

function g = defocus_func(w, sig, S)
    
    r = w/2 - sig*2.355/2; % note the FWHM parameter "width" is equal to 2*r + 2.355*sigma
    
    if r<=0 % if the annulus is so small as to not have a hole anymore, turn it into a small gaussian
        r = 0;
        sig = w./2.355; % gaussian sigma with the regular conversion from FWHM
    end
    
    g = util.shapes.defocus_annulus('r', r, 'sigma', sig, 'size', S, 'norm', 2); 
    
end
    




