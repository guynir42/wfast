function [ps, ps_cov, covariance] = welch_ps(data, varargin)
% Usage: [ps, ps_cov, covariance] = welch_ps(data, varargin)
% Calculate the Power Spectrum (PS) of the input data, using overlapping
% windows and averaging the PS of each window. 
% 
% Optional arguments:
% 
% 

    import util.text.cs;
    
    if nargin==0, help('util.series.welch_ps'); return; end
    
    input = util.text.InputVars;
    input.input_var('window', 256); 
    input.input_var('shape', 'welch'); 
    input.input_var('overlap', 0.5); 
    input.input_var('median', true, 'use_median'); 
    input.input_var('points', [], 'num_points'); 
    input.input_var('covariance', []); % estimate the covariance of "data" against this other time series
    input.scan_vars(varargin{:}); 
    
    if isempty(input.overlap)
        input.overlap = 0.5;
    end
    
    if isempty(input.points)
        input.points = size(data,1);
    end
    
    S = size(data); % keep track of original size of data
    
    data = data(:,:); % linearize any higher dimensions
%     data = permute(data, [1,3,2]); % push dimensions from 2 to 3
    
    % make the window mask
    step = floor(input.window.*(1-input.overlap));
    
    N = floor(size(data,1)./step) - 1;
    if N<1
        error('Window size %d is larger than length of data %d...', input.window, S(1));
    end
    
    if cs(input.shape, 'squared', 'boxcar', 'dirichlet')
        window_shape = ones(input.window, 1, 'like', data); 
    elseif cs(input.shape, 'triangle', 'bartlett')
        window_shape = ones(input.window, 1, 'like', data); 
        window_shape = util.vec.convolution(window_shape, window_shape); 
    elseif cs(input.shape, 'parabolic', 'welch')
        n = (1:input.window)';
        window_shape = 1 - (n-floor(input.window/2)).^2./floor(input.window/2).^2;        
    else
        error('Unknown "shape" parameter "%s". Use "square" or "triangle" or "welch" (default).', input.shape); 
    end
    
%     norm = sum(window_shape.^2); % divide by S(1)??
%     window_shape = window_shape ./ sum(window_shape.^2);
    data_shifted = zeros(length(window_shape), size(data,2), N, 'like', data); 
    
    idx = 1;
    for ii = 1:N
        end_idx = idx+input.window-1;
        data_shifted(:,:,ii) = data(idx:end_idx,:); 
        idx = idx + step;
    end
    
    data_windowed = data_shifted.*window_shape;
    
%     mask = zeros(S(1), N, 'like', data); 
%     
%     idx = 1;
%     
%     for ii = 1:N
%         
%         mask(idx:idx+input.window-1,ii) = window_shape; 
%         
%         idx = idx + step; 
%         
%     end
%     
%     mask = mask(1:S(1),:); 
%     
%     data_windowed = data.*mask; 
    
    data_f = fft(data_windowed, input.points); 
    
    if input.median
        ps = nanmedian(abs(data_f).^2, 3); 
    else
        ps = nanmean(abs(data_f).^2, 3); 
    end
    
    S(1) = input.points;
    ps = reshape(ps, S) * N; 
    
%     if ~isempty(input.covariance)
%         
%         if iscolumn(input.covariance)==0
%             error('We currently only support 1D covariance/auxiliary input!'); 
%         end
%         
%         cov_windowed = input.covariance.*mask; 
%         cov_f = fft(cov_windowed); 
%         
%         if input.median
%             ps_cov = nanmedian(abs(cov_f).^2,2)./norm; 
%         else
%             ps_cov = nanmean(abs(cov_f).^2,2)./norm; 
%         end
% 
%         cov_sections = conj(cov_f).*data_f;
%         
%         % no median: it isn't defined for complex numbers...
%         covariance = nanmean(cov_sections,2)./norm;
%         
%         covariance = reshape(covariance, S); 
%         
%     end
    
end