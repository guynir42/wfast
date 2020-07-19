function data_out = wiener(data, aux, varargin)
% Usage: data_out = wiener(data, aux, varargin)
% Clean up the time series "data" using the auxiliary measurements "aux". 
    
    if nargin==0, help('util.series.wiener'); return; end
    
    if nargin<2 || isempty(aux)
        error('Must supply a valid time-series and a valid auxiliary measurement'); 
    end
    
    if size(data,1)~=size(aux,1)
        error('Size mismatch! size(data)= %s | size(aux)= %s', util.text.print_vec(size(data), 'x'), util.text.print_vec(size(aux), 'x'));
    end
    
    input = util.text.InputVars;
    input.input_var('window', 256); 
    input.input_var('shape', 'welch'); 
    input.input_var('overlap', 0.5); 
    input.input_var('median', true, 'use_median');
    input.input_var('blue', false); % use "double whitening" (or "bluing") to prepare for also doing a matched-filter later
    input.scan_vars(varargin{:}); 
    
    % start by calculating the power spectrum (ps) for data and aux, and also the covariance between them in Fourier space
    [ps, ps_aux, covar] = util.series.welch_ps(data, 'covariance', aux, ...
        'window', input.window, 'shape', input.shape, 'overlap', input.overlap, 'median', input.median');
    
    data_f = fft(data);
    aux_f = fft(aux); 
    
    numer = ps_aux.*data_f - conj(covar).*aux_f;
    denom = ps.*ps_aux - abs(covar).^2; 
    
    if ~input.blue
        denom = sqrt(denom); 
    end
    
    data_filtered = numer./denom; 
    
%     data_filtered(isinf(data_filtered)) = 0; 
    
    data_out = ifft(data_filtered); 
    
end