function data_out = remove_aux(data, aux, varargin)
% Usage: data_out = remove_aux(data, aux, varargin)
% Remove the contribution of the auxiliary "aux" from the time series "data". 
% The two series must have the same length. "data" can be any dimension, 
% while "aux" should be a vector. 
% (support for multiple auxiliaries may be added later)
%
% Will remove the auxiliary by calculating the power-spectrum for data, aux 
% and their covariance, using the util.series.welch_ps() function.
%
% Optional Arguments are passed to the welch_ps function, including:
% window, shape, overlap, median. 
% Additional optional arguements include:
%

    if nargin==0, help('util.series.remove_aux.m'); return; end
    
    if isempty(data), data_out = []; return; end
    
    if isrow(aux)
        aux = aux';
    elseif size(aux,1)>1 && size(aux,2)>2
        error('Input "aux" must be a vector!'); 
    elseif ~ismatrix(aux)
        error('Input "aux" must be a vector!'); 
    end
    
    if nargin<2 || isempty(aux) || size(data,1)~=size(aux,1)
        error('Size mismatch between "data" (%s) and "aux" (%s).', util.text.print_vec(size(data), 'x'), util.text.print_vec(aux, 'x')); 
    end

    input = util.text.InputVars;
    input.input_var('window', 256); 
    input.input_var('shape', 'welch'); 
    input.input_var('overlap', 0.5); 
    input.input_var('median', true, 'use_median');
    input.input_var('blue', false); % use "double whitening" (or "bluing") to prepare for also doing a matched-filter later
    input.scan_vars(varargin{:}); 
    
    data = fillmissing(data, 'linear'); 
    M_data = mean(data, 1); 
    
    data = data - M_data;
    
    % start by calculating the power spectrum (ps) for data and aux, and also the covariance between them in Fourier space
    [ps, ps_aux, covar] = util.series.welch_ps(data, 'covariance', aux, ...
        'window', input.window, 'shape', input.shape, 'overlap', input.overlap, 'median', input.median);
    
    data_f = fft(data);
    aux_f = fft(aux); 
    
    data_filtered = data_f - conj(covar)./ps_aux.*aux_f;
    
    data_out = real(ifft(data_filtered));

    data_out = data_out + M_data; 
    
end

    
    