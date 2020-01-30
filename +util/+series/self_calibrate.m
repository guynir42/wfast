function [flux_cal] = self_calibrate(flux, varargin)
% Usage: [flux_cal] = self_calibrate(flux, varargin)
% Take fluxes from one or several sources from the same epochs and remove 
% large trends, red noise, zero point correction, etc. 
%
% Input: a 2D matrix of fluxes, dim 1 is time, dim 2 is stars. 
%         
% Outputs: a matrix of the same size, with trends removed. 
%
% OPTIONAL ARGUMENTS: (the first 2 arguments can be given without keywords)
%   -timestamps: give the epoch timing (default is just 1:num_epochs). 
%   -errors: give the measured rms for each epoch and for each star. 
%   -use_zp: use a different zero-point for each epoch. Only for multiple stars. 
%            Default is false. 
%   -iterations: how many iterations on zero point removal to do (default 1). 
%                The 1st iteration just uses each star's flux as a measure
%                of the weight it should be given. The next iterations use
%                the cleaned-up flux variance across the whole run to get
%                the relative weight between stars. 
%   -use_poly: use an Nth order polynomial to fit simultaneously to all the
%              stars' lightcurves. Each iteration uses the mean flux of each
%              lightcurve to match them to the same level and then average
%              out the fit coeffs. The average polynomial is subtracted, then
%              another iteration improves the fit. Default true.
%   -order: highest order polynomial for the fit (default 5). 
%   -use_welch: apply Welch's method to calculating the power spectral density, 
%               then use it to rescale the noise at different frequencies 
%               (make it white noise). This can be applied to a single star's 
%               flux, or to multiple stars individually. Default false. 
%   -window: number of samples in the welch window. 
% 

    if nargin==0, help('util.series.self_calibrate'); return; end
    
    input = util.text.InputVars;
    input.use_ordered_numeric = 1;
    input.input_var('timestamps', []); 
    input.input_var('errors', []); 
    input.input_var('use_zp', false, 'zero point'); 
    input.input_var('iterations', 1); 
    input.input_var('use_poly', true, 'use_polynomial', 'polynomial'); 
    input.input_var('order', 5); 
    input.input_var('use_welch', false, 'welch'); 
    input.input_var('window', 128); 
    input.scan_vars(varargin{:}); 
    
    if isempty(input.timestamps)
        input.timestamps = 1:size(flux,1); 
    end
    
    f = flux; % temporary variable to be adjusted in each iteration
    t = input.timestamps; % shorthand for the times
    
    if input.use_zp
        
        for ii = 1:input.iterations
            
        end
        
    end
    
    if input.use_poly
        
        flog = abs(log10(f)); % abs gets rid of any complex numbers brought on by negative fluxes (there may be better ways to do this, maybe use Luptitudes?)  
        
        r = util.fit.polyfit(t, flog, 'order', input.order, 'iterations', 3, 'sigma', 3, 'double', 1); % fit each log-flux independently 

        C = [r.coeffs]; 
        M = 10.^C(1,:); % overall scaling of each star
        C = sum(M.*C,2)./sum(M,2); % get the (flux weighted) average polynomial coeffs
        C(1) = 0; % do not try to make any changes to the overall scaling of each flux
        C = C'; 

        T = ones(size(t,1),1); % construct the time axis in different powers

        for jj = 1:input.order
            T = [T t.^jj]; 
        end

        model = sum(C.*T,2); % this polynomial model represents the best fit to all fluxes

        flog = flog - model; 
       
        f = 10.^flog;
        
    end
        
    if input.use_welch
        % to be continued...
    end
    
    flux_cal = f; % make sure to output something even if all methods are false! 
        
end