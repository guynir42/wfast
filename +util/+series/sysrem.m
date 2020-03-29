function [f_det, a, c, stats] = sysrem(f, varargin)
% Usage: [f_det, a, c, stats] = sysrem(f, varargin)
% Apply the Tamuz et al. 2005 Sys-Rem detrending on fluxes 'f'. 
% Assumes the time axis is the first dimension of 'f' and the different 
% stars are on the second dimension. 
%
% Iteratively find the best a_j ("airmass") for each frame and best c_i 
% (coefficient) for each star, that minimize the rms of r_ij-a_j*c_i
% (r_ij are the residuals after subtracting the mean of f for each frame j 
% and each star i). 
% 
% OUTPUTS: -'f_det' is the detrended fluxes (same size as f). 
%          -'a' and 'c' are the coefficients that were subtracted. 
%          -'stats' contains statistics on the Sys-Rem process. 
%           (currently this is empty, we need to finish that)
% 
% OPTIONAL ARGUMENTS:
% 


    if nargin==0, help('util.series.sysrem'); return; end
    
    input = util.text.InputVars;
    input.use_ordered_numeric = 1;
    input.input_var('errors', [], 'sigma'); % error on each flux measurement, (use vector/scalar to mark errors per frame/per star/global error)
    input.input_var('a', [], 'airmass'); % already have a good initial guess for the 'a' vector
    input.input_var('c', [], 'coeffs', 'coefficients'); % already have a good inital guess for the 'c' vector\
    input.input_var('iterations', 0, 'num_iterations'); % how many times to 'criss-cross' the coefficients 'a' and 'c' after the inital guess
    input.input_var('magnitudes', false); % assume input is given as magnitudes, with negative values taken care of
    input.input_var('pedestal', true, 'use_pedestal', 'constant'); % add some constant to each flux measurement before turning it into magnitudes (instead of putting NaN in negative values)
    input.input_var('subtract', true, 'subtract_mean', 'use_subtract_mean', 'use_mean', 'mean'); % if true, will remove the mean of each flux f_ij to get residual r_ij
    input.input_var('stars', [], 'star_black_list'); % which stars should be avoided when calculating 'a'
    input.scan_vars(varargin{:}); 
    
    
    if input.magnitudes
        m = f;
    else
        if input.pedestal
            P = util.stat.min2(f)*2; 
            if P<0
                P = -P;
            else
                P = 0;
            end
            m = -2.5.*log10(f+P);         
        else
            P = 0;
            f(f<=0) = NaN;
            m = -2.5.*log10(f); 
%             m = asinh(f); 
        end
        
        if ~isempty(input.errors) % errors are given on the flux, not on the magnitudes! 
            input.errors = input.errors./abs(nanmean(f,1)); % translate flux errors to mag errors (approximation of 2.5 to e)
        end
        
    end
    
    if isempty(input.errors)
        input.errors = ones(size(m), 'like', m); % uniform errors! 
    end
    
    if input.subtract
        M = nanmean(m,1); 
        r = m - M;
    else
        M = 0;
        r = m;
    end
    
    if isempty(input.stars)
        input.stars = false(1, size(f,2)); 
    end
    
    input.errors(input.errors==0) = Inf;
    
%     if size(input.errors,1)==1
%         input.errors = repmat(input.errors, [size(f,1),1]); 
%     end
%     
%     if size(input.errors,2)==1
%         input.errors = repmat(input.errors, [1, size(f,2)]); 
%     end
    
    if isempty(input.a) && isempty(input.c) % no 'a' or 'c' given, do a preliminary calculation to get both
        c = ones(1, size(f,2), 'like', f); % assume all stars have the same response coefficient! 
        a = nansum(r(:,~input.stars).*c(~input.stars)./input.errors(:,~input.stars).^2, 2)./nansum(c(~input.stars).^2./input.errors(:,(~input.stars)).^2, 2); % sum over 'i' is sum over stars
        c = nansum(r.*a./input.errors.^2, 1)./nansum(a.^2./input.errors.^2, 1); % sum over 'j' is sum over frames
    elseif isempty(input.a) % only 'a' is empty, use 'c' to find 'a'
        c = input.c;
        a = nansum(r(:,~input.stars).*c(~input.stars)./input.errors(:,~input.stars).^2, 2)./nansum(c(~input.stars).^2./input.errors(:,(~input.stars)).^2, 2); % sum over 'i' is sum over stars
    elseif isempty(input.c) % only 'c' is empty, use 'a' to find 'c' (this can happen if 'a' is given by airmass)
        a = input.a;
        c = nansum(r.*a./input.errors.^2, 1)./nansum(a.^2./input.errors.^2, 1); % sum over 'j' is sum over frames
    else
        a = input.a;
        c = input.c;
    end
    
    for ii = 1:input.iterations
        a = nansum(r(:,~input.stars).*c(~input.stars)./input.errors(:,~input.stars).^2, 2)./nansum(c(~input.stars).^2./input.errors(:,(~input.stars)).^2, 2); % sum over 'i' is sum over stars
        c = nansum(r.*a./input.errors.^2, 1)./nansum(a.^2./input.errors.^2, 1); % sum over 'j' is sum over frames
    end
    
    r_new = r - a.*c;
    
    m_det = r_new + M; 
    
    f_det = 10.^(-0.4.*m_det) - P;
%     f_det = sinh(m_det);
    stats = [];
    
end

