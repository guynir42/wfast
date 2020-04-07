function f_out = binning(f, varargin)
% Usage: f_out = binning(f, bin_factor=2, varargin)
% Take a multi dimensional matrix of time series (with time along 1st dim)
% and bin them by some factor (default 2). 
%
% Optional Arguments:
%   -factor: How many data points to use for each bin (default: 2). 
%   -function: What function to use to combine each bin (default: mean).
%              Can use any of these: mean, median, sum, std, var, min, max. 
%              All functions used will ignore NaNs. 
%   -residual: Do we want to use the last few data points left over after
%              dividing into full bins? 
%   -repmat: Use repmat to stretch the binned data back to the original size,  
%            e.g., for subtracting from the original lightcurve. 
%            Note: using "repmat" automatically also enables "residual". 
%
% NOTE: To get the correct timestamps you have two options:
%       (a) If the timestamps measure the beginning of the measurement then
%           you need to do t2 = t(1:bin_factor:end). That's easy. 
%       (b) If the timestamps are for the middle of the exposure then you
%           need to call this function again: t2 = binning(t, bin_factor). 

    import util.text.cs;

    if nargin==0, help('util.series.binning'); return; end
    
    input = util.text.InputVars;
    input.use_ordered_numeric = 1;
    input.input_var('factor', 2, 'bin_factor', 'binning'); 
    input.input_var('func', 'mean', 'function'); 
    input.input_var('residual', false); % keep the non-dividable residual samples at the end?
    input.input_var('repmat', false); 
    input.scan_vars(varargin{:}); 

    if isempty(input.factor)
        input.factor = 2;
    elseif ~isscalar(input.factor)
        error('Must input an integer scalar binning factor!');
    end
    
    if input.factor<=1 && ~cs(input.func, 'std', 'var')
        f_out = f; 
        return; % short circuit in case factor==1
    elseif input.factor<=1 && cs(input.func, 'std', 'var')
        error('Cannot calculate var/std on binning factor less than 2!'); 
    elseif input.factor~=floor(input.factor)
        error('Must use integer binning factor!');
    end
    
    if isvector(f) && size(f,1)==1
        f = f'; % make sure f is a column vector! 
    end
    
    S = size(f); % keep track of the original size of the input
    N = floor(size(f,1)./input.factor).*input.factor; % the number of samples that divide cleanly by input.factor
    
    if input.residual || input.repmat % add Nans
        
        f = f(:,:); % linearize all higher dimensions... 
        
        if N~=size(f,1)
            N = N+input.factor; % adjust N to the new size of the corrected dataset
            f = vertcat(f,NaN(N - size(f,1),size(f,2)));
        end
        
    else
        f = f(1:N,:); % truncate the residuals and linearize all higher dimensions... 
    end
    
    f = reshape(f, [input.factor size(f,1)./input.factor size(f,2)]);
    
    if cs(input.func, 'mean')
        f = nanmean(f,1); 
    elseif cs(input.func, 'median')
        f = nanmedian(f,1);
    elseif cs(input.func, 'sum')
        f = nansum(f,1);
    elseif cs(input.func, 'std')
        f = nanstd(f,[],1);
    elseif cs(input.func, 'var')
       f = nanvar(f, [], 1);
    elseif cs(input.func, 'min')
        f = nanmin(f, [], 1);
    elseif cs(input.func, 'max')
        f = nanmax(f, [], 1);
    else
        error('Unknown function option: "%s". Use "mean", "median", "sum", "std", "var", "min" or "max"...', input.func); 
    end
    
    if input.repmat
        f = repmat(f, [input.factor, 1, 1]); 
        f_out = reshape(f, [N, S(2:end)]);
        f_out = f_out(1:S(1), :); 
    else
        f = permute(f, [2,3,1]); % get rid of dim 1 which we just integrated on
        f_out = reshape(f, [size(f,1), S(2:end)]);
    end
    
    
    
    
    
    
    