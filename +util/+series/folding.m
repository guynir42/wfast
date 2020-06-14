function [f_folded, var_folded] = folding(f, period, varargin)
% Usage: [f_folded, var_folded] = folding(f, period, varargin)
% Fold the flux data f along a period with "period" bins. 
% The data must be aligned with time on the first dimension, other dimensions
% can be anything else. 
% The flux should be evenly sampled, with any gaps replaced by NaNs. 
% It is usually also a good idea to remove slow trends and to subtract the 
% average of each flux before folding. 
%
% Optional arguments:
%   -method: "fast" (default) assumes uniformly sampled data, and just folds
%            on the "period" by manipulating the matrix dimensions and
%            then summing ther aligned bins. 
%            In the "slow" mode, the "period" input is used to find the right
%            bin for each flux value, using the mod() function on the times. 
%            If the timestamps are not given, they are replaced by linsapce
%            with interval of 1 between each measurement. 
%
%   -timestamps: give explicite timestamps, and give "period" in the same
%                units. This is only used in the "slow" method. 
%
%   -variances: apply the same folding to the variance data. This data can 
%               be used for estimating the noise in each folded bin. 
%
%   -median: use median instead of mean of each phase bin.
%
%   -oversample: use util.series.resample() to increase the sampling rate
%                before folding. This helps for using periods that are more 
%                finely separated than the frame rate. If you really want
%                accuracy, though, you should just use the "slow" mode. 
%
%
% Output: The folded flux, the mean (or median) in each phase bin.  
%         If an additional output argument is given, will output the var in 
%         each phase bin as well. 


    if nargin==0, help('util.series.folding'); return; end
    
    if period==1 || period>size(f,1)
        error('Must give "period" to be larger than 1 and smaller than the length of f');
    end
    
    period = round(period);
    
    input = util.text.InputVars;
    input.use_ordered_numeric = 1; 
    input.input_var('timestamps', []); % can give time stamps for non-equally spaced sampling (assume timestamps are for start-of-acquistion)
    input.input_var('variances', []); % also keep track of the variances
    input.input_var('method', 'fast'); % can choose "fast" (default) or "slow" (go over a range of periods explicitely)
    input.input_var('median', false, 'use_median');
    input.input_var('oversample', []); 
    input.scan_vars(varargin{:});
    
    if input.oversample
        f = util.series.resample(f, input.oversample); 
        period = period*input.oversample;
    end
    
    S = size(f); % keep track of the original size of the input
    N = floor(size(f,1)./period).*period; % the number of samples that divide cleanly by period
    
    if ~isempty(input.timestamps) && size(input.timestamps,1)~=S(1)
        error('Timestamps data must be the same size as flux data!'); 
    end
    
    if ~isempty(input.variances) && ~isscalar(input.variances) && size(input.variances,1)~=S(1)
        error('Variance data must be the same size as flux data!'); 
    end
    
    f = f(:,:); % linearize all higher dimensions... 
    
    if util.text.cs(input.method, 'fast')
    
        if N~=size(f,1)
            f = vertcat(f,NaN(N+period - size(f,1), size(f,2))); % add NaNs to the end of the lightcurve to make it divisible by period
        end

        f = reshape(f, [period size(f,1)./period size(f,2)]); % now the 1st dimension has the length period, the 2nd has different sections at different periods and 3rd all the original other dims. 

        if input.median
            f_folded = nanmedian(f,2);
        else
            f_folded = nanmean(f,2);
        end
        
        num_bins = size(f,2); 

    %     f_folded = permute(f_folded, [1,3,2]); % get rid of dim 2 which we just integrated over
        f_folded = reshape(f_folded, [size(f_folded,1), S(2:end)]);
        
        
        if nargout>1% fold the variance as well

            if isempty(input.variances)
                var_folded = 1./num_bins; % the averaging reduces the variance
            elseif isscalar(input.variances)
                var_folded = input.variances./num_bins; % the averaging reduces the variance
            else

                v = input.variances; 

                if N~=size(f,1)
                    v = vertcat(vf,NaN(N+period - size(v,1), size(v,2))); % add NaNs to the end of the lightcurve to make it divisible by period
                end

                v = reshape(v, [period size(v,1)./period size(v,2)]); % now the 1st dimension has the length period, the 2nd has different sections at different periods and 3rd all the original other dims.

                var_folded = nanmean(v,2); 

                var_folded = reshape(var_folded, [size(var_folded,1), S(2:end)]);

            end

        end
        
    elseif util.text.cs(input.method, 'slow')
        
        if isempty(input.timestamps)
            t = 0:S(1)-1; 
        else
            t = input.timestamps;
        end
        
        dt = nanmedian(diff(t)); % the average time gap? this only works when most of the data is closely bunched with fairly constant period
        
        bins = 0:dt:period; 
        N_phase_bins = length(bins)-1; 
        
        t_folded = mod(t,period); 
        
        assignments = t_folded>=bins(1:end-1) & t_folded<bins(2:end); 
        
%         [mx, idx] = max(assignements,[],2); 
        
        f1 = permute(f, [1,3,2]); % move the different sources index to 3rd dimension
        Z = zeros(size(assignments), 'like', f1); 
        Z(~assignments) = NaN; 
        
        f2 = f1.*assignments + Z; % only have the flux value at the row/column that correspond to the correct period for that flux measurement
        
        if util.text.cs(input.method, 'median')
            f_folded = nanmedian(f2,1); 
        else
            f_folded = nanmean(f2,1); 
        end
        
        f_folded = permute(f_folded, [2,3,1]); 
        
        f_folded = reshape(f_folded, [N_phase_bins, S(2:end)]); % reshape back to original dimensions, except 1st, which is folded
        
    end
    
    if input.oversample
        f_folded = util.series.resample(f_folded, 1./input.oversample); 
    end
    
    if nargout>1
        
    end

