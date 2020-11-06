function corr = correlation(v1, v2, timescale, varargin)
% Usage: corr = correlation(v1, v2, timescale=10, varargin)
% Correlate the two vectors v1 and v2, to see if any part of the series has 
% the same kinds of variations. 
% 
% Calculates v1.*v2 (after removing their mean) and divides by the sum of 
% squares of each series, over the given timescale. 
% The third input is to decide how many frames should be summed when 
% calculating the correlation (in numerator and denominator). 
% The default is 10 frames.  
%
% If v1 or v2 are not single column vectors, they are expanded to matching 
% dimensions just like you get by multiplying the two. 
%
% Output: the correlation coeff summed over a moving window along the time-
%         series, giving the estimate of correlation in each point. 

    if nargin==0, help('util.series.correlation'); return; end
    
    if nargin<3 || isempty(timescale)
        timescale = 10;
    end
    
    input = util.text.InputVars;
    input.use_ordered_numeric = 1; 
    input.input_var('indices', []); 
    input.input_var('movmean', true); % I can probably get rid of the old method with convolutions!
    input.scan_vars(varargin{:}); 
    
    v1 = v1-nanmean(v1);
    v2 = v2-nanmean(v2); 
    
    if input.movmean
        
        numer = movmean(v1.*v2, timescale, 'omitnan');
        sum1 = movmean(v1.^2, timescale, 'omitnan');
        sum2 = movmean(v2.^2, timescale, 'omitnan');
        
    else % old method with convolutions... 
        
        kernel = ones(timescale,1);

        if isa(v1, 'single') && isa(v2, 'single')
            kernel = single(kernel); 
        end

        numer = util.vec.convolution(kernel, v1.*v2); 
        sum1 = util.vec.convolution(kernel, v1.^2); 
        sum2 = util.vec.convolution(kernel, v2.^2);
        
    end
    
    bad_idx = sum1<0 | sum2<0;
   
    corr = numer./sqrt(sum1.*sum2); 
    corr(bad_idx) = 0; 
    
    if ~isempty(input.indices)
        s = substruct('()', repmat({':'}, [1 ndims(corr)]));
        s.subs{1} = input.indices; 
        corr = subsref(corr, s); 
    end
    
end






