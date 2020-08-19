function corr = correlation(v1, v2, timescale, varargin)
% Usage: corr = correlation(v1, v2, timescale=[], varargin)
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
    input.scan_vars(varargin{:}); 
    
    v1 = v1-nanmean(v1);
    v2 = v2-nanmean(v2); 
    
    kernel = ones(timescale,1);
    
    if isa(v1, 'single') && isa(v2, 'single')
        kernel = single(kernel); 
    end
    
    % these are replaced with util.vec.convolution because filter2 can't handle 3D inputs properly... 
%     numer = filter2(kernel,  v1.*v2); 
%     sum1 = filter2(kernel, v1.^2); 
%     sum2 = filter2(kernel, v2.^2); 
    
    numer = util.vec.convolution(kernel, v1.*v2); 
    sum1 = util.vec.convolution(kernel, v1.^2); 
    sum2 = util.vec.convolution(kernel, v2.^2);
    
    bad_idx = sum1<0 | sum2<0;
   
    corr = numer./sqrt(sum1.*sum2); 
    corr(bad_idx) = 0; 
    
    if ~isempty(input.indices)
        s = substruct('()', repmat({':'}, [1 ndims(corr)]));
        s.subs{1} = input.indices; 
        corr = subsref(corr, s); 
    end
    
end






