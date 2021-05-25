function m = mad2(I, varargin)
% Usage: m = mad2(I, varargin)
% Finds the Median Absolute Deviation of an image (MAD). 
% More robust than variance/standard deviation. 
%
% Input: A matrix of arbitrary dimensions. Calculates MAD for each 2D image. 
% Output: The median (or mean) absolute deviation (3rd and higher dims are 
% the same size as input).
% OPTIONAL ARGUMENTS:
%   -type: choose 'median' (default) or 'mean'
%   -multiply or sigma: use this flag to multiply the return values to
%    match the normal distribution sigma parameter (default is "off").
%

    import util.text.cs;
    import util.text.parse_bool;

    if nargin==0, help('util.stat.mad2'); return; end

    type = 1; % this is for 'median', 0 is for 'mean'
    sigma = 0;
    
    for ii = 1:2:length(varargin)
        key = varargin{ii};
        val = varargin{ii+1};
        
        if cs(key, 'type')
            if cs(val, 'median')
                type = 1;
            elseif cs(val, 'mean')
                type = 0;
            else
                error(['Unknown type option "' val '". Use median or mean...']);
            end
        elseif cs(key, 'sigma', 'multiply')
            sigma = parse_bool(val);
        end            
        
    end
    
    S = size(I); 
    
    if ismatrix(I)
        m = mad(I(:), type); 
    else
        I = reshape(I,[S(1).*S(2), 1, S(3:end)]); 
        m = mad(I, type,1); 
    end
    
    if sigma % match to standard deviation (for normal dist. only!)
        if type
            m = m*1.4826;
        else
            m = m*1.253;
        end
    end
    
end