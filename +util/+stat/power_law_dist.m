function values = power_law_dist(index, varargin)
% Usage: values = power_law_dist(index, varargin)
% Get random numbers from a decreasing power law distribution. 
% The index input will be treated as negative either way it is given. 
% 
% OPTIONAL ARGUMENTS:
%   -size: specify output size as a 2-element vector at least. This is the 
%          same strange behavior as in the zeros() or rand() functions. 
%   -min/max: the ends of the distribution to draw from. The "min" should 
%             be chosen as a reasonable value, while "max" can be left as
%             Inf which is fine for descending distributions, as long as 
%             you don't mind occasionally getting very high values. 
%   
    
    if nargin==0, help('util.stat.power_law_dist'); return; end
    
    input = util.text.InputVars;
    input.input_var('size', [1 1]); 
    input.input_var('min', 1);
    input.input_var('max', Inf); 
    input.scan_vars(varargin{:}); 
    
    if index==0
        error('Cannot produce a power law with index=0...'); 
    end
    
    index = abs(index); 
    
    r = rand(input.size); % a random point 0<r<1, we can use inverse sampling agains the CDF
    
    
    % the PDF is proportional to x^-q
    % the CDF is proportional to x^(1-q)
    % the normalization of the CDF is N. 
    % Then we integrate the PDF and solve for x. 
    if index==1
        
        N = log(input.max)-log(input.min);
        
        values = input.min.*exp(r.*N); 
        
    else
        
        N = (input.min.^(1-index)-input.max.^(1-index));
        
        values = input.min + (input.min.^(1-index) - r.*N).^(1./(1-index)); 
        
    end
    
    
end