function [output, categories] = inverseSampling(dist, varargin)
% Usage: [output, categories] = inverseSampling(dist, varargin)
% Calculate a random number with the same distribution as the input values
% in "dist". It calculates the CDF of the inputs and then draws from a 
% uniform random distribution, getting a value for the corresponding value
% category (bin). 
% 
% OPTIONAL ARGUMENTS:
%   -size: the size of the output matrix, given as a vector. Default [1 1]. 
%          Note: if a scalar is given, it will produce an NxN matrix, so
%          make sure to always give at least two-elements. This is the same
%          odd behavior as in functions like zeros(). 
%   -width: size of bin width to be used when histogramming the "dist". 
%           If not given (or empty) will use histcounts() default width. 
%   -min: minimal values of "dist" to use when histogramming the "dist". 
%         Use this to remove low value outliers. Default is unlimited. 
%   -maxn: maximal values of "dist" to use when histogramming the "dist". 
%         Use this to remove high value outliers. Default is unlimited. 
%   -input: use this to specify is the input "dist" is already given as a 
%           distribution, PDF or CDF. Default is to run histcounts(dist). 
%           If giving a PDF or CDF, make sure to also give the "categories"
%           input. In these cases the categories must ordered (it only
%           makes sense to have a CDF if the PDF is for monotonically
%           growing x values). 
%   -categories: if the input is given as PDF or CDF, must give another 
%                vector of values (of equal length as "dist") to specify
%                what underlying value (category) each probability refers
%                to (the x axis to the PDF/CDF y axis). 
%
% EXAMPLE: inverseSampling(dist); // is the same as
%          [N,E]=histcounts(dist); 
%          inverseSampling(N, 'input', 'PSF', 'categories', E(1:end-1)); 
%
% refs: https://en.wikipedia.org/wiki/Inverse_transform_sampling
% refs: https://www.mathworks.com/matlabcentral/answers/152442-generating-random-numbers-from-histogram-data

    import util.text.cs;

    if nargin==0, help('util.stat.inverseSampling'); return; end
    
    if isempty(dist)
        error('Got empty distribution!'); 
    end
    
    input = util.text.InputVars; 
    input.input_var('size', [1,1]); % size of output 
    input.input_var('width', []); % bin width
    input.input_var('min'); % clip the distribution values below this
    input.input_var('max'); % clip the distribution values above this
    input.input_var('input', 'values'); % choose "values" or "PDF" or "CDF"
    input.input_var('categories', []); % choose the x values for PDF/CDF
    input.scan_vars(varargin{:}); 
    
    dist = dist(:); % linearize multiple dimensions
    
    if cs(input.input, 'values')
        
        if ~isempty(input.min)
            dist = dist(dist>=input.min);
        end
        
        if ~isempty(input.max)
            dist = dist(dist<=input.max);
        end
        
        dist = dist(~isnan(dist)); % remove NaNs
        dist = dist(~isinf(dist)); % remove infinities
        
        if isempty(input.width)
            [N,E] = histcounts(dist); 
        else
            [N,E] = histcounts(dist, 'BinWidth', input.width); 
        end
        
        dE = mean(diff(E)); % the step size
        
        x = E(2:end)-dE; % centers of bins
        
        pdf = N./sum(N); 
        
        cdf = cumsum(pdf); 
        
    elseif cs(input.input, 'PDF', 'CDF')
        
        x = input.categories;
        
        if isempty(x)
            error('Must supply categories (x axis) input when giving PDF/CDF data!'); 
        end
        
        if cs(input.input, 'CDF')
            cdf = dist;
        else
            pdf = dist./nansum(dist); 
            cdf = cumsum(pdf); 
        end
        
    else
        error('Unknown "%s" option to "input". Use "values" or "PSF" or "CDF". ', input.input); 
    end

    x = util.vec.tocolumn(x); 
    cdf = util.vec.tocolumn(cdf); 
    
    % now that we have a CDF it should be easy
    r = rand(input.size); 
    
    y = permute(cdf, [2:ndims(r)+1 1]); % change dimensions so the CDF has its vector on a dim bigger by one than ndims(r)
    idx = sum(r>y, ndims(r)+1)+1; % which bin of the CDF each value lives in. 
    output = x(idx); 
    
    categories = x; % make sure to also output the categories/bin centers
    
end










