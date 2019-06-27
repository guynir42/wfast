function [values, radii] = radialStats(M, func, varargin)
% usage: [values, radii] = radialStats(M, func, varargin)
% run annulusPixels on all radii in M, do a function on the pixels or give
% them back as a cell array. 
% if func= "cell" will give the results. can choose any function without parameters... (std, var, mean, max, min...)

    import util.text.cs;
    import util.img.annulusPixels;
    
    if nargin==0, help('util.img.radialStats'); return; end
    
    if nargin<2 || isempty(func)
        func = 'cell';
    end
    
    r_delta = [];
    r_min = [];
    r_max = [];
    edges = 1;
    quadrants = [1 1 1 1];
    
    for ii = 1:2:length(varargin)
    
        if cs(varargin{ii}, {'r_delta', 'dr', 'r_step'})
            r_delta = varargin{ii+1};
        elseif cs(varargin{ii}, 'r_min')
            r_min = varargin{ii+1};
        elseif cs(varargin{ii}, 'r_max')
            r_max = varargin{ii+1};
        elseif cs(varargin{ii}, 'edges')
            edges = varargin{ii+1};
        elseif cs(varargin{ii}, 'quadrants')
            quadrants = varargin{ii+1};
        end
        
    end
    
    if isempty(r_min)
        r_min = 3;
    end
    
    if isempty(r_delta)
        r_delta = 2;
    end
    
    if isempty(r_max)
        r_max = floor(min(size(M,1))/2); 
    end
    
    radii = r_min:r_delta:r_max;
    
    N = length(radii)-1;
    
    if ischar(func) && cs(func, 'cell')
        values{N} = [];
    else
        values = zeros(1, N);
    end
    
    for ii = 1:N
        
        p = annulusPixels(M, radii(ii), radii(ii+1), quadrants, edges);
        
        if ischar(func) && cs(func, 'cell')
            values{ii} = p;
        else
            values(ii) = feval(func, p);
        end
        
    end
    
    radii = radii(1:end-1);
    