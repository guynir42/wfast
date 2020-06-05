function M = weighted_average(values, weights, dim)
    
    if nargin<3 || isempty(dim)
        dim = 1;
    end
    
    M = nansum(values.*weights, dim)./nansum(weights.*~isnan(values), dim);  
    
end