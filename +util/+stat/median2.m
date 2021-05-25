function m = median2(I)
% Usage: m = median2(I)
% Calculate the median of each 2D image in the data input I (ignoring NaNs!). 
% The output "m" will have the same dimensions as the input for all
% dimensions higher than 2, 1st and 2nd dimensions are 1x1.  

    if nargin==0, help('util.stat.median2'); return; end
    
    S = size(I); 
    
    if ismatrix(I)
        m = nanmedian(I(:)); % just find the global median
    else
        I = reshape(I, [S(1).*S(2), 1, S(3:end)]); 
        m = nanmedian(I,1); 
    end
    
end