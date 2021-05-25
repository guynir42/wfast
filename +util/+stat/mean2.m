function m = mean2(I)
% Usage: m = mean2(I)
% Calculate the mean of each 2D image in the data input I (ignoring NaNs!). 
% The output "m" will have the same dimensions as the input for all
% dimensions higher than 2, 1st and 2nd dimensions are 1x1.  

    if nargin==0, help('util.stat.mean2'); return; end
    
    S = size(I); 
    
    if ismatrix(I)
        m = nanmean(I(:)); % just find the global mean
    else
        I = reshape(I, [S(1).*S(2), 1, S(3:end)]); 
        m = nanmean(I,1); 
    end
    
end