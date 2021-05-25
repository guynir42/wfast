function [mx, idx] = max2(I)
% Usage: [mx, idx] = max2(I)
% Calculate the maximum of each 2D image in the data input I (ignoring NaNs!)
% The outputs mx and idx will have the same dimensions as the input for all
% dimensions higher than 2. The 1st and 2nd dimensions for the outputs will
% be 1x1 for mx and 1x2 for idx (the indices of the y and x position of the 
% maximum. 
% 

    if nargin==0, help('util.stat.max2'); return; end
    
    S = size(I); 
    
    if ismatrix(I)
        [mx, idx] = nanmax(I(:)); % just find the global maximum
    else
        I = reshape(I, [S(1).*S(2), 1, S(3:end)]); 
        [mx, idx] = nanmax(I,[],1); 
    end
    
    [row, col] = ind2sub(S(1:2), idx); 
    
    idx = [row, col]; 
    
end