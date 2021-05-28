function [mn, idx] = min2(I)
% Usage: [mn, idx] = max2(I)
% Calculate the minimum of each 2D image in the data input I (ignoring NaNs!)
% The outputs mn and idx will have the same dimensions as the input for all
% dimensions higher than 2. The 1st and 2nd dimensions for the outputs will
% be 1x1 for mn and 1x2 for idx (the indices of the y and x position of the 
% maximum. 
% 

    if nargin==0, help('util.stat.min2'); return; end
    
    [mn, idx] = util.stat.max2(-I); 
    
    mn = - mn; 