function V = var2(I)
% calculates the variance of images in the input matrix. 
% output can be 3D or higher. 

    import util.stat.mean2;
    
    if nargin==0
        help('util.stat.var2');
        return;
    end

    M = mean2(I); % could have 3rd or 4th non scalar dimensions... 

    V = mean2((I-M).^2);
   
end