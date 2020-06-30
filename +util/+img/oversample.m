function I_new = oversample(I, factor)
% Usage: I_new = oversample(I, factor=2)
% Use fft-based sinc interpolation to oversample an image or data cube. 
% 
% Input: -an image or 3D or 4D data cube. 
%        -the factor you would like to use to increase the resolution. 
% NOTE: the factor*image size is rounded to nearest pixel. 
%

    if nargin==0, help('util.img.oversample'); return; end
    
    if nargin<2 || isempty(factor)
        factor = 2;
    end
    
    if ~isscalar(factor) || ~isnumeric(factor)
        error('The second input "factor" must be a numeric scalar!'); 
    end
    
    if factor<=1
        I_new = I;
        return;
    end
    
    N = isnan(I); % keep track of all the NaNs
    
    idx = find(N); 
    [x,y,f,s] = ind2sub(size(N), idx); 
    
    x = round((x-1).*factor) + 1;
    y = round((y-1).*factor) + 1;
    
    I(N) = 0;
    
    S = size(I); 
    if length(S)<3
        S = [S 1]; 
    end
    
    IF = fftshift(fftshift(fft2(I),1),2); 
    
    IF = util.img.pad2size(IF, [round(S(1).*factor), round(S(2).*factor), S(3:end)]); 

    I_new = ifft2(fftshift(fftshift(IF,1),2)); 
    
    idx = sub2ind(size(I_new), x, y, f, s); 

    mask = false(size(I_new)); 
    mask(idx) = 1;
    mask = imdilate(mask, ones(2*factor-1)); 

    I_new(mask) = NaN;

end