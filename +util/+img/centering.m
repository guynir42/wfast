function I_out = centering(I, varargin)
% We need to finish this later... 

    if ndims(I)>4
        error('Cannot handle images with more than 4 dimensions... size(I)= %s.', util.text.print_vec(size(I), 'x')); 
    end

    [~, idx] = util.stat.max2(I); 
    
    cen = floor(size(I)/2) + 1; % image center
    cen = cen(1:2); % keep only 2 first dimensions
    
    I_out = zeros(size(I), 'like', I); 
    
    for ii = 1:size(I,3)
        
        for jj = 1:size(I,4)
    
            I_out(:,:,ii,jj) = circshift(I(:,:,ii,jj), cen-idx(:,:,ii,jj)); 
            
        end
        
    end
    
end