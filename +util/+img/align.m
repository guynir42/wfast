function I_new = align(I, x, y, varargin)
% Usage: I_new = align(I, x, y, varargin)
% Shift each image/cutout in a 2D, 3D or 4D matrix, based on the x,y, values
% given (in units of pixels). The x,y are the shift needed (minus of the 
% image centeroids, if you want to have all images centered). 

    if nargin==0, help('util.img.align'); return; end

    x(isnan(x)) = 0;
    y(isnan(y)) = 0; 
    
    I_new = I; 
    
    for ii = 1:size(I,4)
                
        for jj = 1:size(I,3)

            if nnz(isnan(I(:,:,jj,ii)) | I(:,:,jj,ii)==0)<0.2*numel(I(:,:,jj,ii)) % less than 20% of the pixels are NaN or zero
                try
                    I2 = I(:,:,jj,ii); 
                    I2(isnan(I2)) = 0; 
                    I2 = real(util.img.FourierShift2D(I2,[x(jj,ii),y(jj,ii)])); 
                    I_new(:,:,jj,ii) = I2; 
                catch ME
                    if ~strcmp(ME.identifier, 'MATLAB:regionfill:expectedNonNaN') % ignore these errors silently 
                        rethrow(ME);
                    end
                end
            end

        end

    end
    
end