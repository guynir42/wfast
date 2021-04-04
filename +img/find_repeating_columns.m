function num_columns = find_repeating_columns(Im, varargin)
% Usage: num_columns = find_repeating_columns(Im, varargin)
% Scan a 3D data cube of images, or a 4D matrix of cutouts, and return the 
% number of columns that have at least 50% of the pixel values identical to
% the pixels 10 frames earlier (50% means top half or bottom half). 
% The output is for the number of columns that are copies of PREVIOUSLY 
% captured frames (lower index in the matrix), as it is understood that the
% first copy is the real data and the second copy is the corrupted data. 
%
% For example, if there are two columns in frame 22 that have either top or
% bottom half which is identical to the same pixels in frame 12, then the 
% output "num_columns" would have the value 2 in position 22. 
%
% If the input is 4D (cutouts), then each column in the output represents
% the result for a different cutout. I.e., the output dim1 is the same size
% as the input dim3 and the output dim2 is the same size as the input dim4. 
%
% OPTIONAL ARGUMENTS:
%   -gap: how many frames to look back to find identical columns. If gap=10
%         then frame 37 would be checked against frame 27. Default is 10. 
%   -margins: number of frames/images that are skipped from the beginning/end
%           of the dataset. You can enter a scalar or a 2-element vector, 
%           so that the elements specify the number of frames to skip from
%           the start and end, independently of each other. 
%           The default is zero, but whatever the value is, at least "gap" 
%           number of frames are skipped in the beginning of the frame, 
%           because those are frames that have no matching columns to 
%           compare against. E.g., if gap=10 and margins=5, the first 10 
%           frames and the final 5 frames would be skipped. 
%           Skipped frames just have zero output values, but the size of 
%           the output does not change. 
%   -fraction: the fraction of the frame size, either continuous from the 
%              top or from the bottom, that needs to be identical to count 
%              as a repeated column. Default is 0.5. 
%   -saturation: pixel values above this threshold are not counted as being
%                equal among rows. Default is 5e4. 

    if nargin==0, help('img.find_repeating_columns'); return; end
    
    input = util.text.InputVars; 
    input.input_var('gap', 10); 
    input.input_var('margins', 0); 
    input.input_var('fraction', 0.5); 
    input.input_var('saturation', 5e4);
    input.scan_vars(varargin{:}); 
    
    % first frame to scan (after gap or margins, whichever is bigger)
    start = input.margins(1); 
    
    if start<input.gap
        start = input.gap;
    end
    
    start = start + 1; 
    
    % last frame to scan
    if isscalar(input.margins)
        finish = size(Im,3)-input.margins(1);
    else
        finish = size(Im,3)-input.margins(2);
    end
    
    % pixel index for top and bottom part of the column 
    top = ceil(input.fraction.*size(Im,1));
    bottom = floor((1-input.fraction).*size(Im,1));
    
    % preallocate output matrix
    if islogical(Im)
        num_columns = zeros(size(Im,3),size(Im,4), 'uint16'); 
    else 
        num_columns = zeros(size(Im,3),size(Im,4), 'like', Im); 
    end        
    
    % loop over frames and cutouts
    for ii = start:finish
        
        cols_top = all(squeeze(Im(1:top,:,ii,:))==squeeze(Im(1:top,:,ii-input.gap,:)));
        cols_bottom = all(squeeze(Im(bottom:end,:,ii,:))==squeeze(Im(bottom:end,:,ii-input.gap,:)));

        for jj = 1:size(Im,4) % verify that those columns are not all saturated
            
            cols_top_index = find(cols_top(1,:,jj)); 
            for kk = cols_top_index
                if all(Im(1:top,kk,ii,jj)>=input.saturation)
                    cols_top(1,kk,jj) = 0; 
                end
            end
            
            cols_bottom_index = find(cols_bottom(1,:,jj)); 
            for kk = cols_bottom_index
                if all(Im(bottom:end,kk,ii,jj)>=input.saturation)
                    cols_bottom(1,kk,jj) = 0; 
                end
            end
            
        end % for jj (cutout number)
        
        num_columns(ii,:) = permute(sum(cols_top | cols_bottom, 2), [1,3,2]); 
        
    end % for ii (frame number)
    
end









