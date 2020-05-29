function [cutouts, positions, dims] = jigsaw(Im, varargin)
% Usage: [cutouts, positions, dims] = jigsaw(Im, varargin)
% Cuts an image up into squares (blocks or tiles). 
% Can be useful for making statistical tests on small regions of the image.
% OPTIONAL ARGUMENTS
%   *tile: height and width of tiles (Default 128).
%   *pad_value: how to pad the last tiles in a row/column (Default 0)
%   *overlap: allow N pixels of overlap between tiles (Default no overlap)
%   *squeeze: The cutouts convention is that 3rd dimension is for frames in
%             the batch, e.g., a 3D input matrix of 2048x2048x100 will give
%             an output cutouts of 128x128x100x256 where 4th dimension is
%             the index of the cutout. For a single image this will give a
%             singlton 3rd dimension. So using "squeeze" this dimension is
%             removed. Default is true for convenience. If you expect to
%             get 2D and 3D inputs, use squeeze=0 to get consistent output.
%   *partial: keep the tiles on the edges of the image that don't fully fit
%             into the size of the image (and get filled). Default true. 
%             
% OUTPUTS: the cutouts, the positions matrix (2 columns, one for x and the 
%          other for y, and the dimensions of the x and y of the positions. 

    if nargin==0, help('util.img.jigsaw'); return; end

    input = util.text.InputVars;
    input.use_ordered_numeric = 1;
    input.input_var('tile', 128, 'tile_size', 'size');
    input.input_var('pad_value', 0);
    input.input_var('overlap', 0);
    input.input_var('squeeze', 1);
    input.input_var('partial', true, 'keep_partial'); 
    input.scan_vars(varargin{:});
    
    if isempty(Im)
        cutouts = [];
        positions = []; 
        return;
    end
    
    t = input.tile;
    if isscalar(t)
        t = [t t];
    end
    
    t = ceil(t); % make sure tile sizes are integer! 
    
    S = util.vec.imsize(Im);
    
    if all(t>=S) % if the tile is larger than image return single cutout/position
        cutouts = Im;
        positions = flip(floor(S/2)+1);
        return;
    end
    
    Nx = ceil(S(2)./(t(2)-input.overlap));
    Ny = ceil(S(1)./(t(1)-input.overlap));
    
    positions = zeros(Nx.*Ny,2);
    
    if input.pad_value==0
        cutouts = zeros(t(1),t(2),size(Im,3),Nx.*Ny);
    elseif isnan(input.pad_value) || (ischar(input.pad_value) && isnan(input.pad_value))
        cutouts = nan(t(1),t(2),size(Im,3),Nx.*Ny);
    elseif isnumeric(input.pad_value)
        cutouts = ones(t(1),t(2),size(Im,3),Nx.*Ny).*input.pad_value;
    end
        
    need_clip_x = 0;
    need_clip_y = 0;
    break2 = 0; 
    counter = 1;
    
    for ii = 1:Ny
        
        if break2, break; end
        
        for jj = 1:Nx
            
            start_y = (ii-1)*(t(1)-input.overlap)+1;
            end_y = start_y + t(1) - 1; 
            
            start_x = (jj-1)*(t(2)-input.overlap)+1;
            end_x = start_x + t(2) - 1;
            
            if input.partial
                if end_y>S(1), end_y = S(1); end
                if end_x>S(2), end_x = S(2); end
            else
                if end_x>S(2)
%                     need_clip_x = 1;
                    Nx = jj-1;
                    break; % in this case just skip the partial tiles
                elseif end_y>S(1)
%                     need_clip_y = 1;
                    Ny = ii-1;
                    break2 = 1;
                    break;
                end
            end
            
            C = Im(start_y:end_y,start_x:end_x,:);
            cutouts(1:size(C,1),1:size(C,2),:,counter) = C;
            positions(counter,1) = (end_x+start_x)/2;
            positions(counter,2) = (end_y+start_y)/2;
                        
            counter = counter + 1;
            
        end
    end

    cutouts = cutouts(:,:,:,1:counter-1); 
    
    positions = positions(1:counter-1,:); 
    
    if input.squeeze
        cutouts = squeeze(cutouts);
    end
    
%     if need_clip_x
%         Nx = Nx - 1;
%     end
%     
%     if need_clip_y
%         Ny = Ny - 1;
%     end
    
    dims = [Nx, Ny]; 
    
end