function M_out = pad2size(M_in, size_needed, pad_value)
% usage: pad2size(M_in, size_needed)
% pads the given array M_in to size size_needed, putting it in the middle. 
% Won't shrink the array... can handle 3D or 4D matrices (expands only first 2 dims)
    
    if nargin==0, help('util.img.pad2size'); return; end

    if nargin<3 || isempty(pad_value)
        pad_value = 0;
    end
    
    size_needed = util.vec.imsize(size_needed);
    
    if isempty(M_in)
%         error('Cannot pad an empty matrix...');        
        M_out = zeros(size_needed);
        return;
    end
    
    S_in = size(M_in);
    S_in = S_in(1:2);
    
    if size_needed(1)>S_in(1) || size_needed(2)>S_in(2)
        
        if pad_value==0
            M_out = zeros(max(size_needed(1), S_in(1)), max(size_needed(2), S_in(2)), size(M_in,3), size(M_in,4));
        elseif isnan(pad_value)
            M_out = nan(max(size_needed(1), S_in(1)), max(size_needed(2), S_in(2)), size(M_in,3), size(M_in,4));
        elseif isnumeric(pad_value)
            M_out = pad_value.*ones(max(size_needed(1), S_in(1)), max(size_needed(2), S_in(2)), size(M_in,3), size(M_in,4));
        else
            error('"pad_value" must be numeric or NaN!');
        end
        
        gap = (size_needed-S_in)/2;
        
        y1 = 1+ceil(gap(1));
        y2 = size_needed(1)-floor(gap(1));
        x1 = 1+ceil(gap(2));
        x2 = size_needed(2)-floor(gap(2));
        
        M_out(y1:y2, x1:x2,:,:) = M_in;
        
    else
        M_out = M_in;
    end
    
end