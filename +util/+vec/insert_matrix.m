function [M_out, idx_end] = insert_matrix(M_old, M_new, idx, fill_value, double_up)
% usage: [M_out, idx_end] = insert_matrix(M_old, M_new, idx, fill_value=0, double_up=0)
%
% Insert matrix M_new into matrix M_old at position idx (upper left corner). 
% If M_old is too small it will expand based on the following rules:
% 1) The empty elements will be filled by "fill_value" (default 0). 
% 2) For each dimension, if the new matrix is bigger than the old, the old 
%    just expands to contain it. If the new matrix is smaller then the old 
%    matrix can expand either to fit the new data or to double its size, if
%    doubling in size includes enough room for the new data. 
% 3) Turning off "double_up" (default 1) will cause the old matrix to expand
%    only enough to contain the new data. 
%
% The idea is to reduce the number of re-allocations (like in C++ vectors). 
% The user can use output "idx_end" to know what part of M_out contains data
% and what part is just preallocated space. 

    if nargin==0, help('util.vec.insert_matrix'); return; end

    if isempty(M_new)
        M_out = M_old;
        idx_end = idx;
        return;
    end
    
    if nargin<4 || isempty(fill_value)
        fill_value = 0;
    end
    
    if ~isscalar(fill_value)
        error('fill_value must be scalar! size(fill_value)= %s', util.text.print_vec(size(fill_value)));
    end
    
    if nargin<5 || isempty(double_up)
        double_up = 0;
    end
    
    if ~isscalar(double_up)
        error('double_up must be scalar! size(double_up)= %s', util.text.print_vec(size(double_up)));
    end
    
    if length(idx)<ndims(M_new)
        error('length(idx)= %d is smaller than ndims(M_new)= %d', length(idx), ndims(M_new));
    end
    
    if length(idx)<ndims(M_new)
        error('length(idx)= %d is smaller than ndims(M_new)= %d', length(idx), ndims(M_new));
    end
    
    dims = max([ndims(M_old), ndims(M_new), length(idx)]); % maximum dimension number
    
    S_old = size(M_old);
    if length(S_old)<dims
        S_old = [S_old ones(1, dims-length(S_old))]; % make sure the size vector is as long as the index
    end
    
    S_new = size(M_new);
    if length(S_new)<dims
        S_new = [S_new ones(1, dims-length(S_new))]; % make sure the size vector is as long as the index
    end
    
    S_expand = S_old;
    
    idx_end = idx + S_new - 1; % the index of the last positions of the new matrix inside the old. 
    
    % now idx and idx_end are the same length as S_old and S_new (all equal to dims)
    
    for ii = 1:dims
        
        if idx_end(ii)>S_old(ii) % need to expand on this dimension
            
            if double_up && S_new(ii)<S_old(ii) && S_old(ii).*2>idx_end(ii) % only double up if new data is smaller than old data (otherwise expand to contain just the new data)
                S_expand(ii) = S_old(ii).*2; % double up
            else
                S_expand(ii) = idx_end(ii); % just expand to fit the new data
            end
            
        end
        
    end
    
    % now make a new matrix if needed
    if any(S_expand>S_old)
        if fill_value==0
            M_out = zeros(S_expand, 'like', M_old);
        elseif isnan(fill_value)
            M_out = NaN(S_expand, 'like', M_old);
        elseif isinf(fill_value)
            M_out = Inf(S_expand, 'like', M_old);
        elseif isnumeric(fill_value)
            M_out = fill_value.*ones(S_expand, 'like', M_old);
        else
            error('Must input "fill_value" as a scalar numerical value (or NaN/Inf)... class(fill_value)= %s', class(fill_value));
        end
        
        % copy the old data into new matrix
        c = {};
        for ii = 1:dims
            c{ii} = 1:S_old(ii);
        end
        
        M_out(c{:}) = M_old;
        
    else
        M_out = M_old;
    end
    
    c = {};
    for ii = 1:dims
        c{ii} = idx(ii):idx_end(ii);
    end

    M_out(c{:}) = M_new;

end






