function v_row = topages(vector)
% if input is not a 3D vector (where all elements are in dim 3), 
% topages turns it to a 3D vector. if input vector is a matrix, topages 
% does nothing but displays a warning. 

    if nargin==0
        help('util.vec.topages');
        return;
    end
    
    if isempty(vector)
        v_row = vector;
        return;
    end
    
    vector = squeeze(vector);

    if ~isvector(vector)
        disp(['warning: input to "topages" is a matrix of size ' num2str(size(vector))]);
        v_row = vector;
        return;
    end
    
    if isrow(vector)
        vector = vector';
    end
    
    v_row = permute(vector, [3,2,1]);

end