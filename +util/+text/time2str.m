function str = time2str(time)
% turns a datetime object into a string (FITS format compliant)
    
    if nargin==0
        help('util.text.time2str');
        return;
    end

    if builtin('isempty', time) || isempty(time)
        str = '';
        return;
    end

    vec = datevec(time);

    if isscalar(time)
        str = sprintf('%4d-%02d-%02dT%02d:%02d:%06.3f', vec(1), vec(2), vec(3), vec(4), vec(5), vec(6)); 
    else
        
        for ii = 1:size(vec,1)
            str{ii} = sprintf('%4d-%02d-%02dT%02d:%02d:%06.3f', vec(ii,1), vec(ii,2), vec(ii,3), vec(ii,4), vec(ii,5), vec(ii,6)); 
        end
        
        str = util.vec.tocolumn(str);
        
    end
    
end