function str = print_vec(vector, delimiter, pre_str, post_str)
% prints a row vector into a string, like num2str but with single spaces.
% usage: print_vec(vector, delimiter= ' ', pre_str='', post_str='')

    if nargin==0, help('util.text.print_vec'); return; end

    if nargin<2 || isempty(delimiter)
        delimiter = ' ';
    end

    if nargin<3 || isempty(pre_str)
        pre_str = '';
    end
    
    if nargin<4 || isempty(post_str)
        post_str = '';
    end
    
    if ~isempty(vector)
        str = regexprep(num2str(util.vec.torow(vector)),'\s+',delimiter);
    else
        str = '';
    end
    
    str = [pre_str str post_str];
    
end