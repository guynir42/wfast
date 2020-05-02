function str = print_value(value, delimiter, logical_to_numeric)
% Usage: str = print_value(value, delimiter=' ', logical_to_numeric=1)
% Produce a formated string from a scalar or vector or cell array of values
% separated by a space (or other delimiter). 
% The value or values can be numeric, logical, or strings. 

    if nargin==0, help('util.text.print_value'); return; end
    
    if nargin<2 || isempty(delimiter)
        delimiter = ' ';
    end
    
    if nargin<3 || isempty(logical_to_numeric)
        logical_to_numeric = 1;
    end
    
    if isempty(value)
        str = '';
    elseif iscell(value)
        
        str = util.text.print_value(value{1}, delimiter);
        
        for ii = 2:length(value)
            str = [str delimiter util.text.print_value(value{ii}, delimiter)]; % recursively append these results
        end
        
    elseif isnumeric(value)
        str = util.text.print_vec(value, delimiter); 
    elseif islogical(value) && logical_to_numeric
        str = util.text.print_vec(value, delimiter); 
    elseif ischar(value)
        str = value;
    end
        
end