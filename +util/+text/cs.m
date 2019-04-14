function val = cs(str, varargin)
% Compare strings. Usage: cs(str, varargin)
% If the string "str" is a match for any of the character vectors in
% varargin, then output is TRUE. 
% arguments can be given as a list of strings or as a cell array of strings. 
%
% Doesn't match the whole word, only the initial letters.
% if any argument is numeric, that is used as the minimal number of letters
% needed to compare the string. If none is given (or zero) then the number
% of letters in the given string is used. 
%
% Will remove spaces and underscores before comparing. Ignores case. 

    if nargin==0
        help('util.text.cs');
        return;
    end
    
    val = 0;

    if ~ischar(str)
        return;
    end

    cell_array = {};
    num_letters = [];

    for ii = 1:length(varargin) % make sure all values are in one, continuous cell

        if iscell(varargin{ii})
            cell_array = [cell_array, varargin{ii}{:}];
        elseif ischar(varargin{ii}) || isnumeric(varargin{ii})
            cell_array = [cell_array, varargin{ii}];        
        end

    end

    if isempty(str) && iscell(cell_array) && ~any(cellfun(@isempty, cell_array)) % if str is empty and any of the varargin inputs are not empty! 
        return;
    end

    idx_num = find(~cellfun(@ischar, cell_array));
    if ~isempty(idx_num)
        num_letters = cell_array{idx_num(end)};
        cell_array(idx_num) = [];
    end
    
    if isempty(str) &&  ~isempty(cell_array)
        return;
    end

    str = strrep(str, '_','');
    str = strrep(str, ' ','');

    if iscell(cell_array)
        for ii = 1:length(cell_array)
            cell_array{ii} = strrep(cell_array{ii}, '_','');
            cell_array{ii} = strrep(cell_array{ii}, ' ','');
        end
    else
        cell_array = strrep(cell_array,'_','');    
        cell_array = strrep(cell_array,' ','');
    end

    if isempty(num_letters) || num_letters==0
        num_letters = length(str);
    end

    val = any(strncmpi(str, cell_array, num_letters));

end