function result = parse_value(value)
% Usage: result = parse_value(value)
% Take an input string or numeric or logical and figure out what it is by
% context: e.g., a string '2' will be turned into a numeric 2. 
% 
% String values that are parsed as numeric are NaN, Inf, true, false and []. 

    if nargin==0, help('util.text.parse_value'); return; end
    
    if isempty(value)
        result = [];
        return;
    elseif isnumeric(value) || islogical(value) || isa(value, 'datetime') % just pass numeric/logicals. Also datetime objects get through (presumably to be parsed by head.Epehemeris.parseTime)
        result = value;
        return;
    elseif ischar(value) % this is where the interesting part goes
       % pass! 
    elseif iscell(value)
        
        result = value; % match the length
        
        for ii = 1:length(value)
            result{ii} = util.text.parse_value(value{ii}); % recursively call this function on each word
        end
        
        if all(cellfun(@isnumeric, result))
            result = cell2mat(result); % if we got numeric values return a vector instead of cell
        end
        
        return; 
        
    else
        error('Cannot parse values of class "%s". Try using strings or numeric values instead...', class(value)); 
    end
    
    %%%%%%%%%%%%%% from here on we assume a string input value %%%%%%%%%%%%
    
    value = strtrim(value); % remove trailing whitespace
    words = strsplit(value, {',', ' '}); % split it into cells
    
    if length(words)>1
        
        result = words; % match the length
        
        for ii = 1:length(words)
            result{ii} = util.text.parse_value(words{ii}); % recursively call this function on each word
        end
        
        if all(cellfun(@isnumeric, result))
            result = cell2mat(result); % if we got numeric values return a vector instead of cell
        end
        
        return; % the parsing has been handled recursively, can skip the rest

    elseif length(words)==1
        value = words{1};
    else
        result = [];
        return;
    end
    
    %%%%%%%%%%%%%% from here on we assume a scalar input value %%%%%%%%%%%%
        
    if ~isnan(str2double(value))
        result = str2double(value); 
    elseif length(value)==2 && strcmp(value, '[]')
        result = [];
    elseif strcmpi(value, 'true')
        result = true;
    elseif strcmpi(value, 'false')
        result = false;
    elseif strcmpi(value, 'nan')
        result = NaN;
    else 
        result = value; 
    end
    
end