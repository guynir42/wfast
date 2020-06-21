function values = parse_inputs(str, delims, remove)
% Usage: values = parse_inputs(str, delims={',','|','='}, remove={''', '"'})
% Turn a string of keyword=values into a cell array that can be parsed as 
% varargin to a function. 
% Use "delims" (default {',','|','='}) to separate keys and values. 
% Use "remove" (default {'''', '"'}) to get rid of unwanted characters like
% quote marks that the user may have accidentally inserted for strings. 
% Use remove={} to not remove any characters...

    if nargin==0, help('util.text.parse_inputs'); return; end
    
    if nargin<2 || isempty(delims)
        delims = {',','|','='};
    end
    
    if nargin<3 || (~iscell(remove) && isempty(remove))
        remove = {'''','"'};
    end
    
    values = strsplit(str, delims); 
    
    for ii=1:length(values)
        values{ii} = strtrim(values{ii});
    end
    
    values = values(~cellfun('isempty',values));

    for ii=1:length(values)
        
        for jj = 1:length(remove)
            values{ii} = strrep(values{ii}, remove{jj}, '');
        end
        
        values{ii} = strtrim(values{ii});
        
        if mod(ii,2)==0 % only for values (not for keywords
            values{ii} = util.text.parse_value(values{ii}); 
        end
        
    end
    
    
    
end