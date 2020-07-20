function name = legalize(name, preamble)
% Usage: str_out = legalize(name, preamble='obj')
% Make an object name legal to use in folder names and struct fields. 
% Will replace spaces and illegal characters with '_', and if the name 
% starts with a non-letter character, it will be pre-appended with the 
% string in "preamble" (the default is 'obj'). 

    if nargin==0, help('util.text.legalize'); return; end
    
    if nargin<2 || isempty(preamble)
        preamble = 'obj'; 
    end
    
    name = strip(name); 
    name = replace(name, {' ', '!', '@', '#', '$', '%', '^', '&', '*', '(', ')', '`', '~', '/', '\', ';', ':', '-', '+', '|'}, '_'); 
    
    if isempty(name) || ~isletter(name(1))
        name = [preamble name]; 
    end
    
end