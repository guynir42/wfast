function val = parse_bool(input)
% converts boolean, numbers or text to numeric zero or non-zero.
% If given a string, the values 'yes', 'on', 'true' will output 1.
% All other strings output 0.
% Numerical inputs are returned as is. 
% strings containing numerals are turned to numbers and returned 
% e.g. the input '1' is true. 

    if nargin==0
        help('util.text.parse_bool');
        return;
    end

    if isnumeric(input) || islogical(input)
        val = input;
        return;
    end
    
    val = str2num(input);
    if ~isempty(val)
        return;
    end
    
    if util.text.cs(input, {'yes','on','true'})
        val = 1;
    else
        val = 0;
    end
    

end