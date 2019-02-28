function str_out = strPlusPlus(str_in, increase)
% usage: str_out = strPlusPlus(str_in, increase=1)
% parses the number at the end of the string "str_in" and increases it by 
% "increase" (default 1). 
% If there is no number, just add a number (the value of "increase"). 

    if nargin<2 || isempty(increase)
        increase = 1;
    end
    
    if ~isnumeric(increase)
        error(['Can only accept numeric values of "increase". Got ' num2str(increase)]);
    end
    
    [a,b] = regexp(str_in, '\d+$');
    
    if ~isempty(a)
        numeral = str_in(a:b);
        num = str2double(numeral); 
        num = num + increase;
        str_out = [str_in(1:a-1) num2str(num)];
    else
        str_out = [str_in num2str(increase)];
    end
    
end