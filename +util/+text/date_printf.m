function str = date_printf(str, varargin)
% Usage: str = date_printf(str, varargin)
% Print the string "str", as though given to fprintf(), but also add a date
% and time at the beginning of the string. 
% Example: 
% >> date_printf('value= %d', 10) 
%      2021-04-10 20:30:32.567: value= 10
% The date and time are from the current time, in UTC. 
% The function gives the string without the date, to allow the same string
% to be sent into other functions. It always prints to terminal. 
% 

    if nargin==0, help('util.text.disp_printf'); return; end
    
    str = sprintf(str, varargin{:}); 
    
    now = datetime('now', 'TimeZone', 'UTC', 'Format', 'uuuu-MM-dd HH:mm:ss.sss'); 
    
    disp([char(now) ': ' str]); 
    
    if nargout==0
        clear str;
    end
    
end