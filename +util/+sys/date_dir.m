function val = date_dir(date)
    
    if nargin<1 || isempty(date)
        date = datetime('now', 'TimeZone', 'UTC');
    elseif ~isa(date, 'datetime')
        % ... convert to datetime object! 
    end

    % make sure each night is in the same folder, even after midnight.
    if date.Hour<12 % this is 14:00 local time, pretty safe to say it is after midnight of the previous night but observations still haven't started on the new night
        date.Day = date.Day - 1;
    end

    val = datestr(date, 'yyyy-mm-dd');
    
end