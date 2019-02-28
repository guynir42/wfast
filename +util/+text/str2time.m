function time = str2time(str, timezone)
% parses a time string (FITS format compliant) into a datetime object
% usage: str2time(str, timezone='UTC')
% Optional argument defines the timezone (default: UTC); 

    if nargin==0
        help('util.text.str2time');
        return;
    end

    if isempty(str) || all(strcmp(str, '[]'))
        time = datetime.empty;
        return;
    end
    
    if nargin<2 || isempty(timezone)
        timezone = 'UTC';
    end
    
    if ischar(str)
        vec = sscanf(str, '%4d-%2d-%2dT%2d:%2d:%f')';
    elseif iscell(str)
        for ii = 1:length(str)
            vec(ii,:) = sscanf(str{ii}, '%4d-%2d-%2dT%2d:%2d:%f')';
        end
    else
        error('Must supply a string or cell of strings to "str2time"! Instead got %s.', class(str));
    end
    
    time = datetime(vec, 'TimeZone', timezone);
    
end