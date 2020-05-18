function datestrings = stamps2dates(varargin)
% Usage: datestrings = stamps2dates(varargin)
% Convert timestamp data from files into a cell array of date strings. 
%
% ARGUMENTS:
% -filenames: A string with a filename or cell array of filenames
% -datanames: What name is used for storing timestamps in the HDF5 file
%             (default is "timestamps"). 
% -timestamps: Skip the filenames and just give the timestamps. Must supply
%              the "end_string" and "end_stamp" as well. 
% -end_string: The date string of the file end time. 
% -end_stamp: The time stamp when "end_string" was set. 
% -output: Can specify which format to get the output, 
%          *string (default) formatted as FITS date strings.
%          *datetime objects (MATLAB class). 
%          *seconds, starting from 0 at the beginning of the list. 

% OUTPUT: 
% By default outputs a cell array of date strings. 
% Can also be datetime vector. 

    if isempty(varargin), help('img.stamps2dates'); return; end
    
    input = util.text.InputVars;
    input.input_var('filenames', [], 'names');
    input.input_var('datanames', 'timestamps');
    input.input_var('timestamps', []);
    input.input_var('end_string', []);
    input.input_var('end_stamp', []);
    input.input_var('output', 'string');
    input.scan_vars(varargin{:});
    
    if isempty(input.filenames) % make sure we got timestamps and end_string and end_stamp
        
        if isempty(input.timestamps)
            error('Cannot use "stamps2dates" without filenames or timestamps!');
        elseif isempty(input.end_string)
            error('Cannot use "stamps2dates" without filenames or end_string!');
        elseif isempty(input.end_stamp)
            error('Cannot use "stamps2dates" without filenames or end_stamp!');
        end        
        
    else % use the filenames given
        
        input.datanames = util.text.sa('/',input.datanames);
        
        if iscell(input.filenames)
            
            input.timestamps = [];
            input.end_stamp = [];
            input.end_string = [];
            
            for ii = 1:length(input.filenames)
                input.timestamps = horzcat(input.timestamps,  h5read(input.filenames{ii}, input.datanames)); % TODO: expand to cases where datanames is a cell
                input.end_stamp = horzcat(input.end_stamp, h5readatt(input.filenames{ii}, input.datanames, 't_end_stamp')); % TODO: expand to search for other names
                input.end_string = horzcat(input.end_string, h5readatt(input.filenames{ii}, input.datanames, 't_end')); % TODO: expand to search for other names
            end
            
        elseif ischar(input.filenames)
            input.timestamps = h5read(input.filenames, input.datanames); % TODO: expand to cases where datanames is a cell
            input.end_stamp = h5readatt(input.filenames, input.datanames, 't_end_stamp'); % TODO: expand to search for other names
            input.end_string = h5readatt(input.filenames, input.datanames, 't_end'); % TODO: expand to search for other names
        else
            error('Input "filenames" must be string or cell of strings! Instead got %s...', class(input.filenames));
        end
        
    end
    
    % now do the conversion. input.timestamps can be a column vector or
    % matrix with a column for each file. 
    
    if isscalar(input.end_stamp) 
        input.end_stamp = repmat(input.end_stamp, [1, size(input.timestamps,2)]); 
    end
    
    if ~iscell(input.end_string) 
        str = input.end_string;
        input.end_string = cell(size(input.end_stamp)); 
        input.end_string{:} = str;
    end
    
    datestrings = {};
    
    for ii = 1:size(input.timestamps,2)
        
        dates = util.text.str2time(input.end_string{ii});
        dates = repmat(dates, [size(input.timestamps,1),1]);
        dates.Second = dates.Second - input.end_stamp(ii) + input.timestamps(:,ii);
        
%         datestrings(:,ii) = util.text.time2str(dates);
        datestrings = vertcat(datestrings, dates);
        
    end

    if util.text.cs(input.output, 'string')
        datestrings =  util.text.time2str(datestrings);
    elseif util.text.cs(input.output, 'datetime')
       % pass
    elseif util.text.cs(input.output, 'seconds')
        N = datenum(datestrings);
        datestrings = (N-N(1))*24*3600;
    else
        error('Unknown output request "%s", use "string" or "datetime" instead...');
    end

end