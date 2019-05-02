classdef Logger < handle
% Keeps text file with logging commands and errors from hardware.
%
% Constructor can be called with single text argument to set <name>. 
%
% Use obj.input(text, errflag=0) to input normal operation logs, that will
% be saved to file, with timestamp automatically appended. 
% 
% Use obj.error(text) to produce an error log. This error is saved with a 
% timestamp and also a header showing the error. It will optionally also be
% saved to a separate error file. 
%
% Will generate a new file every day (day starts at 12:00 UTC). 
%
% Filename will be <base_dir>/YYYY-MM-DD_<name> where properties are:
%   *name: the name of hardware/component/sensor (e.g., "mount")
%   *base_dir: the folder where all logfiles are kept. 
%    If left empty, will write to fullfile(getenv('DATA'),'WFAST/logfiles')
%   *use_error_file: if non-zero, will write errors to main file and to 
%    additional error file with same name, appended with _ERROR. 
%
% EXAMPLE: 
% L = util.sys.Logger('mount-ASA'); 
% L.input('Connecting to mount');  % create logfile at C:\Dropbox\data\WFAST\logfiles 
%                                  % with name 2019-04-29_mount-ASA.txt
%                                  % with text: 15:27:34.123 Connecting to mount
% L.error('Mount not responding'); % write this error message (with timestamp
%                                  % and error banner to same file, and also
%                                  % (optionally) to logfile 2019-04-29_mount-ASA_ERROR.txt
% L.base_dir = 'D:\logifles';      % will overwrite the environmental variable 'DATA' 
%
%
    
    properties(Transient=true)
        
    end
    
    properties % objects
        
        hndl; % file handle (using fopen)
        hndl_err; % secondary handle to additional log file (errors only)
        time; % datetime object
        
    end
    
    properties % inputs/outputs
        
        report; % combination of datestr+message, that is written to file. 
        
    end
    
    properties % switches/controls
        
        dev_name = 'log'; % name of device/sensor generating the log (e.g., "mount")
        base_dir = ''; % default is fullfile(getenv('DATA'), 'WFAST/logfiles')
        filename = ''; % file name of log file (e.g., 2019-04-29_mount_ASA.txt)
        use_error_file = 0; % to keep an extra file for errors only
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        fullname; % full path to text file (<basedir>/<filename>)
        
    end
    
    properties(Hidden=true)
       
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = Logger(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'util.sys.Logger')
                if obj.debug_bit
                    fprintf('Logger copy-constructor v%4.2f\n', obj.version);
                end
                obj = util.oop.full_copy(varargin{1});
                
            elseif ~isempty(varargin) && ischar(varargin{1})
                
                obj.dev_name = varargin{1};
                if obj.debug_bit
                    fprintf('Logger constructor v%4.2f | name: %s\n', obj.version, obj.dev_name);
                end
                
            else
                if obj.debug_bit
                    fprintf('Logger constructor v%4.2f\n', obj.version);
                end
            end
            
        end
        
        function delete(obj)
            
            try 
            
                if ~isempty(obj.hndl)
                    fclose(obj.hndl);
                end
            
            catch ME
                warning(ME.getReport);
            end
            
        end
        
    end
    
    methods % reset /clear
        
        function reset(obj)
            
            obj.time = datetime.empty;
            
            try 
                
                if ~isempty(obj.hndl)
                    fclose(obj.hndl);
                end
                
                if ~isempty(obj.hndl_err)
                    fclose(obj.hndl_err);
                end
                
            catch ME
                warning(ME.getReport);
            end
            
            obj.hndl = [];
            obj.hndl_err = [];
            
        end
        
    end
    
    methods % getters
        
        function val = get.base_dir(obj) % allow user to override the <base_dir> property
        
            if isempty(obj.base_dir)
                val = fullfile(util.def.data_folder, 'WFAST/logfiles');
            else
                val = obj.base_dir;
            end
            
        end
        
        function val = get.filename(obj) % allow user to override the <filename> property
            
            if isempty(obj.filename)
                val = [util.sys.date_dir(obj.time) '_' obj.dev_name '.txt'];
            else
                val = obj.filename;
            end
            
        end
        
        function val = get.fullname(obj) % combine <basedir>/<filename>
            
            val = fullfile(obj.base_dir, obj.filename);
            
        end
        
        function val = err_filename(obj) % error filename is the same as <filename> with _ERROR appended
            
            [~, name, ext] = fileparts(obj.filename);
            
            val = [name '_ERROR' ext];
            
        end
        
        function val = err_fullname(obj) % combine <basedir>/<filename>
            
            val = fullfile(obj.base_dir, obj.err_filename);
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % commands
        
        function input(obj, text, errflag) % give bit of text to write to log files. Use errflag=1 to denote errors (default=0) 
            
            if nargin<3 || isempty(errflag)
                errflag = 0;
            end
            
            text = util.text.eraseTags(text);
            
            new_time = datetime('now', 'timezone', 'UTC'); % update to current time
            
            % if we passed noon UTC we should generate new log files
            if ~isempty(obj.time) && new_time.Day>=obj.time.Day && new_time.Hour>=12 && obj.time.Hour<12
                obj.reset;
            end
            
            obj.time = new_time;
            
            timestamp = datestr(obj.time, 'hh:MM:ss.FFF');
            
            if errflag==0
                obj.report = sprintf('%s: %s', timestamp, text); % timestamp and input text 
            else
                obj.report = sprintf('************ EXCEPTION *****************\n%s: %s', timestamp, text); % add banner
            end
            
            if isempty(obj.hndl) || obj.hndl<0
                obj.makeFile; % make sure the file exists or create a new one if needed
            end
            
            if obj.use_error_file && (isempty(obj.hndl_err) || obj.hndl_err<0)
                obj.makeErrFile;
            end
            
            
            try 
                
                fprintf(obj.hndl, '%s\n', obj.report); % print line to normal file
                
                if obj.use_error_file && errflag
                    fprintf(obj.hndl_err, '%s\n', obj.report); % print error line to additional error file
                end
                
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function error(obj, text)
            
            obj.input(text, 1); % use the errflag! 
            
        end
        
        function makeFile(obj)
            
            dir = obj.base_dir;
            
            if ~exist(dir, 'dir')
                mkdir(dir)
            end
            
            try 
                
                obj.hndl = fopen(obj.fullname, 'at');
                
                if obj.hndl<0
                    obj.hndl = [];
                    error('Cannot open the file %s', obj.fullname);
                end
                
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function makeErrFile(obj)
            
            dir = obj.base_dir;
            
            if ~exist(dir, 'dir')
                mkdir(dir)
            end
            
            try 

                obj.hndl_err = fopen(obj.err_fullname, 'at');

                if obj.hndl_err<0
                    obj.hndl_err = [];
                    error('Cannot open the file %s', obj.err_fullname);
                end

            catch ME
                warning(ME.getReport);
            end 
            
        end
        
    end
    
end

