classdef Logger < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        hndl; % file handle
        time@datetime;
        
    end
    
    properties % inputs/outputs
        
        filename; % full path to text file
        
        report; % combination of datestr+message, that is written to file. 
        
    end
    
    properties % switches/controls
        
        name = 'log';
        
        debug_bit = 1;
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Logger(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'util.sys.Logger')
                if obj.debug_bit, fprintf('Logger copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
                return;
            elseif ~isempty(varargin) && ischar(varargin{1})
                
                obj.name = varargin{1};
                if obj.debug_bit, fprintf('Logger constructor v%4.2f | name: %s\n', obj.version, obj.name); end
                
            else
                if obj.debug_bit, fprintf('Logger constructor v%4.2f\n', obj.version); end
                
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
            
            obj.filename = '';
            
            obj.time = datetime.empty;
            
            try 
                
                if ~isempty(obj.hndl)
                    fclose(obj.hndl);
                end
                
            catch ME
                warning(ME.getReport);
            end
            
            obj.hndl = [];
            
        end
        
    end
    
    methods % calculations
        
        function input(obj, text)
            
            obj.time = datetime('now', 'timezone', 'UTC');
            timestamp = datestr(obj.time, 'hh:MM:ss.FFF');
            obj.report = [timestamp ': ' text];
            
            if isempty(obj.hndl) || obj.hndl<0
                obj.makeFile;
            end
            
            try 
                fprintf(obj.hndl, '%s\n', obj.report);
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function makeFile(obj)
            
            dir = fullfile(getenv('DATA'), 'WFAST/logfiles');
            
            if ~exist(dir, 'dir')
                makedir(dir);
            end
            
            if isempty(obj.filename)
                obj.filename = [obj.name '_' util.sys.date_dir(obj.time) '.txt'];
            end
            
            try 
                
                fullname = fullfile(dir, obj.filename);
                
                obj.hndl = fopen(fullname, 'at');
                
                if obj.hndl<0
                    obj.hndl = [];
                    error('Cannot open the file %s', fullname);
                end
                
            catch ME
                warning(ME.getReport);
            end
            
        end
        
    end
    
end

