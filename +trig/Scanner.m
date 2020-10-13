classdef Scanner < handle
% This class overlooks analysis campaigns on multiple runs. 
% There are a few different functions that can be run on different objects:
% (a) Continuous analysis program: define a folder and possible start/end
%     times and let the object push all runs into Analysis on parallel
%     workers until all the runs in range were processed. 
%     This is done using a timer and can be accomplished without graphics
%     so it can be launched on a server with nohup/screen and keep working.
%     The results of each Analysis run is saved in a sub-folder to be
%     picked up later by this or other instances of Scanner. 
%     Use startCampaign() to run. 
%
% (b) Calculate the number of star-hours/coverage/number of detections from
%     the results of a previous campaign. It will go over all folders in
%     range and pick up the RunSummary objects. You can then use plotting
%     tools to show the statistics:
%
% (c) Scanning event candidates for occultations: use showNextCandidates()
%     to open a figure with the next batch of unclassified candidates.
%     This can be done with a different Scanner object, on an instance with
%     graphics, and the saved (classified) candidates should be saved to
%     file when done. 
%     
% Note: to choose the data folder Use chooseFolder(), chooseStartFolder() 
%       and chooseEndFolder() to define the runs for analysis (default is  
%       $DATA/WFAST/ and the range is all folders). 
    
    
    properties(Transient=true)
        
        a@img.Analysis;
        
    end
    
    properties % objects
        
        folders@trig.RunFolder; % output of trig.RunFolder.scan when getting a summary
        summaries@trig.RunSummary; % a list of summary objects for all runs in range
        
        candidates@trig.Candidate; % candidates loaded when running getNextCandidates() or showNextCandidates()
        
        timer; % timer used for continuous analysis
        
    end
    
    properties % inputs/outputs
        
    end
    
    properties % switches/controls
        
        root_folder = ''; % defaults to $DATA/WFAST
        date_start = ''; % default is Jan 1s of this year
        date_end = ''; % default is Dec 31st of this year
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Scanner(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.Scanner')
                if obj.debug_bit>1, fprintf('Scanner copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('Scanner constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
    end
    
    methods % getters
        
        function val = get.root_folder(obj)
            
            if isempty(obj.root_folder)
                t = datetime('now', 'TimeZone', 'UTC'); 
                val = fullfile(getenv('DATA'), sprintf('/WFAST/%04d', t.Year)); 
            else
                val = obj.root_folder; 
            end
            
        end
        
        function val = get.date_start(obj)
            
            if isempty(obj.date_start)
                t = datetime('now', 'TimeZone', 'UTC'); 
                val = sprintf('%04d-%02d-%02d', t.Year, 1, 1); 
            else
                val = obj.date_start;
            end
            
        end
        
        function val = get.date_end(obj)
            
            if isempty(obj.date_end)
                t = datetime('now', 'TimeZone', 'UTC'); 
                val = sprintf('%04d-%02d-%02d', t.Year, 12, 31); 
            else
                val = obj.date_end;
            end
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % API
        
        function startCampaign(obj)
            
            if obj.debug_bit, fprintf('Starting analysis campaign from %s to %s\n', obj.date_start, obj.date_end); end
            
            obj.setup_timer; 
            
        end
        
        function runAnalysis(obj)
            
            [success, report] = obj.callback_timer; 
            
            if success
                if obj.debug_bit, disp(report); end
            else
                error(report); 
            end
            
        end
        
        function chooseFolder(obj, folder)
            
            if nargin<2 || isempty(folder)
                folder = uigetdir(fullfile(getenv('DATA'), 'WFAST')); 
                if ~isempty(folder) && ischar(folder)
                    obj.root_folder = folder;
                end
            else
                obj.root_folder = folder; % what about relative paths?
            end
            
        end
        
        function chooseStartFolder(obj, date_or_folder)
            
            if nargin<2 || isempty(date_or_folder)
                folder = uigetdir(obj.root_folder); 
                if ~isempty(folder) && ischar(folder)
                    [~,folder] = fileparts(folder);
                    obj.date_start = folder(1:10); 
                end
            else
                obj.date_start = date_or_folder; % what about relative paths?
            end
            
        end
        
        function chooseEndFolder(obj, date_or_folder)
            
            if nargin<2 || isempty(date_or_folder)
                folder = uigetdir(obj.root_folder); 
                if ~isempty(folder) && ischar(folder)
                    [~,folder] = fileparts(folder);
                    obj.date_end = folder(1:10); 
                end
            else
                obj.date_end = date_or_folder; % what about relative paths?
            end
            
        end
        
        function getSummary(obj)
           
            % to be continued
            
        end
        
        function [success, report] = getNextCandidates(obj)
            
            
            % get the next folder that needs analysis
            r = trig.RunFolder.scan('folder', obj.root_folder, 'start', obj.date_start, ...
                'end', obj.date_end, 'next', 'unclassified');
            
            if isempty(r) || r.has_candidates==0
                success = 0; 
                report = 'Could not find a folder with unclassified candidates'; 
            else
                load(fullfile(r.folder, 'candidates.mat'))
                obj.candidates = cand; 
                success = 1; 
                report = sprintf('Found %d candidates in %s', length(obj.candidates), util.text.run_id(r.folder)); 
            end
            
        end
        
    end
    
    methods % internal utilities
        
        function setup_timer(obj, ~, ~)
            
            if ~isempty(obj.timer) && isa(obj.timer, 'timer') && isvalid(obj.timer)
                if strcmp(obj.timer.Running, 'on')
                    stop(obj.timer);
                    delete(obj.timer);
                    obj.timer = [];
                end
            end
            
            delete(timerfind('name', 'scanner-timer'));
            
            
            obj.timer = timer('BusyMode', 'drop', 'ExecutionMode', 'fixedRate', 'Name', 'scanner-timer', ...
                'Period', 300, 'StartDelay', 5, 'TimerFcn', @obj.callback_timer, 'ErrorFcn', @obj.setup_timer);
            
            start(obj.timer);
            
        end
        
        function [success, report] = callback_timer(obj, ~, ~)
            
            success = 0; 
            
            % get the next folder that needs analysis
            r = trig.RunFolder.scan('folder', obj.root_folder, 'start', obj.date_start, ...
                'end', obj.date_end, 'next', 'unprocessed');
            
            if isempty(r)
                report = 'Could not find a folder to run'; 
                return;
            end
            
            if isempty(obj.a)
                obj.a = img.Analysis;
            end
            
            worker_idx = obj.a.findWorkerUnread; % get a worker even if it was not read out
            
            if isempty(worker_idx) % if we could find a free worker
                report = 'Could not find a free worker!'; 
            else % we can run this folder now! 
                obj.a.reader.dir.cd(r.folder); 
                obj.a.reader.loadFiles; 
                obj.a.async_run('worker', worker_idx, 'reset', 1, 'logging', 1, 'save', 1); 
                run_id = util.text.run_id(obj.a.reader.current_dir); 
                report = sprintf('Started new run on worker %d for folder %s', worker_idx, run_id);
                success = 1; 
            end
            
        end
        
        function stop_timer(obj, ~, ~)
            
            stop(obj.timer); 
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

