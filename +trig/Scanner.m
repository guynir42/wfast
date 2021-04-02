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
        
        overview@trig.Overview;
        
        timer; % timer used for continuous analysis
        
    end
    
    properties % inputs/outputs
        
    end
    
    properties % switches/controls
        
        root_folder = ''; % defaults to $DATA/WFAST
        date_start = ''; % default is Jan 1s of this year
        date_end = ''; % default is Dec 31st of this year
        date_process = ''; % analysis folders before this date are ignored (defaults to today)
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
                obj.makeAnalysisObject; 
            end
            
        end
        
        function makeAnalysisObject(obj)
            
            obj.a = img.Analysis;
            obj.a.use_save_batched_lightcurves = 0; 

            % TODO: make sure this object has all the correct
            % parametrers, e.g., not to save lightcurves... 
            
        end
        
    end
    
    methods % reset/clear
        
    end
    
    methods % getters
        
        function val = get.root_folder(obj)
            
            if isempty(obj.root_folder)
                y = trig.RunFolder.guess_year(obj.date_start, obj.date_end);
%                 t = datetime('now', 'TimeZone', 'UTC'); 
                val = fullfile(getenv('DATA'), sprintf('/WFAST/%04d', y)); 
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
        
        function all_runs = calcOverview(obj, varargin)
           
            input = util.text.InputVars;
            input.input_var('runs', []); % optionally input the runs from previous calculations...
            input.input_var('classified', false);
            input.scan_vars(varargin{:}); 
            
            if isempty(obj.overview)
                obj.overview = trig.Overview;
            else
                obj.overview.reset;
            end
            
            if isempty(input.runs)  
                t0 = tic;
                all_runs = trig.RunFolder.scan('folder', obj.root_folder, 'start', obj.date_start, ...
                    'end', obj.date_end, 'next', [], 'process_date', obj.date_process); % get all folders
                if obj.debug_bit, fprintf('Time to load run folders is %s\n', util.text.secs2hms(toc(t0))); end
            else
                all_runs = input.runs; 
            end
            
            all_runs = all_runs(logical([all_runs.has_summary])); 
            all_runs = all_runs(logical([all_runs.has_candidates])); 
            
            if input.classified
                all_runs = all_runs(logical([all_runs.has_classifieds])); 
            end
            
            obj.overview.folders = all_runs; 
            
            for ii = 1:length(all_runs)
                
                if obj.debug_bit, fprintf('ii= %d / %d. Loading summary from folder: %s\n', ii, length(all_runs), all_runs(ii).folder); end
                
                try
                    all_runs(ii).loadSummary;
%                     L = load(fullfile(all_runs(ii).folder, all_runs(ii).analysis_folder, 'summary.mat'));
                    obj.overview.input(all_runs(ii).summary); 
                catch ME
                    warning(ME.getReport); 
                end
                
            end
            
        end
        
        function [success, report] = getNextCandidates(obj)
            
            % get the next folder that needs analysis
            r = trig.RunFolder.scan('folder', obj.root_folder, 'start', obj.date_start, ...
                'end', obj.date_end, 'next', 'unclassified', 'process_date', obj.date_process);
            
            if isempty(r) || r.has_candidates==0
                success = 0; 
                report = 'Could not find a folder with unclassified candidates'; 
            else
                
                load(fullfile(r.folder, r.analysis_folder, 'candidates.mat'))
                
                if isempty(cand)
                    
                    obj.candidates = cand;
                    
                    candidates = cand; % the variable name needs to match Candidate.saveClassified()
                    
                    filename = fullfile(r.folder, r.analysis_folder, 'classified.mat');
                    
                    save(filename, 'candidates', '-v7.3');
                    
                    success = 1; 
                    report = sprintf('Found no candidates. Saving empty classified file. '); 
                    
                else
                
                    for ii = 1:length(cand)
                        cand(ii).folder = fullfile(r.folder, r.analysis_folder); % make sure to update each candidate to know what folder it was loaded from! 
                    end

                    obj.candidates = cand; 
                    success = 1; 
                    report = sprintf('Found %d candidates.', length(obj.candidates)); 

                end
                
            end
            
        end
        
        function cand = collectOccultations(obj)
            
            if ~isempty(obj.overview)
                obj.calcOverview;
            end
            
            cand = trig.Candidate.empty; 
            
            for ii = 1:length(obj.overview.folders)
                
                f = fullfile(obj.overview.folders(ii).folder, ...
                       obj.overview.folders(ii).analysis_folder, 'classified.mat');
                
                if exist(f, 'file')
                    L = load(f);

                    new_cand = L.candidates([L.candidates.is_simulated]==0); 

                    cand = vertcat(cand, new_cand); 
                end
                
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
            
            if isempty(obj.date_process)
                obj.date_process = datestr(t, 'yyyy-mm-dd', datetime('today')); 
            end
            
            if isempty(obj.a)
                obj.makeAnalysisObject;
            end
            
            for ii = 1:length(obj.a.futures)
                
                if isa(obj.a.futures{ii}, 'parallel.FevalFuture') && ...
                        strcmp(obj.a.futures{ii}.State, 'finished') && ...
                        obj.a.futures{ii}.Read==0 && ...
                        isempty(obj.a.futures{ii}.Error) % unread, finished runs
                    
                    if ~isempty(obj.a.futures_analysis_folder{ii}) % we know the analysis folder, we can add a diary file to it
                        
                        fid = fopen(fullfile(obj.a.futures_analysis_folder{ii}, 'diary.txt'), 'at');
                        
                        if fid>0
                            on_cleanup = onCleanup(@() fclose(fid));
                            fprintf(fid, '%s', obj.a.futures{ii}.Diary); 
                            if ~isempty(obj.a.futures{ii}.Error)
                                fprintf(fid, '%s', obj.a.futures{ii}.Error.getReport('extended','hyperlinks','off'));
                            end
                        end
                        
                    end

                    obj.a.futures{ii}.fetchOutputs; % make this unread (we can also dump this into a variable, containing a copy of the Analysis object)
                    
                end
                
            end
            
            % get the next folder that needs analysis
            r = trig.RunFolder.scan('folder', obj.root_folder, 'start', obj.date_start, ...
                'end', obj.date_end, 'next', 'unprocessed', 'process_date', obj.date_process);
            
            if isempty(r)
                report = 'Could not find a folder to run'; 
                return;
            end
            
            worker_idx = obj.a.findWorkerUnread; % get a worker even if it was not read out
%             worker_idx = obj.a.findWorker; % get a worker even if it was not read out
            
            if isempty(worker_idx) % if we could find a free worker
                report = 'Could not find a free worker!'; 
            else % we can run this folder now! 
                obj.a.reader.dir.cd(r.folder); 
                obj.a.reader.loadFiles; 
                obj.a.async_run('worker', worker_idx, 'reset', 1, 'logging', 1, 'save', 1); 
                run_id = util.text.run_id(obj.a.reader.current_dir); 
                report = sprintf('Started new run on worker %d for folder %s', worker_idx, run_id);
                if obj.debug_bit, fprintf('%s: %s\n', datetime('now', 'TimeZone', 'UTC'), report); end
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

