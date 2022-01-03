classdef Folder < dynamicprops
% Use this class to scan data folders and get a summary of which folders 
% were analyzed and so on. 
% The main method to use is the static scan() method:
% >> r = run.Folder.scan(...); 
% This returns a vector of Folder objects, one for each run in the 
% root folder (by default scans <DATA>/WFAST/<year> for the current year).
%
% Optional Arguments to the scan() method:
% -folder: specify which folder to scan. This must be a root folder with 
%          subfolders for each date. See default above. 
% -next: Specify if you want all folders to be returned (leave empty, default)
%        or if you want only the next folder that is "unprocessed" or the 
%        next folder that is "unclassified". 
% -start_date and end_date: these define the start and end date (inclusive)
%                           for the scan. Give as string <YYYY-MM-DD> or as
%                           datetime objects. Leave empty to not limit. 
% -process_date: only analysis folders that were generated at or after this
%                date are considered. Older folders are obsolete and do not
%                make the run "processed". Default is defined in the static
%                function default_process_date(). 
% -num_files: only load runs that have at least this many imaging files. 
% -name: only load runs that match this run name exactly. 
% -regexp: only load runs with run names that match this regexp. 
% -glob: only load runs with run names that match this glob expression. 
% -catalog: load the catalog file if it exists in the run folder. Default 
%           is false (this slows down the scan!). 
% -debug_bit: control the verbosity of printouts. Default 0 is no output. 
%
% Thus there are two ways to use this scan. 
% (1) get all the folders with their properties and run summaries. To get 
%     this information leave the "next" argument empty (or don't give it). 
%     This is useful for e.g., summing the total star hours in a given time
%     range/observing season. 
% (2) get the next folder that needs to be processed or the next folder that
%     needs to have its candidates classified. Use the "unprocessed" or 
%     "unclassified" options to the "next" argument. 
%     This is useful for automatically processing all folders or for pulling
%     the unclassified candidates up for scanning. 
%
% Each Folder object will tell you the following information:
%   -folder: full path to the run. 
%   -identifier: each run has a unique identifier no matter what root folder
%                it is saved into. This is the last two folders of the full
%                path, the date and run name folders: e.g., 
%                "2020-06-01/ecliptic_run1". 
%   -analysis_folder: inside each run folder there are multiple analysis 
%                     folders, one for each date when we did analysis. This
%                     date is usually much later than the date the run was
%                     acquired. We always load only the most recent analysis
%                     folder. Folders that were processed before some given
%                     time (the "process_date") are ignored (assuming they
%                     were generated using obsolete analysis code). 
%   -analysis_date: a datetime object (date only) for the date when the 
%                   anaysis was performed. 
%   -num_files: how many HDF5 files are in this folder. 
%   -expT: the exposure time from the README text file (or from the HDF5 if 
%          the README was unparsable). This is usually either 0.03-0.04s, 
%          for fast mode files (that we can process for KBOs) or it will be 
%          longer exposures (typically 3 seconds) that are skipped by the 
%          KBO pipeline. 
%   -frame_rate: the frame rate, also from the README or header. This is
%                usually 25 Hz. 
%   -batch_size: how many frames per file. This is usually either 1 or 100.
%   -is_calibration: mark calibration folders (dark/flat). 
%   -is_full_frame: mark folders where full images were taken (not stack+
%                   cutouts, which are used for the KBO pipeline). 
%                   These are usually tests or images of asteroids etc. 
%   -was_processed: did we do any analysis of this run, after the minimal
%                   "process_date" defined in the scan. 
%   -has_summary: has a RunSummary object saved in a MAT-file. 
%   -has_candidates: has some Candidate objects saved in a MAT-file. 
%   -has_classified: has some Candidate object that have already been 
%                    scanned/vetted by a human and only the good ones saved
%                    separately into a candidates.mat file. 
%   -has_lightcurves_mat: if the processing also saved a lightcurves.mat 
%                         file with the entire lightcurves for the full run.
%                         This is only for older analysis folders. 
%   -has_lightcurves_hdf5: the processing also saved individual batches' 
%                          lightcurves into HDF5 files. There should be 
%                          as many such files as the original stack/cutout
%                          HDF5 files in the run. This is the new way to 
%                          save the photometric results. 
%   -has_microflares: only show folders where some micro flares were
%                     detected and saved to file. 
%   -isFastMode: method that tells which of a vector of Folder was 
%                recorded in fast-mode (i.e., expT<0.05). 
%   -readyForAnalysis: method that tells which of a vector of Folder 
%                      is ready for a new analysis run. The optional arg
%                      is used to define the minimal date before which 
%                      we consider analysis runs to be obsolete. 
%   -readyForClassifying: method that tells which of a vector of Folder 
%                         is ready for classifying/scanning. This means that
%                         there are candidates available but no classified
%                         file is available. 
%   -readyForGettingEvents: method that tells which of a vector of Folder
%                           has classified candidates (saved events) that
%                           should be loaded into a database of events. 
%   -summary: if the folder has an analysis folder (from the new pipeline)
%             then there should be a Summary object saved with the run. This
%             is loaded into the Folder object for easy access to the 
%             summary data, which includes runtime, star hours, header info
%             and so on. See the definition of Summary. 
% 
% 

    
    properties % objects
        
        summary@tno.Summary; 
        cat@head.Catalog; 
        head@head.Header;
        classifications@struct;
        
    end
    
    properties % inputs/outputs
       
        folder = ''; % full path to the folder where this run is saved (not including analysis folder!)
        identifier = ''; % short date and run name combination uniquely identifying this run (e.g., 2020-06-01/ecliptic_run1)
        analysis_folder = ''; % name of the analysis folder, if it exists (e.g., analysis_2020-09-10). Note the date is for when it was processed! 
        analysis_date = []; % datetime object for the time when the analysis was done (date only)
        
        num_files = 0; % how many HDF5 files are in this folder
        expT = NaN; % exposure time loaded from the text file [s]
        frame_rate = NaN; % frame rate loaded from the text file [Hz]
        batch_size = NaN; % how many frames in each batch/file
        RA_deg = NaN; % right ascention of field (degrees)
        Dec_deg = NaN; % declination of field (degrees)
        
        is_calibration = 0; % this is true if the run folder starts with "dark" or "flat"
        is_full_frame = 0; % this is true if the individual files are larger than 100 MB, meaning there are multiple full-frame images saved (not stack+cutout, not single image)
        
        was_processed = 0; % if this run was processed (not including processing dates earlier than the "process_time" given)
        has_summary = 0; % has a RunSummary object saved
        has_candidates = 0; % has some Candidates saved
        has_classified = 0; % has some classified Candidates saved (was scanned and classifications were given)
        
        has_lightcurves_mat = 0; % has a single, giant lightcurves MAT file (old method)
        has_lightcurves_hdf5 = 0; % has individual batches in lightcurve HDF5 files (new method)
        has_microflares = 0; % have there been any flare detections saved in this folder
    end
    
    properties % switches/controls
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        obs_hours; 
        
    end
    
    properties(Hidden=true)
       
        process_date; % test against this time to see if analysis folders are new enough (this is the minimal date given in the scan, not the date the analys was actually made!)
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj_vec = Folder(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'run.Folder')
                if obj_vec.debug_bit>1, fprintf('Folder copy-constructor v%4.2f\n', obj_vec(1).version); end
                obj_vec = util.oop.full_copy(varargin{1});
            else
                if obj_vec.debug_bit>1, fprintf('Folder constructor v%4.2f\n', obj.version); end
            end
            
        end
        
    end
    
    methods % getters
        
        function val = get.obs_hours(obj)
            
            if isempty(obj.frame_rate) || isempty(obj.batch_size) || isempty(obj.num_files)
                val = [];
            else
                val = obj.num_files.*obj.batch_size./obj.frame_rate./3600; 
            end
            
        end
        
        function val = isFastMode(obj_vec)
            
            val = false(size(obj_vec)); 
            
            for ii = 1:numel(obj_vec)
                
                val(ii) = ~isempty(obj_vec(ii).expT) && obj_vec(ii).expT<=0.15; 
                
            end
            
        end
        
        function val = isFullFrame(obj_vec)
            
            val = false(size(obj_vec)); 
            
            for ii = 1:numel(obj_vec)
                
                d = util.sys.WorkingDirectory(obj_vec(ii).folder);
                
                files = d.match('*.h5*'); 
                
                if isempty(files)
                    val = 0;
                else
                
                    file_details = dir(files{1}); 
                    
                    val(ii) = file_details.bytes > 5e8; % files with more than 0.5 GB are probably full-frame! 
                
                end
                
            end
            
        end
        
        function val = readyForAnalysis(obj_vec, minimal_date)
            
            if nargin<2 || ~isempty(minimal_date)
                minimal_date = datetime(obj_vec(1).start_date); % anything before this date is considered obsolete
            end
            
            fast = isFastMode(obj_vec); 
            
            val = fast & minimal_date<[obj_vec.analysis_date]; 
            
        end
        
        function val = readyForClassifying(obj_vec)
            
            val = [obj_vec.has_candidates] & [obj_vec.has_classified]==0; 
            
        end
        
        function val = readyForGettingEvents(obj_vec)
            
            val = [obj_vec.has_classified]; 
            
        end
        
        function val = getClassified(obj_vec)
            
            val = [];
            
            for obj = obj_vec
                L = load(fullfile(obj.folder, obj.analysis_folder, 'classified.mat')); 
                val = vertcat(val, L.candidates); 
            end
            
        end
        
    end
    
    methods % calculations
        
        function str_out = printout(obj)
            
            if obj.was_processed
                str = sprintf(' %-60s summary= %d | candidates= %d | classified= %d', ...
                    fullfile(obj.identifier, obj.analysis_folder), obj.has_summary, obj.has_candidates, obj.has_classified); 
            elseif obj.is_calibration
                str = sprintf(' %-60s calibration run...', obj.identifier); 
            elseif obj.is_full_frame
                str = sprintf(' %-60s full-frame run...', obj.identifier); 
            elseif obj.isFastMode==0
                str = sprintf(' %-60s slow-mode run...', obj.identifier); 
            else
                str = sprintf(' %-60s not processed after %s (fast= %d)', obj.identifier, obj.process_date, obj.isFastMode);
            end
            
            if nargout==0
                disp(str);
            else
                str_out = str;
            end
            
        end
        
        function rsync(obj_vec, varargin) % copy files for all these folders into destination
            
            % NOTE: for this to work on windows you must add this to the path:
            % D:\Dropbox (Weizmann Institute)\software\cwrsync\bin
            % (choose the dropbox folder from the actual machine)
            % for linux just make sure you have rsync on the path
            
            input =  util.text.InputVars;
            input.input_var('destination', ''); % root folder where run identifiers folders are kept (e.g., <Dropbox>/DATA/WFAST/2021)
            input.input_var('full', false); % if true, copy all content of of folders, if false copy only latest analysis folder
            input.input_var('dry_run', false); % don't copy anything, just show the command
            input.input_var('echo', false, 'verbose'); % print out the command output from each rsync call
            input.input_var('ignore_errors', false); % convert errors to warnings and keep going
            input.scan_vars(varargin{:}); 
            
            if isempty(input.destination)
                error('Must provide a destination root folder'); 
            end
            
            for ii = 1:length(obj_vec)
                
                obj = obj_vec(ii); 
                origin = obj.folder;
                target = fullfile(input.destination, obj.identifier);
                
                if ~exist(input.destination, 'dir')
                    error('Could not find the destination folder %s', input.destination);
                end
                
                if ~input.full % copy only the latest analysis folder
                    if isempty(obj.analysis_folder)
                        continue; 
                    else
                        origin = fullfile(origin, obj.analysis_folder); 
                    end
                end
                
                if length(dir(origin))<=2
                    if input.echo
                        disp('origin is an empty folder');
                    end
                    continue; 
                end
                
                if ~exist(target, 'dir')
                    if input.ignore_errors
                        warning('target "%s" does not exist!', target); 
                        continue; 
                    else
                        error('target "%s" does not exist!', target); 
                    end
                end
                
                if ispc
                    
                    origin = strrep(origin, '\', '/'); 
                    origin = strrep(origin, ':', ''); 
                    origin(1) = lower(origin(1)); 
                    origin = ['/cygdrive/', origin]; 
                    
                    target = strrep(target, '\', '/'); 
                    target = strrep(target, ':', '');
                    target(1) = lower(target(1)); 
                    target = ['/cygdrive/', target '/']; 
                    
                end
                
                comm = sprintf('rsync -vzza "%s" "%s"', origin, target);
                
                if input.dry_run
                    comm = [comm ' --dry-run']; 
                end
                
                [status, cmdout] = system(comm); 
                
                if input.echo
                    disp(cmdout);
                end
                
                if status
                    if input.ignore_errors
                        warning('Failed to copy %s:\n %s', origin, cmdout);
                    else
                        error('Failed to copy %s:\n %s', origin, cmdout);
                    end
                end
                
            end
            
        end
                
        function getCandidateSummary(obj) % load classifications summary text file

            obj.classifications = struct('all', [], 'real', [], 'sim', []); 

            cls = tno.Candidate.getListOfClasses; 

            for ii = 1:length(cls)
                cls{ii} = strrep(cls{ii}, ' ', '_'); 
            end
            new_cell = cell(1,length(cls));
            new_cell(:) = {0}; % put zeros everywhere
            cls_arr = [cls;new_cell]; 
            cls_struct = struct(cls_arr{:}); % each classname is followed by a zero
            obj.classifications.all = cls_struct;
            obj.classifications.real = cls_struct;
            obj.classifications.sim = cls_struct;

            if exist(fullfile(obj.folder, obj.analysis_folder, 'classified.txt'), 'file') % read classifications from file
                
                % read the text file
                fid = fopen(fullfile(obj.folder, obj.analysis_folder, 'classified.txt'));
                file_close = onCleanup(@() fclose(fid)); % make sure to close the file at the end
                
                for ii = 1:10000

                    line = fgetl(fid); 

                    if isnumeric(line)
                        break;
                    end
                    
                    line = strip(line); 
                    if ~isempty(line)
                        
                        idx = regexp(line, ':', 'once'); % find the colon
                        if ~isempty(idx) && length(line) > idx
                            this_class = strip(line(1:idx-1));
                            this_class = strrep(this_class, ' ', '_'); 
                            if ismember(this_class, cls)
                                numbers = str2num(line(idx+1:end));
                                if ~isempty(numbers)
                                    
                                    obj.classifications.all.(this_class) = numbers(1); % first number is all events
                                    
                                    if isscalar(numbers) % no simulated events here
                                        obj.classifications.real.(this_class) = numbers(1); 
                                        obj.classifications.sim.(this_class) = 0; 
                                    else % get the number of real vs. simualated
                                        obj.classifications.real.(this_class) = numbers(1)-numbers(2); 
                                        obj.classifications.sim.(this_class) = numbers(2); 
                                    end
                                    
                                end
                                
                            end
                        end
                        
                    end
                    
                end
                
            elseif exist(fullfile(obj.folder, obj.analysis_folder, 'classified.mat'), 'file')
                
                load(fullfile(obj.folder, obj.analysis_folder, 'classified.mat'));
                
                if ~isempty(candidates)
                    
                    for ii = 1:length(candidates)
                        
                        this_class = candidates(ii).classification;
                        this_class = strrep(this_class, ' ', '_'); 
                        obj.classifications.all.(this_class) = obj.classifications.all.(this_class) + 1;
                        if candidates(ii).is_simulated
                            obj.classifications.sim.(this_class) = obj.classifications.sim.(this_class) + 1;
                        else
                            obj.classifications.real.(this_class) = obj.classifications.real.(this_class) + 1;
                        end
                        
                    end
                    
                end
                
                % write a text file for future reference
                fid = fopen(fullfile(obj.folder, obj.analysis_folder, 'classified.txt'), 'wt');
                file_close = onCleanup(@() fclose(fid)); % make sure to close the file at the end
                
                for ii = 1:length(cls)
                    
                    fprintf(fid, '%s : %d', strrep(cls{ii}, '_', ' '), obj.classifications.all.(cls{ii})); 
                    
                    if obj.classifications.sim.(cls{ii})
                        fprintf(fid, ' (%d)', obj.classifications.sim.(cls{ii}));
                    end
                    
                    fprintf(fid, '\n');
                    
                end
                
            else
                error('Could not find any type of classified candidates file!'); 
            end
            
        end
        
    end
    
    methods(Static=true) % scan method lives here!
            
        function obj_vec = scan(varargin) % go over data folders and return a vector of Folder objects
        % Usage: obj_vec = run.Folder.scan(varargin)
        % Goes over all folders in a root data folder and extracts the directory
        % and some additional info on all the relevant folders. 
        %
        % 
        %
        % OPTIONAL ARGUMENTS:
        % -folder: specify which folder to scan. This must be a root folder with 
        %          subfolders for each date. See default above. 
        % -next: Specify if you want all folders to be returned (leave empty, default)
        %        or if you want only the next folder that is "unprocessed" or the 
        %        next folder that is "unclassified". 
        % -start_date and end_date: these define the start and end date (inclusive)
        %                           for the scan. Give as string <YYYY-MM-DD> or as
        %                           datetime objects. Leave empty to not limit. 
        % -process_date: only analysis folders that were generated at or after this
        %                date are considered. Older folders are obsolete and do not
        %                make the run "processed". Default is defined in the static
        %                function default_process_date(). 
        % -num_files: only load runs that have at least this many imaging files. 
        % -name: only load runs that match this run name exactly. 
        % -regexp: only load runs with run names that match this regexp. 
        % -glob: only load runs with run names that match this glob expression. 
        % -catalog: load the catalog file if it exists in the run folder. Default 
        %           is false (this slows down the scan!). 
        % -debug_bit: control the verbosity of printouts. Default 0 is no output. 
        %
        % There are two ways to use this scan. 
        % (1) get all the folders with their properties and run summaries. To get 
        %     this information leave the "next" argument empty (or don't give it). 
        %     This is useful for e.g., summing the total star hours in a given time
        %     range/observing season. 
        % (2) get the next folder that need to be processed or the next folder that
        %     needs to have its candidates classified. Use the "unprocessed" or 
        %     "unclassified" options to the "next" argument. 
        %     This is useful for automatically processing all folders or for pulling
        %     the unclassified candidates up for scanning. 
        
            import util.text.cs;
            
            input = util.text.InputVars;
            input.input_var('folder', '', 'root'); % root folder to scan (must have date folders inside)
            input.input_var('next', []); % leave empty for all folders, use "unprocessed" or "unclassified" to get the next folder that needs processing or classification (respectively)
            input.input_var('start_date', []); % scan starting from runs taken on this date (inclusive). Can be string <YYYY-MM-DD> or datetime object
            input.input_var('end_date', []); % scan upto runs taken on this date (inclusive). Can be string <YYYY-MM-DD> or datetime object
            input.input_var('process_date', run.Folder.default_process_date); % only consider as processed runs that were processed on or after this date
            input.input_var('files', 100, 'min_files', 'minimal_files', 'min_num_files', 'minimal_number_files', 'num_files', 'number_files'); % must have this many HDF5 files for analysis
            input.input_var('name', '', 'run_name'); % match only run names that are equal to this
            input.input_var('regexp', '', 'regular_expression'); % match run names to this regular expression
            input.input_var('glob', '', 'glob_expression', 'wildcard'); % match run names to this wildcard (glob) expression
            input.input_var('catalog', false); % pull the catalog file into each found object
            input.input_var('header', false); % save the full header for each folder 
            input.input_var('debug_bit', 0); % verbosity of printouts
            input.scan_vars(varargin{:}); 
            
            % use only one of the run-name matches 
            if (~isempty(input.name) && ~isempty(input.regexp)) ||...
                    (~isempty(input.name) && ~isempty(input.glob)) || ... 
                    (~isempty(input.regexp) && ~isempty(input.glob))
                
                error('Must supply only one of these optional arguments: "name", "regexp", "glob". ');                 
            end
            
            if ischar(input.start_date)
                input.start_date = datetime(input.start_date); % convert to datetime
            end
            
            if ischar(input.end_date)
                input.end_date = datetime(input.end_date); % convert to datetime
            end
            
            if ischar(input.process_date)
                input.process_date = datetime(input.process_date); % convert to datetime
            end
            
            % make sure start and end date are compatible
            if ~isempty(input.start_date) && ~isempty(input.end_date) && input.start_date>input.end_date
                error('Start date %s is after end date %s', input.start_date, input.end_date); 
            end
            
            %%%%%%%%%%%%%%%% get the root folder %%%%%%%%%%%%%%%%%%%%%
            
            if isempty(input.folder) % default is to load the relevant year folder from the dropbox (must define environmental DATA)
                y = run.Folder.guess_year(input.start_date, input.end_date); 
                d = util.sys.WorkingDirectory(fullfile(getenv('DATA'), sprintf('WFAST/%d',y))); 
            elseif ischar(input.folder) % a string input can mean a few things:

                if ismember(input.folder, {'DATA', 'DATA_TEMP', 'DATA_EXTRA'}) % must define a year also! 
                    d = getenv(folder); 
                    d = util.sys.WorkingDirectory(d); 
                    y = run.Folder.guess_year(input.start_date, input.end_date); 
                    d.cd(sprintf('%d', y));                
                elseif exist(input.folder, 'dir') % just get the full path to the folder
                    d = util.sys.WorkingDirectory(input.folder);
                else
                    error('Could not find the folder "%s". ', input.folder); 
                end

            elseif isa(input.folder, 'util.sys.WorkingDirectory') % give the folder in the form of a WorkingDirectory object
                d = input.folder;
            else
                error('Wrong input to Folder.scan() with class "%s".', class(input.folder)); 
            end
            
            if ~exist(d.pwd, 'dir')
                error('Could not find the folder "%s". ', d.pwd); 
            end
            
            %%%%%%%%%%%%% start scanning %%%%%%%%%%%%%%%%%
            
            obj_vec = run.Folder.empty; % this output is filled out in the loop
            
            list = d.dir; % list of sub-folders (should be date folders)
            root = d.pwd; % the root folder where we come back every time. 
            
            for ii = 1:length(list) % go over dates
                
                if ~regexp(list{ii}, '\d{4}-\d{2}-\d{2}')
                    continue; % skip folders that are not date-folders
                end
                
                if length(list{ii})<10
                    continue;
                end
                
                this_date = datetime(list{ii}(1:10)); % translate the date folder to a datetime object
                
                if ~isempty(input.start_date) % limits on the earliest dates
                    if input.start_date>this_date
                        continue;
                    end
                end
                
                if ~isempty(input.end_date) % limits on the latest dates
                    if input.end_date<this_date
                        continue;
                    end
                end
                
                d.cd(root); % go back to root, then move down
                d.cd(list{ii}); % move into the date folder
                
                if input.debug_bit, fprintf('folder: "%s"\n', d.pwd); end
                
                run_folders = d.dir; % if this is empty, we will just continue to the next date folder
                
                for jj = 1:length(run_folders) % go over run folders
                    
                    d.cd(root); % go back to root (this makes sure we don't lose our place if there is an error in one loop iteration)
                    d.cd(list{ii}); 
                    d.cd(run_folders{jj}); 
                    
                    idx = regexp(run_folders{jj}, '_run\d+$'); 
                    
                    if isempty(idx)
                        run_name = run_folders{jj};
                    else
                        run_name = run_folders{jj}(1:idx-1); 
                    end
                    
                    if ~isempty(input.name)
                        if ~strcmp(input.name, run_name), continue; end
                    elseif ~isempty(input.regexp)
                        if isempty(regexp(run_name, input.regexp, 'once')), continue; end
                    elseif ~isempty(input.glob)
                        if isempty(regexp(run_name, regexptranslate('wildcard', input.glob), 'once')), continue; end
                    end
                    
                    files = d.match('*.h5*'); % cell array with all HDF5 files in this folder
                    
                    if length(files) < input.files
                        continue;
                    end
                    
                    new_obj = run.Folder; % new object added to the list, for this specific run
                    
                    % fill the identification info for this object
                    new_obj.folder = d.pwd;
                    new_obj.identifier = d.two_tail; 
                    
                    if ismember(new_obj.identifier(end), {'\', '/'}) % remove trailing slashes from the identifier
                        new_obj.identifier = new_obj.identifier(1:end-1);
                    end
                    
                    new_obj.num_files = length(files); 
                    
                    if isempty(files) % check if it is an empty folder
                        % do nothing? 
                    elseif regexp(run_folders{jj}, '^(dark|flat).*') % check if it is calibration
                        new_obj.is_calibration = 1; 
                        new_obj.is_full_frame = 1; 
                    elseif new_obj.isFullFrame % check if the files contain multiple full-frame images (as opposed to single images or stack+cutout files)
                        new_obj.is_full_frame = 1; 
                    else % not calibration and not full-frame
                        
                        %%%%%%%% get the parameters %%%%%%%%%%
                    
                        if exist(fullfile(d.pwd, 'A_README.txt'), 'file') % if there is a text file at all... 
                            [new_obj.expT, new_obj.frame_rate, new_obj.batch_size, new_obj.RA_deg, new_obj.Dec_deg] = new_obj.getParameters(fullfile(d.pwd, 'A_README.txt')); 
                        end
                        
                        if isempty(new_obj.expT) || isempty(new_obj.frame_rate) || isempty(new_obj.batch_size) || ...
                            isempty(new_obj.RA_deg) || isempty(new_obj.Dec_deg) % no result, try getting the expT by loading the header

                            new_obj.head = new_obj.getHeader(files); % try to find the header from the README or HDF5 file

                            if ~isempty(new_obj.head)
                                new_obj.expT = new_obj.head.EXPTIME; 
                                new_obj.frame_rate = new_obj.head.FRAMERATE;
                                new_obj.batch_size = new_obj.head.NAXIS3; 
                                new_obj.RA_deg = new_obj.head.RA_DEG; 
                                new_obj.Dec_deg = new_obj.head.DEC_DEG; 
                            end
                            
                        end
                        
                        if isempty(new_obj.expT), new_obj.expT = NaN; end
                        if isempty(new_obj.frame_rate), new_obj.frame_rate = NaN; end
                        if isempty(new_obj.batch_size), new_obj.batch_size = NaN; end
                        if isempty(new_obj.RA_deg), new_obj.RA_deg = NaN; end
                        if isempty(new_obj.Dec_deg), new_obj.Dec_deg = NaN; end
                        
                        %%%%%%%% load the catalog file
                        
                        if input.catalog
                            
                            if exist(fullfile(d.pwd,'catalog.mat'), 'file')
                                new_obj.cat = head.Catalog; 
                                new_obj.cat.loadMAT(fullfile(d.pwd,'catalog.mat')); 
                            end
                                
                        end
                        
                        %%%%%%% load the header if it wasn't loaded already
                        
                        
                        if input.header
                            if isempty(new_obj.head)
                                new_obj.head = new_obj.getHeader(files); % try to find the header from the README or HDF5 file
                            end
                        else
                            new_obj.head = head.Header.empty; % get rid of header if we were not asked to fetch it
                        end
                        
                        new_obj.has_microflares = exist(fullfile(d.pwd, 'micro_flares.mat'), 'file') > 0; 
                        
                        %%%%%%%% go into the analysis folders %%%%%%%%%

                        analysis_folders = sort(d.match_folders('analysis_*')); % get the folders in ascending time order
                        
                        % no analysis has been found for this folder yet
                        new_obj.analysis_folder = ''; 
                        new_obj.analysis_date = ''; 
                        new_obj.was_processed = false;
                        new_obj.has_summary = false;
                        new_obj.has_candidates = false;
                        new_obj.has_classified = false; 
                        new_obj.has_lightcurves_mat = false;
                        new_obj.has_lightcurves_hdf5 = false;
                        
                        if ~isempty(analysis_folders) % we found an analysis folder! 

                            [~, new_obj.analysis_folder] = fileparts(analysis_folders{end});
                            
                            new_obj.analysis_date = datetime(new_obj.analysis_folder(10:end)); % chop off the word "analysis_" and turn it into a datetime

                            new_obj.process_date = input.process_date; % remember the limiting date for considering analysis folders (not the actual date when it was analyized!)
                            
                            if new_obj.analysis_date>=input.process_date % there is an analysis folder that is up-to-date enough to use it:
                                
                                new_obj.was_processed = true; 
                                
                                d.cd(analysis_folders{end}); % now we are in the most recent analysis folder! 

                                if exist(fullfile(d.pwd, 'summary.mat'), 'file') % load the RunSummary from the MAT-file
                                    try 
                                        load(fullfile(d.pwd, 'summary.mat'));                             
                                        obj.summary = summary; 
                                        new_obj.has_summary = true;
                                    catch ME
                                        warning(ME.getReport); 
                                    end
                                end

                                if exist(fullfile(d.pwd, 'candidates.mat'), 'file') % candidate file exists
                                    new_obj.has_candidates = true;
                                end

                                if exist(fullfile(d.pwd, 'classified.mat'), 'file') % load the classified Candidates from the MAT-file
                                    new_obj.has_classified = true;
                                    % load the Candidate objects from the
                                    % MAT-file or from summary text file
                                    new_obj.getCandidateSummary; 
                                end

                                if exist(fullfile(d.pwd, 'lightcurves.mat'), 'file') % check if there is a lightcurves MAT file
                                    new_obj.has_lightcurves_mat = true;
                                end

                                LCs = d.match('*Lightcurves.h5*'); % check if there are lightcurve HDF5 files

                                if ~isempty(LCs)
                                    new_obj.has_lightcurves_hdf5 = true;
                                end

                            end % processing time
                            
                        end % has analysis folders

                    end  % check if there are files and if the folder is for calibration
                    
                    if input.debug_bit, new_obj.printout; end

                    if isempty(input.next) % the default is to get back all the run folders
                        obj_vec(end+1) = new_obj; 
                    elseif cs(input.next, 'unprocessed') % only interested in the first instance that is unprocessed
                        if new_obj.isFastMode && new_obj.was_processed==0 && new_obj.num_files>=input.files
                            obj_vec = new_obj;
                            return;
                        end                            
                    elseif cs(input.next, 'unclassified') % only interested in the first instance that has not yet been classified
                        if new_obj.was_processed && new_obj.has_candidates && ~new_obj.has_classified
                            obj_vec = new_obj;
                            return; 
                        end
                    end

                end % go over run folders
                
            end % go over date folders
            
        end
        
        function val = guess_year(start_date, end_date)
            
            if isempty(start_date) && isempty(end_date)
                val = year(datetime('now'));
            elseif isempty(start_date) % end date is NOT empty
                val = year(end_date);
            elseif isempty(end_date) % start date is NOT empty
                val = year(start_date);
            else % both dates are given! 
                y1 = year(start_date);
                y2 = year(end_date);
                if y1~=y2
                    warning('Start year (%d) and end year (%d) are not the same!', y1, y2); 
                end
                val = y1; % prefer start year
            end
            
        end
        
        function [expT, frame_rate, batch_size, RA, DE] = getParameters(filename) % scan a text file until finding all these parameters
            
            expT = []; % if we can't find the exposure time, return empty
            frame_rate = [];
            batch_size = [];
            RA = [];
            DE = [];
            
            try 
                
                fid = fopen(filename); 
                on_cleanup = onCleanup(@() fclose(fid)); 

                for ii = 1:10000

                    line = fgetl(fid); 

                    if isnumeric(line)
                        break;
                    end

                    if isempty(expT) 
                        idx = regexp(line, '(expT|EXPTIME)', 'end');

                        if ~isempty(idx)
                            expT = util.text.parse_value(line(idx+2:end));       
                            continue;
                        end
                    end
                    
                    if isempty(batch_size)
                        idx = regexp(line, 'NAXIS3', 'end');

                        if ~isempty(idx)
                            batch_size = util.text.parse_value(line(idx+2:end));  
                            continue;
                        end
                    end
                    
                    if isempty(frame_rate)
                        idx = regexp(line, 'FRAMERATE', 'end');

                        if ~isempty(idx)
                            frame_rate = util.text.parse_value(line(idx+2:end));   
                            continue;
                        end
                    end
                    
                    if isempty(RA)
                        idx = regexp(line, 'RA_DEG', 'end');

                        if ~isempty(idx)
                            frame_rate = util.text.parse_value(line(idx+2:end));   
                            continue;
                        end
                    end
                    
                    if isempty(DE)
                        idx = regexp(line, 'DEC_DEG', 'end');

                        if ~isempty(idx)
                            frame_rate = util.text.parse_value(line(idx+2:end));   
                            continue;
                        end
                    end
                    
                    if ~isempty(expT) && ~isempty(frame_rate) && ...
                            ~isempty(batch_size) && ~isempty(RA) && ~isempty(DE)
                        break;
                    end
                        
                end
                
            catch ME
                warning(ME.getReport); 
            end
            
        end
        
        function h = getHeader(files)
            
            h = [];
            
            folder = fileparts(files{1}); 
            
            locations = {'/acquisition/head', '/acquisition/pars', '/camera/head', '/camera/pars'}; 
            
            for kk = 1:length(locations)
                try 
                    h = util.oop.load(fullfile(folder, 'A_README.txt'), 'location', locations{kk}); 
                    break;
                end
            end
            
            if isempty(h)
                try
                    h = util.oop.load(files{1}, 'location', '/header');
                catch 
                    h = util.oop.load(files{1}, 'location', '/pars');
                end

            end

        end
        
        function val = default_process_date % the default minimal processing date for the new pipeline
            val = '2020-10-15'; 
        end
         
    end
    
end

