classdef (CaseInsensitiveProperties) Target < handle
% Contains all the info that is read from a target list file, and can make
% informed decisions on which field should be observed at each given time. 
% 
% Each Target object contains a head.Ephemeris object that does most of the 
% heavy lifting in terms of figuring out airmass and other constraints. 
% The idea is to generate a list of Target objects and update them to the 
% current time (or project them into the future). Then find the best target
% based on priority (in this object) and airmass (in the Ephemeris object). 
%
% Use the static readFile(filename='target_list.txt') function to populate
% a vector of Target objects. Use parse(str) to turn a string formatted as:
% "object_name, RA_Dec_str, key1=val1, key2=val2, ..." into a the Target 
% object, by filling its own properties, the Ephemeris object is has, and 
% that object's "constraints" object. 
%
% Use observable{varargin) to check if the target can be observed based on 
% the internal constraints and time of the Ephemeris object, all of which 
% can be temporarily overriden by the varargin pairs. 
% 
% Use the best_target(obj_vec, varargin) on a vector of Target objects to 
% find the best one, out of those that are observable. Here too we can use
% the varargin to override some constraints or define a different time. 
% Note that giving a "time" parameter would also change the internal state
% of the Ephemeris objects, so make sure to always update them to the right 
% time. Use "now" to get them to current time. 
% This update to the new time will also move around the dynamically allocated
% fields like "ecliptic" or "moon". This is what we want from a scheduler! 
% If best_target() can not find any observable target, it returns an empty
% Target vector. 

    properties(Transient=true)
        
    end
    
    properties % objects
        
        ephem@head.Ephemeris;
        
    end
    
    properties % inputs/outputs
        
        obs_history = []; % struct array with a single struct for each separate run on this target, containing RA_deg, Dec_deg, and start_time and end_times
        start_time = ''; % keep the start time of the current observation (in HH:MM format). 
        start_side = '';
        start_RA_deg = [];
        start_Dec_deg = [];
        
        list_table = []; % a table with survey field coordinates and other details like number of stars, how many hours it was observed, etc.
%         last_list_update = NaT; % datetime with the last time the list file was updated from memory
        
    end
    
    properties % switches/controls
        
        index = [];
        
        use_resolver = 0; % this is true only if we are not given coordinates explicitely!
        
        % scheduling priority and timing
        priority = 1; % higher is more important (fractions are also acceptible). 
        decay_rate = []; % how fast priority decays per hour of observing time.
        
        % instructions to mount/camera
        tracking = 1; % most targets need the telescope to be tracking 
        cam_mode = 'fast'; % choose "fast" for KBO type survey at 25Hz, or "slow" for deeper targets with 3 second exposures
        exp_time = []; % change this only when using "slow" camera mode. Should not be larger than 5 or 10 seconds for guiding to work. 
        
        % additional info
        list_filename = ''; % should we load the target from a list in a text file? 
        sort_columns = {'TotalHours', 'NumStars*'}; % sort first on what was observed least, then on what has the most stars
        repeat_max = 2; % do not observe fields more than two nights in a row
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        % target name and coordinates 
        name; % must give a name for the field/object
        RA; % numeric hours or sexagesimal hours
        Dec; % numeric degrees or sexagesimal degrees
        
        alt;
        airmass;
        side;
        moon; 
        
    end
    
    properties(Hidden=true)
       
        list_filename_loaded = ''; 
        
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = Target(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.sched.Target')
                if obj.debug_bit>1, fprintf('Target copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                
                if obj.debug_bit>1, fprintf('Target constructor v%4.2f\n', obj.version); end
                
                obj.ephem = head.Ephemeris;
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj) % returns the object to what it was when initialized
            
            obj.ephem.reset;
            
            obj.use_resolver = 0; 
            
            % scheduling priority and timing
            obj.priority = 1; % higher is more important (fractions are also acceptible). 
            obj.decay_rate = []; % how fast priority decays per hour of observing time.
            
            % instructions to mount/camera
            obj.tracking = 1; % most targets need the telescope to be tracking 
            obj.cam_mode = 'fast'; % choose "fast" for KBO type survey at 25Hz, or "slow" for deeper targets with 3 second exposures
            obj.exp_time = []; % change this only when using "slow" camera mode. Should not be larger than 5 or 10 seconds for guiding to work. 
            
            obj.ephem.reset;
            
            obj.clear;
            
        end
        
        function clear(obj) % return the target to the state before it was observed at all
            
            obj.ephem.clear;
            
            obj.obs_history = []; 
            obj.start_time = ''; 
            obj.start_RA_deg = [];
            obj.start_Dec_deg = [];
        
        end
        
    end
    
    methods % getters
        
        function val = get.name(obj)
            
            val = obj.ephem.name;
            
        end
        
        function val = get.RA(obj)
            
            val = obj.ephem.RA;
            
        end
        
        function val = get.Dec(obj)
            
            val = obj.ephem.Dec;
            
        end
        
        function val = getCurrentPriority(obj)
            
            if isempty(obj.decay_rate)
                val = obj.priority;
            else
                val = obj.priority - obj.ephem.getTotalRuntimeMinutes/60*obj.decay_rate; 
            end
            
        end
        
        function val = get.alt(obj)
            
            val = obj.ephem.Alt_deg;
            
        end
        
        function val = get.airmass(obj)
            
            val = obj.ephem.AIRMASS;
            
        end
        
        function val = get.side(obj)
            
            if isempty(obj.ephem.HA_deg) || isnan(obj.ephem.HA_deg)
                val = '';
            elseif obj.ephem.HA_deg<0
                val = 'East';
            elseif obj.ephem.HA_deg>=0
                val = 'West';
            end
            
        end
        
        function val = get.moon(obj)
            
            val = obj.ephem.moon_dist; 
            
        end
        
        function val = summary(obj)
            
%             val = sprintf('%s at %s%s with HA= %s and airmass= %f', obj.name, obj.RA, obj.Dec, obj.ephem.HA, obj.ephem.airmass); 
            val = sprintf('%s at %s%s', obj.name, obj.RA, obj.Dec); 
            
        end
        
        function val = details(obj)
            
            val = sprintf('HA= %11s (%s) | Alt= %d | A.M.= %4.2f | moon= %3d | ecl/gal= %2d/%2d', ...
                obj.ephem.HA, obj.side, round(obj.ephem.Alt), obj.airmass, obj.moon, round(obj.ephem.ECL_lat), round(obj.ephem.GAL_lat)); 
            
        end
        
    end
    
    methods % setters
        
        function set.name(obj, val)
            
            obj.ephem.name = val;
            
        end
        
        function set.RA(obj, val)
            
            obj.ephem.RA = val;
            
        end
        
        function set.Dec(obj, val)
            
            obj.ephem.Dec = val;
            
        end
        
    end
    
    methods % calculations
        
        function parse(obj, str, use_resolve) % get target parameters from a line of text
        % Usage: parseText(obj, str, use_resolve=0)
        % The basic target syntax is kept simple based on comma separators. 
        % Each line MUST begin with the name of the object/field. 
        % This is followed by coordinates, e.g., 18:30:00.3+22:40:52.8. 
        % Coordinates must be given as RA (hours) and Dec (degrees) with 
        % a plus or minus separating them. 
        % Second argument is to use resolver after loading (default false):
        % If coordinates are left empty, the object name is given to the 
        % head.Ephemeris object's name interpreter. 
        %
        % The remaining arguments are given as "keyword=value", each pair
        % separated using commas. 
        % 
        % Examples:
        % 1) name only targets: "betelgeuse", "ecliptic". 
        % 2) name and coordinates: "WD0012, 10:20:30.4-55:12:34.5".
        % 3) additional parameters: "galactic, , duration=2.5".
        % NOTE: in case (3) we left the coordinates empty, but you can give
        %       coordinates in addition to optional parameters. 
        %
        % NOTE: there is no need to use quotes of any kind. 
        
            import util.text.cs;
            import util.text.parse_bool; 
            
            if nargin==1, help('obs.sched.Target.parseText'); return; end
        
            if nargin<3 || isempty(use_resolve)
                use_resolve = 0;
            end
            
            if isempty(str)
                error('dome_pc:scheduler:target:parse:empty_string', 'Must give a non-empty string');
            end
            
            % remove all quotes
            str = strrep(str, '"', ''); 
            str = strrep(str, '''', ''); 
            
            c = strsplit(str, ','); % split the input by commas
            
            c = cellfun(@strtrim, c, 'UniformOutput', false); % remove white space
            
            if length(c)>1 && ~isempty(c{2}) && isletter(c{2}(1)) % this means we have skipped giving the RA/Dec
                c = {c{1} '' c{2:end}};
            end
            
            % read the target name
            if isempty(c{1})
                error('dome_pc:scheduler:target:parse:no_name', 'Must give an object name in the first field of the string');
            else
                obj.name = c{1};
            end
            
            % read the coordinates
            if length(c)>1 && ~isempty(c{2})
                
                idx = strfind(c{2}, '+'); % try to find a + separator
                if isempty(idx), idx = strfind(c{2}, '-'); end % try instead to find a - separator
                if isempty(idx)
                    error('dome_pc:shceduler:target:parse:no_coordinates', 'Must provide a coordinates field with a + or - separating the RA and Dec!\n Instead got: %s', c{2}); 
                end
                
                new_RA = strtrim(c{2}(1:idx-1));
                new_Dec = strtrim(c{2}(idx:end)); 
                
                obj.use_resolver = 0; 
                
            else
                obj.use_resolver = 1; 
                new_RA = [];
                new_Dec = [];
            end
            
            % parse additional parameters
            for ii = 3:length(c)
                
                if isempty(c{ii}), continue; end
                
                idx = strfind(c{ii}, '='); 
                if isempty(idx)
                    error('dome_pc:scheduler:target:parse:key_val_bad_format', 'Must provide optional arguments in "keyword=value" format.\n Instead got: %s', c{ii}); 
                end
                
                key = strtrim(c{ii}(1:idx-1));
                val = strtrim(c{ii}(idx+1:end)); 
                val = util.text.parse_value(val); 
                
                if cs(key, 'priority')
                    obj.priority = val;
                elseif cs(key, 'decay rate')
                    obj.decay_rate = val;
                elseif cs(key, 'tracking')
                    obj.tracking = parse_bool(val);
                elseif cs(key, 'cam_mode', 'camera mode', 'mode')
                    obj.cam_mode = val;
                elseif cs(key, 'exp time', 'exposure time')
                    obj.exp_time = val;
                elseif cs(key, 'keyword')
                    obj.ephem.keyword = val;
                elseif cs(key, 'list_filename')
                    obj.list_filename = val;
                elseif cs(key, 'sort_columns')
                    obj.sort_columns = val;
                elseif cs(key, 'repeat_maximum')
                    obj.repeat_max = val; 
                else % assume all other inputs are constraints on the target field
                    obj.ephem.constraints.scan_vars(key, val);
                end
                
            end
            
            if ~isempty(new_RA) && ~isempty(new_Dec)

%                 obj.ephem.input(new_RA, new_Dec); % just give the coordinates
                obj.RA = new_RA;
                obj.Dec = new_Dec;
                obj.ephem.updateSecondaryCoords;

            elseif use_resolve

                obj.ephem.resolve; % use the name resolver. Internally Ephemeris looks first for "keyword", then for "name" to resolve

                if isnan(obj.ephem.RA_deg)
                    warning('Could not parse the name "%s"!', obj.name); 
                end

            end
                
            obj.clear; % clear the observed time and so on
            
        end
        
        function copy_pars(obj, other) % get the parameters from another Target object, without reseting observation history
            
            if nargin<2 || isempty(other) || ~isa(other, 'obs.sched.Target')
                error('dome_pc:scheduler:target:copy_pars:bad_input', 'Must supply a valid Target object!'); 
            end
            
            % copy only the constraints
            obj.ephem.constraints = util.oop.full_copy(other.ephem.constraints); 
            
            obj.use_resolver = other.use_resolver; 
            
            obj.priority = other.priority; 
            obj.decay_rate = other.decay_rate;
            obj.tracking = other.tracking;
            obj.cam_mode = other.cam_mode;
            obj.exp_time = other.exp_time;
            
            
        end
        
        function val = observable(obj, varargin) % check if this object is inside the observational constraints
             
            import util.text.cs;
            
            val = obj.ephem.observable(varargin{:}); 
            
        end
        
        function best = best_target(obj_vec, varargin) % choose the best target out of the list
        % Usage: best = best_target(obj_vec, varargin)
        % Choose the best target from the list based on some criteria. 
        %
        % This is done by comparing the "current_priority" of objects. 
        % If two (observable) targets have the same "current_priority", 
        % then we choose one by setting the time of each target's Ephemeris 
        % and comparing them until the best target is chosen. 
        % 
        % The optional arguments are passed to Epehemris object, and used
        % to compare them until the best target is found. 
        % For example use 'side', 'West' to only take Western targets. 
        % Other options passed to Ephem.better_than() are:
        %  alt_limit=25 deg, duration=0.5 hour, side='both', wind_speed=[] 
        % 
        % Additional OPTIONAL ARGUMENTS that can be given to this function:
        %   -time: choose a time at which to compare targets. Default 'now'. 
        %          This can be used to simulate the decisions made in the 
        %          future or to reconstruct last night's priorities. 
        %   -alt_limit: minimal altitude (angle) above horizon for a target
        %               to be viable/observable. Default is 25 degrees.
        %   -duration: minimal duration (in hours) that a target must be
        %              observable, until meridian or horizon, to be valid. 
        %              Default is 1 hour. 
        %   -airmass: the maximum airmass that is considered observable. 
        %             Default is Inf, meaning we prefer lowest airmass but 
        %             will agree to take anything within the other limits. 
        %   -side: can only choose viable targets on the given side, that is
        %          "East", "West", or "both" (=default, same as empty). 
        %          Note that this overwrides the "wind_speed" input. 
        %   -wind_speed: if wind is above some limit (20km/h) then this 
        %                function will prefer the Eastern target regardless
        %                of airmass. If both are on the same side, this is 
        %                ignored. If not given (or empty) will skip this test.    
        % 
            
            import util.text.cs;
            
            if isempty(obj_vec)
                best = obs.sched.Target.empty; return;
            end
            
            input = util.text.InputVars;
            input.input_var('time', []); % time at which to compare these targets
            input.scan_vars(varargin{:}); 
            
%             if ~isempty(input.time) % get all targets up to date to a specific time (can use "now") make sure dynamic fields are updated! 
%                 
%                 for ii = 1:length(obj_vec)
%                     
%                     obj_vec(ii).ephem.time = head.Ephemeris.parseTime(input.time); 
%                     
%                     if obj_vec(ii).use_resolver
%                         obj_vec(ii).ephem.resolve; % also updates secondaries
%                     end
%                     
%                 end
%                 
%             end
            
            best = obs.sched.Target.empty; % no target is chosen yet
            
            for ii = 1:length(obj_vec)
                
                % first check that new target is inside given constraints
                if obj_vec(ii).observable(varargin{:})
                   
                    if isempty(best) % if we still didn't assign a best target
                        best = obj_vec(ii); 
                    elseif obj_vec(ii).getCurrentPriority>best.getCurrentPriority % compare to another target, check if the priority is higher
                        best = obj_vec(ii);
                    elseif obj_vec(ii).getCurrentPriority==best.getCurrentPriority % if the priority is the same, then check if the Ephemeris is better
                        if obj_vec(ii).ephem.better_than(best.ephem, varargin{:}) % better than best is the new best 
                            best = obj_vec(ii); 
                        end
                    end
                    
                end
                
            end % for ii
            
        end
        
        function start_observation(obj, time)
            
            if nargin<2 || isempty(time)
                time = util.text.time2str(datetime('now', 'TimeZone', 'UTC'));
            end
            
            if isempty(obj.start_time) % if this object is not currently observing!
                
                obj.start_time = time;
                obj.start_side = obj.side;
                obj.start_RA_deg = obj.ephem.RA_deg;
                obj.start_Dec_deg = obj.ephem.Dec_deg;
                
                obj.ephem.time = time;  
                obj.ephem.start_observing;
                
                s = struct('RA_deg', obj.start_RA_deg, 'Dec_deg', obj.start_Dec_deg, ...
                    'side', obj.start_side, 'start_time', obj.start_time, 'end_time', '', 'runtime', 0); 

                if isempty(obj.obs_history)
                    obj.obs_history = s;
                else
                    obj.obs_history(end+1) = s;
                end

            end
            
            % if the object is already observing we will simply ignore this command! 
            
        end
        
        function finish_observation(obj, time) % add a struct to obs_history, tracking the time this object was observed
            
            if nargin<2 || isempty(time)
                time = util.text.time2str(datetime('now', 'TimeZone', 'UTC'));
            end
            
            if isempty(obj.start_time)
                return;
%                 error('Can not finish an observation that has not been started yet! Use start_observation()'); 
            end
            
            if isempty(obj.obs_history(end).end_time)
                obj.obs_history(end).end_time = time; 
                obj.obs_history(end).runtime = obj.compare_times_hours(obj.start_time, time);                 
                obj.ephem.time = time;
                obj.ephem.finish_observing;
            else
                obj.ephem.now_observing = 0; % this should already be 0, but need to verify it
            end
            
            obj.start_time = ''; % when finished we no longer need to keep this info (it is stored in obs_history)
            
        end
        
        function loadFromTargetList(obj)
            
            if isempty(obj.list_filename)
                error('Cannot load a target list without a filename!'); 
            end
            
            if isempty(obj.list_table)
            % first get the most up-to-date target list into a table            
            d = util.sys.WorkingDirectory(fullfile(getenv('DATA'), 'WFAST/target_lists/')); 
            files = d.match([obj.list_filename '*']); % all files matching the filename + additional info like creation date
            files = sort(files); % get the latest file last
            
            if isempty(files)
                error('Could not find any files matching "%s". ', obj.list_filename); 
            end
            
            obj.list_table = readtable(files{end}); 
            obj.list_filename_loaded = files{end}; % keep track of the full name of the file we just loaded
            
            end
            
            T = obj.list_table; % shorthand
            
            time = datetime(T.LastObserved, 'TimeZone', 'UTC', 'Format', 'uuuu-MM-dd HH:mm:ss.sss'); % convert the text into datetime objects

            T.LastObserved = time; 

            %%%%%% choose the best target from this list %%%%%%
            
            % find the column names for sorting
            column_names = T.Properties.VariableNames; 
            cols = [];
            order = {}; 
            for jj = 1:length(obj.sort_columns)
            
                col_name = obj.sort_columns{jj};
                
                if col_name(end)=='*'
                    col_name = col_name(1:end-1); 
                    order{jj} = 'descend'; 
                else
                    order{jj} = 'ascend'; 
                end
                
                col_index = find(strcmp(column_names, col_name), 1, 'first'); 
                
                if isempty(col_index)
                    error('Unknown column name to sort on: "%s". Use one of the column names on the field table.', col_name); 
                end
                
                cols(jj) = col_index;
                
            end
            
            % create ephemeris objects for all fields
            e = util.oop.full_copy(obj.ephem); 
            for jj = height(T):-1:1 % reverse order make sure the allocation happens at the beginning
                fields(jj) = head.Ephemeris; % copy the constraints but keep it separate from the main ephemeris object
                fields(jj).time = e.time;
                fields(jj).constraints = e.constraints; 
                fields(jj).RA_deg = T{jj,'RightAscension'}; 
                fields(jj).Dec_deg = T{jj,'Declination'}; 
                fields(jj).field_id = T{jj,'Index'}; 
                if ~isempty(obj.ephem.field_id) && obj.ephem.field_id==fields(jj).field_id && obj.ephem.now_observing 
                    fields(jj).now_observing = 1; 
                else
                    fields(jj).now_observing = 0;
                end
            end
            
            % clear out data from last night, if this is a new night
            tonight = e.dateNight(e.time); % the datetime for the beginning of this night
            last_observed_nights = e.dateNight(T.LastObserved); 
            
            if all(last_observed_nights<tonight) % not yet observed any of these fields tonight
                
                for jj = 1:height(T) % clear consecutive nights for fields that were not observed last night
                    
                    T{jj, 'TotalHours'} = T{jj, 'TotalHours'} + T{jj, 'TonightHours'}; % accumulate last night's hours into the total count
                    T{jj, 'TonightHours'} = 0; % clear out any hours taken last night
                    
                    if last_observed_nights(jj) < tonight - days(1)
                        T{jj,'MaxConsecutive'} = max(T{jj,'MaxConsecutive'}, T{jj,'ConsecutiveNights'}); % log the longest stretch of consecutive nights before clearing them
                        T{jj,'ConsecutiveNights'} = 0; 
                    end
                end
                
            end
            
            % choose a new field! 
            idx = []; 

            select_idx = find(abs(e.angleDifference(e.LST_deg, T.RightAscension)/15)<2.5 & T.ConsecutiveNights<obj.repeat_max); % narrow down the fields a little
            if isempty(select_idx)
                obj.list_table = T; 
                return;
            end

            Tsel = T(select_idx,:); 

            [Tsort, sort_idx] = sortrows(Tsel, cols, order); 

            options_idx = []; 

            for nn = 0:1
                options_idx = vertcat(options_idx, find(e.dateNight(Tsort.LastObserved)==tonight-days(nn))); % find options for fields that were observed tonight or last night
            end

            options_idx = vertcat(options_idx, (1:height(Tsort))'); % there will be repeated indices, but not too many... 

            target = head.Ephemeris.empty; % fill this if we find a good field

            for jj = 1:length(options_idx) % go over our options, sorted according to: observed tonight, observed last night, everything else 

                idx = select_idx(sort_idx(options_idx(jj))); % convert to original indexing in T instead of the sorted and down-selected Tsort

                fields(idx).updateMoon; % maybe update secondaries? 

                if fields(idx).observable
                    target = fields(idx); 
                    break;
                end

            end

            if isempty(target)
                obj.ephem.RA = ''; 
                obj.ephem.Dec = ''; 
                obj.ephem.field_id = []; 
            else
                % copy the best field parameters
                obj.ephem.RA = target.RA; 
                obj.ephem.Dec = target.Dec; 
                obj.ephem.field_id = T{idx, 'Index'}; 
            end
            
            obj.list_table = T; 
            
        end
        
        function saveListTable(obj)
            
            if isempty(obj.list_filename) || isempty(obj.list_filename_loaded)
                error('Cannot load a target list without a filename!'); 
            end
            
            writetable(obj.list_table, obj.list_filename_loaded); % write the data back into the same file
            
        end
        
    end
    
    methods (Static=true)
        
        function obj_vec = readFile(filename, constraints) % read the text file line by line and parse() each line to a new Target object, returning a vector
        % Usage: obj_vec = readFile(filename='target_list.txt')
        
            if nargin<1 || isempty(filename)
                filename = 'target_list.txt';
            end

            if ~exist(filename, 'file')
                error('dome_pc:scheduler:target:read_file:bad_filename', 'Could not find the file "%s". ', filename);
            end

            if nargin<2 || isempty(constraints)
                constraints = [];
            else
                if ~isa(constraints, 'util.text.InputVars')
                    error('dome_pc:scheduler:target:read_file:bad_constraints', 'Second argument (constraints) to readFile must be a util.text.InputVars setup using a head.Ephemeris object! Instead got a %s', class(constraints)); 
                end
            end
            
            fid = fopen(filename);
            on_cleanup = onCleanup(@() fclose(fid)); % make sure this is called no matter how the function exits
            
            obj_vec = obs.sched.Target.empty; 
            
            for ii = 1:1e4 % arbitrary 
                
                tline = fgetl(fid);
                
                if ~ischar(tline)
                    return;
                end
                
                tline = strtrim(tline);
                
                if isempty(tline)
                    continue;
                end
                
                new_obj = obs.sched.Target; 
                
                if ~isempty(constraints)
                    new_obj.ephem.constraints = util.oop.full_copy(constraints);
                end
                
                new_obj.parse(tline); 
                
                obj_vec = [obj_vec; new_obj]; 
                obj_vec(end).index = length(obj_vec); 
                
            end
            
        end
        
        function val = compare_times_hours(time1, time2)
            
            if ischar(time1)
                time1 = util.text.str2time(time1);
            end
            
            if ischar(time2)
                time2 = util.text.str2time(time2);
            end
            
            d = time2-time1;
            
            val = hours(d); 
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

