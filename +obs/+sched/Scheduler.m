classdef Scheduler < handle
% Wrapper for a vector of obs.sched.Target objects, that also has a GUI
% and can save some constraints regarding the current observations. 
% The list of Targets can be used to choose the best target at any given 
% time, so it is very simple to use it as a scheduler. 
% Currently the fastest way to fill the list is using readFile(filename). 
%
% The important functions are update() which pushes the object to the current
% time, the choose() function that checks if a new target should be observed, 
% which also internally  sets the "continue_run" flag to 0 if the observation
% is over, telling the observer that a new slew/focus/run start needs to 
% be called. 
% Finally, after slewing to a new target, use start_current() to signal to 
% the scheduler that a new run has begun. 
%
% NOTE: each of these functions can be given a time (string or datetime obj)
%       that is parsed and used instead of current time (the default). 
%       You should also provide a second argument set to true to signal the 
%       scheduler to apply all calculations to the simulation list of targets. 
%       In real life, just don't give these arguments. 
%
% Additional functionality is the ability to plot the target positions and 
% observation times on a map of the sky (using the graphics from SkyMap). 
% The simulation mode allows you to move time forward and find what fields
% would be observed at each time, showing how the night would look in advance. 
% Other considerations can be controlled like wind speed, that will put 
% additional constraints on the target selection. So the simulation can run
% with different weather conditions, producing different outcomes. 
%
% Every time you use choose() the scheduler provides a "report" string which
% tells you what the next step should be. 
% There is also a "rationale" string showing all the targets and why that 
% one was selected. 
% In "rationale", the parenthesis includes a reason why a target was excluded, 
% or if it wasn't excluded, it shows the priority and airmass. 
% A target currently being observed has a star next to the airmass, to show
% that it has a slight preference to other targets with similar airmass. 
% The chosen target is highlighted with ***** target ******. 
%
% If you need to get a play-by-play of the reports or rationales, look at 
% the "report_log" or "rationale_log", and their simulation counterparts the
% "report_log_sim" and "rationale_log_sim". 
    
    properties(Transient=true)
        
        gui;
        
    end
    
    properties % objects
        
        map@util.ast.SkyMap; % load this from file, use it to plot the results
        
        targets@obs.sched.Target; % a vector of targets we get from the "target_list.txt" file
        current@obs.sched.Target; % a single Target object that is currently being observed
        
        targets_sim@obs.sched.Target; 
        current_sim@obs.sched.Target;
        
        ephem@head.Ephemeris; % keep track of time and superposed constraints
        
        log@util.sys.Logger; 
        
    end
    
    properties % inputs/outputs
        
        wind_speed = []; % current wind speed given from sensors
%         current_RA = []; % current right ascention, given from last start of new target (numeric degrees)
%         current_Dec = []; % current declination, given from last start of new target (numeric degrees)
        current_side = ''; % current observing side, given from last start of new target
%         need_to_start = 1; % when true, must give telescope the new target, slew and refocus before calling start_current()
        
        % keep track of observation time
        total_night_time = 0; % how much total dark time we had tonight
        total_observed_time = 0; % how much was actually observed (should be equal to night time if there are no conflicting constraints or bad weather)
        obs_history = []; % a strcut array for runs we have done, with fields: name, index (from targets), RA_deg, Dec_deg, start_time and end_time
        
        continue_run = 0; % if this is true, the call to "choose" has concluded with us needing to continue the same target
        report = ''; % report what is the next move for the telescope
        rationale = ''; % save a short description why this target is chosen
        
        report_log = {}; % save a copy of all reports
        rationale_log = {}; % save a copy of all rationales
        
    end
    
    properties % simulation parameters and outputs
        
        sim_time_step = 30; % how many minutes go by before testing for a new target
        sim_wind_start_time = []; % time at which the wind picks up at night, for simulation mode only!         
        sim_wind_end_time = []; % time at which the wind calms down , for simulation mode only! 
        sim_starting_side = 'East'; % it is much easier to just simulate this, no wind considerations. Only for simulations and only when use_stay_on_side=1
        
        sim_slew_time = 0; % simulate the time (in minutes) it takes to slew to target
        sim_flip_time = 0; % simulate the time (in minutes) it takes to do meridian flip
        sim_focus_time = 0; % time to do focus between runs (in minutes). assume focus is done for each new run
        
        sim_plot_pause = 0.0; % length of pause for each simulation sampling point when plotting the simulation
        
        sim_time = [];
        current_side_sim = '';
        obs_history_sim = []; 
        report_log_sim = {}; % save a copy of all reports
        rationale_log_sim = {}; % save a copy of all rationales
        total_night_time_sim = 0;
        total_observed_time_sim = 0;
        
        showing_sim = 0; 
        
    end
    
    properties % switches/controls
        
        filename = fullfile(getenv('DATA'), 'WFAST/target_lists/target_list.txt');
%         filename = 'target_list.txt'; % use the GUI to select a better default file
        
        max_sun_elevation = -10; % consider it night time when sun is this far below horizon
        max_wind = 30; % above this wind speed, we must only observe the East side
        use_stay_on_side = 0; % if this is true, must not choose any target on the flip side (user may choose to start on East side to avoid getting stuck on West when it is windy)
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
    end
    
    properties(Hidden=true)
       
        brake_bit = 1; % use this to stop a simulation in progress... 
        
        prev_datetime = []; % keep track of last time update() was called, so we can sum up the total observing time
        prev_target@obs.sched.Target; % keep track of the last thing we observed (do we even need this? maybe just for reference/debugging)
        
        prev_datetime_sim = []; 
        prev_target_sim@obs.sched.Target; 

        default_filename;
        
        version = 1.02;
        
    end
    
    methods % constructor
        
        function obj = Scheduler(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.sched.Scheduler')
                if obj.debug_bit>1, fprintf('Scheduler copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                
                if obj.debug_bit>1, fprintf('Scheduler constructor v%4.2f\n', obj.version); end
                
                util.oop.save_defaults(obj); % save the value in "filename" to "default_filename" as well (and maybe other properties)
                
                load(fullfile(getenv('DATA'), '/WFAST/saved/sky_map')); 
                
                obj.map = sky_map;
                
                obj.ephem = obj.map.ephem; 
                
                obj.log = util.sys.Logger('Scheduler'); 
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj) % this is called before reading a new target list
            
            obj.targets = obs.sched.Target.empty; 
            obj.current = obs.sched.Target.empty; 
            
            obj.targets_sim = obs.sched.Target.empty; 
            obj.current_sim = obs.sched.Target.empty; 
            
            obj.clear;
            
        end
        
        function clear(obj, use_sim) % call this at the start of the night/simulation
            
            if nargin<2 || isempty(use_sim)
                use_sim = 0;
            end
            
            
            if use_sim==0

                for ii = 1:length(obj.targets)
                    obj.targets(ii).clear;
                end
                
                obj.current = obs.sched.Target.empty; 
                obj.total_night_time = 0;
                obj.total_observed_time = 0;

                obj.prev_datetime = [];
                obj.prev_target = obs.sched.Target.empty;

                obj.obs_history = [];

                obj.report_log = {};
                obj.rationale_log = {}; 

%                 obj.current_RA = [];
%                 obj.current_Dec = [];
                obj.current_side = '';
                
            else

                for ii = 1:length(obj.targets_sim)
                    obj.targets_sim(ii).clear;
                end
                
                obj.current_sim = obs.sched.Target.empty; 
                obj.total_night_time_sim = 0;
                obj.total_observed_time_sim = 0;

                obj.prev_datetime_sim = [];
                obj.prev_target_sim = obs.sched.Target.empty;

                obj.obs_history_sim = [];

                obj.report_log_sim = {};
                obj.rationale_log_sim = {}; 
                
                obj.current_side_sim = '';
                
            end
            
        end
        
    end
    
    methods % getters
        
        function val = find(obj, name, use_sim)
            
            if nargin<3 || isempty(use_sim)
                use_sim = 0;
            end
            
            if ~use_sim
                idx = find(contains({obj.targets.name}, name));
            else
                idx = find(contains({obj.targets_sim.name}, name));
            end
            
            if isempty(idx)
                val = obs.sched.Target.empty;
            elseif ~use_sim
                val = obj.targets(idx); 
            else
                val = obj.targets_sim(idx);
            end
            
        end
        
        function val = wind_state(obj)
            
            if obj.wind_speed>obj.max_wind
                val = 1;
            else
                val = 0; % this includes when wind is empty! 
            end
            
        end
        
        function val = wind_string(obj)
            
            if obj.wind_state
                val = 'windy!';
            else
                val = 'calm';
            end
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function readFile(obj, filename)
            
            if nargin<2 || isempty(filename)
                filename = obj.filename;
            end
            
            obj.reset;
            
            obj.targets = obs.sched.Target.readFile(filename, obj.ephem.constraints);
            
            obj.targets_sim = util.oop.full_copy(obj.targets); 
            
        end
        
        function update(obj, time, use_sim)
            
            if nargin<2 || isempty(time)
                time = datetime('now', 'TimeZone', 'UTC');
            elseif ischar(time)
                time = util.text.str2time(time); 
            end
            
            if nargin<3 || isempty(use_sim)
                use_sim = 0;
            end
            
            if ~use_sim

                if ~isempty(obj.prev_datetime)
                    interval_sec = double(time-obj.prev_datetime); % how long has it been since last update
                    obj.total_night_time = obj.total_night_time + interval_sec/60; % convert to minutes! 
                end

                obj.prev_datetime = time; 

            else

                if ~isempty(obj.prev_datetime_sim)
                    interval_sec = double(time-obj.prev_datetime_sim); % how long has it been since last update
                    obj.total_night_time_sim = obj.total_night_time_sim + interval_sec/60; % convert to minutes! 
                end

                obj.prev_datetime_sim = time; 
            
            end
            
            obj.ephem.time = time;
            obj.ephem.updateSecondaryCoords;
            
        end
        
        function record_rational(obj, time, target_list, new_target, dur_flag)
            
            if nargin<5 || isempty(dur_flag)
                dur_flag = true(length(target_list),1); 
            end
            
            if isa(time, 'datetime')
                time = util.text.time2str(time);
            end
            
            obj.rationale = [time ': '];

            for ii = 1:length(target_list)

                new_str = target_list(ii).name;

                if ~dur_flag(ii) 
                    new_str = sprintf('%s (Duration exceeded)', new_str);
                elseif ~isempty(target_list(ii).ephem.unobservable_reason) % ephem object already recorded why we can't observe this
                    new_str = sprintf('%s (%s)', new_str, target_list(ii).ephem.unobservable_reason);
                else % target is observable, should choose it based on priority and airmass
                    
                    if target_list(ii).ephem.now_observing % currently observed target has some advantage even if it has slightly lower airmass
                        star = '*';
                    else
                        star = '';
                    end
                    
                    new_str = sprintf('%s (P= %4.2f, AM= %4.2f%s)', new_str, target_list(ii).getCurrentPriority, target_list(ii).ephem.AIRMASS, star);
                
                end

                if isequal(target_list(ii), new_target) % highlight the chosen target
                    new_str = sprintf('**** %s ****', new_str); 
                end

                if ii==1
                    obj.rationale = [obj.rationale new_str]; 
                else
                    obj.rationale = [obj.rationale ', ' new_str]; 
                end
                
            end
            
        end
        
        function matchRuntimes(obj, obs_log) % use obs_log to update the target list with the true runtimes and start times
            
            for ii = 1:length(obj.targets)
                
                name = obj.targets(ii).name;
                
                if isfield(obs_log, name)
                    
                    st = obs_log.(name);
                    dur = 0; % total duration of previous runs (not including current run)
                    
                    for jj = 1:length(st)

                        if ~isempty(st(jj).runtime)
                            dur = dur + st(jj).runtime;
                        end
                            
                        if jj==length(st) % last iteration, update the target ephemeris
                            
                            obj.targets(ii).ephem.prev_runtime_minutes = dur/60; % convert to minutes
                            obj.targets(ii).ephem.update;
                            
                            if ~isempty(st(jj).start)
                                obj.targets(ii).ephem.STARTTIME = st(jj).start; 
                            end

                            if ~isempty(st(jj).end) % the last object in the log has finished observing!
                                obj.targets(ii).ephem.now_observing = 0;
                            else % the last object in the log is still running! 
                                obj.targets(ii).ephem.now_observing = 1;
                            end

                        end
                        
                    end
                    
                else
                    % if it doesn't appear on the obs log, should we reset the runtime? 
                end
                
            end
            
            
        end
        
        function new_target = choose(obj, time, varargin)
            
            if nargin<2 || isempty(time)
                time = util.text.time2str(datetime('now', 'TimeZone', 'UTC'));
            elseif isa(time, 'datetime')
                time = util.text.time2str(time); 
            elseif ischar(time)
                if util.text.cs(time, 'now')
                    time = util.text.time2str(datetime('now', 'TimeZone', 'UTC'));
                end
            end
            
            idx_sim = find(strcmpi('use_sim', varargin)); 
            
            if isempty(idx_sim)
                use_sim = 0;
            else
                use_sim = 1;
                varargin(idx_sim) = []; 
            end
            
            arguments = obj.ephem.constraints.output_vars;
            
            % add to arguments the chosen side if "stay_on_side" is true, or because of wind
            % ...
            
            if use_sim
                this_side = obj.current_side_sim;
            else
                this_side = obj.current_side;
            end
            
            if obj.wind_state
                arguments = [arguments 'side', 'East'];
            elseif obj.use_stay_on_side 
                arguments = [arguments 'side', this_side];
            end
            
            if ~isempty(varargin)
                arguments = [arguments varargin];
            end
            
%             fprintf('%s: side= %s\n', time, arguments{end}); 
            
            if ~use_sim
                target_list = obj.targets; % make a copy of the list, that can get shorter because of various external constraints
            else
                target_list = obj.targets_sim; % make a copy of the list, that can get shorter because of various external constraints
            end
            
            new_target = obs.sched.Target.empty;
            
            % see if any of the targets is now observing and is still inside the "continuous" constraint
            for ii = 1:length(target_list)
                
                e = target_list(ii).ephem; % shorthand
                e.time = time; % make sure all targets are updated to current time
                e.updateSecondaryCoords;
                
                if e.now_observing && strcmpi(e.side, this_side) % make sure this object is also on the same side as telescope: if we need to flip, we need to choose a new target
                    
                    if e.getRuntimeMinutes + e.constraints.fudge_time < e.constraints.continuous*60
                        new_target = target_list(ii); % this target must be observed for some time before a new target can be observed... 
                        obj.rationale = sprintf('%s: Target %s has been observed for only %d minutes! Continuing observations... ', time, new_target.name, e.getRuntimeMinutes); 
                        obj.continue_run = 1; % by default we do not need to continue the run, and will instead start a new one
                        break; % skip the other targets (assume there is only one target being observed at each time)
                    end
                    
                end
                
            end
            
            % if we didn't find any target that is still constrained by "continuous", look for the best available target
            if isempty(new_target) % this happens if there is no target currently under "continuous" condition
                
                for ii = 1:length(target_list) % go over the list an re-resolve any targets that need to be updated
                    if target_list(ii).use_resolver
                        target_list(ii).ephem.resolve([], arguments{:}); % also updates secondaries
                    end 
                end

                dur_flag = true(length(target_list),1); % indices of the targets that meet the duration prerequisite

                for ii = 1:length(target_list)
                    
                    e = target_list(ii).ephem; % shorthand
                    
                    if e.getTotalRuntimeMinutes + e.constraints.fudge_time > e.constraints.duration*60
                        dur_flag(ii) = false; % flag as false targets that have been observed more than "duration"
                    end
                    
                end

                target_list_duration = target_list(dur_flag); % narrow down to a subset of legal targets...
                
                new_target = best_target(target_list_duration, 'time', time, arguments{:}); % try to find a good target under the "duration" constraint
                
                if ~isempty(new_target) % we managed to find a good target below "duration" condition
                    obj.record_rational(time, target_list, new_target, dur_flag); % keep a record for each target why it was or wasn't picked
                else % this happens if we didn't find any targets that haven't been observed "duration" hours -> choose a target that has surpassed the "duration" condition
                    
                    new_target = best_target(target_list(~dur_flag), 'time', time, arguments{:}); % run only the targets that have longer than duration
                    
                    if ~isempty(new_target)
                        obj.record_rational(time, target_list, new_target); % keep a record for each target why it was or wasn't picked
                    end
                    
                end
                
            end
            
            obj.continue_run = 0; % in most cases we do not need to continue the run. This can change if the best target chosen is the same one we are observing
            
            % load the target we are now observing
            if ~use_sim
                this_target = obj.current;
                this_side = obj.current_side;
            else
                this_target = obj.current_sim;
                this_side = obj.current_side_sim;
            end
            
            %%%%%% now we have made the selection of the best target, lets see what we can do with it %%%%%%% 
            if isempty(new_target) % no targets are currently observable! 
                
                obj.record_rational(time, target_list, new_target); % keep a record for each target why it was or wasn't picked
                
                obj.report = sprintf('%s: No available targets. Going to idle mode...', time); 
                
            elseif ~isequal(this_target, new_target) || ~obj.compare_coordinates(new_target, use_sim) % different target
            % new target is different from the current one, or current target is empty, 
            % or targets are the same but have different coordinates 
            % (dynamic fields can have different coords for the same name)
                
                obj.report = sprintf('%s: Moving to new object: %50s', time, new_target.summary); 
                
            elseif ~isempty(this_side) && ~isequal(this_side, new_target.side) % same target, new side
                
                obj.report = sprintf('%s: Flip to same object:  %50s', time, new_target.summary); 
                
            else % same target, same side -> continue observing
                obj.report = sprintf('%s: Continue observing:   %50s', time, new_target.summary); 
                obj.continue_run = 1; % this is the only case where we do not need to start a new run! 
            end
            
%             if ~use_sim
%                 obj.report_log{end+1,1} = obj.report;
%             else
%                 obj.report_log_sim{end+1,1} = obj.report;
%             end
            
            if obj.debug_bit>1
                
                if isempty(new_target)
                    disp(obj.report); 
                else
                    disp([obj.report ' | ' new_target.details]); 
                end
                
            end
            
        end
        
        function start_current(obj, time, use_sim)
            
            if nargin<2 || isempty(time)
                time = util.text.time2str(datetime('now', 'TimeZone', 'UTC'));
            elseif isa(time, 'datetime')
                time = util.text.time2str(time); 
            end
            
            if nargin<3 || isempty(use_sim)
                use_sim = 0;
            end
            
            if ~use_sim
            
                if ~isempty(obj.current)

                    obj.current.start_observation(time); 
                    obj.current_side = obj.current.side; 

                    % add the obs_history from the Target object to the Scheduler list
                    s = obj.current.obs_history(end); 
                    s.index = obj.current.index;
                    s.name = obj.current.name;
                    if isempty(obj.obs_history)
                       obj.obs_history = s;
                    else
                        obj.obs_history(end+1) = s;
                    end

                end
            
            else
                
                if ~isempty(obj.current_sim)

                    obj.current_sim.start_observation(time); 
                    obj.current_side_sim = obj.current_sim.side; 

                    % add the obs_history from the Target object to the Scheduler list
                    s = obj.current_sim.obs_history(end); 
                    s.index = obj.current_sim.index;
                    s.name = obj.current_sim.name;
                    if isempty(obj.obs_history_sim)
                       obj.obs_history_sim = s;
                    else
                        obj.obs_history_sim(end+1) = s;
                    end

                end
                
            end
            
            % silently ignore this command if the current target is empty (idle mode)
            
        end
        
        function finish_current(obj, time, use_sim)
            
            if nargin<2 || isempty(time)
                time = util.text.time2str(datetime('now', 'TimeZone', 'UTC'));
            elseif isa(time, 'datetime')
                time = util.text.time2str(time); 
            end
            
            if nargin<3 || isempty(use_sim)
                use_sim = 0;
            end
            
            if ~use_sim

                if isempty(obj.current)
                    return; % silently ignore this command if the current target is empty (idle mode)
                end

                obj.current.finish_observation(time);

                if ~isempty(obj.current) && ~isempty(obj.current.obs_history)
                    obj.obs_history(end).end_time = obj.current.obs_history(end).end_time;
                    obj.obs_history(end).runtime = obj.current.obs_history(end).runtime;
                end

%                 obj.current = obs.sched.Target.empty;
                
            else
                
                if isempty(obj.current_sim)
                    return; % silently ignore this command if the current target is empty (idle mode)
                end

                obj.current_sim.finish_observation(time);
                
                if ~isempty(obj.obs_history_sim)
                    obj.obs_history_sim(end).end_time = obj.current_sim.obs_history(end).end_time;
                    obj.obs_history_sim(end).runtime = obj.current_sim.obs_history(end).runtime;
                end
                
            end
            
        end
        
        function val = compare_coordinates(obj, target, use_sim) % true if current coordinates are close enough to the new target
            
            if nargin<3 || isempty(use_sim)
                use_sim = 0;
            end
            
%             thresh = 30/60; % let's take 30 arcminutes as the threshold for moving to a new field?
            thresh = 1; % one degree threshold, on either RA or Dec, to consider this field as a separate field...
            
            if ~use_sim

                if isempty(obj.obs_history) || isempty(obj.obs_history(end).RA_deg) || isempty(obj.obs_history(end).Dec_deg)
                    val = 0; % no current coordinates, so it can't be close enough
                elseif abs(target.ephem.RA_deg-obj.obs_history(end).RA_deg)<thresh && ...
                        abs(target.ephem.Dec_deg-obj.obs_history(end).Dec_deg)<thresh % both coordinates are close enough
                    val = 1; 
                else
                    val = 0; 
                end

            else
                
                if isempty(obj.obs_history_sim(end).RA_deg) || isempty(obj.obs_history_sim(end).Dec_deg)
                    val = 0; % no current coordinates, so it can't be close enough
                elseif abs(target.ephem.RA_deg-obj.obs_history_sim(end).RA_deg)<thresh && ...
                        abs(target.ephem.Dec_deg-obj.obs_history_sim(end).Dec_deg)<thresh % both coordinates are close enough
                    val = 1; 
                else
                    val = 0; 
                end
                
            end
            
        end
        
        function write_log(obj)
            
            if ~isempty(obj.report)
                obj.log.input(obj.report);
            end
            
            if ~isempty(obj.rationale)
                obj.log.input(obj.rationale);
            end
            
        end
        
    end
    
    methods % simulations
        
        function run_simulation(obj, use_clear, start_time, end_time)
            
            if nargin<2 || isempty(use_clear)
                use_clear = 0; 
            end
            
            if nargin<3 || isempty(start_time)
                start_time = datetime('today', 'TimeZone', 'UTC'); 
                start_time.Hour = 13; % set the time to 16:00 Israel time (in winter??) 
            end
            
            if ischar(start_time)
                start_time = head.Ephemeris.parseTime(start_time);
            end
            
            if nargin<4 || isempty(end_time)
                end_time = datetime('today', 'TimeZone', 'UTC'); 
                end_time = end_time + days(1); % tomorrow morning! 
                end_time.Hour = 04; % set the time to 7:00 Israel time (in winter??) 
            end
            
            if ischar(end_time)
                end_time = head.Ephemeris.parseTime(end_time);
            end
            
            obj.brake_bit = 0; % start running
            
            if use_clear % only start a new simulation if requested specifically to do so... 
                
                obj.sim_time = start_time; % this is the virtual clock we will use throughout
                obj.clear(1); % start a new night (the argument is for use_sim=1)
                
                % move forward in time without calling update() until reaching night time
                for ii = 1:1e4 % arbitrary timeout
                
                    obj.ephem.time = obj.sim_time;

                    if obj.brake_bit
                        return;
                    end

                    obj.ephem.updateSun;

                    if obj.ephem.sun.Alt<obj.max_sun_elevation % sunset! 
                        break; % go to next loop with observations... 
                    end

                    obj.sim_time = obj.sim_time + minutes(obj.sim_time_step); 

                end

                obj.report = sprintf('Sunset at %s. Starting observations...', obj.sim_time); 
                obj.report_log_sim{end+1,1} = obj.report;
                if obj.debug_bit>1, disp(obj.report); end

                obj.current_side_sim = obj.sim_starting_side; % simulations will start with telescope on this side

            end
            
            % move forward in time while observing targets, until sun comes up
            for ii = 1:1e4 % arbitrary timeout
                
                obj.ephem.time = obj.sim_time;
                
                if obj.brake_bit
%                     obj.finish_current(sim_time);
                    break;
                end
                
                obj.ephem.updateSun;
                
                if obj.ephem.sun.Alt>obj.max_sun_elevation % sunrise! 
                    
                    obj.finish_current(obj.sim_time, 1); % second argument is for use_sim=1 
                    obj.brake_bit = 1;

                    obj.report = sprintf('Sunrise at %s. Finished observations...', obj.sim_time);
                    obj.report_log_sim{end+1,1} = obj.report;
                    if obj.debug_bit>1, disp(obj.report); end

                    % need to call anything else to wrap up?

                    break; % finish observing for tonight
                    
                end
                
                new_target = obj.choose(obj.sim_time, 'use_sim'); 
                
                if obj.continue_run==0 % what happens if we continue the same target after a flip? the prev_target will be same as current! is this the right thing to do??
                    obj.prev_target_sim = obj.current_sim;
                    obj.finish_current(obj.sim_time, 1); % second argument is for use_sim=1 
                    obj.current_sim = new_target;
                end
                
                obj.report_log_sim{end+1,1} = obj.report; 
                obj.rationale_log_sim{end+1,1} = obj.rationale;
            
                if isempty(obj.current_sim)
%                     if obj.debug_bit>1, fprintf('No targets available... remain in idle mode\n'); end
                    % do nothing but wait
                else % if obj.need_to_start 
                    
                    if ~isequal(obj.current_side_sim, obj.current_sim.side)
                        obj.sim_time = obj.sim_time + minutes(obj.sim_flip_time);
                    end
                    
                    obj.sim_time = obj.sim_time + minutes(obj.sim_slew_time);
                    obj.sim_time = obj.sim_time + minutes(obj.sim_focus_time);
                    
                    obj.start_current(obj.sim_time, 1); % second argument is for use_sim=1  

                    if ~isempty(obj.gui)
                        obj.gui.update(1); % the argument is for use_sim=1 
                        obj.show('use_sim', 1); 
                    end
                    
                end
                
                obj.sim_time = obj.sim_time + minutes(obj.sim_time_step); 
                
                pause(obj.sim_plot_pause); % pause to make the simulation go a little slower for plotting/visualization
                
            end
            
            if ~isempty(obj.gui)
                obj.gui.update(1); % the argument is for use_sim=1 
                obj.show('use_sim', 1);
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis'); 
            input.input_var('font_size', 20); 
            input.input_var('marker_size', 18);
            input.input_var('use_sim', false);
            input.scan_vars(varargin{:}); 
            
            if isempty(input.ax)
                
                if ~isempty(obj.gui) && obj.gui.check
                    input.ax = obj.gui.axes_image;
                else
                    input.ax = gca;
                end
                
            end
                        
            obj.map.show('ax', input.ax, 'vector', [-5 5], 'LST', obj.ephem.LST_deg/15,...
                'log', 1, 'galactic', 0, 'zenith', 1, 'horizon', 1, ...
                'grid', 0, 'units', 'hours', 'font size', input.font_size); 
            
            title(input.ax, ''); 
            input.ax.FontSize = input.font_size;
            
            input.ax.NextPlot = 'add';
            
            if ~input.use_sim
                history = obj.obs_history;
                this_target = obj.current;
            else
                history = obj.obs_history_sim;
                this_target = obj.current_sim;
            end
            
            % plot current position
            if ~isempty(this_target) && isvalid(this_target) && ~isempty(history(end).RA_deg) && ...
                    ~isempty(history(end).Dec_deg) && ~isempty(history(end).side)
                
                obj.plotObservation(history(end), input.ax, input.marker_size); 

                plot(input.ax, history(end).RA_deg/15, history(end).Dec_deg, 'yo', 'MarkerSize', input.marker_size); 
                
            end
            
            % plot previous targets
            for ii = 1:length(history)
                obj.plotObservation(history(ii), input.ax, input.marker_size); 
            end
            
            input.ax.NextPlot = 'replace';
            
            obj.showing_sim = input.use_sim;
            
        end
        
        function plotObservation(obj, s, ax, marker_size, font_size)
            
            if nargin<5 || isempty(font_size)
                font_size = 16;
            end
            
            if strcmp(s.side, 'West')
                plot(ax, s.RA_deg/15, s.Dec_deg, 'gx', 'MarkerSize', marker_size); 
            else
                plot(ax, s.RA_deg/15, s.Dec_deg, 'g+', 'MarkerSize', marker_size); 
            end
            
            text(ax, s.RA_deg/15+0.2, s.Dec_deg+5, strrep(s.name, '_', ' '), 'Color', 'g', 'FontSize', font_size); 
            
            % at some point we will have to figure out a way to put multiple printouts for the same field, beyond East/West...             
            shift_down = 0;
            if strcmp(s.side, 'West')
                shift_down = 7;
            end
            
            if ~isempty(s.start_time) && ~isempty(s.end_time)
                text(ax, s.RA_deg/15+0.3, s.Dec_deg-shift_down, [s.start_time(12:16) '-' s.end_time(12:16)], 'Color', 'g', 'FontSize', font_size); 
            elseif ~isempty(s.start_time)
                text(ax, s.RA_deg/15+0.3, s.Dec_deg-shift_down, s.start_time(12:16), 'Color', 'y', 'FontSize', font_size); 
            end
            
        end
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = obs.sched.gui.SchedGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
        function constraintsGUI(obj)
            
            obj.ephem.constraints.makeGUI;
            
        end
        
    end    
    
end

