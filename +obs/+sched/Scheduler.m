classdef Scheduler < handle
% Wrapper for a vector of obs.sched.Target objects, that also has a GUI
% and can save some constraints regarding the current observations. 
% The list of Targets can be used to choose the best target at any given 
% time, so it is very simple to use it as a scheduler. 
% Currently the fastest way to fill the list is using readFile(filename). 
%
% The important functions are update() which pushes the object to the current
% time, the choose() function that checks if a new target should be observed, 
% which also internally calls finish_current() and sets the "need_to_start"
% flag to 1 (telling the observer that a new slew/focus/run start needs to 
% be called). Finally, after slewing to a new target, use start_current()
% to signal to the scheduler that a new run has begun. 
% NOTE: each of these functions can be given a time (string or datetime obj)
%       that is parsed and used instead of current time (the default). This
%       is useful for simulations. In real life, just don't give this arg. 
%
% Additional functionality is the ability to plot the target positions and 
% observation times on a map of the sky (using the graphics from SkyMap). 
% The simulation mode allows you to move time forward and find what fields
% would be observed at each time, showing how the night would look in advance. 
% Other considerations can be controlled like wind speed, that will put 
% additional constraints on the target selection. So the simulation can run
% with different weather conditions, producing different outcomes. 
    
    properties(Transient=true)
        
        gui;
        
    end
    
    properties % objects
        
        map@util.ast.SkyMap; % load this from file, use it to plot the results
        
        targets@obs.sched.Target; % a vector of targets we get from the "target_list.txt" file
        current@obs.sched.Target; % a single Target object that is currently being observed
        
        ephem@head.Ephemeris; % keep track of time and superposed constraints
        
        
    end
    
    properties % inputs/outputs
        
        wind_speed = []; % current wind speed given from sensors
        current_RA = []; % current right ascention, given from last start of new target (numeric degrees)
        current_Dec = []; % current declination, given from last start of new target (numeric degrees)
        current_side = ''; % current observing side, given from last start of new target
        need_to_start = 1; % when true, must give telescope the new target, slew and refocus before calling start_current()
        
        % keep track of observation time
        total_night_time = 0; % how much total dark time we had tonight
        total_observed_time = 0; % how much was actually observed (should be equal to night time if there are no conflicting constraints or bad weather)
        obs_history = []; % a strcut array for runs we have done, with fields: name, index (from targets), RA_deg, Dec_deg, start_time and end_time
        
        report = ''; 
        
    end
    
    properties % simulation parameters
        
        sim_time_step = 30; % how many minutes go by before testing for a new target
        sim_wind_start_time = []; % time at which the wind picks up at night, for simulation mode only!         
        sim_wind_end_time = []; % time at which the wind calms down , for simulation mode only! 
        sim_starting_side = 'East'; % it is much easier to just simulate this, no wind considerations. Only for simulations and only when use_stay_on_side=1
        
        sim_slew_time = 0; % simulate the time (in minutes) it takes to slew to target
        sim_flip_time = 0; % simulate the time (in minutes) it takes to do meridian flip
        sim_focus_time = 0; % time to do focus between runs (in minutes). assume focus is done for each new run
        
        sim_plot_pause = 0.5; % length of pause for each simulation sampling point when plotting the simulation
        
    end
    
    properties % switches/controls
        
        filename = 'target_list.txt'; % use the GUI to select a better default file
        
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
        prev_side = ''; % I think we don't need this! 
        
        default_filename;
        
        version = 1.00;
        
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
                obj.map.show_galactic = 0;
                obj.map.show_ra_units = 'hours';
                obj.map.show_log = 1;
                
                obj.ephem = obj.map.ephem; 
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj) % this is called before reading a new target list
            
            obj.targets = obs.sched.Target.empty; 
            obj.current = obs.sched.Target.empty; 
            
            obj.clear;
            
        end
        
        function clear(obj) % call this at the start of the night/simulation
            
            for ii = 1:length(obj.targets)
                obj.targets(ii).clear;
            end
            
            obj.current = obs.sched.Target.empty; 

            obj.current_RA = [];
            obj.current_Dec = [];
            obj.current_side = '';
            
            obj.total_night_time = 0;
            obj.total_observed_time = 0;
            
            obj.prev_datetime = [];
            obj.prev_target = obs.sched.Target.empty;
            obj.prev_side = ''; 
        
            obj.obs_history = [];
            
        end
        
    end
    
    methods % getters
        
        function val = find(obj, name)
            
            idx = find(contains({obj.targets.name}, name));
            
            if isempty(idx)
                val = obs.sched.Target.empty;
            else
                val = obj.targets(idx); 
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
            
            obj.targets = obs.sched.Target.readFile(filename);
            
        end
        
        function update(obj, time)
            
            if nargin<2 || isempty(time)
                time = datetime('now', 'TimeZone', 'UTC');
            elseif ischar(time)
                time = util.text.str2time(time); 
            end
            
            if ~isemprt(obj.prev_datetime)
                interval_sec = double(time-obj.prev_datetime); % how long has it been since last update
                obj.total_night_time = obj.total_night_time + interval_sec/60; % convert to minutes! 
            end
            
            obj.prev_datetime = obj.ephem.time; 
            
            obj.ephem.time = time;
            
        end
        
        function choose(obj, time)
            
            if nargin<2 || isempty(time)
                time = util.text.time2str(datetime('now', 'TimeZone', 'UTC'));
            elseif isa(time, 'datetime')
                time = util.text.time2str(time); 
            end
            
            arguments = obj.ephem.constraints.output_vars;
            
            % add to arguments the chosen side if "stay_on_side" is true, or because of wind
            % ...
            
            new_target = best_target(obj.targets, 'time', time, arguments{:});
            
            if isempty(new_target) % no targets are currently observable! 
                
                obj.report = sprintf('%s: No available targets. Going to idle mode...', time); 
                if obj.debug_bit>1, disp(obj.report); end
                
                obj.prev_target = obj.current;
                obj.finish_current(time);
                obj.current = new_target;
                
            elseif ~isequal(obj.current, new_target) || ~obj.compare_coordinates(new_target) 
            % new target is different from the current one, or current target is empty, 
            % or targets are the same but have different coordinates 
            % (dynamic fields can have different coords for the same name)
                
                obj.report = sprintf('%s: Moving to new object: %50s', time, new_target.summary); 
                if obj.debug_bit>1, disp([obj.report ' | ' new_target.details]); end
                
                obj.prev_target = obj.current;
                obj.finish_current(time);
                obj.current = new_target;
                
            elseif ~isempty(obj.current_side) && ~isequal(obj.current_side, new_target.side)
                
                obj.report = sprintf('%s: Flip to same object:  %50s', time, new_target.summary); 
                if obj.debug_bit>1, disp([obj.report ' | ' new_target.details]); end
                
                obj.finish_current(time);
                % stay with current target, do not update prev_target (is this the right thing to do??)
            else
                obj.report = sprintf('%s: Continue observing:   %50s', time, new_target.summary); 
                if obj.debug_bit>1, disp([obj.report ' | ' new_target.details]); end
            end
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.update;
            end
            
        end
        
        function start_current(obj, time)
            
            if nargin<2 || isempty(time)
                time = util.text.time2str(datetime('now', 'TimeZone', 'UTC'));
            elseif isa(time, 'datetime')
                time = util.text.time2str(time); 
            end
            
            if ~isempty(obj.current)
                
                obj.current.start_observation(time); 
                
                % update the coordinates and side of the telescope at beginning of run! 
                obj.current_RA = obj.current.ephem.RA_deg;
                obj.current_Dec = obj.current.ephem.Dec_deg;
                obj.current_side = obj.current.side; 
                
            end
            
            % silently ignore this command if the current target is empty (idle mode)
            
        end
        
        function finish_current(obj, time)
            
            if nargin<2 || isempty(time)
                time = util.text.time2str(datetime('now', 'TimeZone', 'UTC'));
            elseif isa(time, 'datetime')
                time = util.text.time2str(time); 
            end
            
            if isempty(obj.current)
                return; % silently ignore this command if the current target is empty (idle mode)
            end
            
            obj.current.finish_observation(time);
            
            s = obj.current.obs_history(end); 
            s.index = obj.current.index;
            s.name = obj.current.name;
            
            if isempty(obj.obs_history)
               obj.obs_history = s;
            else
                obj.obs_history(end+1) = s;
            end
            
            
            
        end
        
        function val = compare_coordinates(obj, target) % true if current coordinates are close enough to the new target
            
            thresh = 30/60; % let's take 30 arcminutes as the threshold for moving to a new field?
            
            if isempty(obj.current_RA) || isempty(obj.current_Dec)
                val = 0; % no current coordinates, so it can't be close enough
            elseif abs(target.ephem.RA_deg-obj.current_RA)<thresh && abs(target.ephem.Dec_deg-obj.current_Dec)<thresh % both coordinates are close enough
                val = 1; 
            else
                val = 0; 
            end
            
        end
        
    end
    
    methods % simulations
        
        function run_simulation(obj, start_time, end_time)
            
            if nargin<2 || isempty(start_time)
                start_time = datetime('today', 'TimeZone', 'UTC'); 
                start_time.Hour = 13; % set the time to 16:00 Israel time (in winter??) 
            end
            
            if nargin<3 || isempty(end_time)
                end_time = datetime('today', 'TimeZone', 'UTC'); 
                end_time = end_time + days(1); % tomorrow morning! 
                end_time.Hour = 04; % set the time to 7:00 Israel time (in winter??) 
            end
            
            sim_time = start_time; % this is the virtual clock we will use throughout
            
            obj.clear; % start a new night
            
            obj.brake_bit = 0; % start running
            
            % move forward in time without calling update() until reaching night time
            for ii = 1:1e4 % arbitrary timeout
                
                if obj.brake_bit
                    return;
                end
                
                obj.ephem.time = sim_time;
                obj.ephem.updateSun;
                
                if obj.ephem.sun.Alt<obj.max_sun_elevation % sunset! 
                    break; % go to next loop with observations... 
                end
                
                sim_time = sim_time + minutes(obj.sim_time_step); 
                
            end
            
            obj.report = sprintf('Sunset at %s. Starting observations...', sim_time); 
            if obj.debug_bit>1, disp(obj.report); end
            
            obj.current_side = obj.sim_starting_side; % simulations will start with telescope on this side
            
            % move forward in time while observing targets, until sun comes up
            for ii = 1:1e4 % arbitrary timeout
                
                if obj.brake_bit
                    obj.finish_current(sim_time);
                    return;
                end
                
                obj.ephem.time = sim_time;
                obj.ephem.updateSun;
                
                if obj.ephem.sun.Alt>obj.max_sun_elevation % sunrise! 
                    break; % finish observing for tonight
                end
                
                obj.choose(sim_time); 
                
                if isempty(obj.current)
%                     if obj.debug_bit>1, fprintf('No targets available... remain in idle mode\n'); end
                    % do nothing but wait
                elseif obj.need_to_start 
                    
                    if ~isequal(obj.current_side, obj.current.side)
                        sim_time = sim_time + minutes(obj.sim_flip_time);
                    end
                    
                    sim_time = sim_time + minutes(obj.sim_slew_time);
                    sim_time = sim_time + minutes(obj.sim_focus_time);
                    
                    obj.start_current(sim_time); 

                    if ~isempty(obj.gui)
                        obj.gui.update;
                    end
                    
                end
                
                sim_time = sim_time + minutes(obj.sim_time_step); 
                
                pause(obj.sim_plot_pause); % pause to make the simulation go a little slower for plotting/visualization
                
            end
            
            obj.finish_current(sim_time); 
            
            obj.brake_bit = 1;
            
            if ~isempty(obj.gui)
                obj.gui.update;
            end
            
            % need to call anything else to wrap up?
            
            obj.report = sprintf('Sunrise at %s. Finished observations...', sim_time);
            if obj.debug_bit>1, disp(obj.report); end

        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis'); 
            input.input_var('font_size', 20); 
            input.input_var('marker_size', 18); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.ax)
                
                if ~isempty(obj.gui) && obj.gui.check
                    input.ax = obj.gui.axes_image;
                else
                    input.ax = gca;
                end
                
            end
            
            obj.map.LST = obj.ephem.LST_deg/15;
            
            obj.map.show('vector', [-5 5]); 
            
            title(input.ax, ''); 
            input.ax.FontSize = input.font_size;
            
            input.ax.NextPlot = 'add';
            
            % plot current position
            if ~isempty(obj.current) && isvalid(obj.current) && ~isempty(obj.current_RA) && ~isempty(obj.current_Dec) && ~isempty(obj.current_side)
                s = struct('RA_deg', obj.current_RA, 'Dec_deg', obj.current_Dec, 'side', obj.current_side, ...
                    'start_time', obj.current.start_time, 'end_time', '', 'name', obj.current.name); 
                
                obj.plotObservation(s, input.ax, input.marker_size); 

                plot(input.ax, obj.current_RA/15, obj.current_Dec, 'yo', 'MarkerSize', input.marker_size); 
            end
            
            % plot previous targets
            for ii = 1:length(obj.obs_history)
                obj.plotObservation(obj.obs_history(ii), input.ax, input.marker_size); 
            end
            
            
            input.ax.NextPlot = 'replace';
            
        end
        
        function plotObservation(obj, s, ax, marker_size)
            
            if strcmp(s.side, 'West')
                plot(ax, s.RA_deg/15, s.Dec_deg, 'gx', 'MarkerSize', marker_size); 
            else
                plot(ax, s.RA_deg/15, s.Dec_deg, 'g+', 'MarkerSize', marker_size); 
            end
            
            text(ax, s.RA_deg/15+0.2, s.Dec_deg+5, strrep(s.name, '_', ' '), 'Color', 'g', 'FontSize', 12); 
            
            % at some point we will have to figure out a way to put multiple printouts for the same field, beyond East/West...             
            shift_down = 0;
            if strcmp(s.side, 'East')
                shift_down = 5;
            end
            
            if ~isempty(s.start_time) && ~isempty(s.end_time)
                text(ax, s.RA_deg/15+0.3, s.Dec_deg-shift_down, [s.start_time(12:16) '-' s.end_time(12:16)], 'Color', 'g', 'FontSize', 12); 
            elseif ~isempty(s.start_time)
                text(ax, s.RA_deg/15+0.3, s.Dec_deg-shift_down, s.start_time(12:16), 'Color', 'y', 'FontSize', 12); 
            end
            
        end
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = obs.sched.gui.SchedGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
    end    
    
end

