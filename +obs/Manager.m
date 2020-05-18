classdef (CaseInsensitiveProperties, TruncatedProperties) Manager < handle
% Top level class to control observatory. 
% This class has 3 main roles: 
% (1) To contain all hardware objects (e.g., dome, mount, weather station).
% (2) To make various checks and shutdown the observatory if needed. 
% (3) To open the observatory and start operations via scheduler
%       (THIS IS NOT YET IMPLEMENTED!)
%
% Main features:
%   -dome, and mount are objects connected to hardware for main operations. 
%    If one of these devices fails, the observatory must shut down (or make 
%    a call to get help). 
%   -weather, wind, etc: sensors that check the conditions are viable for 
%    observation. 
%   -checker: a SensorChecker object that collects all the weather sensor 
%    information and makes a summary of the results for the manager. 
%   -timers: each timer has a different time scale for checking the state 
%            of the observatory. Timer t1 is each minute, collecting the 
%            weather data only. Timer t2 runs every five minutes and 
%            decides if the weather conditions are ok and closes/opens 
%            accordingly (also checks t1 is alive). 
%            Time t3 only checks that t2 is alive every half an hour.  
% 
% PLEASE READ THE COMMENTS ON PROPERTIES FOR MORE DETAIL!
    
    
    properties(Transient=true)
        
        t1; % quick timer (every minute or so) just to update sensors/devices
        t2; % check that everything is connected and that weather is good, print to log file
        t3; % verify that the other two are still running (every half an hour or so)
        
        gui;
        
    end
    
    properties % objects
        
        log@util.sys.Logger;
        
        checker@obs.SensorChecker;
        
        dome; % AstroHaven dome
        mount; % ASA mount
        
        weather; % Boltwood weather station
        wind; % windETH sensor
        humidity; % humidity/temperature dog
        temperature; % additional temperature meters
        
        email@obs.comm.Email; 
        
        cam_pc@obs.comm.PcSync; % communications object to camera PC
        
        obs_log; % get this from cam_pc
        
    end
    
    properties % switches/controls
        
        use_maintenance_mode = 0; % this disables the function in t3 that re-enables "use_shutdown" every 30 minutes (use with care!!!)
        use_shutdown = 1; % when this is enabled, observatory shuts down on bad weather/device failure
        use_startup = 0; % when this is enabled, observatory opens up and starts working by itself! 
        
        % use these to override these devices/sensors
        use_dome = 1; % override AstroHaven dome
        use_mount = 1; % ovrride ASA mount
        use_weather = 1; % override Boltwood weather station
        use_wind = 1; % override windETH
        use_humidity = 0; % override humidity dog
        use_temperature = 0; % override other temp sensors
        
        period1 = 60; % time between updates of all devices/sensors
        period2 = 300; % time for equipment/weather check and log file
        period3 = 1800; % time for verifying shorter timers are working (and other tests?)
        
        autostart_min_weather_samples = 4;
        autostart_min_weather_minutes = 20;
        
        brake_bit = 1; % also set the brake bit for mount/dome
        debug_bit = 1;
        
    end
    
    properties % inputs/outputs
                
        devices_ok = 1; % no critical failures in mount/dome/etc
        devices_report = 'OK'; % if failure happens, specify the reason
        
    end
    
    properties(Dependent=true)
        
        sensors_ok; % all weather sensors give good results
        sensors_report; % if bad weather, report it here
        report; % general report: devices, sensors, shutdown state
        
        RA;  % shortcut to mount telRA
        DEC; % shortcut to mount telDEC
        ALT; % shortcut to mount telALT
        LST; % shortcut to mount LST
        
        tracking; % shortcut to mount tracking
        
        is_shutdown; % will be 1 when dome is closed and mount not tracking (add mount in parking at some point)
        
    end
    
    properties(Hidden=true)
       
        latest_email_report_date = ''; % keep track of the last time we sent this, so we don't send multiple emails each day
        
        sensor_ok_history; % a vector of datetime objects, for each time we have measured good weather (this gets reset upon a single bad weather measurement)
        latest_email_autostart_date = '';
        
        
        version = 1.03;
        
    end
    
    methods % constructor
        
        function obj = Manager(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.Manager')
                if obj.debug_bit, fprintf('Manager copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Manager constructor v%4.2f\n', obj.version); end
                
                obj.log = util.sys.Logger('Top_level_manager'); % keep track of commands given and errors received... 
                
                obj.connect; % connect to all hardware and objects
                
            end
            
        end
        
        function connect(obj) % connect to all hardware and objects
            
            if obj.use_dome
                obj.connectDome; % AstroHaven dome
            end
            
            if obj.use_mount
                try
                obj.connectMount; % ASA mount
                catch ME
                    warning(ME.getReport); 
                end
            end
            
            if obj.use_weather
                obj.connectBoltwood; 
            end
            
            if obj.use_wind
                obj.connectWindETH;
            end
            
            % add additional devices
            % ...
            
            obj.connectSensorChecker; % create checker object that collects sensor data
            
            obj.constructPcSync;
            obj.constructEmail;
            
            % start the 3 layers of timers
            obj.setup_t3; % check t2 is alive (half hour period)
            obj.setup_t2; % check devices, collect weather, decide on closing dome, make log report (5 minute period)
            obj.setup_t1; % get sensors to measure weather data for averaging (1 minute period)
                
            
        end
        
        function delete(obj) % destructor
            
            % make sure not to leave hanging timers
            delete(obj.t3);
            delete(obj.t2);
            delete(obj.t1);
            
        end
        
        function connectDome(obj)
            
            obj.log.input('Connecting to dome.');
            
            try 
                
                obj.dome = obs.dome.AstroHaven;
                obj.dome.owner = obj;
                
            catch ME
                
                obj.log.error(ME.getReport);
                obj.log.input('Connecting dome simulator.');
                
                warning(ME.getReport);
                
                disp('Cannot connect to AstroHaven dome. Using simulator instead...');
                
                try
                    obj.dome = obs.dome.Simulator;
                catch ME
                    obj.log.error(ME.getReport);
                    warning(ME.getReport);
                end
                
            end
            
        end
        
        function connectMount(obj)
            
            obj.log.input('Connecting to mount.');
            
            try 
                
                obj.mount = obs.mount.ASA;
                obj.mount.owner = obj;
                
            catch ME
                
                obj.log.error(ME.getReport);
                rethrow(ME);
                
            end
            
        end
        
        function connectBoltwood(obj)
            
            obj.log.input('Connecting to Boltwood weather station.');
            
            try 
                obj.weather = obs.sens.Boltwood;
            catch ME
                
                obj.log.error(ME.getReport);
                obj.log.input('Connecting to weather simulator.');
                warning(ME.getReport);
                
                disp('Cannot connect to Boltwood weather station. Using simulator instead...');
                
                try
                    obj.weather = obs.sens.Simulator;
                catch ME
                    obj.log.error(ME.getReport);
                    warning(ME.getReport);
                end
                
            end
            
        end
        
        function connectWindETH(obj)
            
            obj.log.input('Connecting to WindETH');
            
            try 
                obj.wind = obs.sens.WindETH;
            catch ME
                
                
                obj.log.error(ME.getReport);
                
                warning(ME.getReport);
                
                disp('Cannot connect to WindETH sensor.');
                
            end
            
        end
        
        function connectSensorChecker(obj)
            
            try
                
                obj.checker = obs.SensorChecker(obj); % checks the weather using all sensors
                
            catch ME
                obj.log.error(ME.getReport);
                warning(ME.getReport);
                disp('Cannot initialize SensorChecker!');
            end

        end
        
        function constructPcSync(obj)
            
            try
                
                obj.cam_pc = obs.comm.PcSync('client');
                
            catch ME
                obj.log.error(ME.getReport);
                warning(ME.getReport);
                disp('Cannot create a PcSync tcp/ip object')
            end
            
        end
        
        function constructEmail(obj)
            
            obj.email = obs.comm.Email;
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            % restart the timers
            obj.setup_t3;
            obj.setup_t2;
            obj.setup_t1;
            
        end
        
    end
    
    methods % getters
        
        function val = get.sensors_ok(obj) % shortcut to checker, that is updated on t2
            
            if isempty(obj.checker)
                val = [];
            else
                val = obj.checker.sensors_ok;
            end
            
        end
        
        function val = get.sensors_report(obj) % shortcut to checker, that is updated on t2
            
            if isempty(obj.checker)
                val = '';
            else
                val = obj.checker.report;
            end
            
        end
        
        function val = get.report(obj) % composite report: sensors, devices, shutdown state
            
            val = sprintf('Sensors: %s | Devices: %s | state: %s', obj.sensors_report, obj.devices_report, obj.observatory_state);
            
        end
        
        function val = get.is_shutdown(obj) % true when dome is closed and mount is not tracking
            
            val = 0;
            
            try 
            
                if ~isempty(obj.mount) && obj.mount.status && obj.mount.tracking
                    val = 0;
                    return;
                end
                
                if obj.dome.is_closed==0
                    val = 0;
                    return;
                end
                    
                val = 1; % if all checks are passed, value can be 1
                
            catch ME
                obj.log.error(ME.getReport);
                warning(ME.getReport);
            end
            
        end
        
        function val = observatory_state(obj) % string: SHUT or OPEN
            
            if obj.is_shutdown
                val = 'SHUT';
            else
                val = 'OPEN';
            end
            
        end
        
        function val = get.RA(obj) % shortcut to telescope 
            
            if ~isempty(obj.mount)
                val = obj.mount.telRA;
            else
                val = [];
            end
            
        end
        
        function val = get.DEC(obj) % shortcut to telescope 
            
            if ~isempty(obj.mount)
                val = obj.mount.telDEC;
            else
                val = [];
            end
            
        end
        
        function val = get.ALT(obj) % shortcut to telescope 
            
            if ~isempty(obj.mount)
                val = round(obj.mount.telALT);
            else
                val = [];
            end
            
        end
        
        function val = get.LST(obj) % shortcut to telescope 
            
            if ~isempty(obj.mount)
                val = obj.mount.LST;
            else
                val = [];
            end
            
        end
        
        function val = get.tracking(obj) % shortcut to telescope
            
            if isempty(obj.mount)
                val = [];
            else
                val = obj.mount.tracking;
            end
            
        end
        
        function val = average_temperature(obj) % average of all sensors that can measure this
            
            if ~isfield(obj.checker.temperature, 'func') || isempty(obj.checker.temperature.func)
                val = nanmean(obj.checker.temperature.now);
            else
                val = util.stat.stat_eval(obj.checker.temperature.func, obj.checker.temperature.now, 2); 
            end
            
        end
        
        function val = average_clouds(obj) % average of all sensors that can measure this
            
            if ~isfield(obj.checker.clouds, 'func') || isempty(obj.checker.clouds.func)
                val = nanmean(obj.checker.clouds.now);
            else
                val = util.stat.stat_eval(obj.checker.clouds.func, obj.checker.clouds.now, 2); 
            end
            
        end
        
        function val = average_light(obj) % average of all sensors that can measure this
            
            if ~isfield(obj.checker.light, 'func') || isempty(obj.checker.light.func)
                val = nanmean(obj.checker.light.now);
            else
                val = util.stat.stat_eval(obj.checker.light.func, obj.checker.light.now, 2); 
            end
            
        end
        
        function val = average_wind_speed(obj) % average of all sensors that can measure this
            
            if ~isfield(obj.checker.wind_speed, 'func') || isempty(obj.checker.wind_speed.func)
                val = nanmean(obj.checker.wind_speed.now);
            else
                val = util.stat.stat_eval(obj.checker.wind_speed.func, obj.checker.wind_speed.now, 2); 
            end
            
        end
        
        function val = average_wind_dir(obj) % average of all sensors that can measure this
            
            if ~isfield(obj.checker.wind_dir, 'func') || isempty(obj.checker.wind_dir.func)
                val = nanmean(obj.checker.wind_dir.now);
            else
                val = util.stat.stat_eval(obj.checker.wind_dir.func, obj.checker.wind_dir.now, 2); 
            end
            
        end
        
        function val = average_humidity(obj) % average of all sensors that can measure this
            
            if ~isfield(obj.checker.humidity, 'func') || isempty(obj.checker.humidity.func)
                val = nanmean(obj.checker.humidity.now);
            else
                val = util.stat.stat_eval(obj.checker.humidity.func, obj.checker.humidity.now, 2); 
            end
            
        end
        
        function val = average_pressure(obj) % average of all sensors that can measure this
            
            if ~isfield(obj.checker.pressure, 'func') || isempty(obj.checker.pressure.func)
                val = nanmean(obj.checker.pressure.now);
            else
                val = util.stat.stat_eval(obj.checker.pressure.func, obj.checker.pressure.now, 2); 
            end
        
        end
        
        function val = any_rain(obj) % check if any rain sensor is reporting rain
            
            val = any(obj.checker.rain.now); 
            
        end
        
        function val = areTimersRunning(obj) % quick way to verify timers are still on... 
            
            val(1) = strcmp(obj.t1.Running, 'on');
            val(2) = strcmp(obj.t2.Running, 'on');
            val(3) = strcmp(obj.t3.Running, 'on');
            
        end
        
    end
    
    methods % setters
        
        function set.cam_pc(obj, val) % share the sync object with the telescope
            
            obj.cam_pc = val;
            
            if ~isempty(obj.mount) && isa(obj.mount, 'obs.mount.ASA')
                obj.mount.cam_pc = val;
            end
            
            if ~isempty(obj.dome) && isa(obj.dome, 'obs.dome.AstroHaven')
                obj.dome.cam_pc = val;
            end
            
        end
        
        function set.tracking(obj, val) % shortcut to mount
            
            obj.mount.tracking = val;
            
        end
        
    end
    
    methods % timer related
        
        function stop_timers(obj) % stop all three timers (for debugging only!) make sure to turn them back on! 
            
            obj.stop_t3;
            obj.stop_t2;
            obj.stop_t1;
            
        end
        
        function start_timers(obj) % start all timers
            
            obj.setup_t3;
            obj.setup_t2;
            obj.setup_t1;
            
        end
        
        function callback_t1(obj, ~, ~) % update sensors and GUI
            
            try 
            
                if isempty(obj.checker)
                    obj.connectSensorChecker;
                end
                
                obj.checker.update; % go over all sensors and only tell them to collect data. It's reported back in t2
                
                obj.updateCameraComputer;
                
                if ~isempty(obj.gui) && obj.gui.check
                    obj.gui.update;
                end
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function setup_t1(obj, ~, ~) % start the timer t1 with period1
            
            if ~isempty(obj.t1) && isa(obj.t1, 'timer') && isvalid(obj.t1)
                if strcmp(obj.t1.Running, 'on')
                    stop(obj.t1);
                    delete(obj.t1);
                    obj.t1 = [];
                end
            end
            
            delete(timerfind('name', 'Status-check-t1'));
            
            obj.t1 = timer('BusyMode', 'queue', 'ExecutionMode', 'fixedRate', 'Name', 'Status-check-t1', ...
                'Period', obj.period1, 'StartDelay', obj.period1, ...
                'TimerFcn', @obj.callback_t1, 'ErrorFcn', @obj.setup_t1);
            
            start(obj.t1);
            
        end
        
        function stop_t1(obj, ~, ~) % stop the timer t1
            
            stop(obj.t1);
            
        end 
        
        function callback_t2(obj, ~, ~) % collect (averaged) sensor data and check devices are all ok
            
            try

                % make sure t1 is running! 
                if isempty(obj.t1) || ~isvalid(obj.t1) || strcmp(obj.t1.Running, 'off')
                    obj.setup_t1;
                end

                obj.update; % check devices are functioning
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function setup_t2(obj, ~, ~) % start the timer t2 with period2
            
            if ~isempty(obj.t2) && isa(obj.t2, 'timer') && isvalid(obj.t2)
                if strcmp(obj.t2.Running, 'on')
                    stop(obj.t2);
                    delete(obj.t2);
                    obj.t2 = [];
                end
            end
            
            delete(timerfind('name', 'Status-check-t2'));
            
            obj.t2 = timer('BusyMode', 'queue', 'ExecutionMode', 'fixedRate', 'Name', 'Status-check-t2', ...
                'Period', obj.period2, 'StartDelay', obj.period2, ...
                'TimerFcn', @obj.callback_t2, 'ErrorFcn', @obj.setup_t2);
            
            start(obj.t2);
            
        end
        
        function stop_t2(obj, ~, ~) % stop timer t2
            
            stop(obj.t2);
            
        end 
        
        function callback_t3(obj, ~, ~) % make sure t2 is running and re-enables "use_shutdown"
            
            try

                if obj.use_maintenance_mode==0
                    obj.use_shutdown = 1;
                    obj.checker.use_twilight_mode = 0; 
                end
                
                % make sure t2 is running! 
                if isempty(obj.t2) || ~isvalid(obj.t2) || strcmp(obj.t2.Running, 'off')
                    obj.setup_t2;
                end
                
                % morning report! 
                t = datetime('now', 'TimeZone', 'UTC');
                
                if t.Hour>4 % this is in UTC, so we can translate it to 6 or 6 local time

                    d = datestr(t - days(1), 'yyyy-mm-dd');

                    if isempty(obj.latest_email_report_date) || ~strcmp(obj.latest_email_report_date, d)
                        obj.morning_report;
                    end

                end

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function setup_t3(obj, ~, ~) % start the timer t3 with period3
            
            if ~isempty(obj.t3) && isa(obj.t3, 'timer') && isvalid(obj.t3)
                if strcmp(obj.t3.Running, 'on')
                    stop(obj.t3);
                    delete(obj.t3);
                    obj.t3 = [];
                end
            end
            
            delete(timerfind('name', 'Status-check-t3'));
            
            obj.t3 = timer('BusyMode', 'queue', 'ExecutionMode', 'fixedRate', 'Name', 'Status-check-t3', ...
                'Period', obj.period3, 'StartDelay', obj.period3, ...
                'TimerFcn', @obj.callback_t3, 'ErrorFcn', @obj.setup_t3);
            
            start(obj.t3);
            
        end
        
        function stop_t3(obj, ~, ~) % stop timer t3
            
            stop(obj.t3);
            
        end
        
        function morning_report(obj)
            
            t = datetime('yesterday', 'TimeZone', 'UTC');
            d = datestr(t, 'yyyy-mm-dd');
            
            fprintf('%s: sending morning report by Email!\n', t); 
            
            str = '';
            
            str = sprintf('%s\n%s', str, obj.email.html(sprintf('This is a morning report for the night of %s.', d), 'p', 'font-size:14px')); 
            
%             str = sprintf('%s\n<p style="font-size:14px">', str); 
%             str = sprintf('%s\nThis is a morning report for the night of %s.', str, d);
%             str = sprintf('%s\n</p>', str); 
            
            str = sprintf('%s\n%s', str, obj.email.html(sprintf('Dome status is: %s ', obj.observatory_state), 'p', 'font-size:18px; font-weight:bold')); 

%             str = sprintf('%s\n<p style="font-size:18px;font-weight:bold;">', str); 
%             str = sprintf('%s\n Dome status is: %s ', str, obj.observatory_state);
%             str = sprintf('%s\n</p>', str); 
            
            if ~isempty(obj.obs_log) && isfield(obj.obs_log, 'date') && strcmp(obj.obs_log.date, d) && length(fields(obj.obs_log))>1
                
                list = fields(obj.obs_log); 
                
                obs_str = sprintf('Runs overview for tonight: \n <table style="width:60%%">');
                
                for ii = 1:length(list) % each target has a few runs inside a struct array with this name
                    
                    if strcmp(list{ii}, 'date'), continue; end
                    
                    s = obj.obs_log.(list{ii}); % struct with some info on the run
                    files = 0;
                    time = 0; 
                    
                    for jj = 1:length(s) % go over the different structs for each run 
                        
                        files = files + s(jj).num_files;
                        time = time + s(jj).runtime;
                        
                    end
                    
%                     obs_str = sprintf('%s\n%10s (%d) | files: %4d | runtime: %7.1f ', obs_str, list{ii}, length(s), files, time); 
                    obs_str = sprintf('%s\n <tr> <td> %s (%d) </td><td> files: %4d </td><td> runtime: %7.1f </td></tr>', obs_str, list{ii}, length(s), files, time); 
                    
                end
                
                
                
            else
                obs_str = sprintf('Could not find any information on observations run tonight...'); 
            end
            
%             str = sprintf('%s\n<p> %s \n</p>\n', str, obs_str); % add the observation report to the string
            
            str = [str obj.email.html(obs_str)]; 

            if isfield(obj.cam_pc.incoming, 'drives') && ~isempty(obj.cam_pc.incoming.drives)
            
                drive_str = sprintf('Hard drive space overview:\n'); 
                
                s = obj.cam_pc.incoming.drives;
                f = fields(s); 
                
                gb_limit = 3000; % set this to warn when HDDs are running low...  
                
                drive_str = sprintf('%s\n<table style="width:60%%;font-family:Courier;font-size:12px">\n <tr>', drive_str);
                
                for ii = 1:length(f)
                    drive_str = sprintf('%s\n      <td style="font-family:Courier"> %s:\\ </td> ', drive_str, f{ii}); 
                end
                
                drive_str = sprintf('%s\n  </tr><tr>\n', drive_str);
                
                for ii = 1:length(f)
                    if ( strcmp(f{ii}, 'E') || strcmp(f{ii}, 'F') ) && ~isempty(s.(f{ii})) && s.(f{ii})<gb_limit
                        drive_str = sprintf('%s\n    <td style="background-color:red"> %d Gb </td> ', drive_str, round(s.(f{ii}))); 
                    else
                        drive_str = sprintf('%s\n    <td style="background-color:none"> %d Gb </td> ', drive_str, round(s.(f{ii}))); 
                    end
                end
                
                drive_str = sprintf('%s\n</tr><br>', drive_str);
                
                if isfield(s, 'E') && ~isempty(s.E) && s.E<gb_limit
                    drive_str = sprintf('%s\nWARNING: Drive E is has less than %d Gb left!', drive_str, gb_limit);
                end
                
                if isfield(s, 'F') && ~isempty(s.F) && s.F<gb_limit
                    drive_str = sprintf('%s\nWARNING: Drive F is has less than %d Gb left!', drive_str, gb_limit);
                end
                
%                 str = sprintf('%s\n <p style="font-familiy:Courier;font-size:12px"> %s </p>\n', str, drive_str); 
                str = [str obj.email.html(drive_str, 'p', 'font-family:Courier; font-size:12px;')]; 
                
            end
            
%             obj.email.sendToList('subject', sprintf('[WFAST] Morning report %s %d', d, randi(1000)), 'text', str, 'header', 1, 'footer', 1, 'html', 1); 
            obj.email.sendToList('subject', sprintf('[WFAST] Morning report %s', d), 'text', str, 'header', 1, 'footer', 1, 'html', 1); 
            
            obj.latest_email_report_date = d; % make sure we don't resend this email today (after a successful send!)
            
        end
        
    end
    
    methods % calculations / commands
        
        function update(obj) % check devices and sensors and make a decision if to shut down the observatory

            obj.updateDevices; % runs update() for each critical device (mount, dome) and checks its status

            obj.checker.decision_all; % collect weather data and make a decision

            obj.log.input(obj.report); % summary of observatory status

            if ~isempty(obj.cam_pc) 
                
                if ~isempty(obj.mount.tracking) && obj.mount.tracking && obj.dome.is_closed==0
                    obj.cam_pc.outgoing.stop_camera = 0; % if everything is cool, let the camera keep going
                end
                
                obj.cam_pc.update;
                
            end
            
            % the following conditions for shutting down
            if obj.use_shutdown && obj.devices_ok==0 % critical device failure, must shut down
                if obj.is_shutdown==0 % if already shut down, don't need to do it again
                    fprintf('%s: Device problems... %s \n', datestr(obj.log.time), obj.devices_report); 
                    obj.shutdown;
                end
            end

            if obj.use_shutdown && obj.sensors_ok==0 % one of the sensors reports bad weather, must shut down
                if obj.is_shutdown==0 % if already shut down, don't need to do it again
                    fprintf('%s: Bad weather... %s \n', datestr(obj.log.time), obj.checker.report); 
                    obj.shutdown;
                end
            end
            
            if obj.use_shutdown && obj.checkDayTime % check if the system clock says it is day time
                if obj.is_shutdown==0 % if already shut down, don't need to do it again
                    fprintf('%s: System clock says it is day time... %s \n', datestr(obj.log.time), obj.checker.report); 
                    obj.shutdown;
                end
            end
            
            if obj.is_shutdown
                obj.latest_email_autostart_date = ''; % once shut down, we need once again to send an email if we are ready to open
            end
            
            if obj.use_startup
                
                if obj.sensors_ok
                    
                    t = datetime('now', 'TimeZone', 'UTC');
                    
                    if obj.checkPrevWeather
                        
                        if isempty(obj.latest_email_autostart_date) 
                            obj.email.sendToList('subject', ['Ready to open dome ' util.text.time2str(t)], ...
                                'text', sprintf('Weather data looks good for at least %d minutes. Consider opening for observations. ', ...
                                round(minutes(t-obj.sensor_ok_history(end)))));
                            
                            obj.latest_email_autostart_date = t;
                            
                        end
                        
                    end
                    
                    obj.sensor_ok_history = vertcat(datetime('now', 'TimeZone', 'UTC')); 
                    
                else
                    obj.sensor_ok_history = [];
                end
                
            end
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.update;
            end
            
        end
                
        function updateDevices(obj) % run "update" for each device and see if the status is still good
        
            obj.devices_ok = 1;
            obj.devices_report = 'OK';
            
            if obj.use_dome
                obj.dome.update;
                if obj.dome.status==0
                    obj.devices_ok = 0;
                    obj.devices_report = 'Dome error!';
                    return;
                end
            end
            
            if obj.use_mount
                obj.mount.update;
                if obj.mount.status==0
                    obj.devices_ok = 0;
                    obj.devices_report = 'Mount error!';
                    return;
                end
            end
            
            % add maybe checks for boltwood if we think it is critical?
            
        end
        
        function val = checkDayTime(obj) % return true if it is day time according to the system clock
            
            time = datetime('now', 'TimeZone', 'Asia/Jerusalem');
            
            if time.Hour>7 && time.Hour<16
                val = 1;
            else
                val = 0;
            end
            
        end
        
        function val = checkPrevWeather(obj)
           
            new_hist = datetime('now', 'TimeZone', 'UTC');

            % collect to a new history all the times that are
            % spaced more than 4 minutes from the last accepted entry
            for ii = length(obj.sensor_ok_history):-1:1
                if new_hist(end)-obj.sensor_ok_history(ii)>4*60 % measurements less than 4 minutes apart are ignored
                    new_hist(end+1) = obj.sensor_ok_history(ii); 
                end
            end

            % new_hist should now contain only well spaced times when weather was ok
            if numel(new_hist)>=obj.autostart_min_weather_samples && ...
                minutes(new_hist(end)-new_hist(1))<obj.autostart_min_weather_minutes
                
                val = 1;
                
            else
                val = 0; 
            end
            
        end
        
        function updateCameraComputer(obj) % send updates through the PcSync object to the camera-PC
            
            try 
                
                if ~isfield(obj.cam_pc.incoming, 'time') 
                    % do something like try to reconnect
                    obj.cam_pc.connect;
                end
                
                t_in = util.text.str2time(obj.cam_pc.incoming.time);
                t_now = datetime('now', 'TimeZone', 'UTC'); 
                
                if minutes(t_now-t_in)>10
                    % do something like try to reconnect
                    disp('input from "cam_pc" is out of date by more than 10 minutes!'); 
                    obj.cam_pc.connect;
                end
                
                % check that commands are being echoed back
                if isfield(obj.cam_pc.outgoing, 'command_str') && isfield(obj.cam_pc.outgoing, 'command_time')
                    
                    % if command is successfully echoed, no need to stop sending it
                    if isfield(obj.cam_pc.incoming, 'echo_str') && isfield(obj.cam_pc.incoming, 'echo_time')
                        
                        if strcmp(obj.cam_pc.outgoing.command_str, obj.cam_pc.incoming.echo_str) && ...
                            strcmp(obj.cam_pc.outgoing.command_time, obj.cam_pc.incoming.echo_time)
                            
                            obj.cam_pc.outgoing.command_str = ''; 
                            obj.cam_pc.outgoing.command_time = ''; 
                            % on the receiving side, empty command is ignored
                            
                        end
                        
                    end
                    
                end
                
                if ~obj.cam_pc.is_connected || ~obj.cam_pc.status
%                     obj.cam_pc.connect;
                end

            catch 
                t = datetime('now', 'TimeZone', 'UTC'); 
                fprintf('%s: Failed to connect to camera computer\n', t); 
                % do nothing, as we can be waiting for ever for server to connect
            end
            
            % obj.cam_pc.outgoing.OBJECT = ???
            if ~isempty(obj.mount) && obj.use_mount
                
                if isempty(obj.mount.cam_pc)
                    obj.mount.cam_pc = obj.cam_pc; % share the handle to this object
                end
                
                obj.mount.updateCamera; % does mount have any useful data to share with camera? 
                
            end
            
            obj.cam_pc.outgoing.TEMP_OUT = obj.average_temperature;
            obj.cam_pc.outgoing.WIND_DIR = obj.average_wind_dir;
            obj.cam_pc.outgoing.WIND_SPEED = obj.average_wind_speed;
            obj.cam_pc.outgoing.HUMID_OUT = obj.average_humidity;
            obj.cam_pc.outgoing.LIGHT = obj.average_light; 
            obj.cam_pc.outgoing.PRESSURE = obj.average_pressure;             
            
%             if obj.dome.is_closed
%                 
%                 if ~isfield(obj.cam_pc.outgoing, 'stop_camera') || obj.cam_pc.outgoing.stop_camera==0
%                     obj.log.input('Dome closed, sending camera stop command');
%                     disp(obj.log.report);
%                 end
%                 
%                 obj.cam_pc.outgoing.stop_camera = 1;
%             else
%                 obj.cam_pc.outgoing.stop_camera = 0;
%             end
            
            
            % add additional parameters and some commands like "start run"
            
            obj.cam_pc.update;
            
            % get the observation log from the camera, including how long each target was observed
            if ~isempty(obj.cam_pc.incoming) && isfield(obj.cam_pc.incoming, 'obs_log')
                obj.obs_log = obj.cam_pc.incoming.obs_log; 
            end
            
        end
        
        function commandCamPC(obj, str, cam_mode, exp_time)
            
            if nargin<3 || isempty(cam_mode)
                cam_mode = '';
            end
            
            if nargin<4 || isempty(exp_time)
                exp_time = [];
            end
            
            obj.cam_pc.outgoing.command_str = str;
            obj.cam_pc.outgoing.command_time = util.text.time2str(datetime('now', 'TimeZone', 'UTC')); 
            obj.cam_pc.outgoing.cam_mode = cam_mode;
            obj.cam_pc.outgoing.exp_time = exp_time;
            obj.cam_pc.update; 
            
        end
        
        function closeDome(obj) % shortcut to closing both sides of the dome
            
            obj.dome.closeBothFull;
            
        end
        
        function shutdown(obj) % command to shut down observatory (close dome, stop tracking)
            
            obj.log.input('Shutting down observatory!');
            
            disp([char(obj.log.time) ': Shutting down observatory!']);
            
            try 

                obj.stop; % stop any slews or shutter motion

                obj.mount.tracking = 0; % later add command to park the telescope? 

                obj.closeDome;

                obj.cam_pc.outgoing.stop_camera = 1; % make sure camera stops running also
                obj.cam_pc.update;
                
                % anything else we can do to put the dome to shutdown mode?
                % ...

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function stop(obj) % stop dome and mount motion
            
            obj.brake_bit = 1;
            
            obj.dome.stop;
            obj.mount.stop; 
            
            % verify other objects that may need a brake...
            list = properties(obj);
            
            for ii = 1:length(list)
                
                if isobject(obj.(list{ii})) && isprop(obj.(list{ii}), 'brake_bit')
                    obj.(list{ii}).brake_bit = 1;
                end
                
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = obs.gui.ManagerGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
    end    
    
end

