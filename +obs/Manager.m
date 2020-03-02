classdef (CaseInsensitiveProperties, TruncatedProperties) Manager < handle
% Top level class to control observatory. 
% This class has 2 main roles: 
% (1) To contain all hardware objects (e.g., dome, mount, weather station).
% (2) To make various checks and shutdown the observatory if needed. 
%
% Main features:
%   -dome, and mount are objects connected to hardware for main operations. 
%    If one of these devices fails, the observatory must shut down (or make 
%    a call to get help. 
%   -weather, wind, etc: sensors that check the conditions are viable for 
%    observation. 
%   -checker: a SensorChecker object that collects all the weather sensor 
%    information and makes a summary of the results for the manager. 
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
        
        sync@obs.comm.PcSync;
        
    end
    
    properties % switches/controls
        
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
        
        brake_bit = 1; % also set the brake bit for mount/dome
        debug_bit = 1;
        
    end
    
    properties % inputs/outputs
                
        devices_ok = 1; % no critical failures in mount/dome/etc
        devices_report = 'OK'; % if failure happens, specify it here
        
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
       
        version = 1.01;
        
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
%                 obj.connectMount; % ASA mount
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
                
                obj.sync = obs.comm.PcSync('client');
                
            catch ME
                obj.log.error(ME.getReport);
                warning(ME.getReport);
                disp('Cannot create a PcSync tcp/ip object')
            end
            
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
            
                if isempty(obj.mount) || isempty(obj.mount.tracking) || obj.mount.tracking
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
        
        function val = get.tracking(obj)
            
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
        
        function val = average_pressure(obj)
            
            if ~isfield(obj.checker.pressure, 'func') || isempty(obj.checker.pressure.func)
                val = nanmean(obj.checker.pressure.now);
            else
                val = util.stat.stat_eval(obj.checker.pressure.func, obj.checker.pressure.now, 2); 
            end
        
        end
        
        function val = any_rain(obj)
            
            val = any(obj.checker.rain.now); 
            
        end
        
        function val = areTimersRunning(obj) % quick way to verify timers are still on... 
            
            val(1) = strcmp(obj.t1.Running, 'on');
            val(2) = strcmp(obj.t2.Running, 'on');
            val(3) = strcmp(obj.t3.Running, 'on');
            
        end
        
    end
    
    methods % setters
        
        function set.sync(obj, val)
            
            obj.sync = val;
            obj.mount.sync = val;
            obj.dome.sync = val;
            
        end
        
        function set.tracking(obj, val)
            
            obj.mount.tracking = val;
            
        end
        
    end
    
    methods % timer related
        
        function stop_timers(obj)
            
            obj.stop_t3;
            obj.stop_t2;
            obj.stop_t1;
            
        end
        
        function start_timers(obj)
            
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
        
        function setup_t1(obj, ~, ~)
            
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
        
        function stop_t1(obj, ~, ~)
            
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
        
        function setup_t2(obj, ~, ~)
            
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
        
        function stop_t2(obj, ~, ~)
            
            stop(obj.t2);
            
        end 
        
        function callback_t3(obj, ~, ~) % make sure t2 is running (are there any other checks?)
            
            try

                % make sure t2 is running! 
                if isempty(obj.t2) || ~isvalid(obj.t2) || strcmp(obj.t2.Running, 'off')
                    obj.setup_t2;
                end

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function setup_t3(obj, ~, ~)
            
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
        
        function stop_t3(obj, ~, ~)
            
            stop(obj.t3);
            
        end
        
    end
    
    methods % calculations / commands
        
        function update(obj) % check devices and sensors and make a decision if to shut down the observatory

            obj.updateDevices;

            obj.checker.decision_all; % collect weather data and make a decision

            obj.log.input(obj.report); % summary of observatory status

            if ~isempty(obj.sync) 
                
                if ~isempty(obj.mount.tracking) && obj.mount.tracking && obj.dome.is_closed==0
                    obj.sync.outgoing.stop_camera = 0; % if everything is cool, let the camera keep going
                end
                
                obj.sync.update;
                
            end
            
            if obj.use_shutdown && obj.devices_ok==0 % critical device failure, must shut down
                if obj.is_shutdown==0 % if already shut down, don't need to do it again
                    fprintf('%s: Device problems... %s \n', datestr(obj.log.time), obj.checker.report); 
                    obj.shutdown;
                end
            end

            if obj.use_shutdown && obj.sensors_ok==0 % critical device failure, must shut down
                if obj.is_shutdown==0 % if already shut down, don't need to do it again
                    fprintf('%s: Bad weather... %s \n', datestr(obj.log.time), obj.checker.report); 
                    obj.shutdown;
                end
            end
            
            if obj.use_shutdown && obj.checkDayTime % check if the system clock says it is day time
                if obj.is_shutdown==0 % if already shut down, don't need to do it again
                    fprintf('%s: Day time... %s \n', datestr(obj.log.time), obj.checker.report); 
                    obj.shutdown;
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
        
        function updateCameraComputer(obj)
            
            try 
                
                if ~obj.sync.is_connected % || ~obj.sync.status
                    obj.sync.connect;
                end

            catch 
                t = datetime('now', 'TimeZone', 'UTC'); 
                fprintf('%s: Failed to connect to camera computer\n', t); 
                % do nothing, as we can be waiting for ever for server to connect
            end
            
            % obj.sync.outgoing.OBJECT = ???
            if ~isempty(obj.mount) && obj.use_mount
                
                if isempty(obj.mount.sync)
                    obj.mount.sync = obj.sync; % share the handle to this object
                end
                
                obj.mount.updateCamera;
                
            end
            
            obj.sync.outgoing.TEMP_OUT = obj.average_temperature;
            obj.sync.outgoing.WIND_DIR = obj.average_wind_dir;
            obj.sync.outgoing.WIND_SPEED = obj.average_wind_speed;
            obj.sync.outgoing.HUMID_OUT = obj.average_humidity;
            obj.sync.outgoing.LIGHT = obj.average_light; 
            obj.sync.outgoing.PRESSURE = obj.average_pressure; 
            
            if obj.dome.is_closed
                
                if ~isfield(obj.sync.outgoing, 'stop_camera') || obj.sync.outgoing.stop_camera==0
                    obj.log.input('Dome closed, sending camera stop command');
                    disp(obj.log.report);
                end
                
                obj.sync.outgoing.stop_camera = 1;
            else
%                 obj.sync.outgoing.stop_camera = 0;
            end
            
            
            % add additional parameters and some commands like "start run"
            
            obj.sync.update;
            
        end
        
        function closeDome(obj)
            
            obj.dome.closeBothFull;
            
        end
        
        function shutdown(obj) % command to shut down observatory (close dome, stop tracking)
            
            obj.log.input('Shutting down observatory!');
            
            disp([char(obj.log.time) ': Shutting down observatory!']);
            
            try 

                obj.stop; % stop any slews or shutter motion

                obj.mount.tracking = 0; % later add command to park the telescope? 

                obj.closeDome;

                obj.sync.outgoing.stop_camera = 1; % make sure camera stops running also
                obj.sync.update;
                
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

