classdef (CaseInsensitiveProperties, TruncatedProperties) Manager < handle

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
        
    end
    
    properties % switches/controls
        
        % use these to override these devices/sensors
        use_dome = 1;
        use_mount = 1;
        use_weather = 1;
        use_wind = 1;
        use_humidity = 0;
        use_temperature = 0; 
        
        period1 = 60; % time between updates of all devices/sensors
        period2 = 300; % time for equipment/weather check and log file
        period3 = 1800; % time for verifying shorter timers are working (and other tests?)
        
        brake_bit = 1;
        debug_bit = 1;
        
    end
    
    properties % inputs/outputs
                
        devices_ok = 1;
        devices_report = 'OK';
        
    end
    
    properties(Dependent=true)
        
        sensors_ok;
        sensors_report;
        report; 
        
        RA;
        DEC;
        LST;
        ALT;
        
        is_shutdown;
        
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
                
                obj.connect; % connect to all hardware
                
            end
            
        end
        
        function connect(obj)
            
            if obj.use_dome
                obj.connectDome;
            end
            
            if obj.use_mount
                obj.connectMount;
            end
            
            if obj.use_weather
                obj.connectBoltwood;
            end
            
            if obj.use_wind
                obj.connectWindETH;
            end
            
            % add additional devices
            % ...
            
            
            obj.setup_t3;
            obj.setup_t2;
            obj.setup_t1;
                
            
        end
        
        function delete(obj) % destructor
            
            delete(obj.t3);
            delete(obj.t2);
            delete(obj.t1);
            
        end
        
        function connectDome(obj)
            
            obj.log.input('Connecting to dome.');
            
            try 
                obj.dome = obs.dome.AstroHaven;
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
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.setup_t3;
            obj.setup_t2;
            obj.setup_t1;
            
        end
        
    end
    
    methods % getters
        
        function val = get.sensors_ok(obj)
            
            if isempty(obj.checker)
                val = [];
            else
                val = obj.checker.sensors_ok;
            end
            
        end
        
        function val = get.sensors_report(obj)
            
            if isempty(obj.checker)
                val = '';
            else
                val = obj.checker.report;
            end
            
        end
        
        function val = get.report(obj)
            
            val = sprintf('Sensors: %s | Devices: %s | state: %s', obj.sensors_report, obj.devices_report, obj.observatory_state);
            
        end
        
        function val = get.is_shutdown(obj)
            
            val = 0;
            
            try 
            
                if obj.mount.tracking
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
        
        function val = observatory_state(obj)
            
            if obj.is_shutdown
                val = 'SHUT';
            else
                val = 'OPEN';
            end
            
        end
        
        function val = get.RA(obj)
            
            if ~isempty(obj.mount)
                val = obj.mount.telRA;
            else
                val = [];
            end
            
        end
        
        function val = get.DEC(obj)
            
            if ~isempty(obj.mount)
                val = obj.mount.telDEC;
            else
                val = [];
            end
            
        end
        
        function val = get.LST(obj)
            
            if ~isempty(obj.mount)
                val = obj.mount.LST;
            else
                val = [];
            end
            
        end
        
        function val = get.ALT(obj)
            
            if ~isempty(obj.mount)
                val = round(obj.mount.telALT);
            else
                val = [];
            end
            
        end
        
        function val = average_temp(obj)
            
            val = mean(obj.checker.temp_now, 'omitnan');
            
        end
        
        function val = average_clouds(obj)
            
            val = mean(obj.checker.clouds_now, 'omitnan');
            
        end
        
        function val = average_light(obj)
            
            val = mean(obj.checker.light_now, 'omitnan');
            
        end
        
        function val = average_wind(obj)
            
            val = mean(obj.checker.wind_now, 'omitnan');
            
        end
        
        function val = average_wind_az(obj)
            
            val = mean(obj.checker.wind_az_now, 'omitnan');
            
        end
        
        function val = average_humid(obj)
            
            val = mean(obj.checker.humid_now, 'omitnan');
            
        end
        
        function val = areTimersRunning(obj)
            
            val(1) = strcmp(obj.t1.Running, 'on');
            val(2) = strcmp(obj.t2.Running, 'on');
            val(3) = strcmp(obj.t3.Running, 'on');
            
        end
        
    end
    
    methods % setters
        
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
        
        function callback_t1(obj, ~, ~) % update sensors 
            
            try 
            
                if isempty(obj.checker)
                    obj.connectSensorChecker;
                end
                
                obj.checker.update; % go over all sensors and only tell them to collect data. It's reported back in t2
                
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

                obj.update;
                
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
        
        function callback_t3(obj, ~, ~)
            
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
        
        function updateDevices(obj)
        
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
            
    end
    
    methods % calculations / commands
        
        function update(obj)

            obj.updateDevices;

            obj.checker.decision_all; % collect weather data and make a decision

            obj.log.input(obj.report);

            if obj.devices_ok==0
                if obj.is_shutdown==0
                    obj.shutdown;
                end
            end

            if obj.sensors_ok==0
                if obj.is_shutdown==0
                    obj.shutdown;
                end
            end
            
        end
        
        function shutdown(obj)
            
            obj.log.input('Shutting down observatory!');
            
            try 

                obj.stop; % stop any slews or shutter motion

                obj.mount.tracking = 0; % later add command to park the telescope? 

                obj.dome.closeBothFull;

                % anything else we can do to put the dome to shutdown mode?
                % ...

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function stop(obj)
            
            obj.brake_bit = 1;
            
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

