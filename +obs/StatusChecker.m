classdef StatusChecker < handle

    properties(Transient=true)
        
        t1; % quick timer (every minute or so) just to update sensors/devices
        t2; % check that everything is connected and that weather is good, print to log file
        t3; % verify that the other two are still running (every half an hour or so)
        
    end
    
    properties % objects
        
        owner; % loop back to manager class
        
        devices = {}; % all critical devices that must have status==1
        
        sensors = {}; % all the sensors we have... 
%         sensors_light = {};
%         sensors_clouds = {};
%         sensors_temp = {};
%         sensors_wind = {};
%         sensors_humid = {};
        
        status_log@util.sys.Logger;
        weather_log@util.sys.Logger;
        
    end
    
    properties % inputs/outputs
        
        light_now;
        light_str;
        light_all;
        light_jul;
        
        clouds_now;
        clouds_str;
        clouds_all;
        clouds_jul;
        
        temp_now;
        temp_str;
        temp_all;
        temp_jul;
        
        wind_now;
        wind_str;
        wind_all;
        wind_jul;
        
        wind_az_now;
        wind_az_str;
        wind_az_all;
        wind_az_jul;
        
        humid_now;
        humid_str;
        humid_all;
        humid_jul;
        
        dev_str;
        
        status = 1;
        report = 'OK';
        
    end
    
    properties % switches/controls
        
        period1 = 60
        period2 = 300; % time for equipment/weather check and log file
        period3 = 1800; % time for verifying shorter timers are working (and other tests?)
        
        % light and clouds may be just binary, so we can skip plotting and thresholding them...?
        max_light = 300; % units??
        min_light = -Inf;
        
        max_clouds = Inf; % degree difference to sky??
        min_clouds = -10; % negative --> cloudy??
        
        max_temp = 30;
        min_temp = 0;
        
        max_wind = 30;
        min_wind = -Inf;
        
        max_humid = 60;
        min_humid = -Inf;
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
        
        % list all the classes that status checker is following
        sensor_classes = {'obs.sens.Simulator', 'obs.sens.Boltwood'}; 
        device_classes = {'obs.dome.Simulator', 'obs.dome.AstroHaven', 'obs.mount.Simulator', 'obs.mount.ASA'};
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = StatusChecker(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.StatusChecker')
                if obj.debug_bit, fprintf('StatusChecker copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            
            elseif ~isempty(varargin) && isa(varargin{1}, 'obs.Manager')
                
                if obj.debug_bit, fprintf('StatusChecker constructor v%4.2f\n', obj.version); end
            
                obj.owner = varargin{1};
                
                obj.status_log = util.sys.Logger('Observatory_status');
                obj.weather_log = util.sys.Logger('Weather_report');
                
                obj.connect;
                
                obj.setup_t3;
                obj.setup_t2;
                obj.setup_t1;
                
            else
                
                error('Must supply a "Manager" object to StatusChecker constructor...');

            end
            
        end
        
        function connect(obj)
            
            if isempty(obj.owner)
                error('Must supply a top-level object as owner of StatusChecker!');
            end
            
            list = properties(obj.owner);
            
            obj.devices = {};
            obj.sensors = {};
            
            for ii = 1:length(list)
                
                name = list{ii};
                
                if isobject(obj.owner.(name))
                    
                    if any(strcmp(class(obj.owner.(name)), obj.device_classes))
                        obj.devices{end+1} = obj.owner.(name);
                    end
                    
                    if any(strcmp(class(obj.owner.(name)), obj.sensor_classes))
                        obj.sensors{end+1} = obj.owner.(name);
                    end
                    
                end
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.setup_t3;
            obj.setup_t2;
            obj.setup_t1;
            
            obj.reset_light;
            obj.reset_clouds;
            obj.reset_temp;
            obj.reset_wind;
            obj.reset_humid;
            
            obj.weather_log.reset;
            obj.status_log.reset;
            
        end
        
        function reset_light(obj)
            
            obj.light_all = [];
            obj.light_jul = [];
            
        end
        
        function reset_clouds(obj)
            
            obj.clouds_all = [];
            obj.clouds_jul = [];
            
        end
        
        function reset_temp(obj)
            
            obj.temp_all = [];
            obj.temp_jul = [];
            
        end
        
        function reset_wind(obj)
            
            obj.wind_all = [];
            obj.wind_jul = [];
            
        end
        
        function reset_wind_az(obj)
            
            obj.wind_az_all = [];
            obj.wind_az_jul = [];
            
        end
        
        function reset_humid(obj)
            
            obj.humid_all = [];
            obj.humid_all = [];
            
        end
        
    end
    
    methods % getters
        
        function val = getInstrID(~, other)
            
            if isprop(other, 'id') || ismethod(other, 'id')
                val = other.id;
            else
                c = split(class(other), '.');
                val = c{end};
            end
            
        end
        
        function val = is_running(obj)
            
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
        
        function callback_t1(obj, ~, ~)
            
            for ii = 1:length(obj.devices)
                if ismethod(obj.devices{ii}, 'update')
                    obj.devices{ii}.update;
                end
            end
            
            for ii = 1:length(obj.sensors)
                if ismethod(obj.sensors{ii}, 'update')
                    obj.sensors{ii}.update;
                end
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
        
        function callback_t2(obj, ~, ~)
            
            % make sure t1 is running! 
            if isempty(obj.t1) || ~isvalid(obj.t1) || strcmp(obj.t1.Running, 'off')
                obj.setup_t1;
            end

            str = '';
            obj.collect_light;
            str = [str 'LIGHT: ' obj.light_str];
            
            obj.collect_clouds;
            str = [str 'CLOUDS: ' obj.clouds_str];
            
            obj.collect_temp;
            str = [str 'TEMP: ' obj.temp_str];
            
            obj.collect_wind;
            str = [str 'WIND: ' obj.wind_str];
            
            obj.collect_humid;
            str = [str 'HUMID: ' obj.humid_str];
            
            obj.weather_log.input(str);
            
            obj.decision_all;
            obj.status_log.input(obj.report);
            
            if obj.status==0
                % do something to owner? close dome maybe?
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
            
            % make sure t2 is running! 
            if isempty(obj.t2) || ~isvalid(obj.t2) || strcmp(obj.t2.Running, 'off')
                obj.setup_t2;
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
    
    methods % collect data and make decisions
        
        function collect_light(obj)
            
            val = [];
            str = '';
            
            for ii = 1:length(obj.sensors)
                
                sens = obj.sensors{ii};
                new_val = [];
                
                if isprop(sens, 'light_value_average')
                    new_val = sens.light_value_average;
                elseif isprop(sens, 'light_value')
                    new_val = sens.light_value;
                elseif isprop(sens, 'light')
                    new_val = sens.light;
                elseif isprop(sens, 'Light')
                    new_val = sens.Light;
                elseif isprop(sens, 'daylight')
                    new_val = sens.daylight;
                end
                
                if sens.status==0
                    new_val = NaN;
                end
                
                if ~isempty(new_val)
                    
                    val = [val new_val];
                    
                    new_str = sprintf('%s=%4.2f ', obj.getInstrID(sens), new_val);
                    
                    str = [str new_str];
                    
                end
                
            end
            
            obj.light_now = val; 
            j = juliandate(datetime('now', 'timezone', 'UTC'));
            obj.light_str = str;
            
            if size(val,2)==size(obj.light_all,2)
                obj.light_all = [obj.light_all; val];
                obj.light_jul = [obj.light_jul; j]; 
            elseif isempty(val)
                % pass
            else
                obj.reset_light;
                obj.light_all = val;
                obj.light_jul = j;
            end
            
        end
        
        function collect_clouds(obj)
            
            val = [];
            str = '';
            
            for ii = 1:length(obj.sensors)
                
                sens = obj.sensors{ii};
                new_val = [];
                if isprop(sens, 'temp_sky_average')
                    new_val = sens.temp_sky_average;
                elseif isprop(sens, 'temp_sky')
                    new_val = sens.temp_sky;
                elseif isprop(sens, 'clouds')
                    new_val = sens.clouds;
                elseif isprop(sens, 'Clouds')
                    new_val = sens.Clouds;
                
                end
                
                if sens.status==0
                    new_val = NaN;
                end
                
                if ~isempty(new_val)
                    
                    val = [val new_val];
                    
                    new_str = sprintf('%s=%4.2f ', obj.getInstrID(sens), new_val);
                    
                    str = [str new_str];
                    
                end
                
            end
            
            obj.clouds_now = val; 
            j = juliandate(datetime('now', 'timezone', 'UTC'));
            obj.clouds_str = str;
            
            if size(val,2)==size(obj.clouds_all,2)
                obj.clouds_all = [obj.clouds_all; val];
                obj.clouds_jul = [obj.clouds_jul; j]; 
            elseif isempty(val)
                % pass
            else
                obj.reset_clouds;
                obj.clouds_all = val;
                obj.clouds_jul = j;
            end
            
        end
        
        function collect_temp(obj)
            
            val = [];
            str = '';
            
            for ii = 1:length(obj.sensors)
                
                sens = obj.sensors{ii};
                new_val = [];
                
                if isprop(sens, 'temperature_average')
                    new_val = sens.temperature_average;
                elseif isprop(sens, 'temperature')
                    new_val = sens.temperature;
                elseif isprop(sens, 'Temperature')
                    new_val = sens.Temperature;
                    
                elseif isprop(sens, 'temp')
                    new_val = sens.temp;
                elseif isprop(sens, 'Temp')
                    new_val = sens.Temp;
                end
                
                if sens.status==0
                    new_val = NaN;
                end
                
                if ~isempty(new_val)
                    
                    val = [val new_val];
                    
                    new_str = sprintf('%s=%4.2f ', obj.getInstrID(sens), new_val);
                    
                    str = [str new_str];
                    
                end
                
            end
            
            obj.temp_now = val; 
            j = juliandate(datetime('now', 'timezone', 'UTC'));
            obj.temp_str = str;
            
            if size(val,2)==size(obj.temp_all,2)
                obj.temp_all = [obj.temp_all; val];
                obj.temp_jul = [obj.temp_jul; j]; 
            elseif isempty(val)
                % pass
            else
                obj.reset_temp;
                obj.temp_all = val;
                obj.temp_jul = j;
            end
            
        end
        
        function collect_wind(obj)
            
            val = [];
            str = '';
            
            for ii = 1:length(obj.sensors)
                
                sens = obj.sensors{ii};
                new_val = [];
                if isprop(sens, 'wind_speed_average')
                    new_val = sens.wind_speed_average;
                elseif isprop(sens, 'wind_speed')
                    new_val = sens.wind_speed;
                elseif isprop(sens, 'WindSpeed')
                    new_val = sens.WindSpeed;
                elseif isprop(sens, 'wind')
                    new_val = sens.wind;
                elseif isprop(sens, 'Wind')
                    new_val = sens.Wind;
                end
                
                if sens.status==0
                    new_val = NaN;
                end
                
                if ~isempty(new_val)
                    
                    val = [val new_val];
                    
                    new_str = sprintf('%s=%4.2f ', obj.getInstrID(sens), new_val);
                    
                    str = [str new_str];
                    
                end
                
            end
            
            obj.wind_now = val; 
            j = juliandate(datetime('now', 'timezone', 'UTC'));
            obj.wind_str = str;
            
            if size(val,2)==size(obj.wind_all,2)
                obj.wind_all = [obj.wind_all; val];
                obj.wind_jul = [obj.wind_jul; j]; 
            elseif isempty(val)
                % pass
            else
                obj.reset_wind;
                obj.wind_all = val;
                obj.wind_jul = j;
            end
            
        end
        
        function collect_wind_az(obj)
            
            val = [];
            str = '';
            
            for ii = 1:length(obj.sensors)
                
                sens = obj.sensors{ii};
                new_val = [];
                
                if isprop(sens, 'wind_az_average')
                    new_val = sens.wind_az_average;
                elseif isprop(sens, 'wind_az')
                    new_val = sens.wind_az;
                elseif isprop(sens, 'WindAz')
                    new_val = sens.WindAz;
                elseif isprop(sens, 'wind_direction')
                    new_val = sens.wind_direction;
                elseif isprop(sens, 'WindDirection')
                    new_val = sens.WindDirection;
                end
                
                if sens.status==0
                    new_val = NaN;
                end
                
                if ~isempty(new_val)
                    
                    val = [val new_val];
                    
                    new_str = sprintf('%s=%4.2f ', obj.getInstrID(sens), new_val);
                    
                    str = [str new_str];
                    
                end
                
            end
            
            obj.wind_az_now = val; 
            j = juliandate(datetime('now', 'timezone', 'UTC'));
            obj.wind_az_str = str;
            
            if size(val,2)==size(obj.wind_az_all,2)
                obj.wind_az_all = [obj.wind_az_all; val];
                obj.wind_az_jul = [obj.wind_jul; j]; 
            elseif isempty(val)
                % pass
            else
                obj.reset_wind_az;
                obj.wind_az_all = val;
                obj.wind_az_jul = j;
            end
            
        end
        
        function collect_humid(obj)
            
            val = [];
            str = '';
            
            for ii = 1:length(obj.sensors)
                
                sens = obj.sensors{ii};
                new_val = [];
                
                if isprop(sens, 'humidity_average')
                    new_val = sens.humidity_average;
                elseif isprop(sens, 'humidity')
                    new_val = sens.humidity;
                elseif isprop(sens, 'Humidity')
                    new_val = sens.Humidity;
                elseif isprop(sens, 'humid')
                    new_val = sens.humid;
                elseif isprop(sens, 'Humid')
                    new_val = sens.Humid;
                
                end
                
                if sens.status==0
                    new_val = NaN;
                end
                
                if ~isempty(new_val)
                    
                    val = [val new_val];
                    
                    new_str = sprintf('%s=%4.2f ', obj.getInstrID(sens), new_val);
                    
                    str = [str new_str];
                    
                end
                
            end
            
            obj.humid_now = val; 
            j = juliandate(datetime('now', 'timezone', 'UTC'));
            obj.humid_str = str;
            
            if size(val,2)==size(obj.humid_all,2)
                obj.humid_all = [obj.humid_all; val];
                obj.humid_jul = [obj.humid_jul; j]; 
            elseif isempty(val)
                % pass
            else
                obj.reset_humid;
                obj.humid_all = val;
                obj.humid_jul = j;
            end
            
        end
        
        function status_all = check_devices(obj)
            
            obj.dev_str = '';
            
            status_all = [];
            
            for ii = 1:length(obj.devices)
                
                if isprop(obj.devices{ii}, 'status') 
                    status = obj.devices{ii}.status;
                elseif isprop(obj.devices{ii}, 'Status') 
                    status = obj.devices{ii}.Status;
                else
                    status = NaN;
                end
                
                if status==0
                    
                    if obj.debug_bit, fprintf('Device %s is malfunctioning. Attempting to reconnect...\n', class(obj.devices{ii})); end
                    
                    obj.logger.input(sprintf('Device %s is malfunctioning. Attempting to reconnect...\n', class(obj.devices{ii})));
                    
                    try
                        
                        if ismethod(obj.devices{ii}, 'disconnect')
                            obj.devices{ii}.disconnect;
                        end

                        pause(0.05);

                        if ismethod(obj.devices{ii}, 'connect')
                            obj.devices{ii}.connect;
                        end

                        pause(0.05);

                    catch ME
                        warning(ME.getReport);
                    end
                    
                    if isprop(obj.devices{ii}, 'status')
                        status = obj.devices{ii}.status;
                    elseif isprop(obj.devices{ii}, 'Status') 
                        status = obj.devices{ii}.Status;
                    else
                        status = NaN;
                    end

                end
                
                status_all = [status_all status];
                
                obj.dev_str = [obj.dev_str obj.getInstrID(obj.devices{ii}) ': ' num2str(status)];
                
            end
            
        end
        
        function val = decision_light(obj)
            
            if any(obj.light_now>obj.max_light) || any(obj.light_now<obj.min_light)
                val = 0;
            else
                val = 1;
            end
            
        end
        
        function val = decision_clouds(obj)
            
            if any(obj.clouds_now>obj.max_clouds) || any(obj.clouds_now<obj.min_clouds)
                val = 0;
            else
                val = 1;
            end
            
        end
        
        function val = decision_temp(obj)
            
            if any(obj.temp_now>obj.max_temp) || any(obj.temp_now<obj.min_temp)
                val = 0;
            else
                val = 1;
            end
            
        end
        
        function val = decision_wind(obj)
            
            if any(obj.wind_now>obj.max_wind) || any(obj.wind_now<obj.min_wind)
                val = 0;
            else
                val = 1;
            end
            
        end
        
        function val = decision_humid(obj)
            
            if any(obj.humid_now>obj.max_humid) || any(obj.humid_now<obj.min_humid)
                val = 0;
            else
                val = 1;
            end
            
        end
        
        function val = decision_devices(obj)
            
            status_all = obj.check_devices;
            
            if any(status_all==0)
                val = 0;
            else
                val = 1;
            end
            
        end
        
        function decision_all(obj)
            
            % should I add a rain-checker also??
            
            obj.status = 1;
            obj.report = 'OK';
            
            if obj.decision_devices==0
                obj.status = 0;
                obj.report = ['Critical device failure: ' obj.dev_str];
                return;
            end
            
            if obj.decision_light==0
                obj.status = 0;
                obj.report = ['Light too bright! ' obj.light_str];
                return;
            end
            
            if obj.decision_clouds==0
                obj.status = 0;
                obj.report = ['Sky is cloudy! ' obj.clouds_str];
                return;
            end
            
            if obj.decision_temp==0
                obj.status = 0;
                obj.report = ['Temperature out of range! ' obj.temp_str];
                return;
            end
            
            if obj.decision_wind==0
                obj.status = 0;
                obj.report = ['Wind too strong! ' obj.wind_str];
                return;
            end
            
            if obj.decision_humid==0
                obj.status = 0;
                obj.report = ['Humidity too high! ' obj.humid_str];
                return;
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

