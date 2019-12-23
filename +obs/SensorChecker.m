classdef SensorChecker < handle
% 
%
%
    properties(Transient=true)
        
    end
    
    properties % objects
        
        owner; % loop back to manager class
        
        sensors = {}; % all the sensors we have... 
        
        log@util.sys.Logger;
        
    end
    
    properties % inputs/outputs
        
        data_types = {'light', 'clouds', 'temperature', 'wind_speed', 'wind_dir', 'humidity', 'pressure', 'rain'}; 
        light = struct([]); 
        clouds  = struct([]); 
        temperature = struct([]); 
        wind_speed = struct([]); 
        wind_dir = struct([]); 
        humidity = struct([]); 
        pressure = struct([]); 
        rain = struct([]); 
        
        jd = [];
        
%         light_now;
%         light_ids;
%         light_str;
%         light_all;
%         light_jul;
%         
%         clouds_now;
%         clouds_ids;
%         clouds_str;
%         clouds_all;
%         clouds_jul;
%         
%         temp_now;
%         temp_ids;
%         temp_str;
%         temp_all;
%         temp_jul;
%         
%         wind_now;
%         wind_ids;
%         wind_str;
%         wind_all;
%         wind_jul;
%         
%         wind_az_now;
%         wind_az_ids;
%         wind_az_str;
%         wind_az_all;
%         wind_az_jul;
%         
%         humid_now;
%         humid_ids;
%         humid_str;
%         humid_all;
%         humid_jul;
        
        wise_data_raw;
        wise_data_struct;
        
        status = 1;
        sensors_ok = 1;
        report = 'OK';
        
    end
    
    properties % switches/controls
        
        use_wise_data = 1;
        
        show_day_frac = 0.2; % what fraction of a day to plot back on GUI
        
        use_only_plot_mean = 0; % only display the mean of each weather data type
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
    end
    
    properties(Hidden=true)
        
        % list all the classes that status checker is following
        sensor_classes = {'obs.sens.Simulator', 'obs.sens.Boltwood', 'obs.sens.WindETH', 'obs.sens.VirtualSensor'}; 
        
        version = 1.02;
        % 1.02 (2019/12/16) added virtual sensors, communications with Wise, and put all data into structures
        
    end
    
    methods % constructor
        
        function obj = SensorChecker(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.SensorChecker')
                if obj.debug_bit, fprintf('SensorChecker copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            
            elseif ~isempty(varargin) && isa(varargin{1}, 'obs.Manager')
                
                if obj.debug_bit, fprintf('SensorChecker constructor v%4.2f\n', obj.version); end
            
                obj.owner = varargin{1};
                
                obj.log = util.sys.Logger('Weather_report');
                
                obj.reset;
                
                obj.connect;
                
            else
                
                error('Must supply a "Manager" object to StatusChecker constructor...');

            end
            
        end
        
        function connect(obj)
            
            if isempty(obj.owner)
                error('Must supply a top-level object as owner of StatusChecker!');
            end
            
            obj.reset;
            
            list = properties(obj.owner);
            
            obj.sensors = {};
            
            for ii = 1:length(list)
                
                name = list{ii};
                
                if isobject(obj.owner.(name))
                    
                    if any(strcmp(class(obj.owner.(name)), obj.sensor_classes))
                        obj.sensors{end+1} = obj.owner.(name);
                    end
                    
                end
                
            end
            
            if obj.use_wise_data
                obj.connectVirtualSensors;
            end
            
        end
        
        function connectVirtualSensors(obj)
            
            obj.getWiseData;
            
            virtual_devices = {};
            
            % connect to Wise Boltwood 1 and 2
            for ii = 1:2
            
                device = obs.sens.VirtualSensor;
                device.data_path = {'Boltwood', ii, 'SensorData'};
                device.id = sprintf('BWW%d', ii); 
                virtual_devices{end+1} = device; 
                
            end
            
            % connect to Wise VantagePro
            device = obs.sens.VirtualSensor;
            device.data_path = {'VantagePro2', 'SensorData'};
            device.id = 'VPro'; 
            virtual_devices{end+1} = device; 
            
            % connect to Korean OWL/AWS
            device = obs.sens.VirtualSensor;
            device.data_path = {'OWL', 'AWS'};
            device.id = 'AWS'; 
            virtual_devices{end+1} = device; 

            % connect to Korean OWL/WDS X3
            for ii = 1:3
                
                device = obs.sens.VirtualSensor;
                device.data_path = {'OWL', 'WDS', ii};
                device.id = sprintf('WDS%d', ii); 
                virtual_devices{end+1} = device; 

            end
                        
            % connect to Korean OWL/THS X5
            for ii = 1:3
                
                device = obs.sens.VirtualSensor;
                device.data_path = {'OWL', 'THS', ii};
                device.id = sprintf('THS%d', ii); 
                virtual_devices{end+1} = device; 

            end     
            
            % connect to Korean OWL/CLS X2 (we are going to ignore CLS 3 as it is malfunctioning 
            for ii = 1:2
                
                device = obs.sens.VirtualSensor;
                device.data_path = {'OWL', 'CLS', ii};
                device.id = sprintf('CLS%d', ii); 
                virtual_devices{end+1} = device; 

            end
            
            for ii = 1:length(virtual_devices)
                virtual_devices{ii}.owner = obj;
                virtual_devices{ii}.connect;
            end
            
            obj.sensors = [obj.sensors virtual_devices]; 
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            for ii = 1:length(obj.data_types)
                obj.(['reset_', obj.data_types{ii}]); 
            end
            
%             obj.reset_light;
%             obj.reset_clouds;
%             obj.reset_temp;
%             obj.reset_wind_speed;
%             obj.reset_wind_dir;
%             obj.reset_humid;
%             obj.reset_pressure;
%             obj.reset_rain;
            
            obj.jd = [];
            
            obj.log.reset;
            
        end
        
        function reset_light(obj)
            
            obj.light = struct(); 
            obj.light.names = {'light', 'light_value', 'daylight'};
            obj.light.data = [];
            obj.light.sensors = {}; 
            obj.light.index = [];
            obj.light.now = [];
            obj.light.string = '';
            obj.light.jd = [];
            obj.light.err_str = 'Light is too bright!'; 
            obj.light.units = '';
            
            obj.light.max = 250; % units??
            
        end
        
        function reset_clouds(obj)
            
            obj.clouds = struct();
            obj.clouds.names = {'clouds', 'temp_sky', 'sky_temp', 'skyAmbientTemp', 'skyAmbientTemperature'};
            obj.clouds.data = [];
            obj.clouds.sensors = {};
            obj.clouds.index = [];
            obj.clouds.now = [];
            obj.clouds.string = '';
            obj.clouds.jd = [];
            obj.clouds.units = 'C';
            obj.clouds.err_str = 'Sky is cloudy!'; 
            
            obj.clouds.max = -10; % less negative --> more cloudy (clear days have -30 or -20, very cloudy is -10)
            
        end
        
        function reset_temperature(obj)

            obj.temperature = struct();
            obj.temperature.names = {'temperature', 'ambientTemp', 'ambientTemperature', 'outsideTemp', 'outsideTemperature'}; 
            obj.temperature.data = [];
            obj.temperature.sensors = {};
            obj.temperature.index = [];
            obj.temperature.now = [];
            obj.temperature.string = '';
            obj.temperature.jd = [];
            obj.temperature.units = 'C';
            obj.temperature.err_str = 'Temperature out of range!';
        
            obj.temperature.max = 30; % night time temperature
            obj.temperature.min = 0; % in C obviously
            
            
        end
        
        function reset_wind_speed(obj)

            obj.wind_speed = struct(); 
            obj.wind_speed.names = {'wind', 'wind_speed', 'windSpeed'};
            obj.wind_speed.data = [];
            obj.wind_speed.sensors = {};
            obj.wind_speed.index = [];
            obj.wind_speed.now = [];
            obj.wind_speed.string = '';
            obj.wind_speed.jd = [];
            obj.wind_speed.units = 'km/h';
            obj.wind_speed.err_str = 'Wind too strong!';
        
            obj.wind_speed.max = 50; % km/h
            
        end
        
        function reset_wind_dir(obj)

            obj.wind_dir = struct(); 
            obj.wind_dir.names = {'wind_dir', 'wind_az', 'windDir'};
            obj.wind_dir.data = [];
            obj.wind_dir.sensors = {};
            obj.wind_dir.index  = [];
            obj.wind_dir.now = [];
            obj.wind_dir.string = '';
            obj.wind_dir.jd = [];
            obj.wind_dir.units = 'deg';
            
        end
        
        function reset_humidity(obj)
            
            obj.humidity = struct();
            obj.humidity.names = {'humid', 'humidity', 'outsideHumidity'}; 
            obj.humidity.data = [];
            obj.humidity.sensors = {};
            obj.humidity.index = [];
            obj.humidity.now = [];
            obj.humidity.string = '';
            obj.humidity.jd = [];
            obj.humidity.units = 'percent';
            obj.humidity.err_str = 'Humidity too high!'; 
            
            obj.humidity.max = 85; % percent
            
        end
        
        function reset_pressure(obj)
            
            obj.pressure = struct();
            obj.pressure.names = {'pressure', 'airPressure', 'barometer'}; 
            obj.pressure.data = [];
            obj.pressure.sensors = {};
            obj.pressure.index = [];
            obj.pressure.now = [];
            obj.pressure.string = '';
            obj.pressure.jd = [];
            obj.pressure.units = 'mbar'; 
            obj.pressure.err_str = 'Pressure out of range!'; 
            
        end
        
        function reset_rain(obj)
            
            obj.rain = struct();
            obj.rain.names = {'rain', 'rainRate'}; 
            obj.rain.data = [];
            obj.rain.sensors = {};
            obj.rain.index = [];
            obj.rain.now = [];
            obj.rain.string = '';
            obj.rain.jd = [];
            obj.rain.good_if = false; % it is only OK to observe if this property is false! 
            obj.rain.err_str = 'It is raining!'; 
            
        end
        
    end
    
    methods % getters
        
        function val = getInstrID(~, device) % get a sensor and ask for its identifying string
            
            if isprop(device, 'id') || ismethod(device, 'id')
                val = device.id;
            else
                c = split(class(device), '.');
                val = c{end};
            end
            
        end
        
        function val = getStatus(~, device) % check if device is connected
            
            if isprop(device, 'status') || ismethod(device, 'status')
                val = device.status;
            elseif isprop(device, 'status') || ismethod(device, 'status')
                val = device.Status;
            else
                val = NaN;
            end
            
        end
        
        function val = getInstrTime(~, device)
            
            if isprop(device, 'jd') || ismethod(device, 'status')
                val = datetime(device.jd, 'ConvertFrom', 'juliandate'); 
            elseif isprop(device, 'time')
                val = device.time;
            end
                
        end
        
        function val = getInstrJD(~, device)
            
            if isprop(device, 'jd') || ismethod(device, 'jd')
                val = device.jd;
            elseif isprop(device, 'time')
                val = juliandate(device.time);
            elseif isprop(device, 'updatedAtUT')
                val = julianedate(util.text.str2time(device.updatedAtUT)); 
            end
                
        end
        
        function val = capital(~, name)
            
            name(1) = upper(name(1)); 
            
            val = name;
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % this is called on t1, only update sensors without checking results
        
        function update(obj)
            
            for ii = 1:length(obj.sensors)
                if ismethod(obj.sensors{ii}, 'update')
                    try
                        obj.sensors{ii}.update;
                    catch ME
                        obj.sensors{ii}.connect;
                        warning(ME.getReport);
                    end
                end
            end
            
        end
        
    end
    
    methods % these are all called on t2
           
        function collect(obj, type)
            
            if isprop(obj, type) || ~isa(obj.(type), 'struct')
            
                values = [];
                ids = {};
                str = '';
                idx = [];
%                 time_now = datetime('now', 'TimeZone', 'UTC'); 

                jd_now = juliandate(datetime('now', 'TimeZone', 'UTC'));
                
                for ii = 1:length(obj.sensors)

                    sensor = obj.sensors{ii}; % handle to the sensor
                    if isempty(sensor) % add extra tests here... don't test for status, as it can change (put NaN if bad status)
                        continue;
                    end
                    
                    new_val = []; % put the value we found into this
                    new_id = ''; % id of the sensor
                    new_str = '';
                    new_idx = [];
                    
                    for jj = 1:length(obj.(type).names)
                        
                        name = obj.(type).names{jj};
                        
                        if isprop(sensor, [name '_average']) % prefer average measurements
                            new_val = sensor.([name '_average']);                            
                        elseif isprop(sensor, [obj.capital(name) '_average']) % try with a capital first letter
                            new_val = sensor.([obj.capital(name) '_average']);
                        elseif isprop(sensor, name)
                            new_val = sensor.(name);
                        elseif isprop(sensor, obj.capital(name)) % try with a capital first letter
                            new_val = sensor.(obj.capital(name));
                        end
                        
                        if ~isempty(new_val)
                            
                            new_id = obj.getInstrID(sensor); 
                            new_idx = ii;
%                             new_time = obj.getInstrTime(sensor); 
                            new_jd = obj.getInstrJD(sensor);
                            
                            if ~obj.getStatus(sensor) || jd_now>new_jd + 15/60/24 % time_now>new_time - minutes(15) % if last update of sensor is more than 15 minutes ago it is not used
                                new_val = NaN;
                            end
                            
                            if ischar(new_val)
                                new_val = str2double(new_val);
                            end
                            
                            new_str = sprintf('%s= %4.2f', new_id, new_val); 
                            
                            break; % go to the next sensor! 
                            
                        end
                        
                    end % for jj (names)
            
                    if ~isempty(new_val) 
                        values(end+1) = new_val;
                        ids{end+1} = new_id;
                        idx(end+1) = new_idx; 
                        if isempty(str)
                            str = new_str;
                        else
                            str = sprintf('%s, %s', str, new_str); 
                        end
                    end
                    
                end % for ii (snesors)
                
                if ~isequal(obj.(type).sensors, ids) % added or removed or rearranged sensor list! 
                    obj.(['reset_', type]); % reset only the 
                end
                
                if ~isempty(values)
                    obj.(type).now = values;
                    obj.(type).sensors = ids; 
                    obj.(type).index = idx; 
                    obj.(type).jd(end+1) = jd_now;
                    obj.(type).data(end+1,:) = values;
                    obj.(type).string = str;
                end
                
            else
                error('Cannot find "%s" struct in SensorChecker!', type); 
            end
            
        end
              
        function collect_all(obj)
            
            if obj.use_wise_data
                obj.getWiseData; % first make sure all virtual sensors are updated
            end
            
            str = ''; % this is filled logged in weather-log
            
            for ii = 1:length(obj.data_types)
                
                name = obj.data_types{ii};
                
                obj.collect(name);
                
                str = [str upper(name) ': ' obj.(name).string ' | ']; 
                
            end
            
            obj.log.input(str);
            
        end
        
        function getWiseData(obj) % call the Wise computer and ask for all the weather data
            
            % I got this string from Arie, and had to install cURL for windows but now it works. Note the use of double quotes
            [rc,rv] = system('curl --connect-timeout 2 --silent -X PUT --header "Content-Type: application/x-www-form-urlencoded" --header "Accept: application/json" --data "Action=raw-weather-data&Parameters=" http://132.66.65.9:11111/api/v1/safetymonitor/0/action');
            
            if rc==0
                
                obj.wise_data_raw = rv;
                
                value = jsondecode(rv); 
                
                obj.wise_data_struct = jsondecode(value.Value); 
                
            end
            
        end
           
        function val = checkValueOK(obj, type) % check one data type if it is within bounds
            
            if isprop(obj, type) || ~isa(obj.(type), 'struct')
                
                st = obj.(type);
                value_now = nanmean(st.now); 
                
                if isempty(st.now)
                    val = 1;
                elseif isnan(value_now)
                    val = 0;
                elseif isfield(st, 'max') && value_now>st.max
                    val = 0;
                elseif isfield(st, 'min') && value_now<st.min
                    val = 0;
                else 
                    val = 1;
                end
                
            else
                error('Cannot find "%s" struct in SensorChecker!', type); 
            end
            
        end
        
        function val = checkBoolOK(obj, type) % check one data type with boolean values is not triggering
            
            if isprop(obj, type) || ~isa(obj.(type), 'struct')
                
                val = 1;
                
                st = obj.(type);
                if isfield(st, 'good_if')
                    failues = st.now == ~st.good_if; % locate the real failures, don't count NaNs
                    
                    if any(failues) || all(isnan(st.now)) % if any measurement fails (not on NaN) or all measurements give NaN
                        val = 0; 
                    end
                    
                end
                
            end
            
        end
        
        function decision_all(obj) % summary of all sensor info, make a decision if it is ok to continue
            
            % should I add a rain-checker also??
            obj.collect_all; % first get the data from the sensors
            
%             obj.status = 1; % all devices are responding
            obj.sensors_ok = 1; % if sensors give bad results this turns to zero
            obj.report = 'OK'; % summary of conditions
            
            for ii = 1:length(obj.data_types)
                
                name = obj.data_types{ii}; 
                
                if util.text.cs(name, 'wind_dir', 'pressure'), continue; end
                
                if ~obj.checkValueOK(name) || ~obj.checkBoolOK(name)
                    obj.sensors_ok = 0;
                    obj.report = obj.(name).err_str; 
                    return; 
                end
                
            end
            
            % check the wise general safety flag
            if obj.use_wise_data
                [rc,rv] = system('curl --connect-timeout 2 --silent -X PUT --header "Content-Type: application/x-www-form-urlencoded" --header "Accept: application/json" --data "Action=wise-issafe&Parameters=" http://132.66.65.9:11111/api/v1/safetymonitor/0/action');
                if(rc==0) % check the call succeeded
                    value = jsondecode(rv);
                    if strcmp(value.Value, 'False')
                        obj.sensors_ok = 0;
                        obj.report = 'Wise-unsafe! '; 
                        [rc,rv] = system('curl --connect-timeout 2 --silent -X PUT --header "Content-Type: application/x-www-form-urlencoded" --header "Accept: application/json" --data "Action=wise-unsafereasons&Parameters=" http://132.66.65.9:11111/api/v1/safetymonitor/0/action');
                        if(rc==0)
                            value = jsondecode(rv);
                            obj.report = [obj.report value.Value]; 
                        end
                    end
                end
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function plotWeather(obj, varargin)
            
            import util.text.cs;
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('day_frac', obj.show_day_frac);
            input.scan_vars(varargin{:});
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            input.ax.NextPlot = 'replace'; 
            
            ax_max = 50;
            ax_min = 0;
            
            markers = {'o', '+', '*', 's', 'v', 'x', 'd', 'p', '^', '>', '<', 'h'}; 
            Nmark = length(markers); 
            lines = {'-', '--', ':', '.-'}; 
            colors = {'yellow', 'cyan', 'red', 'magenta', 'black', 'blue', 'green', 'black'}; 
            
            for ii = 1:length(obj.data_types)
                
                name = obj.data_types{ii};
                
                if cs(name, 'rain', 'wind_dir') % don't plot these
                    continue;
                end
                
                [v,t] = obj.getWeatherTimeData(obj.(name).data, obj.(name).jd, input.day_frac);
                
                idx = obj.(name).index;
                unit = obj.(name).units;
                
                if cs(name, 'light') % scale to the axes
                    v = v/16;
                    unit = [unit '/16'];
                end
                
                if cs(name, 'clouds') % get the values on the same scale
                    v = v + 60; 
                    unit = ['+60 ' unit];
                end
                
                if cs(name, 'pressure') % scale the pressure by 20
                    v = v - 975;
                    unit = ['-975 ' unit];
                end
                
                if obj.use_only_plot_mean
                    h = plot(input.ax, t, nanmean(v,2), '-', 'LineWidth', 2, 'DisplayName', sprintf('mean %s [%s]', strrep(name, '_', ' '), unit), 'Color', colors{ii}); 
                    input.ax.NextPlot = 'add';
                else

                    h = plot(input.ax, t, v, '-');

                    for jj = 1:length(h) 
                        h(jj).DisplayName = sprintf('%s [%s]: %s', strrep(name, '_', ' '), unit, obj.(name).sensors{jj}); 
                        h(jj).Marker = markers{mod(idx(jj)-1, Nmark)+1};
                        h(jj).MarkerSize = 8;
                        h(jj).LineWidth = 0.5;
                        h(jj).LineStyle = lines{floor((idx(jj)-1)/Nmark)+1};
                        h(jj).Color = colors{ii}; 
                    end

                    input.ax.NextPlot = 'add';

                    plot(input.ax, t, nanmean(v,2), '-', 'LineWidth', 2, 'HandleVisibility', 'off', 'Color', h(1).Color); 
                
                end
                
                if max(v(:))>ax_max, ax_max = max(v(:))*1.1; end
                if min(v(:))<ax_min, ax_min = min(v(:))*1.1; end
                
            end
            
            input.ax.NextPlot = 'replace'; 
            
            input.ax.YLim = [ax_min, ax_max];
            
            hl = legend(input.ax, 'Location', 'SouthEastOutside');
            hl.FontSize = 12;
            
        end
        
        function [v,t] = getWeatherTimeData(~, v, t, frac_day)
            
            if isempty(t) || isempty(v)
                v = [];
                t = [];
                return;
            end
            
            idx = find(t>t(end)-frac_day, 1, 'first'); 
            if ~isempty(idx) && idx>0
                v = v(idx:end,:);
                t = t(idx:end);
            end
            
            t = datetime(t, 'ConvertFrom', 'juliandate');
            
        end
        
    end    
    
end

