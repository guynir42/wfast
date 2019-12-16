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
        
%         max_light = 250; % units??
%         min_light = -Inf;
        
%         max_clouds = -15; % negative --> cloudy (clear days have -30 or -20, very cloudy is -10)
%         min_clouds = -Inf; % degree difference to sky
        
%         max_temp = 30; % night time temperature
%         min_temp = 0; % in C obviously
        
%         max_wind = 50; % km/h
%         min_wind = -Inf;
        
%         max_humid = 85; % percent
%         min_humid = -Inf;
        
%         max_pressure = Inf; % millibar
%         min_pressure = -Inf; 
        
        show_day_frac = 0.2; % what fraction of a day to plot back on GUI
        
        use_wise_data = 1;
        
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
            obj.light.now = [];
            obj.light.string = '';
            obj.light.jd = [];
            obj.light.err_str = 'Light is too bright!'; 
            obj.light.units = '';
            
            obj.light.max = 250; % units??
            
        end
        
        function reset_clouds(obj)
            
            obj.clouds = struct();
            obj.clouds.names = {'clouds', 'temp_sky', 'sky_temp', 'skyAmbientTemperature'};
            obj.clouds.data = [];
            obj.clouds.sensors = {};
            obj.clouds.now = [];
            obj.clouds.string = '';
            obj.clouds.jd = [];
            obj.clouds.units = 'C';
            obj.clouds.err_str = 'Sky is cloudy!'; 
            
            obj.clouds.max = -15; % negative --> cloudy (clear days have -30 or -20, very cloudy is -10)
            
        end
        
        function reset_temperature(obj)

            obj.temperature = struct();
            obj.temperature.names = {'temperature', 'ambientTemperature'}; 
            obj.temperature.data = [];
            obj.temperature.sensors = {};
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
            obj.wind_speed.now = [];
            obj.wind_speed.string = '';
            obj.wind_speed.jd = [];
            obj.wind_speed.units = 'km/h';
            obj.wind_speed.err_str = 'Wind to strong!';
        
            obj.wind_speed.max = 50; % km/h
            
        end
        
        function reset_wind_dir(obj)

            obj.wind_dir = struct(); 
            obj.wind_dir.names = {'wind_dir', 'wind_az', 'windDir'};
            obj.wind_dir.data = [];
            obj.wind_dir.sensors = {};
            obj.wind_dir.now = [];
            obj.wind_dir.string = '';
            obj.wind_dir.jd = [];
            obj.wind_dir.units = 'deg';
            
        end
        
        function reset_humidity(obj)
            
            obj.humidity = struct();
            obj.humidity.names = {'humid', 'humidity'}; 
            obj.humidity.data = [];
            obj.humidity.sensors = {};
            obj.humidity.now = [];
            obj.humidity.string = '';
            obj.humidity.jd = [];
            obj.humidity.units = 'percent';
            obj.humidity.err_str = 'Humidity too high!'; 
            
            obj.humidity.max = 85; % percent
            
        end
        
        function reset_pressure(obj)
            
            obj.pressure = struct();
            obj.pressure.names = {'pressure', 'airPressure'}; 
            obj.pressure.data = [];
            obj.pressure.sensors = {};
            obj.pressure.now = [];
            obj.pressure.string = '';
            obj.pressure.jd = [];
            obj.pressure.units = 'mbar'; 
            obj.pressure.err_str = 'Pressure out of range!'; 
            
        end
        
        function reset_rain(obj)
            
            obj.rain = struct();
            obj.rain.names = {'rain'}; 
            obj.rain.data = [];
            obj.rain.sensors = {};
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
            
            if isprop(device, 'jd') || ismethod(device, 'status')
                val = device.jd;
            elseif isprop(device, 'time')
                val = juliandate(device.time);
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
                        
                        if ~isempty(new_val) % we found a sensor with the correct
                            
                            new_id = obj.getInstrID(sensor); 
%                             new_time = obj.getInstrTime(sensor); 
                            new_jd = obj.getInstrJD(sensor);
                            
                            if ~obj.getStatus(sensor) || jd_now>new_jd + 15/60/24 % time_now>new_time - minutes(15) % if last update of sensor is more than 15 minutes ago it is not used
                                new_val = NaN;
                            end
                            
                            new_str = sprintf('%s= %4.2f', new_id, new_val); 
                            
                            break; % go to the next sensor! 
                            
                        end
                                                
                    end % for jj (names)
            
                    if ~isempty(new_val) 
                        values(end+1) = new_val;
                        ids{end+1} = new_id;
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
                    obj.(type).jd(end+1) = jd_now;
                    obj.(type).data(end+1,:) = values;
                end
                
            else
                error('Cannot find "%s" struct in SensorChecker!', type); 
            end
            
        end
        
        function collect_light(obj) % to be depricated
            
            val = [];
            str = '';
            obj.light_ids = {};
            
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
                
                if ~isempty(new_val)
                    
                    if sens.status==0
                        new_val = NaN;
                    end
                
                    val = [val new_val];
                    
                    new_str = sprintf('%s=%4.2f ', obj.getInstrID(sens), new_val);
                    
                    str = [str new_str];
                    obj.light_ids{end+1} = obj.getInstrID(sens);
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
        
        function collect_clouds(obj) % to be depricated
            
            val = [];
            str = '';
            obj.clouds_ids = {};
            
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
                
                
                if ~isempty(new_val)
                    
                    if sens.status==0
                        new_val = NaN;
                    end
                
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
                obj.clouds_ids{end+1} = obj.getInstrID(sens);
            end
            
        end
        
        function collect_temp(obj) % to be depricated
            
            val = [];
            str = '';
            obj.temp_ids = {};
            
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
                
                
                
                if ~isempty(new_val)
                    
                    if sens.status==0
                        new_val = NaN;
                    end
                    
                    val = [val new_val];
                    
                    new_str = sprintf('%s=%4.2f ', obj.getInstrID(sens), new_val);
                    
                    str = [str new_str];
                    
                    obj.temp_ids{end+1} = obj.getInstrID(sens);
                    
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
        
        function collect_wind(obj) % to be depricated
            
            val = [];
            str = '';
            obj.wind_ids = {};
            
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
                
                if ~isempty(new_val)
                    
                    if sens.status==0
                        new_val = NaN;
                    end
                
                    val = [val new_val];
                    
                    new_str = sprintf('%s=%4.2f ', obj.getInstrID(sens), new_val);
                    
                    str = [str new_str];
                    
                    obj.wind_ids{end+1} = obj.getInstrID(sens);
                    
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
        
        function collect_wind_az(obj) % to be depricated
            
            val = [];
            str = '';
            obj.wind_az_ids = {};
            
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
                
                
                
                if ~isempty(new_val)
                    
                    if sens.status==0
                        new_val = NaN;
                    end
                    
                    val = [val new_val];
                    
                    new_str = sprintf('%s=%4.2f ', obj.getInstrID(sens), new_val);
                    
                    str = [str new_str];
                    
                    obj.wind_az_ids{end+1} = obj.getInstrID(sens);
                
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
        
        function collect_humid(obj) % to be depricated
            
            val = [];
            str = '';
            obj.humid_ids = {};
            
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
                
                if ~isempty(new_val)
                    
                    if sens.status==0
                        new_val = NaN;
                    end
                
                    val = [val new_val];
                    
                    new_str = sprintf('%s=%4.2f ', obj.getInstrID(sens), new_val);
                    
                    str = [str new_str];
                    
                    obj.humid_ids{end+1} = obj.getInstrID(sens);
                    
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
        
        function getWiseData(obj) % call the Wise computer and ask for all the weather data
            
            % I got this string from Arie, and had to install cURL for windows but now it works. Note the use of double quotes
            [rc,rv] = system('curl --connect-timeout 2 --silent -X PUT --header "Content-Type: application/x-www-form-urlencoded" --header "Accept: application/json" --data "Action=raw-weather-data&Parameters=" http://132.66.65.9:11111/api/v1/safetymonitor/0/action');
            
            if rc==0
                
                obj.wise_data_raw = rv;
                
                value = jsondecode(rv); 
                
                obj.wise_data_struct = jsondecode(value.Value); 
                
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
                
                str = [str upper(name) ': ' obj.(name).string]; 
                
            end
            
%             obj.collect('light');
%             str = [str 'LIGHT: ' obj.light.string];
%             
%             obj.collect('clouds');
%             str = [str 'CLOUDS: ' obj.clouds.string];
%             
%             obj.collect('temperature');
%             str = [str 'TEMP: ' obj.temperature.string];
%             
%             obj.collect('wind_speed');
%             str = [str 'WIND: ' obj.wind_speed.string];
%             
%             obj.collect('wind dir');
%             str = [str 'WIND_AZ: ' obj.wind_dir.string];
%             
%             obj.collect('humidity');
%             str = [str 'HUMID: ' obj.humidity.string];
            
            obj.log.input(str);
            
        end
        
        function val = checkValueOK(obj, type) % check one data type if it is within bounds
            
            if isprop(obj, type) || ~isa(obj.(type), 'struct')
                
                val = 1;
                
                st = obj.(type);
                value_now = nanmean(st.now); 
                
                if isnan(value_now)
                    val = 0;
                end
                
                if isfield(st, 'max') && value_now>st.max
                    val = 0;
                end
                
                if isfield(st, 'min') && value_now<st.min
                    val = 0;
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
        
        function val = decision_light(obj) % to be depricated
            
            if any(obj.light_now>obj.max_light) || any(obj.light_now<obj.min_light)
                val = 0;
            else
                val = 1;
            end
            
        end
        
        function val = decision_clouds(obj) % to be depricated
            
            if any(obj.clouds_now>obj.max_clouds) || any(obj.clouds_now<obj.min_clouds)
                val = 0;
            else
                val = 1;
            end
            
        end
        
        function val = decision_temp(obj) % to be depricated
            
            if any(obj.temp_now>obj.max_temp) || any(obj.temp_now<obj.min_temp)
                val = 0;
            else
                val = 1;
            end
            
        end
        
        function val = decision_wind(obj) % to be depricated
            
            if any(obj.wind_now>obj.max_wind) || any(obj.wind_now<obj.min_wind)
                val = 0;
            else
                val = 1;
            end
            
        end
        
        function val = decision_humid(obj) % to be depricated
            
            if any(obj.humid_now>obj.max_humid) || any(obj.humid_now<obj.min_humid)
                val = 0;
            else
                val = 1;
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
                
                if ~obj.checkValueOK(name) || ~obj.checkBoolOK(name)
                    obj.sensors_ok = 0;
                    obj.report = obj.(name).err_str; 
                end
                
            end
            
%             if obj.decision_light==0
%                 obj.sensors_ok = 0;
%                 obj.report = ['Light too bright! ' obj.light_str];
% %                 obj.owner.dome.closeBothFull; % can we put this somewhere better??
%                 return;
%             end
%             
%             if obj.decision_clouds==0
%                 obj.sensors_ok = 0;
%                 obj.report = ['Sky is cloudy! ' obj.clouds_str];
%                 return;
%             end
%             
%             if obj.decision_temp==0
%                 obj.sensors_ok = 0;
%                 obj.report = ['Temperature out of range! ' obj.temp_str];
%                 return;
%             end
%             
%             if obj.decision_wind==0
%                 obj.sensors_ok = 0;
%                 obj.report = ['Wind too strong! ' obj.wind_str];
%                 return;
%             end
%             
%             if obj.decision_humid==0
%                 obj.sensors_ok = 0;
%                 obj.report = ['Humidity too high! ' obj.humid_str];
%                 return;
%             end
            
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
            
            h_list = [];
            ax_max = 50;
            markers = {'o', '+', '*', 's', 'v', '^','x', 'd', 'p'}; 
            
            for ii = 1:length(obj.data_types)
                
                name = obj.data_types{ii};
                
                if cs(name, 'rain', 'wind_dir') % don't plot these
                    continue;
                end
                
                [v,t] = obj.getWeatherTimeData(obj.(name).data, obj.(name).jd, input.day_frac);
                
                if cs(name, 'light') % scale to the axes
                    v = v/20;
                end
                
                if cs(name, 'clouds') % flip the values
                    v = -v; 
                end
                
                input.ax.ColorOrderIndex = 1;
                
                h = plot(input.ax, t, v, '-', 'Marker', markers{ii});
                
                for jj = 1:length(h) 
                    h(jj).DisplayName = sprintf('%s [%s]: %s', strrep(name, '_', ' '), obj.(name).units, obj.(name).sensors{jj}); 
                end
                
                input.ax.NextPlot = 'add';
                
                h_list = [h_list; h];
                
                if max(v)>ax_max, ax_max = max(v); end
                
            end
            
%             
%             % get temperatures
%             [v,t] = obj.getWeatherTimeData(obj.temp_all, obj.temp_jul, input.day_frac);
%             h = plot(input.ax, t, v, '-');
%             for ii = 1:length(h), h(ii).DisplayName = ['temperature (C) ' obj.temp_ids{ii}]; end 
%             h_list = [h_list; h];
%             if max(v)>ax_max, ax_max = max(v); end
%             hold(input.ax, 'on');
%             
%             % get wind data
%             [v,t] = obj.getWeatherTimeData(obj.wind_all, obj.wind_jul, input.day_frac);
%             h = plot(input.ax, t, v, '-o');
%             for ii = 1:length(h), h(ii).DisplayName = ['wind (km/h) ' obj.wind_ids{ii}]; end
%             h_list = [h_list; h];
%             if max(v)>ax_max, ax_max = max(v); end
%             hold(input.ax, 'on');
%             
%             % get humidity data
%             [v,t] = obj.getWeatherTimeData(obj.humid_all, obj.humid_jul, input.day_frac);
%             h = plot(input.ax, t, v, '-+');
%             for ii = 1:length(h), h(ii).DisplayName = ['Humidity (%) ' obj.humid_ids{ii}]; end
%             h_list = [h_list; h];
%             if max(v)>ax_max, ax_max = max(v); end
%             hold(input.ax, 'on');
%             
%             % get light data
%             [v,t] = obj.getWeatherTimeData(obj.light_all, obj.light_jul, input.day_frac);
%             v = v/20;
%             h = plot(input.ax, t, v, '-*');
%             for ii = 1:length(h), h(ii).DisplayName = ['daylight /20 ' obj.light_ids{ii}]; end
%             h_list = [h_list; h];
%             if max(v)>ax_max, ax_max = max(v); end
%             hold(input.ax, 'on');
            
            input.ax.NextPlot = 'replace'; 
            
            input.ax.YLim = [-10, ax_max];
            
            legend(input.ax, h_list, 'Location', 'SouthWest');
            
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

