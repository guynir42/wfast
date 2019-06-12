classdef SensorChecker < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        owner; % loop back to manager class
        
        sensors = {}; % all the sensors we have... 
        
        log@util.sys.Logger;
        
    end
    
    properties % inputs/outputs
        
        light_now;
        light_ids;
        light_str;
        light_all;
        light_jul;
        
        clouds_now;
        clouds_ids;
        clouds_str;
        clouds_all;
        clouds_jul;
        
        temp_now;
        temp_ids;
        temp_str;
        temp_all;
        temp_jul;
        
        wind_now;
        wind_ids;
        wind_str;
        wind_all;
        wind_jul;
        
        wind_az_now;
        wind_az_ids;
        wind_az_str;
        wind_az_all;
        wind_az_jul;
        
        humid_now;
        humid_ids;
        humid_str;
        humid_all;
        humid_jul;
        
        status = 1;
        sensors_ok = 1;
        report = 'OK';
        
    end
    
    properties % switches/controls
        
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
        
        show_day_frac = 0.2; % what fraction of a day to plot back
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
    end
    
    properties(Hidden=true)
        
        % list all the classes that status checker is following
        sensor_classes = {'obs.sens.Simulator', 'obs.sens.Boltwood', 'obs.sens.WindETH'}; 
        
        version = 1.00;
        
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
            
            obj.reset_light;
            obj.reset_clouds;
            obj.reset_temp;
            obj.reset_wind;
            obj.reset_humid;
            
            obj.log.reset;
            
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
        
        function val = getStatus(obj, device)
            
            if isprop(device, 'status') || ismethod(device, 'status')
                val = device.status;
            elseif isprop(device, 'status') || ismethod(device, 'status')
                val = device.Status;
            else
                val = NaN;
            end
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % decisions on sensor data
        
        function check_devices(obj) % to be depricated (device checks moved to manager)
            
            obj.dev_str = '';
            obj.dev_status = [];
            obj.dev_critical = [];
            
            for ii = 1:length(obj.dev_all)
                
                status_temp = obj.getStatus(obj.dev_all{ii});
                
                if status_temp==0 && ~obj.isIgnored(obj.dev_all{ii})
                    
                    if obj.debug_bit, fprintf('Device %s is malfunctioning. Attempting to reconnect...\n', class(obj.dev_all{ii})); end
                    
                    obj.status_log.input(sprintf('Device %s is malfunctioning. Attempting to reconnect...\n', class(obj.dev_all{ii})));
                    
                    try
                        
                        if ismethod(obj.dev_all{ii}, 'disconnect')
                            obj.dev_all{ii}.disconnect;
                        end

                        pause(0.05);

                        if ismethod(obj.dev_all{ii}, 'connect')
                            obj.dev_all{ii}.connect;
                        end

                        pause(0.05);

                    catch ME
                        warning(ME.getReport);
                    end
                    
                    status_temp = obj.getStatus(obj.dev_all{ii});

                end
                
                obj.dev_status(ii) = status_temp;
                obj.dev_critical(ii) = obj.isCritical(obj.dev_all{ii});
                
                obj.dev_str = [obj.dev_str ' ' obj.getInstrID(obj.dev_all{ii}) ': ' num2str(status_temp)];
                
            end
            
        end
        
    end
    
    methods % this is called on t1, only update sensors without checking results
        
        function update(obj)
            
            for ii = 1:length(obj.sensors)
                if ismethod(obj.sensors{ii}, 'update')
                    obj.sensors{ii}.update;
                end
            end
             
        end
        
    end
    
    methods % these are all called on t2
        
        function collect_light(obj)
            
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
                
                if sens.status==0
                    new_val = NaN;
                end
                
                if ~isempty(new_val)
                    
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
        
        function collect_clouds(obj)
            
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
                obj.clouds_ids{end+1} = obj.getInstrID(sens);
            end
            
        end
        
        function collect_temp(obj)
            
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
                
                if sens.status==0
                    new_val = NaN;
                end
                
                if ~isempty(new_val)
                    
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
        
        function collect_wind(obj)
            
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
                
                if sens.status==0
                    new_val = NaN;
                end
                
                if ~isempty(new_val)
                    
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
        
        function collect_wind_az(obj)
            
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
                
                if sens.status==0
                    new_val = NaN;
                end
                
                if ~isempty(new_val)
                    
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
        
        function collect_humid(obj)
            
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
                
                if sens.status==0
                    new_val = NaN;
                end
                
                if ~isempty(new_val)
                    
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
        
        function collect_all(obj)
            
            str = ''; % this is filled logged in weather-log
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
            
            obj.log.input(str);
            
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
        
        function decision_all(obj) % summary of all sensor info, make a decision if it is ok to continue
            
            % should I add a rain-checker also??
            obj.collect_all; % first get the data from the sensors
            
            obj.status = 1; % all devices are responding
            obj.sensors_ok = 1; % if sensors give bad results this turns to zero
            obj.report = 'OK'; % summary of conditions
            
            if obj.decision_light==0
                obj.sensors_ok = 0;
                obj.report = ['Light too bright! ' obj.light_str];
%                 obj.owner.dome.closeBothFull; % can we put this somewhere better??
                return;
            end
            
            if obj.decision_clouds==0
                obj.sensors_ok = 0;
                obj.report = ['Sky is cloudy! ' obj.clouds_str];
                return;
            end
            
            if obj.decision_temp==0
                obj.sensors_ok = 0;
                obj.report = ['Temperature out of range! ' obj.temp_str];
                return;
            end
            
            if obj.decision_wind==0
                obj.sensors_ok = 0;
                obj.report = ['Wind too strong! ' obj.wind_str];
                return;
            end
            
            if obj.decision_humid==0
                obj.sensors_ok = 0;
                obj.report = ['Humidity too high! ' obj.humid_str];
                return;
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function plotWeather(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('day_frac', obj.show_day_frac);
            input.scan_vars(varargin{:});
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            hold(input.ax, 'off');
            
            h_list = [];
            
            % get temperatures
            [v,t] = obj.getWeatherTimeData(obj.temp_all, obj.temp_jul, input.day_frac);
            h = plot(input.ax, t, v, '-');
            for ii = 1:length(h), h(ii).DisplayName = ['temperature (C) ' obj.temp_ids{ii}]; end 
            h_list = [h_list; h];
            hold(input.ax, 'on');
            
            % get wind data
            [v,t] = obj.getWeatherTimeData(obj.wind_all, obj.wind_jul, input.day_frac);
            h = plot(input.ax, t, v, '-o');
            for ii = 1:length(h), h(ii).DisplayName = ['wind (km/h) ' obj.wind_ids{ii}]; end
            h_list = [h_list; h];
            hold(input.ax, 'on');
            
            % get humidity data
            [v,t] = obj.getWeatherTimeData(obj.humid_all, obj.humid_jul, input.day_frac);
            h = plot(input.ax, t, v, '-+');
            for ii = 1:length(h), h(ii).DisplayName = ['Humidity (%) ' obj.humid_ids{ii}]; end
            h_list = [h_list; h];
            hold(input.ax, 'on');
            
            % get light data
            [v,t] = obj.getWeatherTimeData(obj.light_all, obj.light_jul, input.day_frac);
            h = plot(input.ax, t, v/20, '-*');
            for ii = 1:length(h), h(ii).DisplayName = ['daylight /20 ' obj.light_ids{ii}]; end
            h_list = [h_list; h];
            hold(input.ax, 'on');
            
            hold(input.ax, 'off');
            input.ax.YLim = [-10, 50];
            
            legend(input.ax, h_list, 'Location', 'SouthWest');
            
        end
        
        function [v,t] = getWeatherTimeData(obj, v, t, frac_day)
            
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

