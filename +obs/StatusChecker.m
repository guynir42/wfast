classdef StatusChecker < handle

    properties(Transient=true)
        
        short@timer;
        long@timer;
        
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
        
        wind_dir_now;
        wind_dir_str;
        wind_dir_all;
        wind_dir_jul;
        
        humid_now;
        humid_str;
        humid_all;
        humid_jul;
        
        dev_str;
        
        status = 1;
        report = 'OK';
        
    end
    
    properties % switches/controls
        
        period_short = 10; % time for equipment check
        period_long = 600; % time for verifying short timer is working (and other tests?)
        
        % light and clouds may be just binary, so we can skip plotting and thresholding them...?
        max_light = 1; % units??
        min_light = -Inf;
        
        max_clouds = 1; % units??
        min_clouds = -Inf; 
        
        max_temp = 30;
        min_temp = 0;
        
        max_wind = 20;
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
            else
                
                if obj.debug_bit, fprintf('StatusChecker constructor v%4.2f\n', obj.version); end
            
                obj.status_log = util.sys.Logger('Observatory_status');
                obj.weather_log = util.sys.Logger('Weather_report');
                
%                 obj.connect;
                
                obj.setup_long;
                obj.setup_short;

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
                    
                    if strcmp(class(obj.owner.(name)), obj.device_classes)
                        obj.devices{end+1} = obj.owner.(name);
                    end
                    
                    if strcmp(class(obj.owner.(name)), obj.sensor_classes)
                        obj.sensors{end+1} = obj.owner.(name);
                    end
                    
                end
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.setup_long;
            obj.setup_short;
            
            obj.reset_light;
            obj.reset_clouds;
            obj.reset_temp;
            obj.reset_wind;
            obj.reset_humid;
            
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
        
    end
    
    methods % setters
        
    end
    
    methods % timer related
        
        
        function callback_short(obj, ~, ~)
            
            disp('Running short timer now!');
            
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
        
        function setup_short(obj, ~, ~)
            
            if ~isempty(obj.short) && isvalid(obj.short)
                if strcmp(obj.short.Running, 'on')
                    obj.short.stop;
                end
            end
            
            obj.short = timer('BusyMode', 'queue', 'ExecutionMode', 'fixedRate', 'Name', 'Status-check-short', ...
                'Period', obj.period_short, 'StartDelay', obj.period_short, ...
                'TimerFcn', @obj.callback_short, 'ErrorFcn', @obj.setup_short);
            
            obj.short.start;
            
        end
        
        function val = callback_long(obj, ~, ~)
            
        end
        
        function setup_long(obj, ~, ~)
            
        end
        
    end
    
    methods % collect data and make decisions
        
        function collect_light(obj)
            
        end
        
        function collect_clouds(obj)
            
        end
        
        function collect_temp(obj)
            
            val = [];
            str = '';
            
            for ii = 1:length(obj.sensors)
                
                sens = obj.sensors{ii};
                new_val = [];
                
                if isprop(sens, 'temp')
                    new_val = sens.temp;
                elseif isprop(sens, 'Temp')
                    new_val = sens.Temp;
                elseif isprop(sens, 'temperature')
                    new_val = sens.temperature;
                elseif isprop(sens, 'Temperature')
                    new_val = sens.Temperature;
                end
                
                if sens.status==0
                    new_val = NaN;
                end
                
                if ~isempty(new_val)
                    
                    val = [val new_val];
                    
                    new_str = sprintf('%s: %4.2f ', obj.getInstrID(sens), new_val);
                    
                    str = [str new_str];
                    
                end
                
            end
            
            obj.temp_now = val; 
            j = juliandate(datetime('now', 'timezone', 'UTC'));
            
            if size(val,2)==size(obj.temp_all,2)
                obj.temp_all = [obj.temp_all; val];
                obj.temp_jul = [obj.temp_jul j]; 
            elseif isempty(val)
                % pass
            else
                obj.reset_temp;
                obj.temp_all = val;
                obj.temp_jul = j;
            end
            
        end
        
        function collect_wind(obj)
            
        end
        
        function collect_humid(obj)
            
        end
        
        function check_devices(obj)
            
        end
        
        function val = decision_light(obj)
            
            val = 1;
            
        end
        
        function val = decision_clouds(obj)
            
            val = 1;
            
        end
        
        function val = decision_temp(obj)
            
            if any(obj.temp_now>obj.max_temp) || any(obj.temp_now<obj.min_temp)
                val = 0;
            else
                val = 1;
            end
            
        end
        
        function val = decision_wind(obj)
            
            val = 1;
            
        end
        
        function val = decision_humid(obj)
            
            val = 1;
            
        end
        
        function val = decision_devices(obj)
            
            val = 1;
            
            for ii = 1:length(obj.devices)
                
                if isprop(obj.devices{ii}, 'status') 
                    val = obj.devices{ii}.status;
                elseif isprop(obj.devices{ii}, 'Status') 
                    val = obj.devices{ii}.Status;
                else
                    val = 1;
                end
                
                if val==0
                    return;
                end
                
            end
            
        end
        
        function decision_all(obj)
            
            % should I add a rain-checker also??
            
            obj.status = 1;
            
            if obj.decision_devices==0
                obj.status = 0;
                obj.report = ['Critical device failure: ' obj.dev_str];
                return;
            end
            
            if obj.decision_light==0
                obj.status = 0;
                obj.report = ['Light: ' obj.temp_str];
                return;
            end
            
            if obj.decision_clouds==0
                obj.status = 0;
                obj.report = ['Clouds: ' obj.temp_str];
                return;
            end
            
            if obj.decision_temp==0
                obj.status = 0;
                obj.report = ['Temperature: ' obj.temp_str];
                return;
            end
            
            if obj.decision_wind==0
                obj.status = 0;
                obj.report = ['Wind: ' obj.temp_str];
                return;
            end
            
            if obj.decision_humid==0
                obj.status = 0;
                obj.report = ['Humidity: ' obj.temp_str];
                return;
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

