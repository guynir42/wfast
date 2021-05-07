% Boltwood class  
% Package: +obs/+sens/
% Description: A class for controlling the Boltwood meteorological sensor
%              (from http://diffractionlimited.com/product/boltwood-cloud-sensor-ii/).
% Tested : Matlab R2018a
%     By : David Polishook                    Dec 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: % Create a Boltwood object
%          WB=actxserver('ClarityII.CloudSensorII');
%          % update weather station values in WB object properties
%          update(WB) - update all parameters from theweather program
%          install ASCOM;
%          
% Reliable: 2
%--------------------------------------------------------------------------

classdef Boltwood < handle
    
    
    properties % objects
        
        hndl; % ascom object loaded using actxserver 
        
    end
    
    properties
    
        % generic fields
        status = 0; % false - readings are unreliable, true - ok
        use_this = 1; % is true get the data from this sensor. If not, ignore it and move on
        id = 'BWW'; 
        data@util.vec.CircularBuffer; % Stack object containing data history
        data_col = {'JD', 'DayLightV', 'AmbientT', 'RelSkyT', 'SensorT', 'DewPointT', 'WindSpeed', 'Humidity', 'Rain'}; % Stack object columns
        data_units = {'day', 'relative', 'deg C', 'deg C', 'deg C', 'deg C', 'km/h', 'percentage', 'mm'}; % Stack object column units

        % specific fields
        jd = NaN; % latest date when data was successfully updated
        light = NaN; % Day light value
        temperature = NaN; % Ambient temperature
        temp_sky = NaN; % Sky minus ambient temperature; an indicator for the clouds condition
        temp_sensor = NaN; % Sensor temperature
        temp_dew_point = NaN; % Dew point temperature
        wind_speed = NaN; % Wind speed
        humidity = NaN; % Humiditiy
        rain = NaN; % Rain (logical 0/1)
        
        % Text conditions from Boltwood
        light_condition   = '';                        % Day condition
        cloud_condition = '';                        % Cloud condition
        wind_condition  = '';                        % Wind condition
        rain_condition  = '';                        % Rain condition 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        light_average;
        temperature_average;
        temp_sky_average;
        temp_sensor_average;
        temp_dew_point_average;
        wind_speed_average;
        humidity_average;
        rain_average;
        
    end
    
    properties(Hidden=true)
        
        version = 1.01;
        
    end
    
    methods % Constructor
        
        function obj = Boltwood(varargin)
            
            % Boltwood class constractor
            % Example: obj = obs.sens.Boltwood
            
            if isempty(varargin)
            
                if obj.debug_bit>1, fprintf('Boltwood default constructor v%4.2f\n', obj.version); end
                
                obj.data = util.vec.CircularBuffer;
                obj.connect; % Open Boltwood weather application. Also open connection to application if not already opened.
                
            end
            
        end
        
        function connect(obj) % Open Boltwood weather application. Also open connection to application if not already opened.
            
            [ret, str] = system('taskkill /fi "imagename eq Clarity.exe" /f'); % force killing of the server.

            obj.hndl = actxserver('ClarityII.CloudSensorII');
            
            obj.update;
            
        end

        function disconnect(obj)
            
            [ret, str] = system('taskkill /fi "imagename eq ClarityII.exe" /f'); % force killing of the server.

            delete(obj.hndl);
            obj.hndl = [];
            
        end
        
        function delete(obj)
            
            obj.disconnect;
            
        end
        
    end
    
    methods % getters
        
        function val = get.light_average(obj)
            
            if isempty(obj.data)
                val = [];
            else
                val = obj.data.median;
                if length(val)>=2, val = val(2); end
            end
            
        end
        
        function val = get.temperature_average(obj)
            
            if isempty(obj.data)
                val = [];
            else
                val = obj.data.median;
                if length(val)>=3, val = val(3); end
            end
            
        end
        
        function val = get.temp_sky_average(obj)
            
            if isempty(obj.data)
                val = [];
            else
                val = obj.data.median;
                if length(val)>=4, val = val(4); end
            end
            
        end
        
        function val = get.temp_sensor_average(obj)
            
            if isempty(obj.data)
                val = [];
            else
                val = obj.data.median;
                if length(val)>=5, val = val(5); end
            end
            
        end
        
        function val = get.temp_dew_point_average(obj)
            
            if isempty(obj.data)
                val = [];
            else
                val = obj.data.median;
                val = val(6);
            end
            
        end
        
        function val = get.wind_speed_average(obj)
            
            if isempty(obj.data)
                val = [];
            else
                val = obj.data.median;
                if length(val)>=7, val = val(7); end
            end
            
        end
        
        function val = get.humidity_average(obj)
            
            if isempty(obj.data)
                val = [];
            else
                val = obj.data.median;
                if length(val)>=8, val = val(8); end
            end
            
        end
        
        function val = get.rain_average(obj)
            % this is logical values!!
            
            if isempty(obj.data)
                val = [];
            else
                val = obj.data.median;
                if length(val)>=9, val = val(9); end
            end
            
        end
        
    end
    
    methods
        
        function reset(obj)
            
            obj.data.reset;
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            obj.status = 0;
            obj.jd = NaN;
            obj.light = NaN;
            obj.temperature = NaN;
            obj.temp_sky = NaN;
            obj.temp_sensor = NaN;
            obj.temp_dew_point = NaN;
            obj.wind_speed = NaN;
            obj.humidity = NaN;
            obj.rain = NaN;

            obj.light_condition   = '';
            obj.cloud_condition = '';
            obj.wind_condition  = '';
            obj.rain_condition  = '';
            
        end
        
        function update(obj) % Read weather parameters
            
            obj.status = 0;
            
            % if status is ok than set the time of the last query
            for ii = 1:100
                if obj.hndl.DataReady
                    obj.status = 1;
                    break;
                end
                
                pause(0.05);
                
            end
            
            if obj.status==0
                return;
            end
            
            % If the weather application is down then reinitiate it
            try
                
                for ii = 1:100
                    
                    if obj.hndl.DataReady
                        break;
                    end
                    
                    pause(0.05);
                    
                end
                
            catch

                pause(2)
                obj.status = 0;
                obj.disconnect;
                obj.connect;
                
            end

            % get current time [JD]
            obj.jd = juliandate(datetime('now', 'timezone', 'UTC'));
            
            obj.light = obj.hndl.DayLightV; % Day light value
            
            obj.temperature = obj.hndl.AmbientT; % Ambient temperature
            
            obj.temp_sky = obj.hndl.RelSkyT; % Sky minus ambient temperature; an indicator for the clouds condition
            
            obj.temp_sensor = obj.hndl.SensorT; % Sensor temperature
            
            obj.temp_dew_point = obj.hndl.DewPointT; % Dew point temperature
            obj.wind_speed = obj.hndl.wind; % Wind speed
            obj.humidity = obj.hndl.HumidityPercent; % Humiditiy percentage
            obj.rain = obj.hndl.RainF; % Rain (in mm?)

            % Text conditions from Boltwood
            obj.light_condition = obj.hndl.DayCondition; % Day condition
            obj.cloud_condition = obj.hndl.CloudCondition; % Cloud condition
            obj.wind_condition = obj.hndl.WindCondition; % Wind condition
            obj.rain_condition = obj.hndl.RainCondition; % Rain condition
            
            obj.data.input([obj.jd, obj.light, obj.temperature, obj.temp_sky, obj.temp_sensor, obj.temp_dew_point, ...
                         obj.wind_speed, obj.humidity, obj.rain]);
            
        end
        
    end
    
end

            
