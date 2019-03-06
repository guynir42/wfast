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
        data@util.vec.CircularBuffer; % Stack object containing data history
        data_col = {'JD', 'DayLightV', 'AmbientT', 'RelSkyT', 'SensorT', 'DewPointT', 'WindSpeed', 'Humidity', 'Rain'}; % Stack object columns
        data_units = {'day', 'relative', 'deg C', 'deg C', 'deg C', 'deg C', 'km/h', 'percentage', 'mm'}; % Stack object column units

        % specific fields
        jd = NaN; % latest date when data was successfully updated
        light_value = NaN; % Day light value
        temperature = NaN; % Ambient temperature
        temp_sky = NaN; % Sky minus ambient temperature; an indicator for the clouds condition
        temp_sensor = NaN; % Sensor temperature
        temp_dew_point = NaN; % Dew point temperature
        wind_speed = NaN; % Wind speed
        humidity = NaN; % Humiditiy
        rain = NaN; % Rain
        
        % Text conditions from Boltwood
        day_condition   = '';                        % Day condition
        cloud_condition = '';                        % Cloud condition
        wind_condition  = '';                        % Wind condition
        rain_condition  = '';                        % Rain condition 
        
        debug_bit = 1;
        
    end
    
    properties(Hidden=true)
        
        version = 1.00;
        
    end
    
    
    methods % Constructor
        
        function obj = Boltwood(varargin)
            
            % Boltwood class constractor
            % Example: obj = obs.sens.Boltwood
            
            if isempty(varargin)
            
                if obj.debug_bit, fprintf('Boltwood default constructor v%4.2f\n', obj.version); end
                
                obj.data = util.vec.CircularBuffer;
                obj.connect; % Open Boltwood weather application. Also open connection to application if not already opened.
                
            end
            
        end
        

    end
    
    
    methods % getters/setters

      
    end
    
    methods
        function connect(obj) % Open Boltwood weather application. Also open connection to application if not already opened.
            
            obj.hndl = actxserver('ClarityII.CloudSensorII');
            
        end

        function update(obj) % Read weather parameters
            
            % if status is ok than set the time of the last query
            Dummy = obj.hndl.dataReady;
            if (Dummy)
               % Update status
               obj.status = 1;
            else
               obj.status = 0;
               return;
            end
            
            % get current time [JD]
            obj.jd = celestial.time.julday;
            
            % If the weather application is down then reinitiate it
            try
                obj.hndl.dataReady;
            catch
%                obs.sens.Boltwood(false); % call constructor without update
                pause(5)
                obj.status = false;
                obj.open;
            end

            % send query to Boltwood and get answer
            obj.light_value = obj.hndl.DayLightV; % Day light value
            
            obj.temperature = obj.hndl.AmbientT; % Ambient temperature
            
            obj.temp_sky = obj.hndl.RelSkyT; % Sky minus ambient temperature; an indicator for the clouds condition
            
            obj.temp_sensor = obj.hndl.SensorT; % Sensor temperature
            
            obj.temp_dew_point = obj.hndl.DewPointT; % Dew point temperature
            obj.wind_speed = obj.hndl.wind; % Wind speed
            obj.humidity = obj.hndl.HumidityPercent; % Humiditiy percentage
            obj.rain = obj.hndl.RainF; % Rain (in mm?)

            % Text conditions from Boltwood
            obj.day_condition = obj.hndl.DayCondition; % Day condition
            obj.cloud_condition = obj.hndl.CloudCondition; % Cloud condition
            obj.wind_condition = obj.hndl.WindCondition; % Wind condition
            obj.rain_condition = obj.hndl.RainCondition; % Rain condition
            
            obj.data.input([obj.jd, obj.light_value, obj.temperature, obj.temp_sky, obj.temp_sky, obj.temp_dew_point, ...
                         obj.wind_speed, obj.humidity, obj.rain]);
            
        end
        
    end
    
    
    
end

            
