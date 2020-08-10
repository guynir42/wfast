classdef DomeAssistant < handle
% Class to control the "dome assistant", an Arduino board that has both a 
% set of 8 relays and a temperature-humidity-pressure sensor. 
% The object behaves like any sensor, with a CircularBuffer that has the data
% and some dependent properties like temperature_average that integrate over 
% the last 10 measurements. 
% The SensorChecker should call update() every minute or so, and collect the 
% average measurements when testing the weather. 
%
% To control devices, either use the send() command directly to the Arduino
% parser, or use some of the shortcuts added along with the devices. 
% E.g., the dome lights are controlled using the dependent property "lights". 
% When true, the light is on. You can set it on or off or set it to 'timer' 
% to make it shut off on its own after some time (20 minutes is the default) 
%
% For more information, see the sketch for this Arduino named DomeAssistant, 
% and also the parsing functions in MyArduino/OutputPin.cpp and Timer.cpp.
    
    properties(Transient=true)
       
        time;
        hndl; % serial port
        
    end
    
    properties % objects
        
        data@util.vec.CircularBuffer; % object containing data history for the last few measurements
        
    end
    
    properties % inputs/outputs
        
        data_col = {'JD', 'Temperature', 'Humidity', 'Pressure'}; 
        data_units = {'day', 'deg C', 'percentage', 'mbar'}; % Stack object column units
        
        reply = '';
        
        status = 0; % false - readings are unreliable, true - ok
        
        id = 'assist'; 
        
        JD;
        temperature;
        humidity;
        pressure;
        
        relays;
        
    end
    
    properties % switches/controls
        
        port_number = 14; % e.g., COM7
        timeout = 10; % seconds
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        temperature_average;
        humidity_average;
        pressure_average;
        
        lights; % dome LED lights using the relays
        
    end
    
    properties(Hidden=true)
       
        lights_relay = 'relay1';
        lights_switch = 4; 
        
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = DomeAssistant(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.sens.DomeAssistant')
                if obj.debug_bit>1, fprintf('DomeAssistant copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('DomeAssistant constructor v%4.2f\n', obj.version); end
                
                obj.data = util.vec.CircularBuffer;
                obj.connect; % Open serial port
                obj.lights = 0; 
                
            end
            
        end
        
        function connect(obj)
            
            obj.disconnect;
            
            obj.hndl = serial(sprintf('COM%d', obj.port_number), 'Timeout', obj.timeout, 'Terminator','LF'); 
            
            obj.hndl.BytesAvailableFcnMode = "terminator";
            
            obj.hndl.BytesAvailableFcn = @obj.read_data;
            
            fopen(obj.hndl); 
            
            obj.status = 1;
            
            obj.update; 
            
        end
        
        function disconnect(obj)
            
            delete(obj.hndl);
            obj.hndl = [];
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.reply = '';
            obj.time = [];
            
            obj.JD = [];
            obj.temperature = [];
            obj.humidity = [];
            obj.pressure = [];
            obj.relays = [];
            
            obj.data.reset;
            
        end
        
    end
    
    methods % getters
        
        function val = get.temperature_average(obj)
            
            val = nanmean(obj.data.data(:,2)); 
            
        end
        
        function val = get.humidity_average(obj)
            
            val = nanmean(obj.data.data(:,3)); 
            
        end
        
        function val = get.pressure_average(obj)
            
            val = nanmean(obj.data.data(:,4)); 
            
        end
        
        function val = get.lights(obj)
            
            if isempty(obj.relays)
                val = [];
            else
                val = obj.relays(1); 
            end
            
        end
        
    end
    
    methods % setters
        
        function set.lights(obj, val)
            
            if isempty(val)
                obj.send([obj.lights_relay, ', off']); 
            elseif ischar(val) && util.text.cs(val, 'timer')
                obj.send([obj.lights_relay, ', watch']); 
                obj.send([obj.lights_relay, ', mode, expire']); 
            elseif util.text.parse_bool(val)
                obj.send([obj.lights_relay, ', on']); 
            elseif util.text.parse_bool(val)==0
                obj.send([obj.lights_relay, ', off']); 
            else
                error('Unknown option "%s". Use "on" or "off" or "timer"', val); 
            end
            
            obj.update;
            
        end
        
    end
    
    methods % calculations
        
        function update(obj)
            
            obj.send('measure'); 
            % need to add some way to figure out if status is ok... 
            
        end
        
        function send(obj, val)
            
            fwrite(obj.hndl, sprintf('%s\n', val)); 
            
        end
        
        function read_data(obj, ~, ~)
            
            obj.reply = strtrim(fgetl(obj.hndl));
            
            obj.time = datetime('now', 'TimeZone', 'UTC'); 
            
            if ~isempty(regexp(obj.reply, 'P= \d+ | T= \d+ | H= \d+ | relays: \d+', 'once')) % check the format of the respone is what we expect
            
                vec = sscanf(obj.reply, 'P= %f | T= %f | H= %f');
                
                obj.pressure = vec(1);
                obj.temperature = vec(2);
                obj.humidity = vec(3); 
                
                obj.data.input([juliandate(obj.time), obj.temperature, obj.humidity, obj.pressure]);
                
                idx = regexp(obj.reply, 'relays: \d+'); 
                
                obj.relays = double(obj.reply(idx+8:end))-double('0'); 
                
            end
            
        end
        
    end
    
end

