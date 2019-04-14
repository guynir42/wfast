classdef FocusActuator < handle

% INSTALLATION NOTES:
% In the dropbox HARDWARE/ActuatorDisk library you should install:
% PI_Mercury.CD_Setup.exe
% PI_MATLAB_Driver_GCS2_Setup.exe
% After the first install you should be able to turn on "MikroMove" and
% play with the actuators using their GUI. 
% You MUST make sure each actuator is loaded and calibrated in MikroMove
% After the second install make sure there is a library: 
% C:\Users\Public\PI\PI_MATLAB_Driver_GCS2 (this is hardcoded for laziness)
% 
% HARDWARE SETUP: the controller switches must all be ON except 7 and 8. 
    
    properties % objects and informaion
        
        hndl;
        id = 1;
        
    end
    
    properties % inputs/switches
        
        min_value = 0;
        max_value = 30;
        
    end
    
    properties(Dependent=true) % info outputs
        
        status;
        
        pos;
        
    end
    
    properties(Hidden=true)
        
        axis = '1'; % this is the actuator 1 out of 1 in the single controller (not to be confused with multiple controllers "id")
        controller_serial;
        debug_bit = 1;
        version = 1.02;
        
    end
    
    methods % constructor
        
        function obj = FocusActuator(serial_str)
                       
            if nargin<1 || isempty(serial_str)
                serial_str = '0165500035';
            end
            
            obj.controller_serial = serial_str;
            
            if obj.debug_bit, fprintf('FocusActuator constructor v%4.2f |  axis: %d | S/N: %ld\n', obj.version, obj.axis, obj.controller_serial); end
                        
            addpath('C:\Users\Public\PI\PI_MATLAB_Driver_GCS2');

            obj.connect;
            
        end
        
        function delete(obj)
           
            if obj.debug_bit, disp(['disconnecting from actuator ' obj.controller_serial]); end
            obj.hndl.CloseConnection;
            obj.hndl.Destroy;
            
        end
       
        function connect(obj)
                       
            if isempty(obj.hndl) && ~isa(obj.hndl, 'PI_GCS_Controller'), obj.hndl = PI_GCS_Controller(); end
            
            try
                if ~obj.isConnected, obj.hndl = obj.hndl.ConnectUSB(obj.controller_serial); end
            catch 
                if ~obj.isConnected, obj.hndl = obj.hndl.ConnectUSB(obj.controller_serial); end
            end
            
            obj.hndl = obj.hndl.InitializeController();
            obj.hndl.SVO(obj.axis, 1); % turn on servo...
            obj.hndl.FNL(obj.axis); % measure the reference point
            
        end
        
    end
    
    methods % reset methods
        
    end 
    
    methods % getters
        
        function val = isMoving(obj)
           
            val = obj.hndl.IsMoving(obj.axis);
            
        end
        
        function val = get.pos(obj)
           
            val = obj.hndl.qPOS(obj.axis);
            
        end
        
        function val = getMaxPos(obj)
           
            val = obj.hndl.qTMX(obj.axis);
            
        end     
        
        function val = getMinPos(obj)
           
            val = obj.hndl.qTMN(obj.axis);
            
        end
        
        function val = isConnected(obj)
            
            val = obj.hndl.IsConnected;
            
        end
        
        function val = get.status
            
            val = obj.isConnected; % can we add additional tests on this??
            
        end
        
    end
    
    methods % setters
        
        function set.pos(obj, val)
           
            obj.move(val);
            
        end
        
    end
    
    methods % controls
        
        function rel_move(obj, position)
            
            obj.move(obj.pos + position);
            
        end
        
        function move(obj, position)
           
            if position>obj.max_value
                warning(['cannot move actuator ' num2str(obj.id) ' above max_value= ' num2str(obj.max_value)]);
                position = obj.max_value;
            elseif position<obj.min_value
                warning(['cannot move actuator ' num2str(obj.id) ' below min_value= ' num2str(obj.min_value)]);
                position = obj.min_value;
            end
            
            obj.hndl.MOV(obj.axis, position);
            
        end
        
        function stop(obj)
           
            obj.hndl.STP;
            
        end
        
    end
    
end