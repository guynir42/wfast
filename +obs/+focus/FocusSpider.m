classdef FocusSpider < handle

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

    properties(Transient=true)
        
        gui@obs.focus.gui.SpiderGUI;
        gui_cam; % handle to the camera GUI for feedback
        
    end

    properties % internal info and objects
        
        actuators@obs.focus.FocusActuator;
        
        serial_numbers = {'0165500035', '0165500253', '0165500280'};
        
        % direction vectors for tip/tilt action
        tip_vec = [-0.5 1 -0.5];
        tilt_vec = [0.5 0 -0.5];
        
    end
    
    properties % switches and controls
    
        num_act = 3;
        
        min_pos = 0;
        max_pos = 20;
        max_tip_tilt = 5; % this is absolute tip/tilt        
        
        step = 0.05;
        step_tip = 0.01;
        step_tilt = 0.01;
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
       
        status;
        
        pos;
        tip;
        tilt;
        
    end
    
    properties(Hidden = true)
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = FocusSpider
            
            if obj.debug_bit, fprintf('FocusSpider constructor v%4.2f\n', obj.version); end
            
            obj.connect;
            
        end
        
        function connect(obj)
            
            for ii = 1:obj.num_act
               
                obj.actuators(ii) = obs.focus.FocusActuator(obj.serial_numbers{ii});
                obj.actuators(ii).id = ii;
                
            end
            
        end
        
    end
    
    methods % reset
    
    end
    
    methods % getters
        
        function val = get.status(obj)
            
            val = all([obj.actuators.status]);
            
        end
        
        function val = get.pos(obj)
            
            val = mean([obj.actuators.pos]);
            
        end
        
        function val = get.tip(obj)
            
            val = 0;
            
            p = obj.pos;
            dp = [obj.actuators.pos] - p;
            n = 0; 
            for ii = 1:obj.num_act
                if obj.tip_vec(ii)~=0 
                    val = val + dp(ii)./obj.tip_vec(ii); 
                    n = n+1;
                end
            end
            
            val = val./n;

        end
        
        function val = get.tilt(obj)
            
            val = 0;
            
            p = obj.pos;
            dp = [obj.actuators.pos] - p;
            n = 0; 
            for ii = 1:obj.num_act
                if obj.tilt_vec(ii)~=0 
                    val = val + dp(ii)./obj.tilt_vec(ii); 
                    n = n+1;
                end
            end
            
            val = val./n;  
            
        end
        
    end
    
    methods % setters
        
        function set.pos(obj, val)
            
            if val>obj.max_pos
                warning(['cannot set position to ' num2str(val) ' it is above max_pos= ' num2str(obj.max_pos)]);
                val = obj.max_pos;
            elseif val<obj.min_pos
                warning(['cannot set position to ' num2str(val) ' it is below min_pos= ' num2str(obj.min_pos)]);
                val = obj.min_pos;
            end
            
            p = obj.pos;
            dp = val-p;
            
            for ii = 1:obj.num_act
                obj.actuators(ii).rel_move(dp);
            end
            
        end
        
        function set.tip(obj, val)
            
            if val>obj.max_tip_tilt
                warning(['cannot set tip to ' num2str(val) ' max is ' num2str(obj.max_tip_tilt)]);
                val = obj.max_tip_tilt;
            elseif val<-obj.max_tip_tilt                
                warning(['cannot set tip to ' num2str(val) ' min is ' num2str(-obj.max_tip_tilt)]);
                val = -obj.max_tip_tilt;
            end
            
            t = obj.tip;
            dt = val-t;
            
            for ii = 1:obj.num_act
                obj.actuators(ii).rel_move(dt*obj.tip_vec(ii));
            end
            
        end
        
        function set.tilt(obj, val)
            
            if val>obj.max_tip_tilt
                warning(['cannot set tip to ' num2str(val) ' max is ' num2str(obj.max_tip_tilt)]);
                val = obj.max_tip_tilt;
            elseif val<-obj.max_tip_tilt                
                warning(['cannot set tip to ' num2str(val) ' min is ' num2str(-obj.max_tip_tilt)]);
                val = -obj.max_tip_tilt;
            end
            
            t = obj.tilt;
            dt = val-t;
            
            for ii = 1:obj.num_act
                obj.actuators(ii).rel_move(dt*obj.tilt_vec(ii));
            end
                        
        end
        
    end
    
    methods % commands
        
        function posRelativeMove(obj, val)
            
            obj.pos = obj.pos + val;
            
        end
        
        function tipRelativeMove(obj, val)
            
            obj.tip = obj.tip + val;
            
        end
        
        function tiltRelativeMove(obj, val)
            
            obj.tilt = obj.tilt + val;
            
        end
        
        function home(obj, pos)
           
            if nargin<2 || isempty(pos)
                pos = 0;
            end
            
            for ii = 1:obj.num_act
                
                try
                    obj.actuators(ii).move(pos);
                catch ME
                    disp(['error with actuator number ' num2str(ii)]);
                    rethrow(ME);
                end
                
            end
            
        end
        
        function up(obj)
            
            obj.pos = obj.pos + obj.step;
            
        end
        
        function down(obj)
            
            obj.pos = obj.pos - obj.step;
            
        end
        
        function tip_up(obj)
            
            obj.pos = obj.tip + obj.step_tip;
            
        end
        
        function tip_down(obj)
            
            obj.pos = obj.tip - obj.step_tip;
            
        end
        
        function tilt_up(obj)
            
            obj.pos = obj.tilt + obj.step_tilt;
            
        end
        
        function tilt_down(obj)
            
            obj.pos = obj.tilt - obj.step_tilt;
            
        end
        
        
        
        function demo(obj)
           
            obj.home(30);
            pause(2);
            
            obj.home(0);
            pause(2);
            
            obj.home(20);
            
            pause(2);
            
            theta = 0:10:360;
            
            for t = theta
               
                obj.tip  = 5*cosd(t);
                obj.tilt = 5*sind(t);
                
                pause(0.3);
                
            end
            
            pause(2);
            
            obj.home;
            
        end
        
    end
    
    methods % plotting / printouts
       
        function printout(obj)
           
            str = sprintf('pos= %4.2f | tip= %4.2f | tilt= %4.2f', obj.pos, obj.tip, obj.tilt);
            
            for ii = 1:obj.num_act
                str = [str sprintf(' | a%d.pos= %4.2f', ii, obj.actuators(ii).pos)];
            end
            
%             str = [str sprintf('\n')];
            
            if nargout<1
                disp(str);
            end
            
        end
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = obs.focus.gui.SpiderGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
    end
    
end