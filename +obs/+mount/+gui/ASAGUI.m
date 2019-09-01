classdef ASAGUI < handle
    
    properties 
        
        owner@obs.mount.ASA; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 12;
        big_font_size = 16;
        edit_font_size = 11;
        small_font_size = 10;
        
        color_on = [0 0.3 1];
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_status;
        panel_pointing;
        panel_motion;
        
        panel_rates;
        plot_axes;
        
        panel_object;
        panel_manual;
        panel_engineering;
        panel_limits;
        panel_arduino;
        
        panel_slew;
        
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = ASAGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'obs.mount.ASA')
                
                if obj.debug_bit, fprintf('ASAGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
            else
                error('Input an obs.mount.ASA to constructor of ASAGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            
            obj.fig = util.plot.FigHandler('Mount ASA');
            obj.fig.clear;
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 25;
            movegui(obj.fig.fig, 'center');
            
            N = 9; % number of button lines
            
            pos = N;
            
            %%%%%%%%%%% panel status %%%%%%%%%%%%%%%
            
            % Add buttons using obj.addButton(button_name, var_name='', type='', str1='', str2='', font_size='', split=1, color_on=[], color_off=[], tooltip)
            
            num_buttons = 6;
            pos = pos-1; % horizontal panel
            obj.panel_status = GraphicPanel(obj.owner, [0 pos/N 1 1/N], 'status', 0); % last input is for horizontal
            obj.panel_status.number = num_buttons;
            obj.panel_status.addButton('button_status', 'status', 'info', 'status= ', '', '', [], '', '', 'Is mount responsive and connected');
            obj.panel_status.addButton('button_lst', 'LST', 'info', 'LST= ', '', '', [], '', '', 'Local Sidereal Time');
            obj.panel_status.addButton('button_placeholder1', '', 'custom');
            obj.panel_status.addButton('button_placeholder2', '', 'custom');
            obj.panel_status.addButton('button_placeholder3', '', 'custom');
            obj.panel_status.addButton('button_connect', 'connect', 'push', 'Connect', '', '', [], '', '', 'Reload Autoslew software and connect to telescope');
            obj.panel_status.margin = [0.005 0.2];
            obj.panel_status.make;
            
            %%%%%%%%%%% panel pointing %%%%%%%%%%%%%%%
            
            num_buttons = 6;
            pos = pos-1; % horizontal panel
            obj.panel_pointing = GraphicPanel(obj.owner, [0 pos/N 1 1/N], 'pointing', 0); % last input is for horizontal
            obj.panel_pointing.number = num_buttons;
            obj.panel_pointing.addButton('button_ra', 'telRA', 'info', 'RA= ', '', '', [], '', '', 'Current pointing Right Ascention');
            obj.panel_pointing.addButton('button_dec', 'telDEC', 'info', 'Dec= ', '', '', [], '', '', 'Current pointing Declination');
            obj.panel_pointing.addButton('button_ha', 'telHA', 'info', 'HA= ', '', '', [], '', '', 'Current pointing Hour Angle');
            obj.panel_pointing.addButton('button_alt', 'telALT', 'info', 'ALT= ', '', '', [], '', '', 'Current pointing Altitude');
            obj.panel_pointing.addButton('button_time', 'tel_time_to_limit', 'info', '', ' min. to limit', '', [], '', '', 'Current time to reach limit');
            obj.panel_pointing.addButton('button_pierside', 'pier_side', 'info', '', '', '', [], '', '', 'Current telescope side of pier');
            obj.panel_pointing.margin = [0.005 0.2];
            obj.panel_pointing.make;
            
            %%%%%%%%%%% panel motion %%%%%%%%%%%%%%%
            
            num_buttons = 3;
            pos = pos-2; % horizontal panel
            obj.panel_motion = GraphicPanel(obj.owner, [0 pos/N 0.4 2/N], 'motion', 1); % last input is for vertical
            obj.panel_motion.addButton('button_morors', '', 'custom', 'Motor', '', '', 0.5, '', '', 'Need to implement motor on/off and feedback from mount!');
            obj.panel_motion.addButton('button_tracking', 'tracking', 'toggle', '', '', '', 0.5, obj.color_on, 'red', 'telescope tracking');
            obj.panel_motion.addButton('button_guiding', 'use_guiding', 'toggle', 'guiding off', 'guiding on', '', 0.5, obj.color_on, 'red', 'Apply rate corrections based on inputs from camera');
            obj.panel_motion.addButton('button_reset_rates', 'resetRate', 'push', 'reset rates', '', '', 0.5, '', '', 'Reset guiding rates');
            obj.panel_motion.number = num_buttons;
            obj.panel_motion.margin = [0.02 0.02];
            obj.panel_motion.make;
            
            %%%%%%%%%%% panel rates %%%%%%%%%%%%%%%
            
            obj.panel_rates = uipanel('Title', 'rates', 'Position', [0.4 pos/N 0.6 2/N]);
            
            %%%%%%%%%%% panel object %%%%%%%%%%%%%%%
            
            num_buttons = 6;
            pos = pos - 4;
            obj.panel_object = GraphicPanel(obj.owner, [0 pos/N 0.3 4/N], 'object', 1); % last input is for vertical 
            obj.panel_object.number = num_buttons;
            obj.panel_object.addButton('input_name', 'objName', 'input text', 'obj= ', '', '', 0.7, '', '', 'Input the name of object/field');
            obj.panel_object.addButton('button_resolve', 'inputTarget', 'push', 'resolve', '', '', 0.3, '', '', 'Use object name to auto fill RA/Dec');
            obj.panel_object.addButton('input_ra', 'objRA', 'input', 'RA= ', '', '', 0.5, '', '', 'Input object Right Ascention or use name resolver');
            obj.panel_object.addButton('input_dec', 'objDEC', 'input', 'Dec= ', '', '', 0.5, '', '', 'Input object Declination or use name resolver');
            obj.panel_object.addButton('button_ha', 'objHA', 'info', 'HA= ', '', '', 0.5, '', '', 'Calculated object Hour Angle'); 
            obj.panel_object.addButton('button_alt', 'objALT', 'info', 'ALT= ', '', '', 0.5, '', '', 'Calculated object Altitude'); 
            obj.panel_object.addButton('button_time', 'obj_time_to_limit', 'info', '', ' min. to limit', '', 0.5, '', '', 'Calculated time object has until reaching limit'); 
            obj.panel_object.addButton('button_pierside', 'obj_pier_side', 'info', '', '', '', 0.5); 
            obj.panel_object.margin = [0.02 0.02];
            obj.panel_object.make;
            
            %%%%%%%%%%% panel manual %%%%%%%%%%%%%%%
            
            num_buttons = 4;
            obj.panel_manual = GraphicPanel(obj.owner, [0.3 pos/N 0.3 num_buttons/N], 'manual move', 1); % last input is for vertical 
            obj.panel_manual.addButton('button_NW', '', 'custom', 'NW', '', '', 1/3, '', '', 'Move the telescope to the North West');
            obj.panel_manual.addButton('button_N', '', 'custom', 'N', '', '', 1/3, '', '', 'Move the telescope to the North');
            obj.panel_manual.addButton('button_NE', '', 'custom', 'NE', '', '', 1/3, '', '', 'Move the telescope to the North East');
            obj.panel_manual.addButton('button_W', '', 'custom', 'W', '', '', 1/3, '', '', 'Move the telescope to the West');
            obj.panel_manual.addButton('input_rate', 'move_rate', 'input', ' ', 'deg/s', '', 1/3, '', '', 'Control the manual slew rate (deg/sec)');
            obj.panel_manual.addButton('button_E', '', 'custom', 'E', '', '', 1/3, '', '', 'Move the telescope to the East');
            obj.panel_manual.addButton('button_SW', '', 'custom', 'SW', '', '', 1/3, '', '', 'Move the telescope to the South West');
            obj.panel_manual.addButton('button_S', '', 'custom', 'S', '', '', 1/3, '', '', 'Move the telescope to the South');
            obj.panel_manual.addButton('button_SE', '', 'custom', 'SE', '', '', 1/3, '', '', 'Move the telescope to the South East');
            obj.panel_manual.number = num_buttons;
            obj.panel_manual.margin = [0.05 0.05];
            obj.panel_manual.make;
            
            obj.panel_manual.button_NW.Enable = 'inactive';
            obj.panel_manual.button_NW.control.ButtonDownFcn = @obj.callback_NW;
            obj.panel_manual.button_N.Enable = 'inactive';
            obj.panel_manual.button_N.control.ButtonDownFcn = @obj.callback_N;
            obj.panel_manual.button_NE.Enable = 'inactive';
            obj.panel_manual.button_NE.control.ButtonDownFcn = @obj.callback_NE;
            obj.panel_manual.button_W.Enable = 'inactive';
            obj.panel_manual.button_W.control.ButtonDownFcn = @obj.callback_W;
            obj.panel_manual.button_E.Enable = 'inactive';
            obj.panel_manual.button_E.control.ButtonDownFcn = @obj.callback_E;
            obj.panel_manual.button_SW.Enable = 'inactive';
            obj.panel_manual.button_SW.control.ButtonDownFcn = @obj.callback_SW;
            obj.panel_manual.button_S.Enable = 'inactive';            
            obj.panel_manual.button_S.control.ButtonDownFcn = @obj.callback_S;
            obj.panel_manual.button_SE.Enable = 'inactive';
            obj.panel_manual.button_SE.control.ButtonDownFcn = @obj.callback_SE;
                        
            %%%%%%%%%%% panel limits %%%%%%%%%%%%%%%
            
            num_buttons = 1;
            obj.panel_limits = GraphicPanel(obj.owner, [0.6 (pos+3)/N 0.2 num_buttons/N], 'limits', 1); % last input is for vertical 
            obj.panel_limits.number = num_buttons;
            obj.panel_limits.addButton('button_alt', 'limit_alt', 'info', 'alt limit= ', '', 'small', 0.5, '', '', 'Minimal angle above horizon for telescope motion (degrees)');
            obj.panel_limits.addButton('button_flip', 'limit_flip', 'info', 'flip limit= ', '', 'small', 0.5, '', '', 'Maximum angle beyond meridian after which telescope must flip (degrees)');
            obj.panel_limits.margin = [0.03 0.08];
            obj.panel_limits.make;
            
            %%%%%%%%%%% panel engineering %%%%%%%%%%%%%%%
            
            num_buttons = 3;
            obj.panel_engineering = GraphicPanel(obj.owner, [0.6 (pos)/N 0.2 num_buttons/N], 'engineering slews', 1); % last input is for vertical 
            obj.panel_engineering.number = num_buttons;
            obj.panel_engineering.addButton('button_park', 'park', 'push', 'Park', '', '', [], '', '', 'Send the telescope to park 1 position'); 
            obj.panel_engineering.margin = [0.02 0.02];
            obj.panel_engineering.make;
            
            %%%%%%%%%%% panel arduino %%%%%%%%%%%%%%%
            
            num_buttons = 6;
            obj.panel_arduino = GraphicPanel(obj.owner, [0.8 pos/N 0.2 4/N], 'arduino', 1); % last input is for vertical 
            obj.panel_arduino.number = num_buttons;
            obj.panel_arduino.addButton('button_status', 'ard.status', 'info', 'Status= ', '', '', 0.5, '', '', 'Status of communication with arduino');
            obj.panel_arduino.addButton('button_connect', 'connectArduino', 'push', 'Connect', '', '', 0.5, '', '', 'Attempt to reconnect with arduino');
            obj.panel_arduino.addButton('button_use_accel', 'use_accelerometer', 'toggle', 'accel. off ', 'use accel.', '', 0.5, obj.color_on, 'red', 'Use arduino accelerometer to stop telescope at low altitude');
            obj.panel_arduino.addButton('button_angle', 'ard.ALT', 'info', 'ALT= ', '', 'small', 0.5, '', '', 'Current measured Altitude angle (degrees)');
            obj.panel_arduino.addButton('button_use_ultra', 'use_ultrasonic', 'toggle', 'ultra. off ', 'use ultra.', '', 0.5, obj.color_on, 'red', 'Use arduino ultrasonic sensor to warn agains obstacles in front of telescope');
            obj.panel_arduino.addButton('button_distance', 'ard.distance', 'info', 'dist= ', '', 'small', 0.5, '', '', 'Current measured distance to obstructions (cm)');
            obj.panel_arduino.margin = [0.02 0.02];
            obj.panel_arduino.make;
            
            %%%%%%%%%%% panel slew %%%%%%%%%%%%%%%
            
            pos = pos - 1;
            obj.panel_slew = GraphicPanel(obj.owner, [0 pos/N 1 1/N], '', 1);
            obj.panel_slew.number = 1;
            obj.panel_slew.addButton('button_slew', 'slew', 'push', 'Slew', '', '', 0.2, '', '', 'Immediately start moving telescope to given target');
            obj.panel_slew.addButton('button_stop', 'stop', 'push', 'STOP', '', '', 0.6, '', '', 'Stop current slew and reset guiding rates');
            obj.panel_slew.addButton('button_sync', 'syncToTarget', 'push', 'Sync', '', '', 0.2, '', '', 'Sync current telescope position to current object');
            obj.panel_slew.margin = [0.005 0.1];
            obj.panel_slew.make;
            
            obj.panel_slew.button_slew.Callback = @obj.callback_slew;
            
            obj.fig.fig.WindowButtonUpFcn = @obj.callback_button_up;
            obj.fig.fig.CloseRequestFcn = @obj.callback_close_fig;
            
            obj.update;
            
        end
        
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            if isempty(obj.plot_axes) || ~isvalid(obj.plot_axes)
                delete(obj.panel_rates.Children);
                obj.plot_axes = axes('Parent', obj.panel_rates);
            end
            
            obj.owner.plot_rate(obj.plot_axes);
            
            if obj.owner.objALT<15
                obj.panel_object.button_alt.BackgroundColor = 'red';
            else
                obj.panel_object.button_alt.BackgroundColor = util.plot.GraphicButton.defaultColor;
            end
            
            if strcmp(obj.owner.obj_pier_side, 'pierEast')
                obj.panel_object.button_pierside.Tooltip = 'Object is on West side. To view it, telescope must be on East side';
            elseif strcmp(obj.owner.obj_pier_side, 'pierWest')
                obj.panel_object.button_pierside.Tooltip = 'Object is on East side. To view it, telescope must be on West side';
            elseif strcmp(obj.owner.obj_pier_side, 'pierUnknown')
                obj.panel_object.button_pierside.Tooltip = 'Object is below horizon (or other error)';
            else
                obj.panel_object.button_pierside.Tooltip = 'Unknown pier side. Must be an error'; 
            end
            
            
            if strcmp(obj.owner.pier_side, obj.owner.obj_pier_side)
                obj.panel_object.button_pierside.BackgroundColor = util.plot.GraphicButton.defaultColor;
            elseif ~strcmp(obj.owner.obj_pier_side, 'pierUnknown')
                obj.panel_object.button_pierside.BackgroundColor = 'red';
                obj.panel_object.button_pierside.Tooltip = [obj.panel_object.button_pierside.Tooltip ' (need to flip!)']; 
            end
            
            if obj.owner.hndl.Slewing
                obj.panel_slew.button_stop.BackgroundColor = 'red';
            else
                obj.panel_slew.button_stop.BackgroundColor = util.plot.GraphicButton.defaultColor;
            end
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_slew) && isvalid(obj.panel_slew.panel);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_NW(obj, ~, ~)
            
            if obj.debug_bit>1, disp('Callback: moving NW'); end
            
            obj.owner.hndl.MoveAxis(0,obj.owner.move_rate);
            obj.owner.hndl.MoveAxis(1,obj.owner.move_rate);
            
        end
        
        function callback_N(obj, ~, ~)
            
            if obj.debug_bit>1, disp('Callback: moving NW'); end
            
%             obj.owner.hndl.MoveAxis(0,obj.owner.move_rate);
            obj.owner.hndl.MoveAxis(1,obj.owner.move_rate);
            
        end
        
        function callback_NE(obj, ~, ~)
            
            if obj.debug_bit>1, disp('Callback: moving NW'); end
            
            obj.owner.hndl.MoveAxis(0,-obj.owner.move_rate);
            obj.owner.hndl.MoveAxis(1,obj.owner.move_rate);
            
        end
        
        function callback_W(obj, ~, ~)
            
            if obj.debug_bit>1, disp('Callback: moving NW'); end
            
            obj.owner.hndl.MoveAxis(0,obj.owner.move_rate);
%             obj.owner.hndl.MoveAxis(1,obj.owner.move_rate);
            
        end
        
        function callback_E(obj, ~, ~)
            
            if obj.debug_bit>1, disp('Callback: moving NW'); end
            
            obj.owner.hndl.MoveAxis(0,-obj.owner.move_rate);
%             obj.owner.hndl.MoveAxis(1,obj.owner.move_rate);
            
        end
        
        function callback_SW(obj, ~, ~)
            
            if obj.debug_bit>1, disp('Callback: moving NW'); end
            
            obj.owner.hndl.MoveAxis(0,obj.owner.move_rate);
            obj.owner.hndl.MoveAxis(1,-obj.owner.move_rate);
            
        end
        
        function callback_S(obj, ~, ~)
            
            if obj.debug_bit>1, disp('Callback: moving NW'); end
            
%             obj.owner.hndl.MoveAxis(0,obj.owner.move_rate);
            obj.owner.hndl.MoveAxis(1,-obj.owner.move_rate);
            
        end
        
        function callback_SE(obj, ~, ~)
            
            if obj.debug_bit>1, disp('Callback: moving NW'); end
            
            obj.owner.hndl.MoveAxis(0,-obj.owner.move_rate);
            obj.owner.hndl.MoveAxis(1,-obj.owner.move_rate);
            
        end
        
        function callback_button_up(obj, ~, ~)
            
            if obj.debug_bit>1, disp('Callback: stopping manual move'); end
            
%             obj.owner.hndl.MoveAxis(0,0);
%             obj.owner.hndl.MoveAxis(1,0);
            obj.owner.hndl.AbortSlew;
            obj.owner.hndl.Tracking = obj.owner.was_tracking; % go around the ASA.set.tracking function 
            
        end
        
        function callback_close_fig(obj, ~, ~)
            
            if obj.debug_bit>1, disp('Callback: stopping manual move and closing GUI'); end
            
%             obj.owner.hndl.MoveAxis(0,0);
%             obj.owner.hndl.MoveAxis(1,0);
            
            obj.owner.hndl.AbortSlew;
            obj.owner.hndl.Tracking = obj.owner.was_tracking; % go around the ASA.set.tracking function 
            
            delete(obj.fig.fig);
            
        end
                
        function callback_slew(obj, ~, ~)
            
            if obj.debug_bit, disp('Callback: slew'); end
            
            if obj.owner.check_need_flip
                res = questdlg('Need to flip for this target. Are you sure?', 'Flip needed!', 'Slew', 'Abort', 'Slew');
                if isempty(res) || strcmp(res, 'Abort')
                    return;
                end
            end
            
            obj.owner.slew;
            
            obj.update;
            
        end
        
    end
    
end