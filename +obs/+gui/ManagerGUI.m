classdef ManagerGUI < handle
    
    properties 
        
        owner@obs.Manager; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        menus = {};
        panels = {};
        
        font_size = 12;
        big_font_size = 16;
        edit_font_size = 10;
        small_font_size = 8;
        
        color_on = [0 0.3 1];
        
        color_bg = 0.5*[1 1 1];
        color_info = 0.6*[1 1 1]; 
        color_input = 0.85*[1 1 1]; 
        color_buttons = 0.7*[1 1 1]; 
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_devices;
        
        panel_telescope;
        
        panel_scheduler;
        
        panel_object;
        
        panel_dome;
        
        panel_controls;
        
        panel_report;
        
        panel_weather;
        
        panel_camera;
        
        panel_stop;
        
        panel_image;
        button_reset_axes;
        button_mean_only;
        button_clicker
        axes_image;
        
        menu_options;
        
        menu_objects;
    
    end
    
    properties (Hidden=true)
              
        key_status_shift = 0; 
        
        version = 1.02;
        
    end
            
    methods % constructor
       
        function obj = ManagerGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'obs.Manager')
                
                if obj.debug_bit>1, fprintf('ManagerGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input an obs.Manager to constructor of ManagerGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.MenuItem;

            obj.buttons = {};
            obj.menus = {};
            obj.panels = {}; 
            
            obj.fig = util.plot.FigHandler('Observatory Manager');
            obj.fig.clear;
            obj.fig.maximize
            
            set(obj.fig.fig, 'WindowKeyPressFcn', @obj.callback_key_press);
            set(obj.fig.fig, 'KeyReleaseFcn', @obj.callback_key_released);
            
            N_left = 20;
            pos = N_left;
            
            %%%%%%%%%%% panel devices %%%%%%%%%%%%%%%%
            
            N = 2;
            pos = pos - N;
            
            obj.panel_devices = GraphicPanel(obj.owner, [0.0 pos/N_left 0.2 N/N_left], 'device status');
            obj.panel_devices.number = N;
            obj.panel_devices.addButton('button_dome', 'dome', 'push', 'dome', '', '', 1/3);
            obj.panel_devices.addButton('button_mount', 'mount', 'push', 'mount', '', '', 1/3);
            obj.panel_devices.addButton('button_weather', 'weather', 'push', 'BoltWood', '', '', 1/3);
            obj.panel_devices.addButton('button_wind', 'wind', 'push', 'WindETH', '', '', 1/3);
            obj.panel_devices.margin = [0.01 0.005];
            obj.panel_devices.make;
            
            %%%%%%%%%%% panel telescope %%%%%%%%%%%%%%%
            
            N = 4;
            pos = pos - N;
            obj.panel_telescope = GraphicPanel(obj.owner, [0 pos/N_left 0.2 N/N_left], 'telescope');
            obj.panel_telescope.number = N;
            obj.panel_telescope.addButton('button_RA', 'RA', 'info', 'RA: ', '', '', 0.5);
            obj.panel_telescope.addButton('button_DE', 'DEC', 'info', 'DE: ', '', '', 0.5);            
            obj.panel_telescope.addButton('button_LST', 'LST', 'info', 'LST: ', '', '', 0.5);
            obj.panel_telescope.addButton('button_ALT', 'ALT', 'info', 'ALT: ', ' deg', '', 0.5);
            obj.panel_telescope.addButton('button_tracking', 'tracking', 'toggle', 'tracking is off', 'tracking is on', '', 0.5, obj.color_on, 'red');
            obj.panel_telescope.addButton('button_side', 'mount.telHemisphere', 'info', 'pointing ', '', '', 0.5); 
            obj.panel_telescope.addButton('button_vibrations', '', 'custom', '', '', '', 0.5); 
            obj.panel_telescope.margin = [0.02 0.01];
            obj.panel_telescope.make;
            
            obj.panel_telescope.button_RA.Tooltip = 'Current right ascention of mount';
            obj.panel_telescope.button_DE.Tooltip = 'Current declination of mount';
            obj.panel_telescope.button_LST.Tooltip = 'Local Sidereal Time';
            obj.panel_telescope.button_ALT.Tooltip = 'Current altitidue of mount';
            obj.panel_telescope.button_tracking.Tooltip = 'Turn telescope tracking on and off';
            obj.panel_telescope.button_side.Tooltip = 'The side to which the telescope is currently pointing';
            
            %%%%%%%%%%% panel scheduler %%%%%%%%%%%%%%%
            
            N = 1;
            pos = pos - 1;
            obj.panel_scheduler = GraphicPanel(obj.owner, [0 pos/N_left 0.2 N/N_left], 'scheduler'); 
            obj.panel_scheduler.addButton('button_gui', 'sched', 'push', 'GUI', '', '', 1/3, '', '', 'open the scheduler GUI for more details'); 
            obj.panel_scheduler.addButton('button_constraints', 'sched.constraintsGUI', 'push', 'constraints', '', '', 1/3, '', '', 'open the constraints menu for choosing targets from scheduler'); 
            obj.panel_scheduler.addButton('button_choose', 'chooseNewTarget', 'push', 'choose', '', '', 1/3, '', '', 'choose a target from the scheduler and load it into the "object" field'); 
            obj.panel_scheduler.margin = [0.02 0.01];
            obj.panel_scheduler.make;
            
            %%%%%%%%%%% panel object %%%%%%%%%%%%%%%
            
            N = 5;
            pos = pos - N;
            obj.panel_object = GraphicPanel(obj.owner, [0 pos/N_left 0.2 N/N_left], 'object');            
            obj.panel_object.number = N;
            obj.panel_object.addButton('button_name', 'mount.objName', 'input_text', ' ', '', 'edit', 0.8, '', '', 'Input the name of the object/field'); 
            obj.panel_object.addButton('button_resolve', 'mount.inputTarget', 'push', 'resolve', '', 'edit', 0.2, '', '', 'Resolve name with SIMBAD'); 
            obj.panel_object.addButton('button_ra', 'mount.objRA', 'input text', 'RA= ', '', 'edit', 0.5, '', '', 'Target right ascention');
            obj.panel_object.addButton('button_dec', 'mount.objDEC', 'input text', 'DE= ', '', 'edit', 0.5, '', '', 'Target declination');
            obj.panel_object.addButton('button_prev_objects', '', 'custom', '', '', '', [], '', '', 'List the last objects that were used to for slew');
            obj.panel_object.addButton('button_ALT', 'mount.objALT', 'info', 'ALT= ', '', '', 0.4, '', '', 'Target altitute above horizong (degrees)');
            obj.panel_object.addButton('button_time', 'mount.obj_time_to_limit', 'info', 'lim.= ', ' min', '', 0.4, '', '', 'Calculated time object has until reaching limit'); 
            obj.panel_object.addButton('button_pierside', 'mount.objHemisphere', 'info', ' ', '', '', 0.2, '', '', 'Side of the sky where the object is right now');
            obj.panel_object.addButton('button_slew', 'mount.slewAskFlip', 'push', 'Slew', '', '', 1, '', '', 'Slew the telescope to the given object'); 
            obj.panel_object.margin = [0.02 0.01];
            obj.panel_object.make;
            
            obj.panel_object.button_prev_objects.control.Style = 'popupmenu';
            obj.panel_object.button_prev_objects.Callback = @obj.callback_prev_objects;
            obj.panel_object.button_slew.Callback = @obj.callback_slew;
            
            
            %%%%%%%%%%% panel dome %%%%%%%%%%%%%%%
            
            N = 4;
            pos = pos - N;
            obj.panel_dome = GraphicPanel(obj.owner, [0 pos/N_left 0.2 N/N_left], 'dome');
            obj.panel_dome.number = N;            
            obj.panel_dome.addButton('button_close_dome', 'closeDome', 'push', 'Close Dome');
            obj.panel_dome.addButton('button_shutter_east', '', 'custom', 'East: ', '', 'edit', 0.275);
            obj.panel_dome.addButton('button_tracking', 'dome.use_tracking', 'toggle', 'dome not tracking', 'dome is tracking', 'edit', 0.45, obj.color_on, '', 'dome can be set to slowly open West and close East'); 
            obj.panel_dome.addButton('button_shutter_west', '', 'custom', 'West: ', '', 'edit', 0.275);
            obj.panel_dome.addButton('button_close_east', '', 'custom', 'Close East', '', '', 0.5);
            obj.panel_dome.addButton('button_close_west', '', 'custom', 'Close West', '', '', 0.5);
            obj.panel_dome.addButton('button_open_east', '', 'custom', 'Open East', '', '', 0.5);
            obj.panel_dome.addButton('button_open_west', '', 'custom', 'Open West', '', '', 0.5);
            obj.panel_dome.margin = [0.02 0.01];
            obj.panel_dome.make;
            
            obj.panel_dome.button_close_dome.Tooltip = 'Immediately close both shutters';
            obj.panel_dome.button_shutter_west.Tooltip = 'Position of West shutter';
            obj.panel_dome.button_shutter_east.Tooltip = 'Position of East shutter';
            obj.panel_dome.button_close_west.Callback = @obj.callback_close_west; 
            obj.panel_dome.button_open_west.Callback = @obj.callback_open_west; 
            obj.panel_dome.button_close_east.Callback = @obj.callback_close_east; 
            obj.panel_dome.button_open_east.Callback = @obj.callback_open_east; 
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            N = 4;
            pos = pos - N;
            obj.panel_controls = GraphicPanel(obj.owner, [0 pos/N_left 0.2 N/N_left], 'controls');
            obj.panel_controls.number = N;
            obj.panel_controls.addButton('button_autoshutdown', 'use_shutdown', 'toggle', 'auto shutdown is disabled', 'auto shutdown is enabled', '', 0.7, obj.color_on, 'red'); 
%             obj.panel_controls.addButton('button_twilight', 'checker.use_twilight_mode', 'toggle', 'twilight is off', 'twilight is on', '', 0.3, 'red', obj.color_on);
            obj.panel_controls.addButton('button_twilight', '', 'custom', 'twilight is off', 'twilight is on', '', 0.3, 'red', obj.color_on);
            
            obj.panel_controls.addButton('button_autostartup', 'use_startup', 'toggle', 'auto start up is disabled', 'auto start up is enabled', '', 0.7, obj.color_on, 'red'); 
            obj.panel_controls.addButton('button_lights', 'assist.lights', 'toggle', 'LEDs are off', 'LEDs are on', 'edit', 0.3, 'red', obj.color_on);
            
            obj.panel_controls.addButton('button_weather_check', 'callback_t2', 'push', 'Weather check');
            obj.panel_controls.addButton('button_proceed', 'checkNewTarget', 'push', 'proceed to target'); 
            obj.panel_controls.margin = [0.01 0.01];
            obj.panel_controls.make;
            
            obj.panel_controls.button_autoshutdown.Tooltip = 'Allow manager to close dome and stop tracking if weather is bad or if there is a device failure';
            obj.panel_controls.button_twilight.Tooltip = 'Let the dome stay open during twilight';
            obj.panel_controls.button_twilight.Callback = @obj.callback_twilight_mode;
            obj.panel_controls.button_autostartup.Tooltip = 'Allow manager to open and start observing autonomously. (not implemented yet)';
            obj.panel_controls.button_weather_check.Tooltip = 'Run t2 to check weather and devices'; 
            obj.panel_controls.button_proceed.Tooltip = 'use scheduler to move to new target'; 
            
            %%%%%%%%%%% panel report %%%%%%%%%%%%%%%%%
            
            N_middle = 20;
            pos = N_middle;
            
            N = 1;
            pos = pos - N;
            
            obj.panel_report = GraphicPanel(obj.owner, [0.2, pos/N_middle, 0.8, N/N_middle], 'report');
            obj.panel_report.addButton('button_report', 'report', 'info');
            obj.panel_report.make;
            
            %%%%%%%%%%% panel weather %%%%%%%%%%%%%%%%
            
            N = 2;
            pos = pos - N;
            
            obj.panel_weather = GraphicPanel(obj.owner, [0.2, pos/N_middle, 0.8, N/N_middle], 'weather');
            obj.panel_weather.addButton('button_temperature', 'average_temperature', 'info', 'Amb. Temp= ', 'C', '', 1/4);
            obj.panel_weather.addButton('button_clouds', 'average_clouds', 'info', 'dT= ', 'C', '', 1/4);
            obj.panel_weather.addButton('button_light', 'average_light', 'info', 'Light= ', '', '', 1/4);
            obj.panel_weather.addButton('button_pressure', 'average_pressure', 'info', 'Pres= ', '', '', 1/4); 
            
            obj.panel_weather.addButton('button_wind_speed', 'average_wind_speed', 'info', 'wind= ', ' km/h', '', 1/4);
            obj.panel_weather.addButton('button_wind_dir', 'average_wind_dir', 'info', 'wind dir= ', ' deg', '', 1/4);
            obj.panel_weather.addButton('button_humidity', 'average_humidity', 'info', 'humid= ', '%', '', 1/4);
            obj.panel_weather.addButton('button_rain', 'any_rain', 'info', 'rain= ', '', '', 1/4);
            
            obj.panel_weather.margin = [0.01 0.01];
            obj.panel_weather.number = N;
            
            obj.panel_weather.make;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            N = 15;
            pos = pos - N;
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 pos/N_middle 0.8 N/N_middle]);
            
            obj.makeAxes;
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.0 0.0 0.1 0.05], obj.owner, '', 'custom', 'reset plot');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            
            obj.button_mean_only = GraphicButton(obj.panel_image, [0.0 0.95 0.1 0.05], obj.owner, 'checker.use_only_plot_mean', 'toggle', 'all', 'mean'); 
            
            obj.button_clicker = GraphicButton(obj.panel_image, [0.64 0.45 0.28 0.05], obj.owner, '', 'custom', '', ''); 
            
            %%%%%%%%%%% panel camera %%%%%%%%%%%%%%%%%
            
            obj.panel_camera = GraphicPanel(obj.owner, [0.2 0 0.8 2/N_middle], 'camera');
            obj.panel_camera.addButton('button_stop', 'commandCameraStop', 'push', 'Stop cam', '', '', 0.1, '', '', 'send a stop command to the camera-PC'); 
            obj.panel_camera.addButton('input_arguments', '', 'custom', ' ', '', '', 0.6, '', '', 'use this input to give arguments to the camera-PC'); 
            obj.panel_camera.addButton('button_start', 'commandCameraStart', 'push', 'Start cam', '', '', 0.1, '', '', 'send a start command to the camera-PC'); 
            obj.panel_camera.addButton('button_info', 'camera_info', 'info', ' ', '', '', 0.2, '', '', 'some information from the camera-PC'); 
            
            obj.panel_camera.margin = [0.01 0.1];
            obj.panel_camera.make;
            
            obj.panel_camera.input_arguments.control.Style = 'edit';
            obj.panel_camera.input_arguments.Callback = @obj.callback_arguments;
            
            %%%%%%%%%%% panel stop %%%%%%%%%%%%%%%%%%%
            
%             obj.panel_stop = GraphicPanel(obj.owner, [0.2 0 0.8 2/N_middle]);
            obj.panel_stop = GraphicPanel(obj.owner, [0.7 0.2 0.28 0.2]);
            obj.panel_stop.addButton('button_stop', 'stop', 'push', 'STOP MOUNT AND DOME');
%             obj.panel_stop.margin = [0.01 0.1];
            obj.panel_stop.make;
            
            %%%%%%%%%% Menus %%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.menu_options = MenuItem(obj, '&Options', 'menu');
            obj.menu_options.addButton('menu_dome', '&Dome', 'menu');
            obj.menu_options.menu_dome.addButton('button_east_num', '&East num= %d', 'input', 'dome.number_east', 'How many steps to move the East shutter');
            obj.menu_options.menu_dome.addButton('button_west_num', '&West num= %d', 'input', 'dome.number_west', 'How many steps to move the West shutter');
            
            obj.menu_options.addButton('menu_mount', '&Mount', 'menu');
            obj.menu_options.menu_mount.addButton('button_limit_alt', '&Alt limit= %d', 'input', 'mount.limit_alt', 'Altitude limit (degrees)');
            obj.menu_options.menu_mount.addButton('button_limit_flip', '&Flip limit= %d', 'input', 'mount.limit_flip', 'Meridian flip limit (degrees)');
            obj.menu_options.menu_mount.addButton('button_park', '&Park', 'push', 'mount.park', 'Park telescope');
            
            obj.menu_options.addButton('menu_arduino', '&Arduino', 'menu');
            obj.menu_options.menu_arduino.addButton('button_connect', '&Connect', 'push', 'mount.connectArduino', 'Reconnect bluetooth ScopeAssistant Arduino');
            obj.menu_options.menu_arduino.addButton('button_use_acc', '&Use Acelerometer', 'toggle', 'mount.use_accelerometer', 'Let the Arduino accelerometer stop telescope when it dips too low');
                        
            obj.menu_options.addButton('menu_object', '&Object', 'menu');
            obj.menu_options.menu_object.addButton('button_sync', '&Sync', 'push', 'mount.syncToTarget', 'Sync telescope to current position');
            obj.menu_options.menu_object.addButton('button_reset_prev', '&Reset history', 'push', 'mount.resetPrevObjects', 'Clear the list of previous targets'); 
            
            obj.menu_options.addButton('menu_timers', '&Timers', 'menu');
            obj.menu_options.menu_timers.addButton('button_run_t1', 'run t&1', 'push', 'callback_t1', 'Run the first timer, that updates all devices');
            obj.menu_options.menu_timers.addButton('button_run_t2', 'run t&2', 'push', 'callback_t2', 'Run the second timer, that checks weather and critical device errors (then shuts the dome if needed)');
            obj.menu_options.menu_timers.addButton('button_run_t3', 'run t&3', 'push', 'callback_t3', 'Run the third timer, that makes sure the other two are running');
            obj.menu_options.menu_timers.addButton('button_period1', 'Period1= %ds', 'input', 'period1', 'Adjust the period of the first timer');
            obj.menu_options.menu_timers.addButton('button_period2', 'Period2= %ds', 'input', 'period2', 'Adjust the period of the second timer');
            obj.menu_options.menu_timers.addButton('button_period3', 'Period3= %ds', 'input', 'period3', 'Adjust the period of the third timer');
            
            obj.menu_options.addButton('menu_weather', '&Weather', 'menu'); 
            obj.menu_options.menu_weather.addButton('button_wise', 'Use &Wise', 'toggle', 'checker.use_wise_data', 'Close when the Wise observatory deems it dangerous to operate'); 
            
            obj.menu_objects = MenuItem(obj, '&Device GUIs', 'menu');
            obj.menu_objects.addButton('button_dome', '&Dome', 'push', 'dome', 'dome GUI', 'Open the dome GUI');
            obj.menu_objects.addButton('button_mount', '&Mount', 'push', 'mount', 'mount GUI', 'Open the dome GUI');
            obj.menu_objects.addButton('button_scheduler', '&Scheduler', 'push', 'sched', 'scheduler GUI', 'Open the scheduler GUI');

            for ii = 1:length(obj.panels)
                obj.panels{ii}.panel.BackgroundColor = obj.color_bg; 
            end
            
            for ii = 1:length(obj.buttons)
                
                if strcmp(obj.buttons{ii}.control.Style, 'popupmenu') || util.text.cs(obj.buttons{ii}.type, 'input_text')
                    obj.buttons{ii}.control.BackgroundColor = obj.color_input;
                elseif strcmp(obj.buttons{ii}.type, 'info')
                    obj.buttons{ii}.BackgroundColor = obj.color_info;
                else
                    obj.buttons{ii}.BackgroundColor = obj.color_buttons;
                    obj.buttons{ii}.default_color = obj.color_buttons;
                end
                
            end
            
            obj.panel_camera.input_arguments.BackgroundColor = obj.color_input; 
            
            obj.panel_image.BackgroundColor = obj.color_bg;
            
            obj.update;
            
        end
            
        function makeAxes(obj, ~, ~)
            
            delete(obj.axes_image);
            
            obj.axes_image = axes('Parent', obj.panel_image);
            
            
        end
                
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            for ii = 1:length(obj.menus)
                obj.menus{ii}.update;
            end
            
%             default_color = [1 1 1].*0.94;
            
            % highlight malfunctioning devices! 
            if ~isempty(obj.owner.dome) && obj.owner.dome.status
                obj.panel_devices.button_dome.String = 'dome ok';
                obj.panel_devices.button_dome.BackgroundColor = obj.color_buttons;
            else
                obj.panel_devices.button_dome.String = 'dome error';
                obj.panel_devices.button_dome.BackgroundColor = 'red';
            end
            
            if ~isempty(obj.owner.mount) && obj.owner.mount.status
                obj.panel_devices.button_mount.String = 'mount ok';
                obj.panel_devices.button_mount.BackgroundColor = obj.color_buttons;
            else
                obj.panel_devices.button_mount.String = 'mount error';
                obj.panel_devices.button_mount.BackgroundColor = 'red';
            end
            
            if ~isempty(obj.owner.weather) && obj.owner.weather.status
                obj.panel_devices.button_weather.String = 'BoltWood ok';
                obj.panel_devices.button_weather.BackgroundColor = obj.color_buttons;
            else
                obj.panel_devices.button_weather.String = 'BoltWood error';
                obj.panel_devices.button_weather.BackgroundColor = 'red';
            end
            
            if ~isempty(obj.owner.wind) && obj.owner.wind.status
                obj.panel_devices.button_wind.String = 'WindETH ok';
                obj.panel_devices.button_wind.BackgroundColor = obj.color_buttons;
            else
                obj.panel_devices.button_wind.String = 'WindETH error';
                obj.panel_devices.button_wind.BackgroundColor = 'red';
            end
            
            % update telescope buttons
            if obj.owner.mount.telALT<20
                obj.panel_telescope.button_ALT.control.ForegroundColor = 'red';
            else
                obj.panel_telescope.button_ALT.control.ForegroundColor = 'black';
            end
            
            % update object buttons            
            if obj.owner.mount.objALT<20
                obj.panel_object.button_ALT.control.ForegroundColor = 'red';
            else
                obj.panel_object.button_ALT.control.ForegroundColor = 'black';
            end
            
            if isempty(obj.owner.mount.prev_objects)
                obj.panel_object.button_prev_objects.control.String = {' '};
            else
                obj.panel_object.button_prev_objects.control.String = obj.owner.mount.prev_objects;
            end
            
            if strcmp(obj.owner.mount.obj_pier_side, 'pierUnknown')
                obj.panel_object.button_pierside.ForegroundColor = 'red';
                obj.panel_object.button_pierside.Tooltip = 'Side of sky where object is right now (it is unobservable!)'; 
            elseif strcmp(obj.owner.mount.pier_side, obj.owner.mount.obj_pier_side)
                obj.panel_object.button_pierside.ForegroundColor = 'black';
                obj.panel_object.button_pierside.Tooltip = 'Side of sky where object is right now'; 
            else
                obj.panel_object.button_pierside.ForegroundColor = 'red';
                obj.panel_object.button_pierside.Tooltip = 'Side of sky where object is right now (need to flip!)'; 
            end
            
            obj.updateDomeStatusButtons;
            
            if obj.key_status_shift
                obj.panel_dome.button_close_west.String = 'Close West Full'; 
                obj.panel_dome.button_open_west.String = 'Open West Full'; 
                
                obj.panel_dome.button_close_east.String = 'Close East Full';                 
                obj.panel_dome.button_open_east.String = 'Open East Full'; 
            else
                obj.panel_dome.button_close_west.String = 'Close West'; 
                obj.panel_dome.button_open_west.String = 'Open West'; 

                obj.panel_dome.button_close_east.String = 'Close East';                 
                obj.panel_dome.button_open_east.String = 'Open East'; 
            end
            
            if obj.owner.checker.use_twilight_mode
                obj.panel_controls.button_twilight.String = 'twilight is on'; 
                obj.panel_controls.button_twilight.ForegroundColor = obj.panel_controls.button_twilight.color_on;
            else
                obj.panel_controls.button_twilight.String = 'twilight is off'; 
                obj.panel_controls.button_twilight.ForegroundColor = obj.panel_controls.button_twilight.color_off;
            end
            
            obj.panel_camera.input_arguments.String = obj.owner.camera_args;
            
            obj.updateStopButton;

            obj.owner.checker.plotWeather('ax', obj.axes_image, 'color', obj.color_bg);
                        
        end
        
        function updateDomeStatusButtons(obj)
            
            if obj.owner.dome.status==0
                obj.panel_dome.button_shutter_west.String = 'West: error';
                obj.panel_dome.button_shutter_east.String = 'East: error';
            else
                if obj.owner.dome.shutter_west_deg==0
                    obj.panel_dome.button_shutter_west.String = sprintf('West: open');
                elseif obj.owner.dome.shutter_west_deg==90
                    obj.panel_dome.button_shutter_west.String = sprintf('West: closed');
                else
                    obj.panel_dome.button_shutter_west.String = sprintf('West: %d deg', round(obj.owner.dome.shutter_west_deg));
                end
                
                if obj.owner.dome.shutter_east_deg==0
                    obj.panel_dome.button_shutter_east.String = sprintf('East: open');
                elseif obj.owner.dome.shutter_east_deg==90
                    obj.panel_dome.button_shutter_east.String = sprintf('East: closed');
                else
                    obj.panel_dome.button_shutter_east.String = sprintf('East: %d deg', round(obj.owner.dome.shutter_east_deg));
                end
                
            end
            
        end
        
        function updateStopButton(obj)
            
            if ~isempty(obj.owner.mount) && ~isempty(obj.owner.mount.is_slewing) && obj.owner.mount.is_slewing
                obj.panel_stop.button_stop.BackgroundColor = 'red';
            elseif obj.owner.dome.brake_bit==0
                obj.panel_stop.button_stop.BackgroundColor = 'red';
            else
                obj.panel_stop.button_stop.BackgroundColor = obj.color_buttons;
            end 
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_stop) && isvalid(obj.panel_stop.panel);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_twilight_mode(obj, ~, ~)
            
            if obj.debug_bit>1, disp('Callback: twilight mode'); end
            
            if obj.owner.checker.use_twilight_mode
                obj.owner.checker.use_twilight_mode = 0;
            else
                obj.owner.setup_t3; % first run this so it doesn't turn off twilight mode very shortly after user enables it
                obj.owner.checker.use_twilight_mode = 1;
                obj.owner.callback_t2; % update weather check with new light-level limit
            end
            
            obj.update;
            
        end
        
        function callback_key_press(obj, hndl, event)
            
            if any(contains(event.Modifier,'shift'))
                
                obj.key_status_shift = 1;
                
                obj.panel_dome.button_close_west.String = 'Close West Full'; 
                obj.panel_dome.button_open_west.String = 'Open West Full'; 
                
                obj.panel_dome.button_close_east.String = 'Close East Full';                 
                obj.panel_dome.button_open_east.String = 'Open East Full'; 
                                
            end
            
        end
        
        function callback_key_released(obj, hndl, event)
            
            obj.key_status_shift = 0;
            obj.panel_dome.button_close_west.String = 'Close West'; 
            obj.panel_dome.button_open_west.String = 'Open West'; 

            obj.panel_dome.button_close_east.String = 'Close East';                 
            obj.panel_dome.button_open_east.String = 'Open East'; 
            
        end
        
        function callback_prev_objects(obj, hndl, ~)
            
            if obj.debug_bit>1, disp('Callback: prev_objects'); end

            idx = hndl.Value;
            
            if ~isempty(idx) && ~isempty(hndl.String)
                str = hndl.String{idx};
                if ~isempty(strip(str))
                    obj.owner.mount.parseTargetString(str);
                end
            end
            
            obj.update;
            
        end
        
        function callback_close_west(obj, ~, ~)
            
            if obj.debug_bit>1, disp('callback: close west'); end
            
            if obj.key_status_shift
                obj.owner.dome.closeWestFull;
                obj.key_status_shift = 0;
            else
                obj.owner.dome.closeWest;
            end
            
            obj.update;
            
        end
        
        function callback_open_west(obj, ~, ~)

            if obj.debug_bit>1, disp('callback: open west'); end

            if obj.key_status_shift
                obj.owner.dome.openWestFull;
                obj.key_status_shift = 0;
            else
                obj.owner.dome.openWest;
            end
            
            obj.update;
            
        end
        
        function callback_close_east(obj, ~, ~)

            if obj.debug_bit>1, disp('callback: close east'); end

            if obj.key_status_shift
                obj.owner.dome.closeEastFull;
                obj.key_status_shift = 0;
            else
                obj.owner.dome.closeEast;
            end
            
            obj.update;
            
        end
        
        function callback_open_east(obj, ~, ~)

            if obj.debug_bit>1, disp('callback: open east'); end

            if obj.key_status_shift
                obj.owner.dome.openEastFull;
                obj.key_status_shift = 0;
            else
                obj.owner.dome.openEast;
            end
            
            obj.update;
            
        end
        
        function callback_interval_t1(obj, ~, ~)
            
            if obj.debug_bit>1, disp('callback: interval_t1'); end
            
            rep = util.text.inputdlg('Choose an interval for timer 1', obj.owner.period1);
            
            if ~isempty(rep)
                num = str2num(rep);

                if ~isempty(num) && num>0 && ~isequal(num, obj.owner.period1)
                    obj.owner.period1 = num;
                    obj.owner.setup_t1;
                end
            end
            
            obj.update;
            
        end
        
        function callback_interval_t2(obj, ~, ~)
            
            if obj.debug_bit>1, disp('callback: interval_t2'); end
            
            rep = util.text.inputdlg('Choose an interval for timer 2', obj.owner.period2);
            
            if ~isempty(rep)
                num = str2num(rep);

                if ~isempty(num) && num>0 && ~isequal(num, obj.owner.period2)
                    obj.owner.period2 = num;
                    obj.owner.setup_t2;
                end
            end
            
            obj.update;
            
        end
        
        function callback_interval_t3(obj, ~, ~)
            
            if obj.debug_bit>1, disp('callback: interval_t3'); end
            
            rep = util.text.inputdlg('Choose an interval for timer 3', obj.owner.period3);
            
            if ~isempty(rep)
                num = str2num(rep);

                if ~isempty(num) && num>0 && ~isequal(num, obj.owner.period3)
                    obj.owner.period3 = num;
                    obj.owner.setup_t3;
                end
            end
            
            obj.update;
            
        end
        
        function callback_slew(obj, ~, ~)
            
            if obj.debug_bit>1, disp('Callback: slew'); end
            
            if obj.owner.mount.check_need_flip
                res = questdlg('Need to flip for this target. Are you sure?', 'Flip needed!', 'Slew', 'Abort', 'Slew');
                if isempty(res) || strcmp(res, 'Abort')
                    return;
                end
            end
            
            obj.owner.mount.slew;
            
            obj.update;
            
        end
        
        function callback_arguments(obj, hndl, ~)
            
            if obj.debug_bit>1, disp('Callback: slew'); end
            
            args = hndl.String;
            
            if obj.debug_bit>1, disp(args); end
            
            obj.owner.camera_args = args;
            
            obj.update;
            
        end
        
        function callback_close(obj, ~, ~) % this is unused! 
           
            if obj.debug_bit>1, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end