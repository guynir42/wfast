classdef ManagerGUI < handle
    
    properties 
        
        owner@obs.Manager; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 13;
        big_font_size = 16;
        edit_font_size = 11;
        small_font_size = 9;
        
        color_on = [0 0.3 1];
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_telescope;
        
        panel_dome;
        
        panel_controls;
        
        panel_report;
        
        panel_weather;
        
        panel_objects;
        
        panel_stop;
        
        panel_close;
        button_close;
        
        panel_image;
        button_reset_axes;
        axes_image;
    
    end
    
    properties (Hidden=true)
              
        version = 1.01;
        
    end
            
    methods % constructor
       
        function obj = ManagerGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'obs.Manager')
                
                if obj.debug_bit, fprintf('ManagerGUI constructor v%4.2f\n', obj.version); end
                
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
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            
            obj.fig = util.plot.FigHandler('Observatory Manager');
            obj.fig.clear;
            movegui(obj.fig.fig, 'center');
            
            obj.fig.bottom = -5;
            obj.fig.height = 25;
            obj.fig.width = 36;
            
            
            N_left = 20;
            pos = N_left;
            
            %%%%%%%%%%% panel telescope %%%%%%%%%%%%%%%
            
            N = 5;
            pos = pos - N;
            obj.panel_telescope = GraphicPanel(obj.owner, [0 pos/N_left 0.2 N/N_left], 'telescope');
            obj.panel_telescope.number = N;
            obj.panel_telescope.addButton('button_RA', 'RA', 'info', 'RA: ', '', '', 0.5);
            obj.panel_telescope.addButton('button_DE', 'DEC', 'info', 'DE: ', '', '', 0.5);            
            obj.panel_telescope.addButton('button_LST', 'LST', 'info', 'LST: ', '', '', 0.5);
            obj.panel_telescope.addButton('button_ALT', 'ALT', 'info', 'ALT: ', ' deg', '', 0.5);
            obj.panel_telescope.addButton('button_tracking', 'tracking', 'toggle', 'tracking off', 'tracking on', '', 0.5, obj.color_on, 'red');
            obj.panel_telescope.margin = [0.02 0.01];
            obj.panel_telescope.make;
            
            obj.panel_telescope.button_RA.Tooltip = 'Current right ascention of mount';
            obj.panel_telescope.button_DE.Tooltip = 'Current declination of mount';
            obj.panel_telescope.button_LST.Tooltip = 'Local Sidereal Time';
            obj.panel_telescope.button_ALT.Tooltip = 'Current altitidue of mount';
            obj.panel_telescope.button_tracking.Tooltip = 'Turn telescope tracking on and off';
            
            %%%%%%%%%%% panel dome %%%%%%%%%%%%%%%
            
            N = 4;
            pos = pos - N;
            obj.panel_dome = GraphicPanel(obj.owner, [0 pos/N_left 0.2 N/N_left], 'dome');
            obj.panel_dome.number = N;            
            obj.panel_dome.addButton('button_close_dome', 'closeDome', 'push', 'Close Dome');
            obj.panel_dome.addButton('button_shutter_west', '', 'custom', 'West shutter: ', '', 'small', 0.5);
            obj.panel_dome.addButton('button_shutter_east', '', 'custom', 'East shutter: ', '', 'small', 0.5);
            obj.panel_dome.margin = [0.02 0.01];
            obj.panel_dome.make;
            
            obj.panel_dome.button_close_dome.Tooltip = 'Immediately close both shutters';
            obj.panel_dome.button_shutter_west.Tooltip = 'Position of West shutter';
            obj.panel_dome.button_shutter_east.Tooltip = 'Position of East shutter';
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            N = 9;
            pos = pos - N;
            obj.panel_controls = GraphicPanel(obj.owner, [0 pos/N_left 0.2 N/N_left], 'controls');
            obj.panel_controls.number = N;
            obj.panel_controls.addButton('button_autoshutdown', 'use_shutdown', 'toggle', 'auto shut down disabled', 'auto shutdown enabled', ...
                '', 1, obj.color_on, 'red'); 
            obj.panel_controls.addButton('button_run_t1', '', 'custom', 'run t1', '', '', 0.5);
            obj.panel_controls.addButton('button_interval_t1', '', 'custom', ['P= ' num2str(obj.owner.period1) 's'], '', '', 0.5);
            obj.panel_controls.addButton('button_run_t2', '', 'custom', 'run t2', '', '', 0.5);
            obj.panel_controls.addButton('button_interval_t2', '', 'custom', ['P= ' num2str(obj.owner.period2) 's'], '', '', 0.5);
            obj.panel_controls.addButton('button_run_t3', '', 'custom', 'run t3', '', '', 0.5);
            obj.panel_controls.addButton('button_interval_t3', '', 'custom', ['P= ' num2str(obj.owner.period3) 's'], '', '', 0.5);
            obj.panel_controls.margin = [0.01 0.01];
            obj.panel_controls.make;
            
            obj.panel_controls.button_autoshutdown.Tooltip = 'Allow manager to close dome and stop tracking if weather is bad or if there is a device failure';
            
            obj.panel_controls.button_run_t1.Callback = @obj.callback_run_t1;
            obj.panel_controls.button_run_t1.Tooltip = 'Immediately run the timer 1 callback (update sensors)';
            
            obj.panel_controls.button_interval_t1.Callback = @obj.callback_interval_t1;
            obj.panel_controls.button_interval_t1.Tooltip = 'Change the interval of timer 1';
            
            obj.panel_controls.button_run_t2.Callback = @obj.callback_run_t2;
            obj.panel_controls.button_run_t2.Tooltip = 'Immediately run the time 2 callback (check observatory and shutdown if needed)';
            
            obj.panel_controls.button_interval_t2.Callback = @obj.callback_interval_t2;
            obj.panel_controls.button_interval_t2.Tooltip = 'Change the interval of timer 2';
            
            obj.panel_controls.button_run_t3.Callback = @obj.callback_run_t3;
            obj.panel_controls.button_run_t3.Tooltip = 'Immediately run the timer 3 callback (verify other timers are still running)';
            
            obj.panel_controls.button_interval_t3.Callback = @obj.callback_interval_t3;
            obj.panel_controls.button_interval_t3.Tooltip = 'Change the interval of timer 3';
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Position', [0 0 0.2 2/N_left]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1]+[0.05 0.1 -0.1 -0.2], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            %%%%%%%%%%% panel objects %%%%%%%%%%%%%%%%
            
            N_right = 20;
            pos = N_right;
            
            N = 20;
            pos = pos - N;
            
            obj.panel_objects = GraphicPanel(obj.owner, [0.8 pos/N_right 0.2 N/N_right], 'objects');
            obj.panel_objects.number = N;
%             obj.panel_objects.addButton('button_dome', 'dome', 'push', 'dome');
            obj.panel_objects.addButton('button_dome', 'dome', 'push', 'dome');
            obj.panel_objects.addButton('button_mount', 'mount', 'push', 'mount');
            obj.panel_objects.addButton('button_weather', 'weather', 'push', 'BoltWood');
            obj.panel_objects.addButton('button_wind', 'wind', 'push', 'WindETH');
            obj.panel_objects.margin = [0.01 0.005];
            obj.panel_objects.make;
            
            %%%%%%%%%%% panel report %%%%%%%%%%%%%%%%%
            
            N_middle = 20;
            pos = N_middle;
            
            N = 1;
            pos = pos - N;
            
            obj.panel_report = GraphicPanel(obj.owner, [0.2, pos/N_middle, 0.6, N/N_middle], 'report');
            obj.panel_report.addButton('button_report', 'report', 'info');
            obj.panel_report.make;
            
            %%%%%%%%%%% panel weather %%%%%%%%%%%%%%%%
            
            N = 2;
            pos = pos - N;
            
            obj.panel_weather = GraphicPanel(obj.owner, [0.2, pos/N_middle, 0.6, N/N_middle], 'weather');
            obj.panel_weather.addButton('button_temp', 'average_temp', 'info', 'Amb. Temp= ', 'C', '', 1/3);
            obj.panel_weather.addButton('button_clouds', 'average_clouds', 'info', 'dT= ', 'C', '', 1/3);
            obj.panel_weather.addButton('button_light', 'average_light', 'info', 'Light= ', '', '', 1/3);
            obj.panel_weather.addButton('button_wind', 'average_wind', 'info', 'wind= ', ' km/h', '', 1/3);
            obj.panel_weather.addButton('button_wind_az', 'average_wind_az', 'info', 'wind az= ', ' deg', '', 1/3);
            obj.panel_weather.addButton('button_hummid', 'average_humid', 'info', 'humidity= ', '%', '', 1/3);
            obj.panel_weather.margin = [0.01 0.01];
            obj.panel_weather.number = N;
            
            obj.panel_weather.make;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            N = 15;
            pos = pos - N;
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 pos/N_middle 0.6 N/N_middle]);
            
            obj.makeAxes;
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.9 0.95 0.1 0.05], obj.owner, '', 'custom','reset');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            
            %%%%%%%%%%% panel stop %%%%%%%%%%%%%%%%%%%
            
            obj.panel_stop = GraphicPanel(obj.owner, [0.2 0 0.6 2/N_middle]);
            obj.panel_stop.addButton('button_stop', 'stop', 'push', 'STOP');
            obj.panel_stop.margin = [0.01 0.1];
            obj.panel_stop.make;
            
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
            
            if obj.owner.mount.telRA_deg<obj.owner.mount.LST_deg && strcmp(obj.owner.mount.hndl.SideOfPier, 'pierWest')
                obj.panel_telescope.button_RA.BackgroundColor = 'red';
            else
                obj.panel_telescope.button_RA.BackgroundColor = util.plot.GraphicButton.defaultColor;
            end
            
            if obj.owner.mount.telALT<20
                obj.panel_telescope.button_ALT.control.ForegroundColor = 'red';
            else
                obj.panel_telescope.button_ALT.control.ForegroundColor = 'black';
            end
            
            if obj.owner.dome.status==0
                obj.panel_dome.button_shutter_west.String = 'Shut.West: error';
                obj.panel_dome.button_shutter_east.String = 'Shut.East: error';
            else
                if obj.owner.dome.shutter_west_deg==0
                    obj.panel_dome.button_shutter_west.String = sprintf('Shut.West: open');
                elseif obj.owner.dome.shutter_west_deg==90
                    obj.panel_dome.button_shutter_west.String = sprintf('Shut.West: closed');
                else
                    obj.panel_dome.button_shutter_west.String = sprintf('Shut.West: %d deg', round(obj.owner.dome.shutter_west_deg));
                end
                
                if obj.owner.dome.shutter_east_deg==0
                    obj.panel_dome.button_shutter_east.String = sprintf('Shut.East: open');
                elseif obj.owner.dome.shutter_east_deg==90
                    obj.panel_dome.button_shutter_east.String = sprintf('Shut.East: closed');
                else
                    obj.panel_dome.button_shutter_east.String = sprintf('Shut.East: %d deg', round(obj.owner.dome.shutter_east_deg));
                end
                
            end
            
            obj.panel_controls.button_interval_t1.String = ['P= ' num2str(obj.owner.period1) 's'];
            obj.panel_controls.button_interval_t2.String = ['P= ' num2str(obj.owner.period2) 's'];
            obj.panel_controls.button_interval_t3.String = ['P= ' num2str(obj.owner.period3) 's'];
            
            default_color = [1 1 1].*0.94;
            
            if ~isempty(obj.owner.dome) && obj.owner.dome.status
                obj.panel_objects.button_dome.String = 'dome GUI (status ok)';
                obj.panel_objects.button_dome.BackgroundColor = default_color;
            else
                obj.panel_objects.button_dome.String = 'dome GUI (status: error)';
                obj.panel_objects.button_dome.BackgroundColor = 'red';
            end
            
            if ~isempty(obj.owner.mount) && obj.owner.mount.status
                obj.panel_objects.button_mount.String = 'mount GUI (status: ok)';
                obj.panel_objects.button_mount.BackgroundColor = default_color;
            else
                obj.panel_objects.button_mount.String = 'mount GUI (status: error)';
                obj.panel_objects.button_mount.BackgroundColor = 'red';
            end
            
            if ~isempty(obj.owner.weather) && obj.owner.weather.status
                obj.panel_objects.button_weather.String = 'BoltWood (status: ok)';
                obj.panel_objects.button_weather.BackgroundColor = default_color;
            else
                obj.panel_objects.button_weather.String = 'BoltWood (status: error)';
                obj.panel_objects.button_weather.BackgroundColor = 'red';
            end
            
            if ~isempty(obj.owner.wind) && obj.owner.wind.status
                obj.panel_objects.button_wind.String = 'WindETH (status: ok)';
                obj.panel_objects.button_wind.BackgroundColor = default_color;
            else
                obj.panel_objects.button_wind.String = 'WindETH (status: error)';
                obj.panel_objects.button_wind.BackgroundColor = 'red';
            end
            
            timer_vec = obj.owner.areTimersRunning;
            
            for ii = 1:3
                button = obj.panel_controls.(['button_run_t' num2str(ii)]);
                button2 = obj.panel_controls.(['button_interval_t' num2str(ii)]);
                if timer_vec(ii)
                    button.String = ['run t' num2str(ii)];
                    button.BackgroundColor = util.plot.GraphicButton.defaultColor;
                    button2.BackgroundColor = util.plot.GraphicButton.defaultColor;
                else
                    button.String = ['t' num2str(ii) ' stopped'];
                    button.BackgroundColor = 'red';
                    button2.BackgroundColor = 'red';
                end
            end
            
            obj.owner.checker.plotWeather('ax', obj.axes_image);
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_run_t1(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: run_t1'); end
            
            obj.owner.callback_t1;
            
            obj.update;
            
        end
        
        function callback_interval_t1(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: interval_t1'); end
            
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
        
        
        function callback_run_t2(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: run_t2'); end
            
            obj.owner.callback_t2;
            
            obj.update;
            
        end
        
        function callback_interval_t2(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: interval_t2'); end
            
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
        
        function callback_run_t3(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: run_t3'); end
            
            obj.owner.callback_t3;
            
            obj.update;
            
        end
        
        function callback_interval_t3(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: interval_t3'); end
            
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
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end