classdef ManagerGUI < handle
    
    properties 
        
        owner@obs.Manager; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 16;
        edit_font_size = 14;
        small_font_size = 12;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
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
              
        version = 1.00;
        
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
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            N = 10;
            
            obj.panel_controls = GraphicPanel(obj.owner, [0 (N-9)/N 0.2 9/N], 'controls');
            obj.panel_controls.number = 9;
            obj.panel_controls.addButton('button_run_t1', '', 'custom', 'run t1', '', '', 0.5);
            obj.panel_controls.addButton('button_interval_t1', '', 'custom', ['P= ' num2str(obj.owner.checker.period1)], '', '', 0.5);
            obj.panel_controls.addButton('button_run_t2', '', 'custom', 'run t2', '', '', 0.5);
            obj.panel_controls.addButton('button_interval_t2', '', 'custom', ['P= ' num2str(obj.owner.checker.period2)], '', '', 0.5);
            obj.panel_controls.addButton('button_run_t3', '', 'custom', 'run t3', '', '', 0.5);
            obj.panel_controls.addButton('button_interval_t3', '', 'custom', ['P= ' num2str(obj.owner.checker.period3)], '', '', 0.5);
            
            obj.panel_controls.make;
            obj.panel_controls.button_run_t1.Callback = @obj.callback_run_t1;
            obj.panel_controls.button_interval_t1.Callback = @obj.callback_interval_t1;
            obj.panel_controls.button_run_t2.Callback = @obj.callback_run_t2;
            obj.panel_controls.button_interval_t2.Callback = @obj.callback_interval_t2;
            obj.panel_controls.button_run_t3.Callback = @obj.callback_run_t3;
            obj.panel_controls.button_interval_t3.Callback = @obj.callback_interval_t3;
            
            %%%%%%%%%%% panel objects %%%%%%%%%%%%%%%%
            
            obj.panel_objects = GraphicPanel(obj.owner, [0.8 (N-9)/N 0.2 9/N], 'objects');
            obj.panel_objects.number = 9;
%             obj.panel_objects.addButton('button_dome', 'dome', 'push', 'dome');
            obj.panel_objects.addButton('button_dome', 'dome', 'push', 'dome');
            obj.panel_objects.addButton('button_mount', 'mount', 'push', 'mount');
            obj.panel_objects.addButton('button_weather', 'weather', 'push', 'BoltWood');
            obj.panel_objects.addButton('button_wind', 'wind', 'push', 'WindETH');
            
            obj.panel_objects.make;
            
            %%%%%%%%%%% panel report %%%%%%%%%%%%%%%%%
            
            obj.panel_report = GraphicPanel(obj.owner, [0.2, (N-1)/N, 0.6, 1/N], '');
            obj.panel_report.addButton('button_report', 'report_string', 'info');
            obj.panel_report.make;
            
            %%%%%%%%%%% panel weather %%%%%%%%%%%%%%%%
            
            obj.panel_weather = GraphicPanel(obj.owner, [0.2, (N-3)/N, 0.6, 2/N], 'weather');
            obj.panel_weather.addButton('button_temp', 'average_temp', 'info', 'T= ', 'C', '', 1/3);
            obj.panel_weather.addButton('button_clouds', 'average_clouds', 'info', 'dT= ', 'C', '', 1/3);
            obj.panel_weather.addButton('button_light', 'average_light', 'info', 'L= ', '', '', 1/3);
            obj.panel_weather.addButton('button_wind', 'average_wind', 'info', 'wind= ', 'km/h', '', 1/3);
            obj.panel_weather.addButton('button_wind_az', 'average_wind_az', 'info', 'az= ', 'deg', '', 1/3);
            obj.panel_weather.addButton('button_hummid', 'average_humid', 'info', 'h= ', '%', '', 1/3);
            obj.panel_weather.number = 2;
            
            obj.panel_weather.make;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 1/N 0.6 (N-4)/N]);
            
            obj.makeAxes;
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.9 0.95 0.1 0.05], obj.owner, '', 'custom','reset');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            
            %%%%%%%%%%% panel stop %%%%%%%%%%%%%%%%%%%
            
            obj.panel_stop = GraphicPanel(obj.owner, [0.2 0 0.6 1/N]);
            obj.panel_stop.addButton('button_stop', 'stop', 'push', 'STOP');
            obj.panel_stop.make;
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Position', [0 0 0.2 1/N]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
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
            
            obj.panel_controls.button_interval_t1.String = ['P= ' num2str(obj.owner.checker.period1)];
            obj.panel_controls.button_interval_t2.String = ['P= ' num2str(obj.owner.checker.period2)];
            obj.panel_controls.button_interval_t3.String = ['P= ' num2str(obj.owner.checker.period3)];
            
            default_color = [1 1 1].*0.94;
            
            if ~isempty(obj.owner.dome) && obj.owner.dome.status
                obj.panel_objects.button_dome.String = 'dome: ok';
                obj.panel_objects.button_dome.BackgroundColor = default_color;
            else
                obj.panel_objects.button_dome.String = 'dome: error';
                obj.panel_objects.button_dome.BackgroundColor = 'red';
            end
            
            if ~isempty(obj.owner.mount) && obj.owner.mount.status
                obj.panel_objects.button_mount.String = 'mount: ok';
                obj.panel_objects.button_mount.BackgroundColor = default_color;
            else
                obj.panel_objects.button_mount.String = 'mount: error';
                obj.panel_objects.button_mount.BackgroundColor = 'red';
            end
            
            if ~isempty(obj.owner.weather) && obj.owner.weather.status
                obj.panel_objects.button_weather.String = 'BoltWood: ok';
                obj.panel_objects.button_weather.BackgroundColor = default_color;
            else
                obj.panel_objects.button_weather.String = 'BoltWood: error';
                obj.panel_objects.button_weather.BackgroundColor = 'red';
            end
            
            if ~isempty(obj.owner.wind) && obj.owner.wind.status
                obj.panel_objects.button_wind.String = 'WindETH: ok';
                obj.panel_objects.button_weather.BackgroundColor = default_color;
            else
                obj.panel_objects.button_wind.String = 'WindETH: error';
                obj.panel_objects.button_wind.BackgroundColor = 'red';
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
            
            obj.owner.checker.callback_t1;
            
            obj.update;
            
        end
        
        function callback_interval_t1(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: interval_t1'); end
            
            rep = util.text.inputdlg('Choose an interval for timer 1', obj.owner.checker.period1);
            
            num = str2num(rep);
            
            if ~isempty(num) && num>0 && ~isequal(num, obj.owner.checker.period1)
                obj.owner.checker.period1 = num;
                obj.owner.checker.setup_t1;
            end
            
            obj.update;
            
        end
        
        
        function callback_run_t2(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: run_t2'); end
            
            obj.owner.checker.callback_t2;
            
            obj.update;
            
        end
        
        function callback_interval_t2(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: interval_t2'); end
            
            rep = util.text.inputdlg('Choose an interval for timer 2', obj.owner.checker.period2);
            
            num = str2num(rep);
            
            if ~isempty(num) && num>0 && ~isequal(num, obj.owner.checker.period2)
                obj.owner.checker.period2 = num;
                obj.owner.checker.setup_t2;
            end
            
            obj.update;
            
        end
        
        function callback_run_t3(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: run_t3'); end
            
            obj.owner.checker.callback_t3;
            
            obj.update;
            
        end
        
        function callback_interval_t3(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: interval_t3'); end
            
            rep = util.text.inputdlg('Choose an interval for timer 3', obj.owner.checker.period3);
            
            num = str2num(rep);
            
            if ~isempty(num) && num>0 && ~isequal(num, obj.owner.checker.period3)
                obj.owner.checker.period3 = num;
                obj.owner.checker.setup_t3;
            end
            
            obj.update;
            
        end
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end