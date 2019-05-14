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
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 25;
            movegui(obj.fig.fig, 'center');
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            N = 10;
            
            obj.panel_controls = GraphicPanel(obj.owner, [0 (N-9)/N 0.2 9/N], 'controls');
            obj.panel_controls.number = 9;
            obj.panel_controls.make;
            
            %%%%%%%%%%% panel objects %%%%%%%%%%%%%%%%
            
            obj.panel_objects = GraphicPanel(obj.owner, [0.8 (N-9)/N 0.2 9/N], 'objects');
            obj.panel_objects.number = 9;
%             obj.panel_objects.addButton('button_dome', 'dome', 'push', 'dome');
            obj.panel_objects.addButton('button_dome', 'dome', 'push', 'dome');
            obj.panel_objects.addButton('button_mount', 'mount', 'push', 'mount');
            obj.panel_objects.make;
            
            
            %%%%%%%%%%% panel report %%%%%%%%%%%%%%%%%
            
            obj.panel_report = GraphicPanel(obj.owner, [0.2, (N-1)/N, 0.6, 1/N], '');
            obj.panel_report.addButton('button_report', 'report_string', 'info');
            obj.panel_report.make;
            
            %%%%%%%%%%% panel weather %%%%%%%%%%%%%%%%
            
            obj.panel_weather = GraphicPanel(obj.owner, [0.2, (N-3)/N, 0.6, 2/N], 'weather');
            obj.panel_weather.number = 2;
            
            obj.panel_weather.make;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 1/N 0.6 (N-4)/N]);
            
            obj.makeAxes;
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.9 0.95 0.1 0.05], obj.owner, '', 'custom','reset');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            
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
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end