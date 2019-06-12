classdef AstroHavenGUI < handle
    
    properties 
        
        owner@obs.dome.AstroHaven; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 16;
        edit_font_size = 14;
        small_font_size = 12;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_status;
        
        panel_shutter_west;
        panel_shutter_east;
        panel_shutter_both;
        
        panel_stop1;
        panel_stop2;
        
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = AstroHavenGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'obs.dome.AstroHaven')
                
                if obj.debug_bit, fprintf('AstroHavenGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input an AstroHaven to constructor of AstroHavenGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            
            obj.fig = util.plot.FigHandler('AstroHaven dome');
            obj.fig.clear;
            obj.fig.height = 16;
            obj.fig.width = 12;
            movegui(obj.fig.fig, 'center');
            obj.fig.left = 24;
            obj.fig.bottom = -5;
            
            N = 10; % total number of rows
            pos = N;
            
            %%%%%%%%%%% panel status %%%%%%%%%%%%%%%%%
            
            pos = pos - 3;
            
            obj.panel_status = GraphicPanel(obj.owner, [0 pos/N 1 3/N], 'status');
            obj.panel_status.addButton('button_close', '', 'custom', 'close GUI', '', '', 1/4);
            obj.panel_status.addButton('button_status', 'status', 'info', 'status= ', '', '', 1/4);
            obj.panel_status.addButton('button_reply', 'reply', 'info', 'reply= ', '', '', 1/4);
            obj.panel_status.addButton('button_connect', 'connect', 'push', '', '', '', 1/4);
            obj.panel_status.addButton('button_west_str', 'shutter_west', 'info', 'West: ', '', '', 1/2);
            obj.panel_status.addButton('button_east_str', 'shutter_east', 'info', 'East: ', '', '', 1/2);
            obj.panel_status.addButton('button_west_deg', 'shutter_west_deg', 'info', 'West: ', ' deg', '', 1/2);
            obj.panel_status.addButton('button_east_deg', 'shutter_east_deg', 'info', 'East: ', ' deg', '', 1/2);
            obj.panel_status.make;
            
            obj.panel_status.button_close.Callback = @obj.callback_close;
            
            %%%%%%%%%%% panel stop2 %%%%%%%%%%%%%%%%%%%
            
            pos = pos - 1;
            
            obj.panel_stop1 = GraphicPanel(obj.owner, [0 pos/N 1 1/N]);
            obj.panel_stop1.addButton('button_stop', 'stop', 'push', 'STOP!');
            obj.panel_stop1.make;
            
            %%%%%%%%%%% panel shutters %%%%%%%%%%%%%%%
            
            pos = pos - 5;
            
            obj.panel_shutter_west = GraphicPanel(obj.owner, [0/3 pos/N 1/3 5/N], 'West');
            obj.panel_shutter_west.addButton('button_full_close', 'closeWestFull', 'push', 'closeWestFull');
            obj.panel_shutter_west.addButton('button_close', 'closeWest', 'push', 'closeWest');
            obj.panel_shutter_west.addButton('button_number', 'number_west', 'input', 'N= ');
            obj.panel_shutter_west.addButton('button_open', 'openWest', 'push', 'openWest');
            obj.panel_shutter_west.addButton('button_full_open', 'openWestFull', 'push', 'openWestFull');
            obj.panel_shutter_west.make;
            
            obj.panel_shutter_both = GraphicPanel(obj.owner, [1/3 pos/N 1/3 5/N], 'Both');
            obj.panel_shutter_both.addButton('button_full_close', 'closeBothFull', 'push', 'closeBothFull');
            obj.panel_shutter_both.addButton('button_close', 'closeBoth', 'push', 'closeBoth');
            obj.panel_shutter_both.addButton('button_number', 'number_both', 'input', 'N= ');
            obj.panel_shutter_both.addButton('button_open', 'openBoth', 'push', 'openBoth');
            obj.panel_shutter_both.addButton('button_full_open', 'openBothFull', 'push', 'openBothFull');
            obj.panel_shutter_both.make;
            
            obj.panel_shutter_east = GraphicPanel(obj.owner, [2/3 pos/N 1/3 5/N], 'East');
            obj.panel_shutter_east.addButton('button_full_close', 'closeEastFull', 'push', 'closeEastFull');
            obj.panel_shutter_east.addButton('button_close', 'closeEast', 'push', 'closeEast');
            obj.panel_shutter_east.addButton('button_number', 'number_east', 'input', 'N= ');
            obj.panel_shutter_east.addButton('button_open', 'openEast', 'push', 'openEast');
            obj.panel_shutter_east.addButton('button_full_open', 'openEastFull', 'push', 'openEastFull', '', '', 1, 'blue');
            obj.panel_shutter_east.make;
            
            %%%%%%%%%%% panel stop2 %%%%%%%%%%%%%%%%%%%
            
            obj.panel_stop2 = GraphicPanel(obj.owner, [0 0 1 1/N]);
            obj.panel_stop2.addButton('button_stop', 'stop', 'push', 'STOP!');
            obj.panel_stop2.make;
            
            obj.update;
            
        end
                        
        function update(obj,~,~)
            
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            if obj.owner.brake_bit
                obj.panel_stop1.button_stop.BackgroundColor = util.plot.GraphicButton.defaultColor;
                obj.panel_stop2.button_stop.BackgroundColor = util.plot.GraphicButton.defaultColor;
            else
                obj.panel_stop1.button_stop.BackgroundColor = 'red';
                obj.panel_stop2.button_stop.BackgroundColor = 'red';
            end
            
            if obj.owner.status
                obj.panel_status.button_status.BackgroundColor = util.plot.GraphicButton.defaultColor;
            else
                obj.panel_status.button_status.BackgroundColor = 'red';
            end
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_stop2) && is_valid(obj.panel_stop2);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end