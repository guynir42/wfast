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
        
        panel_shutter1;
        panel_shutter2;
        panel_shutter_both;
        
        panel_stop;
        
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
            N = 8;
            
            %%%%%%%%%%% panel status %%%%%%%%%%%%%%%%%
            
            obj.panel_status = GraphicPanel(obj.owner, [0 (N-1)/N 1 1/N], 'status');
            obj.panel_status.addButton('button_status', 'status', 'info', 'status= ', '', '', 1/3);
            obj.panel_status.addButton('button_reply', 'reply', 'info', 'rep= ', '', '', 1/3);
            obj.panel_status.addButton('button_connect', 'connect', 'push', '', '', '', 1/3);
            obj.panel_status.make;
            
            %%%%%%%%%%% panel shutters %%%%%%%%%%%%%%%
            
            obj.panel_shutter1 = GraphicPanel(obj.owner, [0 (N-3)/N 1 2/N], 'shutter1');
            obj.panel_shutter1.number = 2;
            obj.panel_shutter1.addButton('button_state', 'shutter1_deg', 'info', 'S1= ', ' deg', '', 1/3);
            obj.panel_shutter1.addButton('button_full_close', 'close1Full', 'push', 'closeFull', '', '', 1/3);
            obj.panel_shutter1.addButton('button_full_open', 'open1Full', 'push', 'openFull', '', '', 1/3);
            obj.panel_shutter1.addButton('button_number', 'number1', 'input', 'N= ', '', '', 1/3);
            obj.panel_shutter1.addButton('button_close', 'close1', 'push', 'close', '', '', 1/3);
            obj.panel_shutter1.addButton('button_open', 'open1', 'push', 'open', '', '', 1/3);
            
            obj.panel_shutter1.make;
            
            obj.panel_shutter2 = GraphicPanel(obj.owner, [0 (N-5)/N 1 2/N], 'shutter2');
            obj.panel_shutter2.number = 2;
            obj.panel_shutter2.addButton('button_state', 'shutter2_deg', 'info', 'S2= ', ' deg ', '', 1/3);
            obj.panel_shutter2.addButton('button_full_close', 'close2Full', 'push', 'closeFull', '', '', 1/3);
            obj.panel_shutter2.addButton('button_full_open', 'open2Full', 'push', 'openFull', '', '', 1/3);
            obj.panel_shutter2.addButton('button_number', 'number2', 'input', 'N= ', '', '', 1/3);
            obj.panel_shutter2.addButton('button_close', 'close2', 'push', 'close', '', '', 1/3);
            obj.panel_shutter2.addButton('button_open', 'open2', 'push', 'open', '', '', 1/3);
            
            obj.panel_shutter2.make;
            
            obj.panel_shutter_both = GraphicPanel(obj.owner, [0 (N-7)/N 1 2/N], 'both shutters');
            obj.panel_shutter_both.number = 2;
            obj.panel_shutter_both.addButton('button_state', '', 'custom', 'both: ', '', '', 1/3);
            obj.panel_shutter_both.addButton('button_full_close', 'closeBothFull', 'push', 'closeFull', '', '', 1/3);
            obj.panel_shutter_both.addButton('button_full_open', 'openBothFull', 'push', 'openFull', '', '', 1/3);
            obj.panel_shutter_both.addButton('button_number', 'number_both', 'input', 'N= ', '', '', 1/3);
            obj.panel_shutter_both.addButton('button_close', 'closeBoth', 'push', 'close', '', '', 1/3);
            obj.panel_shutter_both.addButton('button_open', 'openBoth', 'push', 'open', '', '', 1/3);
            
            obj.panel_shutter_both.make;
            
            %%%%%%%%%%% panel stop %%%%%%%%%%%%%%%%%%%
                        
            obj.panel_stop = GraphicPanel(obj.owner, [0 0 1 1/N]);
            obj.panel_stop.addButton('button_close', '', 'custom', 'CLOSE', '', '', 1/3);
            obj.panel_stop.addButton('button_stop', 'stop', 'push', 'STOP', '', '', 2/3);
            obj.panel_stop.make;
            obj.panel_stop.button_close.Callback = @obj.callback_close;
            
            obj.update;
            
        end
                        
        function update(obj,~,~)
            
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            obj.panel_shutter_both.button_state.String = [obj.owner.shutter1 '/' obj.owner.shutter2];
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_stop) && is_valid(obj.panel_stop);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end