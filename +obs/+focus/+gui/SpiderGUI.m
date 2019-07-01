classdef SpiderGUI < handle
    
    properties 
        
        owner@obs.focus.FocusSpider; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 12;
        big_font_size = 14;
        edit_font_size = 11;
        small_font_size = 10;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_status;
        
        panel_pos;
        panel_tip;
        panel_tilt;
        
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
       
        function obj = SpiderGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'obs.focus.FocusSpider')
                
                if obj.debug_bit, fprintf('SpiderGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input an obs.focus.FocusSpider to constructor of SpiderGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            
            obj.fig = util.plot.FigHandler('focus spider');
            
            obj.fig.bottom = 5;
            obj.fig.height = 10;
            obj.fig.width = 8;
            movegui(obj.fig.fig, 'center');
            obj.fig.left = 20;
            
            N = 7;
            pos = N;
            
            %%%%%%%%%%% panel status %%%%%%%%%%%%%%%%%%
            
            pos = pos-2;
            obj.panel_status = GraphicPanel(obj.owner, [0 pos/N 1 2/N], 'connection', 1);
            obj.panel_status.addButton('status', 'status', 'info', 'status= ', '', '', 0.7);
            obj.panel_status.addButton('connect', 'connect', 'custom', 'CONNECT', '', '', 0.3);
            obj.panel_status.addButton('vector', 'pos_vec', 'info', 'pos= ');
            obj.panel_status.margin = [0.01 0.05];
            obj.panel_status.make;
            obj.panel_status.connect.Callback = @obj.callback_connect;
            
            %%%%%%%%%%% panel pos %%%%%%%%%%%%%%%%%%
            
            pos = pos-4;
            obj.panel_pos = GraphicPanel(obj.owner, [0/3 pos/N 1/3 4/N], 'position');
            obj.panel_pos.addButton('value', 'pos', 'input', ' ');
            obj.panel_pos.addButton('up', 'up', 'push', 'UP');
            obj.panel_pos.addButton('step', 'step', 'input', ' ');
            obj.panel_pos.addButton('down', 'down', 'push', 'DOWN');
            obj.panel_pos.margin = [0.1 0.03];
            obj.panel_pos.make;
            
            %%%%%%%%%%% panel tip %%%%%%%%%%%%%%%%%%
            
            obj.panel_tip = GraphicPanel(obj.owner, [1/3 pos/N 1/3 4/N], 'tip');
            obj.panel_tip.addButton('value', 'tip', 'input', ' ');
            obj.panel_tip.addButton('up', 'tip_up', 'push', 'UP');
            obj.panel_tip.addButton('step', 'step_tip', 'input', ' ');
            obj.panel_tip.addButton('down', 'tip_down', 'push', 'DOWN');
            obj.panel_tip.margin = [0.1 0.03];
            obj.panel_tip.make;
            
            %%%%%%%%%%% panel tilt %%%%%%%%%%%%%%%%%%
            
            obj.panel_tilt = GraphicPanel(obj.owner, [2/3 pos/N 1/3 4/N], 'tilt');
            obj.panel_tilt.addButton('value', 'tilt', 'input', ' ');
            obj.panel_tilt.addButton('up', 'tilt_up', 'push', 'UP');
            obj.panel_tilt.addButton('step', 'step_tilt', 'input', ' ');
            obj.panel_tilt.addButton('down', 'tilt_down', 'push', 'DOWN');
            obj.panel_tilt.margin = [0.1 0.03];
            obj.panel_tilt.make;
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            pos = pos - 1;
            obj.panel_close = uipanel('Position', [0 pos 1 1/N]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.font_size = 'big';
            obj.button_close.Callback = @obj.callback_close;
            
            obj.update;
            
        end
        
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            if obj.owner.status
                obj.panel_status.connect.String = 'Disconnect';
            else
                obj.panel_status.connect.String = 'connect';
            end
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_connect(obj, ~, ~)
           
            if obj.owner.status
                if obj.debug_bit, disp('callback: disconnect'); end
                obj.owner.disconnect;
            else
                if obj.debug_bit, disp('callback: connect'); end
                obj.owner.connect;
            end
            
            obj.update;
            
        end
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end