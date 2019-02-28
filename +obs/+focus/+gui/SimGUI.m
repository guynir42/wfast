classdef SimGUI < handle
    
    properties 
        
        owner@obs.focus.Simulator; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 16;
        edit_font_size = 14;
        small_font_size = 12;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_position;
        panel_defocus;
        
        panel_close;
        button_close;
        
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = SimGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'obs.focus.Simulator')
                
                if obj.debug_bit, fprintf('SimGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input an obs.focus.Simulator to constructor of SimGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            
            obj.buttons = {};
            
            obj.fig = util.plot.FigHandler('sim focuser');
            
            obj.fig.reset;
            obj.fig.bottom = 5;
            obj.fig.height = 8;
            obj.fig.width = 4;
            movegui(obj.fig.fig, 'center');
            
            %%%%%%%%%%% panel position %%%%%%%%%%%%%%%
            
            obj.panel_position = GraphicPanel(obj.owner, [0 0.6 1 0.4], 'position');
            obj.panel_position.addButton('button_pos', 'pos', 'input', 'pos= ');
            obj.panel_position.addButton('button_best', 'best_pos', 'input', 'best= ');
            obj.panel_position.make;
            
            %%%%%%%%%%% panel defocus %%%%%%%%%%%%%%%%
            
            obj.panel_position = GraphicPanel(obj.owner, [0 0.2 1 0.4], 'defocus');
            obj.panel_position.addButton('button_parameter', 'defocus_parameter', 'input', 'par= ');
            obj.panel_position.addButton('button_inner', 'inner_annulus', 'input', 'inner= ');
            obj.panel_position.make;
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
                        
            obj.panel_close = uipanel('Position', [0 0 1 0.2]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE');
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