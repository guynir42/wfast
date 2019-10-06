classdef objGUI < handle
    
    properties 
        
        owner@package.OBJ; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        menus = {};
        
        latest_error = '';
        latest_warning = '';
        
        font_size = 12;
        big_font_size = 16;
        edit_font_size = 11;
        small_font_size = 10;
        
        color_on = [0 0.3 1];
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_controls;
        
        panel_contrast;
    
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
       
        function obj = objGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'package.OBJ')
                
                if obj.debug_bit, fprintf('objGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            elseif isa(owner, '...')
                
                
                
            else
                error('Input an ... to constructor of ...!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            obj.menus = {};
            
            obj.fig = util.plot.FigHandler('...');
            obj.fig.clear;
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 25;
            movegui(obj.fig.fig, 'center');
            
            
            %%%%%%%%%%%%%%%%%%%%%%% MENUS %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % MenuItem(parent, text, type, variable, tooltip, separator)
            % obj.addButton(text, type, variable, tooltip, separator)
            % menu types: menu, toggle, push, input, input_text, info, custom
            
            
            
            N_left = 10; % number of buttons on left side
            
            %%%%%%%%%%%%%%%%%%% LEFT SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            pos = N_left;
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            % Add buttons using obj.addButton(button_name, var_name='', type='', str1='', str2='', font_size='', split=1, color_on=[], color_off=[], tooltip)
            
            num_buttons = 5;
            pos = pos-num_buttons;
            obj.panel_controls = GraphicPanel(obj.owner, [0 pos/N_left 0.2 num_buttons/N_left], 'controls', 1); % last input is for vertical (default)
            obj.panel_controls.number = num_buttons;
            
            obj.panel_controls.make;
            
            %%%%%%%%%%% panel contrast %%%%%%%%%%%%%%%
            
            obj.panel_contrast = util.plot.ContrastLimits(obj.axes_image, obj.fig.fig, [0 pos/N_left 0.2 5/N_left], 1); % last input is for vertical (default)
            obj.panel_contrast.font_size = obj.font_size;
            obj.panel_contrast.big_font_size = obj.big_font_size;
            obj.panel_contrast.small_font_size = obj.small_font_size;
            obj.panel_contrast.edit_font_size = obj.edit_font_size;
            
            %%%%%%%%%%%%%%%%%%% RIGHT SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_right = 10; pos = N_right;
                        
            %%%%%%%%%%%%%%%%%%%%%% MIDDLE %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_middle = 10; pos = N_middle;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            N = 9; pos = pos - N;
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 pos/N_middle 0.6 N/N_middle]);
                        
            obj.makeAxes;
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.85 0.95 0.15 0.05], obj.owner, '', 'custom', 'new axes', '');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            obj.button_reset_axes.Tooltip = 'Create a new image axis, zoomed out and with default contrast limits'; 
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            N = 1; pos = pos - N;            
            obj.panel_close = uipanel('Position', [0 pos 0.2 N/N_middle]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            obj.update;
            
        end
            
        function makeAxes(obj, ~, ~)
            
            delete(obj.axes_image);
            
            obj.axes_image = axes('Parent', obj.panel_image);
            
            obj.panel_contrast.ax = obj.axes_image;
%             colorbar(obj.axes_image);
            axis(obj.axes_image, 'image');
            
            obj.panel_contrast.ax = obj.axes_image;
            
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
            
            obj.panel_contrast.update;
            
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