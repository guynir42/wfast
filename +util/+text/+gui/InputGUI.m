classdef InputGUI < handle
    
    properties 
        
        owner@util.text.InputVars; % link back to containg object

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
        
        panel_inputs; 
        
        panel_close;
        button_close;
        
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = InputGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'util.text.InputVars')
                
                if obj.debug_bit>1, fprintf('InputGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input an InputVars object to constructor of InputGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            import util.plot.MenuItem;
            
            obj.buttons = {};
            obj.menus = {};
            
            obj.fig = util.plot.FigHandler(sprintf('input parameters %s', obj.owner.input_object_name));
            obj.fig.clear;
            obj.fig.height = 10;
            obj.fig.width = 6;
%             movegui(obj.fig.fig, 'center');
            obj.fig.left = obj.fig.left + 12;
            
            %%%%%%%%%%%%%%%%%%%%%%% MENUS %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % MenuItem(parent, text, type, variable, tooltip, separator)
            % obj.addButton(name, text, type, variable, tooltip, separator)
            % menu types: menu, toggle, push, input, input_text, info, custom
            
            
            %%%%%%%%%%% panel inputs %%%%%%%%%%%%%%%%%%
            
            list = obj.owner.list_added_properties;
            
            N = length(list); 
            
            obj.panel_inputs = GraphicPanel(obj.owner, [0 1/(N+2) 1 (N+1)/(N+2)], '', 1, 'graphic_user_interface');
            
            for ii = 1:N
                
                obj.panel_inputs.addButton(['button_' list{ii}], '', 'custom', list{ii}, '', '', 0.4, '', '', obj.owner.comment_dictionary.(list{ii}));
                if obj.owner.logical_dictionary.(list{ii})
                    obj.panel_inputs.addButton(['input_' list{ii}], list{ii}, 'toggle', list{ii}, list{ii}, '', 0.6, obj.color_on, 'k', obj.owner.comment_dictionary.(list{ii}));
                elseif isobject(obj.owner.(list{ii})) % obj.owner.nested_dictionary.(list{ii})
                    obj.panel_inputs.addButton(['input_' list{ii}], list{ii}, 'push', 'Object GUI', '', '', 0.6, '', '', obj.owner.comment_dictionary.(list{ii}));
                else
                    obj.panel_inputs.addButton(['input_' list{ii}], list{ii}, 'input generic', ' ', '', '', 0.6, '', '', obj.owner.comment_dictionary.(list{ii}));
                end
                
            end
            
            obj.panel_inputs.addButton('button_restore_defaults', 'return_to_defaults', 'push', 'restore defaults'); 
            
            obj.panel_inputs.make;
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Position', [0 0 1 1/(N+2)]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE', '', '', 'graphic_user_interface');
            obj.button_close.Callback = @obj.callback_close;
            
            obj.fig.height = N+2;
            obj.fig.center;
            obj.fig.left = 30;
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
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end