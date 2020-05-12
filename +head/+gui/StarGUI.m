classdef StarGUI < handle
    
    properties 
        
        owner@head.Star; % link back to containing object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 16;
        edit_font_size = 14;
        small_font_size = 12;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_chooser;
        button_back;
        button_forward;
        input_number;
        
        panel_input;
    
        panel_close;
        button_close;
        
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = StarGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'head.Star')
                
                if obj.debug_bit>1, fprintf('StarGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                                
            else
                error('Input a head.Star to constructor of StarGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            
            if isempty(obj.fig)
                obj.fig = util.plot.FigHandler('Stars');
            end
            
            obj.fig.reset;
            obj.fig.left = 6;
            if ~isempty(obj.owner.head.gui)
                obj.fig.left = obj.owner.head.gui.fig.left+obj.owner.head.gui.fig.width;
            end
            obj.fig.bottom = 4;
            obj.fig.height = 16;
            obj.fig.width = 10;
%             obj.fig.name = '...';
%             movegui(obj.fig.fig, 'center');
            
            %%%%%%%%%%% panel chooser %%%%%%%%%%%%%%%%
            
            obj.panel_chooser = uipanel('Title', 'Choose star', 'Position', [0 0.9 1 0.1]);
            obj.button_back = GraphicButton(obj.panel_chooser, [0/3 0 1/3 1], obj.owner, '', 'custom', 'back');            
            obj.input_number = GraphicButton(obj.panel_chooser, [1/3 0 1/3 1], obj.owner, '', 'input custom', ' ');
            obj.button_forward = GraphicButton(obj.panel_chooser, [2/3 0 1/3 1], obj.owner, '', 'custom', 'forward');
            
            obj.button_back.Callback = @obj.callback_back;
            obj.input_number.Callback = @obj.callback_number;
            obj.button_forward.Callback = @obj.callback_forward;
            
            %%%%%%%%%%% panel input %%%%%%%%%%%%%%%%%%
            
            obj.panel_input = GraphicPanel(obj.owner, [0 0.1 1 0.8], 'input data');
            obj.panel_input.addButton('button_name', '', 'custom', 'name: ','','', 0.2);
            obj.panel_input.addButton('input_name', 'name', 'input text', ' ', '', '', 0.8);
            
            obj.panel_input.addButton('button_mag', '', 'custom', 'mag: ','','', 0.2);
            obj.panel_input.addButton('input_mag', 'mag', 'input', ' ', '', '', 0.8);
                        
            obj.panel_input.addButton('button_primary', '', 'custom', 'primary index: ','','', 0.7);
            obj.panel_input.addButton('input_primary', 'primary_index', 'input', ' ', '', '', 0.3);
                        
            obj.panel_input.addButton('button_coord', '', 'custom', 'x,y:','','', 0.2);
            obj.panel_input.addButton('input_coord', '', 'input custom', ' ', '', '', 0.8);
            
            obj.panel_input.addButton('button_type', '', 'custom', 'type: ','','', 0.2);
            obj.panel_input.addButton('input_type', 'type', 'input text', ' ', '', '', 0.8);
            
            obj.panel_input.number = 6;
            obj.panel_input.make;
            obj.panel_input.input_coord.Callback = @obj.callback_coord;
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
                        
            obj.panel_close = uipanel('Position', [0 0 1 0.1]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            obj.update;
            
        end
        
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
            
            obj.input_number.String = [num2str(obj.owner.index) '/' num2str(length(obj.owner.head.stars))];
            
            if isempty(obj.owner.primary_index)
                obj.panel_input.button_coord.String = 'x,y:';
                obj.panel_input.input_coord.String = num2str([obj.owner.anchor_x obj.owner.anchor_y]);
            else
                obj.panel_input.button_coord.String = 'r,\theta:';
                obj.panel_input.input_coord.String = num2str([obj.owner.offset_sep obj.owner.offset_angle]);
            end
            
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
        end
        
        function c = check(obj)
           
            c = ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
        
        function set.owner(obj, new_owner)
            
            obj.owner = new_owner;
            obj.panel_input.owner = new_owner;
            
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.owner = new_owner;
            end
            
        end
        
        function callback_back(obj, ~, ~)
            
            num = obj.owner.index;
            
            if num==1
                num = length(obj.owner.head.stars);
            else
                num = num - 1;
            end
            
            obj.owner = obj.owner.head.stars(num);
            
            obj.update;
            
        end
        
        function callback_forward(obj, ~, ~)
            
            num = obj.owner.index;
            
            if num==length(obj.owner.head.stars)
                num = 1;
            else
                num = num + 1;
            end
            
            obj.owner = obj.owner.head.stars(num);
            
            obj.update;
            
        end
        
        function callback_number(obj, hndl, ~)
            
            values = util.text.extract_numbers(hndl.String);
            num = values{1};
            
            if num>0 && num<=length(obj.owner.head.stars)
                obj.owner = obj.owner.head.stars(num);
            end
            
            obj.update;
            
        end
        
        function callback_coord(obj, hndl, ~)
            
            num = util.text.extract_numbers(hndl.String);
            
            num = cell2mat(num);
                        
            if ~isempty(num)
                val1 = num(1);
            else
                val1 = [];
            end
            
            if length(num)>1
                val2 = num(2);
            else
                val2 = [];
            end
            
            if isempty(obj.owner.primary_index)
                obj.owner.anchor_x = val1;
                obj.owner.anchor_y = val2;
            else
                obj.owner.offset_sep = val1;
                obj.owner.offset_angle = val2;
            end
            
            obj.update;
            
        end
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: close'); end
            
            obj.owner.head.gui.update;
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end