classdef ContrastLimits < handle
% manipulates the CLim of some axes object. Also has autodyn. 
% replaces a uipanel inside a GUI object or figure. 
%
% TEST PROTOCOL: f1=figure; c=util.ContrastLimits([], f1); delete(f1);

    properties % objects
        
        ax; % axes to manipulate
        panel; % uipanel with all buttons
        parent; % some figure or panel
        gui; % loops to itself!!!
        
        buttons = {};
        
        input_clim;
        button_max;
        button_min;
        
        button_high;
        slider_high;
        button_low;
        slider_low;
        
        button_auto;
        button_reset;
        
    end
    
    properties % switches and controls
        
        max_val = 1e5;
        min_val = 0;
        clim;
        
        font_size = 18;
        edit_font_size = 14;
        small_font_size = 12;
        position = [0 0 1 1];
        vertical_mode = 1;
        
        debug_bit = 1;
        
    end
    
    properties(Hidden=true)
        version = 1.00;
    end
    
    methods % constructor
        
        function obj = ContrastLimits(ax, parent, position, vertical_mode)
           
%             if nargin<1 || isempty(ax)
%                 error('must supply an axes object!')
%             end
            
            if nargin<2 || isempty(parent)
                parent = gcf;            
            end
            
            if nargin<3 || isempty(position)
                position = [];
            end
            
            if nargin<4 || isempty(vertical_mode)
                vertical_mode = [];
            end
            
            obj.ax = ax;
            obj.parent = parent;
            obj.gui = obj;
            
            if ~isempty(position)
                obj.position = position;
            end
            
            if ~isempty(vertical_mode)
                obj.vertical_mode = vertical_mode;
            end
            
            if ~obj.check
                obj.make;
            end
            
        end
        
    end
        
    methods % setters
        
        function set.clim(obj, val)
            
            if isempty(val)
                obj.clim = [];
                return;
            end
            
            if all(val==0) || val(1)>val(2)
                return;
            end
            
            obj.clim = val;
            
            obj.ax.CLim = val;
            
        end
        
        function set.position(obj, val)
            
            if ~all(obj.position==val)
                obj.position = val;
                obj.make;
            end
            
        end
        
        function set.vertical_mode(obj, val)
            
            if obj.vertical_mode~=val
                obj.vertical_mode = val;
                obj.make;
            end
            
        end
        
    end
        
    methods % make, update, check
        
        function p = calcPositions(obj)
            
            num = 5;
            
            if obj.vertical_mode
                
                p = {};
                p{end+1} = [0.0 4/num 1.0 1/num];
                p{end+1} = [0.0 3/num 0.5 1/num];
                p{end+1} = [0.5 3/num 0.5 1/num];
                p{end+1} = [0.0 2/num 0.2 1/num];
                p{end+1} = [0.2 2/num 0.8 1/num];
                p{end+1} = [0.0 1/num 0.2 1/num];
                p{end+1} = [0.2 1/num 0.8 1/num];
                p{end+1} = [0.0 0/num 0.5 1/num];
                p{end+1} = [0.5 0/num 0.5 1/num];
                
            else
                
                p = {};
                p{end+1} = [0.0/num 0.0 1.0/num 1];
                p{end+1} = [1.0/num 0.0 0.5/num 1];
                p{end+1} = [1.5/num 0.0 0.5/num 1];
                p{end+1} = [2.0/num 0.0 0.2/num 1];
                p{end+1} = [2.2/num 0.0 0.8/num 1];
                p{end+1} = [3.0/num 0.0 0.2/num 1];
                p{end+1} = [3.2/num 0.0 0.8/num 1];
                p{end+1} = [4.0/num 0.0 0.5/num 1];
                p{end+1} = [4.5/num 0.0 0.5/num 1];
                
            end
            
        end
        
        function make(obj)
            
            obj.buttons = {};
            
            import util.plot.GraphicButton;
            
            p = obj.calcPositions;
            
            obj.panel = uipanel(obj.parent, 'Title', 'contrast', 'Units', 'Normalized', 'Position', obj.position);
            
            obj.input_clim = GraphicButton(obj.panel, p{1}, obj, 'clim', 'input', 'CLim= ');
            
            obj.button_min = GraphicButton(obj.panel, p{2}, obj, 'changeMin', 'push', 'min');            
            obj.button_max = GraphicButton(obj.panel, p{3}, obj, 'changeMax', 'push', 'max');   
            obj.button_low = GraphicButton(obj.panel, p{4}, obj, '', 'custom', 'low:');
            obj.button_low.font_size = 'small';
            obj.slider_low = uicontrol(obj.panel, 'Style', 'Slider',  'Units', 'Normalized', 'Position', p{5}, 'Callback', @obj.callback_sliders);
            obj.button_high = GraphicButton(obj.panel, p{6}, obj, '', 'custom', 'high:');
            obj.button_high.font_size = 'small';
            obj.slider_high = uicontrol(obj.panel, 'Style', 'Slider',  'Units', 'Normalized', 'Position', p{7}, 'Callback', @obj.callback_sliders);
            
            obj.button_auto = GraphicButton(obj.panel, p{8}, obj, 'autodyn', 'push', 'auto');
            obj.button_reset = GraphicButton(obj.panel, p{9}, obj, 'reset', 'push', 'reset');
            
            obj.update;
            
        end
        
        function update(obj, ~, ~)
            
            if ~obj.check
                return;
            end
            
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            if length(obj.clim)==2
                obj.input_clim.String = sprintf('CLim= %4.2f %4.2f', obj.clim(1), obj.clim(2));
            end
            
            obj.button_min.String = num2str(obj.min_val);
            obj.button_max.String = num2str(obj.max_val);
                        
            obj.slider_low.Min = obj.min_val;
            obj.slider_low.Max = obj.max_val;
            obj.slider_low.SliderStep = [ 0.01 0.1 ];
            
            if isempty(obj.clim)
                obj.slider_low.Value = obj.slider_low.Min;
            else
                if obj.clim(1)>obj.min_val
                    obj.slider_low.Value = obj.clim(1);
                else
                    obj.slider_low.Value = obj.min_val;
                end
            end
            
            obj.slider_high.Min = obj.min_val;
            obj.slider_high.Max = obj.max_val;
            obj.slider_high.SliderStep = [ 0.01 0.1 ];
            
            if isempty(obj.clim)
                obj.slider_high.Value = obj.slider_high.Min;
            else
                if obj.clim(2)<obj.max_val
                    obj.slider_high.Value = obj.clim(2);
                else
                    obj.slider_high.Value = obj.max_val;
                end
            end
            
        end
        
        function val = check(obj)
            
            val = ~isempty(obj.slider_low) && isvalid(obj.slider_low);
            
        end
        
    end
    
    methods % callbacks
        
        function changeMax(obj, ~, ~)
            
            rep = util.text.inputdlg('max value: ', '');
            
            if ~isempty(rep)
                val = util.text.extract_numbers(rep);
                obj.max_val = val{1};
                
                if ~isempty(obj.clim) && obj.max_val<obj.clim(2)
                    if obj.max_val>obj.clim(1)
                        obj.clim = [obj.clim(1) obj.max_val];
                    else
                        obj.clim = [obj.min_val obj.max_val];
                    end
                end
                
            end
            
            obj.update;
            
        end
        
        function changeMin(obj, ~, ~)
            
            rep = util.text.inputdlg('min value: ', '');
            
            if ~isempty(rep)
                val = util.text.extract_numbers(rep);
                obj.min_val = val{1};
                
                if ~isempty(obj.clim) && obj.min_val<obj.clim(1)
                    if obj.min_val<obj.clim(2)
                        obj.clim = [obj.min_val obj.min_val+(obj.clim(2)-obj.clim(1))];
                    else
                        obj.clim = [obj.min_val obj.clim(2)];
                    end
                end
                
            end
            
            obj.update;
            
        end
        
        function autodyn(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: autodyn'); end
            
            im = findobj(obj.ax, 'type', 'Image');
            if isempty(im)
                return;
            end
            
            new_lim = util.img.autodyn(im.CData);
            if new_lim(1)>=new_lim(2)
                new_lim(1) = new_lim(2).*0.8;
            end
            
            obj.clim = new_lim
            
            if obj.clim(1)<0 
                obj.min_val = -2.^ceil(log2(abs(obj.clim(1)))); 
            else
                obj.min_val = 0;
            end
            
            obj.max_val = 2.^ceil(log2(obj.clim(2)));
            
            obj.update;
            
        end
        
        function reset(obj, ~, ~)
                        
            if obj.debug_bit, disp('callback: reset'); end
            
            obj.clim = [];
            
            obj.ax.CLimMode = 'auto';
            
            obj.update;
            
        end
        
        function callback_sliders(obj, ~, ~)
                        
            if obj.debug_bit, disp('callback: sliders'); end
            
            obj.clim = [obj.slider_low.Value obj.slider_high.Value];
            
            obj.update;
            
        end
        
    end
    
end