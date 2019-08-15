classdef GenGUI < handle
    
    properties 
        
        owner@occult.CurveGenerator; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 16;
        edit_font_size = 14;
        small_font_size = 12;
        
        total_delay = 0.05;
        
        brake_bit = 1;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_R;
        panel_r;
        panel_b;
        panel_v;
        panel_t;
        
        panel_noise;
        
        panel_scalars;
        
        panel_display;
        
        panel_image;
        button_reset_axes;
        axes_image;
        
        panel_bottom;
        button_stop;
        button_close;
        button_reset;
        
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = GenGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'occult.CurveGenerator')
                
                if obj.debug_bit, fprintf('GenGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input a CurveGenerator to constructor of GenGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            
            obj.buttons = {};
            
            obj.fig = util.plot.FigHandler('curve generator');
            
            obj.fig.clear;
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 25;
            movegui(obj.fig.fig, 'center');
            
            %%%%%%%%%%% parameters panels %%%%%%%%%%%%
            
            obj.panel_R = GraphicPanel(obj.owner, [0.0 0.2 0.2 0.2], '', 1);
            obj.panel_R.addButton('reset', 'reset_R', 'push', 'R', '', '', 0.1);
            obj.panel_R.addButton('value', 'R', 'input', 'R= ', '', '', 0.9);
            obj.panel_R.addButton('play', '', 'custom', 'P', '', '', 0.1);
            obj.panel_R.addButton('slider', '', 'custom', '', '', '', 0.9);
            obj.panel_R.make;
            obj.panel_R.play.Callback = @obj.callback_play_R;
            obj.panel_R.slider.control.Style = 'slider';
            obj.panel_R.slider.control.Min = obj.owner.min_R;
            obj.panel_R.slider.control.Max = obj.owner.max_R;
%             obj.panel_R.slider.control.SliderStep = [0.01 0.1];
            obj.panel_R.slider.Callback = @obj.callback_slide_R;
            
            obj.panel_r = GraphicPanel(obj.owner, [0.2 0.2 0.2 0.2], '', 1);
            obj.panel_r.addButton('reset', 'reset_r', 'push', 'R', '', '', 0.1);
            obj.panel_r.addButton('value', 'r', 'input', 'r= ', '', '', 0.9);
            obj.panel_r.addButton('play', '', 'custom', 'P', '', '', 0.1);
            obj.panel_r.addButton('slider', '', 'custom', '', '', '', 0.9);
            obj.panel_r.make;
            obj.panel_r.play.Callback = @obj.callback_play_r;
            obj.panel_r.slider.control.Style = 'slider';
            obj.panel_r.slider.control.Min = obj.owner.min_r;
            obj.panel_r.slider.control.Max = obj.owner.max_r;
%             obj.panel_r.slider.control.SliderStep = [0.01 0.1];
            obj.panel_r.slider.Callback = @obj.callback_slide_r;
            
            obj.panel_b = GraphicPanel(obj.owner, [0.4 0.2 0.2 0.2], '', 1);
            obj.panel_b.addButton('reset', 'reset_b', 'push', 'R', '', '', 0.1);
            obj.panel_b.addButton('value', 'b', 'input', 'b= ', '', '', 0.9);
            obj.panel_b.addButton('play', '', 'custom', 'P', '', '', 0.1);
            obj.panel_b.addButton('slider', '', 'custom', '', '', '', 0.9);
            obj.panel_b.make;
            obj.panel_b.play.Callback = @obj.callback_play_b;
            obj.panel_b.slider.control.Style = 'slider';
            obj.panel_b.slider.control.Min = obj.owner.min_b;
            obj.panel_b.slider.control.Max = obj.owner.max_b;
%             obj.panel_b.slider.control.SliderStep = [0.01 0.1];
            obj.panel_b.slider.Callback = @obj.callback_slide_b;
            
            obj.panel_v = GraphicPanel(obj.owner, [0.6 0.2 0.2 0.2], '', 1);
            obj.panel_v.addButton('reset', 'reset_v', 'push', 'R', '', '', 0.1);
            obj.panel_v.addButton('value', 'v', 'input', 'v= ', '', '', 0.9);
            obj.panel_v.addButton('play', '', 'custom', 'P', '', '', 0.1);
            obj.panel_v.addButton('slider', '', 'custom', '', '', '', 0.9);
            obj.panel_v.make;
            obj.panel_v.play.Callback = @obj.callback_play_v;
            obj.panel_v.slider.control.Style = 'slider';
%             obj.panel_v.slider.control.SliderStep = [0.1 1];
            obj.panel_v.slider.control.Min = obj.owner.min_v;
            obj.panel_v.slider.control.Max = obj.owner.max_v;
            obj.panel_v.slider.Callback = @obj.callback_slide_v;
            
            obj.panel_t = GraphicPanel(obj.owner, [0.8 0.2 0.2 0.2], '', 1);
            obj.panel_t.addButton('reset', 'reset_t', 'push', 'R', '', '', 0.1);
            obj.panel_t.addButton('value', 't', 'input', 't= ', 'ms', '', 0.9);
            obj.panel_t.addButton('play', '', 'custom', 'P', '', '', 0.1);
            obj.panel_t.addButton('slider', '', 'custom', '', '', '', 0.9);
            obj.panel_t.make;
            obj.panel_t.play.Callback = @obj.callback_play_t;
            obj.panel_t.slider.control.Style = 'slider';
            obj.panel_t.slider.control.Min = obj.owner.min_t;
            obj.panel_t.slider.control.Max = obj.owner.max_t;
%             obj.panel_t.slider.control.SliderStep = [0.5 5];
            obj.panel_t.slider.Callback = @obj.callback_slide_t;
            
            %%%%%%%%%%% panel noise %%%%%%%%%%%%%%%%%%
            
            obj.panel_noise = GraphicPanel(obj.owner, [0/3 0.1 1/3 0.1], 'noise', 0);
            obj.panel_noise.addButton('generate', 'generateNoise', 'push', 'generate', '', '');
            obj.panel_noise.addButton('snr', 'snr', 'input', 's= ', '', 'small'); 
            obj.panel_noise.addButton('number', 'num_noise_iterations', 'input', 'N= ', '', '');
            obj.panel_noise.make; 
            
            %%%%%%%%%%% panel scalars %%%%%%%%%%%%%%%%
            
            obj.panel_scalars = GraphicPanel(obj.owner, [1/3 0.1 1/3 0.1], 'scalars', 0);
            obj.panel_scalars.addButton('T', 'T', 'input', 'T= ', 'ms', '');
            obj.panel_scalars.addButton('f', 'f', 'input', 'f= ', 'Hz', '');
            obj.panel_scalars.addButton('W', 'W', 'input', 'W= ', 's', '');
            
            obj.panel_scalars.make;
            
            %%%%%%%%%%% panel display %%%%%%%%%%%%%%%%
            
            obj.panel_display = GraphicPanel(obj.owner, [2/3 0.1 1/3 0.1], 'display', 0);
            obj.panel_display.addButton('num_display', 'num_display', 'input', 'N= ', '', '');
            obj.panel_display.addButton('num_display_noise', 'num_display_noise', 'input', 'Nn= ', '', '');
            obj.panel_display.addButton('show_noise', 'show_noise', 'toggle', 'noise', 'noise', '', 1, 'red', 'black');
            
            obj.panel_display.make;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.0 0.4 1 0.6]);
                        
            obj.makeAxes;
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.9 0.95 0.1 0.05], obj.owner, '', 'custom','reset');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            
            %%%%%%%%%%% panel bottom %%%%%%%%%%%%%%%%%
            
            obj.panel_bottom = uipanel('Position', [0.0 0.0 1 0.1]);
            
            obj.button_close = GraphicButton(obj.panel_bottom, [0.0 0 0.1 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            obj.button_stop = GraphicButton(obj.panel_bottom, [0.1 0 0.8 1], obj.owner, '', 'custom', 'STOP');
            obj.button_stop.Callback = @obj.callback_stop;
            
            obj.button_reset = GraphicButton(obj.panel_bottom, [0.9 0 0.1 1], obj.owner, 'reset', 'push', 'RESET');
            
            obj.update;
            
        end
            
        function makeAxes(obj, ~, ~)
            
            delete(obj.axes_image);
            
            obj.axes_image = axes('Parent', obj.panel_image);
            
            axis(obj.axes_image, 'image');
            
        end
                
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            if obj.owner.lc.is_updated==0
                obj.owner.getLightCurves;
            end
            
            if obj.owner.lc.is_noise_updated==0
                obj.owner.generateNoise;
            end
            
            obj.owner.lc.plot('ax', obj.axes_image);
            title(obj.axes_image, '');
            xlabel(obj.axes_image, '');
            ylabel(obj.axes_image, '');
            xtickformat(obj.axes_image, '%gs');
            obj.axes_image.YLim = [0 1.35];
            obj.axes_image.XLim = [-0.6.*obj.owner.W, 0.6.*obj.owner.W];
            
            text(-1.2, 0.2, sprintf('runtime= %5.3fs', obj.owner.runtime_get));
            
            if isscalar(obj.owner.R) || all(obj.owner.R==obj.owner.R(1))
                obj.panel_R.slider.Value = obj.owner.R(1);
                obj.panel_R.value.String = ['R= ' num2str(obj.owner.R(1))];
            end
            
            if isscalar(obj.owner.r) || all(obj.owner.r==obj.owner.r(1))
                obj.panel_r.slider.Value = obj.owner.r(1);
                obj.panel_r.value.String = ['r= ' num2str(obj.owner.r(1))];
            end
            if isscalar(obj.owner.b) || all(obj.owner.b==obj.owner.b(1))
                obj.panel_b.slider.Value = obj.owner.b(1);
                obj.panel_b.value.String = ['b= ' num2str(obj.owner.b(1))];
            end
            if isscalar(obj.owner.v) || all(obj.owner.v==obj.owner.v(1))
                obj.panel_v.slider.Value = obj.owner.v(1);
                obj.panel_v.value.String = ['v= ' num2str(obj.owner.v(1))];
            end
            if isscalar(obj.owner.t) || all(obj.owner.t==obj.owner.t(1))
                obj.panel_t.slider.Value = obj.owner.t(1);
                obj.panel_t.value.String = ['t= ' num2str(obj.owner.t(1))];
            end
            
            h = legend(obj.axes_image, 'Location', 'Southeast');
            h.FontSize = 10;
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_bottom) && isvalid(obj.panel_bottom);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_play_R(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: play R'); end
            
%             mn = obj.panel_R.slider.control.Min;
            mn = obj.panel_R.slider.Value;
            mx = obj.panel_R.slider.control.Max;
            step = obj.panel_R.slider.control.SliderStep(1);
            
            vec = mn:step:mx;
            
            obj.brake_bit = 0;
            
            for ii = 1:length(vec)
                
                if obj.brake_bit, break; end
                
                obj.owner.R = vec(ii);
                obj.update;
                delay = max([0 obj.total_delay-obj.owner.runtime_get]); 
                pause(delay);
                drawnow;
            
            end
            
        end
        
        function callback_play_r(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: play r'); end
            
%             mn = obj.panel_r.slider.control.Min;
            mn = obj.panel_r.slider.Value;
            mx = obj.panel_r.slider.control.Max;
            step = obj.panel_r.slider.control.SliderStep(1);
            
            vec = mn:step:mx;
            
            obj.brake_bit = 0;
            
            for ii = 1:length(vec)
                
                if obj.brake_bit, break; end
                
                obj.owner.r = vec(ii);
                obj.update;
                delay = max([0 obj.total_delay-obj.owner.runtime_get]); 
                pause(delay);
                drawnow;
            
            end
            
        end
        
        function callback_play_b(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: play b'); end
            
%             mn = obj.panel_b.slider.control.Min;
            mn = obj.panel_b.slider.Value;
            mx = obj.panel_b.slider.control.Max;
            step = obj.panel_b.slider.control.SliderStep(1);
            
            vec = mn:step:mx;
            
            obj.brake_bit = 0;
            
            for ii = 1:length(vec)
                
                if obj.brake_bit, break; end
                
                obj.owner.b = vec(ii);
                obj.update;
                delay = max([0 obj.total_delay-obj.owner.runtime_get]); 
                pause(delay);
                drawnow;
                
            end
            
        end
        
        function callback_play_v(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: play v'); end
            
%             mn = obj.panel_v.slider.control.Min;
            mn = obj.panel_v.slider.Value;
            mx = obj.panel_v.slider.control.Max;
            step = obj.panel_v.slider.control.SliderStep(1);
            
            vec = mn:step:mx;
            
            obj.brake_bit = 0;
            
            for ii = 1:length(vec)
                
                if obj.brake_bit, break; end
                
                obj.owner.v = vec(ii);
                obj.update;
                delay = max([0 obj.total_delay-obj.owner.runtime_get]); 
                pause(delay);
                drawnow;
            
            end
            
        end
        
        function callback_play_t(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: play t'); end
            
%             mn = obj.panel_t.slider.control.Min;
            mn = obj.panel_t.slider.Value;
            mx = obj.panel_t.slider.control.Max;
            step = 1;
            
            vec = mn:step:mx;
            
            obj.brake_bit = 0;
            
            for ii = 1:length(vec)
                
                if obj.brake_bit, break; end
                
                obj.owner.t = vec(ii);
                obj.update;
                delay = max([0 obj.total_delay-obj.owner.runtime_get]); 
                pause(delay);
                drawnow;
            
            end
            
        end
        
        function callback_slide_R(obj, hndl, ~)
            
            if obj.debug_bit, disp('callback: slide R'); end
            
            obj.owner.R = hndl.Value;
            
            obj.update;
            
        end
        
        function callback_slide_r(obj, hndl, ~)
            
            if obj.debug_bit, disp('callback: slide r'); end
            
            obj.owner.r = hndl.Value;
            
            obj.update;
            
        end
        
        function callback_slide_b(obj, hndl, ~)
            
            if obj.debug_bit, disp('callback: slide b'); end
            
            obj.owner.b = hndl.Value;
            
            obj.update;
            
        end
        
        function callback_slide_v(obj, hndl, ~)
            
            if obj.debug_bit, disp('callback: slide v'); end
            
            obj.owner.v = hndl.Value;
            
            obj.update;
            
        end
        
        function callback_slide_t(obj, hndl, ~)
            
            if obj.debug_bit, disp('callback: slide t'); end
            
            obj.owner.t = hndl.Value;
            
            obj.update;
            
        end
        
        function callback_stop(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: stop'); end
            
            obj.brake_bit = 1;
            drawnow;
            
            obj.update;
            
        end
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end