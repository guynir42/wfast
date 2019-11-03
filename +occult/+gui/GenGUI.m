classdef GenGUI < handle
    
    properties 
        
        owner@occult.CurveGenerator; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 14;
        big_font_size = 18;
        edit_font_size = 12;
        small_font_size = 10;
        
        color_on = [0 0.3 1]; % apply this color when on
        
        total_delay = 0.05;
        use_uniform_limits = 1;
        show_what = 'lightcurve';
        
        brake_bit = 1;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_controls;
        
        panel_r;
        panel_r2;
        panel_d;
        panel_th;
        
        panel_R;
        panel_b;
        panel_v;
        panel_t;
        
        panel_noise;
        
        panel_scalars;
        
        panel_display;
        
        panel_image;
        button_reset_axes;
        button_show_what;
        button_limits;
        axes_image;
        
        panel_bottom;
        
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
            obj.fig.height = 20;
            obj.fig.width = 25;
            movegui(obj.fig.fig, 'center');
            
            N = 20;
            
            %%%%%%%%%%% control panel %%%%%%%%%%%%%%%%
            
            obj.panel_controls = GraphicPanel(obj.owner, [0/5, 3/N 1/5 6/N], 'controls'); 
            obj.panel_controls.number = 6;
            obj.panel_controls.addButton('button_source', 'use_source_matrix', 'toggle', 'new method', 'source matrix', '', [], '', '', 'Use the old (and fast) source matrix or use the new method that allows binary KBOs'); 
            obj.panel_controls.addButton('button_binary', 'use_binary', 'toggle', 'no binary', 'use binary', '', [], obj.color_on, '', 'Add the secondary component to the occultation. Only works in "new method"'); 
            obj.panel_controls.addButton('button_geometric', 'use_geometric', 'toggle', 'no geometric', 'use geometric', '', [], obj.color_on, '', 'Use geometric approximation for high r and R values. Only works in "new method"'); 
            obj.panel_controls.addButton('input_rho_step', 'rho_step', 'input', 'rho_step= ', '', '', [], '', '', 'Resolution step for new method: how many samples per FSU'); 
            
            obj.panel_controls.make;
            
            
            %%%%%%%%%%% parameters panels %%%%%%%%%%%%
            
            obj.panel_r = GraphicPanel(obj.owner, [1/5 6/N 1/5 3/N], '', 1);
            obj.panel_r.addButton('text', '', 'custom', 'occulter radius'); 
            obj.panel_r.addButton('reset', 'reset_r', 'push', util.text.unicode('recycle'), '', '', 0.15, '', '', 'Reset r to the default value');
            obj.panel_r.addButton('value', 'r', 'input', 'r= ', '', '', 0.85, '', '', 'Input a value for r, the occulter radius');
            obj.panel_r.addButton('play', '', 'custom', util.text.unicode('play'), '', '', 0.15, '', '', 'Play an animation scanning through the range of parameter values');
            obj.panel_r.addButton('slider', '', 'custom', '', '', '', 0.85, '', '', 'Change the value of r, the occulter radius');
            obj.panel_r.margin = [0.01 0.05];
            obj.panel_r.make;
            
            obj.panel_r.text.Enable = 'inactive';
            obj.panel_r.play.Callback = @obj.callback_play_r;
            obj.panel_r.slider.control.Style = 'slider';
            obj.panel_r.slider.control.Min = obj.owner.min_r;
            obj.panel_r.slider.control.Max = obj.owner.max_r;
            obj.panel_r.slider.control.Value = obj.owner.r(1);
            obj.panel_r.slider.control.SliderStep = [0.01 0.1]./(obj.owner.max_r-obj.owner.min_r);
            obj.panel_r.slider.Callback = @obj.callback_slide_r;
            
            obj.panel_r2 = GraphicPanel(obj.owner, [2/5 6/N 1/5 3/N], '', 1);
            obj.panel_r2.addButton('text', '', 'custom', 'companion radius'); 
            obj.panel_r2.addButton('reset', 'reset_r2', 'push', util.text.unicode('recycle'), '', '', 0.15, '', '', 'Reset r2 to the default value');
            obj.panel_r2.addButton('value', 'r2', 'input', 'r_comp= ', '', '', 0.85, '', '', 'Input a value for r2, the companion radius (binary mode only)');
            obj.panel_r2.addButton('play', '', 'custom', util.text.unicode('play'), '', '', 0.15, '', '', 'Play an animation scanning through the range of parameter values');
            obj.panel_r2.addButton('slider', '', 'custom', '', '', '', 0.85, '', '', 'Change the value of r2, the companion radius (binary mode only)');
            obj.panel_r2.margin = [0.01 0.05];
            obj.panel_r2.make;
            
            obj.panel_r2.text.Enable = 'inactive';
            obj.panel_r2.play.Callback = @obj.callback_play_r2;
            obj.panel_r2.slider.control.Style = 'slider';
            obj.panel_r2.slider.control.Min = obj.owner.min_r2;
            obj.panel_r2.slider.control.Max = obj.owner.max_r2;
            obj.panel_r2.slider.control.Value = obj.owner.r2(1);
            obj.panel_r2.slider.control.SliderStep = [0.01 0.1]./(obj.owner.max_r2-obj.owner.min_r2);
            obj.panel_r2.slider.Callback = @obj.callback_slide_r2;
            
            obj.panel_d = GraphicPanel(obj.owner, [3/5 6/N 1/5 3/N], '', 1);
            obj.panel_d.addButton('text', '', 'custom', 'binary separation'); 
            obj.panel_d.addButton('reset', 'reset_d', 'push', util.text.unicode('recycle'), '', '', 0.15, '', '', 'Reset d to the default value');
            obj.panel_d.addButton('value', 'd', 'input', 'd= ', '', '', 0.85, '', '', 'Input a value for d, the binary separation (binary mode only)');
            obj.panel_d.addButton('play', '', 'custom', util.text.unicode('play'), '', '', 0.15, '', '', 'Play an animation scanning through the range of parameter values');
            obj.panel_d.addButton('slider', '', 'custom', '', '', '', 0.85, '', '', 'Change the value of d, the binary separation (binary mode only)');
            obj.panel_d.margin = [0.01 0.05];
            obj.panel_d.make;
            
            obj.panel_d.text.Enable = 'inactive';
            obj.panel_d.play.Callback = @obj.callback_play_d;
            obj.panel_d.slider.control.Style = 'slider';
            obj.panel_d.slider.control.Min = obj.owner.min_d;
            obj.panel_d.slider.control.Max = obj.owner.max_d;
            obj.panel_d.slider.control.Value = obj.owner.d(1);
            obj.panel_d.slider.control.SliderStep = [0.01 0.1]./(obj.owner.max_d-obj.owner.min_d);
            obj.panel_d.slider.Callback = @obj.callback_slide_d;
            
            obj.panel_th = GraphicPanel(obj.owner, [4/5 6/N 1/5 3/N], '', 1);
            obj.panel_th.addButton('text', '', 'custom', 'position angle'); 
            obj.panel_th.addButton('reset', 'reset_th', 'push', util.text.unicode('recycle'), '', '', 0.15, '', '', 'Reset th to the default value');
            obj.panel_th.addButton('value', 'th', 'input', 'th= ', '', '', 0.85, '', '', 'Input a value for th, the binary position angle (binary mode only)');
            obj.panel_th.addButton('play', '', 'custom', util.text.unicode('play'), '', '', 0.15, '', '', 'Play an animation scanning through the range of parameter values');
            obj.panel_th.addButton('slider', '', 'custom', '', '', '', 0.85, '', '', 'Change the value of th, the binary position angle (binary mode only)');
            obj.panel_th.margin = [0.01 0.05];
            obj.panel_th.make;
            
            obj.panel_th.text.Enable = 'inactive';
            obj.panel_th.play.Callback = @obj.callback_play_th;
            obj.panel_th.slider.control.Style = 'slider';
            obj.panel_th.slider.control.Min = obj.owner.min_th;
            obj.panel_th.slider.control.Max = obj.owner.max_th;
            obj.panel_th.slider.control.Value = obj.owner.th(1);
            obj.panel_th.slider.control.SliderStep = [0.01 0.1]./(obj.owner.max_th-obj.owner.min_th);
            obj.panel_th.slider.Callback = @obj.callback_slide_th;
            
            obj.panel_R = GraphicPanel(obj.owner, [1/5 3/N 1/5 3/N], '', 1);
            obj.panel_R.addButton('text', '', 'custom', 'stellar radius'); 
            obj.panel_R.addButton('reset', 'reset_R', 'push', util.text.unicode('recycle'), '', '', 0.15, '', '', 'Reset R to the default value');
            obj.panel_R.addButton('value', 'R', 'input', 'R= ', '', '', 0.85, '', '', 'Input a value for R, the stellar radius');
            obj.panel_R.addButton('play', '', 'custom', util.text.unicode('play'), '', '', 0.15, '', '', 'Play an animation scanning through the range of parameter values');
            obj.panel_R.addButton('slider', '', 'custom', '', '', '', 0.85, '', '', 'Change the value of R, the stellar radius');
            obj.panel_R.margin = [0.01 0.05];
            obj.panel_R.make;
            
            obj.panel_R.text.Enable = 'inactive';
            obj.panel_R.play.Callback = @obj.callback_play_R;
            obj.panel_R.slider.control.Style = 'slider';
            obj.panel_R.slider.control.Min = obj.owner.min_R;
            obj.panel_R.slider.control.Max = obj.owner.max_R;
            obj.panel_R.slider.control.Value = obj.owner.R(1);
            obj.panel_R.slider.control.SliderStep = [0.01 0.1]./(obj.owner.max_R-obj.owner.min_R);
            obj.panel_R.slider.Callback = @obj.callback_slide_R;
            
            obj.panel_b = GraphicPanel(obj.owner, [2/5 3/N 1/5 3/N], '', 1);
            obj.panel_b.addButton('text', '', 'custom', 'impact parameter'); 
            obj.panel_b.addButton('reset', 'reset_b', 'push', util.text.unicode('recycle'), '', '', 0.15, '', '', 'Reset b to the default value');
            obj.panel_b.addButton('value', 'b', 'input', 'b= ', '', '', 0.85, '', '', 'Input a value for b, the impact parameter');
            obj.panel_b.addButton('play', '', 'custom', util.text.unicode('play'), '', '', 0.15, '', '', 'Play an animation scanning through the range of parameter values');
            obj.panel_b.addButton('slider', '', 'custom', '', '', '', 0.85, '', '', 'Change the value of b, the impact parameter');
            obj.panel_b.margin = [0.01 0.05];
            obj.panel_b.make;
            
            obj.panel_b.text.Enable = 'inactive';
            obj.panel_b.play.Callback = @obj.callback_play_b;
            obj.panel_b.slider.control.Style = 'slider';
            obj.panel_b.slider.control.Min = obj.owner.min_b;
            obj.panel_b.slider.control.Max = obj.owner.max_b;
            obj.panel_b.slider.control.Value = obj.owner.b(1);
            obj.panel_b.slider.control.SliderStep = [0.01 0.1]./(obj.owner.max_b-obj.owner.min_b)*2;
            obj.panel_b.slider.Callback = @obj.callback_slide_b;
            
            obj.panel_v = GraphicPanel(obj.owner, [3/5 3/N 1/5 3/N], '', 1);
            obj.panel_v.addButton('text', '', 'custom', 'velocity'); 
            obj.panel_v.addButton('reset', 'reset_v', 'push', util.text.unicode('recycle'), '', '', 0.15, '', '', 'Reset v to the default value');
            obj.panel_v.addButton('value', 'v', 'input', 'v= ', '', '', 0.85, '', '', 'Input a value for v, the occulter velocity');
            obj.panel_v.addButton('play', '', 'custom', util.text.unicode('play'), '', '', 0.15, '', '', 'Play an animation scanning through the range of parameter values');
            obj.panel_v.addButton('slider', '', 'custom', '', '', '', 0.85, '', '', 'Change the value of v, the occulter velocity');
            obj.panel_v.margin = [0.01 0.05];
            obj.panel_v.make;
            
            obj.panel_v.text.Enable = 'inactive';
            obj.panel_v.play.Callback = @obj.callback_play_v;
            obj.panel_v.slider.control.Style = 'slider';
            obj.panel_v.slider.control.Min = obj.owner.min_v;
            obj.panel_v.slider.control.Max = obj.owner.max_v;
            obj.panel_v.slider.control.Value = obj.owner.v(1);
            obj.panel_v.slider.control.SliderStep = [0.01 0.1]./(obj.owner.max_v-obj.owner.min_v)*100;
            obj.panel_v.slider.Callback = @obj.callback_slide_v;
            
            obj.panel_t = GraphicPanel(obj.owner, [4/5 3/N 1/5 3/N], '', 1);
            obj.panel_t.addButton('text', '', 'custom', 'crossing time'); 
            obj.panel_t.addButton('reset', 'reset_t', 'push', util.text.unicode('recycle'), '', '', 0.15, '', '', 'Reset t to the default value');
            obj.panel_t.addButton('value', 't', 'input', 't= ', 'ms', '', 0.85, '', '', 'Input a value for t, the subframe occultation crossing time');
            obj.panel_t.addButton('play', '', 'custom', util.text.unicode('play'), '', '', 0.15, '', '', 'Play an animation scanning through the range of parameter values');
            obj.panel_t.addButton('slider', '', 'custom', '', '', '', 0.85, '', '', 'Change the value of t, the subframe occultation crossing time');
            obj.panel_t.margin = [0.01 0.05];
            obj.panel_t.make;
            
            obj.panel_t.text.Enable = 'inactive';
            obj.panel_t.play.Callback = @obj.callback_play_t;
            obj.panel_t.slider.control.Style = 'slider';
            
            obj.panel_t.slider.control.Min = obj.owner.min_t;
            obj.panel_t.slider.control.Max = obj.owner.max_t;
            obj.panel_t.slider.control.Value = obj.owner.t(1);

            obj.panel_t.slider.control.SliderStep = [0.01 0.1]./(obj.owner.max_t-obj.owner.min_t)*100;
            obj.panel_t.slider.Callback = @obj.callback_slide_t;
            
            %%%%%%%%%%% panel noise %%%%%%%%%%%%%%%%%%
            
            obj.panel_noise = GraphicPanel(obj.owner, [0/3 1/N 1/3 2/N], 'noise', 0);
            obj.panel_noise.addButton('generate', '', 'custom', obj.random_dice, '', 'big', [], '', '', 'Generate new random noise samples');
            obj.panel_noise.addButton('snr', 'snr', 'input', 's= ', '', '', [], '', '', 'Star S/N (flux over rms), determines noise amplitude'); 
            obj.panel_noise.addButton('number', 'num_noise_iterations', 'input', 'N= ', '', '', [], '', '', 'How many noise iterations should be generated');
            obj.panel_noise.margin = [0.02 0.01];
            obj.panel_noise.make; 
            
            obj.panel_noise.generate.Callback = @obj.callback_randomize;
            
            %%%%%%%%%%% panel scalars %%%%%%%%%%%%%%%%
            
            obj.panel_scalars = GraphicPanel(obj.owner, [1/3 1/N 1/3 2/N], 'scalars', 0);
            obj.panel_scalars.addButton('T', 'T', 'input', 'T= ', 'ms', '', [], '', '', 'Exposure time (ms)');
            obj.panel_scalars.addButton('f', 'f', 'input', 'f= ', 'Hz', '', [], '', '', 'Frame rate (Hz)');
            obj.panel_scalars.addButton('W', 'W', 'input', 'W= ', 's', '', [], '', '', 'Time window for lightcurve (seconds)');
            obj.panel_scalars.margin = [0.02 0.01];
            obj.panel_scalars.make;
            
            %%%%%%%%%%% panel display %%%%%%%%%%%%%%%%
            
            obj.panel_display = GraphicPanel(obj.owner, [2/3 1/N 1/3 2/N], 'display', 0);
            obj.panel_display.addButton('num_display', 'num_display', 'input', 'Nlcs= ', '', '', [], '', '', 'Max number of lightcuves to display');
            obj.panel_display.addButton('num_display_noise', 'num_display_noise', 'input', 'Nnoise= ', '', '', [], '', '', 'Max number of noise iterations to display, per lightcurve');
            obj.panel_display.addButton('show_noise', 'show_noise', 'toggle', 'noise', 'noise', '', 1, obj.color_on, 'black', 'Turn on/off display of noise');
            obj.panel_display.margin = [0.02 0.01];
            obj.panel_display.make;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [0 9/N 1 (N-9)/N]);
                        
            obj.makeAxes;
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.9 0.93 0.1 0.07], obj.owner, '', 'custom','new axes');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            obj.button_reset_axes.Tooltip = 'Make a new axes';
            obj.button_show_what = GraphicButton(obj.panel_image, [0.0 0.93 0.1 0.07], obj.owner, '', 'custom', ''); 
            obj.button_show_what.Callback = @obj.callback_show_what;
            obj.button_show_what.control.Style = 'popupmenu';
            obj.button_show_what.control.String = {'lightcurve', 'map'};
            obj.button_limits = GraphicButton(obj.panel_image, [0.0 0.0 0.1 0.07], obj.owner, '', 'custom', 'auto limits', 'uinform limits', 'small'); 
            obj.button_limits.Callback = @obj.callback_limits;
            
            %%%%%%%%%%% panel bottom %%%%%%%%%%%%%%%%%
            
            obj.panel_bottom = GraphicPanel(obj.owner, [0 0 1 1/N], '', 1);
            obj.panel_bottom.addButton('button_close', '', 'custom', 'CLOSE', '', '', 0.1, '', '', 'Close the GUI'); 
            obj.panel_bottom.addButton('button_stop', '', 'custom', 'STOP', '', '', 0.8, '', '', 'Stop all running animations'); 
            obj.panel_bottom.addButton('button_reset', 'reset', 'push', 'RESET', '', '', 0.1, '', '', 'Reset the generator parameters'); 
            obj.panel_bottom.margin = [0.01 0.1];
            obj.panel_bottom.make;
            
%             obj.panel_bottom = uipanel('Position', [0.0 0.0 1 0.1]);
%             
%             obj.button_close = GraphicButton(obj.panel_bottom, [0.0 0 0.1 1], obj.owner, '', 'custom', 'CLOSE');
            obj.panel_bottom.button_close.Callback = @obj.callback_close;
              
%             obj.button_stop = GraphicButton(obj.panel_bottom, [0.1 0 0.8 1], obj.owner, '', 'custom', 'STOP');
            obj.panel_bottom.button_stop.Callback = @obj.callback_stop;
%             
%             obj.button_reset = GraphicButton(obj.panel_bottom, [0.9 0 0.1 1], obj.owner, 'reset', 'push', 'RESET');
            
            obj.update;
            
        end
            
        function makeAxes(obj, ~, ~)
            
            delete(obj.axes_image);
            
            obj.axes_image = axes('Parent', obj.panel_image);
            
            obj.axes_image.Box = 'on';
            
            obj.axes_image.Position(1) = 0.07;
            obj.axes_image.Position(3) = 0.90;
            
            
%             axis(obj.axes_image, 'image');
            
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
            
            if util.text.cs(obj.show_what, 'lightcurve')
                
                obj.makeAxes;
                obj.owner.lc.plot('ax', obj.axes_image);
                title(obj.axes_image, '');
                xlabel(obj.axes_image, '');
                ylabel(obj.axes_image, '');
                xtickformat(obj.axes_image, '%gs');
                
                if obj.use_uniform_limits
                    
                    obj.axes_image.YLim = [0 2];
                    
                    obj.button_limits.String = 'uniform limits';
                    
                    text(-obj.owner.W/2, 0.3, sprintf('runtime= %5.3fs', obj.owner.runtime_get), 'FontSize', obj.font_size);
                    text( obj.owner.W/2, 0.3, sprintf('S/N= %5.3fs', obj.owner.detectionSNR), 'FontSize', obj.font_size, 'HorizontalAlignment', 'right');

                else
                    
                    obj.axes_image.YLim(1) = min(obj.owner.lc.flux(:)).*0.9;
                    obj.axes_image.YLim(2) = max(obj.owner.lc.flux(:)).*1.1;
                    
                    obj.button_limits.String = 'auto limits';
                    
                    text(-obj.owner.W/2, (1+obj.axes_image.YLim(1))/2, sprintf('runtime= %5.3fs', obj.owner.runtime_get), 'FontSize', obj.font_size);
                    text( obj.owner.W/2, (1+obj.axes_image.YLim(1))/2, sprintf('S/N= %5.3fs', obj.owner.detectionSNR), 'FontSize', obj.font_size, 'HorizontalAlignment', 'right');

                end
                    
                obj.axes_image.XLim = [-0.6.*obj.owner.W, 0.6.*obj.owner.W];

            elseif util.text.cs(obj.show_what, 'map')
                
                obj.makeAxes;
                obj.owner.showMap;
                
                text(obj.axes_image.XLim(1)/2, obj.axes_image.YLim(2)/2, sprintf('runtime= %5.3fs', obj.owner.runtime_get), 'FontSize', obj.font_size);
                text(obj.axes_image.XLim(1)/2, obj.axes_image.YLim(2)/1.5, sprintf('S/N= %5.3fs', obj.owner.detectionSNR), 'FontSize', obj.font_size);
                
            else
                error('Unknown "show_what" option: "%s". Use "lightcurve" or "map"', obj.show_what);
            end
            
            if isscalar(obj.owner.r) || all(obj.owner.r==obj.owner.r(1))
                obj.panel_r.slider.Value = obj.owner.r(1);
                obj.panel_r.value.String = ['r= ' num2str(obj.owner.r(1))];
            end
                        
            if isscalar(obj.owner.r2) || all(obj.owner.r2==obj.owner.r2(1))
                obj.panel_r2.slider.Value = obj.owner.r2(1);
                obj.panel_r2.value.String = ['r_comp= ' num2str(obj.owner.r2(1))];
            end
            
            if isscalar(obj.owner.d) || all(obj.owner.d==obj.owner.d(1))
                obj.panel_d.slider.Value = obj.owner.d(1);
                obj.panel_d.value.String = ['d= ' num2str(obj.owner.d(1))];
            end
            
            if isscalar(obj.owner.th) || all(obj.owner.th==obj.owner.th(1))
                obj.panel_th.slider.Value = obj.owner.th(1);
                obj.panel_th.value.String = ['th= ' num2str(obj.owner.th(1))];
            end
            
            if isscalar(obj.owner.R) || all(obj.owner.R==obj.owner.R(1))
                obj.panel_R.slider.Value = obj.owner.R(1);
                obj.panel_R.value.String = ['R= ' num2str(obj.owner.R(1))];
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
            
            if obj.brake_bit
                obj.panel_r.play.String = util.text.unicode('play');
                obj.panel_r2.play.String = util.text.unicode('play');
                obj.panel_d.play.String = util.text.unicode('play');
                obj.panel_th.play.String = util.text.unicode('play');
                obj.panel_R.play.String = util.text.unicode('play');
                obj.panel_b.play.String = util.text.unicode('play');
                obj.panel_v.play.String = util.text.unicode('play');
                obj.panel_t.play.String = util.text.unicode('play');
            else
                obj.panel_r.play.String = util.text.unicode('stop');
                obj.panel_r2.play.String = util.text.unicode('stop');
                obj.panel_d.play.String = util.text.unicode('stop');
                obj.panel_th.play.String = util.text.unicode('stop');
                obj.panel_R.play.String = util.text.unicode('stop');
                obj.panel_b.play.String = util.text.unicode('stop');
                obj.panel_v.play.String = util.text.unicode('stop');
                obj.panel_t.play.String = util.text.unicode('stop');
            end
            
            obj.panel_r.slider.control.Min = obj.owner.r_range(1);
            obj.panel_r2.slider.control.Min = obj.owner.r2_range(1);
            obj.panel_d.slider.control.Min = obj.owner.d_range(1);
            obj.panel_th.slider.control.Min = obj.owner.th_range(1);
            obj.panel_R.slider.control.Min = obj.owner.R_range(1);
            obj.panel_b.slider.control.Min = obj.owner.b_range(1);
            obj.panel_v.slider.control.Min = obj.owner.v_range(1);
            obj.panel_t.slider.control.Min = obj.owner.t_range(1);
            
            obj.panel_r.slider.control.Max = obj.owner.r_range(2);
            obj.panel_r2.slider.control.Max = obj.owner.r2_range(2);
            obj.panel_d.slider.control.Max = obj.owner.d_range(2);
            obj.panel_th.slider.control.Max = obj.owner.th_range(2);
            obj.panel_R.slider.control.Max = obj.owner.R_range(2);
            obj.panel_b.slider.control.Max = obj.owner.b_range(2);
            obj.panel_v.slider.control.Max = obj.owner.v_range(2);
            obj.panel_t.slider.control.Max = obj.owner.t_range(2);
            
            h = legend(obj.axes_image, 'Location', 'Southeast');
            h.FontSize = 10;
            
        end
        
        function str = random_dice(obj, number)
            
            if nargin<2 || isempty(number)
                number = 3;
            end
            
            rnd_num = randi([-5,0],3);
            
            str = '';
            
            for ii = 1:number
                str(ii) = util.text.unicode('dice', rnd_num(ii));
            end
             
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_bottom) && is_valid(obj.panel_bottom);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_play_r(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: play r'); end
            
            if obj.brake_bit==0
                
                obj.brake_bit = 1;
                obj.update;
                
            else
                
                val = obj.panel_r.slider.Value;
                mn = obj.panel_r.slider.control.Min;
                mx = obj.panel_r.slider.control.Max;
                if val==mx
                    val = mn;
                end
                
                step = obj.panel_r.slider.control.SliderStep(1)*(mx-mn);

                vec = val:step:mx;

                obj.brake_bit = 0;

                for ii = 1:length(vec)

                    if obj.brake_bit, break; end

                    obj.owner.r = vec(ii);
                    obj.update;
                    delay = max([0 obj.total_delay-obj.owner.runtime_get]); 
                    pause(delay);
                    drawnow;

                end

                obj.brake_bit = 1;
                
                obj.update;
                
            end
            
        end
        
        function callback_play_r2(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: play r2'); end
            
            if obj.brake_bit==0
                
                obj.brake_bit = 1;
                obj.update;
                
            else
                
                val = obj.panel_r2.slider.Value;
                mn = obj.panel_r2.slider.control.Min;
                mx = obj.panel_r2.slider.control.Max;
                if val==mx
                    val = mn;
                end
                
                step = obj.panel_r2.slider.control.SliderStep(1)*(mx-mn);

                vec = val:step:mx;

                obj.brake_bit = 0;

                for ii = 1:length(vec)

                    if obj.brake_bit, break; end

                    obj.owner.r2 = vec(ii);
                    obj.update;
                    delay = max([0 obj.total_delay-obj.owner.runtime_get]); 
                    pause(delay);
                    drawnow;

                end

                obj.brake_bit = 1;
                
                obj.update;
                
            end
                
        end
        
        function callback_play_d(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: play d'); end
            
            if obj.brake_bit==0
                
                obj.brake_bit = 1;
                obj.update;
                
            else
                
                val = obj.panel_d.slider.Value;
                mn = obj.panel_d.slider.control.Min;
                mx = obj.panel_d.slider.control.Max;
                if val==mx
                    val = mn;
                end
                
                step = obj.panel_d.slider.control.SliderStep(1)*(mx-mn);

                vec = val:step:mx;

                obj.brake_bit = 0;

                for ii = 1:length(vec)

                    if obj.brake_bit, break; end

                    obj.owner.d = vec(ii);
                    obj.update;
                    delay = max([0 obj.total_delay-obj.owner.runtime_get]); 
                    pause(delay);
                    drawnow;

                end

                obj.brake_bit = 1;
                
                obj.update;
                
            end
                
        end
        
        function callback_play_th(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: play th'); end
            
            if obj.brake_bit==0
                
                obj.brake_bit = 1;
                obj.update;
                
            else
                
                val = obj.panel_th.slider.Value;
                mn = obj.panel_th.slider.control.Min;
                mx = obj.panel_th.slider.control.Max;
                if val==mx
                    val = mn;
                end
                
                step = obj.panel_th.slider.control.SliderStep(1)*(mx-mn);

                vec = val:step:mx;

                obj.brake_bit = 0;

                for ii = 1:length(vec)

                    if obj.brake_bit, break; end

                    obj.owner.th = vec(ii);
                    obj.update;
                    delay = max([0 obj.total_delay-obj.owner.runtime_get]); 
                    pause(delay);
                    drawnow;

                end

                obj.brake_bit = 1;
                
                obj.update;
                
            end
                
        end
        
        function callback_play_R(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: play R'); end
            
            if obj.brake_bit==0
                
                obj.brake_bit = 1;
                obj.update;
                
            else
                
                val = obj.panel_R.slider.Value;
                mn = obj.panel_R.slider.control.Min;
                mx = obj.panel_R.slider.control.Max;
                if val==mx
                    val = mn;
                end
                
                step = obj.panel_R.slider.control.SliderStep(1)*(mx-mn);

                vec = val:step:mx;

                obj.brake_bit = 0;

                for ii = 1:length(vec)

                    if obj.brake_bit, break; end

                    obj.owner.R = vec(ii);
                    obj.update;
                    delay = max([0 obj.total_delay-obj.owner.runtime_get]); 
                    pause(delay);
                    drawnow;

                end

                obj.brake_bit = 1;
                
                obj.update;
                
            end
                
        end
        
        function callback_play_b(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: play b'); end
            
            if obj.brake_bit==0
                
                obj.brake_bit = 1;
                obj.update;
                
            else
                
                val = obj.panel_b.slider.Value;
                mn = obj.panel_b.slider.control.Min;
                mx = obj.panel_b.slider.control.Max;
                if val==mx
                    val = mn;
                end
                
                step = obj.panel_b.slider.control.SliderStep(1)*(mx-mn);

                vec = val:step:mx;

                obj.brake_bit = 0;

                for ii = 1:length(vec)

                    if obj.brake_bit, break; end

                    obj.owner.b = vec(ii);
                    obj.update;
                    delay = max([0 obj.total_delay-obj.owner.runtime_get]); 
                    pause(delay);
                    drawnow;

                end

                obj.brake_bit = 1;
                
                obj.update;
                
            end
            
        end
        
        function callback_play_v(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: play v'); end
            
            if obj.brake_bit==0
                
                obj.brake_bit = 1;
                obj.update;
                
            else
                
                val = obj.panel_v.slider.Value;
                mn = obj.panel_v.slider.control.Min;
                mx = obj.panel_v.slider.control.Max;
                if val==mx
                    val = mn;
                end
                
                step = obj.panel_v.slider.control.SliderStep(1)*(mx-mn);

                vec = val:step:mx;

                obj.brake_bit = 0;

                for ii = 1:length(vec)

                    if obj.brake_bit, break; end

                    obj.owner.v = vec(ii);
                    obj.update;
                    delay = max([0 obj.total_delay-obj.owner.runtime_get]); 
                    pause(delay);
                    drawnow;

                end

                obj.brake_bit = 1;
                
                obj.update;
                
            end
            
        end
        
        function callback_play_t(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: play t'); end
            
            if obj.brake_bit==0
                
                obj.brake_bit = 1;
                obj.update;
                
            else
                
                val = obj.panel_t.slider.Value;
                mn = obj.panel_t.slider.control.Min;
                mx = obj.panel_t.slider.control.Max;
                if val==mx
                    val = mn;
                end
                
                step = obj.panel_t.slider.control.SliderStep(1)*(mx-mn);

                vec = val:step:mx;

                obj.brake_bit = 0;

                for ii = 1:length(vec)

                    if obj.brake_bit, break; end

                    obj.owner.t = vec(ii);
                    obj.update;
                    delay = max([0 obj.total_delay-obj.owner.runtime_get]); 
                    pause(delay);
                    drawnow;

                end

                obj.brake_bit = 1;
                
                obj.update;
                
            end
            
        end
        
        function callback_slide_r(obj, hndl, ~)
            
            if obj.debug_bit, disp('callback: slide r'); end
            
            obj.owner.r = hndl.Value;
            
            obj.update;
            
        end
        
        function callback_slide_r2(obj, hndl, ~)
            
            if obj.debug_bit, disp('callback: slide r2'); end
            
            obj.owner.r2 = hndl.Value;
            
            obj.update;
            
        end
        
        function callback_slide_d(obj, hndl, ~)
            
            if obj.debug_bit, disp('callback: slide d'); end
            
            obj.owner.d = hndl.Value;
            
            obj.update;
            
        end
        
        function callback_slide_th(obj, hndl, ~)
            
            if obj.debug_bit, disp('callback: slide th'); end
            
            obj.owner.th = hndl.Value;
            
            obj.update;
            
        end
        
        function callback_slide_R(obj, hndl, ~)
            
            if obj.debug_bit, disp('callback: slide R'); end
            
            obj.owner.R = hndl.Value;
            
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
        
        function callback_show_what(obj, hndl, ~)
           
            if obj.debug_bit, disp('Callback: show what'); end
            
            obj.show_what = hndl.String{hndl.Value};
            
            obj.update;
            
        end
        
        function callback_limits(obj, ~, ~)
            
            if obj.debug_bit, disp('Callback: limits'); end
            
            obj.use_uniform_limits = ~obj.use_uniform_limits;
            
            obj.update;
            
        end
        
        function callback_randomize(obj, ~, ~)
            
            if obj.debug_bit, disp('Callback: generate noise'); end
            
            obj.owner.generateNoise;
            
            obj.panel_noise.generate.String = obj.random_dice;
                        
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