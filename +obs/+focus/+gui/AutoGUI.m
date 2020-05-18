classdef AutoGUI < handle
    
    properties 
        
        owner@obs.focus.AutoFocus; % link back to containg object

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
        panel_photometry;
        panel_tip_tilt;
        
        panel_results;
        
        panel_run; 
        
        panel_close;
        button_close;
        
        panel_image;
        button_reset_axes;
        button_num_plots;
        axes_image;
    
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = AutoGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'obs.focus.AutoFocus')
                
                if obj.debug_bit>1, fprintf('AutoGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input an obs.focuse.AutoFocus object to constructor of AutoGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.MenuItem;
            
            obj.buttons = {};
            obj.menus = {};
            
            obj.fig = util.plot.FigHandler('auto focus');
            obj.fig.clear;
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 36;
            movegui(obj.fig.fig, 'center');
            
            
            %%%%%%%%%%%%%%%%%%%%%%% MENUS %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % MenuItem(parent, text, type, variable, tooltip, separator)
            % obj.addButton(name, text, type, variable, tooltip, separator)
            % menu types: menu, toggle, push, input, input_text, info, custom
            
            
            
            N = 10; % number of buttons on left side
            width = 0.25;
            %%%%%%%%%%%%%%%%%%% LEFT SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            pos = N;
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            % Add buttons using obj.addButton(button_name, var_name='', type='', str1='', str2='', font_size='', split=1, color_on=[], color_off=[], tooltip)
            
            num_buttons = 4;
            pos = pos-num_buttons;
            obj.panel_controls = GraphicPanel(obj.owner, [0 pos/N width num_buttons/N], 'controls', 1); % last input is for vertical (default)
            obj.panel_controls.addButton('input_exp_time', 'expT', 'input', 'expT= ', '', '', 0.5, '', '', 'exposure time (seconds) for each image'); 
            obj.panel_controls.addButton('input_batch_size', 'batch_size', 'input', 'N= ', '', '', 0.5, '', '', 'number of images for each batch / each focus position'); 
            obj.panel_controls.addButton('button_loop_back', 'use_loop_back', 'toggle', 'no loop', 'loop back', '', 0.5, obj.color_on, '', 'scan the focus positions in both directions'); 
            obj.panel_controls.addButton('button_placeholder', 'use_fit_curves', 'toggle', 'no fit', 'use fit', '', 0.5, obj.color_on, '', 'fit each star width curve to a parabola'); 
            obj.panel_controls.addButton('input_range', 'range', 'input', 'range= ', '', '', 0.5, '', '', 'how far in either direction from the focus point we want to scan'); 
            obj.panel_controls.addButton('input_step', 'step', 'input', 'step= ', '', '', 0.5, '', '', 'step size for sampling the focus positions'); 
                        
            obj.panel_controls.number = num_buttons;
            
            obj.panel_controls.make;
            
            %%%%%%%%%%% panel photometry %%%%%%%%%%%%%%%

            num_buttons = 2;
            pos = pos-num_buttons;
            obj.panel_photometry = GraphicPanel(obj.owner, [0 pos/N width num_buttons/N], 'photometry', 1); % last input is for vertical (default)
            obj.panel_photometry.addButton('input_cut_size', 'cut_size', 'input', 'cut_size= ', '', '', 0.5, '', '', 'size of cutouts (in pixels)'); 
            obj.panel_photometry.addButton('input_aperture', 'aperture', 'input', 'ap= ', '', '', 0.5, '', '', 'photometric aperture radius (in pixels)'); 
            obj.panel_photometry.addButton('input_annulus', 'annulus', 'input', 'ann= ', '', '', 0.5, '', '', 'photometric annulus inner and outer radius (in pixels)'); 
            obj.panel_photometry.addButton('input_gaussian', 'gaussian', 'input', 'gauss= ', '', '', 0.5, '', '', 'size of Gaussian sigma '); 
            
            obj.panel_photometry.number = num_buttons;
            
            obj.panel_photometry.make;
            
            %%%%%%%%%%% panel tip tilt %%%%%%%%%%%%%%%

            num_buttons = 3;
            pos = pos-num_buttons;
            obj.panel_tip_tilt = GraphicPanel(obj.owner, [0 pos/N width num_buttons/N], 'tip tilt', 1); % last input is for vertical (default)
            obj.panel_tip_tilt.addButton('input_angle', 'angle', 'input', 'angle= ', '', '', 0.5, '', '', 'angle between camera x axis and positive-tip actuator (East direction)'); 
            obj.panel_tip_tilt.addButton('input_diameter', 'spider_diameter', 'input', 'diameter= ', '', '', 0.5, '', '', 'diameter of focus spider, measuring at the position of actuators (in centimeters)'); 
            obj.panel_tip_tilt.addButton('input_pixel_size', 'pixel_size', 'input', 'pix. size= ', '', '', 0.5, '', '', 'size of pixels (in microns)'); 
            obj.panel_tip_tilt.addButton('button_placeholder', '', 'custom', '', '', '', 0.5); 
            obj.panel_tip_tilt.addButton('button_show', 'show', 'push', 'show surface fit', '', '', 1, '', '', 'show the best 2D fit to the positions of each minimal focus point.'); 
            
            obj.panel_tip_tilt.number = num_buttons;
            
            obj.panel_tip_tilt.make;
            
            
            %%%%%%%%%%% panel results %%%%%%%%%%%%%%%%%%
            
            obj.panel_results = GraphicPanel(obj.owner, [width, 0.9, 1-width, 0.1], 'results'); 
            obj.panel_results.addButton('button_results', 'printout', 'info');
            obj.panel_results.make;
            
            %%%%%%%%%%% panel run %%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_run = GraphicPanel(obj.owner, [width, 0.0, 1-width, 0.1], 'run'); 
            obj.panel_run.addButton('button_autofocus', '', 'custom', 'Start autofocus run', '', '', 0.8, '', '', 'run autofocus routine using the camera'); 
            obj.panel_run.addButton('button_focuser', 'cam.focuser', 'push', 'manual', '', '', 0.2, '', '', 'open the focuser GUI to set the focus/tip/tilt manually'); 
            obj.panel_run.make;
            
            obj.panel_run.button_autofocus.Callback = @obj.callback_autofocus;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [width 0.1 1-width 0.8]);
                        
            obj.makeAxes;
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.85 0.95 0.15 0.05], obj.owner, '', 'custom', 'new axes', '');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            obj.button_reset_axes.Tooltip = 'create a new image axis, zoomed out and with default contrast limits'; 
                        
            obj.button_num_plots = GraphicButton(obj.panel_image, [0.00 0.0 0.15 0.05], obj.owner, 'num_plots', 'input', 'num= ', '');
            obj.button_num_plots.Tooltip = 'how many plots to show at the same time'; 
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Position', [0 0 width 1/N]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            obj.update;
            
        end
            
        function makeAxes(obj, ~, ~)
            
            delete(obj.axes_image);
            
            obj.axes_image = axes('Parent', obj.panel_image);

            axis(obj.axes_image, 'normal');
            
            obj.update;
            
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
           
            if obj.owner.cam.is_running_focus
                obj.panel_run.button_autofocus.String = 'STOP';
            else
                obj.panel_run.button_autofocus.String = 'Start autofocus run';
            end
            
            obj.owner.calculate;
            obj.owner.plot;
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_autofocus(obj, ~, ~)
            
            if obj.debug_bit>1, disp('callback: start/stop'); end
            
            if obj.owner.cam.is_running_focus
                obj.owner.cam.brake_bit = 1;
            else
                obj.owner.cam.autofocus;
            end
            
            obj.update;
            
        end
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end