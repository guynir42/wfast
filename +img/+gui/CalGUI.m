classdef CalGUI < handle
   
    properties % useful stuff
        
        cal@img.Calibration;
        fig@util.plot.FigHandler;
        
        buttons = {};
        menus = {};
        
        font_size = 13;
        big_font_size = 16;
        edit_font_size = 12;
        small_font_size = 10;
        
        color_on = [0 0.3 1];

        debug_bit = 1;
        
    end
    
    properties % gui objects
        
        panel_dark;
%         button_make_dark;
%         button_reset_dark;
%         button_reader_dark;
%         button_browse_dark;
%         button_num_dark;
        
        panel_flat;
%         button_make_flat;
%         button_reset_flat;
%         button_reader_flat;
%         button_browse_flat;
%         button_num_flat;
        
        panel_analysis;
        button_calc_gain;
        button_calc_lightcurves
        button_calc_covariance;
        
        panel_controls;
        button_use_flat;
        button_sub_median;
        button_interp_mask;
        button_conv_mask;
        input_mask_sigma;
        input_mask_var_thresh;
        input_replace_value;
        button_autosave;
                
        panel_actions;
        button_show;
        button_save;
        button_save_browse;
        button_load;
        button_load_browse;
        
        panel_image;
        axes_dark;
        axes_var;
        axes_mask;
        axes_flat;
        
        panel_close;
        button_close;
        
        panel_progress;
        button_progress;
        
        panel_stop;
        button_stop;
        
    end
    
    properties(Hidden=true)
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = CalGUI(calibration_obj)
           
            if nargin<1 || isempty(calibration_obj) || ~isa(calibration_obj, 'img.Calibration')
                error('must input a Calibration object to create GUI!');
            end
            
            obj.cal = calibration_obj;
            
        end
        
    end
    
    methods % make/update GUI

        function makeGUI(obj)
        
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            
            if isempty(obj.fig)
                obj.fig = util.plot.FigHandler('calibration GUI');
            end
            
            obj.buttons = {};
            
            obj.fig.reset;
            obj.fig.height = 20;
            obj.fig.width = 30;
            obj.fig.bottom = 2;
            obj.fig.name = 'Calibration';
            
            %%%%%%%%%%%%%%%%%%%%%%% dark panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_dark = GraphicPanel(obj.cal, [0.0 0.7 0.2 0.3], 'make dark'); 
            obj.panel_dark.addButton('button_num', 'num_darks', 'info', 'Ndark= ', '', '', 1, '', '', 'number of files in current calibration dark'); 
            obj.panel_dark.addButton('button_reader', 'reader_dark', 'push', 'reader', '', '', 0.6, '', '', 'open the GUI for the dark Reader object');
            obj.panel_dark.addButton('button_browse', 'browseDark', 'push', 'browse', '', '', 0.4, '', '', 'browse the for a folder with dark images');
            obj.panel_dark.addButton('button_make', 'makeDark', 'push', 'makeDark', '', '', 0.6, '', '', 'call makeDark to read all files in the folder into the calibration file');
            obj.panel_dark.addButton('button_reset', 'resetDark', 'push', 'reset', '', '', 0.4, '', '', 'reset the dark calibration');
            obj.panel_dark.margin = [0.02 0.05]; 
            obj.panel_dark.make;
            
            
%             obj.panel_dark = uipanel('Title', 'dark', 'Units', 'Normalized', 'Position', [0.0 0.7 0.2 0.3]);
%             num = 3;
%             obj.button_num_dark = GraphicButton(obj.panel_dark, [0 2/num 1 1/num], obj.cal, 'num_darks', 'info', 'Ndark= ');
%             obj.button_reader_dark = GraphicButton(obj.panel_dark, [0 1/num 0.7 1/num], obj.cal, 'reader_dark', 'push', 'reader');
%             obj.button_reader_dark.control.BackgroundColor = [0.0 0.7 1];
%             obj.button_browse_dark = GraphicButton(obj.panel_dark, [0.7 1/num 0.3 1/num], obj.cal, 'browseDark', 'push', 'browse');
%             obj.button_make_dark = GraphicButton(obj.panel_dark, [0 0/num 0.7 1/num], obj.cal, 'makeDark', 'push', 'makeDark');            
%             obj.button_reset_dark = GraphicButton(obj.panel_dark, [0.7 0/num 0.3 1/num], obj.cal, 'resetDark', 'push', 'reset');
                        
            %%%%%%%%%%%%%%%%%%%%%%% flat panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      
            obj.panel_dark = GraphicPanel(obj.cal, [0.0 0.4 0.2 0.3], 'make flat'); 
            obj.panel_dark.addButton('button_num', 'num_flats', 'info', 'Nflat= ', '', '', 1, '', '', 'number of files in current calibration flat'); 
            obj.panel_dark.addButton('button_reader', 'reader_flat', 'push', 'reader', '', '', 0.6, '', '', 'open the GUI for the flat Reader object');
            obj.panel_dark.addButton('button_browse', 'browseFlat', 'push', 'browse', '', '', 0.4, '', '', 'browse the for a folder with flat images');
            obj.panel_dark.addButton('button_make', 'makeFlat', 'push', 'makeFlat', '', '', 0.6, '', '', 'call makeFlat to read all files in the folder into the calibration file');
            obj.panel_dark.addButton('button_reset', 'resetFlat', 'push', 'reset', '', '', 0.4, '', '', 'reset the flat calibration');
            obj.panel_dark.margin = [0.02 0.05]; 
            obj.panel_dark.make;
            
%             obj.panel_flat = uipanel('Title', 'flat', 'Units', 'Normalized', 'Position', [0.0 0.4 0.2 0.3]);
%             num = 3;
%             obj.button_num_flat = GraphicButton(obj.panel_flat, [0 2/num 1 1/num], obj.cal, 'num_flats', 'info', 'Nflat= ');
%             obj.button_reader_flat = GraphicButton(obj.panel_flat, [0 1/num 0.7 1/num], obj.cal, 'reader_flat', 'push', 'reader');
%             obj.button_reader_flat.control.BackgroundColor = [0.0 0.7 1];
%             obj.button_browse_flat = GraphicButton(obj.panel_flat, [0.7 1/num 0.3 1/num], obj.cal, 'browseFlat', 'push', 'browse');
%             obj.button_make_flat = GraphicButton(obj.panel_flat, [0 0/num 0.7 1/num], obj.cal, 'makeFlat', 'push', 'makeFlat');            
%             obj.button_reset_flat = GraphicButton(obj.panel_flat, [0.7 0/num 0.3 1/num], obj.cal, 'resetFlat', 'push', 'reset');
            
            %%%%%%%%%%%%%%%%%%%%%%% controls panel %%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_controls = GraphicPanel(obj.cal, [0.8 0.3 0.2 0.7], 'controls'); 
            obj.panel_controls.addButton('button_use_flat', 'use_flat', 'toggle', 'not using flat', 'using flat', '', 1, obj.color_on, '', 'apply flat field when doing calibration'); 
            obj.panel_controls.addButton('button_sub_median', 'use_subtract_median', 'toggle', 'not using median', 'using median', '', 1, obj.color_on, '', 'subtract the median of each incoming image'); 
            obj.panel_controls.addButton('button_use_interp_mask', 'use_interp_mask', 'toggle', 'interp mask is off', 'interp mask is on', '', 1, obj.color_on, '', 'replace bad pixels with interpolation of surrounding pixels'); 
            obj.panel_controls.addButton('button_use_conv_mask', 'use_conv_mask', 'toggle', 'conv mask is off', 'conv mask is on', '', 1, obj.color_on, '', 'replace bad pixels with the average of surrounding pixels'); 
            obj.panel_controls.addButton('input_mask_sigma', 'dark_mask_sigma', 'input', 'mask sigma= ', '', '', 1, '', '', 'pixels above this number of times the noise RMS are flagged as bad'); 
            obj.panel_controls.addButton('input_mask_var_max', 'dark_mask_var_max', 'input', 'mask var max= ', '', '', 1, '', '', 'pixels with variance above this threshold are flagged as bad'); 
            obj.panel_controls.addButton('input_mask_var_min', 'dark_mask_var_min', 'input', 'mask var min= ', '', '', 1, '', '', 'pixels with variance below this value are flagged as bad'); 
            obj.panel_controls.addButton('input_replace_value', 'replace_value', 'input', 'replace val= ', '', '', 1, '', '', 'replace value of bad pixels with this (unless using conv/interp mask');
            obj.panel_controls.margin = [0.05 0.02];
            obj.panel_controls.make;
            
%             obj.panel_controls = uipanel('Title', 'controls', 'Units', 'Normalized', 'Position', [0.8 0.3 0.2 0.7]);
%             num = 8;
%             obj.button_use_flat = GraphicButton(obj.panel_controls, [0 7/num 1 1/num], obj.cal, 'use_flat', 'toggle', 'no flat', 'use flat');
%             obj.button_sub_median = GraphicButton(obj.panel_controls, [0 6/num 1 1/num], obj.cal, 'use_subtract_median', 'toggle', 'no sub median', 'sub median');
%             obj.button_interp_mask = GraphicButton(obj.panel_controls, [0 5/num 1 1/num], obj.cal, 'use_interp_mask', 'toggle', 'no interp mask', 'use interp mask');
%             obj.button_conv_mask = GraphicButton(obj.panel_controls, [0 4/num 1 1/num], obj.cal, 'use_conv_mask', 'toggle', 'no conv mask', 'use conv mask');
%             obj.input_mask_sigma = GraphicButton(obj.panel_controls, [0 3/num 1 1/num], obj.cal, 'dark_mask_sigma', 'input', 'mask sigma= ');
%             obj.input_mask_var_thresh = GraphicButton(obj.panel_controls, [0 2/num 1 1/num], obj.cal, 'dark_mask_var_thresh', 'input', 'var thresh= ');
%             obj.input_replace_value = GraphicButton(obj.panel_controls, [0 1/num 1 1/num], obj.cal, 'replace_value', 'input', 'replace_value= ');
%             obj.button_autosave = GraphicButton(obj.panel_controls, [0 0/num 1 1/num], obj.cal, 'autosave', 'toggle', 'no autosave', 'autosave');
            
            %%%%%%%%%%%%%%%%%%%%%%% utils panel %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_analysis = GraphicPanel(obj.cal, [0.0 0.1 0.2 0.3], 'analysis'); 
            obj.panel_analysis.addButton('button_calc_gain', 'use_calc_gain', 'toggle', 'gain calc is off', 'gain calc is on', '', 1, obj.color_on, '', 'calculate the gain while running makeFlat'); 
            obj.panel_analysis.addButton('button_calc_lightcurves', 'use_calc_lightcurve', 'toggle', 'lightcurves is off', 'lightcurves is on', '', 1, obj.color_on, '', 'calculate lightcurves while running makeFlat'); 
            obj.panel_analysis.addButton('button_calc_covariance', 'use_calc_covariance', 'toggle', 'calc cov is off', 'calc cov is on', '', 1, obj.color_on, '', 'calculate the covariance matrix while running makeFlat'); 
            obj.panel_analysis.margin = [0.02 0.05]; 
            obj.panel_analysis.make;
            
%             obj.panel_analysis = uipanel('Title', 'analysis', 'Units', 'Normalized', 'Position', [0.0 0.1 0.2 0.3]);
%             num = 3;
%             obj.button_calc_gain = GraphicButton(obj.panel_analysis, [0 2/num 1 1/num], obj.cal, 'use_calc_gain', 'toggle', 'no gain', 'calc gain');
%             obj.button_calc_lightcurves = GraphicButton(obj.panel_analysis, [0 1/num 1 1/num], obj.cal, 'use_calc_lightcurve', 'toggle', 'no lightcurves', 'calc lightcurves');
%             obj.button_calc_covariance = GraphicButton(obj.panel_analysis, [0 0/num 1 1/num], obj.cal, 'use_calc_covariance', 'toggle', 'no covariance', 'calc covariance');
                        
            %%%%%%%%%%%%%%%%%%%%%%% actions panel %%%%%%%%%%%%%%%%%%%%%%%%%

            obj.panel_actions = GraphicPanel(obj.cal, [0.8 0.0 0.2 0.3], 'actions'); 
            obj.panel_actions.addButton('button_show', 'show', 'push', 'Show', '', '', 1, '', '', 'show the dark/flat images');
            obj.panel_actions.addButton('button_load', 'loadByDate', 'push', 'Load', '', '', 0.5, '', '', 'load a calibration file based on today''s date');
            obj.panel_actions.addButton('button_load_browse', 'load_browse', 'push', 'browse', '', '', 0.5, '', '', 'browse for a calibration file to load');
            obj.panel_actions.addButton('button_save', 'save', 'push', 'Save', '', '', 0.5, '', '', 'save calibration files in the calibration folder and along with the raw data files');
            obj.panel_actions.addButton('button_save_browse', 'save_browse', 'push', 'browse', '', '', 0.5, '', '', 'browse the filesystem and save a new calibration file');
            obj.panel_actions.addButton('button_use_autosave', 'use_autosave', 'toggle', 'autosave is off', 'autosave is on', '', 1, obj.color_on, '', 'automatically save calibration after running makeDark or makeFlat'); 
            obj.panel_actions.margin = [0.02 0.02]; 
            obj.panel_actions.make;
            
%             obj.panel_actions = uipanel('Title', 'actions', 'Units', 'Normalized', 'Position', [0.8 0.0 0.2 0.3]);
%             num = 3;
%             obj.button_show = GraphicButton(obj.panel_actions, [0 2/num 1 1/num], obj.cal, 'show', 'push', 'Show');
%             obj.button_load = GraphicButton(obj.panel_actions, [0 1/num 0.5 1/num], obj.cal, 'load', 'push', 'Load');
%             obj.button_load_browse = GraphicButton(obj.panel_actions, [0.5 1/num 0.5 1/num], obj.cal, 'load_browse', 'push', 'browse');
%             obj.button_save = GraphicButton(obj.panel_actions, [0.0 0/num 0.5 1/num], obj.cal, 'save', 'push', 'Save');
%             obj.button_save_browse = GraphicButton(obj.panel_actions, [0.5 0/num 0.5 1/num], obj.cal, 'save_browse', 'push', 'browse');
                        
            %%%%%%%%%%%%%%%%%%%%%%% close panel %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Title', '', 'Units', 'Normalized', 'Position', [0 0 0.2 0.1]);
            
            obj.button_close = GraphicButton(obj.panel_close, [0.05 0.1 0.9 0.8], obj.cal, '', 'custom', 'CLOSE GUI');
            obj.button_close.Callback = @obj.callback_close;
            
            %%%%%%%%%%%%%%%%%%%%%%% progress panel %%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_progress = uipanel('Title', '', 'Units', 'Normalized', 'Position', [0.2 0.9 0.6 0.1]);
            
            obj.button_progress = GraphicButton(obj.panel_progress, [0.01 0.1 0.98 0.8], obj.cal, '', 'custom');
            obj.button_progress.Callback = @obj.update;
            
            %%%%%%%%%%%%%%%%%%%%%%% stop panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_stop = uipanel('Title', '', 'Units', 'Normalized', 'Position', [0.2 0.0 0.6 0.1]);
            
            obj.button_stop = GraphicButton(obj.panel_stop, [0.01 0.1 0.98 0.8], obj.cal, '', 'custom');
            obj.button_stop.Callback = @obj.callback_stop;
            
            %%%%%%%%%%%%%%%%%%%%%%% image panel %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            obj.panel_image = uipanel('Title', 'image', 'Units', 'Normalized', 'Position', [0.2 0.1 0.6 0.8]);
            
            obj.update;
            
        end
        
        function update(obj,~,~)
                        
            if ~obj.checkGUI
                return;
            end
            
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%% progress panel %%%%%%%%%%%%%%%%%%%%%%%%
           
            obj.button_progress.String = obj.cal.prog.show;
            obj.button_progress.FontSize = obj.small_font_size;
            
            %%%%%%%%%%%%%%%%%%%%%%% stop panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if isempty(obj.cal.mode)
                obj.button_stop.String = 'IDLE';
            else
                if obj.cal.brake_bit
                    obj.button_stop.String = 'Continue';
                else
                    obj.button_stop.String = 'STOP';
                end
            end
            
            drawnow;
            
        end
        
        function updateGUI(obj)
            obj.update;
        end
        
        function c = check(obj)
            c = obj.checkGUI;
        end            
        
        function c = checkGUI(obj)
            
            c = ~isempty(obj) && ~isempty(obj.panel_image) && isvalid(obj.panel_image);
            
        end
        
    end
    
    methods % callbacks
               
        %%%%%%%%%%%%%%%%%%%%%%% close panel %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% stop panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function callback_stop(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: stop'); end
            
            if ~isempty(obj.cal.mode) && obj.cal.brake_bit==0
                obj.cal.brake_bit = 1;
            else
                obj.cal.brake_bit = 0;
                obj.cal.loop;
            end
            
            obj.update;
            
        end
        
    end
    
end