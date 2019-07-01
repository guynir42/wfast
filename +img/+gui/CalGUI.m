classdef CalGUI < handle
   
    properties % useful stuff
        
        cal@img.Calibration;
        fig@util.plot.FigHandler;
        
        buttons = {};
        
        font_size = 18;
        edit_font_size = 14;
        small_font_size = 12;

        debug_bit = 1;
        
    end
    
    properties % gui objects
        
        panel_dark;
        button_make_dark;
        button_reset_dark;
        button_reader_dark;
        button_browse_dark;
        button_num_dark;
        
        panel_flat;
        button_make_flat;
        button_reset_flat;
        button_reader_flat;
        button_browse_flat;
        button_num_flat;
        
        panel_utils;
        button_use_mv_points;
        button_calc_gain;
        
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
            
            if isempty(obj.fig)
                obj.fig = util.plot.FigHandler('calibration GUI');
            end
            
            obj.buttons = {};
            
            obj.fig.reset;
            obj.fig.height = 20;
            obj.fig.width = 26;
            obj.fig.bottom = 2;
            obj.fig.name = 'Calibration';
            
            %%%%%%%%%%%%%%%%%%%%%%% dark panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_dark = uipanel('Title', 'dark', 'Units', 'Normalized', 'Position', [0.0 0.7 0.2 0.3]);
            
            num = 3;
            
            obj.button_num_dark = GraphicButton(obj.panel_dark, [0 2/num 1 1/num], obj.cal, 'num_darks', 'info', 'Ndark= ');
            obj.button_reader_dark = GraphicButton(obj.panel_dark, [0 1/num 0.7 1/num], obj.cal, 'reader_dark', 'push', 'reader');
            obj.button_reader_dark.control.BackgroundColor = [0.0 0.7 1];
            obj.button_browse_dark = GraphicButton(obj.panel_dark, [0.7 1/num 0.3 1/num], obj.cal, 'browseDark', 'push', 'browse');
            obj.button_make_dark = GraphicButton(obj.panel_dark, [0 0/num 0.7 1/num], obj.cal, 'makeDark', 'push', 'makeDark');            
            obj.button_reset_dark = GraphicButton(obj.panel_dark, [0.7 0/num 0.3 1/num], obj.cal, 'resetDark', 'push', 'reset');
                        
            %%%%%%%%%%%%%%%%%%%%%%% flat panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
            obj.panel_flat = uipanel('Title', 'flat', 'Units', 'Normalized', 'Position', [0.0 0.4 0.2 0.3]);
            
            num = 3;
                        
            obj.button_num_flat = GraphicButton(obj.panel_flat, [0 2/num 1 1/num], obj.cal, 'num_flats', 'info', 'Nflat= ');
            obj.button_reader_flat = GraphicButton(obj.panel_flat, [0 1/num 0.7 1/num], obj.cal, 'reader_flat', 'push', 'reader');
            obj.button_reader_flat.control.BackgroundColor = [0.0 0.7 1];
            obj.button_browse_flat = GraphicButton(obj.panel_flat, [0.7 1/num 0.3 1/num], obj.cal, 'browseFlat', 'push', 'browse');
            obj.button_make_flat = GraphicButton(obj.panel_flat, [0 0/num 0.7 1/num], obj.cal, 'makeFlat', 'push', 'makeFlat');            
            obj.button_reset_flat = GraphicButton(obj.panel_flat, [0.7 0/num 0.3 1/num], obj.cal, 'resetFlat', 'push', 'reset');
            
            %%%%%%%%%%%%%%%%%%%%%%% controls panel %%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_controls = uipanel('Title', 'controls', 'Units', 'Normalized', 'Position', [0.8 0.3 0.2 0.7]);
            
            num = 8;
            
            obj.button_use_flat = GraphicButton(obj.panel_controls, [0 7/num 1 1/num], obj.cal, 'use_flat', 'toggle', 'no flat', 'use flat');
            obj.button_sub_median = GraphicButton(obj.panel_controls, [0 6/num 1 1/num], obj.cal, 'use_subtract_median', 'toggle', 'no sub median', 'sub median');
            obj.button_interp_mask = GraphicButton(obj.panel_controls, [0 5/num 1 1/num], obj.cal, 'use_interp_mask', 'toggle', 'no interp mask', 'use interp mask');
            obj.button_conv_mask = GraphicButton(obj.panel_controls, [0 4/num 1 1/num], obj.cal, 'use_conv_mask', 'toggle', 'no conv mask', 'use conv mask');
            obj.input_mask_sigma = GraphicButton(obj.panel_controls, [0 3/num 1 1/num], obj.cal, 'dark_mask_sigma', 'input', 'mask sigma= ');
            obj.input_mask_var_thresh = GraphicButton(obj.panel_controls, [0 2/num 1 1/num], obj.cal, 'dark_mask_var_thresh', 'input', 'var thresh= ');
            obj.input_replace_value = GraphicButton(obj.panel_controls, [0 1/num 1 1/num], obj.cal, 'replace_value', 'input', 'replace_value= ');
            obj.button_autosave = GraphicButton(obj.panel_controls, [0 0/num 1 1/num], obj.cal, 'autosave', 'toggle', 'no autosave', 'autosave');
            
            %%%%%%%%%%%%%%%%%%%%%%% actions panel %%%%%%%%%%%%%%%%%%%%%%%%%
                        
            obj.panel_actions = uipanel('Title', 'actions', 'Units', 'Normalized', 'Position', [0.8 0.0 0.2 0.3]);

            num = 3;
            
            obj.button_show = GraphicButton(obj.panel_actions, [0 2/num 1 1/num], obj.cal, 'show', 'push', 'Show');
            obj.button_load = GraphicButton(obj.panel_actions, [0 1/num 0.5 1/num], obj.cal, 'load', 'push', 'Load');
            obj.button_load_browse = GraphicButton(obj.panel_actions, [0.5 1/num 0.5 1/num], obj.cal, 'load_browse', 'push', 'browse');
            obj.button_save = GraphicButton(obj.panel_actions, [0.0 0/num 0.5 1/num], obj.cal, 'save', 'push', 'Save');
            obj.button_save_browse = GraphicButton(obj.panel_actions, [0.5 0/num 0.5 1/num], obj.cal, 'save_browse', 'push', 'browse');
                        
            %%%%%%%%%%%%%%%%%%%%%%% utils panel %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_utils = uipanel('Title', 'utilities', 'Units', 'Normalized', 'Position', [0.0 0.1 0.2 0.3]);
            
            num = 3;
            
            obj.button_use_mv_points = GraphicButton(obj.panel_utils, [0 2/num 1 1/num], obj.cal, 'use_mv_points', 'toggle', 'no save MV', 'save MV');
            obj.button_calc_gain = GraphicButton(obj.panel_utils, [0 1/num 1 1/num], obj.cal, 'calcGain', 'push', 'calcGain');
                        
            %%%%%%%%%%%%%%%%%%%%%%% close panel %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Title', '', 'Units', 'Normalized', 'Position', [0 0 0.2 0.1]);
            
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.cal, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            %%%%%%%%%%%%%%%%%%%%%%% progress panel %%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_progress = uipanel('Title', '', 'Units', 'Normalized', 'Position', [0.2 0.9 0.6 0.1]);
            
            obj.button_progress = GraphicButton(obj.panel_progress, [0 0 1 1], obj.cal, '', 'custom');
            obj.button_progress.Callback = @obj.update;
            
            %%%%%%%%%%%%%%%%%%%%%%% stop panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_stop = uipanel('Title', '', 'Units', 'Normalized', 'Position', [0.2 0.0 0.6 0.1]);
            
            obj.button_stop = GraphicButton(obj.panel_stop, [0 0 1 1], obj.cal, '', 'custom');
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
            
            if ~isempty(obj.cal.gain)
                
                obj.button_calc_gain.String = ['gain= ' num2str(obj.cal.gain)];
                
            else
                
                obj.button_calc_gain.String = 'calcGain';
                
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
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% stop panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function callback_stop(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: stop'); end
            
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