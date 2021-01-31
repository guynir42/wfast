classdef LightGUI < handle
    
    properties 
        
        owner@img.Lightcurves; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 16;
        edit_font_size = 14;
        small_font_size = 12;
        
        color_on = [0 0.3 1];
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_methods;
        
        panel_control;
        
        panel_process;
        
        panel_psd;
        
        panel_display;
        
        panel_stats;
        
        panel_info;
        
        panel_close;
        button_close;
        
        panel_image;
        button_reset_axes;
        button_log_scale
        button_cut_number; 
        button_datetime;
        axes_image;
    
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = LightGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'img.Lightcurves')
                
                if obj.debug_bit>1, fprintf('LightGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input an img.Lightcurves to constructor of LightGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            
            obj.fig = util.plot.FigHandler('lightcurves');
            obj.fig.clear;
            
            obj.fig.bottom = 5;
            obj.fig.height = 20;
            obj.fig.width = 32;
            movegui(obj.fig.fig, 'center');
            
            %%%%%%%%%%% panel methods %%%%%%%%%%%%%%%
            
%             obj.panel_methods = GraphicPanel(obj.owner, [0 0.9 0.2 0.1], 'methods'); 
%             obj.panel_methods.addButton('button_signal_method', 'signal_method', 'info', ' ', '', '', 0.5);
%             obj.panel_methods.addButton('button_background_method', 'background_method', 'info', ' ', '', '', 0.5);
%             obj.panel_methods.make;
            
            
            N = 19;
            pos = N;

            %%%%%%%%%%% panel control %%%%%%%%%%%%%%%
            
            pos = pos - 1;
            obj.panel_control = GraphicPanel(obj.owner, [0 pos/N 0.2 1/N], 'control'); 
            obj.panel_control.number = 1;
            obj.panel_control.addButton('button_preallocate', 'use_preallocate', 'toggle', 'preallocate', 'preallocate', '', 0.5, obj.color_on);
            obj.panel_control.addButton('button_double_up', 'use_double_up', 'toggle', 'double up', 'double up', '', 0.5, obj.color_on);
            obj.panel_control.margin = [0.01 0.02];
            obj.panel_control.make;
            
            %%%%%%%%%%% panel process %%%%%%%%%%%%%%%%
            
            pos = pos - 8;
            obj.panel_process = GraphicPanel(obj.owner, [0 pos/N 0.2 8/N], 'process');
            obj.panel_process.number = 5;
            obj.panel_process.addButton('button_background', 'use_subtract_backgrounds', 'toggle', 'sub b/g', 'sub b/g', '', 0.5, obj.color_on);
            obj.panel_process.addButton('button_back_median', 'use_background_median', 'toggle', 'b/g median', 'b/g median', '', 0.5, obj.color_on);
%             obj.panel_process.addButton('button_flagged', 'use_skip_flagged', 'toggle', 'skip flag', 'skip flag', '', 0.5, obj.color_on);
            obj.panel_process.addButton('button_outliers', 'use_remove_outliers', 'toggle', 'outliers', 'outliers', '', 0.5, obj.color_on);
            obj.panel_process.addButton('input_outlier_sigma', 'outlier_sigma', 'input', 'sig= ', '', '', 0.5);
            obj.panel_process.addButton('button_bad_times', 'use_skip_bad_times', 'toggle', 'bad times', 'bad times', '', 0.5, obj.color_on);
            obj.panel_process.addButton('input_fraction', 'bad_times_fraction', 'input', 'frac= ', '', '', 0.5);
            obj.panel_process.addButton('button_airmass', 'use_airmass_correction', 'toggle', 'airmass', 'airmass', '', 0.5, obj.color_on);
            obj.panel_process.addButton('button_placeholder', '', 'custom', ' ', ' ', '', 0.5, obj.color_on);
            obj.panel_process.addButton('button_width', 'use_psf_correction', 'toggle', 'PSF width', 'PSF width', '', 0.5, obj.color_on);
            obj.panel_process.addButton('button_offset_fit', 'use_offset_fit', 'toggle', 'offsets', 'offsets', '', 0.5, obj.color_on);
            obj.panel_process.addButton('button_zero_point', 'use_zero_point', 'toggle', 'zero point', 'zero point', '', 0.5, obj.color_on);
            obj.panel_process.addButton('button_zp_order', 'zero_point_spatial_order', 'input', 'zp order= ', '', '', 0.5);
            obj.panel_process.addButton('button_sysrem', 'use_sysrem', 'toggle', 'sysrem', 'sysrem', '', 0.5, obj.color_on);
            obj.panel_process.addButton('input_sysrem_iter', 'sysrem_iterations', 'input', 'iter= ', '', '', 0.5);
            obj.panel_process.addButton('button_self_exclude', 'use_self_exclude', 'toggle', 'self', 'self', '', 0.5, obj.color_on); 
            
            
%             obj.panel_process.addButton('button_zavitzky_golay', 'use_savitzky_golay', 'toggle', 'savitzky golay', 'savitzky golay', '', 0.5, obj.color_on); 
%             obj.panel_process.addButton('input_sg_order', 'sg_order', 'input', 'SG order= ', '', '', 0.5); 
%             obj.panel_process.addButton('input_sg_length', 'sg_length', 'input', 'SG len= ', '', '', 0.5);
%             obj.panel_process.addButton('input_missing', 'missing_method', 'input_text', 'miss= ', '', '', 0.5); 

%             obj.panel_process.addButton('button_psf_correction', 'use_psf_correction', 'toggle', 'PSF width', 'PSF width', '', 0.5, obj.color_on);
%             obj.panel_process.addButton('button_polynomial', 'use_polynomial', 'toggle', 'polynomial', 'polynomial', '', 0.5, obj.color_on);
            
            obj.panel_process.margin = [0.01 0.01];
            obj.panel_process.make;

            %%%%%%%%%%% panel display %%%%%%%%%%%%%%%%
            
            pos = pos - 6;
            obj.panel_display = GraphicPanel(obj.owner, [0 pos/N 0.2 6/N], 'display');
            obj.panel_display.number = 6;
            obj.panel_display.addButton('button_index', 'index_flux', 'input', 'idx= ','', '', 0.5, '', '', 'which aperture type to display'); 
            obj.panel_display.addButton('button_show_what', 'show_what', 'picker', obj.owner.show_what, '', '', 0.7); 
            obj.panel_display.addButton('button_flux_type', 'show_flux_type', 'picker', obj.owner.show_flux_type, '', '', 0.3); 
            obj.panel_display.addButton('input_num_stars', 'show_num_stars', 'input', 'Nstars= ');
            obj.panel_display.addButton('input_show_indices', 'show_indices', 'input', 'idx= '); 
            obj.panel_display.addButton('button_smooth', 'use_smooth', 'toggle', 'smooth', 'smooth', '', 0.5, obj.color_on);
            obj.panel_display.addButton('input_smooth', 'smooth_interval', 'input', ' ', '', '', 0.5);
            obj.panel_display.addButton('button_spawn', 'spawnFigure', 'push', 'spawn', '', '', 0.5); 
            obj.panel_display.margin = [0.01 0.02];
            obj.panel_display.make;
            
            obj.panel_display.button_show_what.String = obj.owner.show_what_list;
            obj.panel_display.button_show_what.Callback = @obj.callback_show_what;
            
            obj.panel_display.button_flux_type.String = obj.owner.show_flux_type_list;
            obj.panel_display.button_flux_type.Callback = @obj.callback_flux_type;
            
            %%%%%%%%%%% panel psd %%%%%%%%%%%%%%%%%%%%
                        
            pos = pos - 1;
            obj.panel_psd = GraphicPanel(obj.owner, [0 pos/N 0.2 1/N], 'display');
            obj.panel_psd.number = 1;
            obj.panel_psd.addButton('button_welch', 'use_welch', 'toggle', 'welch', 'welch', '', 0.5, obj.color_on);
            obj.panel_psd.make;
            
            %%%%%%%%%%% panel statistics %%%%%%%%%%%%%
            
            pos = pos - 2;
            obj.panel_stats = GraphicPanel(obj.owner, [0 pos/N 0.2 2/N], 'statistics');
            obj.panel_stats.number = 2;
            obj.panel_stats.addButton('input_jump', 'sampling_jump', 'input', 'jump= ', '', '', 0.5);
            obj.panel_stats.addButton('input_rms_points', 'num_points_rms', 'input', 'rms pts= ', '', '', 0.5); 
%             obj.panel_stats.addButton('button_bins', 'use_show_re_bins', 'toggle', 'show mag', 'show bin', '', 0.5);
            obj.panel_stats.addButton('input_mag_limit', 'show_mag_limit', 'input', 'mag lim= ', '', '', 0.5);
            obj.panel_stats.addButton('button_fits', 'use_show_bin_fits', 'toggle', 'show fits', 'show fits', '', 0.5, obj.color_on);
            obj.panel_stats.make;
            
            %%%%%%%%%%% panel info %%%%%%%%%%%%%%%%%%%
            
            obj.panel_info = GraphicPanel(obj.owner, [0.2 0.9 0.8 0.1], 'info', 0);
            obj.panel_info.number = 1;
            obj.panel_info.addButton('button_info', 'print_pars','info'); 
            obj.panel_info.make;
            obj.panel_info.button_info.font_size = 'edit';
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 0.0 0.8 0.9]);
            
            obj.makeAxes;
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.9 0.95 0.1 0.05], obj.owner, '', 'custom', 'reset axis', '', 'small');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            
            obj.button_log_scale = GraphicButton(obj.panel_image, [0.0 0.95 0.1 0.05], obj.owner, 'use_show_log', 'toggle', 'linear', 'log', 'small');
            
            obj.button_cut_number = GraphicButton(obj.panel_image, [0.0 0.0 0.2 0.05], obj.owner, '', 'custom', 'star: '); 
            
            obj.button_datetime = GraphicButton(obj.panel_image, [0.8 0.0 0.2 0.05], obj.owner, 'use_show_datetime', 'toggle', 'seconds', 'datetime'); 
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
                        
            obj.panel_close = uipanel('Position', [0 0 0.2 1/N]);
            obj.button_close = GraphicButton(obj.panel_close, [0.1 0.1 0.8 0.8], obj.owner, '', 'custom', 'CLOSE GUI');
            obj.button_close.Callback = @obj.callback_close;
            
            obj.update;
            
        end
            
        function makeAxes(obj, ~, ~)
            
            delete(obj.axes_image);
            
            obj.axes_image = axes('Parent', obj.panel_image);
            
            obj.owner.last_shown = '';
            obj.owner.last_xlim = [];
            obj.owner.last_ylim = [];
            
            obj.update;
            
        end
        
        function adjustAxes(obj, ~, ~)
            
            axis(obj.axes_image, 'normal'); 
            
        end
        
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            obj.button_cut_number.String = 'star= ';
            
            obj.owner.show;
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_show_what(obj, hndl, ~)
            
            if obj.debug_bit>1, disp('callback: show_what'); end
            
            obj.owner.show_what = hndl.String{hndl.Value};

            obj.update;
            
        end
        
        function callback_flux_type(obj, hndl, ~)
            
            if obj.debug_bit>1, disp('callback: flux_type'); end
            
            obj.owner.show_flux_type = hndl.String{hndl.Value};

            obj.update;
            
        end
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end