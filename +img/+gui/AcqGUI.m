classdef AcqGUI < handle
    
    properties 
        
        owner@img.Acquisition; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        menus = {};
        
        latest_error = '';
        latest_warning = '';
        latest_message = ''; 
        
        font_size = 13;
        big_font_size = 16;
        edit_font_size = 12;
        small_font_size = 10;
        
        color_on = [0 0.3 1];
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        menu_options;
        
        menu_objects;
        
        panel_controls;
        
        panel_contrast;
        
        panel_sidebar;
    
%         panel_objects;
%         
%         panel_save;
%         
%         panel_sync;
        
        panel_run;
        button_run;
        
        panel_info;
        
        panel_image;
        
        button_reset_axes;
        button_unlock;
        button_batch_num;
        button_time_left;
        button_gb_left; 
        
        button_show_what;
        button_gray;
%         button_src_status;
        button_flip;
        button_show_switch;
        button_num_stars;
        axes_image;
        
%         panel_close;
%         button_close;
        
    end
    
    properties (Hidden=true)
              
        version = 1.01;
        
    end
            
    methods % constructor
       
        function obj = AcqGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'img.Acquisition')
                
                if obj.debug_bit>1, fprintf('AcqGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input an Acquisition object to constructor of AcqGUI!');
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
            
            obj.fig = util.plot.FigHandler('Acquisition');
            
            obj.fig.bottom = 5;
            obj.fig.height = 20;
            obj.fig.width = 30;
            movegui(obj.fig.fig, 'center');
            obj.fig.clear;
            
            
            %%%%%%%%%%%%%%%%%%%%%%% MENUS %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % MenuItem(parent, text, type, variable, tooltip, separator)
            % obj.addButton(text, type, variable, tooltip, separator)
            % menu types: menu, toggle, push, input, input_text, info, custom
            
            obj.menu_options = MenuItem(obj, '&Options', 'menu');
            
            obj.menu_options.addButton('menu_parameters', '&Parameters', 'menu');
            obj.menu_options.menu_parameters.addButton('button_expT', 'exp&T= %6.3f s', 'input', 'expT', 'Exposure time for each frame (seconds)');
            obj.menu_options.menu_parameters.addButton('button_frame_rate', '&Frame_rate= %6.3f Hz', 'input', 'frame_rate', 'frame rate (cadence) given in Hz');
            
            obj.menu_options.addButton('menu_header', '&Header', 'menu');
            obj.menu_options.menu_header.addButton('button_object', '&Object: %s', 'input_text', 'run_name', 'object or field name (used for naming the data folder');
            obj.menu_options.menu_header.addButton('button_ra', '&RA: %s', 'input_text', 'head.RA', 'object right ascention (hours)');
            obj.menu_options.menu_header.addButton('button_dec', '&Dec: %s', 'input_text', 'head.Dec', 'object declination (degrees)');
            
            
            obj.menu_options.addButton('menu_find_stars', '&Find stars', 'menu');
            obj.menu_options.menu_find_stars.addButton('button_astrometry', '&Astrometry', 'toggle', 'use_astrometry', 'use astrometry (with GAIA DR2) to get star positions/magnitudes');
            obj.menu_options.menu_find_stars.addButton('button_mextractor', '&Mextractor', 'toggle', 'use_mextractor', 'use mextractor to find the star positions');
            obj.menu_options.menu_find_stars.addButton('button_arbitrary', '&Arbitrary pos', 'toggle', 'use_arbitrary_pos', 'set random cutout positions (debugging only!)');
            obj.menu_options.menu_find_stars.addButton('button_threshold', '&Threshold', 'input', 'detect_thresh', 'threshold for detection in stack image (noise rms)');
            obj.menu_options.menu_find_stars.addButton('button_bad_pixels', '&Bad pixels', 'toggle', 'use_remove_bad_pixels', 'remove bad pixels before finding stars');
            obj.menu_options.menu_find_stars.addButton('button_saturated', '&Remove saturated', 'toggle', 'use_remove_saturated', 'remove any stars with any saturated pixels');
            obj.menu_options.menu_find_stars.addButton('button_sat_value', '&Saturation value', 'input', 'saturation_value', 'what pixel value is considered saturated');
            obj.menu_options.menu_find_stars.addButton('button_min_temp', '&Minimal temperature', 'input', 'min_star_temp', 'take only stars with temperature above this value (Kelvin). Works only with astrometry');
            
            obj.menu_options.addButton('menu_cutouts', '&Cutouts', 'menu');
            obj.menu_options.menu_cutouts.addButton('button_num_stars', '&Num stars', 'input', 'num_stars', 'Maximum number of cutouts around stars');
            obj.menu_options.menu_cutouts.addButton('button_cut_size', '&Cutout size', 'input', 'cut_size', 'Size of star cutouts (pixels)');
            obj.menu_options.menu_cutouts.addButton('button_egdes', '&Edge distance', 'input', 'avoid_edges', 'Distance from edges to avoid finding stars (pixels)');
            obj.menu_options.menu_cutouts.addButton('button_num_bgs', 'Num &Backgrounds', 'input', 'num_backgrounds', 'Number of background cutouts');
            obj.menu_options.menu_cutouts.addButton('button_size_bgs', '&Background size', 'input', 'cut_size_bg', 'Size of background cutouts (pixels)');
            obj.menu_options.menu_cutouts.addButton('button_check', '&Position check', 'toggle', 'use_check_positions', 'Toggle checking flux and realigning or quitting run if no stars are found');
            obj.menu_options.menu_cutouts.addButton('button_adjust', '&Adjust cutouts', 'toggle', 'use_adjust_cutouts', 'Toggle cutout adjustments: reposition cutouts based on centroids');
            obj.menu_options.menu_cutouts.addButton('button_lock_adjust', '&Lock adjust', 'toggle', 'use_lock_adjust', 'Toggle cutout adjust lock: all cutouts would move together');
            
            obj.menu_options.addButton('menu_photometry', '&Photometry', 'menu');
            obj.menu_options.menu_photometry.addButton('button_use_simple', '&Simple photometry', 'toggle', 'use_simple_photometry', 'just sum the cutouts, not using Photometry/Lightcurve objects');
            obj.menu_options.menu_photometry.addButton('button_use_store_photometry', 'S&tore photometry', 'toggle', 'use_store_photometry', 'keep a copy of the photometric products in the Lightcurve object');
            obj.menu_options.menu_photometry.addButton('button_use_save_photometry', 'Sa&ve photometry', 'toggle', 'use_save_photometry', 'save the flux and other photometric products in the HDF5 files');
            obj.menu_options.menu_photometry.addButton('button_model_psf', '&PSF model', 'toggle', 'use_model_psf', 'stack the cutouts and fit it to a PSF model');
            obj.menu_options.menu_photometry.addButton('button_num_cutouts', '&Num cutouts', 'input', 'num_phot_cutouts', 'How many stars/cutouts to input to photomery and store in lightcurves');
            
            obj.menu_options.addButton('menu_sync', '&Synchronization', 'menu');
            obj.menu_options.menu_sync.addButton('button_use_sync', '&Use sync', 'toggle', 'use_sync', 'use the sync between the two computers');
            obj.menu_options.menu_sync.addButton('button_ignore_manager', '&Ignore manager', 'toggle', 'use_ignore_manager', 'when true, Acquisition will ignore all commands from Manager');
            obj.menu_options.menu_sync.addButton('button_stop', 'Sync &stop', 'toggle', 'use_sync_stop', 'when true, will stop camera when given command from Manager, when dome is close, when telescope is not tracking');
            obj.menu_options.menu_sync.addButton('button_guiding', '&Guiding', 'toggle', 'use_autoguide', 'use this to pass position data back to mount');
            
            obj.menu_options.addButton('menu_deflate', '&Deflate', 'menu');
            obj.menu_options.menu_deflate.addButton('button_autodeflate', '&Autodeflate', 'toggle', 'use_autodeflate', 'turn on the deflate automatically in the morning');


            
            
            % add latest_error/warning to MenuItem
            
            obj.menu_objects = MenuItem(obj, 'O&bjects', 'menu');
            obj.menu_objects.addButton('button_camera', '&Camera', 'push', 'cam', 'Start the camera GUI'); 
            obj.menu_objects.addButton('button_header', '&Header', 'push', 'head', 'Start the header object GUI'); 
            obj.menu_objects.addButton('button_calibration', 'C&alibration', 'push', 'cal', 'Start the Calibration object GUI'); 
            obj.menu_objects.addButton('button_clipper', 'C&lipper', 'push', 'clip', 'Start the Clipper GUI'); 
            obj.menu_objects.addButton('button_photometry', '&Photometry', 'push', 'phot', 'Start the Photometry GUI'); 
            obj.menu_objects.addButton('button_lightcurves', '&Lightcurves', 'push', 'lightcurves', 'Start the Lightcurves GUI'); 
            obj.menu_objects.addButton('button_buffers', '&Buffers', 'push', 'buf', 'Start the BufferWheel GUI'); 
            obj.menu_objects.addButton('button_focuser', '&Focuser', 'custom', '', 'Start the Focuser GUI'); 
            obj.menu_objects.addButton('button_deflator', '&Deflator', 'push', 'deflator', 'Start the Deflator GUI'); 

            obj.menu_objects.button_focuser.Callback = @obj.callback_focuser;
            
            %%%%%%%%%%%%%%% BUTTONS %%%%%%%%%%%%%%%%%
            
            N_left = 15;
            N_right = 15;
            W_left = 0.2;
            W_right = 0.2;
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            pos = N_left;
            N = 10;
            pos = pos - N;
            
            obj.panel_controls = GraphicPanel(obj.owner, [0 pos/N_left W_left N/N_left], 'controls');
            obj.panel_controls.number = N;
            
            obj.panel_controls.addButton('input_runtime', 'total_runtime', 'input', 'run time= ', '', '', 0.6, '', '', 'Total time for the next run (in minutes)');
            obj.panel_controls.addButton('chooser_units', '', 'custom', '', '', 'small', 0.4); 
            obj.panel_controls.addButton('button_slow_mode', 'setupSlowMode', 'push', 'Slow mode', '', '', 0.5, '', '', 'Setup low-cadence observations with T=3s, single image files, 500 stars max');
            obj.panel_controls.addButton('button_fast_mode', 'setupFastMode', 'push', 'Fast mode', '', '', 0.5, '', '', 'Setup high-cadence observations with T=30ms, f=25Hz, 100 images per batch/file');
            obj.panel_controls.addButton('button_single', 'single', 'push', 'Take single batch', '', '', [], '', '', 'Take a single image and show it on screen. Does not save the image to file');
            obj.panel_controls.addButton('button_live', 'startLiveView', 'push', 'Start live view', '', '', [], '', '', 'Open the camera GUI and start the live-view video mode. Does not save any images to file');
            obj.panel_controls.addButton('button_auto_focus', 'cam.autofocus', 'push', 'Auto-focus', '', '', 0.7, '', '', 'Start a focus-run. Does not save any images to file');
            obj.panel_controls.addButton('button_manual_focus', '', 'custom', 'manual', '', '', 0.3, '', '', 'open the focus GUI for manual focusing');
            obj.panel_controls.addButton('button_preview', '', 'custom', 'Preview run', '', '', [], '', '', 'Start a preview run, checking full acquisition pipeline. Does not save any images to file');
            
            obj.panel_controls.margin = [0.02 0.01];
            obj.panel_controls.make;
            
            obj.panel_controls.chooser_units.control.Style = 'popupmenu';
            obj.panel_controls.chooser_units.control.String = {'seconds', 'minutes', 'hours', 'batches', 'frames'};
            
            obj.panel_controls.chooser_units.Callback = @obj.callback_units_picker;
            obj.panel_controls.button_preview.Callback = @obj.callback_preview;
            obj.panel_controls.button_manual_focus.Callback = @obj.callback_focuser;
            
            %%%%%%%%%%% panel contrast %%%%%%%%%%%%%%%
            
            N = 5;
            pos = pos - N;
            obj.panel_contrast = util.plot.ContrastLimits(obj.axes_image, obj.fig.fig, [0 pos/N_left W_left N/N_left], 1, [0.02 0.01]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RIGHT SIDE %%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%% panel sidebar %%%%%%%%%%%%%%%%
            pos = N_right;
            N = 15;
            pos = pos - N;

            obj.panel_sidebar = GraphicPanel(obj.owner, [1-W_right pos/N_right W_right N/N_right], '');
            obj.panel_sidebar.addButton('button_info', 'info_long', 'info', ' ', '', 'small');
            obj.panel_sidebar.margin = [0.02 0.01];
            obj.panel_sidebar.make;
            obj.panel_sidebar.button_info.control.Style = 'text';
            obj.panel_sidebar.button_info.control.HorizontalAlignment = 'left';
            
            
            %%%%%%%%%%% panel run %%%%%%%%%%%%%%%%%%
            
            obj.panel_run = GraphicPanel(obj.owner, [W_left 0 (1-W_left-W_right) 0.1]);
            
            obj.panel_run.addButton('button_run', '', 'custom', 'RUN', '', '', 0.8);
            obj.panel_run.addButton('button_continue', 'run', 'push', 'Continue', '', '', 0.2);
            obj.panel_run.margin = [0.005 0.1];
            obj.panel_run.make;
            obj.panel_run.button_run.Callback = @obj.callback_run;
            obj.panel_run.button_continue.Callback = @obj.callback_continue;
            obj.panel_run.button_run.Tooltip = 'Start/stop the run';
            obj.panel_run.button_continue.Tooltip = 'Continue the current run without reseting';
            
            %%%%%%%%%%% panel info %%%%%%%%%%%%%%%%%%
            
            obj.panel_info = GraphicPanel(obj.owner, [W_left 0.9 1-W_left-W_right 0.1], 'info', 1);
            obj.panel_info.addButton('button_message', '', 'custom', '', '', 'small'); 
            obj.panel_info.addButton('button_info', 'info_short', 'info', ' ', '', 'small'); 
            
            obj.panel_info.margin = [0.005 0.03];
            obj.panel_info.make;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [W_left 0.1 1-W_left-W_right 0.8]);
                        
            obj.makeAxes;
            
            obj.button_show_what = GraphicButton(obj.panel_image, [0 0.00 0.15 0.05], obj.owner, 'show_what', 'picker', 'full', '', 'small');
            obj.button_show_what.Callback = @obj.callback_show_what;
            obj.button_show_what.String = obj.owner.show_what_list;
            
            obj.button_gray = GraphicButton(obj.panel_image, [0.15 0.00 0.1 0.05], obj.owner, 'use_show_gray', 'toggle', 'color', 'gray', 'small'); 
            obj.button_gray.Tooltip = 'show the images in color or in gray-scale'; 
            
            obj.button_flip = GraphicButton(obj.panel_image, [0.9 0.00 0.1 0.05], obj.owner, 'use_flip', 'toggle', 'no flip', 'flip', 'small');
            obj.button_show_switch = GraphicButton(obj.panel_image, [0.8 0.00 0.1 0.05], obj.owner, 'use_show', 'toggle', 'no show', 'show', 'small');
                        
            obj.button_num_stars = GraphicButton(obj.panel_image, [0.00 0.95 0.1 0.05], obj.owner, 'display_num_rect_stars', 'input', 'rect= ', '', 'small');
            
            obj.button_unlock = GraphicButton(obj.panel_image, [0.9 0.95 0.1 0.05], obj.owner, 'unlock', 'push', 'unlock');
            
            obj.button_show_what.Tooltip = 'Display images, stack, or stack processed (dark, flat, background removed)';
            obj.button_flip.Tooltip = 'Flip the image 180 degrees (for viewing beyond the meridian)';
            obj.button_unlock.Tooltip = 'Manually remove all locks left over from run-time errors. Make sure camera is stopped!!';
            
            obj.update;
            
        end
            
        function makeAxes(obj, ~, ~)
            
            delete(obj.axes_image);
            
            obj.axes_image = axes('Parent', obj.panel_image);
            
            obj.panel_contrast.ax = obj.axes_image;
%             colorbar(obj.axes_image);
            axis(obj.axes_image, 'image');
            
            obj.panel_contrast.ax = obj.axes_image;
            
            if ~isempty(obj.owner.images)
                obj.owner.show;
            end
            
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
            
            obj.panel_contrast.font_size = obj.font_size;
            obj.panel_contrast.edit_font_size = obj.edit_font_size;
            obj.panel_contrast.update;
            
            if obj.owner.is_slow_mode
                obj.panel_controls.button_slow_mode.ForegroundColor = obj.color_on;
            else
                obj.panel_controls.button_slow_mode.ForegroundColor = 'black';
            end
            
            if obj.owner.is_fast_mode
                obj.panel_controls.button_fast_mode.ForegroundColor = obj.color_on;
            else
                obj.panel_controls.button_fast_mode.ForegroundColor = 'black';
            end
            
            if obj.owner.brake_bit && (obj.owner.is_running || obj.owner.is_running_single) % started a run but not yet took off the brake bit
                obj.panel_run.button_run.Enable = 'off';
            else

                if obj.owner.brake_bit
                    obj.panel_run.button_run.Enable = 'on';
                   obj.panel_run.button_run.String = 'START NEW RUN';                
                else % still running
                    obj.panel_run.button_run.Enable = 'on';
                    obj.panel_run.button_run.String = 'STOP';
                end

            end
            
            if obj.owner.batch_counter>0 && obj.owner.brake_bit && ~(obj.owner.is_running || obj.owner.is_running_single)
                obj.panel_run.button_continue.Enable = 'on';
            else
                obj.panel_run.button_continue.Enable = 'off';
            end
            
            if obj.owner.is_running || obj.owner.is_running_single
                
                for ii = 1:length(obj.panel_controls.buttons)
                    obj.panel_controls.buttons(ii).Enable = 'off';
                end
                
            else
                   
                for ii = 1:length(obj.panel_controls.buttons)
                    obj.panel_controls.buttons(ii).Enable = 'on';
                end
                
            end
            
            obj.panel_controls.chooser_units.control.Value = find(strcmpi(obj.panel_controls.chooser_units.control.String, obj.owner.runtime_units));
            obj.button_show_what.control.Value = find(strcmpi(obj.button_show_what.control.String, obj.owner.show_what));

            obj.latest_warning = lastwarn;

            display_string = ''; 
            
            if isempty(display_string) && ~isempty(obj.latest_error)
                display_string = regexprep(obj.latest_error, '\n*|\s{2,}', '... ');
            end
            
            % I wanted to stop showing these warnings, for now...
%             if isempty(display_string) && ~isempty(obj.latest_warning)
%                 display_string = ['Warning: ' regexprep(obj.latest_warning, '\n*|\s{2,}', '... ')];
%             end
            
            if isempty(display_string) && ~isempty(obj.latest_message)
                display_string = obj.latest_message; 
            end
            
            obj.panel_info.button_message.String = display_string;    
            
            if obj.owner.use_show
                obj.owner.show;
            end

            drawnow;
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_image) && isvalid(obj.panel_image);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_units_picker(obj, hndl, ~)
            
            if obj.debug_bit>1, disp('callback: units picker'); end
            
            total_seconds = obj.owner.total_runtime.*obj.owner.convertRuntimeToSeconds;
            
            obj.owner.runtime_units = hndl.String{hndl.Value};
            
            obj.owner.total_runtime = total_seconds./obj.owner.convertRuntimeToSeconds;
            
            obj.update;
            
        end
        
        function callback_run(obj, ~, ~)
            
            if obj.owner.brake_bit
                if obj.debug_bit>1, disp('callback: run'); end
                
                obj.latest_error = '';
                obj.latest_warning = '';
                obj.latest_message = '';
                lastwarn('');
                                
                try
                    obj.owner.use_save = 1;
                    obj.owner.run('reset', 1);
                catch ME
                    obj.latest_error = util.text.eraseTags(ME.getReport());
                    obj.update;
                    rethrow(ME);
                end
                
            else
                if obj.debug_bit>1, disp('callback: stop'); end            
                obj.owner.brake_bit = 1;
                obj.update;
            end
            
        end
        
        function callback_preview(obj, ~, ~)
            
            if obj.debug_bit>1, disp('callback: preview'); end

            obj.latest_error = '';
            obj.latest_warning = '';
            lastwarn('');

            try
                obj.owner.use_save = 0;
                obj.owner.run('reset', 1, 'num_batches', 30);
            catch ME
                obj.latest_error = util.text.eraseTags(ME.getReport());
                obj.update;
                rethrow(ME);
            end

        end
        
        function callback_continue(obj, ~, ~)
            
            if obj.debug_bit>1, disp('callback: continue'); end

            obj.latest_error = '';
            obj.latest_warning = '';
            lastwarn('');

            try
                obj.owner.run;
            catch ME
                obj.latest_error = util.text.eraseTags(ME.getReport());
                obj.update;
                rethrow(ME);
            end


        end
        
        function callback_num_files(obj, ~, ~)
            
            if obj.debug_bit>1, disp('callback: num_files'); end
            
            if ~isinf(obj.owner.num_files)
                obj.owner.num_batches = obj.owner.num_files;
            end
            
            obj.update;
        
        end
        
        function callback_show_what(obj, hndl, ~)
            
            if obj.debug_bit>1, disp('callback: show_what'); end
            
%             obj.owner.show_what_index = hndl.Value;
%             obj.owner.show_what = obj.owner.show_what_list{obj.owner.show_what_index};

            obj.owner.show_what = hndl.String{hndl.Value};

            obj.owner.show;
            obj.panel_contrast.autodyn;
            
            obj.update;
            
        end
        
        function callback_focuser(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: focuser'); end
        
            obj.owner.cam.focuser.makeGUI;
            
            obj.update;
            
        end
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end