classdef ProcGUI < handle
    
    properties 
        
        owner@img.Processor; % link back to containg object

        fig@util.plot.FigHandler;
        
        buttons = {};
        menus = {};
        panels = {};
        
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
        
        menu_options;
        menu_coadds;
        menu_stars;
        menu_cutouts;
        menu_save;
        menu_fits;
        
        panel_controls;
        
        panel_objects;
        
        panel_contrast;
        
        panel_info;
        panel_stop;
        
        panel_close;
        button_close;
        
        panel_image;
        button_reset_axes;
        button_display_text;
        input_num_rect;
        button_expt;
        button_indices;
        panel_legend;
        panel_buffers;
        
        axes_image;
    
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = ProcGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'img.Processor')
                
                if obj.debug_bit>1, fprintf('ProcGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input an img.Processor to constructor of ProcGUI!');
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
            obj.panels = {}; 
            
            obj.fig = util.plot.FigHandler('processor');
            obj.fig.clear;
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 30;
            obj.fig.center;
%             obj.fig.maximize;
            
            
            %%%%%%%%%%%%%%%%%%%%%%% MENUS %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % MenuItem(parent, text, type, variable, tooltip, separator)
            % obj.addButton(name, text, type, variable, tooltip, separator)
            % menu types: menu, toggle, push, input, input_text, info, custom
            
            obj.menu_options = MenuItem(obj, '&Options', 'menu'); 
            obj.menu_options.addButton('button_calibration', '&Load calibration', 'toggle', 'pars.use_auto_load_cal', 'find the appropriate calibration file based on observation date'); 
            obj.menu_options.addButton('button_astrometry', '&Astrometry', 'toggle', 'pars.use_astrometry', 'run astrometry and save a new catalog file'); 
                        
            obj.menu_coadds = MenuItem(obj, 'Co-&Adds', 'menu'); 
            obj.menu_coadds.addButton('button_coadds', '&Use coadds', 'toggle', 'pars.use_coadds', 'coadd some images before finding stars, etc.'); 
            obj.menu_coadds.addButton('input_coadd_size', 'Coadd &Size', 'input', 'pars.coadd_size', 'how many images to put into each coadd'); 
            
            
            obj.menu_stars = MenuItem(obj, '&Find Stars', 'menu'); 
            obj.menu_stars.addButton('input_num_stars', '&Num stars', 'input', 'pars.num_stars', 'maximum number of stars to track for photometry'); 
            obj.menu_stars.addButton('input_threshold', '&Threshold', 'input', 'pars.threshold', 'how much above the background noise the stars need to be to be detected'); 
            obj.menu_stars.addButton('input_saturation', '&Saturation', 'input', 'pars.saturation', 'if any pixels in the star region surpass this value, flag the star as saturated'); 
            obj.menu_stars.addButton('button_bg_map', '&Background map', 'toggle', 'pars.use_bg_map', 'interpolate the background to a different value in each point in the image'); 
            
            obj.menu_stars.addButton('menu_bad_matches', 'Bad &Matches', 'menu');
            obj.menu_stars.menu_bad_matches.addButton('button_remove', '&Remove bad', 'toggle', 'pars.use_remove_bad_matches', 'stars without a GAIA match are removed, unless they are very bright');
            obj.menu_stars.menu_bad_matches.addButton('input_snr', '&Min S/N', 'input', 'pars.bad_match_min_snr', 'stars above this S/N value are not removed, even without a GAIA match (e.g., asteroids)'); 
            
            obj.menu_stars.addButton('menu_kernel', '&Kernel estimate', 'menu');
            obj.menu_stars.menu_kernel.addButton('input_initial', '&Initial width', 'input', 'pars.initial_guess_psf_width', 'initial guess of PSF width to find a small sample of stars used to measure the PSF width and filter kernel'); 
            obj.menu_stars.menu_kernel.addButton('input_num_stars', '&Num stars', 'input', 'pars.num_stars_filter_kernel', 'number of stars used for the first detection run, from which we calculate the PSF width and filter kernel');
            
            obj.menu_cutouts = MenuItem(obj, '&Cutouts', 'menu'); 
            obj.menu_cutouts.addButton('input_cut_size', '&Cut size', 'input', 'pars.cut_size', 'size of the cutout stamps (pixels)'); 
            obj.menu_cutouts.addButton('input_lock', '&Lock adjust', 'input_text', 'pars.lock_adjust', 'choose to lock the cutout adjustments to "all", "none" or "stars" only (those with GAIA match)'); 
            obj.menu_cutouts.addButton('button_use_check_flux', 'Check &Flux', 'toggle', 'pars.use_check_flux', 'check the flux each batch and call off analysis if the flux is lost multiple times'); 
            obj.menu_cutouts.addButton('input_num_failed', '&Max failed files', 'input', 'pars.max_failed_files', 'how many files in a row must fail to call off the analysis'); 
            
            obj.menu_save = MenuItem(obj, '&Save', 'menu'); 
            obj.menu_save.addButton('button_save', '&Save results', 'toggle', 'pars.use_save_results', 'automatically save lightcurves and other products to disk at end of run'); 
            obj.menu_save.addButton('button_overwrite', '&Overwrite', 'toggle', 'pars.use_overwrite', 'if processor folder already exists for today, silently overwrite the contents with a new analysis'); 
            obj.menu_save.addButton('input_output', 'Output &Folder', 'input_text', 'pars.output_folder', 'give a different name to the output folder, instead of just processor_YYYY-MM-DD'); 
            obj.menu_save.addButton('button_now', 'Save &Now', 'push', 'saveResults', 'save the results of the analysis into a processor folder right now', 1); 
            
            obj.menu_fits = MenuItem(obj, 'FITS', 'menu'); 
            obj.menu_fits.addButton('button_save', '&Save FITS', 'toggle', 'pars.fits.use_save', 'save the calibrated images into FITS files');
            obj.menu_fits.addButton('button_use_roi', 'Use &ROI', 'toggle', 'pars.fits.use_roi', 'cut out a Region Of Interest before saving FITS files');
            obj.menu_fits.addButton('input_size', 'ROI &Size', 'input', 'pars.fits.roi_size', 'define the size of the Region Of Interest as scalar or [height,width] vector');
            obj.menu_fits.addButton('input_coords', 'ROI &Coordinates', 'input_text', 'pars.fits.roi_coordinates', 'select the position of the Region Of Interest, use "header" to auto-center on object RA/Dec');
            obj.menu_fits.addButton('button_flip', 'Flip &Image', 'toggle', 'pars.fits.use_flip', 'flip the FITS images 180 degrees after getting ROI');
            obj.menu_fits.addButton('input_dir', '&Directory', 'input_text', 'pars.fits.directory', 'choose the name of the FITS files subfolder (or absolute path to any folder)');
            obj.menu_fits.addButton('input_rename', 'File &Name', 'input_text', 'pars.fits.rename', 'rename each file to a new string + zero padded serial number');
            obj.menu_fits.addButton('button_finding', '&Finding chart', 'toggle', 'pars.fits.use_finding', 'save a finding chart PNG along with the FITS images');
            
            
            %%%%%%%%%%%%%%%%%%% LEFT SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            width = 0.25;
            N = 13; % number of buttons on left side
            
            pos = N;
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            % Add buttons using obj.addButton(button_name, var_name='', type='', str1='', str2='', font_size='', split=1, color_on=[], color_off=[], tooltip)
            
            num_buttons = 4;
            pos = pos-num_buttons;
            obj.panel_controls = GraphicPanel(obj.owner, [0 pos/N width num_buttons/N], 'controls', 1); % last input is for vertical (default)            
            obj.panel_controls.addButton('button_browse', 'browse', 'push', 'Browse', '', '', [], '', '', 'use the system dialog to find a folder with data files'); 
            obj.panel_controls.addButton('button_continue', 'cont', 'push', 'Continue run', '', '', [], '', '', 'continue the run from where it stopped, without resetting'); 
            obj.panel_controls.addButton('input_num_files', 'pars.num_files', 'input', 'num_files= ', '', '', [], '', '', 'how many files we want this run to go for'); 
            obj.panel_controls.addButton('button_run', 'run', 'push', 'Start new run', '', '', [], '', '', 'start a new run, resetting the current results'); 
            obj.panel_controls.number = num_buttons;
            obj.panel_controls.margin = [0.02 0.03];
            obj.panel_controls.make;
            
            %%%%%%%%%%% panel objects %%%%%%%%%%%%%%%
            
            num_buttons = 3;
            pos = pos-num_buttons;
            obj.panel_objects = GraphicPanel(obj.owner, [0 pos/N width num_buttons/N], 'objects', 1); % last input is for vertical (default)            
            obj.panel_objects.addButton('button_calibration', 'cal', 'push', 'Calibration', '', '', 0.5, '', '', 'open the GUI for the calibration object'); 
            obj.panel_objects.addButton('button_reader', 'reader', 'push', 'Reader', '', '', 0.5, '', '', 'open the GUI for the file reader'); 
            obj.panel_objects.addButton('button_photometry', 'phot', 'push', 'Photometry', '', '', 0.5, '', '', 'open the GUI for photometry object'); 
            obj.panel_objects.addButton('button_lightcurves', 'light', 'push', 'Lightcurves', '', '', 0.5, '', '', 'open the GUI for the lightcurves object'); 
            obj.panel_objects.addButton('button_pars', 'pars', 'push', 'Parameters', '', '', 0.5, '', '', 'open the GUI for parameters'); 
            obj.panel_objects.addButton('button_head', 'head', 'push', 'Header', '', '', 0.5, '', '', 'open the GUI for the header object'); 
            obj.panel_objects.number = num_buttons;
            obj.panel_objects.margin = [0.02 0.03];
            obj.panel_objects.make;
            
            %%%%%%%%%%% panel contrast %%%%%%%%%%%%%%%
            
            num_buttons = 5;
            pos = pos-num_buttons;
            obj.panel_contrast = util.plot.ContrastLimits(obj.axes_image, obj.fig.fig, [0 pos/N width 5/N], 1, [0.02 0.03]); % 4th input is for vertical (default)
            obj.panel_contrast.font_size = obj.font_size;
            obj.panel_contrast.big_font_size = obj.big_font_size;
            obj.panel_contrast.small_font_size = obj.small_font_size;
            obj.panel_contrast.edit_font_size = obj.edit_font_size;
            obj.panel_contrast.margin = [0.02 0.03];
            
            %%%%%%%%%%%%%%%%%%%%%% MIDDLE %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_info = GraphicPanel(obj.owner, [width (N-1)/N 1-width 1/N], 'info', 1); % last input is for vertical (default)    
            obj.panel_info.addButton('button_info', 'printout', 'custom'); % no need to update this automatically, set to custom
            obj.panel_info.make;
            
            obj.panel_stop = GraphicPanel(obj.owner, [width 0 1-width 1/N], '', 1); % last input is for vertical (default)    
            obj.panel_stop.addButton('button_stop', 'stop', 'push', 'Stop run'); 
            obj.panel_stop.make;
            
            obj.panel_image = uipanel('Title', '', 'Position', [width 1/N 1-width (N-2)/N]);
                        
            obj.makeAxes;
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.85 0.95 0.15 0.05], obj.owner, '', 'custom', 'new axes', '');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            obj.button_reset_axes.Tooltip = 'Create a new image axis, zoomed out and with default contrast limits'; 
            
            obj.button_display_text = GraphicButton(obj.panel_image, [0.0 0.0 0.15 0.05], obj.owner, 'pars.use_display_rect_text', 'toggle', 'text: off', 'text: on'); 
            obj.button_display_text.Tooltip = 'turn on display of cutout number with each square';
            
            obj.input_num_rect = GraphicButton(obj.panel_image, [0.15 0.0 0.2 0.05], obj.owner, 'pars.display_num_rect_stars', 'input', 'num cutouts= '); 
            obj.input_num_rect.Tooltip = 'number of squares around stars to show'; 
            
            obj.button_expt = GraphicButton(obj.panel_image, [0.0 0.95 0.15 0.05], obj.owner, 'head.EXPTIME', 'info', 'T= ', 's'); 
            
            obj.button_indices = GraphicButton(obj.panel_image, [0.6 0.0 0.4 0.05], obj.owner, 'getIndicesString', 'info', ' '); 
            
            obj.panel_legend = uipanel('Title', 'legend', 'Position', [width, 0.2, 0.11, 0.3]); 
            
            obj.panel_buffers = uipanel('Title', 'buffers', 'Position', [0.93 0.2 0.05 0.6]); 
            for ii = 1:obj.owner.pars.coadd_size
                uicontrol(obj.panel_buffers, 'Units', 'Normalized', 'Style', 'pushbutton', ...
                    'Position', [0, (ii-1)/obj.owner.pars.coadd_size, 1, 1/obj.owner.pars.coadd_size]); 
            end
            
            obj.makeLegend;
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Position', [0 0 width 1/N]);
            obj.button_close = GraphicButton(obj.panel_close, [0.02 0.03 0.96 0.95], obj.owner, '', 'custom', 'CLOSE GUI');
            obj.button_close.Callback = @obj.callback_close;
            
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
           
        function makeLegend(obj)
            
            delete(obj.panel_legend.Children);
            
            N = 6;
            color = 0.6;
            
            uicontrol(obj.panel_legend, 'Style', 'text', 'FontSize', obj.edit_font_size, ...
                'ForegroundColor', 'black', 'String', 'normal stars', ...
                'Units', 'Normalized', 'Position', [0.0 (N-1)/N 1.0 1/N], ...
                'BackgroundColor', [1,1,1].*color); 
            
            uicontrol(obj.panel_legend, 'Style', 'text', 'FontSize', obj.edit_font_size, ...
                'ForegroundColor', 'red', 'String', 'saturated', ...
                'Units', 'Normalized', 'Position', [0.0 (N-2)/N 1.0 1/N], ... 
                'BackgroundColor', [1,1,1].*color); 
            
            uicontrol(obj.panel_legend, 'Style', 'text', 'FontSize', obj.edit_font_size, ...
                'ForegroundColor', 'green', 'String', 'extended source', ...
                'Units', 'Normalized', 'Position', [0.0 (N-3)/N 1.0 1/N], ...
                'BackgroundColor', [1,1,1].*color); 
           
            uicontrol(obj.panel_legend, 'Style', 'text', 'FontSize', obj.edit_font_size, ...
                'ForegroundColor', 'yellow', 'String', 'narrow source', ...
                'Units', 'Normalized', 'Position', [0.0 (N-4)/N 1.0 1/N], ...
                'BackgroundColor', [1,1,1].*color); 
            
            uicontrol(obj.panel_legend, 'Style', 'text', 'FontSize', obj.edit_font_size, ...
                'ForegroundColor', 'white', 'String', 'no Gaia match', ...
                'Units', 'Normalized', 'Position', [0.0 (N-5)/N 1.0 1/N], ...
                'BackgroundColor', [1,1,1].*color); 
                        
            uicontrol(obj.panel_legend, 'Style', 'text', 'FontSize', obj.edit_font_size, ...
                'ForegroundColor', 'magenta', 'String', 'forced target', ...
                'Units', 'Normalized', 'Position', [0.0 (N-6)/N 1.0 1/N], ...
                'BackgroundColor', [1,1,1].*color); 
            
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
            
            obj.panel_contrast.update;
            
            obj.panel_controls.input_num_files.String = sprintf('num_files= %d', obj.owner.getNumFiles); 
            
            obj.updateBufferPanel;
            
            obj.owner.show('ax', obj.axes_image); 
            
        end
        
        function updateBufferPanel(obj)
            
            for ii = 1:obj.owner.pars.coadd_size
                
                color = 0.94*[1 1 1];
                
                if ii <= length(obj.owner.buffers)
                    if obj.owner.buffers(ii).is_loaded
                        color = 'green';
                    end

                    if obj.owner.buffers(ii).is_processed
                        color = 'magenta';
                    end
                end
                
                obj.panel_buffers.Children(obj.owner.pars.coadd_size-ii+1).BackgroundColor = color;
                
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