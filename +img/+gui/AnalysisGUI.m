classdef AnalysisGUI < handle
    
    properties 
        
        owner@img.Analysis; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 12;
        big_font_size = 16;
        edit_font_size = 11;
        small_font_size = 10;
        
        color_on = [0 0.3 1];
        
        dialog_fig;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_controls;
        
        panel_contrast;
    
        panel_close;
        button_close;
    
        panel_objects;
        panel_fits;
            
        panel_progress;
        
        panel_info;
        
        panel_image;
        button_reset_axes;
        button_batch_counter;
        button_flip; 
        input_num_rect;
        axes_image;
    
        panel_run;
        
    end
    
    properties (Hidden=true)
              
        version = 1.01;
        
    end
            
    methods % constructor
       
        function obj = AnalysisGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'img.Analysis')
                
                if obj.debug_bit, fprintf('AnalysisGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input an img.Analysis to constructor of AnalysisGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            
            obj.fig = util.plot.FigHandler('analysis');
            obj.fig.clear;
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 25;
            movegui(obj.fig.fig, 'center');
            
            %%%%%%%%%%%%%%%%%%% LEFT SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_left = 20; pos = N_left;
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            % Add buttons using obj.addButton(button_name, var_name='', type='', str1='', str2='', font_size='', split=1, color_on=[], color_off=[])
            
            N = 13; pos = pos - N;
            
            obj.panel_controls = GraphicPanel(obj.owner, [0 pos/N_left 0.25 N/N_left], 'controls', 1); % last input is for vertical (default)
            obj.panel_controls.number = N;
            obj.panel_controls.addButton('button_browse', 'chooseDir', 'push', 'browse', '', '', 0.5, '', '', 'Choose a folder to analize'); 
            obj.panel_controls.addButton('button_async_run', 'async_run', 'custom', 'run async', '', '', 0.5, '', '', 'Run the current object in a parallel worker');
            obj.panel_controls.addButton('button_reset', 'reset', 'push', 'RESET', '', '', 0.5, '', '', 'Start a new run by reseting all events and lightcurves');
            obj.panel_controls.addButton('input_num_batches', 'num_batches', 'input', 'Nbatch= ', '', 'small', 0.5, '', '', 'Maximum batches, limited by user input or by number of files in reader'); 
            obj.panel_controls.addButton('button_bg_stack', 'use_background_stack', 'toggle', 'b/g stack off', 'b/g stack on', 'small', 0.5, obj.color_on, '', 'Subtract background from stack images');
            obj.panel_controls.addButton('button_bg_cutouts', 'use_background_cutouts', 'toggle', 'b/g cutouts off', 'b/g cutouts on', 'small', 0.5, obj.color_on, '', 'Subtract background from cutouts');
            obj.panel_controls.addButton('button_astrometry', 'use_astrometry', 'auto', 'astrometry', '', 'small', 0.5, obj.color_on, '', 'Use mextractor, astrometry and match resulting stars to GAIA');             obj.panel_controls.addButton('button_cutouts', 'use_cutouts', 'toggle', 'no cutouts', 'with cutouts', 'small', 0.5, obj.color_on, '', 'Run cutouts analysis'); 
            obj.panel_controls.addButton('button_photometry', 'use_photometry', 'toggle', 'no photometry', 'with photometry', 'small', 0.5, obj.color_on, '', 'Run photometry on cutotus'); 
            obj.panel_controls.addButton('button_psf_model', 'use_psf_model', 'toggle', 'no PSF modeling', 'with PSF modeling', 'small', 0.5, obj.color_on, '', 'Calculate the model PSF from all cutouts.'); 
            obj.panel_controls.addButton('button_events', 'use_event_finding', 'toggle', 'no event search', 'use event search', 'small', 0.5, obj.color_on, '', 'Run event search using match filtering of lightcurves'); 
            obj.panel_controls.addButton('button_auto_load_cal', 'use_auto_load_cal', 'toggle', 'no load cal', 'auto load cal', 'small', 0.5, obj.color_on, '', 'Always test for most relevant calibration file when starting a new run'); 
            obj.panel_controls.addButton('button_check_flux', 'use_check_flux', 'toggle', 'no check flux', 'use check flux', 'small', 0.5, obj.color_on, '', 'Check if flux is lost a few times in a row, quit the run. '); 
            obj.panel_controls.margin = [0.03 0.01];
            obj.panel_controls.make;
            
            obj.panel_controls.button_async_run.Callback = @obj.callback_async_run;
            
            %%%%%%%%%%% panel contrast %%%%%%%%%%%%%%%
            
            N = 5; pos = pos - N;
            
            obj.panel_contrast = util.plot.ContrastLimits(obj.axes_image, obj.fig.fig, [0 pos/N_left 0.25 5/N_left], 1, [0.01 0.01]); % 4th input is for vertical (default), 5th is margin
            obj.panel_contrast.font_size = obj.font_size;
            obj.panel_contrast.big_font_size = obj.big_font_size;
            obj.panel_contrast.small_font_size = obj.small_font_size;
            obj.panel_contrast.edit_font_size = obj.edit_font_size;
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            N = 2; pos = pos - N;
            
            obj.panel_close = uipanel('Position', [0 pos 0.25 N/N_left]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1]+[0.1 0.1 -0.2 -0.2], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            %%%%%%%%%%%%%%%%%%% RIGHT SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_right = 20; pos = N_right;
            
            %%%%%%%%%%% panel objects %%%%%%%%%%%%%%%%
            
            N = 16; pos = pos - N;
            
            obj.panel_objects = GraphicPanel(obj.owner, [0.8 pos/N_right 0.2 N/N_right], 'objects'); 
            obj.panel_objects.number = N;
            obj.panel_objects.addButton('button_parameters', 'pars', 'push', 'Parameters');
            obj.panel_objects.addButton('button_reader', 'reader', 'push', 'Reader');
            obj.panel_objects.addButton('button_calibration', 'cal', 'push', 'Calibration');
            obj.panel_objects.addButton('button_clipper', 'clip', 'push', 'Clipper');
            obj.panel_objects.addButton('button_background', 'back', 'push', 'Background');
            obj.panel_objects.addButton('button_photometry', 'phot', 'push', 'Photometry');
%             obj.panel_objects.addButton('button_light_basic', 'light_basic', 'push', 'basic', '', '', 1/3);
            obj.panel_objects.addButton('button_lightcurves', 'lightcurves', 'push', 'lightcurves', '', '', 1);
%             obj.panel_objects.addButton('button_light_gauss', 'light_gauss', 'push', 'gauss', '', '', 1/3);
            obj.panel_objects.addButton('button_finder', 'finder', 'push', 'Finder');
            obj.panel_objects.margin = [0.1 0.005];
            obj.panel_objects.make;
            
            %%%%%%%%%%% panel fits %%%%%%%%%%%%%%%%
            
            N = 4; pos = pos - N;
            
            obj.panel_fits = GraphicPanel(obj.owner, [0.8 pos/N_right 0.2 N/N_right], 'fits'); 
            obj.panel_fits.number = N;
            obj.panel_fits.addButton('button_fits_save', 'use_fits_save', 'toggle', 'fits save off', 'fits save on', '', [], obj.color_on);
            obj.panel_fits.addButton('button_fits_flip', 'use_fits_flip', 'toggle', 'fits flip off', 'fits flip on', '', [], obj.color_on);
            obj.panel_fits.addButton('button_fits_roi', 'use_fits_roi', 'toggle', 'fits ROI off', 'fits ROI on', '', [], obj.color_on);
            obj.panel_fits.addButton('input_fits_roi', 'fits_roi', 'input', 'ROI= ');
            
            obj.panel_fits.margin = [0.1 0.005];
            obj.panel_fits.make;
            
            
            %%%%%%%%%%%%%%%%%%%%%% MIDDLE %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_middle = 20; pos = N_middle;
            
            %%%%%%%%%%% panel progres %%%%%%%%%%%%%%%%
            
            N = 1; pos = pos - N;
            obj.panel_progress = GraphicPanel(obj.owner, [0.25 pos/N_middle 0.55 N/N_middle]);
            obj.panel_progress.addButton('button_progress', '', 'custom', ' ', '', 'small');
            obj.panel_progress.margin = [0.0 0.0];
            obj.panel_progress.make;
            
            %%%%%%%%%%% panel info %%%%%%%%%%%%%%%%%%%
            
            N = 3; pos = pos - N;
            obj.panel_info = GraphicPanel(obj.owner, [0.25 pos/N_middle 0.55 N/N_middle], '', 1); 
            obj.panel_info.addButton('button_directory', 'directory', 'info', ' ', '', 'small');
            obj.panel_info.addButton('button_filename', 'filename', 'info', ' ', '', 'small');
%             obj.panel_info.addButton('button_fwhm', 'FWHM', 'info', 'FWHM= ', 'pix', '', 0.2); 
            obj.panel_info.addButton('button_seeing', 'seeing', 'info', 'seeing= ', '"', 'edit', 0.2);
            obj.panel_info.addButton('button_RA', 'pars.RA', 'info', 'RA= ', '', 'edit', 0.2); 
            obj.panel_info.addButton('button_DE', 'pars.DEC', 'info', 'DE= ', '', 'edit', 0.2); 
            obj.panel_info.addButton('button_exptime', 'pars.EXPTIME', 'info', 'EXPTIME= ', 's', 'edit', 0.2); 
            obj.panel_info.margin = 0.00;
            obj.panel_info.make;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            N = 14; pos = pos - N;
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.25 pos/N_middle 0.55 N/N_middle]);
                        
            obj.makeAxes;
            
            obj.button_batch_counter = GraphicButton(obj.panel_image, [0.0 0.0 0.15 0.05], obj.owner, 'batch_counter', 'info', 'N= ');
            obj.button_batch_counter.Tooltip = 'How many batches already finished';
            
            obj.input_num_rect = GraphicButton(obj.panel_image, [0.0 0.95 0.15 0.05], obj.owner, 'display_num_rect_stars', 'input', ' ', ' rect');
            obj.input_num_rect.Tooltip = 'How many rectangles (maximum) to show on screen';
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.85 0.95 0.15 0.05], obj.owner, '', 'custom', 'new axes', '');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            obj.button_reset_axes.Tooltip = 'Create a new image axis, zoomed out and with default contrast limits'; 
            
            obj.button_flip = GraphicButton(obj.panel_image, [0.85 0.00 0.15 0.05], obj.owner, 'use_display_flip', 'toggle', 'no flip', 'flip on');
            
            %%%%%%%%%%% panel run/stop %%%%%%%%%%%%%%%
            
            N = 2; pos = pos - N;
            
            obj.panel_run = GraphicPanel(obj.owner, [0.25 pos/N_middle 0.55 N/N_middle]);
            
            obj.panel_run.addButton('button_run', '', 'custom', 'RUN', '', 'big');
            obj.panel_run.margin = [0.02 0.1];
            obj.panel_run.make;
            obj.panel_run.button_run.Callback = @obj.callback_run;
            
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
                
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            obj.panel_progress.button_progress.String = obj.owner.prog.show;
            
            if obj.owner.brake_bit
                if obj.owner.batch_counter==0
                    obj.panel_run.button_run.String = 'START NEW RUN';
                else
                    obj.panel_run.button_run.String = 'CONTINUE';
                end
            else
                obj.panel_run.button_run.String = 'STOP';
            end
            
            obj.panel_contrast.update;
            
            if ~isempty(obj.owner.stack_proc)
                obj.owner.show('axes', obj.axes_image);
            end
            
            if obj.owner.pars.EXPTIME>=0.1
                obj.panel_info.button_exptime.BackgroundColor = 'red';
            else
                obj.panel_info.button_exptime.BackgroundColor = util.plot.GraphicButton.defaultColor;
            end
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
    
        function callback_run(obj, ~, ~)
            
            if obj.owner.brake_bit
                if obj.debug_bit, disp('callback: run'); end            
                obj.owner.run;
            else
                if obj.debug_bit, disp('callback: stop'); end            
                obj.owner.brake_bit = 1;
            end
            
            obj.update;
            
        end
        
        function callback_async_run(obj, ~, ~)
            
            delete(obj.dialog_fig);
            
            margin = 0.05;
            pos_x = margin;
            pos_y = 0.25;
            width = 0.28;
            height = 0.15;
            font_size = 14;
            font_size_small = 12;
            
            idx = obj.owner.findWorker;
            
            log_time = datetime('now', 'TimeZone', 'UTC');
            log_dir = fullfile(obj.owner.reader.dir.pwd, ['analysis_' char(log_time, 'yyyy-MM-dd')]);
            
            obj.dialog_fig = figure('units','pixels', 'position',[150 500 500 300], 'menubar','none', 'numbertitle','off', 'name', 'verify parallel run', 'resize','on');
            
            if isempty(idx)
                str = 'Cannot find a free worker to run analysis!';
                ok = 0;
%             elseif exist(log_dir, 'dir')
%                 str = 'Analysis folder already exists!'; 
%                 ok = 0;
            else
                str = sprintf('Starting a new run in worker %d!', idx);
                ok = 1;
            end
            
            button_text = uicontrol('Style', 'text', 'String', str, ... 
                'Units', 'Normalized', 'Position', [0 0.9 1 0.1], 'Parent', obj.dialog_fig,...
                'HorizontalAlignment', 'center', 'FontSize', 14);
            
            button_folder = uicontrol('Style', 'text', 'String', sprintf('Folder: %s', obj.owner.reader.dir.two_tail), ... 
                'Units', 'Normalized', 'Position', [0 0.8 1 0.1], 'Parent', obj.dialog_fig,...
                'HorizontalAlignment', 'center', 'FontSize', font_size_small);
            
            if exist(log_dir, 'dir')
                str = 'Analysis folder already exists!'; 
            % ... any other warnings? 
            else
                str =  '';
            end
            
            button_warning = uicontrol('Style', 'text', 'String', str, ... 
                'Units', 'Normalized', 'Position', [0 0.7 1 0.1], 'Parent', obj.dialog_fig,...
                'HorizontalAlignment', 'center', 'FontSize', font_size_small, 'ForegroundColor', 'red');
            
            button_cancel = uicontrol('Parent', obj.dialog_fig, 'Style', 'pushbutton', 'String', 'Cancel', ...
                'Units', 'Normalized', 'Position', [margin, margin, width, height], 'FontSize', font_size, ...
                'Callback', @(~, ~, ~) delete(obj.dialog_fig)); 
            
            uicontrol(button_cancel);
            
            if ok
            
                % first column of buttons
                
                pos_y = button_warning.Position(2) - height; 

                button_reset = uicontrol('Parent', obj.dialog_fig, 'Style', 'checkbox', 'String', 'reset', 'Value', 1, ...
                    'Units', 'Normalized', 'Position', [pos_x pos_y width height], 'FontSize', font_size, ...
                    'Tooltip', 'Use this to make sure Analysis object is "reset" before starting the async-run');

                pos_y = pos_y - height;

                button_logging = uicontrol('Parent', obj.dialog_fig, 'Style', 'checkbox', 'String', 'logging', 'Value', 1, ...
                    'Units', 'Normalized', 'Position', [pos_x pos_y width height], 'FontSize', font_size,...
                    'Tooltip', 'Use this to produce text file outputs during the async-run');

                pos_y = pos_y - height;

                button_save = uicontrol('Parent', obj.dialog_fig, 'Style', 'checkbox', 'String', 'save', 'Value', 1, ...
                    'Units', 'Normalized', 'Position', [pos_x pos_y width height], 'FontSize', font_size, ...
                    'Tooltip', 'Use this to save MAT files of the event Finder and Lightcurve objects');

                pos_y = pos_y - height;

                button_run = uicontrol('Parent', obj.dialog_fig, 'Style', 'pushbutton', 'String', 'Run', ...
                    'Units', 'Normalized', 'Position', [1-margin-width, margin, width, height], 'FontSize', font_size, ...
                    'Callback', @func); 

                % second column of buttons
                pos_x = pos_x + margin + width;
                pos_y = button_warning.Position(2) - height;
                
                button_overwrite = uicontrol('Parent', obj.dialog_fig, 'Style', 'checkbox', 'String', 'overwrite', 'Value', 0, ...
                    'Units', 'Normalized', 'Position', [pos_x pos_y width height], 'FontSize', font_size, ...
                    'Tooltip', 'Use this to overwrite existing analysis folder (careful: another analysis may be ongoing!)');

                
                uicontrol(button_run); % bring the selection to "run" button

            end
            
            function func(~,~)
                
                reset = button_reset.Value;
                logging = button_logging.Value;
                save = button_save.Value;
                over = button_overwrite.Value;
                
                delete(obj.dialog_fig);
                
                if over==0 && exist(log_dir, 'dir')
                    
                    fprintf('Folder "%s " already exists!\n', log_dir);
                    
                else
                
                    fprintf('Running analysis on separate worker. reset: %d, logging: %d, save: %d, overwrite: %d\n', reset, logging, save, over);
                
                    obj.owner.async_run('reset', reset, 'logging', logging, 'save', save, 'overwrite', over);
                    
                end
                
            end
            
        end
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end