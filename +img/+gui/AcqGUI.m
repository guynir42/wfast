classdef AcqGUI < handle
    
    properties 
        
        owner@img.Acquisition; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 13;
        big_font_size = 16;
        edit_font_size = 12;
        small_font_size = 11;
        
        color_on = [0 0.3 1];
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_controls;
        
        panel_contrast;
    
        panel_objects;
        
        panel_save;
        
        panel_sync;
        
        panel_run;
        button_run;
        
        panel_info;
        panel_image;
        button_reset_axes;
        button_batch_num;
        button_time_left;
        button_gb_left; 
        
        button_show_what;
        button_flip;
        axes_image;
    
        panel_close;
        button_close;
        
    end
    
    properties (Hidden=true)
              
        version = 1.01;
        
    end
            
    methods % constructor
       
        function obj = AcqGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'img.Acquisition')
                
                if obj.debug_bit, fprintf('AcqGUI constructor v%4.2f\n', obj.version); end
                
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
            
            obj.buttons = {};
            
            obj.fig = util.plot.FigHandler('Acquisition');
            
            obj.fig.bottom = 5;
            obj.fig.height = 20;
            obj.fig.width = 36;
            movegui(obj.fig.fig, 'center');
            obj.fig.clear;
            
            N_left = 20;
            N_right = 20;
            W_left = 0.25;
            W_right = 0.15;
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            pos = N_left;
            N = 13;
            pos = pos - N;
            obj.panel_controls = GraphicPanel(obj.owner, [0 pos/N_left W_left N/N_left], 'controls');
            obj.panel_controls.number = N;
            obj.panel_controls.addButton('button_source_choose', 'chooseSource', 'push', 'choose src', '', '', 2/3);
            obj.panel_controls.addButton('button_source_gui', 'src', 'push', 'GUI', '', 'big', 1/3);
            
            obj.panel_controls.addButton('button_reset', 'reset', 'push', 'RESET', '', '', 1/3);
            obj.panel_controls.addButton('input_name', 'run_name', 'input_text', 'object= ', '', '', 2/3);
            
            obj.panel_controls.addButton('button_num_files', 'num_files', 'info', 'max= ', '', '', 1/3);
            obj.panel_controls.addButton('input_num_batches', 'num_batches', 'input', 'Nbatch= ', '', '', 1/3);
            obj.panel_controls.addButton('input_batch_size', 'batch_size', 'input', 'Nframe= ', '', '', 1/3);
            
            obj.panel_controls.addButton('button_stars_found', 'num_stars_found', 'info', ' ', '', '', 1/3);
            obj.panel_controls.addButton('input_num_stars', 'num_stars', 'input', 'Nstars= ', '', '', 1/3);
            obj.panel_controls.addButton('input_cut_size', 'cut_size', 'input', 'size= ', '', '', 1/3);
            
            obj.panel_controls.addButton('input_edges', 'avoid_edges', 'input', 'edges= ', '', '', 1/3);
            obj.panel_controls.addButton('input_num_bgs', 'num_backgrounds', 'input', 'Nbgs= ', '', '', 1/3);
            obj.panel_controls.addButton('input_cut_size_bg', 'cut_size_bg', 'input', 'size= ', '', '', 1/3);
            
            obj.panel_controls.addButton('input_expT', 'expT', 'input', 'T= ', 's', '', 0.5);
            obj.panel_controls.addButton('input_frame_rate', 'frame_rate', 'input', 'f= ', 'Hz', '', 0.5);
            
            obj.panel_controls.addButton('button_mextractor', 'use_mextractor', 'toggle', 'mextractor off', 'mextractor on', '', 0.5, obj.color_on);
            obj.panel_controls.addButton('button_arbitrary_pos', 'use_arbitrary_pos', 'toggle', 'find pos', 'arbitrary pos', '', 0.5, obj.color_on, 'red');
            
            obj.panel_controls.addButton('button_adjust', 'use_adjust_cutouts', 'toggle', 'no adjust', 'adjust on', '', 0.5, obj.color_on);
            obj.panel_controls.addButton('button_background', 'use_background', 'toggle', 'sub b/g off', 'sub b/g on', '', 0.5, obj.color_on);
            
            obj.panel_controls.addButton('button_simple_phot', 'use_simple_photometry', 'toggle', 'full phot', 'simple phot', '', 0.5, obj.color_on);
            obj.panel_controls.addButton('button_model_psf', 'use_model_psf', 'toggle', 'no model PSF', 'use model PSF', '', 0.5, obj.color_on);
            
            obj.panel_controls.margin = [0.02 0.01];
            obj.panel_controls.make;
            
            obj.panel_controls.button_num_files.Callback = @obj.callback_num_files;
%             obj.panel_controls.button_num_files.Enable = 'inactive';
%             obj.panel_controls.button_stars_found.Enable = 'inactive';
            
            obj.panel_controls.button_source_choose.Tooltip = 'Choose camera or reader';
            obj.panel_controls.button_source_gui.Tooltip = 'Open the GUI of the camera/reader';
            obj.panel_controls.button_reset.Tooltip = 'Starts a new run. Reset the positions and batch counter';
            obj.panel_controls.input_name.Tooltip = 'Object name (will be used also as folder name)';
            obj.panel_controls.button_num_files.Tooltip = 'Maximum number of batches from reader. For camera it is Inf';
            obj.panel_controls.input_num_batches.Tooltip = 'How many batches should this run be';
            obj.panel_controls.input_batch_size.Tooltip = ['Number of images for each batch (default=' num2str(obj.owner.default_batch_size) ')'];
            obj.panel_controls.button_stars_found.Tooltip = 'How many stars were found by "findStars" or mextractor';
            obj.panel_controls.input_num_stars.Tooltip = 'Maximum number of stars/cutouts in this run';
            obj.panel_controls.input_cut_size.Tooltip = 'Size of square cutouts (in pixels)'; 
            obj.panel_controls.input_edges.Tooltip = 'Avoids finding stars close to the edge of the frame (pixels)';
            obj.panel_controls.input_num_bgs.Tooltip = 'Number of background sampling points';
            obj.panel_controls.input_cut_size_bg.Tooltip = 'Size of square cutouts of background sample points';
            obj.panel_controls.input_expT.Tooltip = 'Exposure time for each frame (seconds)';
            obj.panel_controls.input_frame_rate.Tooltip = 'Nominal frame rate to be uphel by camera (Hz). Use NaN to make camera take images as fast as possible';
            obj.panel_controls.button_mextractor.Tooltip = 'Use mextractor and astrometry to find stars and match them to a catalog';
            obj.panel_controls.button_arbitrary_pos.Tooltip = 'Find stars in the image or choose arbitrary positions (for debugging)';
            obj.panel_controls.button_adjust.Tooltip = 'Adjust the cutout positions as stars drift in the field';
            obj.panel_controls.button_background.Tooltip = 'Use background subtraction on all analysis data';
            obj.panel_controls.button_simple_phot.Tooltip = 'Use simple summing of cutouts <or> full photometry object';
            obj.panel_controls.button_model_psf.Tooltip = 'Use the stack of PSFs to model the average width, etc...';
            
            %%%%%%%%%%% panel contrast %%%%%%%%%%%%%%%
            
            N = 5;
            pos = pos - N;
            obj.panel_contrast = util.plot.ContrastLimits(obj.axes_image, obj.fig.fig, [0 pos/N_left W_left N/N_left], 1, [0.02 0.01]);
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Position', [0 0 W_left 2/N_left]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1]+[0.05 0.1 -0.1 -0.2], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            obj.button_close.Tooltip = 'Close the GUI';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RIGHT SIDE %%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%% panel objects %%%%%%%%%%%%%%%%
            
            pos = N_right;
            N = 14;
            pos = pos - N;
            obj.panel_objects = GraphicPanel(obj.owner, [1-W_right pos/N_right W_right N/N_right], 'objects');
            obj.panel_objects.number = N;
            obj.panel_objects.addButton('button_pars', 'pars', 'push', 'Header GUI');
            obj.panel_objects.addButton('button_camera', 'cam', 'push', 'Camera GUI');
            obj.panel_objects.addButton('button_reader', 'reader', 'push', 'Reader GUI');
            obj.panel_objects.addButton('button_buffers', 'buffers', 'push', 'Buffers GUI');
            obj.panel_objects.addButton('button_calibration', 'cal', 'push', 'Calibration GUI');
            obj.panel_objects.addButton('button_background', 'back', 'push', 'Background GUI');
            obj.panel_objects.addButton('button_clipper', 'clip', 'push', 'Clipper GUI');
            obj.panel_objects.addButton('button_clipper_bg', 'clip_bg', 'push', 'Clipper BG GUI');
            obj.panel_objects.addButton('button_photometry', 'phot', 'push', 'Photometry GUI');
            obj.panel_objects.addButton('button_phot_stack', 'phot_stack', 'push', 'stack Phot GUI');
            obj.panel_objects.addButton('button_lightcurves', 'lightcurves', 'push', 'Lightcurves GUI');
            obj.panel_objects.addButton('button_deflator', 'deflator', 'push', 'Deflator GUI');
            obj.panel_objects.addButton('button_focuser', '', 'custom', 'Focuser: '); 
            obj.panel_objects.margin = [0.02 0.01];
            obj.panel_objects.make;
            obj.panel_objects.button_focuser.Callback = @obj.callback_focuser;
            
            %%%%%%%%%%% panel save %%%%%%%%%%%%%%%%%%%
            
            N = 3;
            pos = pos - N;
            obj.panel_save = GraphicPanel(obj.owner, [1-W_right pos/N_right W_right N/N_right], 'save');
            obj.panel_save.number = N;
            obj.panel_save.addButton('button_save', 'use_save', 'toggle', 'save off', 'save on', '', 1, obj.color_on, 'red');
            obj.panel_save.addButton('button_trig', 'use_triggered_save', 'toggle', 'trig save off', 'trig save on', '', 1, obj.color_on);
            obj.panel_save.margin = [0.02 0.01];
            obj.panel_save.make;
            
            obj.panel_save.button_save.Tooltip = 'Save cutouts and stack images of each batch';
            obj.panel_save.button_trig.Tooltip = 'Save full frame images when triggered on events';
            
            %%%%%%%%%%% panel sync %%%%%%%%%%%%%%%%%%%
            
            N = 3;
            pos = pos - N;
            obj.panel_sync = GraphicPanel(obj.owner, [1-W_right pos/N_right W_right N/N_right], 'sync');
            obj.panel_sync.number = N;
            obj.panel_sync.addButton('button_sync', 'use_sync', 'toggle', 'sync off', 'sync on', '', 1, obj.color_on);
            obj.panel_sync.addButton('button_sync_stop', 'use_sync_stop', 'toggle', 'stop ignored', 'stop enabled', '', 1, obj.color_on);
            obj.panel_sync.addButton('button_autoguide', 'use_autoguide', 'toggle', 'guiding off', 'guiding on', '', 1, obj.color_on);
            obj.panel_sync.make;
            
            obj.panel_sync.button_ignore_manager.Tooltip = 'Allow Manager to update camera with weather and coordinates';
            obj.panel_sync.button_ignore_manager_stop.Tooltip = 'Allow Manager to send stop command when tracking is off or when shutting down observatory';
            
            %%%%%%%%%%% panel run %%%%%%%%%%%%%%%%%%
            
            % obj.panel_run = uipanel('Title','', 'Position', [W_left, 0, 1-W_left-W_right, 0.1]);
            obj.panel_run = GraphicPanel(obj.owner, [W_left 0 (1-W_left-W_right) 0.1]);
            
            obj.panel_run.addButton('button_preview', 'runPreview', 'push', 'PREVIEW', '', '', 0.15);
            obj.panel_run.addButton('button_run', '', 'custom', 'RUN', '', '', 0.7);
            obj.panel_run.addButton('button_focus', 'runFocus', 'push', 'FOCUS', '', '', 0.15);
            obj.panel_run.margin = [0.005 0.1];
            obj.panel_run.make;
            obj.panel_run.button_run.Callback = @obj.callback_run;
            obj.panel_run.button_preview.Tooltip = 'Run a single batch just to see the field';
            obj.panel_run.button_run.Tooltip = 'Start a run / continue a run / stop the run';
            obj.panel_run.button_focus.Tooltip = 'Start an autofocus run';
            
            %%%%%%%%%%% panel info %%%%%%%%%%%%%%%%%%
            
            obj.panel_info = GraphicPanel(obj.owner, [W_left 0.9 1-W_left-W_right 0.1], 'info', 0);
            obj.panel_info.addButton('button_frame_rate', 'frame_rate_average', 'info', 'f= ', ' Hz'); 
            obj.panel_info.addButton('button_width', 'average_width', 'info', 'width= ', ' px', 'small', 0.5);
            obj.panel_info.addButton('button_seeing', 'average_width', 'custom', 'seeing= ', ' "', 'small', 0.5);
            obj.panel_info.addButton('button_min_axis', 'minor_axis', 'info', 'b= ', ' px', 'small', 0.5);
            obj.panel_info.addButton('button_maj_axis', 'major_axis', 'info', 'a= ', ' px', 'small', 0.5);
            obj.panel_info.addButton('button_offset_y', 'average_offsets', 'custom', 'dy= ', ' px', 'small', 0.5);
            obj.panel_info.addButton('button_offset_x', 'average_offsets', 'custom', 'dx= ', ' px', 'small', 0.5);
            obj.panel_info.addButton('button_flux', 'average_flux', 'info', 'flux= ', '', 'small', 0.5);
            obj.panel_info.addButton('button_bg', 'average_background', 'info', 'b/g= ', '', 'small', 0.5);
            obj.panel_info.addButton('button_temperature', 'sensor_temperature', 'info', 's.temp= '); 
            obj.panel_info.margin = [0.005 0.03];
            obj.panel_info.make;
            
            obj.panel_info.button_frame_rate.Tooltip = 'Measured frame rate (averaged over a few batches)';
            obj.panel_info.button_width.Tooltip = '2nd moment of stars, averaged over all cutouts (pixels)';
            obj.panel_info.button_seeing.Tooltip = 'Full width at half maximum (arcsec)';
            obj.panel_info.button_min_axis.Tooltip = 'Minor axis of the stacked PSF';
            obj.panel_info.button_maj_axis.Tooltip = 'Major axis of the stacked PSF';
            obj.panel_info.button_offset_y.Tooltip = 'Average offset/drift/1st moment in the y direction (pixels)';
            obj.panel_info.button_offset_x.Tooltip = 'Average offset/drift/1st moment in the x direction (pixels)';
            obj.panel_info.button_flux.Tooltip = 'Average flux of all stars';
            obj.panel_info.button_bg.Tooltip = 'Average background of all stars';
            obj.panel_info.button_temperature.Tooltip = 'Sensor temperature (C)';
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [W_left 0.1 1-W_left-W_right 0.8]);
                        
            obj.makeAxes;
            
            obj.button_batch_num = GraphicButton(obj.panel_image, [0.00 0.95 0.1 0.05], obj.owner, 'batch_counter', 'info', 'N= ', '', 'small');
            obj.button_time_left = GraphicButton(obj.panel_image, [0.25 0.95 0.35 0.05], obj.owner, 'getTimeLeftHMS', 'info', ' ', '', 'small');
            obj.button_gb_left = GraphicButton(obj.panel_image, [0.60 0.95 0.15 0.05], obj.owner, 'getGbLeft', 'info', ' ', 'Gb', 'small');
            
            obj.button_show_what = GraphicButton(obj.panel_image, [0 0.00 0.15 0.05], obj.owner, 'show_what', 'picker', 'full', '', 'small');
            obj.button_show_what.Callback = @obj.callback_show_what;
            obj.button_show_what.String = obj.owner.show_what_list;
            
            obj.button_flip = GraphicButton(obj.panel_image, [0.9 0.00 0.1 0.05], obj.owner, 'use_flip', 'toggle', 'flip', '', 'small');
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.9 0.95 0.1 0.05], obj.owner, '', 'custom','new axes');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            
            obj.button_batch_num.Tooltip = 'How many batches already taken in this run';
            obj.button_time_left.Tooltip = 'Estimate for remaining run time';
            obj.button_gb_left.Tooltip = 'Estimate for hard drive spaced needed for the remaining batches in this run';
            obj.button_show_what.Tooltip = 'Display images, stack, or stack processed (dark, flat, background removed)';
            obj.button_flip.Tooltip = 'Flip the image 180 degrees (for viewing beyond the meridian)';
            obj.button_reset_axes.Tooltip = 'Generate a new image axes with contrast and zoom initialized';
            
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
            
            obj.panel_contrast.font_size = obj.font_size;
            obj.panel_contrast.edit_font_size = obj.edit_font_size;
            obj.panel_contrast.update;
            
            if isa(obj.owner.src, 'file.Reader')
                obj.panel_controls.button_source_choose.String = 'src: Reader';
            elseif isa(obj.owner.src, 'img.Simulator')
                obj.panel_controls.button_source_choose.String = 'src: Simulator';
            elseif isa(obj.owner.src, 'obs.cam.CameraControl') || isa(obj.owner.src, 'obs.cam.Andor')
                obj.panel_controls.button_source_choose.String = 'src: Camera';
            end
            
            obj.panel_info.button_seeing.String = sprintf('seeing= %4.2f"', obj.owner.average_width.*obj.owner.pars.plate_scale.*2.355);
            
%             obj.panel_info.button_offsets.String = sprintf('dx,dy= %4.2f,%4.2f', obj.owner.average_offsets(2), obj.owner.average_offsets(1));
            if length(obj.owner.average_offsets)==2
                obj.panel_info.button_offset_x.String = sprintf('dx= %4.2f', obj.owner.average_offsets(2));
                obj.panel_info.button_offset_y.String = sprintf('dy= %4.2f', obj.owner.average_offsets(1));
            else
                obj.panel_info.button_offset_x.String = 'dx= ';
                obj.panel_info.button_offset_y.String = 'dy= ';
            end
            
            if obj.owner.brake_bit && (obj.owner.is_running || obj.owner.is_running_single) % started a run but not yet took off the brake bit
                obj.panel_run.button_run.Enable = 'off';
%                 obj.panel_run.button_run.String = 'STARTING UP';
            else

                if obj.owner.brake_bit
                    obj.panel_run.button_run.Enable = 'on';
                    if obj.owner.batch_counter>0
                        obj.panel_run.button_run.String = 'CONTINUE THIS RUN';
                    else
                        obj.panel_run.button_run.String = 'START NEW RUN';
                    end
    %                 
                else % still running
                    obj.panel_run.button_run.Enable = 'on';
                    obj.panel_run.button_run.String = 'STOP';
                end

            end
            
            if obj.owner.is_running || obj.owner.is_running_single
                
                obj.panel_controls.button_source_choose.Enable = 'off';
                obj.panel_controls.button_reset.Enable = 'off';
                obj.panel_controls.input_name.Enable = 'off';
                obj.panel_controls.input_num_batches.Enable = 'off';
                obj.panel_controls.input_batch_size.Enable = 'off';
                obj.panel_controls.input_num_stars.Enable = 'off';
                obj.panel_controls.input_cut_size.Enable = 'off';
                obj.panel_controls.input_edges.Enable = 'off';
                obj.panel_controls.input_num_bgs.Enable = 'off';
                obj.panel_controls.input_cut_size_bg.Enable = 'off';
                obj.panel_controls.input_expT.Enable = 'off';
                obj.panel_controls.input_frame_rate.Enable = 'off';
                obj.panel_controls.button_mextractor.Enable = 'off';
                obj.panel_controls.button_arbitrary_pos.Enable = 'off';
                obj.panel_controls.button_adjust.Enable = 'off';
                obj.panel_controls.button_background.Enable = 'off';
                obj.panel_controls.button_simple_phot.Enable = 'off';
                obj.panel_controls.button_model_psf.Enable = 'off';
                
                obj.panel_run.button_preview.Enable = 'off';
                obj.panel_run.button_focus.Enable = 'off';
                
            else
                
                obj.panel_controls.button_source_choose.Enable = 'on';
                obj.panel_controls.button_reset.Enable = 'on';
                obj.panel_controls.input_name.Enable = 'on';
                obj.panel_controls.input_num_batches.Enable = 'on';
                obj.panel_controls.input_batch_size.Enable = 'on';
                obj.panel_controls.input_num_stars.Enable = 'on';
                obj.panel_controls.input_cut_size.Enable = 'on';
                obj.panel_controls.input_edges.Enable = 'on';
                obj.panel_controls.input_num_bgs.Enable = 'on';
                obj.panel_controls.input_cut_size_bg.Enable = 'on';
                obj.panel_controls.input_expT.Enable = 'on';
                obj.panel_controls.input_frame_rate.Enable = 'on';
                obj.panel_controls.button_mextractor.Enable = 'on';
                obj.panel_controls.button_arbitrary_pos.Enable = 'on';
                obj.panel_controls.button_adjust.Enable = 'on';
                obj.panel_controls.button_background.Enable = 'on';
                obj.panel_controls.button_simple_phot.Enable = 'on';
                obj.panel_controls.button_model_psf.Enable = 'on';
                
                obj.panel_run.button_preview.Enable = 'on';
                obj.panel_run.button_focus.Enable = 'on';
                
            end
            
            for ii = 1:length(obj.button_show_what.String)
                if util.text.cs(obj.owner.show_what, obj.button_show_what.String{ii})
                    obj.button_show_what.Value = ii;
                    break;
                end
            end
            
            if ~isempty(obj.owner.cam) && ~isempty(obj.owner.cam.focuser)
                obj.panel_objects.button_focuser.String = sprintf('Focuser: %6.4f', obj.owner.cam.focuser.pos);
            end
            
            drawnow;
            
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
                obj.update;
            end
            
        end
        
        function callback_num_files(obj, ~, ~)
            
            
            if obj.debug_bit, disp('callback: num_files'); end
            
            if ~isinf(obj.owner.num_files)
                obj.owner.num_batches = obj.owner.num_files;
            end
            
            obj.update;
        
        end
        
        function callback_show_what(obj, hndl, ~)
            
            if obj.debug_bit, disp('callback: show_what'); end
            
%             obj.owner.show_what_index = hndl.Value;
%             obj.owner.show_what = obj.owner.show_what_list{obj.owner.show_what_index};

            obj.owner.show_what = hndl.String{hndl.Value};

            obj.owner.show;
            obj.panel_contrast.autodyn;
            
            obj.update;
            
        end
        
        function callback_focuser(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: focuser'); end
        
            obj.owner.cam.focuser.makeGUI;
            
            obj.update;
            
        end
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end