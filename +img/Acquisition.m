classdef Acquisition < file.AstroData

    properties(Transient=true)
        
        gui;
        
        % utilities
        prog@util.sys.ProgressBar;
        audio@util.sys.AudioControl;
        runtime_buffer@util.vec.CircularBuffer;
        
    end
    
    properties % objects
        
        % general objects
        head@head.Header;
        cat@head.Catalog;
        
        % input objects
        cam; % can be obs.cam.CameraControl or obs.cam.Andor
        reader@file.Reader;
%         sim; % later add the class for the simulator        
        src; % can be camera, reader or simulator
        
        % image processing
        cal@img.Calibration;
        back@img.Background;
        clip@img.Clipper;
        clip_bg@img.Clipper;
        phot_stack@img.Photometry;
        phot@img.Photometry;
        flux_buf@util.vec.CircularBuffer;
        lightcurves@img.Lightcurves;
        light_stack@img.Lightcurves;
        
        model_psf@img.ModelPSF;
        
        % output to file
        buf@file.BufferWheel;
        
        deflator@file.Deflator;
        
        sync@obs.comm.PcSync;
        
        timer;
        slow_timer;
        backup_timer;
        
        log@util.sys.Logger;
        
        obs_log; % struct with observation time and number of files and other metadata for each target
        
        forced_targets = []; % structs with target parameters (at least RA/Dec) to add to the positions/cutouts/catalog
        
        micro_flares = []; % struct with details on flares captured by cosmic ray detector
        
    end
    
    properties % inputs/outputs
        
        num_stars_found;
        cutouts_proc;
        cutouts_bg_proc;
        
        stack_cutouts; 
        stack_cutouts_bg;
        stack_proc;
        
    end
    
    properties % switches/controls
        
        run_name_append = ''; % append this to the run name/folder (e.g., "slow" or "roi")
        
        total_runtime; % how much run time we want for the next run 
        runtime_units = 'minutes';
        
        use_focus_on_start = 1; % when true, will do a focus run every time the command to start a new run is given from PcSync
        
        use_background = 0; % not sure if subtracting the background gives us anything?
        
        % these swithces determine how stars are picked when run begins
        detect_thresh = 15; % minimal S/N of the stack stars for selecting cutouts
        use_remove_bad_pixels = true;
        use_remove_saturated = true; % remove all stars with any pixels above saturation value
        
        use_astrometry = true; % calculate the star positions matched to GAIA DR2 and save catalog
        
        use_cutouts = true;
        use_adjust_cutouts = 1; % use adjustments in software (not by moving the mount)
        use_lock_adjust = 0; % make all cutouts move together based on the average drift (this is not yet implemented!)
        % add switch to allow unlocking only some of the cutouts (e.g.,
        % those without a match to GAIA
        
        use_simple_photometry = 1; % use only sums on the cutouts instead of Photometry object for full cutouts (for now we keep this on, to maintain 25 Hz)
        use_store_photometry = 0; % store the photometric products in the Lightcurve object for entire run
        use_save_photometry = 1; % save the flux and other products in the HDF5 files along with the images
        use_save_stack_lcs = 1; % save the Lightcurves object from phot_stack to disk at end of run
        
        use_dynamic_cutouts = 1; % use find_cosmic_rays to detect bleeps in the full data cube and assign cutouts to them
        num_dynamic_cutouts = 5; % how many additional cutouts we want
        
        use_model_psf = 0;
        
        use_cam_focusing = 1;
        num_focus_iterations = 6;
        
        use_check_positions = 1;
        
        use_save = false; % must change this when we are ready to really start
%         use_triggered_save = false;
        
        use_sync = 1; % if false, do not send or receive messages from PcSync object
        use_ignore_manager = 0; % if true, will not use any data from Manager (via sync object)
        use_sync_stop = 0; % if false, will not respect stop commands from Manager (via sync object)
        use_autoguide = 1; % if true, send back adjustments on drifts to telescope
        use_ignore_sync_object_name = 0; % if true, will not update the head.OBJECT field so it can be changed manually
        use_autodeflate = 1;
        
        % display parameters
        use_show = true;
        
        show_what = 'images'; % can choose "images" or "stack"
        display_num_rect_stars = 30;
        display_num_rect_bg = 30;
        
        use_flip = 0; % flip view by 180 degrees (for meridien flip)
        use_show_gray = 0; % display in gray instead of default colormap
        
        use_audio = 1;
        use_progress = 1;
        use_print_timing = 0; % print a line each batch with the timing data for each part of the processing chain
        
        debug_bit = 1;
        log_level = 1;
        
    end
    
    properties(Dependent=true)
        
        % get these from camera/reader
        num_batches;
        batch_size;
        expT;
        frame_rate;
        use_roi;
        roi_size;
        roi_center;
        
        % get these from Clipper
        num_stars;
        cut_size;
        avoid_edges;
        
        frame_rate_average; % calculated from this object's timing data
        
    end
        
    properties(Dependent=true, Hidden=true)
        
        run_name; % the parameter "OBJECT" is used here to give a name to the whole run
        buf_full; % camera's buffers are used for full-frame dump on triggers
        
        % these are read only camera info
        sensor_temperature;
        
        % get these from background Clipper
        num_backgrounds;
        cut_size_bg;
        
        average_offsets; % latest adjustment to x and y
        average_width; % of all stars in the stack
        average_background; % use this to get an indication of sky brightness
        average_variance; % use this to get an indication of noise in the image
        average_flux; % average flux of all stars, indicative of the sky transparency
        
    end
    
    properties(Hidden=true)
        
        start_index; % if non empty, use this number as initial index for run (e.g., to continue from where we stopped)
        
        use_refine_bg = 0; % use the cutouts_bg of each frame to estimate a different background for each point
        
        saturation_value = 50000; % consider any pixels above this to be saturated
        min_star_temp; % set a lower limit on temperature of stars for findStarsMAAT;
        num_phot_cutouts; % limit the number of cutouts given to photomery (in the fast cadence) to save runtime and RAM in lightcurves
        
        use_quick_find_stars = true; % use new method that runs faster
        use_mextractor = false; % use mextractor to identify stars and find their WCS and catalog mag/temp
        use_arbitrary_pos = false; % don't look for stars (e.g., when testing with the dome closed)
        
        lost_stars_fraction = 0.3;
        lost_flux_threshold = 0.3;
        max_failed_batches = 3; % if star flux is lost for more than this number of batches, quit the run
        
        camera_angle = 60; % degrees between image top and cardinal south/north (when after meridian)
        
        pass_source = {}; % parameters to pass to camera/reader/simulator
        pass_cal = {}; % parameters to pass to calibration object
        pass_back = {}; % parameters to pass to background object
        pass_phot = {}; % parameters to pass to photometry object
        pass_show = {}; % parameters to pass to show function
        
        prev_stack;        
        ref_stack;
        ref_positions;
        
        prev_fluxes; % fluxes measured in previous batch (for triggering)
        
        prev_average_width;
        
        batch_counter = 0;
        
        num_cosmic_rays = 0;
        
        use_verify_gui = 1; % when true (default) the GUI is loaded if it is not open during an acquisition
        
        star_props; % table with the results from quick_find_stars
        
        brake_bit = 1; % when this is set to 1 (using the GUI, for example), the run stops. 
        is_running = 0; % when this is 1, cannot start a new run or anything
        is_running_single = 0; % when this is 1, cannot start a new run or anything
        
        latest_command_str = '';
        latest_command_time = '';
        latest_command_pars = '';
        
        failed_batch_counter = 0;
        
        default_saturation_value;
        
        default_run_name;
        
        default_num_batches;
        default_batch_size;
        default_expT;
        default_frame_rate;
        default_roi_size;
        default_roi_center; 
        
        default_num_stars;
        default_cut_size;
        default_num_backgrounds;
        default_cut_size_bg;
        
        slow_mode_expT = 3;
        slow_mode_frame_rate = 1/3.01; % a little bit lower than 1/3
        slow_mode_batch_size = 1;
        
        fast_mode_expT = 0.03;
        fast_mode_frame_rate = 25;
        fast_mode_batch_size = 100;
        
        start_index_;
        use_background_;
        use_refine_bg_;
        use_remove_saturated_;
        saturation_value_;
        use_mextractor_;
        use_arbitrary_pos_;
        use_cutouts_;
        use_adjust_cutouts_;
        use_simple_photometry_;
        use_model_psf_;
        use_save_;
%         use_triggered_save_;
        use_show_;
        use_audio_;
        use_progress_;
        pass_source_;
        pass_cal_;
        pass_back_;
        pass_phot_;
        pass_show_;
        debug_bit_;
        log_level_;

        run_name_;
        batch_size_;
        total_runtime_;
        num_batches_;
        frame_rate_;
        expT_;
        num_stars_;
        cut_size_;
        num_backgrounds_;
        cut_size_bg_;
            
        use_roi_;
        roi_size_;
        roi_center_;
        
        show_what_list = {'images', 'stack', 'stack_proc'};
        
        drive_space_gb = []; 
        
        version = 1.07;
        
    end
    
    methods % constructor
        
        function obj = Acquisition(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'Acquisition')
                if obj.debug_bit>1, fprintf('Acquisition copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                
                if obj.debug_bit>1, fprintf('Acquisition constructor v%4.2f\n', obj.version); end
                
                obj.log = util.sys.Logger('Acquisition', obj);
                obj.log.heartbeat(600);
                
                obj.reader = file.Reader;
%                 obj.sim; % fill this when we have a simulator
                
                obj.cal = img.Calibration;
                obj.back = img.Background;
                obj.clip = img.Clipper;
                obj.clip.use_adjust = 0; % adjust through photometry object?
                obj.clip_bg = img.Clipper;
                obj.clip_bg.use_adjust = 0;
                
                obj.phot_stack = img.Photometry;
                obj.phot_stack.index = 2;
                
                obj.flux_buf = util.vec.CircularBuffer;
                
                obj.phot = img.Photometry;
                obj.phot.num_threads = 8;
                obj.lightcurves = img.Lightcurves;
                obj.light_stack = img.Lightcurves;
                
                % we don't want to start doing serious calibration on the
                % lightcurves object during acquisition! 
                obj.lightcurves.use_remove_outliers = 0;
                obj.lightcurves.use_zero_point = 0;
                obj.lightcurves.use_airmass_correction = 0;
                obj.lightcurves.use_psf_correction = 0;
                obj.lightcurves.use_sysrem = 0;
                
                obj.model_psf = img.ModelPSF;
                
                obj.buf = file.BufferWheel(10);
%                 obj.buf.product_type = 'Cutouts';
                obj.buf.use_save_raw_images = 0; % do not save the full frame images! 
                
                obj.runtime_buffer = util.vec.CircularBuffer;
                obj.runtime_buffer.titles = {'num_frames', 'time'};
                
                try 
                    obj.prog = util.sys.ProgressBar;
                    obj.audio = util.sys.AudioControl;
                catch ME
                    obj.log.input(['Warning: ' ME.getReport]);
                    warning(ME.getReport); 
                end
                
                obj.src = obj.reader;
                
                obj.deflator = file.Deflator;
                
                obj.sync = obs.comm.PcSync('server'); % no connect command is given, since this is still blocking indefinitely...
                obj.sync.reco.delay_time_minutes = 10; % this prevents the slow_timer from repeatedly trying to reconnect 
                
%                 obj.sync.name = 'Cam-PC';
                obj.setup_timer;
                obj.setup_slow_timer;
                obj.setup_backup_timer;
                
                obj.setupDefaults;
                
                obj.head = head.Header; % this also gives "head" to all sub-objects
                obj.cat = head.Catalog(obj.head);
                obj.cat.RA_scan_step = 2;
                obj.cat.RA_scan_range = 0;
                obj.cat.Dec_scan_step = 2;
                obj.cat.Dec_scan_range = 0;
                
                obj.lightcurves.head = obj.head;
                obj.lightcurves.cat = obj.cat;
                obj.light_stack.head = obj.head;
                obj.light_stack.cat = obj.cat;
                
                util.oop.save_defaults(obj); % make sure each default_XXX property is updated with the current XXX property value. 
                
                obj.stash_parameters;
                
%                 obj.sync.connect;
                
            end
            
        end

        function setupDefaults(obj)

            obj.runtime_units = 'hours';
            obj.total_runtime = 4;
%             obj.num_batches = 500;
            obj.batch_size = 100;
            
            obj.num_stars = 2000;
            obj.cut_size = 15;
            obj.avoid_edges = 50;
            
            obj.num_backgrounds = 50;
            obj.cut_size_bg = 32;

        end

    end
    
    methods % reset/clear
        
        function reset(obj)
            
            if obj.brake_bit==0
                warning('Must stop acquisition before resetting!');
                return;
            end

            list = properties(obj);
            
            for ii = 1:length(list)
                
                if isobject(obj.(list{ii})) && ~isempty(obj.(list{ii})) && ismethod(obj.(list{ii}), 'reset') 
                    if ~util.text.cs(list{ii}, 'sync', 'af')
                        obj.(list{ii}).reset;
                    end
                end
                
            end
            
            obj.num_stars_found = [];
            obj.prev_fluxes = [];
            obj.batch_counter = 0;
            obj.start_index = 1;
            obj.positions = [];
            obj.positions_bg = [];
            
            obj.failed_batch_counter = 0;
            
            obj.sync.outgoing.RA_rate_delta = 0;
            obj.sync.outgoing.DE_rate_delta = 0;

            obj.star_props = [];
            
            obj.forced_indices = [];
            obj.unlocked_indices = [];
            obj.dynamic_indices = []; 
            
            obj.num_cosmic_rays = 0;
            obj.micro_flares = [];
            
            obj.clear;
            
        end
        
        function resetForcedTargets(obj)
            
            obj.forced_targets = {};
            
        end
        
        function clear(obj)
            
            obj.prev_fluxes = obj.fluxes; % keep one batch history
            
            list = properties(obj);
            
            for ii = 1:length(list)
                
                if isobject(obj.(list{ii})) && ~isempty(obj.(list{ii})) && ismethod(obj.(list{ii}), 'clear') 
                    obj.(list{ii}).clear;
                end
                
            end
            
            clear@file.AstroData(obj); % clear everything in this object that is included in file.AstroData (including positions!)
            
            obj.prev_stack = obj.stack_proc; % keep one stack from last batch
            
            % clear some of the additional data products
            obj.stack_proc = [];
            obj.stack_cutouts = [];
            obj.stack_cutouts_bg = [];
            obj.cutouts_proc = [];
            obj.cutouts_bg = [];
            obj.cutouts_bg_proc = [];
            
        end
        
    end
    
    methods % getters

        function val = get.run_name(obj)
            
            if ~isempty(obj.head)
                val = obj.head.OBJECT;
            else
                val = [];
            end
            
        end
        
        function val = get.buf_full(obj)
            
            if isempty(obj.cam)
                val = [];
            else
                val = obj.cam.buffers;
            end 
            
        end
        
        function val = num_files(obj)
            
            if isa(obj.src, 'file.Reader')
                val = obj.src.num_files;
            else
                val = Inf;
            end
            
        end
        
        function val = get.num_batches(obj)
            
            if isprop(obj.src, 'num_batches')
                val = obj.src.num_batches;
             else
                val = [];
            end 
            
        end
        
        function val = get.batch_size(obj)
            
            if isprop(obj.src, 'batch_size')
                val = obj.src.batch_size;
            elseif isprop(obj.src, 'num_frames_per_batch')
                val = obj.src.num_frames_per_batch;
            else
                val = [];
            end 
            
        end
        
        function val = getFrameRateEstimate(obj)
            
            if ~isempty(obj.frame_rate) && ~isnan(obj.frame_rate)
                val = obj.frame_rate;
            elseif ~isempty(obj.frame_rate_average)
                val = obj.frame_rate_average;
            else
                val = [];
            end
            
        end
        
        function val = getTimeLeft(obj)
            
            if isempty(obj.getFrameRateEstimate)
                val = [];
            else
                val = (obj.num_batches-obj.batch_counter).*obj.batch_size./obj.getFrameRateEstimate;
            end
            
        end
        
        function val = getTimeLeftHMS(obj)
            
            t = obj.getTimeLeft;
            
            if isempty(t)
                val = '';
            else
                val = util.text.secs2hms(t);
            end
            
        end
        
        function val = getGbPerBatch(obj)
            
            stack_size = obj.src.ROI(3).*obj.src.ROI(4).*4;
            
            if isempty(obj.num_stars_found)
                cutout_size = obj.cut_size.^2.*obj.batch_size.*obj.num_stars*2;
                pos_size = 2.*obj.num_stars.*4;
                fluxes_size = obj.batch_size.*obj.num_stars.*4;
            else
                cutout_size = obj.cut_size.^2.*obj.batch_size.*obj.num_stars_found*2;
                pos_size = 2.*obj.num_stars_found.*4;
            fluxes_size = obj.batch_size.*obj.num_stars_found.*4;
            end
            
            background_size = obj.cut_size_bg.^2.*obj.batch_size.*obj.num_backgrounds.*2;
            pos_size_background = 2.*obj.num_backgrounds.*4;
            
            timestamps_size = obj.batch_size.*8;
            
            val = stack_size+cutout_size+pos_size+fluxes_size+background_size+pos_size_background+timestamps_size;
            val = val./1024.^3;
            
        end
        
        function val = getGbLeft(obj)
            
            val = obj.getGbPerBatch.*(obj.num_batches-obj.batch_counter);
            
        end
        
        function val = info_short(obj)
            
            val = sprintf('N= %d/%d batches | time left: %s | disk needed= %5.2f GB | frame rate= %4.2f Hz', ...
                obj.batch_counter, obj.num_batches, obj.getTimeLeftHMS, obj.getGbLeft, obj.frame_rate_average);
            
        end
        
        function val = info_long(obj)
            
            val = sprintf('Observation parameters\n--------------------------------------');
            
            if ~isempty(obj.head)
                
                val = sprintf('%s\n object: %s', val, obj.head.OBJECT);
                
                val = sprintf('%s\n RA:    %s hours\n Dec: %s deg', val, obj.head.RA, obj.head.Dec);
            
                val = sprintf('%s\n--------------------------------------', val);
                
                val = sprintf('%s\n telRA:    %s hours\n telDec: %s deg', val, obj.head.telRA, obj.head.telDec); 
                
                if obj.cat.success
                    val = sprintf('%s\n realRA:    %s hours\n realDec: %s  deg', val, head.Ephemeris.deg2hour(obj.cat.central_RA), head.Ephemeris.deg2sex(obj.cat.central_Dec)); 
                else
                    val = sprintf('%s\n Could not solve astrometry', val); 
                end
                
                val = sprintf('%s\n--------------------------------------', val);
            end
            
            val = sprintf('%s\n exp. time= %4.2f s \n frame rate= %4.2f Hz', val, obj.expT, obj.frame_rate);
            
            val = sprintf('%s\n num. batches= %d \n batch size= %d frames', val, obj.num_batches, obj.batch_size);
            
            val = sprintf('%s\n num. stars= %d / %d ', val, obj.num_stars_found, obj.num_stars);
            
            val = sprintf('%s\n cutout size= %dx%d pix \n edges= %d pix', val, obj.cut_size, obj.cut_size, obj.avoid_edges);
            
            val = sprintf('%s\n num.backgrounds= %d \n b/g cut. size= %dx%d pix', val, obj.num_backgrounds, obj.cut_size_bg, obj.cut_size_bg);
            
            val = sprintf('%s\n--------------------------------------', val);
            
            if ~isempty(obj.head) && ~isempty(obj.head.SCALE)
                if obj.use_model_psf
                    val = sprintf('%s\n seeing= %s"', val, util.text.print_vec(round(obj.head.SCALE.*sort([obj.average_width.*2.355 obj.model_psf.fwhm]),2), '-'));
                else
                    val = sprintf('%s\n seeing= %4.2f"', val, obj.head.SCALE.*obj.average_width.*2.355);
                end
            end
            
            if obj.use_model_psf
                val = sprintf('%s\n PSF widths= %4.2f / %4.2f pix', val, obj.model_psf.minor_axis, obj.model_psf.major_axis);
                val = sprintf('%s\n PSF angle= %4.2f deg', val, obj.model_psf.angle);
            end
            
            val = sprintf('%s\n LIMMAG_D= %4.2f', val, obj.head.LIMMAG_DET); 
            
            if length(obj.average_offsets)==2
                val = sprintf('%s\n dx/dy= %4.2f / %4.2f pix', val, obj.average_offsets(2), obj.average_offsets(1));
            end
            
            val = sprintf('%s\n mean flux= %.1f', val, obj.average_flux);
            val = sprintf('%s\n mean background= %.1f', val, obj.average_background);
            
            val = sprintf('%s\n AIRMASS= %5.3f', val, obj.head.AIRMASS); 
            val = sprintf('%s\n MOON DIST= %d', val, obj.head.ephem.moon_dist); 
            
            if isprop(obj.cam, 'focuser') && ~isempty(obj.cam.focus_pos_tip_tilt)
                val = sprintf('%s\n focus: %5.3f (%5.3f / %5.3f)', ...
                    val, obj.cam.focus_pos_tip_tilt(1), obj.cam.focus_pos_tip_tilt(2), obj.cam.focus_pos_tip_tilt(3));
            end
            
            val = sprintf('%s\n--------------------------------------', val);
            
            try
                val = sprintf('%s\n sensor temperature= %4.2f', val, obj.sensor_temperature);
            catch ME
                warning(ME.getReport); 
            end
            
            if isa(obj.src, 'obs.cam.Andor')

                if obj.src.use_async
                    val = sprintf('%s\n source: CAM-async', val);
                else
                    val = sprintf('%s\n source: CAM-sync', val);
                end
                
            elseif isa(obj.src, 'file.Reader')
                val = sprintf('%s\n source: reader', val);
            end
            
        end
        
        function val = get.expT(obj)
            
            if isprop(obj.src, 'expT')
                val = obj.src.expT;
            elseif ~isempty(obj.head)
                val = obj.head.expT;
            else
                val = [];
            end 
            
        end
        
        function val = get.frame_rate(obj)
            
            if isprop(obj.src, 'frame_rate')
                val = obj.src.frame_rate;
            elseif ~isempty(obj.head)
                val = obj.head.frame_rate;
            else
                val = [];
            end 
            
        end
        
        function val = get.use_roi(obj)
            
            if isprop(obj.src, 'use_roi')
                val = obj.src.use_roi;
            else
                val = [];
            end
            
        end
        
        function val = get.roi_size(obj)
            
            if isprop(obj.src, 'roi_size')
                val = obj.src.roi_size;
            elseif isprop(obj.src, 'im_size')
                val = obj.src.im_size;
            else
                val = [];
            end
            
        end
        
        function val = get.roi_center(obj)
            
            if isprop(obj.src, 'roi_center')
                val = obj.src.roi_center;
            elseif isprop(obj.src, 'center_region')
                val = obj.src.center_region;
            else
                val = [];
            end
            
        end
        
        function val = get.num_stars(obj)
            
            if ~isempty(obj.clip)
                val = obj.clip.num_stars;
            else
                val = [];
            end
            
        end
        
        function val = get.cut_size(obj)
            
            if ~isempty(obj.clip)
                val = obj.clip.cut_size;
            else
                val = [];
            end
            
        end
        
        function val = get.avoid_edges(obj)
            
            if ~isempty(obj.clip)
                val = obj.clip.avoid_edges;
            else
                val = [];
            end
            
        end
        
        function val = get.num_backgrounds(obj)
            
            if ~isempty(obj.clip_bg)
                val = obj.clip_bg.num_stars;
            else
                val = [];
            end
            
        end
        
        function val = get.cut_size_bg(obj)
            
            if ~isempty(obj.clip_bg)
                val = obj.clip_bg.cut_size;
            else
                val = [];
            end
            
        end
        
        function val = get.frame_rate_average(obj)
            
            if is_empty(obj.runtime_buffer)
                val = [];
            else
                N = sum(obj.runtime_buffer.data(:,1));
                T = sum(obj.runtime_buffer.data(:,2));
                val = N/T;
            end
        end
        
        function val = get.sensor_temperature(obj)
            
            if isa(obj.src, 'obs.cam.Andor')
                val = obj.src.getTemperatureHW;
            else
                val = [];
            end
            
        end
        
        function val = get.average_flux(obj)
            
            val = obj.phot_stack.average_flux;
            
        end
        
        function val = get.average_background(obj)
            
            val = obj.phot_stack.average_background;
            
        end
        
        function val = get.average_variance(obj)
            
            val = obj.phot_stack.average_variance;
            
        end
        
        function val = get.average_offsets(obj)
            
            val = [obj.phot_stack.average_offset_x obj.phot_stack.average_offset_y];
            
        end
        
        function val = get.average_width(obj)
            
%             if obj.use_model_psf
%                 val = (obj.model_psf.maj_axis+obj.model_psf.min_axis)/2;
%             else
                val = obj.phot_stack.average_width;
%             end
            
        end
        
        function val = getWidthEstimate(obj)
            
            if isempty(obj.prev_average_width) || ~isreal(obj.prev_average_width) || obj.prev_average_width<=0
                val = 1.6;
            elseif obj.prev_average_width>2
                val = 2;
            else
                val = obj.prev_average_width;
            end
            
        end
        
    end
    
    methods % setters
        
        function set.head(obj,val)
            
            obj.head = val;
            
            list = properties(obj);
            
            for ii = 1:length(list)
                
                if isobject(obj.(list{ii})) && ~isempty(obj.(list{ii})) && isprop(obj.(list{ii}), 'head') 
                    obj.(list{ii}).head = val;
                end
                
            end
            
        end
        
        function set.src(obj, val)
            
            if isempty(obj.src) || isempty(val)
                obj.src = val;
            elseif obj.src~=val
                num_batches = obj.src.num_batches;
                obj.src = val;
                obj.src.num_batches = num_batches;
            end
            
        end
        
        function set.run_name(obj, val)
            
            if ~isempty(obj.head)
                obj.head.OBJECT = val;
            end
            
        end
        
        function set.total_runtime(obj, val)
            
            obj.total_runtime = val;
            
            obj.num_batches = ceil(obj.total_runtime.*obj.convertRuntimeToSeconds.*obj.getFrameRateEstimate./obj.batch_size);
            
        end
        
        function set.num_batches(obj, val)
            
            if isprop(obj.src, 'num_batches')
                obj.src.num_batches = val;
            end 
            
        end
        
        function set.batch_size(obj, val)
            
            if isprop(obj.src, 'batch_size')
                obj.src.batch_size = val;
            elseif isprop(obj.src, 'num_frames_per_batch')
                obj.src.num_frames_per_batch = val;
            end 
            
        end
        
        function set.expT(obj, val)
            
            if isprop(obj.src, 'expT')
                obj.src.expT = val;
            end 
            
        end
        
        function set.frame_rate(obj, val)
            
            if isprop(obj.src, 'frame_rate')
                obj.src.frame_rate = val;
            end 
            
        end
        
        function set.use_roi(obj, val)
            
            if isprop(obj.src, 'use_roi')
                obj.src.use_roi = val;
            end
            
        end
        
        function set.roi_size(obj, val)
            
            if isprop(obj.src, 'roi_size')
                obj.src.roi_size = val;
            end
            
        end
        
        function set.roi_center(obj, val)
            
            if isprop(obj.src, 'roi_center')
                obj.src.roi_center = val;
            end
            
        end
        
        function set.num_stars(obj, val)
            
            if ~isempty(obj.clip)
                obj.clip.num_stars = val;
            end
            
        end
        
        function set.cut_size(obj, val)
            
            if ~isempty(obj.clip)
                obj.clip.cut_size = val;
            end
            
        end
        
        function set.avoid_edges(obj, val)
            
            if ~isempty(obj.clip)
                obj.clip.avoid_edges = val;
            end
            
        end
        
        function set.num_backgrounds(obj, val)
            
            if ~isempty(obj.clip_bg)
                obj.clip_bg.num_stars = val;
            end
            
        end
        
        function set.cut_size_bg(obj, val)
            
            if ~isempty(obj.clip_bg)
                obj.clip_bg.cut_size = val;
            end
            
        end
        
        function set.brake_bit(obj, val)
            
            obj.brake_bit = val;
            
            if obj.brake_bit && isprop(obj.src, 'brake_bit')
                obj.src.brake_bit = val;
            end
            
        end
        
        function set.use_show_gray(obj, val)
            
            obj.use_show_gray = val;
            
            if obj.use_show_gray
                colormap(obj.gui.axes_image, 'gray');
            else
                colormap(obj.gui.axes_image, 'default');
            end
            
        end
        
    end
    
    methods % utilities
        
        function chooseSource(obj, source)
            
            import util.text.cs;
            
            if nargin<2 || isempty(source)
                source = questdlg('Which source?', 'choose source', 'Reader', 'Camera', 'Simulator', 'Reader');
                if isempty(source)
                    return;
                end
                
                if cs(source, 'file reader', 'reader')
                    obj.chooseDir;
                end
                
            end
            
            if ischar(source)
                if cs(source, 'Zyla')
                    
                    if isempty(obj.cam) || ~isa(obj.cam.cam, 'obs.cam.ZylaControl')
%                         obj.cam = obs.cam.CameraControl('zyla');
                        obj.cam = obs.cam.Andor;
                        obj.cam.head = obj.head;
                    end
                    
                    obj.src = obj.cam;
                    
                elseif cs(source, 'Dhyana')
                    error('This camera is no longer supported...');
                    if isempty(obj.cam) || ~isa(obj.cam.cam, 'obs.DhyanaControl')
                        obj.cam = obs.CameraControl('dhyana');
                        obj.cam.head = obj.head;
                    end
                    
                    obj.src = obj.cam;
                    
                elseif cs(source, 'simcamera')
                    error('This camera is no longer supported...');
                    if isempty(obj.cam) || ~isa(obj.cam.cam, 'obs.cam.SimCamera')
                        obj.cam = obs.cam.CameraControl('sim');
                        obj.cam.head = obj.head;
                    end
                    obj.src = obj.cam;
                    
                    
%                 elseif cs(source, 'simulator')
%                     
%                     if isempty(obj.sim)
%                         obj.sim = img.Simulator;
%                         obj.sim.head = obj.head;
%                     end
%                     
%                     obj.src = obj.sim;

                elseif cs(source, {'file reader', 'reader'})
                    
                    if isempty(obj.reader)
                        obj.reader = file.Reader;
                        obj.cam.head = obj.head;
                    end
                    
                    obj.src = obj.reader;
                    
                elseif cs(source, 'camera')
                    if isempty(obj.cam)
%                         obj.cam = obs.cam.CameraControl;
                        obj.cam = obs.cam.Andor;
                        obj.cam.head = obj.head;
                    end
                    obj.default_expT = obj.cam.default_expT;
                    obj.default_frame_rate = obj.cam.default_frame_rate;
                    
                    obj.src = obj.cam;
                else
                    warning('unknown source "%s"', source);
                end
            else
                if isa(source, 'obs.cam.CameraControl') || isa(source, 'obs.cam.Andor') || isa(source, 'file.Reader') % add simulator check, too
                    obj.src = source;
                    obj.src.head = obj.head;
                else
                    warning(['unknown source class: ' class(source)]);
                end
            end
            
            obj.update;
            
        end
        
        function chooseDir(obj, dirname) % do I need this?
            
            if isempty(obj.src)
                error('Cannot chooseDir with an empty source');
            end
            
            if nargin<2 || isempty(dirname)
                dirname = '';
            end
            
            if ~isa(obj.src, 'file.Reader')
                error('Cannot chooseDir for a source of type %s', class(obj.src));
            end
            
            if isempty(dirname)
                obj.reader.browseDir;
            else
                if ~obj.reader.dir.cd(dirname)
                    error('cannot find directory %s', dirname);
                end
            end
            
            obj.cal.reader_dark.dir.cd(obj.reader.dir.pwd);
            obj.cal.reader_dark.dir.cd('..');
            if obj.cal.reader_dark.dir.smart_cd('dark')
                obj.cal.reader_dark.loadFiles;
            end
            
            obj.cal.reader_flat.dir.cd(obj.reader.dir.pwd);
            obj.cal.reader_flat.dir.cd('..');
            if obj.cal.reader_flat.dir.smart_cd('flat')
                obj.cal.reader_flat.loadFiles;
            end
            
            obj.cal.load;            
            
        end
        
        function setupSlowMode(obj)
            
            obj.expT = obj.slow_mode_expT;
            obj.frame_rate = obj.slow_mode_frame_rate;
            obj.batch_size = obj.slow_mode_batch_size;
            
            obj.total_runtime = obj.total_runtime; % this triggers the setter and updates num_batches
            
        end
        
        function setupFastMode(obj)
            
            obj.expT = obj.fast_mode_expT;
            obj.frame_rate = obj.fast_mode_frame_rate;
            obj.batch_size = obj.fast_mode_batch_size;
            
            obj.total_runtime = obj.total_runtime; % this triggers the setter and updates num_batches
            
        end
        
        function setupTestMode(obj)
            
            obj.use_arbitrary_pos = 1;
            obj.use_astrometry = 0;
            obj.use_adjust_cutouts = 0;
            obj.use_autodeflate = 0;
            obj.use_check_positions = 0;
            obj.use_sync_stop = 0;
                        
        end
        
        function cancelTestMode(obj)
            
            obj.use_arbitrary_pos = 0;
            obj.use_astrometry = 1;
            obj.use_adjust_cutouts = 1;
            obj.use_autodeflate = 1;
            obj.use_check_positions = 1;
%             obj.use_sync_stop = 1;
            
        end
        
        function createForcedTarget(obj, varargin)
            
            s = struct;
            
            for ii = 1:2:length(varargin)
                
                key = varargin{ii};
                
                if length(varargin)>ii
                    val = varargin{ii+1};
                else
                    val = NaN;
                end
                
                % what about strings??
%                 if ~isscalar(val)
%                     val = NaN;
%                 end
                
                s.(key) = val;

            end
            
            if ~isfield(s, 'RA') || ~isfield(s, 'Dec')
                error('Must input a target with RA and Dec'); 
            end
            
            if ischar(s.RA)
                s.RA = head.Ephemeris.hour2deg(s.RA); 
            end
            
            if ischar(s.Dec)
                s.Dec = head.Ephemeris.sex2deg(s.Dec);
            end
            
            if isempty(obj.forced_targets)
                obj.forced_targets{1} = s;
            else
                obj.forced_targets{end+1} = s;
            end
            
        end
        
    end
    
    methods % commands/calculations
        
        function unlock(obj) % manually remove all locks left over from crashes in mid-recording. Make sure camera is really done! 
            
            obj.is_running_single = 0;
            obj.is_running = 0;
            obj.brake_bit = 1;
            
        end
        
        function run(obj, varargin)
            
            if obj.is_running || obj.is_running_single
                disp('Already running, set is_running and is_running_single to zero...');
                return;
            else
                obj.is_running = 1;
            end
            
            check = obj.startup(varargin);
            
            try 
                
                cleanup = onCleanup(@obj.finishup);
                
                if check==0, return; end
            
                if obj.start_index
                    idx = obj.start_index;
                else
                    idx = 1;
                end

                for ii = idx:obj.num_batches

                    if obj.brake_bit
                        return;
                    end

                    obj.batch;
                    
                end

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function check = single(obj) % get a single batch and calculate stack_proc from it... 
            
            check = 0;
            
            if obj.is_running_single
                disp('Already running "single". Set is_running_single to zero...');
                return;
            else
                obj.is_running_single = 1;
            end
            
            obj.log.input('Getting single batch from source');
            
            try 
                
                obj.reset;
                
                obj.src.single; 
                
                obj.copyFrom(obj.src); % get the data into this object

                % if src is using ROI, must update the calibration object to do the same
                if isprop(obj.src, 'use_roi') && obj.src.use_roi 

                    obj.cal.use_roi = 1;

                    obj.cal.ROI = obj.src.ROI;

                else
                    obj.cal.use_roi = 0;
                end

                obj.calcStack;
                
                if ~isempty(obj.gui)
                    obj.show; 
                    obj.gui.update; 
                end
                
                check = 1;
                obj.is_running_single = 0;
                
            catch ME
                obj.log.error(ME.getReport);
                obj.is_running_single = 0;
                rethrow(ME);
            end
            
        end
        
        function startLiveView(obj)
            
            if ~isempty(obj.cam.gui) && obj.cam.gui.check
                figure(obj.cam.gui.fig.fig)
            else
                obj.cam.makeGUI;
            end
            
            obj.cam.live('autodyn', 1);
            
        end
        
        function success = runFlat(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('number', 20); % how many flat files we want en-total
            input.input_var('iterations', 4); % how many iterations of focus run we want to do
            input.scan_vars(varargin{:}); 
            
            success = 0;
            
            % find out how many flat files we already took
            
            on_cleanup = onCleanup(@obj.setBrake);
            obj.brake_bit = 0;
            
            for ii = 1:input.iterations

                if obj.brake_bit
                    break;
                end
                
                N_flats = 0; 

                if exist(fullfile(getenv('DATA_TEMP'), util.sys.date_dir('now')), 'dir')

                    d = util.sys.WorkingDirectory;
                    d.cd(fullfile(getenv('DATA_TEMP'), util.sys.date_dir('now')));

                    if d.cd('flat')

                        N_flats = length(d.match('*.h5*')); 

                    end

                end

                N = input.number-N_flats;

                if N>0

                    obj.cam.record('mode', 'flat', 'num_batches', N, 'batch_size', 100, 'frame_rate', 10, 'expT', 0.03); 

                else
                    success = 1;
                    break;
                end

            end
            
        end
        
        function runFocus(obj, varargin)
            
            obj.log.input('Running autofocus loop.');
            
            prev_focus_point = NaN;
            
            try
                
                obj.brake_bit = 0;
                
                if ~isempty(obj.gui) && obj.gui.check
                    obj.gui.latest_error = '';
                    obj.gui.latest_warning = '';
                    lastwarn('');
                    obj.gui.update;
                end
                
                success = 0;
                
                if isempty(obj.cam.focuser)
                    error('must be connected to camera and focuser!');
                end
                
                for ii = 1:obj.num_focus_iterations
                    
                    obj.sync.outgoing.report = 'Focusing';
                    obj.sync.update;
                    obj.parseCommands; % allow stopping between focus runs...
                    
                    if obj.brake_bit, return; end
                    
                    if obj.debug_bit, fprintf('Running focus loop attempt %d\n', ii); end
                    
                    if obj.use_cam_focusing
                        val = obj.cam.autofocus('iteration', ii); % run the autofocus loop inside the Andor class (no calibration...)
                    else
                        
                        obj.cam.af.cam = obj.cam;
                        
                        old_pos = obj.cam.focuser.pos;
                        
                        try % this runs a single camera focus loop
                            
                            obj.cam.is_running_focus = 1; % tell focus GUI to display an indicator that focus is running...
                            
                            % make sure finishup is called in the end
                            on_cleanup = onCleanup(@() obj.cam.finish_focus(old_pos));
                            
                            % update the object with the camera parameters
                            obj.cam.af.x_max = obj.cam.ROI(4);
                            obj.cam.af.y_max = obj.cam.ROI(3);
                            
                            obj.cam.af.pixel_size = obj.head.PIXSIZE;
                            
                            if ~isempty(obj.head.FIELDROT)
                                obj.cam.af.angle = obj.head.FIELDROT;
                            elseif util.text.cs(obj.head.cam_name, 'Zyla')
                                obj.cam.af.angle = -15;
                            elseif util.text.cs(obj.head.cam_name, 'Balor')
                                obj.cam.af.angle = -60;
                            end
                            
                            if obj.cam.is_running
                                disp('Camera is already running. Set is_running to zero...');
                                return;
                            end
                            
                            % the focus positions to scan
                            p = obj.cam.af.getPosScanValues(obj.cam.focuser.pos);
                            
                            if obj.cam.af.use_loop_back
                                p = [p flip(p)];
                            end
                            
                            obj.cam.af.reset;
                            obj.cam.af.iteration = ii;
                            
                            if isempty(obj.cam.af.gui) || ~obj.cam.af.gui.check
                                obj.cam.af.makeGUI;
                            end
                            
                            obj.cam.af.gui.update;
                            obj.cam.af.plot;
                            
                            ax = obj.cam.af.gui.axes_image;
                            
                            text(ax, mean(ax.XLim), mean(ax.YLim), 'Finding stars for focus run!', 'FontSize', 36, 'HorizontalAlignment', 'Center');
                            
                            drawnow;
                            
                            figure(obj.cam.af.gui.fig.fig);
                            
                            % find stars, the quick version!
                            obj.cam.single('frame rate', obj.cam.af.frame_rate, 'exp time', obj.cam.af.expT, 'batch size', obj.cam.af.batch_size); % need to first update with observational parameters (using varargin!)
                            obj.num_sum = size(obj.cam.images,3);
                            obj.stack = single(sum(obj.cam.images,3));
                            obj.stack_proc = obj.cal.input(obj.stack, 'num', obj.num_sum); % calibrated sum
                            obj.star_props = util.img.quick_find_stars(obj.stack_proc, 'threshold', obj.cam.af.threshold,...
                                'saturation', 5e4*obj.num_sum, 'unflagged', 1, 'num_stars', obj.cam.af.num_stars);
                            
                            if obj.brake_bit, return; end
                            
                            if isempty(obj.star_props)
                                error('Could not find any stars for doing focus!');
                            end
                            
                            obj.cam.focuser.pos = p(1);
                            
                            % before we start the loop
                            obj.cam.startup('async', 0, 'num_batches', length(p), 'batch_size', obj.cam.af.batch_size, ...
                                'reset', 1, 'frame_rate', obj.cam.af.frame_rate, 'exp time', obj.cam.af.expT, ...
                                'save', 0, 'show', 1, 'log_level', 1, 'progress', 0, 'audio', 1, varargin{:});
                            
                            for jj = 1:obj.num_batches
                                
                                if obj.brake_bit, return; end
                                if obj.cam.brake_bit, return; end
                                
                                obj.cam.focuser.pos = p(jj);
                                
                                obj.cam.batch;
                                obj.images = obj.cam.images;
                                
                                obj.cutouts = util.img.mexCutout(obj.images, obj.star_props.pos, obj.cam.af.cut_size, NaN);
                                obj.cutouts_proc = obj.cal.input(obj.cutouts, 'clip', obj.star_props.pos);
                                
                                C = nansum(obj.cutouts_proc,3);
                                
                                obj.cam.af.phot_struct = util.img.photometry2(C, 'aperture', obj.cam.focus_aperture, 'use_aperture', 0, ...
                                    'use_forced', 1, 'gauss_sigma', 5, 'use_gaussian', 1, 'index', 3, 'threads', 4);
                                
                                A = obj.cam.af.phot_struct.forced_photometry.area;
                                B = obj.cam.af.phot_struct.forced_photometry.background;
                                fluxes = obj.cam.af.phot_struct.forced_photometry.flux - A.*B;
                                %                     widths = phot_struct.apertures_photometry.width;
                                
                                widths = util.img.fwhm(C-permute(B, [1,3,4,2]), 'method', 'filters', 'min_size', 0.25, 'max_size', 10, 'step', 0.2)/2.355;
                                widths(widths>10 | widths<0.1) = NaN;
                                
                                if jj==1
                                    [~,idx] = nanmax(fluxes);
                                    obj.cam.af.star_idx = idx;
                                end
                                
                                obj.cam.af.input(jj, p(jj), widths, (abs(fluxes)), obj.star_props.pos, C);
                                
                                obj.cam.af.plot;
                                obj.cam.af.showCutout;
                                
                                drawnow;
                                
                            end % for jj
                            
                            obj.cam.af.calculate;
                            obj.cam.af.fitSurface;
                            obj.cam.af.findPosTipTilt;
                            obj.cam.af.plot;
                            
                            if ~isempty(obj.cam.af.gui)
                                obj.cam.af.gui.update;
                            end
                            
                            fprintf('FOCUSER RESULTS: width= %f | pos= %f | tip= %f | tilt= %f\n', obj.cam.af.found_width, obj.cam.af.found_pos, obj.cam.af.found_tip, obj.cam.af.found_tilt);
                            
                            if ~isnan(obj.cam.af.found_pos)
                                obj.cam.focuser.pos = obj.cam.af.found_pos;
                                if obj.cam.af.use_fit_tip_tilt
                                    obj.cam.focuser.tipRelativeMove(obj.cam.af.found_tip);
                                    obj.cam.focuser.tiltRelativeMove(obj.cam.af.found_tilt);
                                end
                            else
                                disp('The location of new position is NaN. Choosing original position');
                                obj.cam.focuser.pos = old_pos;
                            end
                            
                            %                 if isprop(obj.cam.focuser, 'tip') && ~isempty(obj.af.found_tip)
                            %                     obj.cam.focuser.tip = obj.cam.focuser.tip + obj.af.found_tip;
                            %                 end
                            %
                            %                 if isprop(obj.cam.focuser, 'tilt') && ~isempty(obj.af.found_tilt)
                            %                     obj.cam.focuser.tilt = obj.cam.focuser.tilt + obj.af.found_tilt;
                            %                 end
                            
                            success = 1;
                            
                        catch ME
                            
                            obj.cam.focuser.pos = old_pos;
                            obj.log.error(ME.getReport);
                            rethrow(ME);
                            
                        end
                        
                        
                    end % focus using Analysis tools, not internally in the camera class
                    
                    if obj.brake_bit || val==0
                        break;
                    end
                    
                    min_pos = obj.cam.af.pos(2);
                    max_pos = obj.cam.af.pos(end-1);
                    
                    if obj.debug_bit, fprintf('Resulting focus point is %4.2f at FWHM of %4.2f"\n', obj.cam.af.found_pos, obj.cam.af.found_width.*2.355.*obj.head.SCALE); end
                    
                    if obj.cam.af.found_width<1 && obj.cam.af.found_pos>=min_pos && obj.cam.af.found_pos<=max_pos
                        success = 1;
                        break; % if focus is good enough and not at the edges of the range, we don't need to repeat it
                    end
                    
                    if abs(obj.cam.af.found_pos - prev_focus_point)<0.05
                        success = 1;
                        break; % if focus returns to the same position each time, we can give up on it
                    end
                    
                end % for ii (iterations)
                
                % could not find a decent focus:
                if ~isempty(obj.cam.af.gui) && obj.cam.af.gui.check
                    
                    obj.cam.af.plot;
                    
                    obj.cam.af.gui.update;
                    
                    if success
                        util.plot.inner_title(sprintf('Focus success after %d iterations!', ii),...
                            'ax', obj.cam.af.gui.axes_image, 'FontSize', 26, 'Position', 'North');
                    else
                        util.plot.inner_title(sprintf('Focus failed after %d iterations!', ii),...
                            'ax', obj.cam.af.gui.axes_image, 'FontSize', 26, 'Position', 'North');
                        pause(0.5);
                    end
                    
                end
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
            obj.brake_bit = 1;
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.update;
            end
            
        end
        
    end
    
    methods % timer related
        
        function stop_timers(obj)
            
            stop(obj.backup_timer);
            stop(obj.slow_timer); 
            stop(obj.timer); 
            
        end
        
        function start_timers(obj)
            
            obj.setup_backup_timer;
            obj.setup_slow_timer;
            obj.setup_timer;
            
        end
        
        function setup_timer(obj, ~, ~)
            
            if ~isempty(obj.timer) && isa(obj.timer, 'timer') && isvalid(obj.timer)
                if strcmp(obj.timer.Running, 'on')
                    stop(obj.timer);
                    delete(obj.timer);
                    obj.timer = [];
                end
            end
            
            delete(timerfind('name', 'acquisition-timer'));
            
            obj.timer = timer('BusyMode', 'drop', 'ExecutionMode', 'fixedRate', 'Name', 'acquisition-timer', ...
                'Period', 30, 'StartDelay', 30, 'TimerFcn', @obj.callback_timer, 'ErrorFcn', @obj.setup_timer, 'BusyMode', 'drop');
            
            start(obj.timer);
            
        end
        
        function callback_timer(obj, ~, ~)
            
%             disp('timer')
            
            if obj.brake_bit
                obj.updateSyncData; % this also updates the obs_log
                if ~isempty(obj.gui) 
                    obj.gui.update;
                end
            elseif obj.is_running==0 && obj.is_running_single==0 && obj.cam.is_running==0 && obj.cam.is_running_focus==0
                obj.sync.outgoing.report = 'idle';
            end
            
        end
        
        function setup_slow_timer(obj, ~, ~)
            
            if ~isempty(obj.slow_timer) && isa(obj.slow_timer, 'timer') && isvalid(obj.slow_timer)
                if strcmp(obj.slow_timer.Running, 'on')
                    stop(obj.slow_timer);
                    delete(obj.slow_timer);
                    obj.slow_timer = [];
                end
            end
            
            delete(timerfind('name', 'acquisition-slow-timer'));
            
            obj.slow_timer = timer('BusyMode', 'drop', 'ExecutionMode', 'fixedRate', 'Name', 'acquisition-slow-timer', ...
                'Period', 300, 'StartDelay', 300, 'TimerFcn', @obj.callback_slow_timer, 'ErrorFcn', @obj.setup_slow_timer);
            
            start(obj.slow_timer);
            
            
        end
        
        function callback_slow_timer(obj, ~, ~)
            
            if ~strcmp(obj.timer.Running, 'on')
                obj.setup_timer;
            end
            
            obj.sync.update;
            pause(0.05); 
            
            if obj.sync.status==0 && obj.brake_bit % do not stop to reconnect during acquisition! 
%                 disp('connecting to PcSync'); 
                obj.connectSync;
                
                if obj.sync.status
                    obj.sync.reco.lock;
                end
                
                pause(0.05);
               
            elseif obj.brake_bit
                obj.sync.setup_callbacks;
            end
            
        end
        
        function setup_backup_timer(obj, ~, ~)
            
            if ~isempty(obj.backup_timer) && isa(obj.backup_timer, 'timer') && isvalid(obj.backup_timer)
                if strcmp(obj.backup_timer.Running, 'on')
                    stop(obj.backup_timer);
                    delete(obj.backup_timer);
                    obj.backup_timer = [];
                end
            end
            
            delete(timerfind('name', 'acquisition-backup-timer'));
            
            obj.backup_timer = timer('BusyMode', 'drop', 'ExecutionMode', 'fixedRate', 'Name', 'acquisition-backup-timer', ...
                'Period', 1800, 'StartDelay', 1800, 'TimerFcn', @obj.callback_backup_timer, 'ErrorFcn', @obj.setup_backup_timer);
            
            start(obj.backup_timer);
            
            
        end
        
        function callback_backup_timer(obj, ~, ~)
            
            if ~strcmp(obj.slow_timer.Running, 'on')
                obj.setup_slow_timer;
            end
            
        end
        
        function updateSyncData(obj)
            
%             disp('updateSyncData'); 
            
            try 
                
                if obj.use_sync && obj.sync.status

                    drawnow;
                    
                    if obj.brake_bit==0
                        obj.sync.read_data_rx;
                        obj.sync.read_data_tx;
                    end
                    
                    obj.sync.update;

                    s = obj.sync.incoming;
                    
                    list = head.Header.makeSyncList; 

                    if obj.use_ignore_sync_object_name
                        list = list(~strcmp(list, 'OBJECT'));
                    end
                    
                    if obj.brake_bit==0 % things we want to NOT update during a run? 
                        list = list(~strcmpi(list, 'OBJECT'));
                        list = list(~strcmpi(list, 'RA'));
                        list = list(~strcmpi(list, 'DEC'));
                        list = list(~strcmpi(list, 'RA_DEG'));
                        list = list(~strcmpi(list, 'DEC_DEG'));
                    end
                    
                    if ~isempty(s) && isstruct(s)

                        for ii = 1:length(list)
                            if isfield(s, list{ii})
                                try
                                    obj.head.(list{ii}) = s.(list{ii});
                                end
                            end
                        end

                    end

                    % we have two ways to update the log: fast and slow
                    if obj.brake_bit % slow mode, just go over all files and check the data

                        obj.obs_log = obj.makeObsLog;

                    else % fast mode (while taking images, only update the latest obs_log, don't go scanning files)

                        if isempty(obj.obs_log) || ~isfield(obj.obs_log, 'date') || ~strcmp(obj.obs_log.date, obj.buf.date_dir) % only do the slow update if there is no obs_log saved or if it is outdated
                            obj.obs_log = obj.makeObsLog;
                        else

                            if ~isfield(obj.obs_log, obj.run_name)
                                obj.obs_log.(obj.run_name) = struct('name', obj.run_name, 'start', '', 'end', '', 'runtime', [], 'num_files', []);
                            end

                            s_obs = obj.obs_log.(obj.run_name); % get the structs for this run name
                            s_obs(end).runtime = obj.t_end_stamp;
                            s_obs(end).end_time = obj.t_end;
                            s_obs(end).num_files = obj.batch_counter;
                            obj.obs_log.(obj.run_name) = s_obs; % structs are not handles! 

                        end

                    end

                    obj.sync.outgoing.obs_log = obj.obs_log; % update the Manager on how much observing time is invested in each target
                    obj.drive_space_gb = obj.getDriveSpace; % calculate how much space is left in each drive (in Gb) 
                    obj.sync.outgoing.drives = obj.drive_space_gb;

                    if obj.cat.success
                        obj.sync.outgoing.obsRA = head.Ephemeris.deg2hour(obj.cat.central_RA); 
                        obj.sync.outgoing.obsDec = head.Ephemeris.deg2sex(obj.cat.central_Dec); 
                    end
                    
                    if ~isempty(obj.cam) && ~isempty(obj.cam.af) && obj.cam.af.success
                        obj.sync.outgoing.focus_pos = obj.cam.af.found_pos;
                        obj.sync.outgoing.focus_width = obj.cam.af.found_width;
                    end
                    
                    if isfield(obj.sync.outgoing, 'report') && strcmp(obj.sync.outgoing.report, 'Running')
                        obj.sync.outgoing.batch_counter = obj.batch_counter;
                        obj.sync.outgoing.total_batches = obj.num_batches;
                        obj.sync.outgoing.runtime = obj.prog.getElapsed; 
                    end
                    
                    if obj.use_sync_stop && isfield(s, 'stop_camera') && s.stop_camera
                        obj.brake_bit = 1; % allow for manager to send stop command to camera... 
                    end
                    
                    if obj.brake_bit==0 && obj.batch_counter>1
                        obj.sync.outgoing.report = 'Running';
                    end
                    
                    obj.parseCommands;
                    
                    obj.sync.update;
                
                end
                
            catch ME
                warning(ME.getReport)
            end
            
        end
        
        function parseCommands(obj)
            
%             disp('parseCommands'); 
            
            import util.text.cs;
            
            try 
                
                if ~isfield(obj.sync.incoming, 'command_str')
                    return;
                end
                
                obj.sync.outgoing.echo_str = obj.sync.incoming.command_str;
                obj.sync.outgoing.echo_time = obj.sync.incoming.command_time;
                
%                 fprintf('command_str= %s | command_time= %s\n', obj.sync.incoming.command_str, obj.sync.incoming.command_time); 
                
                obj.sync.read_data(obj.sync.hndl_rx, 'rx'); 
                obj.sync.read_data(obj.sync.hndl_tx, 'tx'); 
                
                if ~isempty(obj.sync.incoming.command_str) && ~isempty(obj.sync.incoming.command_time)
                    
                    t1 = util.text.str2time(obj.latest_command_time);
                    t2 = util.text.str2time(obj.sync.incoming.command_time);
                    
                    if isempty(obj.latest_command_time)
                        dt = 0;
                    else
                        dt = abs(minutes(t2-t1)); % the time since we got the last command needs to be long enough
                    end

                    if strcmp(obj.sync.incoming.command_str, obj.latest_command_str) && dt<1 % cannot take the same command over and over in such a short interval...
                        return;
                    end 
                    
                    obj.latest_command_str = obj.sync.incoming.command_str;
                    obj.latest_command_time = obj.sync.incoming.command_time;
                    obj.latest_command_pars = obj.sync.incoming.command_pars;
                    
                    if cs(obj.sync.incoming.command_str, 'start') && ... % GOT COMMAND TO START A NEW RUN
                            obj.brake_bit && obj.is_running==0 && obj.is_running_single==0 &&...
                            obj.cam.is_running==0 && obj.cam.is_running_focus==0

                        args = util.text.parse_inputs(obj.latest_command_pars);
                        
                        input = util.text.InputVars;
                        input.input_var('focus', obj.use_focus_on_start, 'use_focus');
                        input.input_var('mode', 'fast', 'cam_mode', 'camera_mode'); 
                        input.input_var('exp_time', [], 'exposure_time');
                        input.input_var('frame_rate', []); 
                        input.scan_vars(args{:}); 
                        
                        obj.sync.outgoing.error = ''; 
                        obj.sync.outgoing.report = 'Starting';
                        obj.sync.outgoing.batch_counter = 0;
                        obj.sync.outgoing.total_batches = 0;
                        obj.sync.outgoing.runtime = 0;
                        obj.sync.update;
                        pause(0.1); % leave time to update
                        
                        if input.focus 
                            disp('Now running focus by order of dome-PC'); % this message will be removed later on...
                            obj.runFocus;
                        end
                        
                        obj.log.input(sprintf('Starting run command from Dome-PC. Args= "%s"', obj.latest_command_pars)); 
                        disp(obj.log.report); 
                        if cs(input.mode, 'fast')
                            obj.setupFastMode;
                        elseif cs(input.mode, 'slow')
                            obj.setupSlowMode;
                        end
                        
                        if ~isempty(input.exp_time)
                            args{end+1} = 'expT'; % this is the only format that is understood by run()
                            args{end+1} = input.exp_time;
                            
                            if isempty(input.frame_rate)
                                input.frame_rate = (1./input.exp_time).*0.99; % make the frame rate a little lower than the expected
                            end
                            
                            args{end+1} = 'frame_rate';
                            args{end+1} = input.frame_rate; 
                        
                        end
                        
                        obj.use_save = 1; % verify that we are saving images, unless asked not to by dome-PC
                        
                        if ~isempty(obj.gui) && obj.gui.check
                            figure(obj.gui.fig.fig); % pop the GUI back on top
                        end
                        
                        obj.run('reset', 1, args{:}); % the rest of the inputs from dome-pc are parsed in the regular way
                        
                    elseif cs(obj.sync.incoming.command_str, 'start') % only get here if the timing is wrong... 
                        warning('Got command to start running, but already running...'); 
                    elseif cs(obj.sync.incoming.command_str, 'stop') 
                        
                        obj.brake_bit = 1;
                        
                        obj.sync.outgoing.report = 'idle';

                    else
                        error('Unknown command: %s! Use "start" or "stop", etc...', obj.sync.incoming.command_str); 
                    end

                    % careful: there is a "return" statement in there
                    
                    obj.sync.outgoing.report = 'idle'; % even if we started a new run, it will have to be over before we can set report back to "idle"
                
                end
                
            catch ME
                obj.sync.outgoing.error = sprintf('error! \n%s', ME.getReport); 
                rethrow(ME); 
            end

            
        end
        
    end
    
    methods(Hidden=true) % internal functions
        
        function input = makeInputVars(obj, varargin)
            
            idx = util.text.InputVars.isInputVars(varargin);
            
            if ~isempty(varargin) && any(idx) % was given an InputVars object
            
                idx = find(idx, 1, 'first');
                input = varargin{idx}; 
                varargin(find(idx, 1, 'first')) = []; % remove the one cell with InputVars in it
            
            else % use default values (load them from Acquisition object)
            
                input = util.text.InputVars;
                input.input_var('use_reset', false, 'reset'); 
                input.input_var('start_index', []);
                input.input_var('use_background', []);
                input.input_var('use_refine_bg', []);
                input.input_var('use_remove_saturated', []);
                input.input_var('saturation_value', []);
                input.input_var('use_mextractor', []);
                input.input_var('use_arbitrary_pos', []);
                input.input_var('use_cutouts', []);
                input.input_var('use_adjust_cutouts', []);
                input.input_var('use_simple_photometry', []);
                input.input_var('use_model_psf', []);
                input.input_var('use_save', [], 'save');
                input.input_var('use_trigger_save', []);
                input.input_var('use_show', [], 'show');
                input.input_var('use_audio', []);
                input.input_var('use_progress', []);
                input.input_var('pass_source', {}, 7); % cell array to pass to camera/reader/simulator
                input.input_var('pass_cal', {}, 7); % cell array to pass to calibration object
                input.input_var('pass_back', {}, 7); % cell array to pass to background object
                input.input_var('pass_phot', {}, 7); % cell array to pass to photometry object
                input.input_var('pass_show', {}, 7) % cell array to pass to show function
                input.input_var('debug_bit', []);
                input.input_var('log_level', []);
                
%                 input.input_var('run_name', '', 'name', 'object', 'objname');
                input.input_var('RA', [], 'right ascention', 'right ascension');
                input.input_var('DE', [], 'declination');
                input.input_var('expT', [], 'T', 'exposure time');
                input.input_var('frame_rate', []); 
                input.input_var('num_batches', [], 'Nbatches');
                input.input_var('total_runtime', [], 'runtime');
                input.input_var('batch_size', [], 'frames');
                input.input_var('num_stars', [], 'Nstars');
                input.input_var('cut_size', []);
                input.input_var('num_backgrounds', [], 'Nbackgrounds');
                input.input_var('cut_size_bg', []);
                
                input.scan_obj(obj); % overwrite defaults using values in Acquisition object

            end
            
            input.scan_vars(varargin{:});
            
        end
        
        function stash_parameters(obj, input) % sets all camera parameters into hidden "stash" parameters. If given an InputVars object, will load parameters from it to the camera object.
            
            obj.start_index_ = obj.start_index;
            obj.use_background_ = obj.use_background;
            obj.use_refine_bg_ = obj.use_refine_bg;
            obj.use_remove_saturated_ = obj.use_remove_saturated;
            obj.saturation_value_ = obj.saturation_value;
            obj.use_mextractor_ = obj.use_mextractor;
            obj.use_arbitrary_pos_ = obj.use_arbitrary_pos;
            obj.use_cutouts_ = obj.use_cutouts;
            obj.use_adjust_cutouts_ = obj.use_adjust_cutouts;
            obj.use_simple_photometry_ = obj.use_simple_photometry;
            obj.use_model_psf_ = obj.use_model_psf;
            obj.use_save_ = obj.use_save;
%             obj.use_triggered_save_ = obj.use_triggered_save;
            obj.use_show_ = obj.use_show;
            obj.use_audio_ = obj.use_audio;
            obj.use_progress_ = obj.use_progress;
            obj.pass_source_ = obj.pass_source;
            obj.pass_cal_ = obj.pass_cal;
            obj.pass_back_ = obj.pass_back;
            obj.pass_phot_ = obj.pass_phot;
            obj.pass_show_ = obj.pass_show;
            obj.debug_bit_ = obj.debug_bit; 
            obj.log_level_ = obj.log_level;
            
            obj.run_name_ = obj.run_name;
            obj.batch_size_ = obj.batch_size;
            obj.total_runtime_ = obj.total_runtime;
            obj.num_batches_ = obj.num_batches;
            obj.frame_rate_ = obj.frame_rate;
            obj.expT_ = obj.expT;
            obj.num_stars_ = obj.num_stars;
            obj.cut_size_ = obj.cut_size;
            obj.num_backgrounds_ = obj.num_backgrounds;
            obj.cut_size_bg_ = obj.cut_size_bg;
            
            obj.use_roi_ = obj.use_roi;
            obj.roi_size_ = obj.roi_size;
            obj.roi_center_ = obj.roi_center;
            
            % optionally, give the "input" object properties into the public properties of "obj"
            if nargin>1 && ~isempty(input) && isa(input, 'util.text.InputVars')
                list = properties(input);
                for ii = 1:length(list)
                    if isprop(obj, list{ii})
                        obj.(list{ii}) = input.(list{ii});
                    end
                end
            end
            
        end
        
        function unstash_parameters(obj) % return "stashed" parameters to camera object after run is done. 
            
            obj.start_index = obj.start_index_;
            obj.use_background = obj.use_background_;
            obj.use_refine_bg = obj.use_refine_bg_;
            obj.use_remove_saturated = obj.use_remove_saturated_;
            obj.saturation_value = obj.saturation_value_;
            obj.use_mextractor = obj.use_mextractor_;
            obj.use_arbitrary_pos = obj.use_arbitrary_pos_;
            obj.use_cutouts = obj.use_cutouts_;
            obj.use_adjust_cutouts = obj.use_adjust_cutouts_;
            obj.use_simple_photometry = obj.use_simple_photometry_;
            obj.use_model_psf = obj.use_model_psf_;
            obj.use_save = obj.use_save_;
%             obj.use_triggered_save = obj.use_triggered_save_;
            obj.use_show = obj.use_show_;
            obj.use_audio = obj.use_audio_;
            obj.use_progress = obj.use_progress_;
            obj.pass_source = obj.pass_source_;
            obj.pass_cal = obj.pass_cal_;
            obj.pass_back = obj.pass_back_;
            obj.pass_phot = obj.pass_phot_;
            obj.pass_show = obj.pass_show_;
            obj.debug_bit = obj.debug_bit_; 
            obj.log_level = obj.log_level_;
            
            obj.run_name = obj.run_name_;
            obj.batch_size = obj.batch_size_;
            obj.total_runtime = obj.total_runtime_;
            obj.num_batches = obj.num_batches_;
            obj.frame_rate = obj.frame_rate_;
            obj.expT = obj.expT_;
            obj.num_stars = obj.num_stars_;
            obj.cut_size = obj.cut_size_;
            obj.num_backgrounds = obj.num_backgrounds_;
            obj.cut_size_bg = obj.cut_size_bg_;
            
            obj.use_roi = obj.use_roi_;
            obj.roi_size = obj.roi_size_;
            obj.roi_center = obj.roi_center_;
            
        end
        
        function update(obj, input)
            
%             disp('update');
            
            if nargin>=2 && ~isempty(input) && isa(input, 'util.text.InputVars')
                
                if ~isempty(input.RA)
                    obj.head.RA = input.RA;
                end

                if ~isempty(input.DE)
                    obj.head.DE = input.DE;
                end

            end
            
            if obj.use_sync && obj.use_ignore_manager==0 && ~isempty(obj.sync)
                obj.updateSyncData;
            end
            
            obj.head.update;
            
            if isa(obj.src, 'obs.cam.Andor')
                obj.head.INST = obs.cam.mex_new.get(obj.src.hndl, 'name'); 
                obj.head.PIXSIZE = obs.cam.mex_new.get(obj.src.hndl, 'pixel width'); 
                obj.head.NAXIS1 = obs.cam.mex_new.get(obj.src.hndl, 'width');
                obj.head.NAXIS2 = obs.cam.mex_new.get(obj.src.hndl, 'height');
                obj.head.NAXIS3 = obj.batch_size;
                obj.head.NAXIS4 = obj.num_stars;
            end
            
            obj.cal.camera_name = obj.head.INST;
            
        end
        
        function check = startup(obj, varargin)
            
            check = 0;
            
            % must make this shorter or break it down to smaller functions
            try 
                
                if obj.brake_bit==0
                    disp('Cannot start a new acquisition while old one is still runnning (turn off brake_bit)');
                    return;
                end
                
                input = obj.makeInputVars(varargin{:});
                
                obj.update(input); % update header object to current time and input run name, RA/DE if given to input.
                
                if obj.use_sync && obj.sync.status==0
                    obj.gui.latest_error = 'Cannot start a run without a PcSync connection to dome-PC!';
                    warning(obj.gui.latest_error); 
                    obj.log.error(obj.gui.latest_error); 
                    return;
                end
                
                if isempty(obj.head.OBJECT)
                    obj.gui.latest_error = 'Cannot start a run without a valid OBJECT name in header!';
                    warning(obj.gui.latest_error); 
                    obj.log.error(obj.gui.latest_error); 
                    return;
                end
                
                obj.stash_parameters(input);
            
                if obj.use_arbitrary_pos && obj.use_astrometry
                    obj.gui.latest_error = 'Cannot use astrometry and arbitrary pos at the same time!';
                    warning(obj.gui.latest_error); 
                    obj.log.error(obj.gui.latest_error); 
                    return;
                end
                
                obj.buf.use_save_photometry = obj.use_save_photometry; 
                
                if isempty(obj.num_batches)
                    obj.gui.latest_error = 'Must input a number of batches!';
                    warning(obj.gui.latest_error);
                    obj.log.error(obj.gui.latest_error);
                    return;
                end
                
                obj.sync.hndl_rx.BytesAvailableFcn = '';
                obj.sync.hndl_tx.BytesAvailableFcn = '';
                
                if ~isempty(obj.gui) && obj.gui.check
                    
                    if input.use_reset
                        obj.gui.latest_message = sprintf('Starting new run "%s" with %d batches. Save is %d', obj.run_name, obj.num_batches, obj.use_save);
                    else
                        obj.gui.latest_message = sprintf('Continuing new run "%s" with %d batches. Save is %d', obj.run_name, obj.num_batches, obj.use_save);
                    end
                    
                    obj.gui.update;
                    
                end
                
                if input.log_level
                    obj.log.input(sprintf('Starting a new run "%s" (saving is %d). ', obj.run_name, obj.use_save)); 
                end
                
                if ~obj.cal.checkDark
                    obj.gui.latest_error = 'Cannot start a new run without loading darks into calibration object!';
                    warning(obj.gui.latest_error);
                    obj.log.error(obj.gui.latest_error);
                    return;
                end
                
                if obj.getTimeLeft>3600*10
                    obj.gui.latest_error = sprintf('Run scheduled to take %4.2f hours with these parameter... aborting!', obj.getTimeLeft/3600);
                    warning(obj.gui.latest_error);
                    obj.log.error(obj.gui.latest_error);
                    return;
                end
                
                if obj.use_save && obj.getGbLeft>util.sys.disk_space(obj.buf.directory)*1.0 % only throw an error if the required disk space is bigger than storage! 
                    obj.gui.latest_error = sprintf('Run scheduled requires an estimated %5.2f Gb of storage. Only %5.2f Gb available on drive!', obj.getGbLeft, util.sys.disk_space(obj.buf.directory));
                    warning(obj.gui.latest_error);
                    obj.log.error(obj.gui.latest_error);
                    return;
                end
                
                if input.use_reset % this parameter is not saved in the object because we only use it here... 
                    
                    obj.reset;
                    
                    if obj.debug_bit, disp(['Starting run "' obj.run_name '" for ' num2str(obj.num_batches) ' batches.']); end

                end
                
                % update the obs_log with the new run
                if obj.use_save && input.use_reset
                    
                end
                
                if obj.display_num_rect_stars>100
                    obj.display_num_rect_stars = 100; % do not plot too many rectangles, as it slows down the display
                end
                
                if obj.use_save && obj.use_autodeflate
                    obj.deflator.setup_timer;
                end
                
                if isempty(obj.positions)
                    
                    str = 'Positions field empty. Calling single then findStars';
                    
                    if ~isempty(obj.gui)
                        pause(1); 
                        obj.gui.latest_message = str;
                        obj.gui.update; 
                    end
                    
                    if obj.debug_bit, disp(str); end
                    
                    obj.sync.outgoing.report = 'Finding stars';
                    obj.sync.update;
                    
                    obj.single;
                    obj.findStars;
                    
                    if obj.use_astrometry
                        
                        obj.sync.outgoing.report = 'Astrometry';
                        obj.sync.update;
                        obj.runAstrometry;
                        
                        % add forced cutouts
                        if obj.cat.success
                            obj.addForcedPositions; 
                        end
                        
                    end
                    
                end
                
                num_frames = obj.num_batches.*obj.batch_size;
                num_stars = size(obj.positions,1);
                num_apertures = length(obj.phot.aperture).*obj.phot.use_aperture + length(obj.phot.aperture).*obj.phot.use_forced; % can later change the second length() to the forced aperture list
                
                obj.lightcurves.startup(num_frames, num_stars, num_apertures); 
                
%                 obj.update(input); % update header object to current time and input run name, RA/DE if given to input.
                
                obj.head.RUNSTART = util.text.time2str(obj.head.ephem.time);
                obj.head.NAXIS4 = obj.num_stars_found;
                
                if obj.src.use_roi
                    obj.buf.product_type_append = 'ROI';
                else
                    obj.buf.product_type_append = strrep(obj.buf.product_type_append, 'ROI', '');
                end
                
                if obj.use_save && input.use_reset % only do these things when startin a real run (with saving)
                    
                    try % get the observation log updated that a new run has begum

                        if ~isempty(obj.obs_log)
                            obj.obs_log = obj.makeObsLog;
                        end

                        start = util.text.time2str(datetime('now', 'TimeZone', 'UTC')); 

                        if ~isfield(obj.obs_log, obj.run_name) % this run name has not been created yet
                            obj.obs_log.(obj.run_name) = struct('name', obj.run_name, 'start', start, 'end', '', 'runtime', 0, 'num_files', 0); % make a new struct for this name
                        else % there are previous runs with this name, need to 
                            obj.obs_log.(obj.run_name) = vertcat(obj.obs_log.(obj.run_name), struct('name', obj.run_name, 'start', start, 'end', '', 'runtime', 0, 'num_files', 0)); 
                        end

                    catch ME
                        warning(ME.getReport); 
                    end
                    
                    try % make a folder for the new files
                        
                        basename = obj.buf.makeFullpath; % the name from the object_name, ignoring the override
                        
                        for ii = 1:100
                            
                            dirname = sprintf('%s%s_run%d', basename, obj.run_name_append, ii);
                            if ~exist(dirname, 'dir')
                                mkdir(dirname);
                                obj.buf.directory_override = dirname;
                                break;
                            end
                            
                        end
                        
                    catch ME
                        warning(ME.getReport);
                    end

                    if ~isempty(obj.cat.success) && obj.cat.success % save the catalog file 

                        try
                            filename = fullfile(obj.buf.directory, 'catalog.mat');
                            obj.cat.saveMAT(filename);
                        catch ME
                            warning(ME.getReport);
                        end

                    end

                end

                if obj.use_save % save the README file for each time the run is started/continued
                    
                    try
                        filename = obj.buf.getReadmeFilename;
                        util.oop.save(obj, filename, 'name', 'acquisition', 'hidden', false); 
                    catch ME
                        warning(ME.getReport); 
                    end
                    
                end
                
                if obj.use_audio, try obj.audio.playTakeForever; catch ME, warning(ME.getReport); end, end

                obj.src.num_batches = obj.num_batches;
                
                if isa(obj.src, 'file.Reader') % we need to change the reader interface to be more inline with the camera...
                    obj.src.num_frames_per_batch = obj.batch_size;
                    obj.src.num_files_per_batch = 1;
                    % what if batch_size is bigger than 100??
                end
                
                obj.src.startup('use_save', 0, 'use_reset', input.use_reset, 'use_async', 1, obj.pass_source{:});

                if obj.use_progress
                    obj.prog.start(obj.num_batches); % maybe use continue if not restarting? 
                end
                
                obj.brake_bit = 0;
                
                check = 1;
                
            catch ME
                obj.log.error(ME.getReport);
                obj.unstash_parameters;
                obj.is_running = 0;
                obj.sync.outgoing.error = ME.getReport; 
                obj.sync.outgoing.report = 'idle';
                rethrow(ME);
            end
            
        end
        
        function finishup(obj)
            
            if obj.use_progress
                obj.prog.finish;
            end
            
            obj.src.finishup;
                        
            if obj.use_save
                
                try % readme file
                    filename = obj.buf.getReadmeFilename('Z');
                    util.oop.save(obj, filename, 'name', 'acquisition', 'hidden', false); 
                catch ME
                    warning(ME.getReport);
                end
                
                try % save lightcurves
                
                    if obj.use_save_stack_lcs
                        obj.light_stack.saveAsMAT(fullfile(obj.buf.directory, 'lightcurves.mat'));
                    end
                
                catch ME
                    warning(ME.getReport);
                end
                
                try % save any microflares
                    
                    if ~isempty(obj.micro_flares)
                        
                        micro_flares = obj.micro_flares;
                        head = obj.head;
                        
                        save(fullfile(obj.buf.directory, 'micro_flares.mat'), 'micro_flares', 'head', '-v7.3'); 
                        
                    end
                    
                catch ME
                    warning(ME.getReport); 
                end
                
            end
            
            if obj.debug_bit, disp(['Finished run "' obj.run_name '" with ' num2str(obj.batch_counter) ' batches.']); end
            
            if obj.use_audio
                try
                    obj.audio.playShowsOver;
                catch ME
                    warning(ME.getReport);
                end
            end
            
            obj.brake_bit = 1;
            
            obj.unstash_parameters;
            
            obj.buf.directory_override = '';
            
            obj.sync.outgoing.RA_rate_delta = 0;
            obj.sync.outgoing.DE_rate_delta = 0;
            
            obj.sync.setup_callbacks;
            
            obj.sync.outgoing.report = 'idle'; % maybe only do this if there was no error? 
            
            obj.lightcurves.finishup;
            
            obj.is_running = 0;
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.update;
            end
            
        end
        
        function batch(obj)
            
            if obj.src.is_finished
                obj.brake_bit = 1;
                return;
            end
            
            obj.head.update;
            
            t0 = tic;
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.update; 
            end
            
            obj.prev_stack = obj.stack_proc; % keep one stack from last batch
            obj.clear;
            obj.positions = obj.clip.positions; % recover the positions from the copy saved in the Clipper
            
            t_batch = tic;
            obj.src.batch; % produce the data (from camera, file, or simulator)
            t_batch = toc(t_batch); 
            
            drawnow; % make sure commands to the GUI and other callbacks are noticed... 
            
            obj.update;
            
            t_copy = tic;
            obj.copyFrom(obj.src); % get the data into this object
            t_copy = toc(t_copy);
            
            if obj.debug_bit>1, fprintf('Starting batch %d. Loaded %d images from "%s" source.\n', obj.batch_counter+1, size(obj.images,3), class(obj.src)); end
            
            obj.head.END_STAMP = obj.t_end_stamp;
                        
            J = juliandate(util.text.str2time(obj.t_end));
            
            obj.juldates = J + (obj.timestamps - obj.t_end_stamp)/24/3600; 
            
            % if src is using ROI, must update the calibration object to do the same
            if isprop(obj.src, 'use_roi') && obj.src.use_roi 
                
                obj.cal.use_roi = 1;
                
                obj.cal.ROI = obj.src.ROI;
                
            else
                obj.cal.use_roi = 0;
            end
            
            t_stack = tic;
            obj.calcStack;
            t_stack = toc(t_stack);
            
            if obj.use_cutouts
                
                t_cut = tic;
                obj.calcCutouts;
                t_cut = toc(t_cut);
            
                t_light = tic;
                obj.calcLightcurves;
                t_light = toc(t_light); 
            
            end
            
            obj.positions = obj.clip.positions;
            obj.positions_bg = obj.clip_bg.positions;
            
            t_save = tic;
            
            if obj.use_save
                
                obj.buf.input(obj);
%                 obj.buf.clearImages; % do we need this if we have set use_save_raw_images=0 in the buffers?
                obj.buf.save;
                obj.buf.nextBuffer;
                % check triggering then call the camera save mode
                
            end
            
            t_save = toc(t_save);
            
            if ismethod(obj.src, 'next')
                obj.src.next;
            end
            
            obj.batch_counter = obj.batch_counter + 1;
            
            if obj.use_progress
                obj.prog.showif(obj.batch_counter);
            end
            
            obj.start_index = obj.start_index + 1;
            
            t_show = tic;
            
            if obj.use_show
                obj.show;
            end
            
            t_show = toc(t_show);
            
            t_gui = tic;
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.update; 
            end
            
            drawnow;
            
            t_gui = toc(t_gui);
            
            obj.runtime_buffer.input([size(obj.images,3), toc(t0)]);
            
            if obj.use_print_timing
                fprintf('TIMING: batch= %4.2f | copy= %4.2f | stack= %4.2f | cut= %4.2f | light= %4.2f | save= %4.2f | show= %4.2f | gui= %4.2f\n', ...
                    t_batch, t_copy, t_stack, t_cut, t_light, t_save, t_show, t_gui); 
            end
            
        end
        
        function calcStack(obj)
            
            % make the basic stack image
            obj.num_sum = size(obj.images,3);
%             obj.stack = util.stat.sum_single(obj.images); % sum along the 3rd dimension directly into single precision
            try
                
                obj.stack = single(sum(obj.images,3)); % sum along the 3rd dimension directly into single precision
                
            catch ME
                if strcmp(ME.identifier, 'MATLAB:nomem')
                    disp('Ran out of memory while stacking images. Trying again with a split array...');
                    pause(0.1);
                    
                    indices = ceil([size(obj.images,1) size(obj.images,2)]./2); 
                
                    obj.stack(1:indices(1), 1:indices(2)) = single(sum(obj.images(1:indices(1), 1:indices(2),:),3)); 
                    obj.stack(1:indices(1), indices(2)+1:end) = single(sum(obj.images(1:indices(1), indices(2)+1:end,:),3)); 
                    obj.stack(indices(1)+1:end,1:indices(2)) = single(sum(obj.images(indices(1)+1:end,1:indices(2),:),3)); 
                    obj.stack(indices(1)+1:end,indices(2)+1:end) = single(sum(obj.images(indices(1)+1:end,indices(2)+1:end,:),3)); 
                    
                else
                    rethrow(ME);
                end
            
            end
            
            obj.stack_proc = obj.cal.input(obj.stack, 'sum', obj.num_sum); % stack after calibration
            
            % make the background cutouts of the stack 
            if isempty(obj.clip_bg.positions) % only if we didn't already assign positions to the bg_cutouts
                obj.clip_bg.arbitraryPositions('im_size', size(obj.stack));
            end
            
            obj.stack_cutouts_bg = obj.clip_bg.input(obj.stack_proc); % dim 1&2 are y&x, dim 3 is scalar, dim 4 is star number.

            if obj.use_background
            
                obj.back.input(obj.stack_cutouts_bg, obj.clip_bg.positions);
                % should also get variance from background object...
                
                B = obj.back.getImage(size(obj.stack_proc));
                obj.stack_proc = obj.stack_proc - B;

            end
            
            if ~isempty(obj.positions) % if no positions are known just skip this part
                
                if obj.use_dynamic_cutouts
                    
                    mask = imdilate(obj.stack_proc>1024*1, ones(5)); % dilate the area around stars
                    mask = logical(obj.cal.dark_mask + mask); % add the bad pixels to the mask
                    pos = util.img.find_cosmic_rays(obj.images, mask, 2*256, 6, 100, 0); % the arguments are: images, mask, threshold, num_threads, max_number, debug_bit
                    
                    if ~isempty(pos)
                        pos = sortrows(pos, 4, 'descend'); 
                    end
                    
                    obj.num_cosmic_rays = obj.num_cosmic_rays + size(pos,1); % keep track of how many such events we found
                    
                    new_pos = ones(obj.num_dynamic_cutouts,1).*floor(size(obj.stack)/2); 
                    
                    for ii = 1:min(size(pos,1), obj.num_dynamic_cutouts)
                        new_pos(ii,1:2) = pos(ii,1:2); 
                    end
                    
                    obj.positions(obj.dynamic_indices,:) = [];
                    obj.dynamic_indices = []; 
                    obj.dynamic_indices = size(obj.positions,1)+1:size(obj.positions,1)+size(new_pos,1); 
                    
                    obj.positions = vertcat(obj.positions, new_pos); 
                    
                    obj.clip.positions = obj.positions; 
                    
                end

                obj.stack_cutouts = obj.clip.input(obj.stack_proc);  
                
                obj.phot_stack.input(obj.stack_cutouts, 'positions', obj.clip.positions, 'timestamps', obj.timestamps(1),'juldates', obj.juldates(1)); % run photometry on the stack to verify flux and adjust positions
                if ~isempty(obj.phot_stack.gui) && obj.phot_stack.gui.check, obj.phot_stack.gui.update; end
                
                obj.light_stack.getData(obj.phot_stack)

                obj.prev_average_width = obj.average_width; % keep track of the average width
                
                if obj.use_check_positions && obj.use_arbitrary_pos==0
                    obj.checkRealign;
                end
                
                % store the latest fluxes from the stack cutouts
                obj.flux_buf.input(obj.phot_stack.fluxes);

                if obj.use_adjust_cutouts
                    offsets = obj.average_offsets;
                    offsets(isnan(offsets)) = 0;
                    obj.clip.positions = double(obj.clip.positions + offsets);
                else
                    
                end
                
                if obj.batch_counter>2 && obj.use_sync && obj.use_autoguide && obj.sync.status && ~isempty(obj.getFrameRateEstimate)
                    % send the average adjustment back to mount controller (should we still adjust the cutouts though??)
                    rot = [cosd(obj.camera_angle) sind(obj.camera_angle); -sind(obj.camera_angle) cosd(obj.camera_angle)];
                    vec = rot*(obj.average_offsets.*obj.head.SCALE./obj.batch_size.*obj.getFrameRateEstimate)'; % units of arcsec/second
                    vec = vec*0.7;
                    dRA = vec(1)/15; % convert from arcsec to sidereal seconds
                    dDE = vec(2);
                    
                    obj.sync.outgoing.RA_rate_delta = dRA;
                    obj.sync.outgoing.DE_rate_delta = dDE;
                    
%                     fprintf('dx= %6.4f | dy= %6.4f | dRA= %6.4f | dDE= %6.4f\n', obj.average_offsets(1), obj.average_offsets(2), dRA, dDE); 
                    
                    obj.sync.update;
                     
                end

                if obj.use_model_psf && ~isempty(obj.stack_cutouts)
                    
                    obj.model_psf.input(obj.stack_cutouts, obj.phot_stack.offsets_x, obj.phot_stack.offsets_y);
                    
%                     obj.prev_average_width = obj.average_width;
                    
                end
                
            end
            
        end
        
        function findStars(obj)
            
            S = obj.stack_proc;

            if obj.use_remove_saturated % can we get rid of this and instead use the saturation avoidance in quick_find_stars?
                mu = median(squeeze(util.stat.corner_mean(util.img.jigsaw(obj.stack_proc))));
                sig = median(squeeze(util.stat.corner_std(util.img.jigsaw(obj.stack_proc))));
                S = util.img.remove_saturated(S, 'saturation', 4.5e4*obj.num_sum, 'threshold', mu+5*sig, 'dilate', 4); % note the saturation value is X100 because we are looking at the stack
            end
            
            if obj.use_remove_bad_pixels
                S = util.img.maskBadPixels(S, NaN); 
            end
            
            if obj.use_arbitrary_pos
                obj.clip.arbitraryPositions('im_size', size(obj.stack)); % maybe add some input parameters?
                obj.positions = obj.clip.positions;
            elseif obj.use_mextractor
                obj.findStarsMAAT;
            elseif obj.use_quick_find_stars
                % replaced the psf width with 1, because getWidthEstimate was returning unreasonable results... 
                T = util.img.quick_find_stars(S, 'psf', 1, 'number', obj.num_stars, 'sigma', obj.detect_thresh, ...
                    'dilate', obj.cut_size-5, 'saturation', obj.saturation_value.*obj.num_sum, 'edges', obj.avoid_edges, 'unflagged', 1); 
                if isempty(T)
                    error('Could not find any stars using quick_find_stars!');
                end
                
                obj.head.THRESH_DETECTION = obj.detect_thresh; 
                
                obj.clip.positions = T.pos;
                obj.positions = T.pos;
                obj.num_stars_found = size(obj.clip.positions,1);
                
                obj.star_props = T; 
                
            else
                obj.clip.findStars(S);
            end
            
            obj.ref_stack = obj.stack_proc;
            obj.ref_positions = obj.clip.positions;
            
        end
        
        function runAstrometry(obj)
            
            if obj.debug_bit, disp('runAstrometry'); end
            
            if ~isempty(obj.gui)
                obj.gui.latest_message = 'running astrometry...';
                obj.gui.update;
            end
            
            obj.cat.detection_threshold = obj.detect_thresh;
            obj.cat.detection_stack_number = obj.num_sum;
            obj.cat.detection_exposure_time = obj.expT;
            
            if isfield(obj.sync.incoming, 'TELRA')
                obj.head.TELRA = obj.sync.incoming.TELRA;
                obj.head.TELRA_DEG = obj.sync.incoming.TELRA_DEG;
            end
            
            if isfield(obj.sync.incoming, 'TELDEC')
                obj.head.TELDEC = obj.sync.incoming.TELDEC;
                obj.head.TELDEC_DEG = obj.sync.incoming.TELDEC_DEG;
            end
            
            if ~isempty(obj.star_props)
                obj.cat.input(obj.star_props); % new method gives additiona data to calculate limiting magnitude etc... 
            else
                obj.cat.inputPositions(obj.positions); % old method
            end
            
            if ~isempty(obj.cat.data) && ~isempty(obj.cat.success) && obj.cat.success % successfully filled the catalog
                
                obj.positions = obj.cat.positions; % if used some filter on the stars we found
                obj.ref_positions = obj.cat.positions;
                obj.clip.positions = obj.cat.positions;
                
                str = sprintf('Successfully solved astrometry! Coordinates: %s %s', ...
                    head.Ephemeris.deg2hour(obj.cat.central_RA), head.Ephemeris.deg2sex(obj.cat.central_Dec)); % need to pull the coordinates in a more serious way
                
                obj.cat.num_stars = obj.num_stars;
                
                obj.positions = obj.cat.positions; % usually we will already have positions so this should do nothing (unless this analysis is on full frame rate images)
                
            else
                str = sprintf('Could not fit astrometric solution...');
            end
            
            if ~isempty(obj.gui)
                obj.gui.latest_message = str;
                obj.gui.update;
            end
            
            if obj.debug_bit, disp(str); end
            obj.head.LIMMAG_DETECTION = obj.cat.detection_limit;
            
            [obj.object_idx, dist] = obj.cat.findNearestObject;
            
            if dist>5/3600
                obj.object_idx = []; % if the closest star found is more than 5" from the required position, it is not really a good match!
            end
            
        end
        
        function addForcedPositions(obj)
            
            % check if object position is already included in the list of targets
            
            found = 0;
            tol = 3; % arcsec
            
            for ii = 1:length(obj.forced_targets)
                
                RA = obj.head.RA_DEG;
                DE = obj.head.DEC_DEG;
            
                if abs(obj.forced_targets{ii}.RA - RA)*3600<tol && abs(obj.forced_targets{ii}.Dec - DE)*3600<tol % if an existing target is very close to this forced target
                    found = ii;
                    break;
                end
                
            end
            
            if found==0 % if not, add it now! 
%                 obj.createForcedTarget('RA', obj.head.RA, 'Dec', obj.head.Dec); 
                xy = obj.cat.coo2xy(obj.head.RA, obj.head.DEC); 
                if ~(xy(1)<1 || xy(1)>obj.head.NAXIS2 || xy(2)<1 || xy(2)>obj.head.NAXIS1)
                    obj.positions = vertcat(obj.positions, xy); % add the position fields
                    obj.forced_indices(end+1) = size(obj.positions,1); % log the indices of the new positions into the forced_indices field
                    if obj.debug_bit, fprintf('adding a forced cutout number %d at x= %f | y= %f\n', obj.forced_indices(end), xy(1), xy(2)); end
                end
            end

%             default_table = table;
%             
%             props = obj.cat.data.Properties.VariableNames;
% 
%             for ii = 1:length(props)
%                 
%                 if iscell(obj.cat.data{1, props{ii}})
%                     default_table.(props{ii}) = {NaN};
%                 else
%                     default_table.(props{ii}) = NaN(size(obj.cat.data{1, props{ii}}), 'Like', obj.cat.data{1, props{ii}});
%                 end
%                 
%             end

            for ii = 1:length(obj.forced_targets)
                
                % first get the X/Y from the RA/Dec
                xy = obj.cat.coo2xy(obj.forced_targets{ii}.RA, obj.forced_targets{ii}.Dec); 
                
                if xy(1)<1 || xy(1)>obj.head.NAXIS2 || xy(2)<1 || xy(2)>obj.head.NAXIS1
                    continue; % if target is outside field of view, skip it
                end
                
                obj.positions = vertcat(obj.positions, xy); % add the position fields
                obj.forced_indices(end+1) = size(obj.positions,1); % log the indices of the new positions into the forced_indices field
                
                % this is commented as I've decided I don't want to update
                % the catalog with the forced/dynamic targets... 
                
                % add the RA/Dec and X/Y to the catalog
%                 t = default_table; 
%                 
%                 list = fields(obj.forced_targets{ii});
%                 
%                 for jj = 1:length(list)
%                     if any(strcmp(list{jj}, obj.cat.data.Properties.VariableNames))
%                         t.(list{jj}) = obj.forced_targets{ii}.(list{jj});
%                     end
%                 end
%                 
%                 t.pos = xy;
%                 
%                 obj.cat.data = vertcat(obj.cat.data, t);
%                 
%                 obj.cat.magnitudes = obj.cat.data.Mag_BP;
%                 obj.cat.coordinates = [obj.cat.data.RA, obj.cat.data.Dec]; 
%                 obj.cat.temperatures = obj.cat.data.Teff;
                
            end
            
            obj.clip.positions = obj.positions; % update the clipper with the new forced positions
            obj.ref_positions = obj.positions; 
            
        end
        
        function findStarsMAAT(obj)
            
            if isempty(which('mextractor'))
                error('Cannot load the MAAT package. Make sure it is on the path...');
            end
             
            % add additional tests to remove irrelvant stars
            
            obj.cat.threshold = obj.detect_thresh;
            obj.cat.input(obj.stack_proc);
            
            T = obj.cat.data;
            
            if obj.min_star_temp
                T = T(T{:,'Teff'}>=obj.min_star_temp,:); % select only stars with temperature above minimal level (hotter stars have smaller angular scale)
            end
            
            T = sortrows(T, 'Mag_G'); % sort stars from brightest to faintest
            
            if height(T)>obj.num_stars
                T = T(1:obj.num_stars,:);
            end
            
            obj.num_stars_found = height(T);
            
            obj.positions = [T.XPEAK_IMAGE T.YPEAK_IMAGE];
            obj.clip.positions = obj.positions;
            
            obj.magnitudes = T{:,'Mag_G'};
            obj.coordinates = [T.RA T.Dec];
            
            obj.obj_idx = obj.cat.findNearestObject;
            
        end
        
        function val = checkFluxes(obj)
            
            if size(obj.flux_buf.data, 2)~=size(obj.phot_stack.fluxes,2)
                obj.flux_buf.reset;
            end
            
            if is_empty(obj.flux_buf)
                val = 1;
            else
                
                mean_fluxes = obj.flux_buf.mean;
                mean_fluxes(mean_fluxes<=0) = NaN;

                new_fluxes = obj.phot_stack.fluxes;
                new_fluxes(isnan(mean_fluxes)) = [];
                mean_fluxes(isnan(mean_fluxes)) = [];

                flux_lost = nnz(new_fluxes<obj.lost_flux_threshold*mean_fluxes)>obj.lost_stars_fraction*numel(mean_fluxes); % lost <fraction> of the flux in more than <fraction> of the stars...
                % want to add more tests...?

                val = ~flux_lost;
                
            end
            
        end
        
        function checkRealign(obj)
            
            if size(obj.flux_buf.data, 2)~=size(obj.phot_stack.fluxes,2)
                obj.flux_buf.reset;
            end
            
            if obj.checkFluxes % check that stars are still aligned properly... 
                
                if obj.use_sync, obj.sync.outgoing.stars_visible = 1;  end
                
            else
                
                if obj.debug_bit, disp('Lost star positions, using quick_align'); end
                    
                [~,shift] = util.img.quick_align(obj.stack_proc, obj.ref_stack);
                obj.clip.positions(1:size(obj.ref_positions,1),:) = double(obj.ref_positions + flip(shift));

                % this shift should also be reported back to mount controller? 

                obj.stack_cutouts = obj.clip.input(obj.stack_proc);

                obj.phot_stack.input(obj.stack_cutouts, 'positions', obj.positions); % run photometry on the stack to verify flux and adjust positions
                if ~isempty(obj.phot_stack.gui) && obj.phot_stack.gui.check, obj.phot_stack.gui.update; end

                % second test, after quick_align, to verify the stars are back!
                if obj.checkFluxes
                    
                    obj.failed_batch_counter = 0;
                    
                    if obj.use_sync, obj.sync.outgoing.stars_visible = 1;  end
                    
                else
                    
                    obj.failed_batch_counter = obj.failed_batch_counter + 1;
                    
                    if obj.failed_batch_counter>obj.max_failed_batches
                        
                        fprintf('Cannot find stars %d times in a row. Quiting run...\n', obj.failed_batch_counter);
                        
                        obj.brake_bit = 1; % finish this batch and then quit the run
                        
                        if obj.use_sync, obj.sync.outgoing.stars_visible = 0;  end
                        
                    end
                    
                end
                
            end
            
        end
        
        function calcCutouts(obj)
            
            obj.cutouts = obj.clip.input(obj.images);
            
            obj.cutouts_proc = obj.cal.input(obj.cutouts, 'clip', obj.clip);
            
            obj.cutouts_bg = obj.clip_bg.input(obj.images);
            
            obj.cutouts_bg_proc = obj.cal.input(obj.cutouts_bg, 'clip', obj.clip_bg);
            
            if obj.use_background
                
                B = obj.back.getPoints(obj.clip.positions); 
                B = permute(B, [4,3,2,1]); % turn the column vector into a 4D vector
                % can also get variance from background object...

                if obj.use_refine_bg
                    % use bg_cutouts to calculate overall differences between
                    % frames to correct for regional results from background
                    % object (based on the stack). 
                end

                obj.cutouts_proc = obj.cutouts_proc - B/obj.num_sum;

            end
            
        end
        
        function calcLightcurves(obj)
            
            if isempty(obj.num_phot_cutouts) || any(isinf(obj.num_phot_cutouts)) || any(isnan(obj.num_phot_cutouts))
                C = obj.cutouts_proc; 
                P = obj.positions;
            elseif ischar(obj.num_phot_cutouts) && util.text.cs(obj.num_phot_cutouts, 'all')
                C = obj.cutouts_proc;
                P = obj.positions;
            elseif isnumeric(obj.num_phot_cutouts) && obj.num_phot_cutouts<=size(obj.cutouts_proc,4)
                C = obj.cutouts_proc(:,:,:,1:obj.num_phot_cutouts); % in this case we only take some of the stars into the photometry/lightcurves pipeline
                P = obj.positions(1:obj.num_phot_cutouts,:); 
            else
                C = obj.cutouts_proc;
                P = obj.positions;
            end
            
            if obj.use_simple_photometry
                
                obj.fluxes = permute(util.stat.sum2(C), [3,4,2,1]); 
                
            else % use extensive photometry method
            
                obj.phot.input('images', C, 'timestamps', obj.timestamps, 'positions', P, 'juldates', obj.juldates); % add variance input? 
                if ~isempty(obj.phot.gui) && obj.phot.gui.check, obj.phot.gui.update; end

                obj.fluxes = obj.phot.fluxes;
                obj.errors = obj.phot.errors;
                obj.areas = obj.phot.areas;
                obj.backgrounds = obj.phot.backgrounds;
                obj.variances = obj.phot.variances;
                obj.offsets_x = obj.phot.offsets_x;
                obj.offsets_y = obj.phot.offsets_y;
                obj.centroids_x = obj.phot.centroids_x;
                obj.centroids_y = obj.phot.centroids_y;
                obj.widths = obj.phot.widths;
                obj.bad_pixels = obj.phot.bad_pixels;
                obj.flags = obj.phot.flags;
                
                if obj.use_store_photometry
                    obj.lightcurves.getData(obj.phot);
                end
                
            end
            
            if obj.use_dynamic_cutouts
                
                for ii = 1:length(obj.dynamic_indices)
                    
                    i2 = obj.dynamic_indices(ii); 
                    
                    f = obj.fluxes(:,i2);
                    
                    % get some statistics on the lightcurve
                    S = nanstd(f); 
                    M = nanmean(f); 
                    [mx, idx] = nanmax(f); 
                    N = nnz((f-M)./S>5);
                    
                    % get the 1st moments
                    C = obj.cutouts_proc(:,:,:,i2); 
                    C2 = C - median(C(:)); 
                    
                    [X,Y] = meshgrid(floor(-size(C,2)/2)+1:floor(size(C,2)/2));
                    
                    cx = squeeze(sum(X.*abs(C2))./sum(abs(C2))); 
                    cy = squeeze(sum(Y.*abs(C2))./sum(abs(C2))); 
                    
                    if N>=3 % at least 3 points above 5 sigma...
                        
                        st = struct('filename', obj.buf.filename, 'batch_index', obj.batch_counter+1, 'frame_index', idx, ...
                            'peak', mx, 'mean', M, 'std', S, 'num_frames', N, 'flux', f, 'cx', cx, 'cy', cy, ...
                            'cutouts', obj.cutouts_proc(:,:,:,i2)); 
                        
                        if isempty(obj.micro_flares)
                            obj.micro_flares = st;
                        else
                            obj.micro_flares(end+1) = st;
                        end
                        
                    end
                    
                end
                
            end
            
            if obj.lightcurves.gui.check
                obj.lightcurves.gui.update;
            end
            
        end
        
        function s = makeObsLog(obj, date) % scan the folders in "data_temp" to get up-to-date info on number of files and runtime for each target
            
            if nargin<2 || isempty(date)
                date = obj.buf.date_dir; 
            end
            
            s = struct;
            s.date = date;
            
            d = util.sys.WorkingDirectory(obj.buf.base_dir);
            if exist(fullfile(d.pwd, date), 'dir')
                d.cd(date);
            elseif exist(fullfile(getenv('DATA_EXTRAS'), date), 'dir')
                d.cd(fullfile(getenv('DATA_EXTRAS'), date))
            else
                s = struct('date', date); 
                return;
            end
            
            list = d.dir; 
            
            for ii = 1:length(list)
                
                idx = regexp(list{ii}, '_run\d+$');
                if isempty(idx)
                    name = list{ii}; 
                else
                    name = list{ii}(1:idx-1); 
                end
                
                if ~isfield(s, name)
                    s.(name) = struct.empty;
                end
                
                new_struct = struct('name', name, 'start', '', 'end', '', 'runtime', [], 'num_files', []); 
                new_struct = obj.getObsLogFromFolder(fullfile(d.pwd, list{ii}), new_struct);
                s.(name) = vertcat(s.(name), new_struct); 
                
            end
            
        end
        
        function log_struct = getObsLogFromFolder(obj, folder, log_struct)
            
            import util.text.parse_value;
            
            if nargin<3 || isempty(log_struct)
                log_struct = struct;
                log_struct.name = '';
                log_struct.start = '';
                log_struct.end = '';
                log_struct.num_files = '';
                log_struct.runtime = [];
            end
            
            d = util.sys.WorkingDirectory(folder);
            
            files = d.match('*.h5*'); 
            readme = d.match('Z_README.txt');
            
            if ~isempty(readme)
                
                fid = fopen(readme{1}); 
                on_cleanup = onCleanup(@() fclose(fid)); 
                
                for ii = 1:1e4
                    
                    tline = fgetl(fid); 
                    
                    if ~ischar(tline)
                        break;
                    end
                    
                    [~, idx] = regexp(tline, 'RUNSTART:'); 
                    
                    if ~isempty(idx) && isempty(log_struct.start)
                        log_struct.start = strtrim(tline(idx+1:end)); 
                    end
                    
                    [~, idx] = regexp(tline, 'ENDTIME:');
                    
                    if ~isempty(idx) && isempty(log_struct.end)
                        log_struct.end = strtrim(tline(idx+1:end)); 
                    end
                    
                    [~, idx] = regexp(tline, 'END_STAMP:'); 
                    
                    if ~isempty(idx) && isempty(log_struct.runtime)
                        log_struct.runtime = parse_value(tline(idx+1:end)); 
                    end
                    
                end % for ii (file lines)
                
            elseif ~isempty(files) % can't find the readme, use the last HDF5 file instead
                
                try 
                    log_struct.start = h5readatt(files{end}, '/head', 'RUNSTART');
                    log_struct.end = h5readatt(files{end}, '/head', 'ENDTIME'); 
                    log_struct.runtime = h5readatt(files{end}, '/head', 'END_STAMP'); 
                catch 
                    log_struct.start = h5readatt(files{end}, '/header', 'RUNSTART');
                    log_struct.end = h5readatt(files{end}, '/header', 'ENDTIME'); 
                    log_struct.runtime = h5readatt(files{end}, '/header', 'END_STAMP'); 
                end
                
            end
            
            log_struct.num_files = length(files); 
            
        end
        
        function s = getDriveSpace(obj)
            
            s = struct;
            
            for ii = double('A'):double('Z')
                
                name = [char(ii) ':\'];
                if exist(name, 'dir')
                    s.(char(ii)) = util.sys.disk_space(name); 
                end
                
            end
            
        end
        
        function connectSync(obj)
            
            obj.log.input('Connecting to PcSync as server');
            
            try
                obj.sync.connect;
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function setBrake(obj)
            
            obj.brake_bit = 1;
            
        end
        
        function val = convertRuntimeToSeconds(obj)
            
            import util.text.cs;
            
            if cs(obj.runtime_units, 'seconds')
                val = 1;
            elseif cs(obj.runtime_units, 'minutes')
                val = 60;
            elseif cs(obj.runtime_units, 'hours')
                val = 3600;
            elseif cs(obj.runtime_units, 'batches')
                val = obj.batch_size./obj.getFrameRateEstimate;
            elseif cs(obj.runtime_units, 'frames')
                val = 1./obj.getFrameRateEstimate;
            else
                error('Unknown runtime units "%s". Use seconds, minutes, hours, batches or frames', obj.runtime_units);
            end
            
        end
        
        function val = is_slow_mode(obj)
            
            val = true;
            val = val && ~isempty(obj.expT) && obj.expT==obj.slow_mode_expT;
            val = val && ~isempty(obj.frame_rate) && obj.frame_rate==obj.slow_mode_frame_rate;
            val = val && ~isempty(obj.batch_size) && obj.batch_size==obj.slow_mode_batch_size;
            
        end
        
        function val = is_fast_mode(obj)
            
            val = true;
            val = val && ~isempty(obj.expT) && obj.expT==obj.fast_mode_expT;
            val = val && ~isempty(obj.frame_rate) && obj.frame_rate==obj.fast_mode_frame_rate;
            val = val && ~isempty(obj.batch_size) && obj.batch_size==obj.fast_mode_batch_size;
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = img.gui.AcqGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
        function show(obj, varargin)
            
            import util.text.cs;
            
            try
                
                input = util.text.InputVars;
                input.input_var('ax', [], 'axes', 'axis');
                input.scan_vars(varargin{:});

                try

                    if obj.use_verify_gui && obj.brake_bit==0
                        if isempty(obj.gui) || obj.gui.check==0
                            obj.makeGUI;
                        end
                    end

                catch ME
                    warning(ME.getReport);
                end

                if isempty(input.ax)
                    if ~isempty(obj.gui) && obj.gui.check
                        input.ax = obj.gui.axes_image;
                    else
                        input.ax = gca;
                    end
                end
                
                if cs(obj.show_what, 'images')
                    I = obj.images(:,:,end);
                elseif cs(obj.show_what, 'stack')
                    I = obj.stack;
                elseif cs(obj.show_what, 'stack_proc')
                    I = obj.stack_proc;
                else
                    error('Unknown option for "display_what". Use "images", "raw", or "stack"');
                end

                if obj.use_flip
                    I = rot90(I,2);
                end

                util.plot.setImage(I, input.ax);
                
                if obj.use_show_gray
%                     colormap(input.ax, 'gray'); % changing colormaps is very slow!
                    rect_color = 'white';
                else
%                     colormap(input.ax, 'default'); % changing colormaps is very slow!
                    rect_color = 'black'; 
                end
                
                obj.clip.showRectangles('num', obj.display_num_rect_stars, 'color', rect_color, 'ax', input.ax, 'flip', obj.use_flip, 'delete', 1, 'text', 0);
                obj.clip_bg.showRectangles('num', obj.display_num_rect_bg, 'color', 'red', 'ax', input.ax, 'flip', obj.use_flip, 'delete', 0, 'text', 0);
                
            catch ME
                warning(ME.getReport);
            end
            
        end
        
    end    
    
    methods(Access=protected) % bullshit functions to override setter/getter in AstroData (is this not cancelled??)
        
        function val = getPositions(obj)
            
            val = obj.clip.positions;
            
        end
        
        function setPositions(obj, val)
            
            obj.clip.positions = val;
            
        end
        
        function val = getPositionsBG(obj)
            
            val = obj.clip_bg.positions;
            
        end
        
        function setPositionsBG(obj, val)
            
            obj.clip_bg.positions = val;
            
        end
        
    end
    
end

