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
        pars@head.Parameters;
        
        % input objects
%         cam@obs.cam.CameraControl;
        cam; % can be obs.cam.CameraControl or obs.cam.Andor
        reader@file.Reader;
        sim; % later add the class for the simulator        
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
        
        af@obs.focus.AutoFocus;
        
        % output to file
        buf@file.BufferWheel;
        
        deflator@file.Deflator;
        
        log@util.sys.Logger;
        
    end
    
    properties % inputs/outputs
        
        num_stars_found;
        cutouts_proc;
        cutouts_bg_proc;
        
        stack_cutouts; 
        stack_cutouts_bg;
        stack_proc;
        prev_stack;        
        ref_stack;
        ref_positions;
        
        prev_fluxes; % fluxes measured in previous batch (for triggering)
        
        batch_counter = 0;
        
    end
    
    properties % switches/controls
        
        start_index; % if non empty, use this number as initial index for run (e.g., to continue from where we stopped)
        
        use_background = 1;
        use_refine_bg = 0;
        
        % these swithces determine how stars are picked when run begins
        use_remove_saturated = false; % remove all stars with any pixels about saturation value
        saturation_value = 50000; % consider any pixels above this to be saturated
        
        use_quick_find_stars = true; % use new method that runs faster
        use_mextractor = false; % use mextractor to identify stars and find their WCS and catalog mag/temp
        use_arbitrary_pos = false; % don't look for stars (e.g., when testing with the dome closed)
        
        use_cutouts = true;
        use_adjust_cutouts = 1; % use adjustments in software (not by moving the mount)
        use_simple_photometry = 1; % use only sums on the cutouts instead of Photometry object for full cutouts
        
        use_save = false; % must change this when we are ready to really start
        use_triggered_save = false;
        
        % display parameters
        use_show = true;
        
        show_what = 'images'; % can choose "images" or "stack"
        num_display_stars = 30;
        num_display_bg = 30;
        
        use_flip = 0; % flip view by 180 degrees (for meridien flip)
        
        use_audio = 1;
        use_progress = 1;
        
        pass_source = {}; % parameters to pass to camera/reader/simulator
        pass_cal = {}; % parameters to pass to calibration object
        pass_back = {}; % parameters to pass to background object
        pass_phot = {}; % parameters to pass to photometry object
        pass_show = {}; % parameters to pass to show function
        
        debug_bit = 1;
        log_level = 1;
        
    end
    
    properties(Dependent=true)
        
        run_name; % the parameter "target_name" is used here to give a name to the whole run
        buf_full; % camera's buffers are used for full-frame dump on triggers
        
        frame_rate_average; % calculated from this object's timing data
        
        % these are read only camera info
        frame_rate_measured;
        sensor_temperature;
        
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
       
        brake_bit = 1; % when this is set to 1 (using the GUI, for example), the run stops. 
        is_running = 0; % when this is 1, cannot start a new run or anything
        is_running_single = 0; % when this is 1, cannot start a new run or anything
        
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
        
        start_index_;
        use_background_;
        use_refine_bg_;
        use_remove_saturated_;
        saturation_value_;
        use_mextractor_;
        use_arbitrary_pos_;
        use_cutouts_;
        use_adjust_cutouts_;
        use_save_;
        use_triggered_save_;
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
        
        version = 1.05;
        
    end
    
    methods % constructor
        
        function obj = Acquisition(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'Acquisition')
                if obj.debug_bit, fprintf('Acquisition copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                
                if obj.debug_bit, fprintf('Acquisition constructor v%4.2f\n', obj.version); end
                
                obj.log = util.sys.Logger('Acquisition', obj);
                obj.log.heartbeat(600);
                
                obj.reader = file.Reader;
                obj.sim; % fill this when we have a simulator
                
                obj.cal = img.Calibration;
                obj.back = img.Background;
                obj.clip = img.Clipper;
                obj.clip.use_adjust = 0; % adjust through photometry object?
                obj.clip_bg = img.Clipper;
                obj.clip_bg.use_adjust = 0;
                
                obj.phot_stack = img.Photometry;
                obj.phot_stack.use_aperture = 1;
                obj.phot_stack.use_gaussian = 0;
                obj.flux_buf = util.vec.CircularBuffer;
                
                obj.phot = img.Photometry;
                obj.phot.use_aperture = 0;
                obj.phot.use_gaussian = 0;
                obj.lightcurves = img.Lightcurves;
                
                obj.af = obs.focus.AutoFocus;
                
                obj.buf = file.BufferWheel;
                obj.buf.product_type = 'Cutouts';
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
                
                obj.setupDefaults;
                
                obj.pars = head.Parameters; % this also gives "pars" to all sub-objects
                
                util.oop.save_defaults(obj); % make sure each default_XXX property is updated with the current XXX property value. 
                
            end
            
        end

        function setupDefaults(obj)

            obj.num_batches = 2;
            obj.batch_size = 100;
            
            obj.num_stars = 250;
            obj.cut_size = 20;
            
            obj.num_backgrounds = 50;
            obj.cut_size_bg = 32;

        end

    end
    
    methods % reset/clear
        
        function reset(obj)
            
%             disp('resetting');
            
            if obj.brake_bit==0
                warning('Must stop acquisition before resetting!');
                return;
            end

            list = properties(obj);
            
            for ii = 1:length(list)
                
                if isobject(obj.(list{ii})) && ~isempty(obj.(list{ii})) && ismethod(obj.(list{ii}), 'reset') 
                    obj.(list{ii}).reset;
                end
                
            end

            obj.clear;
            
            obj.num_stars_found = [];
            obj.prev_fluxes = [];
            obj.batch_counter = 0;
            obj.start_index = 1;
            obj.positions = [];

        end
        
        function clear(obj)
            
            obj.prev_fluxes = obj.fluxes; % keep one batch history
            
            list = properties(obj);
            
            for ii = 1:length(list)
                
                if isobject(obj.(list{ii})) && ~isempty(obj.(list{ii})) && ismethod(obj.(list{ii}), 'clear') 
                    obj.(list{ii}).clear;
                end
                
            end
            
        end
        
    end
    
    methods % getters

        function val = get.run_name(obj)
            
            if ~isempty(obj.pars)
                val = obj.pars.target_name;
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
        
        function val = getTimeLeft(obj)
            
            if ~isempty(obj.frame_rate_average)
                val = (obj.num_batches-obj.batch_counter).*obj.batch_size./obj.frame_rate_average;
            elseif ~isempty(obj.frame_rate) && ~isnan(obj.frame_rate)
                val = (obj.num_batches-obj.batch_counter).*obj.batch_size./obj.frame_rate;
            else
                val = [];
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
        
        function val = get.expT(obj)
            
            if isprop(obj.src, 'expT')
                val = obj.src.expT;
            elseif ~isempty(obj.pars)
                val = obj.pars.expT;
            else
                val = [];
            end 
            
        end
        
        function val = get.frame_rate(obj)
            
            if isprop(obj.src, 'frame_rate')
                val = obj.src.frame_rate;
            elseif ~isempty(obj.pars)
                val = obj.pars.frame_rate;
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
        
        function val = get.frame_rate_measured(obj)
            
            if isa(obj.src, 'obs.cam.Andor')
                val = obj.src.frame_rate_measured;
            else
                val = [];
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
            
            val = obj.phot_stack.average_width;
            
        end
        
    end
    
    methods % setters
        
        function set.pars(obj,val)
            
            obj.pars = val;
            
            list = properties(obj);
            
            for ii = 1:length(list)
                
                if isobject(obj.(list{ii})) && ~isempty(obj.(list{ii})) && isprop(obj.(list{ii}), 'pars') 
                    obj.(list{ii}).pars = val;
                end
                
            end
            
        end
        
        function set.run_name(obj, val)
            
            if ~isempty(obj.pars)
                obj.pars.target_name = val;
            end
            
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
        
    end
    
    methods % utilities
        
        function input = makeInputVars(obj, varargin)
            
            idx = util.text.InputVars.isInputVars(varargin);
            
            if ~isempty(varargin) && any(idx) % was given an InputVars object
            
                idx = find(idx, 1, 'first');
                input = varargin{idx}; 
                varargin(find(idx, 1, 'first')) = []; % remove the one cell with InputVars in it
            
            else % use default values (load them from Acquisition object)
            
                input = util.text.InputVars;
                input.input_var('use_reset', 0, 'reset'); 
                input.input_var('start_index', []);
                input.input_var('use_background', []);
                input.input_var('use_refine_bg', []);
                input.input_var('use_remove_saturated', []);
                input.input_var('saturation_value', []);
                input.input_var('use_mextractor', []);
                input.input_var('use_arbitrary_pos', []);
                input.input_var('use_cutouts', []);
                input.input_var('use_adjust_cutouts', []);
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
                
                input.input_var('run_name', '', 'name', 'object', 'objname');
                input.input_var('RA', [], 'right ascention', 'right ascension');
                input.input_var('DE', [], 'declination');
                input.input_var('expT', [], 'T', 'exposure time');
                input.input_var('frame_rate', []); 
                input.input_var('num_batches', [], 'Nbatches');
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
            obj.use_save_ = obj.use_save;
            obj.use_triggered_save_ = obj.use_triggered_save;
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
            obj.use_save = obj.use_save_;
            obj.use_triggered_save = obj.use_triggered_save_;
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
                        obj.cam.pars = obj.pars;
                    end
                    
                    obj.src = obj.cam;
                    
                elseif cs(source, 'Dhyana')
                    error('This camera is no longer supported...');
                    if isempty(obj.cam) || ~isa(obj.cam.cam, 'obs.DhyanaControl')
                        obj.cam = obs.CameraControl('dhyana');
                        obj.cam.pars = obj.pars;
                    end
                    
                    obj.src = obj.cam;
                    
                elseif cs(source, 'simcamera')
                    error('This camera is no longer supported...');
                    if isempty(obj.cam) || ~isa(obj.cam.cam, 'obs.cam.SimCamera')
                        obj.cam = obs.cam.CameraControl('sim');
                        obj.cam.pars = obj.pars;
                    end
                    obj.src = obj.cam;
                    
                    
                elseif cs(source, 'simulator')
                    
                    if isempty(obj.sim)
                        obj.sim = img.Simulator;
                        obj.sim.pars = obj.pars;
                    end
                    
                    obj.src = obj.sim;

                elseif cs(source, {'file reader', 'reader'})
                    
                    if isempty(obj.reader)
                        obj.reader = file.Reader;
                        obj.cam.pars = obj.pars;
                    end
                    
                    obj.src = obj.reader;
                    
                elseif cs(source, 'camera')
                    if isempty(obj.cam)
%                         obj.cam = obs.cam.CameraControl;
                        obj.cam = obs.cam.Andor;
                        obj.cam.pars = obj.pars;
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
                    obj.src.pars = obj.pars;
                else
                    warning(['unknown source class: ' class(source)]);
                end
            end
           
        end
        
        function chooseDir(obj, dirname)
            
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
        
    end
    
    methods % commands/calculations
        
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

                    obj.prog.showif(obj.batch_counter);

                end

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function update(obj, input)
            
            if ~isempty(input.RA)
                obj.pars.RA = input.RA;
            end
            
            if ~isempty(input.DE)
                obj.pars.DE = input.DE;
            end
            
            if ~isempty(input.run_name)
                obj.pars.target_name = input.run_name;
            end
            
            obj.pars.update;
            
        end
        
        function check = startup(obj, varargin)
            
            check = 0;
            
            try 
                
                if obj.brake_bit==0
                    disp('Cannot start a new acquisition while old one is still runnning (turn off brake_bit)');
                    return;
                end

                if ~isempty(obj.gui) && obj.gui.check
                    obj.gui.update;
                end
                
                input = obj.makeInputVars(varargin{:});
            
                if input.log_level
                    obj.log.input('Starting a new run with Acquisition.'); % maybe add some more info here...?
                end
                
                if ~obj.cal.checkDark
                    error('Cannot start a new run without loading darks into calibration object!');
                end

                obj.stash_parameters(input);
                
                if obj.getTimeLeft>3600*10
                    error('Run scheduled to take %4.2f hours with these parameter... aborting!', obj.getTimeLeft/3600); 
                end
                
                if obj.getGbLeft>util.sys.disk_space(obj.buf.directory)
                    error('Run scheduled requires an estimated %5.2f Gb of storage. Only %5.2f Gb available on drive!', obj.getGbLeft, util.sys.disk_space(obj.buf.directory));
                end
                    
                if input.use_reset % this parameter is not saved in the object because we only use it here... 
                    
                    obj.reset;
                    
                    if obj.debug_bit, disp(['Starting run "' input.run_name '" for ' num2str(input.num_batches) ' batches.']); end

                end
                
                if isempty(obj.positions)
                    if obj.debug_bit, disp('Positions field empty. Calling single then findStars'); end
                    obj.single;
                    obj.findStars;
                end
                
                obj.update(input); % update pars object to current time and input run name, RA/DE if given to input.

                if obj.use_save
                    try
                        filename = obj.buf.getReadmeFilename;
                        util.oop.save(obj, filename, 'name', 'acquisition'); 
                    catch ME
                        warning(ME.getReport);
                    end
                end

                if obj.use_audio
                    try obj.audio.playTakeForever; catch ME, warning(ME.getReport); end
                end

                obj.src.num_batches = obj.num_batches;
                
                if isa(obj.src, 'file.Reader') % we need to change the reader interface to be more inline with the camera...
                    obj.src.num_frames_per_batch = obj.batch_size;
                    obj.src.num_files_per_batch = 1;
                    % what if batch_size is bigger than 100??
                end

                obj.src.startup('use_save', 0, obj.pass_source{:});

                if obj.use_progress
                    obj.prog.start(obj.num_batches);
                end
                
                obj.brake_bit = 0;
                
                check = 1;
                
            catch ME
                obj.log.error(ME.getReport);
                obj.unstash_parameters;
                obj.is_running = 0;
                rethrow(ME);
            end
            
        end
        
        function finishup(obj)
            
            if obj.use_progress
                obj.prog.finish;
            end
            
            obj.src.finishup;
            
            if obj.use_save
                try
                    filename = obj.buf.getReadmeFilename('Z');
                    util.oop.save(obj, filename, 'name', 'acquisition'); 
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
            
            t = tic;
            
            obj.prev_stack = obj.stack_proc; % keep one stack from last batch
            obj.clear;
            
            obj.src.batch; % produce the data (from camera, file, or simulator)
            
            if obj.debug_bit>1, fprintf('Starting batch %d. Loaded %d images from "%s" source.\n', obj.batch_counter+1, size(obj.images,3), class(obj.src)); end
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.update; 
            end
            
            drawnow; % make sure commands to the GUI and other callbacks are noticed... 
            
            obj.copyFrom(obj.src); % get the data into this object
            
            % if src is using ROI, must update the calibration object to do the same
            if isprop(obj.src, 'use_roi') && obj.src.use_roi 
                
                obj.cal.use_roi = 1;
                
                obj.cal.ROI = obj.src.ROI;
                
            else
                obj.cal.use_roi = 0;
            end
            
            obj.calcStack;
            
            if obj.use_cutouts
                obj.calcCutouts;
                obj.calcLightcurves;
                obj.calcTrigger;
            end
            
            obj.positions = obj.clip.positions;
            obj.positions_bg = obj.clip_bg.positions;
            
            if obj.use_save
                obj.buf.input(obj);
                obj.buf.clearImages; % do we need this if we have set use_save_raw_images=0 in the buffers?
                obj.buf.save;
                obj.buf.nextBuffer;
                % check triggering then call the camera save mode
                
            end
            
            if ismethod(obj.src, 'next')
                obj.src.next;
            end
            
            obj.batch_counter = obj.batch_counter + 1;
            
            if obj.use_progress
                obj.prog.showif(obj.batch_counter);
            end
            
            obj.start_index = obj.start_index + 1;
            
            if obj.use_show
                obj.show;
            end
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.update; 
            end
            
            drawnow;
            
            obj.runtime_buffer.input([size(obj.images,3), toc(t)]);
            
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
                
                obj.src.single; 
                obj.images = obj.src.images;
                obj.timestamps = obj.src.timestamps;
                obj.stack = sum(obj.src.images,3);
                obj.num_sum = size(obj.src.images,3);

                obj.stack_proc = obj.cal.input(obj.stack, 'sum', obj.num_sum); % stack after calibration

                check = 1;
                obj.is_running_single = 0;
                
            catch ME
                obj.log.error(ME.getReport);
                obj.is_running_single = 0;
                rethrow(ME);
            end
            
%             cleanup = onCleanup(@obj.finishup);
%             obj.startup('num_batches', 1, 'start_index', 1, 'use_save', 0, 'use_audio', 0, 'use_progress', 0);
%             
%             try 
%                 obj.batch;
%             catch ME
%                 obj.log.error(ME.getReport);
%                 rethrow(ME);
%             end
            
        end
        
        function calcStack(obj)
            
            % make the basic stack image
            obj.num_sum = size(obj.images,3);
            obj.stack = util.stat.sum_single(obj.images); % sum along the 3rd dimension directly into single precision
            
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
            
%             if obj.batch_counter==0
%                 obj.findStars; % this replaces the clipper findStars with somethinig better, or just runs it explicitely... 
%             end
            
            obj.stack_cutouts = obj.clip.input(obj.stack_proc); % if no positions are known, it will call "findStars"
            
            obj.phot_stack.input(obj.stack_cutouts, 'positions', obj.clip.positions); % run photometry on the stack to verify flux and adjust positions
            if ~isempty(obj.phot_stack.gui) && obj.phot_stack.gui.check, obj.phot_stack.gui.update; end
            
            obj.checkRealign;
            
            % store the latest fluxes from the stack cutouts
            obj.flux_buf.input(obj.phot_stack.fluxes);
            
            if obj.use_adjust_cutouts
                obj.clip.positions = double(obj.clip.positions + obj.average_offsets);
            else
                % must send the average adjustment back to mount controller (should we still adjust the cutouts though??)
            end
            
        end
        
        function findStars(obj)
            
%             disp('findStars');
            
            if obj.use_remove_saturated
                mu = median(squeeze(util.stat.corner_mean(util.img.jigsaw(obj.stack_proc))));
                sig = median(squeeze(util.stat.corner_std(util.img.jigsaw(obj.stack_proc))));
                obj.stack_proc = util.img.remove_saturated(obj.stack_proc, 'saturation', 4.5e4*obj.num_sum, 'threshold', mu+5*sig, 'dilate', 4); % note the saturation value is X100 because we are looking at the stack
            end
            
            if obj.use_arbitrary_pos
                obj.clip.arbitraryPositions; % maybe add some input parameters?
            elseif obj.use_mextractor
                % add the code for mextractor+astrometry here
            elseif obj.use_quick_find_stars
                T = util.img.quick_find_stars(obj.stack_proc, 'psf', 1.6, 'number', obj.num_stars,...
                    'dilate', obj.cut_size-5, 'saturation', obj.saturation_value.*obj.num_sum); 
                if isempty(T)
                    error('Could not find any stars using quick_find_stars!');
                end
                
                obj.clip.positions = T.pos;
                obj.positions = T.pos;
                obj.num_stars_found = size(obj.clip.positions,1);
                
            else
                obj.clip.findStars(obj.stack_proc);                
            end
            
            obj.ref_stack = obj.stack_proc;
            obj.ref_positions = obj.clip.positions;
            
        end
        
        function checkRealign(obj)
            
            if size(obj.flux_buf.data, 2)~=size(obj.phot_stack.fluxes,2)
                obj.flux_buf.reset;
            end
            
            if ~is_empty(obj.flux_buf) % check that stars are still aligned properly... 
                
                mean_fluxes = obj.flux_buf.mean;
                mean_fluxes(mean_fluxes<=0) = NaN;

                new_fluxes = obj.phot_stack.fluxes;
                new_fluxes(isnan(mean_fluxes)) = [];
                mean_fluxes(isnan(mean_fluxes)) = [];

                % maybe 0.5 is arbitrary and should be turned into parameters?
                if sum(new_fluxes<0.5*mean_fluxes)>0.5*numel(mean_fluxes) % lost half the flux in more than half the stars...

                    disp('Lost star positions, using quick_align');
                    
                    [~,shift] = util.img.quick_align(obj.stack_proc, obj.ref_stack);
                    obj.clip.positions = double(obj.ref_positions + flip(shift));
                    
                    % this shift should also be reported back to mount controller? 

                    obj.stack_cutouts = obj.clip.input(obj.stack_proc);

                    obj.phot_stack.input(obj.stack_cutouts, 'positions', obj.positions); % run photometry on the stack to verify flux and adjust positions
                    if ~isempty(obj.phot_stack.gui) && obj.phot_stack.gui.check, obj.phot_stack.gui.update; end
                    
                end

                % add second test and maybe quit the run if it fails...
                
            end
            
        end
        
        function calcCutouts(obj)
            
            obj.cutouts = obj.clip.input(obj.images);
            
            obj.cutouts_proc = obj.cal.input(obj.cutouts, 'clip', obj.clip);
            
            obj.cutouts_bg = obj.clip_bg.input(obj.images);
            
            obj.cutouts_bg_proc = obj.cal.input(obj.cutouts_bg, 'clip', obj.clip_bg);
            
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
        
        function calcLightcurves(obj)
            
            if obj.use_simple_photometry
                
                obj.fluxes = permute(util.stat.sum2(obj.cutouts_proc), [3,4,2,1]); 
                
%                 bad_pixels = util.stat.sum2(isnan(obj.cutouts_proc));
%                 obj.lightcurves.input(obj.fluxes, 'bad_pixels', bad_pixels);
                
            else % use extensive photometry method
            
                obj.phot.input('images', obj.cutouts_proc, 'timestamps', obj.timestamps, 'positions', obj.positions); % add variance input? 
                if ~isempty(obj.phot.gui) && obj.phot.gui.check, obj.phot.gui.update; end

                obj.fluxes = obj.phot.fluxes;

%                 obj.lightcurves.input(obj.phot.fluxes, obj.timestamps, obj.phot.weights, ...
%                     obj.phot.backgrounds, obj.phot.variances, ...
%                     obj.phot.offsets_x, obj.phot.offsets_y, obj.phot.centroids_x, obj.phot.centroids_y, ...
%                     obj.phot.widths, obj.phot.bad_pixels); % store the full lightcurves...
                
            end
            
            if obj.lightcurves.gui.check
                obj.lightcurves.gui.update;
            end
            
        end
        
        function calcTrigger(obj)
            
            % to be implemented!
            
        end
        
        function roi(obj, position) % this is cool but do we need this??
            
            import util.text.cs;
            
            if cs(position, 'none')
                obj.use_roi = 0;
            elseif cs(position, 'center')
                obj.use_roi = 1;
                obj.roi_x1 = 2160-256; % must change this to width/height of source
                obj.roi_x2 = 2160+255;
                obj.roi_y1 = 2560-256;
                obj.roi_y2 = 2560+255;
            elseif cs(position, 'northwest')
                obj.use_roi = 1;
                obj.roi_x1 = 1;
                obj.roi_x2 = 512;
                obj.roi_y1 = 1;
                obj.roi_y2 = 512;
            elseif cs(position, 'northeast')
                obj.use_roi = 1;
                obj.roi_x1 = 2160-512+1;
                obj.roi_x2 = 2160;
                obj.roi_y1 = 1;
                obj.roi_y2 = 512;
            elseif cs(position, 'southeast')
                obj.use_roi = 1;
                obj.roi_x1 = 2160-512+1;
                obj.roi_x2 = 2160;
                obj.roi_y1 = 2560-512+1;
                obj.roi_y2 = 2560;
            elseif cs(position, 'southwest')
               obj.use_roi = 1;
               obj.roi_x1 = 1;
                obj.roi_x2 = 512;
                obj.roi_y1 = 2560-512+1;
                obj.roi_y2 = 2560;
            else
                error('Unknown ROI position "%s". Use "center", "none", "NorthEast" etc...', position);
            end
            
        end
        
    end
    
    methods % optional run commands
        
        function runPreview(obj, varargin) % I'm not sure we need this... 
            
            if obj.is_running
                disp('Cannot run a preview during another run');
                return;
            else
                obj.is_running = 1;
            end
            
            try 
            obj.single;
            obj.findStars;
            obj.show;
            obj.is_running = 0;
             
            catch ME
                obj.is_running = 0;
                rethrow(ME);
            end
            
%             if ~isempty(varargin) && isa(varargin{1}, 'util.text.InputVars')
%                 input = varargin{1};
%                 input.scan_vars(varargin{2:end});
%             else
%                 input = obj.makeInputVars('num_batches', 1, 'batch_size', 1, 'expT', 1, 'num_stars', 0, ...
%                     'use_save', 0, 'use_audio', 0, varargin{:});
%             end
%             
%             obj.run(input); % take the same run-loop but with different input parameters
            
        end
        
        function runLive(obj, varargin) % I'm not sure we need this... 
            
            if obj.is_running
                return;
            else
                obj.is_running = 1;
            end
            
            if ~isempty(varargin) && isa(varargin{1}, 'util.text.InputVars')
                input = varargin{1};
                input.scan_vars(varargin{2:end});
            else
                input = obj.makeInputVars('num_batches', 1e6, 'batch_size', 1, 'num_stars', 0, ...
                    'use_save', 0, 'use_audio', 0, 'use_show', 1, varargin{:});
            end
            
            obj.run(input); % take the same run-loop but with different input parameters
            
        end
        
        function runFocus(obj, varargin)
            
            if obj.is_running || obj.is_running_single
                disp('Already running, set is_running and is_running_single to zero...');
                return;
            else
                obj.is_running = 1;
            end
            
            obj.log.input('Running autofocus');
            
            try
               
                if ~isempty(obj.gui) && obj.gui.check
                    obj.gui.update;
                end
                
                old_pos = [];
                
                if isempty(obj.cam) || isempty(obj.cam.focuser)
                    error('must be connected to camera and focuser!');
                end

                input = obj.makeInputVars('reset', 0, 'batch_size', 100, 'expT', 0.025, 'frame_rate', 25, ...
                    'num_stars', 25, 'use audio', 0, 'use_save', 0, 'run name', 'focus', 'prog', 0, ...
                    'pass_source', {'async', 0}, varargin{:});

                obj.reset;
                
                check = obj.single;
                if check==0, return; end
                
                obj.findStarsFocus;
                
                old_pos = obj.cam.focuser.pos; % keep this in case of error/failure
                
                obj.cam.focuser.pos = obj.cam.focuser.pos - obj.af.range; % starting point for scan... 
            
                obj.af.startup(old_pos, obj.positions, input.num_stars);
            
                pause(0.1);
            
                input.num_batches = length(obj.af.pos);
            
                obj.stash_parameters(input);
                obj.brake_bit = 0;
                cleanup = onCleanup(@obj.finishupFocus);
%                 obj.startup(input);

                p = [obj.af.pos flip(obj.af.pos)];

                for ii = 1:length(p)
                    if obj.brake_bit
                        return;
                    end

                    try 
                        
                        obj.cam.focuser.pos = p(ii);
                        
                        pause(0.1);
            
                    catch ME
                        warning(ME.getReport);
                    end
                    
                    check = obj.single;
                    if check==0, return; end
                    
                    obj.phot_stack.input(obj.clip.input(obj.stack_proc), 'positions', obj.clip.positions); 

                    obj.checkRealign;

                    if obj.use_adjust_cutouts
                        obj.clip.positions = double(obj.clip.positions + obj.average_offsets);
                    end
                    
                    if ~isempty(obj.gui) && obj.gui.check
                        obj.show;
                        obj.gui.update;
                    end
                    
                    obj.af.input(ii, obj.cam.focuser.pos, obj.phot_stack.widths, obj.phot_stack.fluxes);
                    
                    obj.af.plot;

                    drawnow;

                end
            
                obj.af.calculate;
                obj.brake_bit = 1;
                
                fprintf('FOCUSER RESULTS: pos= %f | tip= %f | tilt= %f\n', obj.af.found_pos, obj.af.found_tip, obj.af.found_tilt);
            
            catch ME
                
                disp(['Focus has failed, returning focuser to previous position= ' num2str(old_pos)]); 
                
                try 
                    obj.cam.focuser.pos = old_pos; 
                end
                
                obj.clip.reset; % don't save these star positions! 
                obj.positions = [];
                rethrow(ME);
                
            end
            
%             obj.cam.focuser.pos = obj.af.pos(1); % go back to lower position, then go back up (like the tank cannon)
            
            if ~isnan(obj.af.found_pos)
                obj.cam.focuser.pos = obj.af.found_pos;
            else
                disp('The location of new position is NaN. Choosing original position');
                obj.cam.focuser.pos = old_pos;
            end
            
            if isprop(obj.cam.focuser, 'tip') && ~isempty(obj.af.found_tip)
                obj.cam.focuser.tip = obj.cam.focuser.tip + obj.af.found_tip;
            end
            
            if isprop(obj.cam.focuser, 'tilt') && ~isempty(obj.af.found_tilt)
                obj.cam.focuser.tilt = obj.cam.focuser.tilt + obj.af.found_tilt;
            end
            
            obj.af.plot;
            
            obj.clip.reset; % don't save these star positions! 
            obj.positions = [];
            
        end
        
        function findStarsFocus(obj, varargin)
            
%             disp('findStarsFocus');
            
            I = obj.stack_proc;
            S = size(I);
            C = obj.cut_size;
            
            [M,V] = util.img.im_stats(I);
            
            T_all = table;
            
            markers = round(S'.*[1/3 2/3]); % divide the sensor to 1/3rds 
            
            unmask{1} = false(S);
            unmask{1}(markers(1,1):markers(1,2), markers(2,1):markers(2,2)) = 1; % only the central part is unmasked
            
            unmask{2} = false(S);
            unmask{2}(C:markers(1,1), C:markers(2,1)) = 1; % upper left corner
            
            unmask{3} = false(S);
            unmask{3}(markers(1,2):end-C+1, C:markers(2,1)) = 1; % lower left corner
            
            unmask{4} = false(S);
            unmask{4}(C:markers(1,1), markers(2,2):end-C+1) = 1; % upper right corner
            
            unmask{5} = false(S);
            unmask{5}(markers(1,2):end-C+1, markers(2,2):end-C+1) = 1; % lower right corner
            
            for ii = 1:5
            
                I_masked = I;
                I_masked(~unmask{ii}) = M;
                
                T = util.img.quick_find_stars(I_masked, 'mean', M, 'std', sqrt(V), 'number', 5, varargin{:});
            
                if ~isempty(T)
                    T_all = vertcat(T_all, T);
                end
                
            end
            
            obj.clip.positions = T_all.pos;
            obj.positions = double(T_all.pos);
            
            if obj.gui.check
                obj.show;
            end
            
        end
        
        function finishupFocus(obj)
            
            if obj.debug_bit, disp(['Finished run "' obj.run_name '" with ' num2str(obj.batch_counter) ' batches.']); end
            
            obj.brake_bit = 1;
            
            obj.unstash_parameters;
            
            obj.is_running = 0;
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.update;
            end
            
        end
        
        function findStarsFocusOld(obj) % find stars in order in five locations around the sensor
            
            I = obj.stack_proc;
            S = size(I);
            C = obj.cut_size;
            
            I = util.img.maskBadPixels(I);
            
            if obj.use_remove_saturated
                mu = median(squeeze(util.stat.corner_mean(util.img.jigsaw(I))));
                sig = median(squeeze(util.stat.corner_std(util.img.jigsaw(I))));
                I = util.img.remove_saturated(I, 'saturation', 4.5e4, 'threshold', mu+5*sig, 'dilate', 4);
            end
            
            I = conv2(I, util.img.gaussian2(2), 'same'); % smoothing filter
            
            markers = round(S'.*[1/3 2/3]); % divide the sensor to 1/3rds 
            
            mask{1} = false(S);
            mask{1}(markers(1,1):markers(1,2), markers(2,1):markers(2,2)) = 1; % only the central part is unmasked
            
            mask{2} = false(S);
            mask{2}(C:markers(1,1), C:markers(2,1)) = 1; % upper left corner
            
            mask{3} = false(S);
            mask{3}(markers(1,2):end-C+1, C:markers(2,1)) = 1; % lower left corner
            
            mask{4} = false(S);
            mask{4}(C:markers(1,1), markers(2,2):end-C+1) = 1; % upper right corner
            
            mask{5} = false(S);
            mask{5}(markers(1,2):end-C+1, markers(2,2):end-C+1) = 1; % lower right corner
            
            pos = zeros(obj.num_stars, 2);
            
            for ii = 1:obj.num_stars
                
                for jj = 1:100
                
                    [mx,idx] = util.stat.max2(I.*mask{mod(ii-1,5)+1}); % find the maximum in each masked area
                    
                    if any(idx-floor(C/2)<1), break; end
                    
                    I(idx(1)-floor(C/2):idx(1)+floor(C/2)+1, idx(2)-floor(C/2):idx(2)+floor(C/2)+1) = NaN; % remove found stars
                
                    if mx<obj.saturation_value*obj.num_sum % found a good star
                        pos(ii,:) = flip(idx); % x then y!
                        break; % pick this star and keep going
                    else % continue to look for stars in this quadrant
                        if obj.debug_bit, disp(['Star is saturated max= ' num2str(mx)]); end
                    end
                    
                end
                
            end
            
            obj.clip.positions = pos;
            obj.positions = double(pos);
            
            if obj.gui.check
                obj.show;
            end
            
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

                obj.clip.showRectangles('num', obj.num_display_stars, 'color', 'black', 'ax', input.ax, 'flip', obj.use_flip, 'delete', 1, 'text', 0);
                obj.clip_bg.showRectangles('num', obj.num_display_bg, 'color', 'red', 'ax', input.ax, 'flip', obj.use_flip, 'delete', 0, 'text', 0);

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

