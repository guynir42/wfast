classdef CameraControl < file.AstroData

    properties(Transient=true)
        
        gui@obs.cam.gui.CamGUI;
        audio@util.sys.AudioControl;       
        
    end
    
    properties % objects and useful resources
        
        pars@head.Parameters; % parameters used in the observations
        cam; % this can be either a ZylaControl, a DhyanaControl or a SimCamera object
        buffers@file.BufferWheel;
        buf_live@file.BufferWheel;
        loupe@obs.cam.LoupeMagnifier;
        focuser; % connects to focus and has a focus screen (add constructor in getter)
%         af@obs.AutoFocus; 

        stats = struct('max', [], 'mean', [], 'std', [], 'cmean', [], 'cstd', [], 't1', [], 't2', []);
        
        v_curve_width;
        v_curve_position;
        
    end
     
    properties(Dependent=true) % outputs/switches that are stored in the specific camera
        
        hndl;
        cam_name;
        im_size;
        
        ROI; % gets the Region Of Interest position (left, top, width, height)
        
    end
    
    properties(Dependent=true) % outputs/switches that are stored in the current buffer
               
%         images;
%                 
%         timestamps;        
%         t_start;
%         t_end;
%         
%         psfs;
%         psf_sampling;
        
        latest_filename;
        
    end
    
    properties % inputs and switches
        
        mode = 'stars'; % can also choose "dark" or "flat"
        batch_size = 100; % how many frames in each batch
        num_batches = 2;
        
        use_flip_image = 0;
             
        height; % ROI height 
        width; % ROI width
        left; % ROI left corner
        top; % ROI top corner
        
        f = 100; % frame rate (Hz)
        T = 0.007; % exposure time (seconds)
        
        preview_time = 1; % seconds (used when calling preview action)
        num_frames_record = 2; % default number of frames to record
        mean_frame_rate;
        
        use_verify_frame_rate = 0;
        use_async = 1; % must remove this after changin the mex file for "capture"
        
        use_error_log = 1;
%         frame_delay_ms = 0;
        delay_batch_seconds = 0;
        
        debug_bit = 1;
        
    end
        
    properties(Hidden=true)
        
        start_time; % saved when batch started
        end_times; % end times of all batches
        num_frames_in_last_batch = 1;
        
        use_hardware_settings_warnings = 1; % this can be used to turn off warnings when the user settings are not correctly set in the hardware
        use_hardware_settings_errors = 1; % this can be used to turn off errors when the user settings are not correctly set in the hardware
        
        back_average = 100; % how many frames you want to average in "live view" to get the mean/max/std stats
        
        mex_flag = [0 0 0 0]; % for starting and stopping the camera...
        
        mode_list = {'stars', 'dark', 'flat'};
        mode_index = 1;
        
        default_f; % frame rate
        default_T; % exposure time
        default_batch_size; % number of frames per file
        default_num_batches; 
        default_camera = 'sim';
        
        current_ROI_size = 512; % this is for width and height
        current_ROI_size_list = [512, 1024, 256];
        current_ROI_size_index = 1;
        
        default_preview_time;
        default_num_frames_record;
        
        brake_bit = 1;
        now_recording = 0;
        now_live = 0;
        previous_folder = ''; % this is to send a warning when trying to write twice to the same folder
        use_prev_folder_warning = 1; % ask if the user is sure about recording twice to the same folder
        
        batch_counter = 0; % how many batches were recorded in this run
        clip_limit = 16000; % any pixels above this value will display a "clipping" warning... 
        is_clipping = 0; % update this manually in record/live view
        
        float_epsilon = 1e-3;
        error_log;
        
        version = 5.03;
        
    end
    
    methods % constructor
        
        function obj = CameraControl(type, varargin)
            
            if nargin<1 || isempty(type)
                type = '';
            end
            
            try % setup parameters
               
                obj.pars = head.Parameters;
                
            catch ME
                rethrow(ME);
            end
            
            try % setup buffers
                obj.setupBuffers;
            catch ME
                warning(ME.getReport);
            end
            
            try % setup camera
                obj.setupCamera(type, varargin{:});
            catch ME
                warning(ME.getReport)
            end
            
            try % setup all other objects
                obj.setupFocuser;
                obj.setupAuxObjects;
            catch ME
                warning(ME.getReport);
            end
            
            util.oop.save_defaults(obj);
            
            obj.initialize; % initialize all default values
            
        end
        
        function setupCamera(obj, type, varargin)
           
            import util.text.cs;
            
            if nargin<2 || isempty(type)
                type = questdlg('choose camera', 'choose camera', 'Simulator', 'Zyla', 'DSLR', 'Simulator');
            end
            
            if cs(type, {'simulator','simcamera'})
                obj.cam = obs.cam.SimCamera(varargin{:});
            elseif cs(type, 'zylacontrol')
                obj.cam = obs.cam.ZylaControl(varargin{:});
            elseif cs(type, 'dhyanacontrol')
                obj.cam = obs.cam.DhyanaControl(varargin{:});
            elseif cs(type, 'DSLR')
                obj.cam = obs.cam.DSLRController(varargin{:});
                obj.T = 1;
            else
                error(['unknown camera type: ' type ' use sim, zyla, dhyana or DSLR']);
            end
            
            obj.unzoom;
                                    
        end
                
        function setupBuffers(obj)

            obj.buffers = file.BufferWheel(5, obj.pars);
            obj.buffers.camera_mex_flag = obj.mex_flag;
            obj.buf_live = file.BufferWheel(2, obj.pars);
            obj.buf_live.camera_mex_flag = obj.mex_flag;
            
        end
        
        function setupFocuser(obj)
            
            try
                obj.focuser = obs.focus.FocusControl;
            catch ME
                disp('Cannot connect to focuser...');
                obj.focuser = obs.focus.FocusControl.empty;
            end
            
            if obj.gui.check
                obj.gui.makeFocusPanel;
            end
            
        end
        
        function setupAuxObjects(obj)
            
            obj.loupe = obs.cam.LoupeMagnifier(obj);
            obj.audio = util.sys.AudioControl;
            
        end
        
    end
    
    methods % reset methods
        
        function reconnect(obj)
            
            obj.cam.shutdown;
            obj.cam.initialize;
            
        end
        
        function initialize(obj) % this sets all the camera features to the default working settings
        
            util.oop.load_defaults(obj);
            
        end
        
        function shutdown(obj)
           
            obj.cam.shutdown;
            
        end
        
        function reset(obj) % get ready for a new run

            obj.buffers.reset; 
            util.vec.mex_change(obj.mex_flag, 4, 0); % reset the counter on the mex_flag
            
        end
        
        function resetStats(obj, Nbatches)
            
            obj.stats.max = zeros(Nbatches,1);
            obj.stats.mean = zeros(Nbatches,1);
            obj.stats.std = zeros(Nbatches,1);
            obj.stats.cmean = zeros(Nbatches,1);
            obj.stats.cstd = zeros(Nbatches,1);
            obj.stats.t_end = zeros(Nbatches,1);
            
        end
        
        function clear(obj) % removes outputs and intermediary data
            
            % maybe clear the buffers??
            
        end
        
        function clearStars(obj)
           
            obj.pars.clearStars;
            
        end
        
    end

    methods % getters
              
        function val = get.hndl(obj)
            
            if isempty(obj.cam)
                val = [];
            else
                val = obj.cam.hndl;
            end
            
        end
        
        function val = get.cam_name(obj)
            
            if isempty(obj.cam)
                val = '';
            else
                val = obj.cam.name;
            end
            
        end
        
        function val = get.im_size(obj)
           
            val = [obj.height, obj.width];
            
        end
        
        function val = get.ROI(obj)
            
            val = [obj.left, obj.top, obj.width, obj.height];
            
        end
        
        function val = getROItext(obj)
            
            val = util.text.print_vec(obj.im_size, 'x');
            
        end
        
        function val = getAverageFrameRate(obj)
            
            val = [];
            
        end
        
        function c = checkIsZoomed(obj)
           
            c = obj.left>1 || obj.top>1 || obj.width<obj.maxWidth || obj.height<obj.maxHeight;
            
        end
                        
        function val = is_capturing(obj)
            
            val = obj.mex_flag(1);
            
        end
        
        function val = is_finished(obj)
            
            val = obj.brake_bit;
            
        end
        
        function val = get.latest_filename(obj)
            
            val = obj.buffers.prev_buf.latest_filename;
            
        end
        
    end
    
    methods % setters
                
        function set.hndl(obj, val)
            
            if ~isempty(obj.cam)
                obj.cam.hndl = val;
            end
            
        end
        
        function set.pars(obj, val)
           
            obj.pars = val;
            
            if ~isempty(obj.buffers)
                obj.buffers.pars = val;
            end
            
            if ~isempty(obj.buf_live)
                obj.buf_live.pars = val;
            end
            
        end
        
        function set.T(obj, val)
            
            obj.T = val;
            
            if obj.setExpTimeHW(val) % if you succeed in setting this parameter, update the "pars" object
                obj.pars.T = val;
            end
            
        end
        
        function set.f(obj, val)
            
            obj.f = val;
            
            if obj.setFrameRateHW(val) % if you succeed in setting this parameter, update the "pars" object
                obj.pars.f = val;
            end
            
        end
        
        function set.height(obj, val)
            
            obj.height = val;
            
            if obj.setHeightHW(val) % if you succeed in setting this parameter, update the "pars" object
                obj.pars.AOI_height = val;
            end
            
        end
        
        function set.width(obj, val)
            
            obj.width = val;
            
            if obj.setWidthHW(val) % if you succeed in setting this parameter, update the "pars" object
                obj.pars.AOI_width = val;
            end
            
        end
        
        function set.top(obj, val)
            
            obj.top = val;
            
            if obj.setTopHW(val) % if you succeed in setting this parameter, update the "pars" object
                obj.pars.AOI_top = val;
            end
            
        end
        
        function set.left(obj, val)
            
            obj.left = val;
            
            if obj.setLeftHW(val) % if you succeed in setting this parameter, update the "pars" object
                obj.pars.AOI_left = val;
            end
            
        end
        
        function setCoordinates(obj, RA, DE)
           
            obj.pars.RA = RA;
            obj.pars.DE = DE;
            
        end
        
        function setTargetName(obj, name)
            
            obj.pars.target_name = name;
            
        end
        
        function addStar(obj,varargin)
           
            obj.pars.addStar(varargin{:});
            
        end
        
        function set.mex_flag(~, ~)
           
            disp('You must not change "mex_flag" directly. Use util.vec.mex_change or obj.stop or obj.unlockMexFlag');
            
        end
        
    end
    
    methods % activate the camera (API)
        
        function live(obj, varargin) % capture single-frame batches but do not save them
            
            input = obj.makeInputVars(varargin{:});
            input.num_batches = 1e6; % arbitrary timeout
            input.batch_size = 1;
            input.f = min(10, 0.95./input.T);
            
            obj.startup(input, obj.buf_live);
            
            if ~obj.gui.check
                obj.makeGUI;
            end
                        
            obj.now_recording = 0;
            obj.now_live = 1;
                                    
            finishup_cleanup = onCleanup(@obj.finishup); % make sure "finishup" is called when exiting this function

            for ii = 1:input.num_batches
                
                if obj.brake_bit || ~obj.gui.check 
                    disp('braking out of live view...');
                    obj.stop;
                    break;
                end
                
                ok = obj.batch(obj.buf_live);
                
                if ok==0, obj.stop; break; end
                
                obj.next(obj.buf_live);
                
            end
            
        end
        
        function record(obj, varargin) % record some batches
            
            input = obj.makeInputVars(varargin{:});
                        
            if obj.mex_flag(1)
                disp('Camera is already running...');
                return;
            end
            
            if ~obj.readyRecording(input)
                return; % if we fail some tests inside "readyRecording" 
            end
                        
            finishup_cleanup = onCleanup(@obj.finishup); % make sure "finishup" is called when exiting this function

            try
                                
                for ii = 1:input.num_batches
                    
                    if obj.brake_bit, break; end
                    
                    ok = obj.batch;                    
                    if ok==0, obj.stop; break; end
%                     pause(1.5); % simulate some pipeline processing time... (debugging only)
%                     disp(['index= ' num2str(obj.buffers.index)]);
                    obj.save;
                    obj.buffers.nextBuffer;
                    
                end
                
            catch ME
                obj.now_recording = 0;
                rethrow(ME);
            end
            
            obj.endRecording;
            
        end
        
        function autofocus(obj, varargin)
            
            % add a check that focuser is connected...
            
            input = obj.makeInputVars(varargin{:});
            input.input_var('range', 500);
            input.T = 1;
            input.num_batches = 100; % arbitrary timeout
            input.batch_size = 1;
            input.f = min(10, 0.95./input.T);
            input.scan_vars(varargin{:});
            
            obj.startup(input, obj.buf_live);
            
            if ~obj.gui.check
                obj.makeGUI;
            end
                        
            obj.now_recording = 0;
            obj.now_live = 1;
                        
            finishup_cleanup = onCleanup(@obj.finishup); % make sure "finishup" is called when exiting this function

            obj.v_curve_width = [];
            obj.v_curve_position = linspace(obj.focuser.pos-input.range, obj.focuser.pos+input.range, input.num_batches);
            obj.focuser.pos = obj.focuser.pos-input.range;
            
            for ii = 1:input.num_batches
                
                if obj.brake_bit 
                    disp('braking out of autofocus');
                    obj.stop;
                    break;
                end
                
                obj.focuser.pos = obj.v_curve_position(ii);
                
                ok = obj.batch(obj.buf_live);
                
                if ok==0, obj.stop; break; end
                
                I = obj.loupe.image_zoomed;
                I = util.img.centering(I);
%                 w = util.img.fwhm2(I);
                [~, w] = util.img.moments(double(I));
                w = mean(w(1:2));
                
                obj.v_curve_width(ii) = w;
                
                if obj.gui.check
                    delete(findobj(obj.gui.axes_image, 'Type', 'text'));
                    util.plot.inner_title(sprintf('pos= %d | width= %g', obj.focuser.pos, w), 'ax', obj.gui.axes_image, 'Position', 'South');
                end
                
                obj.next(obj.buf_live);
                
            end
            
            if length(obj.v_curve_position)>length(obj.v_curve_width)
                obj.v_curve_position(length(obj.v_curve_width):end) = [];
            end
            
        end
        
        function single(obj, varargin) % do a single exposure using "live"
            
            obj.live(varargin{:}, 'num_batches', 1);
            
        end
        
        function preview(obj, varargin) % do a single exposure using "live", but with expt=preview_time
            
            obj.live(varargin{:}, 'expT', obj.preview_time, 'num_batches', 1);
            
        end
        
        function recordLater(obj, delay_time, delay_units, varargin)
        % Setup a delayed recording
        % Usage: obj.recordLater(delay_time, delay_units='seconds', varargin)
        
            import util.text.cs;
            input = util.oop.full_copy(obj.rec_inputs);
            input.scan_vars(varargin{:});
            
            if nargin<2, help('obs.cam.CameraControl.recordLater'); return; end
            
            if nargin<3 || isempty(delay_units)
                delay_units = 'seconds';
            end
            
            if obj.debug_bit, disp(['setting up delayed recording of ' num2str(input.num_batches) ' batches... delayed by ' num2str(delay_time) ' ' delay_units]); end
            
            if cs(delay_units, 'seconds')
                dividor = 10^floor(log10(delay_time)-0.5);
                multiplier = 1;
            elseif cs(delay_units, 'minutes')
                dividor = 1;
                multiplier = 60;
            elseif cs(delay_units, 'hours')
                dividor = 3600;
                multiplier = 3600;
            else
                error(['unknown time units: ' delay_units]);
            end
            
            for t = delay_time:-1:1
                if mod(t-1, dividor)==0
                    disp(['t= ' num2str(t) ' ' delay_units]);
                end
                pause(multiplier);
            end
            
            obj.record(varargin{:});
            
        end
        
    end
    
    methods % access to hardware
        
        function val = floatMatch(obj, val1, val2, epsilon)
        % usage: val = obj.floatMatch(val1, val2, epsilon=obj.float_epsilon(1e-3))
        % compare two values and check if their relative size is smaller than
        % epsilon (relative size is abs(v1-v2)/(abs(v1)+abs(v2))).
        
            if nargin==0, help('obs.cam.CameraControl.floatMatch'); return; end
            
            if nargin<4 || isempty(epsilon)
                epsilon = obj.float_epsilon;
            end
            
            if isempty(val1) || isempty(val2)
                check_result = 0; % empty results are always a mis-match
            else
                check_result = abs(val1-val2)/(abs(val1)+abs(val2))<epsilon; % check if the results are close enough
            end
            
            if nargout>0
                val = check_result;
            else
                fprintf('VALUE MISMATCH: %s= %f | %s= %f\n', inputname(val1), val1, inputname(val2), val2);
            end            
            
        end
        
        function val = getTemperatureHW(obj)
            
            if obj.is_capturing
                val = [];
            else

                try
                    val = obj.cam.getTemperature;
                catch ME
                    val = [];
                    warning(ME.getReport);
                end

            end
            
        end
        
        function check = setTemperatureHW(obj, val)
            
            if obj.is_capturing
               check = 0; % if camera is running this check fails
            else
                check = 1; % assume the value is passed successfully, but change to 0 if there is an error or value mis-match
                
                try 
                    obj.cam.setTemperature(val);
                    hw_val = obj.cam.getTemperature; % check the HW actually got updated (directly to HW, so no checks or try/catch
                    if ~obj.floatMatch(val, hw_val)
                        error('User value: %f | HW value: %f',  val, hw_val); % this error (as well as error in set/get) is caught right after this...
                    end
                catch ME
                    check = 0;
                    if obj.use_hardware_settings_warnings
                        warning(ME.getReport);
                    end
                end
                
            end
            
        end
        
        function val = getExpTimeHW(obj)
            
            if obj.is_capturing
                val = [];
            else

                try
                    val = obj.cam.getExpTime;
                catch ME
                    val = [];
                    warning(ME.getReport);
                end

            end
        end
        
        function check = setExpTimeHW(obj, val)
            
            if nargin<2 || isempty(val)
                val = obj.T;
            end
            
            if obj.is_capturing
               check = 0; % if camera is running this check fails
            else
                check = 1; % assume the value is passed successfully, but change to 0 if there is an error or value mis-match
                
                try 
                    obj.cam.setExpTime(val);
                    hw_val = obj.cam.getExpTime; % check the HW actually got updated (directly to HW, so no checks or try/catch
                    if ~obj.floatMatch(val, hw_val)
                        error('User value: %f | HW value: %f',  val, hw_val); % this error (as well as error in set/get) is caught right after this...
                    end
                catch ME
                    check = 0;
                    if obj.use_hardware_settings_warnings
                        warning(ME.getReport);
                    end
                end
                
            end
            
        end
        
        function val = getExpTimeMaxHW(obj)
            
            
            if obj.is_capturing
                val = [];
            else

                try
                   val = obj.cam.getExpTimeMax;
                catch ME
                   val = [];
                   warning(ME.getReport);
                end
                
            end
            
        end
        
        function val = getExpTimeMinHW(obj)
            
            if obj.is_capturing
                val = [];
            else
                
                try
                    val = obj.cam.getExpTimeMin;
                catch ME
                    val = [];
                    warning(ME.getReport);
                end
                
            end
            
        end
                
        function val = getFrameRateHW(obj)
            
            if obj.is_capturing
                val = [];
            else

                try
                    val = obj.cam.getFrameRate;
                catch ME
                    val = [];
                    warning(ME.getReport);
                end

            end
        end
        
        function check = setFrameRateHW(obj, val)
            
            if nargin<2 || isempty(val)
                val = obj.f;
            end
            
            if obj.is_capturing
               check = 0; % if camera is running this check fails
            else
                check = 1; % assume the value is passed successfully, but change to 0 if there is an error or value mis-match
                
                try 
                    obj.cam.setFrameRate(val);
                    hw_val = obj.cam.getFrameRate; % check the HW actually got updated (directly to HW, so no checks or try/catch
                    if obj.use_verify_frame_rate
                        if ~obj.floatMatch(val, hw_val) % this may cause trouble... 
                            error('User value: %f | HW value: %f',  val, hw_val); % this error (as well as error in set/get) is caught right after this...
                        end
                    end
                catch ME
                    check = 0;
                    if obj.use_hardware_settings_warnings
                        warning(ME.getReport);
                    end
                end
                
            end
            
        end
        
        function val = getFrameRateMaxHW(obj)
           
            if obj.is_capturing
                val = [];
            else

                try
                    val = obj.cam.getFrameRateMax;
                catch ME
                    val = [];
                    warning(ME.getReport);
                end

            end
            
        end
        
        function val = getFrameRateMinHW(obj)
            
            if obj.is_capturing
                val = [];
            else

               try
                   val = obj.cam.getFrameRateMin;
               catch ME
                   val = [];
                   warning(ME.getReport);
               end
               
            end
            
        end
        
        function val = getHeightHW(obj)
            
            if obj.is_capturing
                val = [];
            else

                try
                    val = obj.cam.getHeight;
                catch ME
                    val = [];
                    warning(ME.getReport);
                end
            
            end
            
        end
        
        function check = setHeightHW(obj, val)
            
            if nargin<2 || isempty(val)
                val = obj.height;
            end
            
            if obj.is_capturing
               check = 0; % if camera is running this check fails
            else
                check = 1; % assume the value is passed successfully, but change to 0 if there is an error or value mis-match
                
                try 
                    obj.cam.setHeight(val);
                    hw_val = obj.cam.getHeight; % check the HW actually got updated (directly to HW, so no checks or try/catch
                    if val~=hw_val
                        error('User value: %d | HW value: %d',  val, hw_val); % this error (as well as error in set/get) is caught right after this...
                    end
                catch ME
                    check = 0;
                    if obj.use_hardware_settings_warnings
                        warning(ME.getReport);
                    end
                end
                
            end
            
        end
        
        function val = getWidthHW(obj)
               
            if obj.is_capturing
                val = [];
            else
         
                try
                    val = obj.cam.getWidth;
                catch ME
                    val = [];
                    warning(ME.getReport);
                end

            end
            
        end
        
        function check = setWidthHW(obj, val)
            
            if nargin<2 || isempty(val)
                val = obj.width;
            end
            
            if obj.is_capturing
               check = 0; % if camera is running this check fails
            else
                check = 1; % assume the value is passed successfully, but change to 0 if there is an error or value mis-match
                
                try 
                    obj.cam.setWidth(val);
                    hw_val = obj.cam.getWidth; % check the HW actually got updated (directly to HW, so no checks or try/catch
                    if val~=hw_val
                        error('User value: %d | HW value: %d',  val, hw_val); % this error (as well as error in set/get) is caught right after this...
                    end
                catch ME
                    check = 0;
                    if obj.use_hardware_settings_warnings
                        warning(ME.getReport);
                    end
                end
                
            end
            
        end
        
        function val = getTopHW(obj)
            
            if obj.is_capturing
                val = [];
            else

                try
                    val = obj.cam.getTop;
                catch ME
                    val = [];
                    warning(ME.getReport);
                end
                
            end
            
        end
        
        function check = setTopHW(obj, val)
            
            if nargin<2 || isempty(val)
                val = obj.top;
            end
            
            if obj.is_capturing
               check = 0; % if camera is running this check fails
            else
                check = 1; % assume the value is passed successfully, but change to 0 if there is an error or value mis-match
                
                try 
                    obj.cam.setTop(val);
                    hw_val = obj.cam.getTop; % check the HW actually got updated (directly to HW, so no checks or try/catch
                    if val~=hw_val
                        error('User value: %d | HW value: %d',  val, hw_val); % this error (as well as error in set/get) is caught right after this...
                    end
                catch ME
                    check = 0;
                    if obj.use_hardware_settings_warnings
                        warning(ME.getReport);
                    end
                end
                
            end
            
        end
        
        function val = getLeftHW(obj)
            
            if obj.is_capturing
                val = [];
            else

                try
                    val = obj.cam.getLeft;
                catch ME
                    val = [];
                    warning(ME.getReport);
                end
            end
            
        end
        
        function check = setLeftHW(obj, val)
            
            if nargin<2 || isempty(val)
                val = obj.left;
            end
            
            if obj.is_capturing
               check = 0; % if camera is running this check fails
            else
                check = 1; % assume the value is passed successfully, but change to 0 if there is an error or value mis-match
                
                try 
                    obj.cam.setLeft(val);
                    hw_val = obj.cam.getLeft; % check the HW actually got updated (directly to HW, so no checks or try/catch
                    if val~=hw_val
                        error('User value: %d | HW value: %d',  val, hw_val); % this error (as well as error in set/get) is caught right after this...
                    end
                catch ME
                    check = 0;
                    if obj.use_hardware_settings_warnings
                        warning(ME.getReport);
                    end
                end
                
            end
            
        end
        
        function val = maxHeight(obj)
            val = obj.cam.maxHeight;
        end
        
        function val = maxWidth(obj)
            val = obj.cam.maxWidth;
        end
        
    end
    
    methods % inner working 
        
        function input = makeInputVars(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('num_batches', obj.num_batches, 'Nbatches');
            input.input_var('batch_size', obj.batch_size);
            input.input_var('f', obj.f, 'frame rate', 'frequency');
            input.input_var('T', obj.T, 'exposure time', 'expT');
            input.input_var('height', obj.height);
            input.input_var('width', obj.width);
            input.input_var('top', obj.top);
            input.input_var('left', obj.left);
            input.scan_vars(varargin{:});
            
        end
        
        function cycleModes(obj)
            
            obj.mode_index = obj.mode_index + 1;
            if obj.mode_index>length(obj.mode_list)
                obj.mode_index = 1;
            end
            
            obj.mode = obj.mode_list{obj.mode_index};
            
        end
         
        function cycleROIsize(obj)
            
            obj.current_ROI_size_index = obj.current_ROI_size_index + 1;
            if obj.current_ROI_size_index>length(obj.current_ROI_size_list)
                obj.current_ROI_size_index = 1;
            end
            
            obj.current_ROI_size = obj.current_ROI_size_list(obj.current_ROI_size_index);
            
        end
        
        function zoom(obj, height, width, center_y, center_x)
            
            if nargin==1 && obj.checkIsZoomed % use the default call to "zoom" as "unzoom" as well
                obj.unzoom;
                return;
            end
            
            if nargin<2 || isempty(height)
                height = obj.current_ROI_size;
            end
            
            if nargin<3 || isempty(width)
                width = height;
            end
            
            if nargin<4 || isempty(center_y)
                center_y = floor((obj.maxHeight)/2)+1;
            end
            
            if nargin<5 || isempty(center_x)
                center_x = floor((obj.maxWidth)/2)+1;
            end
            
            % check the inputs make sense            
            if height>obj.maxHeight, error('Cannot set height to %d, maxHeight is %d', height, obj.maxHeight); end
            
            if width>obj.maxWidth, error('Cannot set width to %d, maxWidth is %d', width, obj.maxWidth); end                
            
            top = center_y - floor(height/2);            
            if top<1, error('Cannot set center_y to %d, the height (%d) is too big!', center_y, height); end
            
            left = center_x - floor(width/2);            
            if left<1, error('Cannot set center_x to %d, the width (%d) is too big!', center_x, width); end
            
            % if it is all ok, update the settings
            obj.height = height;            
            obj.width = width;
            obj.top = top;
            obj.left = left;
            
            if obj.debug_bit, disp(['Zoomed to ' util.text.print_vec(obj.ROI)]); end
            
        end
        
        function unzoom(obj)
            
            obj.height = obj.maxHeight;
            obj.width = obj.maxWidth;
            obj.top = 1;
            obj.left = 1;
            
            if obj.debug_bit, disp('Unzoomed'); end
            
        end
        
        function update(obj, input) % make sure the user inputs are in fact set in the camera... 
            
            import util.text.cs;
            
            % update hardware of any new settings
            obj.setExpTimeHW(input.T);
            obj.setFrameRateHW(input.f);
            obj.setHeightHW(input.height);
            obj.setWidthHW(input.width);
            obj.setTopHW(input.top);
            obj.setLeftHW(input.left);
            
            obj.update_pars(input);
            
            if cs(obj.mode, 'stars')
                if obj.checkIsZoomed
                    obj.buffers.product_type = 'RawAOI';
                    obj.cam.product_type = 'RawAOI';
                else
                    obj.buffers.product_type = 'RawFull';
                    obj.cam.product_type = 'RawFull';
                end
            elseif cs(obj.mode, 'dark')
                obj.buffers.product_type = 'dark';
                obj.cam.product_type = 'dark';
                obj.pars.type = 'dark';
                obj.pars.is_dark = 1;
            elseif cs(obj.mode, 'flat')
                obj.buffers.product_type = 'flat';
                obj.cam.product_type = 'flat';
                obj.pars.type = 'flat';
                obj.pars.is_flat = 1;
            else
                error('Unknown camera mode: %s. Use stars, dark or flat...', obj.mode);
            end
            
            if cs(obj.cam_name, 'sim camera')
                obj.buffers.product_type = [obj.buffers.product_type 'Sim'];
                obj.pars.is_sim = 1;
            end
            
        end
        
        function update_pars(obj, input) % verify HW settings are updated in parameter object
            
            import util.text.cs;
            
            obj.pars.datapath = obj.buffers.base_dir;
            
            obj.pars.T = obj.getExpTimeHW;
            obj.pars.f = obj.getFrameRateHW;
            
            obj.pars.AOI_height = obj.getHeightHW;
            obj.pars.AOI_top = obj.getTopHW;
            obj.pars.AOI_width = obj.getWidthHW;
            obj.pars.AOI_left = obj.getLeftHW;
            
            obj.pars.batch_size = input.batch_size;
            obj.pars.gain = obj.cam.gain;
            obj.pars.instrument = obj.cam.name;
            
            obj.pars.update;
            
            obj.pars.type = 'science';
            obj.pars.is_dark = cs(obj.mode, 'dark');
            obj.pars.is_flat = cs(obj.mode, 'flat');
            obj.pars.is_sim = isa(obj.cam, 'obs.cam.SimCamera');
            
        end
        
        function startup(obj, input, buf) % get ready to take images
            
            import util.text.cs;
            
            if obj.mex_flag(1)
                disp('Camera is already running. Stopping and restarting...');
                obj.stop;
                return;
            end
            
            if nargin<2 || isempty(input) || ~isa(input, 'util.text.InputVars')
                input1 = obj.makeInputVars;
                if isnumeric(input) && isscalar(input)
                    input1.num_batches = input;
                end
                input = input1;
            end
            
            if nargin<3 || isempty(buf)
                buf = obj.buffers;
            end
            
            obj.brake_bit = 0;
            
            obj.update(input);
            
            buf.reset;
            
            obj.resetStats(input.num_batches); % mean, max, std etc are calculated over a bunch of frames (not used right now)
            obj.end_times = zeros(input.num_batches,1); % preallocate
            
            obj.unlockMexFlag; % this should be replaced by a more carefull unlocking inside buffer>reset
            
            obj.batch_counter = 0;
            
%             if cs(obj.cam.name, 'Zyla')
%                 obj.captureMex(input.num_batches, input.batch_size, buf); % start rolling camera!
%             else
% %                 util.vec.mex_change(obj.mex_flag, 1, 1);
%             end
            
            obj.cam.startup(obj, input, buf);
               
            obj.start_time = tic;
                        
            obj.gui.update;
            
        end
        
        function finishup(obj) % after taking images is done
            
            import util.text.cs;
            
            if obj.debug_bit>5
                disp('finishup');
            end
            
            obj.stop;
            
            obj.now_recording = 0;
            obj.now_live = 0;
            
            obj.cam.finishup(obj);
            
            % this is currently unused (stats are too slow to calculate)
            f_names = fields(obj.stats);
            for jj = 1:length(f_names)
                obj.stats.(f_names{jj})(obj.batch_counter+1:end) = []; % clip the trailing zeros...
            end
            
            obj.end_times(obj.batch_counter+1:end) = []; % clip the trailing zeros...
            
            drawnow;
            
            obj.gui.update;
            
%             if cs(obj.cam.name, 'Zyla')
%                 
%             else
%                 util.vec.mex_change(obj.mex_flag, 1, 0);
%             end
            
        end
        
        function ok = readyRecording(obj, input) % only call this before recording, calls "startup" internally
            
            ok = 1;
            
            obj.update(input);
            
            if obj.use_prev_folder_warning && strcmp(obj.previous_folder, obj.buffers.directory)
                val = questdlg('Folder has already been recorded to... Continue anyway?', 'warning', 'Yes', 'No', 'Yes');
                if ~util.text.cs(val, 'Yes')
                    obj.brake_bit=1;
                    ok = 0;
                    return;
                end
            end
            
            obj.startup(input); % update parameters, start the camera... 
            
            filename = obj.buffers.getReadmeFilename;
            util.oop.save(obj, filename, 'name', 'camera'); 
            
            try obj.audio.playTakeForever; catch ME, warning(ME.getReport); end
            
            obj.now_recording = 1;
            
        end
        
        function endRecording(obj) % only call this after recording, calls "finishup" internally
                                    
            obj.now_recording = 0;
            
%             obj.finishup;
                        
            obj.previous_folder = obj.buffers.directory;
            
            try % sound effects
                
                obj.audio.playShowsOver;
                
            catch ME
                warning(ME.getReport)
            end
            
            filename = obj.buffers.getReadmeFilename('Z');
            
            util.oop.save(obj, filename, 'name', 'camera'); 
            
        end
        
%         function captureMex(obj, num_batches, batch_size, buffer_wheel)
%                         
%             if nargin<3 || isempty(batch_size)
%                 batch_size = obj.batch_size;
%             end
%             
%             if nargin<4 || isempty(buffer_wheel)
%                 buffer_wheel = obj.buffers;
%             end
%             
%             if ~isa(buffer_wheel, 'file.BufferWheel')
%                 error('Must give a file.BufferWheel object as input 4. Was given a %s...', class(buffer_wheel));
%             end
%             
%             if obj.use_error_log
%                 obj.error_log = zeros(num_batches,1); 
%             else
%                 obj.error_log = [];
%             end
%             
%             obs.cam.mex.capture(obj, obj.mex_flag, buffer_wheel.buf, obj.buffers.index_rec_vec, num_batches, batch_size, obj.error_log);
%             
%         end
        
        function ok = batch(obj, buf)
            
            ok = 1;
            
            if nargin<2 || isempty(buf) % use second argument input to read batches into buf_live or something...
                buf = obj.buffers;
            end
            
            if obj.brake_bit
                obj.stop;
                ok = 0;
                return;
            end
            
            obj.batch_counter = obj.batch_counter + 1; % displayed in the GUI
            
            if obj.delay_batch_seconds
                pause(obj.delay_batch_seconds);
            end
            
            obj.cam.batch(buf.this_buf, obj);
            
            buf.waitForRecording; % waits for the camera to fill the buffers
            buf.loadDataFromBuffer;
            
            if obj.debug_bit>5
                disp(['reading out batch ' num2str(obj.batch_counter) ' from buffer ' num2str(buf.index) ' | read_flag: ' util.text.print_vec(buf.this_buf.mex_flag_read)]);
            end
            
            obj.takeFrom(buf); % copies the pointers to the data in "buf"
            
%             I = buf.this_buf.images_raw;
%             if ~isempty(I), I(1) = I(1); end % do we need to make a copy of the data...?
%             t = buf.this_buf.timestamps; % decide later what to do with this...
%             if ~isempty(t), t(1) = t(1); end % do we need to make a copy of the data...?
            
            I = obj.images; % shorthand
            obj.loupe.input(I); % used for magnifying the image but also for autofocus
            
            % this code is way too slow for running full-frame
%             obj.setStats(I); 
%             [MX,MN,SD,CM,CS,FR] = obj.getStats;
%             obj.mean_frame_rate = FR.*size(I,3);
            
            buf.vec2times; % update the timing data from t_vec to t_start, t_end and t_end_stamp
            obj.pars.ephem.time = util.text.str2time(buf.t_start);
            obj.end_times(obj.batch_counter) = toc(obj.start_time);
            obj.mean_frame_rate = obj.getMeanFrameRate;
            obj.num_frames_in_last_batch = size(I,3);
            
            % check for clipping 
            MX = util.stat.max2(I(:,:,end)); % checking whole batch is too slow... need to get frame max from camera maybe?

            if MX>obj.clip_limit
                obj.is_clipping = 1;
            else
                obj.is_clipping = 0;
            end
            
            try % try to display the data on the GUI

                if obj.gui.check

                    h = findobj(obj.gui.axes_image, 'type','image');

                    if isempty(h) && obj.gui.check
                        disp('making new axes');
                        obj.gui.makeAxes;
                        h = findobj(obj.gui.axes_image, 'type','image');
                    end

                    if ~isempty(h)

                        if obj.use_flip_image
                            I = rot90(I(:,:,end), 2); % use this when viewing after the meridien flip
                        end

                        h.CData = I(:,:,end);

                    end
                    
%                     title(sprintf('max= %4.2f | mean= %4.2f | std= %4.2f | c.mean= %4.2f | c.std= %4.2f', MX, MN, SD, CM, CS), 'FontSize', 12, 'parent', obj.gui.axes_image);
                    title(sprintf('max= %4.2f | t-end(%d)= %4.1f | ROI: %s', MX, obj.batch_counter, obj.end_times(obj.batch_counter), util.text.print_vec(obj.ROI)), 'FontSize', 12, 'parent', obj.gui.axes_image);
                    
                    if obj.now_live && obj.batch_counter>1
                        obj.quickUpdateGUI;
                    else
                        obj.gui.update; % this is way too slow for live view!
                    end
                    
                end
                
            catch ME
                warning(ME.getReport);
            end
            
            drawnow;
            
        end
        
        function quickUpdateGUI(obj)
            
            obj.gui.panel_controls.button_frame_rate.String = ['f= ' num2str(obj.mean_frame_rate)];
            obj.gui.button_number.String = ['N= ' num2str(obj.batch_counter)];
            
            if obj.is_clipping
                obj.gui.button_clipping.control.BackgroundColor = 'yellow';
                obj.gui.button_clipping.String = 'Clipping!';
            else
                obj.gui.button_clipping.control.BackgroundColor = util.plot.GraphicButton.defaultColor;
                obj.gui.button_clipping.String = 'ready';
                if obj.now_live
                    obj.gui.button_clipping.String = 'live';
                elseif obj.now_recording
                    obj.gui.button_clipping.String = 'recording';
                end
            end
            
            if obj.brake_bit
                obj.gui.panel_run.button_run.String = 'RUN';
            else
                obj.gui.panel_run.button_run.String = 'STOP';
            end
            
        end
        
        function next(obj, buf)
           
            if nargin<2 || isempty(buf)
                buf = obj.buffers;
            end
            
            buf.nextBuffer;
                
        end
        
        function save(obj, buf)
           
            if nargin<2 || isempty(buf) % use this input to read batches into buf_live or something...
                buf = obj.buffers;
            end
            
            buf.gui.update;
            buf.save;
            
        end
        
        function stop(obj)
            
            util.vec.mex_change(obj.mex_flag, 2, 1);

            obj.brake_bit = 1;
            
        end
        
        function unlockMexFlag(obj)
            
            for ii = 1:length(obj.mex_flag)
                util.vec.mex_change(obj.mex_flag, ii, 0);
            end
            
        end
        
    end
    
    methods % statistics
        
        function setStats(obj, I) % save the averages of the current frame (too slow. not used right now)
           
            MX = util.stat.max2(I);
            MN = util.stat.mean2(I);
            SD = util.stat.std2(I);
            CM = util.stat.corner_mean(I);
            CS = util.stat.corner_std(I);
            T1 = toc(obj.start_time);
            
            obj.stats.max(obj.batch_counter) = mean(MX,3);
            obj.stats.mean(obj.batch_counter) = mean(MN,3);
            obj.stats.std(obj.batch_counter) = mean(SD,3);
            obj.stats.cmean(obj.batch_counter) = mean(CM,3);
            obj.stats.cstd(obj.batch_counter) = mean(CS,3);
            obj.stats.t_end(obj.batch_counter) = T1;
             
        end
        
        function [MX,MN,SD,CM,CS,FR] = getStats(obj) % calculate the averages of previous frames (too slow. not used right now)
            
            idx = max(1, obj.batch_counter-obj.back_average);
            
            MX = mean(obj.stats.max(idx:obj.batch_counter));
            MN = mean(obj.stats.mean(idx:obj.batch_counter));
            SD = mean(obj.stats.std(idx:obj.batch_counter));
            CM = mean(obj.stats.cmean(idx:obj.batch_counter));
            CS = mean(obj.stats.cstd(idx:obj.batch_counter));
            FR = (obj.batch_counter-idx+1)/(obj.stats.t_end(obj.batch_counter)-obj.stats.t_end(idx));
            
        end
        
        function val = getMeanFrameRate(obj) % we need to get this for monitoring...

            idx = max(1, obj.batch_counter-obj.back_average); % how many batches back we want to use to calculate this

            val = obj.num_frames_in_last_batch*(obj.batch_counter-idx)/(obj.end_times(obj.batch_counter)-obj.end_times(idx));
            
        end
        
    end
       
    methods % plot/GUI stuff
     
        function makeGUI(obj)
           
            if isempty(obj.gui)
                obj.gui = obs.cam.gui.CamGUI(obj);                
            end
            
            if ~isempty(obj.focuser) && isempty(obj.focuser.gui)
                obj.focuser.gui = obj.gui;
            end
            
            obj.gui.make;
            
        end
        
    end

end