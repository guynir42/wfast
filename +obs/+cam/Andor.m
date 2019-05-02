classdef Andor < file.AstroData

    properties(Transient=true)
        
        gui;
        audio@util.sys.AudioControl;
        prog@util.sys.ProgressBar;
        time_buffer@util.vec.CircularBuffer;
        latest_input@util.text.InputVars;
        
    end
    
    properties % objects
        
        pars@head.Parameters; % parameters used in the observations
        buffers@file.BufferWheel;
        focuser;
        log@util.sys.Logger;
        hndl; % pointer to the Andor SDK handle of the camera
        
    end
    
    properties % inputs/outputs
        
        status = 0;
        frame_rate_measured;
        frame_rate_camera;
        
    end
    
    properties % switches/controls
        
        mode = 'science'; % can also choose "dark" or "flat"
        batch_size = 100; % how many frames in each batch
        num_batches = 2;
        
        frame_rate = 28; % frame rate (Hz)
        expT = 0.025; % exposure time (seconds)
        expT_deep = 1; % for previews
        
        use_async = 1; % must remove this after changin the mex file for "capture"
        
        % Region Of Interest
        use_roi = 0; % when enabled the camera will use im_size and center_region
        
        use_show_flipped = 0; % flip the display for North is up after meridian flip
        
        debug_bit = 1;
        log_level = 1;
        
    end
    
    properties(Dependent=true)
        
        im_size; % size of the ROI (if scalar then width=height). If use_roi=0 will output maximum height/width
        center_region; % (y then x) of the center of the ROI. If empty, choose the center of the field
        ROI; % [top, left, height, width] as given to camera
        
    end
    
    properties(Hidden=true)
       
        mex_flag = [0 0 0 0]; % for starting and stopping the camera...
        
        max_height;
        max_width;
        
        mode_list = {'science', 'dark', 'flat'};
        mode_index = 1;
        
        default_frame_rate; 
        default_expT; 
        default_expT_deep;
        default_batch_size; 
        default_num_batches; 
        
        im_size_zoomed = 512; % default zoom size
        center_region_zoomed = []; % default zoom center (empty=centered)
        
        previous_folder = ''; % this is to send a warning when trying to write twice to the same folder
        use_prev_folder_warning = 1; % ask if the user is sure about recording twice to the same folder
        
        batch_counter = 0; % how many batches were recorded in this run
        clip_limit = 16000; % any pixels above this value will display a "clipping" warning... 
        is_clipping = 0; % update this manually in record/live view
        
        % these are hardware related parameters (used for synchronous record)
        imageSizeBytes; % need this to setup buffers for the camera
        AOIwidth_c;  % this helps understand how to parse the raw buffers
        AOIheight_c; % this helps understand how to parse the raw buffers
        AOIstride; % this helps understand how to parse the raw buffers
        clockFreq; % translate clock ticks to seconds
        
        num_restarts; % how many times did the camera get stuck (synchronous mode only!)
        
        brake_bit = 1;
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Andor(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.cam.Andor')
                if obj.debug_bit, fprintf('Andor camera copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Andor camera constructor v%4.2f\n', obj.version); end
                
                obj.log = util.sys.Logger('Andor_camera');
                
                util.oop.save_defaults(obj);
                
                obj.connect;
               
                obj.log.input('Loading auxiliary objects');
                
                try 
                    
                    obj.pars = head.Parameters;

                    obj.setupBuffers;
                    
                    obj.setupFocuser;
                    
                    obj.setupAudio;
                    
                    obj.setupProgressBar;
                    
                    obj.setupTimeBuffer;
                    
                    obj.update;

                catch ME
                    obj.log.error(ME.getReport);
                    rethrow(ME);
                end
                
            end
            
        end

        function connect(obj)
            
            obj.log.input('Connecting to Andor camera...');
            
            try 
                    
                rc = AT_InitialiseLibrary; AT_CheckError(rc);
                    
                [rc, obj.hndl] = AT_Open(0); AT_CheckError(rc);
%                 obj.hndl = obs.cam.mex_new.connect; % the new mex code is better than the matlab SDK because...?
                
                [rc] = AT_SetBool(obj.hndl, 'SensorCooling', 1); AT_CheckWarning(rc);        
              
                rc = AT_SetBool(obj.hndl, 'SpuriousNoiseFilter', 0); AT_CheckWarning(rc);
                rc = AT_SetBool(obj.hndl, 'StaticBlemishCorrection', 0); AT_CheckWarning(rc);
                
                [rc] = AT_SetEnumString(obj.hndl,'ElectronicShutteringMode','Rolling'); AT_CheckWarning(rc);
                               
                [rc] = AT_SetEnumString(obj.hndl,'SimplePreAmpGainControl','16-bit (low noise & high well capacity)'); AT_CheckWarning(rc);
                [rc] = AT_SetEnumString(obj.hndl,'PixelEncoding','Mono16'); AT_CheckWarning(rc);
                                
                % Enable Metadata
                [rc] = AT_SetBool(obj.hndl,'MetadataEnable',1); AT_CheckWarning(rc);
                [rc] = AT_SetBool(obj.hndl,'MetadataTimestamp',1); AT_CheckWarning(rc);
                
                [rc] = AT_SetEnumString(obj.hndl,'CycleMode','Continuous'); AT_CheckWarning(rc); % I think we don't need to limit the number of frames (there are explicit stops in the code).
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function disconnect(obj)
            
            obj.log.input('Disconnecting Andor camera...');
            
            try 
                obj.hndl = obs.cam.mex_new.disconnect;
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function reconnect(obj)
            
            obj.disconnect;
            pause(0.05);
            obj.connect;
            
        end
        
        function setupBuffers(obj)
            
            obj.buffers = file.BufferWheel(5, obj.pars);
            obj.buffers.camera_mex_flag = obj.mex_flag;

        end
        
        function setupFocuser(obj)
           
            try
                obj.focuser = obs.focus.FocusSpider;
            catch ME
                disp('Cannot connect to focuser, using simulator instead');
                obj.log.input('Cannot connect to focuser, using simulator instead');
                obj.focuser = obs.focus.Simulator;
            end
 
        end
        
        function setupAudio(obj)
            
            try
                obj.audio = util.sys.AudioControl;
            catch ME
                disp('Cannot connect to audio driver. Continuing without it!');
                obj.log.input('Cannot connect to audio driver. Continuing without it!');
            end
            
        end
        
        function setupProgressBar(obj)
            
            obj.prog = util.sys.ProgressBar;
            
        end
        
        function setupTimeBuffer(obj)
            
            obj.time_buffer = util.vec.CircularBuffer;
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj) % get ready for a new run

            obj.buffers.reset; 
            util.vec.mex_change(obj.mex_flag, 4, 0); % reset the counter on the mex_flag
            
            obj.time_buffer.reset;
            
            obj.frame_rate_camera = [];
            obj.frame_rate_measured = [];
            
            obj.batch_counter = 0;
            
        end
        
        function clear(obj) % removes outputs and intermediary data
            
            % maybe clear the buffers??
            
        end
        
    end
    
    methods % getters
        
        function val = is_finished(obj)
            
            val = obj.brake_bit;
            
        end
        
        function val = get.max_height(obj)
            
            if isempty(obj.max_height)
                obj.max_height = obj.maxHeightHW;
            end
            
            val = obj.max_height;
            
        end
        
        function val = get.max_width(obj)
            
            if isempty(obj.max_width)
                obj.max_width = obj.maxWidthHW;
            end
            
            val = obj.max_width;
            
        end
        
        function val = get.im_size(obj) % if empty, gives the max frame size
            
            if obj.use_roi==0
                val = [obj.max_height, obj.max_width];
            elseif isempty(obj.im_size_zoomed)
                val = [obj.max_height, obj.max_width];
            elseif isscalar(obj.im_size_zoomed)
                val = obj.im_size_zoomed.*[1 1];
            else
                val = obj.im_size_zoomed;
            end
            
        end
        
        function val = get.center_region(obj)
            
            if obj.use_roi==0
                val = []; % empty means centered! 
            else
                val = obj.center_region_zoomed;
            end
            
        end
        
        function val = get.ROI(obj)
            
            if obj.use_roi==0
                val = [1 1 obj.max_height obj.max_width];
            else
                val = obj.zoom2roi(obj.im_size, obj.center_region);
            end
            
        end
        
        function val = zoom2roi(obj, im_size, center)
            
            if nargin<2
                im_size = obj.im_size_zoomed;
            end
            
            if isscalar(im_size)
                im_size = im_size.*[1 1];
            end
            
            if im_size(1)>obj.max_height
                im_size(1) = obj.max_height;
            end
            
            if im_size(2)>obj.max_width
                im_size(2) = obj.max_width;
            end
            
            if nargin<3
                center = obj.center_region;
            end
            
            if isscalar(center)
                center = center.*[1 1];
            end
            
            if isempty(center)
                center = floor([obj.max_height, obj.max_width]/2)+1;
            end
            
            h = im_size(1);
            w = im_size(2);
            t = center(1) - ceil((h-1)/2);
            b = center(1) + ceil(h/2) - 1;
            l = center(2) - ceil((w-1)/2);
            r = center(2) + ceil(2/2) - 1;
            
            if t<1, h = h + t - 1; t = 1; end
            if b>obj.max_height, h = h - b + obj.max_height; b = obj.max_height; end
            if l<1, w = w + l - 1; l = 1; end
            if r>obj.max_width, h = h - r + obj.max_width; r = obj.max_height; end
            
            val = [t,l,h,w];
            
        end
        
        function [im_size, center] = roi2zoom(obj, roi)
            
            im_size(1) = roi(3);
            im_size(2) = roi(4);
            center(1) = roi(1) + floor(roi(3)/2);
            center(2) = roi(2) + floor(roi(4)/2);
            
        end
        
        function val = getZoomStr(obj, im_size, center) % a string describing the zoom state as given by camera defaults or by the two inputs
            
            if nargin<2
                im_size = obj.im_size;
            end
            
            if isscalar(im_size)
                im_size = im_size.*[1 1];
            end
            
            if nargin<3
                center = obj.center_region;
            end
            
            if isscalar(center)
                center = center.*[1 1];
            end
            
            % verify that the center+size is not outside the bounds, which makes the ROI smaller than you think...
            roi = obj.zoom2roi(im_size, center);
            [im_size, center] = obj.roi2zoom(roi); 
            
            if isempty(im_size) || (roi(3)==obj.max_height && roi(4)==obj.max_width)
                val = 'full frame';
                return;
            else
                val = util.text.print_vec(im_size, 'x');
            end
            
            if ~isempty(center) && (center(1)~=floor(obj.max_height/2)+1 || center(2)~=floor(obj.max_width/2)+1)
                val = sprintf('%s at %s', val, util.text.print_vec(center, ','));
            end
            
        end
        
    end
    
    methods % setters
        
        function set.pars(obj, val)
           
            obj.pars = val;
            
            if ~isempty(obj.buffers)
                obj.buffers.pars = val;
            end
        
        end
        
        function set.mex_flag(~, ~)
           
            disp('You must not change "mex_flag" directly. Use util.vec.mex_change or obj.stop or obj.unlockMexFlag');
            
        end
        
        function set.im_size(obj, val)
            
            if isempty(val) || (val(1)==obj.max_height && val(2)==obj.max_width) % input size is used to unzoom
                obj.use_roi = 0;
            else % change the zoomed size
                obj.im_size_zoomed = val;
                obj.use_roi = 1;
            end
            
        end
        
        function set.center_region(obj, val)
            
            obj.center_region_zoomed = val; 
            
        end
        
        function set.center_region_zoomed(obj, val)
            
            obj.center_region_zoomed = round(val); 
            
        end
        
        function zoom(obj, varargin)
            
            obj.use_roi = 0;
            
            if isempty(varargin)
                
            elseif isscalar(varargin)
                if isnumeric(varargin{1})
                    if length(varargin{1})<=2
                        obj.im_size = varargin{1};
                    elseif length(varargin{1})==3
                        obj.im_size = varargin{1}(1:2);
                        obj.center_region = varargin{1}(3).*[1 1];
                    elseif length(varargin{1})==4
                        obj.im_size = varargin{1}(1:2);
                        obj.center_region = varargin{1}(3:4);
                    else
                        error('must input a vector with 4 or less elements! size(varargin{1})= %s', util.text.print_vec(size(varargin{1})));
                    end
                        
                elseif ischar(varargin{1})
                    if util.text.cs(varargin{1}, 'unzoom')
                        obj.use_roi = 0;
                    elseif util.text.cs(varargin{1}, 'zoom')
                        obj.use_roi = 1;
                    else
                        error('Unknown zoom command "%s". Use "zoom" or "unzoom" instead...', varargin{1});
                    end
                end
            elseif length(varargin)==2
                obj.im_size(1) = varargin{1};
                obj.im_size(2) = varargin{2};
            elseif length(varargin)==3
                obj.im_size(1) = varargin{1};
                obj.im_size(2) = varargin{2};
                obj.center_region(1) = varargin{3};
                obj.center_region(2) = varargin{3};
            elseif length(varargin)==4
                obj.im_size(1) = varargin{1};
                obj.im_size(2) = varargin{2};
                obj.center_region(1) = varargin{3};
                obj.center_region(2) = varargin{4};
            end
            
        end
        
        function unzoom(obj)
            
            obj.use_roi = 0;
            
        end
        
    end
    
    methods % high level run commands
        
        function record(obj, varargin)
            
            obj.run('save', 1, 'audio', 1, 'progress', 1, 'log_level', 1, 'async', 1, 'show', 1, varargin{:});
            
        end
        
        function live(obj, varargin)
            
            obj.run('num_batches', 1e6, 'batch_size', 1, 'show', 1, 'save', 0,...
                'audio', 0, 'async', 0, 'progress', 0, 'log_level', 0, ...
                'frame rate', [], varargin{:}); 
            
        end
        
        function preview(obj, varargin)
            
            obj.run('num_batches', 1, 'batch_size', 1, 'expT', obj.expT_deep, 'frame_rate', [], ...
                'show', 1, 'save', 0, 'audio', 0, 'async', 0,...
                'progress', 0, 'log_level', 0, varargin{:}); 
            
        end
        
        function autofocus(obj, varargin)
            
            error('This is not yet implemented!');
            
        end
        
    end
    
    methods % interface commands to camera
        
        function update(obj) % check that connection is still ok
            
            disp('updating camera');
            
            obj.status = 1;
            
            % add some tests to see if camera is alive
            rc = AT_GetFloat(obj.hndl, 'ExposureTime');
            
            if rc
                obj.status = 0; 
            end
            
        end
        
        function run(obj, varargin)
            
            try
                input = obj.makeInputVars(varargin{:});
                obj.latest_input = input;
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
            if input.log_level
                obj.log.input(sprintf('Starting a new run with %d batches of %d images each. expT= %f, async= %d, save= %d, ROI= %s',...
                    input.num_batches, input.batch_size, input.expT, input.async, input.save, obj.getZoomStr(input.im_size, input.center))); 
            end
            
            try
                
                obj.startup(input);
                on_cleanup = onCleanup(@() obj.finishup(input));
                
                for ii = 1:input.num_batches
                    
                    if obj.brake_bit, break; end
                    
                    obj.batch(input);
                    
                end
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function startup(obj, varargin)
            
            input = obj.makeInputVars(varargin{:});
            
            obj.reset;
            
            % make sure hardware is updated... 
            obj.setROI_HW(obj.zoom2roi(input.im_size, input.center));
            obj.setExpTimeHW(input.expT);
            obj.setFrameRateHW(input.frame_rate);
            
            if isempty(input.frame_rate) % in this mode the camera takes an image as soon as it gets a command to "software trigger"
%                 [rc] = AT_SetEnumString(obj.hndl,'CycleMode','Continuous'); AT_CheckWarning(rc); 
                [rc] = AT_SetEnumString(obj.hndl,'TriggerMode','Software'); AT_CheckWarning(rc);
            else % in this mode there is a fixed frame rate, so the frame rate may be lower than the maximum 
%                 [rc] = AT_SetEnumString(obj.hndl,'CycleMode','Fixed'); AT_CheckWarning(rc);
                [rc] = AT_SetEnumString(obj.hndl,'TriggerMode','Internal'); AT_CheckWarning(rc);
            end
            
            obj.check_inputs(input); % check the hardware/input configuration is compatible
            obj.update_pars(input); % update the Parameter object
            obj.update_buffers(input); % upodate the Buffer object
            
            obj.num_restarts = 0; % how many times did the camera have to be restarted (synchronous mode only!)
            
            obj.batch_counter = 0; % how many batches were already taken (for GUI display, etc)
            
            if ~isempty(obj.gui)
                obj.gui.update; % update before turning off brake bit (so hardware info gets updated before the run)
            end
            
            obj.brake_bit = 0; % this allows all the loops to continue. It becomes 1 if the GUI stop button is pressed.
            
            if input.save
                filename = obj.buffers.getReadmeFilename;
                util.oop.save(obj, filename, 'name', 'camera');  
            end
            
            if input.audio
                try 
                    obj.audio.playTakeForever;
                catch ME
                    warning(ME.getReport);
                end
            end
            
            if input.progress
                obj.prog.start(input.num_batches);
            end
            
            if input.async
                
                if obj.mex_flag(1)
                    obj.stop;
                    error('Camera is already running. Stopping and restarting...');
                end
                
                obj.unlockMexFlag; % this should be replaced by a more carefull unlocking inside buffer>reset

                for ii = 1:length(obj.buffers.buf)
                    util.vec.mex_change(obj.buffers.buf(ii).mex_flag_record, 1, 1); % lock all the buffers for recording, they will unlock once the camera fills them... 
                end
                
                obs.cam.mex_new.startup(obj, obj.mex_flag, obj.buffers.buf, obj.buffers.index_rec_vec, input.num_batches, input.batch_size); % call the mex file for async recording
                
            else
                
                if strcmp(obj.getCycleModeHW, 'Fixed')
                    [rc] = AT_SetInt(obj.hndl,'FrameCount',input.batch_size); AT_CheckWarning(rc);
                end
                
                [rc, obj.imageSizeBytes] = AT_GetInt(obj.hndl, 'ImageSizeBytes'); AT_CheckWarning(rc)

                % note width and height are flipped due to C -> matlab conventions
                [rc, obj.AOIwidth_c] = AT_GetInt(obj.hndl, 'AOIWidth'); AT_CheckWarning(rc)
                [rc, obj.AOIheight_c] = AT_GetInt(obj.hndl, 'AOIHeight'); AT_CheckWarning(rc)

                [rc, obj.AOIstride] = AT_GetInt(obj.hndl, 'AOIStride'); AT_CheckWarning(rc)

                [rc, obj.clockFreq] = AT_GetInt(obj.hndl, 'TimestampClockFrequency'); AT_CheckWarning(rc)
                
                [rc] = AT_Flush(obj.hndl); AT_CheckWarning(rc);
                
                for X = 1:10 % setup buffers
                    [rc] = AT_QueueBuffer(obj.hndl, obj.imageSizeBytes); AT_CheckWarning(rc);
                end

                [rc] = AT_Command(obj.hndl, 'TimestampClockReset'); AT_CheckWarning(rc); % reset the time stamp clock
                
                [rc] = AT_Command(obj.hndl, 'AcquisitionStart'); AT_CheckWarning(rc);

            end
            
            if ~isempty(obj.gui)
                obj.gui.update; % update after everything is set up
            end
            
        end
        
        function finishup(obj, varargin)
            
            input = obj.makeInputVars(varargin{:});
            
            if input.async
                
                obj.stop;
                
            else
                
                [rc] = AT_Command(obj.hndl,'AcquisitionStop'); AT_CheckWarning(rc);
            
                [rc] = AT_Flush(obj.hndl); AT_CheckWarning(rc);
            
            end
            
            if input.save
                filename = obj.buffers.getReadmeFilename('Z');
                util.oop.save(obj, filename, 'name', 'camera');  
            end
            
            if input.audio
                try
                    obj.audio.playShowsOver;
                catch ME
                    warning(ME.getReport);
                end
            end
            
            if input.progress
                obj.prog.finish(obj.batch_counter);
            end
            
            obj.brake_bit = 1;
            
        end
        
        function batch(obj, varargin)
            
            input = obj.makeInputVars(varargin{:});
            
            if input.async
                obj.batch_async(input)
            else
                obj.batch_sync(input)
            end
            
            obj.batch_counter = obj.batch_counter + 1;
            
            if input.show
                obj.show(input.pass_show{:});
            end

            obj.frame_rate_camera = size(obj.images,3)./(obj.t_end_stamp - obj.timestamps(1)); % how many frames in what time
            
            if ~isempty(obj.buffers.t_start)
                
                time = util.text.str2time(obj.buffers.t_start);
                dt = seconds(time - obj.pars.ephem.time); % how much time passed since last update
                obj.pars.ephem.time = time;
                
                obj.time_buffer.input([size(obj.images,3), dt]); % how many images, how much time passed (total)
                
                obj.frame_rate_measured = sum(obj.time_buffer.data(:,1))./sum(obj.time_buffer.data(:,2)); 
                
            end
            
            if input.save
                obj.buffers.save;
            end
            
            obj.buffers.nextBuffer;
            
            if ~isempty(obj.gui)
                obj.gui.update;
            end
            
            if input.progress
                obj.prog.showif(obj.batch_counter)
            end
            
            drawnow;
            
        end
        
    end
    
    methods % internal functions 
        
        function input = makeInputVars(obj, varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'util.text.InputVars')
                input = varargin{1};
                varargin = varargin(2:end);
            else
                input = util.text.InputVars;
                input.input_var('mode', obj.mode); % science, dark, flat
                input.input_var('num_batches', obj.num_batches, 'Nbatches');
                input.input_var('batch_size', obj.batch_size);
                input.input_var('frame_rate', obj.frame_rate, 'frequency');
                input.input_var('expT', obj.expT, 'exposure time', 'expT'); 
                input.input_var('im_size', obj.im_size); % set the height and width of the ROI
                input.input_var('center', obj.center_region, 'center_region', 'center_ROI'); % set the center of the ROI (set to center of field if empty)
                input.input_var('async', true, 'use async', 5); % use mex-async code or simple batch-by-batch interface
                input.input_var('show', true, 'use_show', 5); % if you want to show each batch's images/stack
                input.input_var('save', true, 'use_save', 5); % if you want to save each batch 
                input.input_var('audio', false, 'use_audio', 5); % turn on/off audio signals
                input.input_var('progress', true, 'use_progress', 5); % display a progress bar on screen
                input.input_var('log_level', 1); % choose if and how much logging you want for this run (1 is only start of run). Errors are always logged. 
                input.input_var('pass_show', {}, 6); % parameters to pass to "show" function
            end
            
            input.scan_vars(varargin{:});
            
        end
        
        function cycleModes(obj)
            
            obj.mode_index = obj.mode_index + 1;
            if obj.mode_index>length(obj.mode_list)
                obj.mode_index = 1;
            end
            
            obj.mode = obj.mode_list{obj.mode_index};
            
        end
        
        function check_inputs(obj, input)
            
            if util.text.cs(input.mode, 'dark', 'flat')
            
                if obj.is_zoomed_HW
                    error('We should never take dark/flat images in ROI!');
                end
                
            end
            
            [rc, cooling] = AT_GetBool(obj.hndl, 'SensorCooling'); AT_CheckWarning(rc);
            if ~cooling, error('Camera sensor cooling is off!'); end
            
            % add other checks for consistency in parameters... 
            
        end
        
        function update_pars(obj, input)
            
            obj.pars.datapath = obj.buffers.base_dir;
            
            obj.pars.expT = obj.getExpTimeHW;
            obj.pars.frame_rate = obj.getFrameRateHW;
            
            if ~isempty(obj.focuser)
                obj.pars.focus = obj.focuser.pos;
            end
            
            roi = obj.getROI_HW;
            obj.pars.AOI_top = roi(1);
            obj.pars.AOI_left = roi(2);
            obj.pars.AOI_height = roi(3);
            obj.pars.AOI_width = roi(4);
            
            obj.pars.im_size = [roi(3) roi(4)];
            
            obj.pars.batch_size = input.batch_size;
            obj.pars.gain = obj.getGain;
            obj.pars.instrument = obj.getCameraNameHW;
            
            obj.pars.update; % get the current time etc
            
            obj.pars.type = input.mode;
            obj.pars.is_dark = util.text.cs(input.mode, 'dark');
            obj.pars.is_flat = util.text.cs(input.mode, 'flat');
            
        end
        
        function update_buffers(obj, input)
            
            import util.text.cs;
            
            if cs(input.mode, 'science')
                
            elseif cs(input.mode, 'dark')
                obj.buffers.product_type = 'Dark';
            elseif cs(input.mode, 'flat')
                obj.buffers.product_type = 'Flat';
            end
            
        end
        
        function batch_async(obj, input)
            
            obj.buffers.waitForRecording; % waits for the camera to fill the buffers
            obj.buffers.vec2times; % update the timing data from t_vec to t_start, t_end and t_end_stamp
            obj.buffers.loadDataFromBuffer;
            obj.copyFrom(obj.buffers); % copies the pointers to the data in "buf"
            
            if obj.debug_bit>5
                disp(['reading out batch ' num2str(obj.batch_counter) ' from buffer ' num2str(buf.index) ' | read_flag: ' util.text.print_vec(buf.this_buf.mex_flag_read)]);
            end
            
        end
        
        function batch_sync(obj, input)
            
            temp_images = zeros(obj.AOIwidth_c, obj.AOIheight_c, input.batch_size, 'uint16'); % notice the axes are inverted because this is C vs. Matlab! 
            obj.timestamps = zeros(input.batch_size,1);
            
            timeout = max([1, input.expT.*2]); % if exposure time is short, set timeout to 1 second, otherwise set it to 2*expT
            
            [rc] = AT_QueueBuffer(obj.hndl,obj.imageSizeBytes); AT_CheckWarning(rc);

            if strcmp(obj.getTriggerModeHW, 'Software')
                [rc] = AT_Command(obj.hndl,'SoftwareTrigger'); AT_CheckWarning(rc);
            end

            obj.t_start = util.text.time2str(datetime('now', 'TimeZone', 'UTC'));
            
            for ii = 1:input.batch_size
                
                buf = [];
                
                [rc, buf] = AT_WaitBuffer(obj.hndl, timeout*1000); % timeout in milliseconds! 
                
                if rc==0 % if we did not timeout or other error
                    
                    if ii<input.batch_size
                        if strcmp(obj.getTriggerModeHW, 'Software')
                            [rc] = AT_Command(obj.hndl,'SoftwareTrigger'); AT_CheckWarning(rc);
                            [rc] = AT_QueueBuffer(obj.hndl, obj.imageSizeBytes); AT_CheckWarning(rc); 
                        end
                    end
                    
                    [rc,buf2] = AT_ConvertMono16ToMatrix(buf, obj.AOIheight_c, obj.AOIwidth_c, obj.AOIstride); AT_CheckWarning(rc);
                    
                    temp_images(:,:,ii) = buf2;
                    
                    % Get timestamp and convert it into seconds
                    [rc,ticks] = AT_GetTimeStamp(buf, obj.imageSizeBytes); AT_CheckWarning(rc);
                    obj.timestamps(ii) = double(ticks)./double(obj.clockFreq);
                    
                else % if the acquisition got stuck... 
                    if rc==13
                        disp('ERROR 13: timeout...');
                    end
                    
                    obj.num_restarts = obj.num_restarts + 1;
                    disp(['Restarting acquisition... num_restart= ' num2str(obj.num_restarts)]);
                    obj.restart_sync;
                    
                end
                
            end
            
            obj.t_end = util.text.time2str(datetime('now', 'TimeZone', 'UTC'));
            obj.t_end_stamp = obj.getTimestampHW;
            
            temp_images(:,:,ii+1:end) = [];
            obj.timestamps(ii+1:end) = [];
            
            obj.images = temp_images;
            
            obj.buffers.input(obj); % pass the timestamps, images, t_start, t_end, t_end_stamp to buffers for saving
            
        end
        
        function restart_sync(obj)
            
            [rc] = AT_Command(obj.hndl,'AcquisitionStop'); AT_CheckWarning(rc);
            
            [rc] = AT_Flush(obj.hndl); AT_CheckWarning(rc);
            
            for X = 1:10
                [rc] = AT_QueueBuffer(obj.hndl, obj.imageSize); AT_CheckWarning(rc);
            end
                        
            [rc] = AT_Command(obj.hndl, 'AcquisitionStart'); AT_CheckWarning(rc);
            
            if strcmp(obj.getTriggerModeHW, 'Software')
                [rc] = AT_Command(obj.hndl,'SoftwareTrigger'); AT_CheckWarning(rc);
            end
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
    
    methods % hardware access
        
        function val = getTemperatureHW(obj)
           
            [rc, val] = AT_GetFloat(obj.hndl, 'SensorTemperature'); AT_CheckWarning(rc);
             
        end
        
        function val = getExpTimeHW(obj)
           
            [rc, val] = AT_GetFloat(obj.hndl, 'ExposureTime'); AT_CheckWarning(rc);
                        
        end
        
        function [val_min, val_max] = getExpTimeLimitsHW(obj)
           
            [rc, val_min] = AT_GetFloatMin(obj.hndl, 'ExposureTime'); AT_CheckWarning(rc);
            [rc, val_max] = AT_GetFloatMax(obj.hndl, 'ExposureTime'); AT_CheckWarning(rc);
            
        end
        
        function setExpTimeHW(obj, val) % sets the lowest possible frame rate. User must then set the frame rate using "setFrameRateHW" based on the new expT
            
            if isempty(val)
                error('Must input a valid exposure time');
            end
            
            rc = AT_SetFloat(obj.hndl, 'ExposureTime', val); AT_CheckWarning(rc);
            
            [~, f_max] = obj.getFrameRateLimitsHW;
            
            rc = AT_SetFloat(obj.hndl, 'FrameRate', f_max*0.95); AT_CheckWarning(rc);
            
        end
        
        function val = getFrameRateHW(obj)
            
            if ~util.text.cs(obj.getTriggerModeHW, 'internal') % if we are not in "internal" trigger mode, there is no meaning to "frame rate"
                val = [];
            else
                [rc, val] = AT_GetFloat(obj.hndl, 'FrameRate'); AT_CheckWarning(rc);
            end
            
        end
        
        function [val_min, val_max] = getFrameRateLimitsHW(obj)
           
            [rc, val_min] = AT_GetFloatMin(obj.hndl, 'FrameRate'); AT_CheckWarning(rc);
            [rc, val_max] = AT_GetFloatMax(obj.hndl, 'FrameRate'); AT_CheckWarning(rc);
            
        end
        
        function setFrameRateHW(obj, val)
            
            [~,f_max] = obj.getFrameRateLimitsHW;
            
            if val>=f_max % cannot set frame rate to higher than top frame rate (set by expT)
                str = sprintf('Cannot set frame rate to %g. Setting to %g instead', val, f_max.*0.95);
                disp(str);
                obj.log.input(str);
                val = f_max.*0.95;
            end 
            
            if isempty(val)
                obj.setTriggerModeHW('software');
            else
                obj.setTriggerModeHW('internal');
                rc = AT_SetFloat(obj.hndl, 'FrameRate', val); AT_CheckWarning(rc);
            end
            
        end
        
        function val = getROI_HW(obj)
            
            [rc, t] = AT_GetInt(obj.hndl, 'AOITop'); AT_CheckWarning(rc);
            [rc, l] = AT_GetInt(obj.hndl, 'AOILeft'); AT_CheckWarning(rc);
            [rc, h] = AT_GetInt(obj.hndl, 'AOIHeight'); AT_CheckWarning(rc);
            [rc, w] = AT_GetInt(obj.hndl, 'AOIWidth'); AT_CheckWarning(rc);
            
            val = [l, t, w, h]; % flip axis from C to matlab
            
        end
        
        function val = is_zoomed_HW(obj)
            
            % flip axis from C to matlab
            [rc, w] = AT_GetInt(obj.hndl, 'AOIHeight'); AT_CheckWarning(rc);
            [rc, h] = AT_GetInt(obj.hndl, 'AOIWidth'); AT_CheckWarning(rc);
            
            if h==obj.max_height && w==obj.max_width
                val = 0;
            else
                val = 1;
            end
            
        end
        
        function val = maxWidthHW(obj)
            
            [rc, val] = AT_GetIntMax(obj.hndl, 'AOIHeight'); AT_CheckWarning(rc); % flip axis from C to matlab
            
        end
        
        function val = maxHeightHW(obj)
            
            [rc, val] = AT_GetIntMax(obj.hndl, 'AOIWidth'); AT_CheckWarning(rc); % flip axis from C to matlab
            
        end
        
        function setROI_HW(obj, varargin)
            
            if isscalar(varargin) && isnumeric(varargin{1}) && length(varargin{1})==4
                t = varargin{1}(1);
                l = varargin{1}(2);
                h = varargin{1}(3);
                w = varargin{1}(4);
            elseif length(varargin)==4
                t = varargin{1};
                l = varargin{2};
                h = varargin{3};
                w = varargin{4};
            else
                error('Supply the ROI to be set in the format <left, top, width, height>.');
            end
            
%             fprintf('l= %f | t= %f | w= %f | h= %f\n', l, t, w, h);
            
            % flip axis from C to matlab
            rc = AT_SetInt(obj.hndl, 'AOIHeight', w); AT_CheckError(rc);
            rc = AT_SetInt(obj.hndl, 'AOIWidth', h); AT_CheckError(rc);
            rc = AT_SetInt(obj.hndl, 'AOITop', l); AT_CheckError(rc);
            rc = AT_SetInt(obj.hndl, 'AOILeft', t); AT_CheckError(rc);
            
        end
        
        function val = getCycleModeHW(obj)
            
            [rc, ind] = AT_GetEnumIndex(obj.hndl, 'CycleMode'); AT_CheckWarning(rc);
            [rc, val] = AT_GetEnumStringByIndex(obj.hndl, 'CycleMode', ind , 100); AT_CheckWarning(rc);
            
        end
        
        function setCycleModeHW(obj, mode)
           
            import util.text.cs;
            
            if cs(mode, 'continuous')
                mode = 'Continuous';
            elseif cs(mode, 'fixed')
                mode = 'Fixed';
            else
                error(['unknown cycle mode: ' mode ]);
            end
            
            rc = AT_SetEnumString(obj.hndl, 'CycleMode', mode); AT_CheckWarning(rc);
            
        end 
        
        function val = getTriggerModeHW(obj)
            
            [rc, ind] = AT_GetEnumIndex(obj.hndl, 'TriggerMode'); AT_CheckWarning(rc);
            [rc, val] = AT_GetEnumStringByIndex(obj.hndl, 'TriggerMode', ind , 100); AT_CheckWarning(rc);
            
        end
        
        function setTriggerModeHW(obj, mode)
           
            import util.text.cs;
            
            if cs(mode, 'internal')
                mode = 'Internal';
            elseif cs(mode, 'software')
                mode = 'Software';
            else
                error(['unknown trigger mode: ' mode ]);
            end
            
            rc = AT_SetEnumString(obj.hndl, 'TriggerMode', mode); AT_CheckWarning(rc);
            
        end
        
        function val = getTimestampHW(obj)
            
%             [rc, obj.clockFreq] = AT_GetInt(obj.hndl, 'TimestampClockFrequency'); AT_CheckWarning(rc)
            
            [rc, ticks] = AT_GetInt(obj.hndl, 'TimestampClock'); AT_CheckWarning(rc)
            
            val = double(ticks)/obj.clockFreq;
            
        end
        
        function val = getGain(obj)
            
            val = 0.6; % there must be a better way to do this??
            val = []; % I don't think 0.6 is reliable enough... 
            
        end
        
        function val = getCameraNameHW(obj)
            
            [rc, val] = AT_GetString(obj.hndl, 'CameraName', 100); AT_CheckWarning(rc);
            
        end
        
        function str_out = printout_HW(obj, full)
            
            if nargin<2 || isempty(full)
                full = 0;
            end
            
            % cycle = obj.getCycleModeHW;
            trigger = obj.getTriggerModeHW;
            
            roi = obj.getROI_HW;
            
            temp = obj.getTemperatureHW;
            f = obj.getFrameRateHW;
            T = obj.getExpTimeHW;
            
            str = sprintf('trig= %s | T= %5.3f | f= %5.2f | temp= %4.2f | ROI= %s\n', trigger, T, f, temp, util.text.print_vec(roi));
    
            if full
                
                [rc, baseline] = AT_GetInt(obj.hndl, 'Baseline'); AT_CheckWarning(rc);
                [rc, pixwidth] = AT_GetFloat(obj.hndl, 'PixelWidth'); AT_CheckWarning(rc);  
                [rc, readtime] = AT_GetFloat(obj.hndl, 'ReadoutTime'); AT_CheckWarning(rc);
                [rc, noisefilter] = AT_GetBool(obj.hndl, 'SpuriousNoiseFilter'); AT_CheckWarning(rc);
                [rc, blemish] = AT_GetBool(obj.hndl, 'StaticBlemishCorrection'); AT_CheckWarning(rc);
                [rc, cooling] = AT_GetBool(obj.hndl, 'SensorCooling'); AT_CheckWarning(rc);

                str2 = sprintf('baseline= %d | pixwidth= %f | readtime= %f | noise_filter= %d | blemish remove= %d | cooling= %d\n', baseline, pixwidth, readtime, noisefilter, blemish, cooling);
                
                str = [str str2];
                
            end
    
            if nargout>0
                str_out = str;
            else
                disp(str);
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function makeGUI(obj)
           
            if isempty(obj.gui)
                obj.gui = obs.cam.gui.AndorGUI(obj);                
            end
            
            obj.gui.make;
            
        end
        
        function show(obj, varargin)
            
            import util.text.cs;
            
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
            
            I = obj.images(:,:,end);    
            
            delete(findobj(input.ax, 'type', 'rectangle'));
            
            h = findobj(input.ax, 'type','image');
            
            if isempty(h) 
                if ~isempty(obj.gui) && obj.gui.check
                    obj.gui.makeAxes;
                    input.ax = obj.gui.axes_image;
                    h = findobj(input.ax, 'type','image');
                end
                
                if isempty(h)
                    imagesc(input.ax, I);
                    axis(input.ax, 'image');
                    colorbar(input.ax);
                end
                
                h = findobj(input.ax, 'type','image');
                
            end

            if ~isempty(h)

                if obj.use_show_flipped
                    I = rot90(I, 2); % use this when viewing after the meridian flip
                end

                h.CData = I;

            end
            
        end
        
        function showROI(obj, roi, ax)
            
            if nargin<2 || isempty(roi)
                roi = obj.zoom2roi(obj.im_size_zoomed, obj.center_region_zoomed);
            end
            
            if nargin<3 || isempty(ax)
                
                if ~isempty(obj.gui) && obj.gui.check
                    ax = obj.gui.axes_image;
                end
                
            end
            
            h = findobj(ax, 'type','image');
            
            if isempty(h) || size(h.CData,1)~=obj.max_height || size(h.CData,2)~=obj.max_width
                return; % don't show a rectangle unless it is on the full frame image... 
            end
                
            delete(findobj(ax, 'type', 'rectangle'));
            
            rectangle(ax, 'Position', [roi(2) roi(1), roi(4), roi(3)]);
            
        end
        
    end    
    
end

