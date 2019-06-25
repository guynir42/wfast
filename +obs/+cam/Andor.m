classdef Andor < file.AstroData
% Control an Andor Camera, e.g., Zyla, Marana, Balor. 
% If you only want to record some images, use run(...) with optional arguments.
% To run this object under some analysis code, use startup(...) to give 
% parameters, then use "batch" to generate some images. Finish the run by 
% calling "finishup". 
% Each batch consists of filling properties "images" and "timestamps" and 
% a few timing parameters ("t_start", "t_end", "t_end_stamp"). See the 
% description of file.AstroData for more details. 
%hutter
% The camera can run in synchronous/async modes. 
% Synchronuous mode just holds up the main thread until all images are captured. 
% Async mode lets the camera run on a separate thread using C++/mex. 
% In this case the BufferWheel is filled one buffer at a time while the user 
% is free to run analysis code or send the buffers to be saved. 
% Note: in synchronuos mode the buffers are also being filled, to allow saving. 
%
% CONTROL PARAMETERS: Each switch/control parameter has a private version 
% marked with a trailing underscore (e.g., num_batches_). These are filled 
% with the user defined parameters by calling "stash_parameters" which is 
% called in "startup". Then a util.text.InputVars is created to parse any
% user inputs to "run" or "startup". These override the parameters of the 
% object, which are then used for the entire run. When "finishup" is called, 
% the stashed values are returned using "unstash_parameters". 
% This lets the user define different, temporary parameter setups without 
% changing the default values the camera uses. 
% Example: change the <batch_size> property when running live view by calling
% live('batch_size', 1), but when the live run is done, <batch_size> will 
% return to its original value (not necessarily 100, which is the global 
% default, but any value given to the camera, e.g., using the GUI). 
%
% LIST OF CONTROL PARAMETERS:
%   *mode: can choose "science" (default), "dark" or "flat". 
%   *num_batches: how many batches to capture (default=2).
%   *batch_size: how many frames in each batch (default=100).
%   *expT: exposure time, in seconds (default=0.025).
%   *expT_deep: used for preview images (default is 1 second). 
%   *frame_rate: when this is NaN, camera will run in TriggerMode='software'
%                which means a frame is taken as soon as camera is ready. 
%                If frame_rate is a number, then camera will trigger on this
%                frame rate, even if it leaves some dead time. 
%   *use_async: use the C++/mex interface for multithreaded capture (default=true).
%   *use_show: show the last image from each batch on GUI axis or on current axis (default=true).
%   *use_save: save full frame images to file. Used for "record" command (default=false). 
%   *use_roi: Take Region Of Interest (ROI) images instead of full frame. 
%             The actual position and size of the ROI is given by the next
%             two parameters. The default is false. 
%   *im_size: define the size of the image when using ROI (default: 512). 
%             When use_roi=0 the im_size shows the maximum height/width. 
%             To access the zoomed value even if camera is not in ROI, use 
%             the hidden variable <im_size_zoomed>. 
%   *center_region: the y and x position of the center of the ROI. 
%                   When empty (default) it will give the center of the image. 
%                   When use_roi=0 the center_region shows the center of the image. 
%                   To access the zoomed value even if camera is not in ROI, use 
%                   the hidden variable <center_region_zoomed>.
%   *debug_bit: determines the level of debugging info printed onscreen (default=1).
%   *log_level: determines the level of info output to text file (default=1).
%   *use_audio: turn on/off the audio queues when starting/finishing a run (default=false). 
%   *use_progress: turn on/off the progress bar printing onscreen (default=false). 
%   *pass_show: cell array with parameters to pass to "show" function (default={}). 
%
% NOTE: There are 3 levels of input parameters for running the camera:
%       1) Hardware parameters, as the camera was setup when it last took 
%          images. These is only relevant when a run has started, after many
%          changes to parameters are given from object/user parameters.
%       2) Object defined parameters, as given by setting object properties 
%          or by changing them via GUI. These override (1) at the beginning
%          of each run. 
%       3) User defined parameters, as given by varargin pairs to "run" or 
%          "startup" commands. These will override what is given in (2) at
%          the beginning of each run. 
%
% EXAMPLE USE CASES:
% 1) Running a recording session using object defined parameters.
%   >> a = obs.cam.Andor;
%   >> a.num_batches = 20; % number of batches for all following runs. 
%   >> a.record; % now proceed to record 20 batches. 
%
% 2) Recording some images using the varargin interface: 
%   >> a.record('num_batches', 20, 'batch_size', 100, 'use_save', 1); 
%
% 3) Use camera inside of analysis loop. 
%   >> a.startup('num_batches', 20, 'use_save', 0); 
%   >> for ii = 1:20
%   >>     a.batch;
%   >>     some_processing_function(a.images, a.timestamps); 
%   >> end
%   >> a.finishup; 
% 
% NOTE: use the syntax 
%       >> on_cleanup = onCleanup(@a.finishup);
%       right after the call to "startup" to make sure finishup is called 
%       even if there was an error or user ctrl+C break during the loop. 
%
    
    properties(Transient=true)
        
        gui;
        
        audio@util.sys.AudioControl;
        prog@util.sys.ProgressBar;
        runtime_buffer@util.vec.CircularBuffer; % keeps track of the last few batches' runtime and number of frames, to calculate frame_rate_measured
        
    end
    
    properties % objects
        
        pars@head.Parameters; % parameters used in the observations
        
        buffers@file.BufferWheel;
        
        focuser;
        
        log@util.sys.Logger;
        
        hndl; % pointer to the Andor SDK handle of the camera
        
    end
    
    properties % inputs/outputs
        
        status = 0; % if 0 there is some problem communicating with hardware
        frame_rate_measured; % average of a few batches
        frame_rate_camera; % frame rate from timestamps as given from camera (not including dead time and read time)
        
    end
    
    properties % switches/controls
        
        mode = 'science'; % can also choose "dark" or "flat"
        batch_size = 100; % how many frames in each batch
        num_batches = 2; % how many batches before finishing the run
        
        expT = 0.03; % exposure time (seconds)
        expT_deep = 1; % for previews (this is never stashed or updated from InputVars)
        frame_rate = 25; % frame rate (Hz) used for general capures. Use NaN to let camera take images as soon as it can, regardless of timing. 
        frame_rate_live = 6; % use this frame rate as the default for "live". 
        
        use_async = true; % use the new C++/mex interface for multithreaded capture
        use_show = true; % show the last image of each batch
        use_save = false; % save full frame images
        
        % Region Of Interest
        use_roi = false; % when enabled the camera will use im_size and center_region
        
        pass_show = {}; % parameters to pass to "show" function
        
        %...remember to stash any additional parameters! 
        
        % display/gui parameters 
        use_show_flipped = 0; % flip the display for North is up after meridian flip
        
        debug_bit = 1; % level of debugging output on screen
        log_level = 1; % level of output to text file
        use_audio = false; % audio queues for start/finish runs
        use_progress = false; % show progress bar onscreen
        
        % add new switches also to makeInputVars, stash_parameters and unstash_parameters...
        
    end
    
    properties(Dependent=true)
        
        % these parameters are updated with respect to use_roi...
        im_size; % size of the ROI (if scalar then width=height). If use_roi=0 will output maximum height/width
        center_region; % (y then x) of the center of the ROI. If empty, choose the center of the field
        ROI; % [top, left, height, width] as given to camera
        
    end
    
    properties(Hidden=true)
       
        mex_flag = [0 0 0 0]; % for starting and stopping the camera...
        
        max_height; % lazy loaded once each run from camera hardware
        max_width; % lazy loaded once each run from camera hardware
        
        mode_list = {'science', 'dark', 'flat'}; % possible modes for camera
        mode_index = 1; % which mode on the list we are in right now
        
        % these are used when GUI is given an empty value
        default_num_batches; 
        default_batch_size; 
        default_expT; 
        default_expT_deep;
        default_frame_rate; 
        default_frame_rate_live;
        
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
        
        is_running = 0;
        brake_bit = 1; % when 1 the camera is stopped. Is set to 0 on "startup". 
        
        version = 1.01;
        
    end
    
    properties(Hidden=true, Transient=true) % stashed parameters...
        
        % these are stashed each new run, so user can override the parameters
        % of the object. When run is finished, these will be returned to the
        % visible, object parameters for the next runs. 
        mode_;
        batch_size_;
        num_batches_;
        frame_rate_;
        expT_;
        use_async_;
        use_show_;
        use_save_;
        use_roi_;
        im_size_zoomed_;
        center_region_zoomed_;
        debug_bit_;
        log_level_;
        use_audio_;
        use_progress_;
        pass_show_;
        
    end
    
    methods % constructor
        
        function obj = Andor(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.cam.Andor')
                if obj.debug_bit, fprintf('Andor camera copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Andor camera constructor v%4.2f\n', obj.version); end
                
                obj.log = util.sys.Logger('Andor_camera'); % start the logger before anything else
                
                util.oop.save_defaults(obj); % store default values from class definitions
                
                obj.connect; % initialize library and get a camera handle
               
                obj.log.input('Loading camera auxiliary objects');
                
                try 
                    
                    obj.pars = head.Parameters;

                    obj.setupBuffers;
                    
                    obj.setupFocuser;
                    
                    obj.setupAudio;
                    
                    obj.setupProgressBar;
                    
                    obj.setupTimeBuffer;
                    
                    obj.update;
                    
                    try 
                        pause(5);
                        obj.setupDefaultFocusPosition;
                    catch ME
                        warning(ME.getReport);
                    end

                catch ME
                    obj.log.error(ME.getReport);
                    rethrow(ME);
                end
                
            end
            
        end

        function connect(obj) % initialize library and get a camera handle 
            
            obj.log.input('Connecting to Andor camera...');
            
            try 
                    
                rc = obs.cam.sdk.AT_InitialiseLibrary; obs.cam.sdk.AT_CheckError(rc);
                    
                [rc, obj.hndl] = obs.cam.sdk.AT_Open(0); obs.cam.sdk.AT_CheckError(rc);
%                 obj.hndl = obs.cam.mex_new.connect; % the new mex code is better than the matlab SDK because...?
                
                [rc] = obs.cam.sdk.AT_SetBool(obj.hndl, 'SensorCooling', 1); obs.cam.sdk.AT_CheckWarning(rc);        
              
                rc = obs.cam.sdk.AT_SetBool(obj.hndl, 'SpuriousNoiseFilter', 0); obs.cam.sdk.AT_CheckWarning(rc);
                rc = obs.cam.sdk.AT_SetBool(obj.hndl, 'StaticBlemishCorrection', 0); obs.cam.sdk.AT_CheckWarning(rc);
                
%                 [rc] = obs.cam.sdk.AT_SetEnumString(obj.hndl,'ElectronicShutteringMode','Rolling'); obs.cam.sdk.AT_CheckWarning(rc);
                
                [rc, val] = obs.cam.sdk.AT_IsEnumIndexImplemented(obj.hndl, 'SimplePreAmpGainControl',2); obs.cam.sdk.AT_CheckWarning(rc);
                if(val)
                    [rc] = obs.cam.sdk.AT_SetEnumString(obj.hndl,'SimplePreAmpGainControl','16-bit (low noise & high well capacity)'); obs.cam.sdk.AT_CheckWarning(rc);
                end
                [rc] = obs.cam.sdk.AT_SetEnumString(obj.hndl,'PixelEncoding','Mono16'); obs.cam.sdk.AT_CheckWarning(rc);
                                
                % Enable Metadata
                [rc] = obs.cam.sdk.AT_SetBool(obj.hndl,'MetadataEnable',1); obs.cam.sdk.AT_CheckWarning(rc);
                [rc] = obs.cam.sdk.AT_SetBool(obj.hndl,'MetadataTimestamp',1); obs.cam.sdk.AT_CheckWarning(rc);
                
                [rc] = obs.cam.sdk.AT_SetEnumString(obj.hndl,'CycleMode','Continuous'); obs.cam.sdk.AT_CheckWarning(rc); % I think we don't need to limit the number of frames (there are explicit stops in the code).
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function disconnect(obj) % disconnect the current camera handle
            
            obj.log.input('Disconnecting Andor camera...');
            
            try 
                obs.cam.mex_new.disconnect(obj.hndl);
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function reconnect(obj) % disconnect then connect
            
            obj.disconnect;
            pause(0.05);
            obj.connect;
            
        end
        
        function setupBuffers(obj) % creat a buffer wheel and give it the camera's mex_flag
            
            obj.buffers = file.BufferWheel(5, obj.pars);
            obj.buffers.camera_mex_flag = obj.mex_flag;

        end
        
        function setupFocuser(obj) % try to connect to focuser
           
            try
                obj.focuser = obs.focus.FocusSpider;
            catch ME
                disp('Cannot connect to focuser, using simulator instead');
                warning(ME.getReport);
                obj.log.input('Cannot connect to focuser, using simulator instead');
                obj.log.input(ME.getReport);
                obj.focuser = obs.focus.Simulator;
                
            end
 
        end
        
        function setupDefaultFocusPosition(obj)
            
%             obj.focuser.pos = 4.69; % from our focus test on 26/5/19
            obj.focuser.pos = 4.5857; % updated at 28/5/19
            obj.focuser.pos = 4.658; % updated at 03/6/19
%             obj.focuser.tip = 3.25;
            obj.focuser.tip = 0;
%             obj.focuser.tilt = 1.05;
            obj.focuser.tilt = 0;
            
        end
        
        function setupAudio(obj) % try to load the audio player
            
            try
                obj.audio = util.sys.AudioControl;
            catch ME
                disp('Cannot connect to audio driver. Continuing without it!');
                obj.log.input('Cannot connect to audio driver. Continuing without it!');
            end
            
        end
        
        function setupProgressBar(obj) % make a progress bar object
            
            obj.prog = util.sys.ProgressBar;
            
        end
        
        function setupTimeBuffer(obj) % make a circular buffer for number of frames and runtime
            
            obj.runtime_buffer = util.vec.CircularBuffer;
            obj.runtime_buffer.titles = {'num_frames', 'time'};
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj) % get ready for a new run

            if obj.brake_bit==0
                warning('Cannot reset camera while it is running! (set brake_bit=1)'); 
                return;
            end
            
            obj.buffers.reset; 
            util.vec.mex_change(obj.mex_flag, 4, 0); % reset the counter on the mex_flag
            
            obj.runtime_buffer.reset;
            
            obj.frame_rate_camera = [];
            obj.frame_rate_measured = [];
            
            obj.max_height = [];
            obj.max_width = [];
            
            obj.batch_counter = 0;
            
        end
        
        function clear(obj) % removes outputs and intermediary data
            
            % maybe clear the buffers??
            % clearing the buffers will probably crash the camera mex code
            
        end
        
    end
    
    methods % getters
        
        function val = is_finished(obj) % when run is finished (brake_bit is 1). Used as indicator for controlling objects (e.g., Acquisition)
            
            val = obj.brake_bit;
            
        end
        
        function val = get.max_height(obj) % lazy load the maximum height from camera hardware
            
            if isempty(obj.max_height)
                obj.max_height = obj.maxHeightHW;
            end
            
            val = obj.max_height;
            
        end
        
        function val = get.max_width(obj) % lazy load the maximum width from camera hardware
            
            if isempty(obj.max_width)
                obj.max_width = obj.maxWidthHW;
            end
            
            val = obj.max_width;
            
        end
        
        function val = get.im_size(obj) % if empty (or when use_roi=0) gives the max frame size
            
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
        
        function val = get.center_region(obj) % if empty (or when use_roi=0) gives the center of the full-frame image
            
            if obj.use_roi==0
                val = floor([obj.max_height, obj.max_width]/2)+1;
            elseif isempty(obj.center_region_zoomed)
                val = floor([obj.max_height, obj.max_width]/2)+1;
            else
                val = obj.center_region_zoomed;
            end
            
        end
        
        function val = get.ROI(obj) % if use_roi=1, get user-defined ROI. Otherwise give the full-frame region
            
            if obj.use_roi==0
                val = [1 1 obj.max_height obj.max_width];
            else
                val = obj.zoom2roi(obj.im_size, obj.center_region);
            end
            
        end
        
        function val = zoom2roi(obj, im_size, center) % translate im_size and center_region to ROI 4-vector, including shrinking the ROI if center is too close to the edges
        % if im_size is empty, use the predefined im_size_zoomed. 
        % if im_size is scalar, expand it to 2-element vector with identical elements. 
        % if center is empty, use center of full frame image. 
        % if center is scalar, expand it to 2-element vector with identical elements. 
        
            
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
        
        function [im_size, center] = roi2zoom(obj, roi) % translate the ROI into im_size and center_region values
            
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
            
                if ~isempty(center) && (center(1)~=floor(obj.max_height/2)+1 || center(2)~=floor(obj.max_width/2)+1)
                    val = sprintf('%s at %s', val, util.text.print_vec(center, ','));
                end
            
            end
            
            
        end
        
    end
    
    methods % setters
        
        function set.pars(obj, val) % make sure <pars> is given to <buffers> as well.
           
            obj.pars = val;
            
            if ~isempty(obj.buffers)
                obj.buffers.pars = val;
            end
        
        end
        
        function set.mex_flag(~, ~) % prevents users from altering the value of mex_flag... 
           
            disp('You must not change "mex_flag" directly. Use util.vec.mex_change or obj.stop or obj.unlockMexFlag');
            
        end
        
        function set.im_size(obj, val) % set the im_size used when zoomed in. If setting a value smaller than full-frame, will also set use_roi=1.
            
            if isempty(val) || (val(1)==obj.max_height && val(2)==obj.max_width) % input size is used to unzoom
                obj.use_roi = 0;
            else % change the zoomed size
                obj.im_size_zoomed = val;
                obj.use_roi = 1;
            end
            
        end
        
        function set.center_region(obj, val) % set the center_region used when zoomed in. 
            
            obj.center_region_zoomed = val; 
            
        end
        
        function set.center_region_zoomed(obj, val) % make sure any values entered to center_region_zoomed are rounded to integers. 
            
            obj.center_region_zoomed = round(val); 
            
        end
        
        function zoom(obj, varargin) % quick access to zoom in or out and set the zoom parameters. 
        % can accept multiple lengths of inputs, in the form of a string, 
        % a vector, or a comma separated list of numerical values (with 1 to 4 inputs): 
        %   *input is "unzoom" it will set use_roi=0
        %   *input is "zoom" or is empty, set use_roi=1 and keep existing zoom parameters
        %   *input is numeric vector of 1-4 elements or a list of 1-4 numeric arguments:
        %           1- set im_size height and width to given value. 
        %           2- set im_size height to 1st value and width to 2nd. 
        %           3- set im_size as before, set both center_region elements to 3rd value. 
        %           4- set im_size as before, set center(1) to 3rd value and center(2) to 4th value. 
        
            if isempty(varargin)
                obj.use_roi = 1;
            elseif isscalar(varargin)
                if isempty(varargin{1})
                    obj.use_roi = 1;
                elseif isnumeric(varargin{1})
                    if length(varargin{1})<=2
                        obj.im_size = varargin{1};
                        obj.use_roi = 1;
                    elseif length(varargin{1})==3
                        obj.im_size = varargin{1}(1:2);
                        obj.center_region = varargin{1}(3).*[1 1];
                        obj.use_roi = 1;
                    elseif length(varargin{1})==4
                        obj.im_size = varargin{1}(1:2);
                        obj.center_region = varargin{1}(3:4);
                        obj.use_roi = 1;
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
                obj.use_roi = 1;
            elseif length(varargin)==3
                obj.im_size(1) = varargin{1};
                obj.im_size(2) = varargin{2};
                obj.center_region(1) = varargin{3};
                obj.center_region(2) = varargin{3};
                obj.use_roi = 1;
            elseif length(varargin)==4
                obj.im_size(1) = varargin{1};
                obj.im_size(2) = varargin{2};
                obj.center_region(1) = varargin{3};
                obj.center_region(2) = varargin{4};
                obj.use_roi = 1;
            end
            
        end
        
        function unzoom(obj) % set use_roi=0 and leave zoom parameters as they are
            
            obj.use_roi = 0;
            
        end
        
    end
    
    methods % high level run commands
        
        function record(obj, varargin) % capture some batches and save them to file
            
            obj.run('save', 1, 'audio', 1, 'progress', 1, 'log_level', 1, 'async', 1, 'show', 1, varargin{:});
            
        end
        
        function live(obj, varargin) % run for up to 1e6 batches of sinlge image, showing each one, at a constant rate (without saving images)
            
            obj.run('num_batches', 1e6, 'batch_size', 1, 'show', 1, 'save', 0,...
                'audio', 0, 'async', 0, 'progress', 0, 'log_level', 0, ...
                'frame rate', obj.frame_rate_live, varargin{:}); 
            
        end
        
        function preview(obj, varargin) % take a single, deep exposure, show it, and don't save it
            
            obj.run('num_batches', 1, 'batch_size', 1, 'expT', obj.expT_deep, 'frame_rate', NaN, ...
                'show', 1, 'save', 0, 'audio', 0, 'async', 0,...
                'progress', 0, 'log_level', 0, varargin{:}); 
            
        end
        
        function single(obj, varargin)
            
            obj.run('num_batches', 1, 'show', 1, 'save', 0, 'audio', 0, 'async', 0,...
                'progress', 0, 'log_level', 0, varargin{:}); 
            
        end
        
        function autofocus(obj, varargin) % not yet implemented
            
            error('This is not yet implemented!');
            
        end
        
    end
    
    methods % interface commands to camera
        
        function update(obj) % check that connection is still ok
            
            if obj.debug_bit>3, disp('updating camera'); end
            
            rc = obs.cam.sdk.AT_GetFloat(obj.hndl, 'ExposureTime'); % contact hardware to see if it is connected
            
            if rc % non-zero value indicates an error
                obj.status = 0; 
            else
                obj.status = 1;
            end
            
            % maybe add some more tests to see if camera is alive??
            
        end
        
        function run(obj, varargin) % calls "startup", loops over calls to "batch", then calls "finishup"
            
            if obj.is_running
                disp('Camera is already running. Set is_running to zero...');
                return;
            else
                obj.is_running = 1;
            end
            
            obj.startup(varargin);
            
            try
                
                on_cleanup = onCleanup(@obj.finishup);
                
                for ii = 1:obj.num_batches
                    
                    if obj.brake_bit, break; end
                    
                    obj.batch;
                    
                end
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function startup(obj, varargin)
            
            try 
                
                if obj.brake_bit==0
                    disp('Cannot startup camera while it is still runnning (turn off brake_bit)');
                end
                
                input = obj.makeInputVars(varargin{:});
                
                if input.log_level
                    obj.log.input(sprintf('Starting a new run with %d batches of %d images each. expT= %f, async= %d, save= %d, ROI= %s',...
                        input.num_batches, input.batch_size, input.expT, input.use_async, input.use_save, obj.getZoomStr(input.im_size, input.center_region))); 
                end

                obj.stash_parameters(input); % put the user-defined parameters into hidden "stashed" parameters to be loaded at end of run
                
                obj.reset;
                
                if isempty(obj.audio)
                    obj.setupAudio;
                end
                
                if isempty(obj.prog)
                    obj.setupProgressBar;
                end

                if isempty(obj.runtime_buffer)
                    obj.setupTimeBuffer;
                end
                
                % make sure hardware is updated... 
                obj.setROI_HW(obj.zoom2roi(obj.im_size, obj.center_region)); % if use_roi=0 these parameters just give full-frame
                obj.setExpTimeHW(obj.expT); 
                obj.setFrameRateHW(obj.frame_rate);

                obj.setShutterModeHW('rolling'); % maybe add this as an optional argument?
                
                if isempty(obj.frame_rate) || isnan(obj.frame_rate) % in this mode the camera takes an image as soon as it gets a command to "software trigger"
                    [rc] = obs.cam.sdk.AT_SetEnumString(obj.hndl,'TriggerMode','Software'); obs.cam.sdk.AT_CheckWarning(rc);
                else % in this mode there is a fixed frame rate, so the frame rate may be lower than the maximum 
                    [rc] = obs.cam.sdk.AT_SetEnumString(obj.hndl,'TriggerMode','Internal'); obs.cam.sdk.AT_CheckWarning(rc);
                end

                obj.check_inputs; % check the hardware/input configuration is compatible
                obj.update_pars; % update the Parameter object
                obj.update_buffers; % upodate the Buffer object

                obj.num_restarts = 0; % how many times did the camera have to be restarted (synchronous mode only!)

                obj.batch_counter = 0; % how many batches were already taken (for GUI display, etc)

                if obj.use_save
                    filename = obj.buffers.getReadmeFilename;
                    util.oop.save(obj, filename, 'name', 'camera');  
                end

                if obj.use_audio
                    try 
                        obj.audio.playTakeForever;
                    catch ME
                        warning(ME.getReport);
                    end
                end

                if obj.use_progress
                    obj.prog.start(input.num_batches);
                end

                if obj.use_async

                    if obj.mex_flag(1)
                        obj.stop;
                        error('Camera is already running. Stopping and restarting...');
                    end

                    obj.unlockMexFlag; % this should be replaced by a more carefull unlocking inside buffer>reset

                    for ii = 1:length(obj.buffers.buf)
                        util.vec.mex_change(obj.buffers.buf(ii).mex_flag_record, 1, 1); % lock all the buffers for recording, they will unlock once the camera fills them... 
                    end

                    obs.cam.mex_new.startup(obj, obj.mex_flag, obj.buffers.buf, obj.buffers.index_rec_vec, obj.num_batches, obj.batch_size); % call the mex file for async recording

                else % synchronous startup option

                    if strcmp(obj.getCycleModeHW, 'Fixed') % this shouldn't happen (we don't use this mode anymore)
                        [rc] = obs.cam.sdk.AT_SetInt(obj.hndl,'FrameCount',obj.batch_size); obs.cam.sdk.AT_CheckWarning(rc);
                    end

                    % save these parameters from hardware to save calls whenever we need to parse the raw data
                    [rc, obj.imageSizeBytes] = obs.cam.sdk.AT_GetInt(obj.hndl, 'ImageSizeBytes'); obs.cam.sdk.AT_CheckWarning(rc); % size of raw data buffer for single image

                    % note width and height are flipped due to C -> matlab conventions
                    [rc, obj.AOIwidth_c] = obs.cam.sdk.AT_GetInt(obj.hndl, 'AOIWidth'); obs.cam.sdk.AT_CheckWarning(rc);
                    [rc, obj.AOIheight_c] = obs.cam.sdk.AT_GetInt(obj.hndl, 'AOIHeight'); obs.cam.sdk.AT_CheckWarning(rc);

                    [rc, obj.AOIstride] = obs.cam.sdk.AT_GetInt(obj.hndl, 'AOIStride'); obs.cam.sdk.AT_CheckWarning(rc);

                    [rc, obj.clockFreq] = obs.cam.sdk.AT_GetInt(obj.hndl, 'TimestampClockFrequency'); obs.cam.sdk.AT_CheckWarning(rc);

                    [rc] = obs.cam.sdk.AT_Flush(obj.hndl); obs.cam.sdk.AT_CheckWarning(rc); % clear any leftover hardware buffers from last run

                    for X = 1:10 % setup new hardware buffers
                        [rc] = obs.cam.sdk.AT_QueueBuffer(obj.hndl, obj.imageSizeBytes); obs.cam.sdk.AT_CheckWarning(rc);
                    end

                    [rc] = obs.cam.sdk.AT_Command(obj.hndl, 'TimestampClockReset'); obs.cam.sdk.AT_CheckWarning(rc); % reset the time stamp clock to now

                    [rc] = obs.cam.sdk.AT_Command(obj.hndl, 'AcquisitionStart'); obs.cam.sdk.AT_CheckWarning(rc); % start rolling the camera

                end

                obj.brake_bit = 0; % this allows all the loops to continue. It becomes 1 if the GUI stop button is pressed.

                if ~isempty(obj.gui)
                    obj.gui.update; % update after everything is set up
                end

            catch ME
                obj.log.error(ME.getReport);
                obj.is_running = 0;
                rethrow(ME);
            end
                
        end
        
        function finishup(obj) % called at end of each run. Shuts down everything (when use_save=1 will also add a finish readme file)
            
            try 
                if obj.use_async

                    obj.stop;

                else

                    [rc] = obs.cam.sdk.AT_Command(obj.hndl,'AcquisitionStop'); obs.cam.sdk.AT_CheckWarning(rc);

                    [rc] = obs.cam.sdk.AT_Flush(obj.hndl); obs.cam.sdk.AT_CheckWarning(rc);

                end

                if obj.use_save
                    filename = obj.buffers.getReadmeFilename('Z');
                    util.oop.save(obj, filename, 'name', 'camera');  
                end

                if obj.use_audio
                    try
                        obj.audio.playShowsOver;
                    catch ME
                        warning(ME.getReport);
                    end
                end

                if obj.use_progress
                    obj.prog.finish(obj.batch_counter);
                end

                obj.brake_bit = 1;

                obj.unstash_parameters; % return all parameters to values they were in before this run
                
                obj.is_running = 0;
                
                if ~isempty(obj.gui)
                    obj.gui.update;
                end
            
            catch ME
                obj.log.error(ME.getReport);
                obj.is_running = 0;
                rethrow(ME);
            end
        end
        
        function batch(obj)
            
            if obj.use_async
                obj.batch_async;
            else
                obj.batch_sync;
            end
            
            obj.batch_counter = obj.batch_counter + 1;
            
            if obj.use_show && ~isempty(obj.gui) && obj.gui.check
                obj.show(obj.pass_show{:});
            end

            obj.frame_rate_camera = size(obj.images,3)./(obj.t_end_stamp - obj.timestamps(1)); % how many frames in what time (not including dead time and read time)
            
            if ~isempty(obj.buffers.t_start)
                
                time = util.text.str2time(obj.buffers.t_start);
                dt = seconds(time - obj.pars.ephem.time); % how much time passed since last update
                obj.pars.ephem.time = time; % update the "pars" object with current start time
                
                obj.runtime_buffer.input([size(obj.images,3), dt]); % how many images, how much time passed (total)
                
                obj.frame_rate_measured = sum(obj.runtime_buffer.data(:,1))./sum(obj.runtime_buffer.data(:,2)); % average over the last few (10) batches
                
            end
            
            if obj.use_save
                obj.buffers.save;
            end
            
            obj.buffers.nextBuffer; % make sure to rotate the buffers
            
            if ~isempty(obj.gui)
                obj.gui.update;
            end
            
            if obj.use_progress
                obj.prog.showif(obj.batch_counter); % show progress bar onscreen
            end
            
            drawnow; % make sure to update graphic objects and get callback interrupts
            
        end
        
    end
    
    methods % internal functions 
        
        function input = makeInputVars(obj, varargin) % parse a varargin for keyword-value pairs or an InputVars object, outputs an InputVars with defaults and overrides by user inputs
            
            idx = util.text.InputVars.isInputVars(varargin);
            
            if ~isempty(varargin) && any(idx)
                idx = find(idx, 1, 'first');
                input = varargin{idx}; 
                varargin(find(idx, 1, 'first')) = []; % remove the one cell with InputVars in it
            else
                input = util.text.InputVars;
                input.input_var('mode', obj.mode); % science, dark, flat
                input.input_var('num_batches', obj.num_batches, 'Nbatches');
                input.input_var('batch_size', obj.batch_size);
                input.input_var('frame_rate', obj.frame_rate, 'frequency');
                input.input_var('expT', obj.expT, 'exposure time', 'expT'); 
                input.input_var('use_roi', obj.use_roi); % turn on/off the ROI mode
                input.input_var('im_size', obj.im_size); % set the height and width of the ROI
                input.input_var('center_region', obj.center_region, 'center', 'center_ROI'); % set the center of the ROI (set to center of field if empty)
                input.input_var('use_async', obj.use_async, 'async', 5); % use mex-async code or simple batch-by-batch interface
                input.input_var('use_show', obj.use_show, 'show', 5); % if you want to show each batch's images/stack
                input.input_var('use_save', obj.use_save, 'save', 5); % if you want to save each batch 
                input.input_var('use_audio', obj.use_audio, 'audio', 5); % turn on/off audio signals
                input.input_var('use_progress', obj.use_progress, 'progress', 5); % display a progress bar on screen
                input.input_var('log_level', obj.log_level); % choose if and how much logging you want for this run (1 is only start of run). Errors are always logged. 
                input.input_var('debug_bit', obj.debug_bit, 'debug');
                input.input_var('pass_show', obj.pass_show, 6); % parameters to pass to "show" function
            end
            
            input.scan_vars(varargin{:});
            
        end
        
        function stash_parameters(obj, input) % sets all camera parameters into hidden "stash" parameters. If given an InputVars object, will load parameters from it to the camera object. 
            
            obj.mode_ = obj.mode;
            obj.batch_size_ = obj.batch_size;
            obj.num_batches_ = obj.num_batches;
            obj.frame_rate_ = obj.frame_rate;
            obj.expT_ = obj.expT;
            obj.use_async_ = obj.use_async;
            obj.use_show_ = obj.use_show;
            obj.use_save_ = obj.use_save;
            obj.use_roi_ = obj.use_roi;
            obj.im_size_zoomed_ = obj.im_size_zoomed;
            obj.center_region_zoomed_ = obj.center_region_zoomed;
            obj.debug_bit_ = obj.debug_bit; 
            obj.log_level_ = obj.log_level;
            obj.use_audio_ = obj.use_audio;
            obj.use_progress_ = obj.use_progress;
            obj.pass_show_ = obj.pass_show;
            
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
            
            obj.mode = obj.mode_;
            obj.batch_size = obj.batch_size_;
            obj.num_batches = obj.num_batches_;
            obj.frame_rate = obj.frame_rate_;
            obj.expT = obj.expT_;
            obj.use_async = obj.use_async_;
            obj.use_show = obj.use_show_;
            obj.use_save = obj.use_save_;
            obj.use_roi = obj.use_roi_;
            obj.im_size_zoomed = obj.im_size_zoomed_;
            obj.center_region_zoomed = obj.center_region_zoomed_;
            obj.debug_bit = obj.debug_bit_;
            obj.log_level = obj.log_level_;
            obj.use_audio = obj.use_audio_;
            obj.use_progress = obj.use_progress_;
            obj.pass_show = obj.pass_show_;
            
        end
        
        function cycleModes(obj) % cycle through "science", "dark" and "flat" modes
            
            obj.mode_index = obj.mode_index + 1;
            if obj.mode_index>length(obj.mode_list)
                obj.mode_index = 1;
            end
            
            obj.mode = obj.mode_list{obj.mode_index};
            
        end
        
        function check_inputs(obj) % check that input parameters make sense
            
            if util.text.cs(obj.mode, 'dark', 'flat')
            
                if obj.is_zoomed_HW
                    error('We should never take dark/flat images in ROI!');
                end
                
            end
            
            [rc, cooling] = obs.cam.sdk.AT_GetBool(obj.hndl, 'SensorCooling'); obs.cam.sdk.AT_CheckWarning(rc);
            if ~cooling, error('Camera sensor cooling is off!'); end
            
            % add other checks for consistency in parameters... 
            
        end
        
        function update_pars(obj) % make sure "pars" object is updated with hardware values
            
%             obj.pars.datapath = obj.buffers.base_dir;
            
            obj.pars.expT = obj.getExpTimeHW;
            obj.pars.frame_rate = obj.getFrameRateHW;
            
            if ~isempty(obj.focuser)
                obj.pars.FOCUS_POS = obj.focuser.pos;
            end
            
%             roi = obj.getROI_HW;
%             obj.pars.AOI_top = roi(1);
%             obj.pars.AOI_left = roi(2);
%             obj.pars.AOI_height = roi(3);
%             obj.pars.AOI_width = roi(4);
            
            obj.pars.ROI = obj.getROI_HW;

            obj.pars.im_size = [obj.pars.ROI(3) obj.pars.ROI(4)];
            
            obj.pars.batch_size = obj.batch_size;
            obj.pars.gain = obj.getGain;
            obj.pars.instrument = obj.getCameraNameHW;
            
            obj.pars.update; % get the current time etc
            
            obj.pars.type = obj.mode;
            obj.pars.is_dark = util.text.cs(obj.mode, 'dark');
            obj.pars.is_flat = util.text.cs(obj.mode, 'flat');
            
        end
        
        function update_buffers(obj) % make sure "buffers" object is setup according to input parameters (e.g., product type)
            
            import util.text.cs;
            
            if cs(obj.mode, 'science')
                if obj.use_roi
                    obj.buffers.product_type = 'ROI';
                else
                    obj.buffers.product_type = 'Raw';
                end
            elseif cs(obj.mode, 'dark')
                obj.buffers.product_type = 'Dark';
            elseif cs(obj.mode, 'flat')
                obj.buffers.product_type = 'Flat';
            end
            
        end
        
        function batch_async(obj) % take a batch using the C++/mex multithread interface
            
            obj.buffers.waitForRecording; % waits for the camera to fill the buffers
            obj.buffers.vec2times; % update the timing data from t_vec to t_start, t_end and t_end_stamp
            obj.buffers.loadDataFromBuffer;
            obj.copyFrom(obj.buffers); % copies the pointers to the data in "buf"
            
            if obj.debug_bit>5
                disp(['reading out batch ' num2str(obj.batch_counter) ' from buffer '...
                    num2str(obj.buffers.index) ' | read_flag: ' util.text.print_vec(obj.buffers.this_buf.mex_flag_read)]);
            end
            
        end
        
        function batch_sync(obj) % take a batch using the basic SDK tools 
            
            temp_images = zeros(obj.AOIwidth_c, obj.AOIheight_c, obj.batch_size, 'uint16'); % notice the axes are inverted because this is C vs. Matlab! 
            obj.timestamps = zeros(obj.batch_size,1);
            
            timeout = max([1, obj.expT.*2]); % if exposure time is short, set timeout to 1 second, otherwise set it to 2*expT
            
            [rc] = obs.cam.sdk.AT_QueueBuffer(obj.hndl,obj.imageSizeBytes); obs.cam.sdk.AT_CheckWarning(rc);

%             if strcmp(obj.getTriggerModeHW, 'Software') % I think we can just give a SoftwareTrigger which is ignored in any other mode
                [rc] = obs.cam.sdk.AT_Command(obj.hndl,'SoftwareTrigger'); obs.cam.sdk.AT_CheckWarning(rc);
%             end

            obj.t_start = util.text.time2str(datetime('now', 'TimeZone', 'UTC')); % get time from system clock
            
            for ii = 1:obj.batch_size
                
                buf = [];
                
                [rc, buf] = obs.cam.sdk.AT_WaitBuffer(obj.hndl, timeout*1000); % timeout in milliseconds! 
                
                if rc==0 % if we did not timeout or other error
                    
                    if ii<obj.batch_size
%                         if strcmp(obj.getTriggerModeHW, 'Software') % I think we can just give a SoftwareTrigger which is ignored in any other mode
                            [rc] = obs.cam.sdk.AT_Command(obj.hndl,'SoftwareTrigger'); obs.cam.sdk.AT_CheckWarning(rc);
                            [rc] = obs.cam.sdk.AT_QueueBuffer(obj.hndl, obj.imageSizeBytes); obs.cam.sdk.AT_CheckWarning(rc); 
%                         end
                    end
                    
                    [rc,buf2] = obs.cam.sdk.AT_ConvertMono16ToMatrix(buf, obj.AOIheight_c, obj.AOIwidth_c, obj.AOIstride); obs.cam.sdk.AT_CheckWarning(rc); % convert the raw buffer to a matrix we can use
                    
                    temp_images(:,:,ii) = buf2; % collect all images in the batch
                    
                    % Get timestamp and convert it into seconds
%                     [rc,ticks] = obs.cam.sdk.AT_GetTimeStamp(buf, obj.imageSizeBytes); obs.cam.sdk.AT_CheckWarning(rc);
%                     obj.timestamps(ii) = double(ticks)./double(obj.clockFreq);
                    
                else % if the acquisition got stuck... 
                    if rc==13
                        disp('ERROR 13: timeout...');
                    end
                    
                    obj.num_restarts = obj.num_restarts + 1;
                    disp(['Restarting acquisition... num_restart= ' num2str(obj.num_restarts)]);
                    obj.restart_sync; % stop the camera and restart the acquisition
                    
                end
                
            end
            
            obj.t_end = util.text.time2str(datetime('now', 'TimeZone', 'UTC')); % get time from system clock
            obj.t_end_stamp = obj.getTimestampHW; % get timestamps from camera clock to match to absolute time
            
            temp_images(:,:,ii+1:end) = []; % if batch ended early, get rid of empty frames
            obj.timestamps(ii+1:end) = []; % if batch ended early, get rid of empty timestamps
            
            obj.images = temp_images;
            
            obj.buffers.input(obj); % pass the timestamps, images, t_start, t_end, t_end_stamp to buffers for saving
            
        end
        
        function restart_sync(obj) % stop the camera and restart the acquisition (only using the synchronuous SDK tools)
            
            [rc] = obs.cam.sdk.AT_Command(obj.hndl,'AcquisitionStop'); obs.cam.sdk.AT_CheckWarning(rc);
            
            [rc] = obs.cam.sdk.AT_Flush(obj.hndl); obs.cam.sdk.AT_CheckWarning(rc); % flush existing hardware buffers
            
            for X = 1:10 % setup new hardware buffers
                [rc] = obs.cam.sdk.AT_QueueBuffer(obj.hndl, obj.imageSizeBytes); obs.cam.sdk.AT_CheckWarning(rc);
            end
                        
            [rc] = obs.cam.sdk.AT_Command(obj.hndl, 'AcquisitionStart'); obs.cam.sdk.AT_CheckWarning(rc);
            
%             if strcmp(obj.getTriggerModeHW, 'Software') % I think we can just give a SoftwareTrigger which is ignored in any other mode
                [rc] = obs.cam.sdk.AT_Command(obj.hndl,'SoftwareTrigger'); obs.cam.sdk.AT_CheckWarning(rc);
%             end
        end
        
        function stop(obj) % stops the acquisition and breaks out of any loops
            
            util.vec.mex_change(obj.mex_flag, 2, 1); % C++ interfaces through this flag

            obj.brake_bit = 1; % matlab loops are stopped with this value
            
        end
        
        function unlockMexFlag(obj) % used after an error in C++/mex interface
            
            for ii = 1:length(obj.mex_flag)
                util.vec.mex_change(obj.mex_flag, ii, 0);
            end
            
        end
        
    end
    
    methods % hardware access
        
        function val = getTemperatureHW(obj) % sensor temperature
           
            [rc, val] = obs.cam.sdk.AT_GetFloat(obj.hndl, 'SensorTemperature'); obs.cam.sdk.AT_CheckWarning(rc);
             
        end
        
        function val = getExpTimeHW(obj) % actual value as given to camera
           
            [rc, val] = obs.cam.sdk.AT_GetFloat(obj.hndl, 'ExposureTime'); obs.cam.sdk.AT_CheckWarning(rc);
                        
        end
        
        function [val_min, val_max] = getExpTimeLimitsHW(obj)
           
            [rc, val_min] = obs.cam.sdk.AT_GetFloatMin(obj.hndl, 'ExposureTime'); obs.cam.sdk.AT_CheckWarning(rc);
            [rc, val_max] = obs.cam.sdk.AT_GetFloatMax(obj.hndl, 'ExposureTime'); obs.cam.sdk.AT_CheckWarning(rc);
            
        end
        
        function setExpTimeHW(obj, val) % sets the frame rate to maximum (x0.95), considering the newly set exposure time
            
            if isempty(val)
                error('Must input a valid exposure time');
            end
            
            rc = obs.cam.sdk.AT_SetFloat(obj.hndl, 'ExposureTime', val); obs.cam.sdk.AT_CheckWarning(rc);
            
            if strcmp(obj.getTriggerModeHW, 'Internal')
                [~, f_max] = obj.getFrameRateLimitsHW;
                rc = obs.cam.sdk.AT_SetFloat(obj.hndl, 'FrameRate', f_max*0.95); obs.cam.sdk.AT_CheckWarning(rc);
            end
        end
        
        function val = getFrameRateHW(obj) % actual value as given to camera (relevnat only in TriggerMode="internal")
            
            if ~util.text.cs(obj.getTriggerModeHW, 'internal') % if we are not in "internal" trigger mode, there is no meaning to "frame rate"
                val = [];
            else
                [rc, val] = obs.cam.sdk.AT_GetFloat(obj.hndl, 'FrameRate'); obs.cam.sdk.AT_CheckWarning(rc);
            end
            
        end
        
        function [val_min, val_max] = getFrameRateLimitsHW(obj)
           
            [rc, val_min] = obs.cam.sdk.AT_GetFloatMin(obj.hndl, 'FrameRate'); obs.cam.sdk.AT_CheckWarning(rc);
            [rc, val_max] = obs.cam.sdk.AT_GetFloatMax(obj.hndl, 'FrameRate'); obs.cam.sdk.AT_CheckWarning(rc);
            
        end
        
        function setFrameRateHW(obj, val) % set to empty or NaN to ignore frame rate and set hardware TriggerMode to "internal"
            
            [~,f_max] = obj.getFrameRateLimitsHW;
            
            if val>=f_max % cannot set frame rate to higher than top frame rate (set by expT)
                str = sprintf('Cannot set frame rate to %g. Setting to %g instead', val, f_max.*0.95);
                disp(str);
                obj.log.input(str);
                val = f_max.*0.95;
            end 
            
            if isempty(val) || isnan(val)
                obj.setTriggerModeHW('software');
            else
                obj.setTriggerModeHW('internal');
                rc = obs.cam.sdk.AT_SetFloat(obj.hndl, 'FrameRate', val); obs.cam.sdk.AT_CheckWarning(rc);
            end
            
        end
        
        function val = getROI_HW(obj) % from hardware. ROI is defined as [top, left, height, width]. 
            
            [rc, t] = obs.cam.sdk.AT_GetInt(obj.hndl, 'AOITop'); obs.cam.sdk.AT_CheckWarning(rc);
            [rc, l] = obs.cam.sdk.AT_GetInt(obj.hndl, 'AOILeft'); obs.cam.sdk.AT_CheckWarning(rc);
            [rc, h] = obs.cam.sdk.AT_GetInt(obj.hndl, 'AOIHeight'); obs.cam.sdk.AT_CheckWarning(rc);
            [rc, w] = obs.cam.sdk.AT_GetInt(obj.hndl, 'AOIWidth'); obs.cam.sdk.AT_CheckWarning(rc);
            
            val = [l, t, w, h]; % flip axis from C to matlab
            
        end
        
        function val = is_zoomed_HW(obj) % check if ROI is smaller than max width/height
            
            % flip axis from C to matlab
            [rc, w] = obs.cam.sdk.AT_GetInt(obj.hndl, 'AOIHeight'); obs.cam.sdk.AT_CheckWarning(rc);
            [rc, h] = obs.cam.sdk.AT_GetInt(obj.hndl, 'AOIWidth'); obs.cam.sdk.AT_CheckWarning(rc);
            
            if h==obj.max_height && w==obj.max_width
                val = 0;
            else
                val = 1;
            end
            
        end
        
        function val = maxWidthHW(obj) % gets maximum width from hardware
            
            [rc, val] = obs.cam.sdk.AT_GetIntMax(obj.hndl, 'AOIHeight'); obs.cam.sdk.AT_CheckWarning(rc); % flip axis from C to matlab
            
        end
        
        function val = maxHeightHW(obj) % gets maximum height from hardware
            
            [rc, val] = obs.cam.sdk.AT_GetIntMax(obj.hndl, 'AOIWidth'); obs.cam.sdk.AT_CheckWarning(rc); % flip axis from C to matlab
            
        end
        
        function setROI_HW(obj, varargin) % input 4-element vector or 4 inputs. ROI is defined as [top, left, height, width]. 
            
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
            rc = obs.cam.sdk.AT_SetInt(obj.hndl, 'AOIHeight', w); obs.cam.sdk.AT_CheckError(rc);
            rc = obs.cam.sdk.AT_SetInt(obj.hndl, 'AOIWidth', h); obs.cam.sdk.AT_CheckError(rc);
            rc = obs.cam.sdk.AT_SetInt(obj.hndl, 'AOITop', l); obs.cam.sdk.AT_CheckError(rc);
            rc = obs.cam.sdk.AT_SetInt(obj.hndl, 'AOILeft', t); obs.cam.sdk.AT_CheckError(rc);
            
        end
        
        function val = getCycleModeHW(obj) % should always be "continuous" but maybe set to "fixed" (fixed number of frames)
            
            [rc, ind] = obs.cam.sdk.AT_GetEnumIndex(obj.hndl, 'CycleMode'); obs.cam.sdk.AT_CheckWarning(rc);
            [rc, val] = obs.cam.sdk.AT_GetEnumStringByIndex(obj.hndl, 'CycleMode', ind , 100); obs.cam.sdk.AT_CheckWarning(rc);
            
        end
        
        function setCycleModeHW(obj, mode) % should always be "continuous" but can also be set to "fixed" (fixed number of frames)
           
            import util.text.cs;
            
            if cs(mode, 'continuous')
                mode = 'Continuous';
            elseif cs(mode, 'fixed')
                mode = 'Fixed';
            else
                error(['unknown cycle mode: ' mode ]);
            end
            
            rc = obs.cam.sdk.AT_SetEnumString(obj.hndl, 'CycleMode', mode); obs.cam.sdk.AT_CheckWarning(rc);
            
        end 
        
        function val = getTriggerModeHW(obj) % when in "software", camera takes images as soon as it can. In "internal", camera maintains constant frame rate
            
            [rc, ind] = obs.cam.sdk.AT_GetEnumIndex(obj.hndl, 'TriggerMode'); obs.cam.sdk.AT_CheckWarning(rc);
            [rc, val] = obs.cam.sdk.AT_GetEnumStringByIndex(obj.hndl, 'TriggerMode', ind , 100); obs.cam.sdk.AT_CheckWarning(rc);
            
        end
        
        function setTriggerModeHW(obj, mode) % set to "software" when frame_rate is empty or NaN. Otherwise set to "internal". 
           
            import util.text.cs;
            
            if cs(mode, 'internal')
                mode = 'Internal';
            elseif cs(mode, 'software')
                mode = 'Software';
            else
                error(['unknown trigger mode: ' mode ]);
            end
            
            rc = obs.cam.sdk.AT_SetEnumString(obj.hndl, 'TriggerMode', mode); obs.cam.sdk.AT_CheckWarning(rc);
            
        end
        
        function val = getTimestampHW(obj) % get timestamp (in seconds), directly from camera clock
            
            [rc, ticks] = obs.cam.sdk.AT_GetInt(obj.hndl, 'TimestampClock'); obs.cam.sdk.AT_CheckWarning(rc)
            
            val = double(ticks)/obj.clockFreq;
            
        end
        
        function val = getGain(obj) % we don't have a way to get this from the camera just yet
            
            val = 0.6; % there must be a better way to do this??
            val = []; % I don't think 0.6 is reliable enough... 
            
        end
        
        function val = getCameraNameHW(obj) % get name from camera directly
            
            [rc, val] = obs.cam.sdk.AT_GetString(obj.hndl, 'CameraName', 100); obs.cam.sdk.AT_CheckWarning(rc);
            
        end
        
        function str_out = printout_HW(obj, full) % get hardware parameters and print them to string or to screen. If full=1 will print second line with more specific parameters. 
            
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
                
                [rc, baseline] = obs.cam.sdk.AT_GetInt(obj.hndl, 'Baseline'); obs.cam.sdk.AT_CheckWarning(rc);
                [rc, pixwidth] = obs.cam.sdk.AT_GetFloat(obj.hndl, 'PixelWidth'); obs.cam.sdk.AT_CheckWarning(rc);  
                [rc, readtime] = obs.cam.sdk.AT_GetFloat(obj.hndl, 'ReadoutTime'); obs.cam.sdk.AT_CheckWarning(rc);
                [rc, noisefilter] = obs.cam.sdk.AT_GetBool(obj.hndl, 'SpuriousNoiseFilter'); obs.cam.sdk.AT_CheckWarning(rc);
                [rc, blemish] = obs.cam.sdk.AT_GetBool(obj.hndl, 'StaticBlemishCorrection'); obs.cam.sdk.AT_CheckWarning(rc);
                [rc, cooling] = obs.cam.sdk.AT_GetBool(obj.hndl, 'SensorCooling'); obs.cam.sdk.AT_CheckWarning(rc);

                str2 = sprintf('baseline= %d | pixwidth= %f | readtime= %f | noise_filter= %d | blemish remove= %d | cooling= %d\n', baseline, pixwidth, readtime, noisefilter, blemish, cooling);
                
                str = [str str2];
                
            end
    
            if nargout>0
                str_out = str;
            else
                disp(str);
            end
            
        end
        
        function val = getShutterModeHW(obj)
            
            [rc, idx] = obs.cam.sdk.AT_GetEnumIndex(obj.hndl,'ElectronicShutteringMode'); obs.cam.sdk.AT_CheckWarning(rc);
            
            [rc, val] = obs.cam.sdk.AT_GetEnumStringByIndex(obj.hndl,'ElectronicShutteringMode', idx, 30); obs.cam.sdk.AT_CheckWarning(rc);

        end
        
        function setShutterModeHW(obj, val)
            
            if isempty(val) || ~ischar(val) 
                error('Input a shutter mode, either "global" or "rolling"'); 
            elseif util.text.cs(val, 'global')
                [rc] = obs.cam.sdk.AT_SetEnumString(obj.hndl,'ElectronicShutteringMode','Global'); obs.cam.sdk.AT_CheckWarning(rc);
            elseif util.text.cs(val, 'rolling')
                [rc] = obs.cam.sdk.AT_SetEnumString(obj.hndl,'ElectronicShutteringMode','Rolling'); obs.cam.sdk.AT_CheckWarning(rc);
            else
                error('Unknown shutter mode: "%s". Use "global" or "rolling"', val);
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function makeGUI(obj) % generate a GUI for this object
           
            if isempty(obj.gui)
                obj.gui = obs.cam.gui.AndorGUI(obj);                
            end
            
            obj.gui.make;
            
        end
        
        function show(obj, varargin) % show last image from last batch onto GUI axes (or current axes). 
            
            import util.text.cs;
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis'); % use this as axes to draw on. Default=[] means use GUI axes (or if it is invalid, use current axes)
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
        
        function showROI(obj, roi, ax) % show a rectangle of the ROI on top of a full frame image
            
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

