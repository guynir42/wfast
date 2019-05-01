classdef Andor < file.AstroData

    properties(Transient=true)
        
        gui;
        audio@util.sys.AudioControl;
        
    end
    
    properties % objects
        
        pars@head.Parameters; % parameters used in the observations
        buffers@file.BufferWheel;
        focuser;
        hndl; % pointer to the Andor SDK handle of the camera
        
    end
    
    properties % inputs/outputs
        
        mean_frame_rate;
        
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
        ROI_size; % scalar size of the ROI 
        ROI_center_xy; % if empty, choose the center of the field
        ROI; % [left, top, width, height] as given to camera
        
        use_show_flipped = 0; % flip the display for North is up after meridian flip
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        mex_flag = [0 0 0 0]; % for starting and stopping the camera...
        
        mode_list = {'science', 'dark', 'flat'};
        mode_index = 1;
        
        default_frame_rate; 
        default_expT; 
        default_expT_deep;
        default_batch_size; 
        default_num_batches; 
        
        now_recording = 0;
        now_live = 0;
        previous_folder = ''; % this is to send a warning when trying to write twice to the same folder
        use_prev_folder_warning = 1; % ask if the user is sure about recording twice to the same folder
        
        batch_counter = 0; % how many batches were recorded in this run
        clip_limit = 16000; % any pixels above this value will display a "clipping" warning... 
        is_clipping = 0; % update this manually in record/live view
        
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
                obj.log.input('Connecting to Andor Camera!');
                
                try 

                    obj.pars = head.Parameters;

                    obj.setupBuffers;
                    
                    obj.setupFocuser;
                    

                catch ME
                    obj.log.error(ME.getReport);
                    rethrow(ME);
                end
                
            end
            
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
        
    end
    
    methods % reset/clear
        
        function connect(obj)
            
        end
        
        function disconnect(obj)
            
        end
        
        function reconnect(obj)
            
            obj.disconnect;
            obj.connect;
            
        end
        
        function reset(obj) % get ready for a new run

            obj.buffers.reset; 
            util.vec.mex_change(obj.mex_flag, 4, 0); % reset the counter on the mex_flag
            
        end
        
        function clear(obj) % removes outputs and intermediary data
            
            % maybe clear the buffers??
            
        end
        
    end
    
    methods % getters
        
        function val = is_capturing(obj)
            
            val = obj.mex_flag(1);
            
        end
        
        function val = is_finished(obj)
            
            val = obj.brake_bit;
            
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
        
    end
    
    methods % actions
        
        function run(obj, varargin)
            
        end
        
        function batch(obj, varargin)
            
        end
        
    end
    
    methods % internal functions 
        
        function input = makeInputVars(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('num_batches', obj.num_batches, 'Nbatches');
            input.input_var('batch_size', obj.batch_size);
            input.input_var('frame_rate', obj.frame_rate, 'frequency');
            input.input_var('expT', obj.expT, 'exposure time', 'expT');
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
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
            
            
        end
        
    end    
    
end

