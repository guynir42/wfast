classdef ZylaControl < handle
   
    properties % outputs
        
        images@uint16;
        timestamps;
        timestamps_matlab;
        
    end
    
    properties % internal info
       
%         gui@obs.gui.ZylaGUI;
        gui;
        
        name = 'Zyla_5.5';
        
        num_restarts = 0;
        
        clockFreq;
        stride;        
        imageSize;
         
        readtime;
                
        height;
        width;
        
        gain = 0.6;
        
        matlab_time;
        
    end
    
    properties % inputs and switches
        
        hndl;
        
        product_type = 'RawFull';
        
        timeout = 1000; % in milliseconds...

        mode = 2;
        
        debug_bit = 1;
        
    end
    
    properties(Hidden=true)
        
        version = 1.01;
        
    end
    
    methods % constructor/initialize/shutdown
        
        function obj = ZylaControl
           
            if obj.debug_bit, fprintf('ZylaControl constructor v%4.2f\n', obj.version); end
            
            obj.initialize;
            
%             obj.gui = obs.gui.ZylaGUI(obj);
            
        end
        
        function reconnect(obj)
            
            obj.initialize;
            
        end
        
        function initialize(obj)
                                   
            try
                
                if obj.debug_bit, disp('initializing camera'); end
                addpath(getenv('ZYLA'));
                rc = AT_InitialiseLibrary; AT_CheckError(rc);
                
                [rc, obj.hndl] = AT_Open(0); AT_CheckError(rc);
                
                [rc] = AT_SetBool(obj.hndl, 'SensorCooling', 1); AT_CheckWarning(rc);        
              
                [rc, cooling] = AT_GetBool(obj.hndl, 'SensorCooling'); AT_CheckWarning(rc);

                if obj.debug_bit, disp(['Camera initialized, sensor cooling is ' num2str(cooling)]); end
                
                rc = AT_SetBool(obj.hndl, 'SpuriousNoiseFilter', 0); AT_CheckWarning(rc);
                rc = AT_SetBool(obj.hndl, 'StaticBlemishCorrection', 0); AT_CheckWarning(rc);
                
                [rc] = AT_SetEnumString(obj.hndl,'ElectronicShutteringMode','Rolling'); AT_CheckWarning(rc);
                               
                [rc] = AT_SetEnumString(obj.hndl,'SimplePreAmpGainControl','16-bit (low noise & high well capacity)'); AT_CheckWarning(rc);
                [rc] = AT_SetEnumString(obj.hndl,'PixelEncoding','Mono16'); AT_CheckWarning(rc);
                                
                % Enable Metadata
                [rc] = AT_SetBool(obj.hndl,'MetadataEnable',1); AT_CheckWarning(rc);
                [rc] = AT_SetBool(obj.hndl,'MetadataTimestamp',1); AT_CheckWarning(rc);
                
                if obj.mode==1 % testing the differences between these modes. in mode 2 need to do command: SoftwareTrigger
                    [rc] = AT_SetEnumString(obj.hndl,'CycleMode','Fixed'); AT_CheckWarning(rc);
                    [rc] = AT_SetEnumString(obj.hndl,'TriggerMode','Internal'); AT_CheckWarning(rc);
                elseif obj.mode==2                    
                    [rc] = AT_SetEnumString(obj.hndl,'CycleMode','Continuous'); AT_CheckWarning(rc);
                    [rc] = AT_SetEnumString(obj.hndl,'TriggerMode','Software'); AT_CheckWarning(rc);
                end
                
                % get some info from camera
                [rc, baseline] = AT_GetInt(obj.hndl, 'Baseline'); AT_CheckWarning(rc);
                [rc, pixwidth] = AT_GetFloat(obj.hndl, 'PixelWidth'); AT_CheckWarning(rc);  
                [rc, obj.readtime] = AT_GetFloat(obj.hndl, 'ReadoutTime'); AT_CheckWarning(rc);
                [rc, noisefilter] = AT_GetBool(obj.hndl, 'SpuriousNoiseFilter'); AT_CheckWarning(rc);
                [rc, blemish] = AT_GetBool(obj.hndl, 'StaticBlemishCorrection'); AT_CheckWarning(rc);
                
                if obj.debug_bit
                    disp(['the baseline is ' num2str(baseline) ' | pixel width= ' num2str(pixwidth)...
                    ' | readout time= ' num2str(obj.readtime) ' | noise reduction= ' num2str(noisefilter),... 
                    ' | blemish removal= ' num2str(blemish)]);
                end
                
            catch ME
                warning(getReport(ME));
                obj.hndl = 0;                
            end
            
        end
        
        function shutdown(obj)
           
            obj.finishup;
            % other things like unload the library?

            rc = AT_Close(obj.hndl); AT_CheckWarning(rc);

            rc = AT_FinaliseLibrary; AT_CheckWarning(rc);
            
        end
        
    end
    
    methods % reset and clear
        
    end
    
    methods % getters
        
        function val = getTemperature(obj)
           
            val = 1000;
%             [rc, val] = AT_GetFloat(obj.hndl, 'Temperature'); AT_CheckWarning(rc);
             
        end
        
        function val = getExpTime(obj)
           
            [rc, val] = AT_GetFloat(obj.hndl, 'ExposureTime'); AT_CheckWarning(rc);
                        
        end
        
        function [val_min, val_max] = getExpTimeLimits(obj)
           
            [rc, val_min] = AT_GetFloatMin(obj.hndl, 'ExposureTime'); AT_CheckWarning(rc);
            [rc, val_max] = AT_GetFloatMax(obj.hndl, 'ExposureTime'); AT_CheckWarning(rc);
            
        end
        
        function val = getExpTimeMax(obj)
                        
            [rc, val] = AT_GetFloatMax(obj.hndl, 'ExposureTime'); AT_CheckWarning(rc);

        end
        
        function val = getExpTimeMin(obj)
                        
            [rc, val] = AT_GetFloatMin(obj.hndl, 'ExposureTime'); AT_CheckWarning(rc);

        end
        
        function val = getFrameRate(obj)
           
            [rc, val] = AT_GetFloat(obj.hndl, 'FrameRate'); AT_CheckWarning(rc);
                        
        end
        
        function [val_min, val_max] = getFrameRateLimits(obj)
           
            [rc, val_min] = AT_GetFloatMin(obj.hndl, 'FrameRate'); AT_CheckWarning(rc);
            [rc, val_max] = AT_GetFloatMax(obj.hndl, 'FrameRate'); AT_CheckWarning(rc);
            
        end
        
        function val = getFrameRateMax(obj)
                   
            [rc, val] = AT_GetFloatMax(obj.hndl, 'FrameRate'); AT_CheckWarning(rc);
     
        end
        
        function val = getFrameRateMin(obj)
           
            [rc, val] = AT_GetFloatMin(obj.hndl, 'FrameRate'); AT_CheckWarning(rc);

        end
        
        function val = getHeight(obj)
                            
            [rc, val] = AT_GetInt(obj.hndl, 'AOIWidth'); AT_CheckWarning(rc);
            
        end
        
        function val = getWidth(obj)
            
            [rc, val] = AT_GetInt(obj.hndl, 'AOIHeight'); AT_CheckWarning(rc);
            
        end
        
        function val = getLeft(obj)
            
            [rc, val] = AT_GetInt(obj.hndl, 'AOILeft'); AT_CheckWarning(rc);
            
        end
        
        function val = getTop(obj)
            
            [rc, val] = AT_GetInt(obj.hndl, 'AOITop'); AT_CheckWarning(rc);
            
        end
        
        function val = maxWidth(obj)
            
            val = 2160;
            
        end
        
        function val = maxHeight(obj)
            
            val = 2560;
            
        end
                
        function val = getCycleMode(obj)
            
            [rc, ind] = AT_GetEnumIndex(obj.hndl, 'CycleMode'); AT_CheckWarning(rc);
            [rc, val] = AT_GetEnumStringByIndex(obj.hndl, 'CycleMode', ind , 100); AT_CheckWarning(rc);
            
        end
        
        function val = getTriggerMode(obj)
            
            [rc, ind] = AT_GetEnumIndex(obj.hndl, 'TriggerMode'); AT_CheckWarning(rc);
            [rc, val] = AT_GetEnumStringByIndex(obj.hndl, 'TriggerMode', ind , 100); AT_CheckWarning(rc);
            
        end
        
        function val = getTimestamp(obj)
            
            [rc, obj.clockFreq] = AT_GetInt(obj.hndl, 'TimestampClockFrequency'); AT_CheckWarning(rc)
                        
            [rc, ticks] = AT_GetInt(obj.hndl, 'TimestampClock'); AT_CheckWarning(rc)
            
            val = double(ticks)/double(obj.clockFreq);
            
        end
        
    end
    
    methods % setters
        
        function setExpTime(obj, T)
            
            assert(~isempty(T), 'cannot accept empty values!');
            
            try 
            
                rc = AT_SetFloat(obj.hndl, 'ExposureTime', T); AT_CheckWarning(rc);
                
                obj.timeout = max([1000 obj.getExpTime.*1200]);
            
            catch ME
                if obj.debug_bit, disp(ME.getReport); end
            end
            
        end
        
        function setFrameRate(obj, f)
                        
            assert(~isempty(f), 'cannot accept empty values!');
            
            try 
                rc = AT_SetFloat(obj.hndl, 'FrameRate', f);
                if rc==6
                    disp(['FrameRate ' num2str(f) ' Hz out of range for these settings. Using ' num2str(obj.getFrameRate)])
                end
            catch ME
                if obj.debug_bit, disp(ME.getReport); end
            end
        end
        
        function setHeight(obj, val)
           
            rc = AT_SetInt(obj.hndl, 'AOIWidth', val); AT_CheckWarning(rc);
            
        end
        
        function setWidth(obj, val)
            
            rc = AT_SetInt(obj.hndl, 'AOIHeight', val); AT_CheckWarning(rc);
            
        end
        
        function setLeft(obj, val)
            
            rc = AT_SetInt(obj.hndl, 'AOITop', val); AT_CheckWarning(rc);
            
        end
        
        function setTop(obj, val)
            
            rc = AT_SetInt(obj.hndl, 'AOILeft', val); AT_CheckWarning(rc);
            
        end
        
        function setCenterX(obj, val)
            
            width = obj.getWidth;
            left = val - floor(width/2) + 1;
            
            assert(left>1, ['cannot set center_x to ' num2str(val) ' the width ' num2str(width) ' is too big!']);
            
            obj.setLeft(left);
            
        end
        
        function setCenterY(obj, val)
            
            height = obj.getHeight;
            top = val - floor(height/2) + 1;
            
            assert(top>1, ['cannot set center_y to ' num2str(val) ' the height ' num2str(height) ' is too big!']);
            
            obj.setTop(top);
            
        end
        
        function zoom(obj, width, height, center_x, center_y)
           
            if nargin<2 || isempty(width)
                width = 512;
            end
            
            if nargin<3 || isempty(height)
                height = width;
            end
            
            if nargin<4 || isempty(center_x)
                center_x = floor((obj.maxWidth)/2)+1;
            end
            
            if nargin<5 || isempty(center_y)
                center_y = floor((obj.maxHeight)/2)+1;
            end
            
            obj.setWidth(width);
            obj.setHeight(height);
            obj.setCenterX(center_x);
            obj.setCenterY(center_y);
            
        end
        
        function full_frame(obj)
           
            obj.setLeft(1);
            obj.setTop(1);
            obj.setWidth(obj.maxWidth);
            obj.setHeight(obj.maxHeight);
            
        end
        
        function setCycleMode(obj, mode)
           
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
        
        function setTriggerMode(obj, mode)
           
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
        
    end
    
    methods % actions
        
        function startup(obj, owner, input, buffer_wheel)
            
            for ii = 1:length(buffer_wheel.buf)
                util.vec.mex_change(buffer_wheel.buf(ii).mex_flag_record, 1, 1); % lock all the buffers for recording, they will unlock once the camera fills them... 
            end
            
            obj.captureMex(owner, input.num_batches, input.batch_size, buffer_wheel); % start rolling camera!
            
%             % setup buffers
%             [rc] = AT_Flush(obj.hndl); AT_CheckWarning(rc);
%             
%             [rc, obj.imageSize] = AT_GetInt(obj.hndl, 'ImageSizeBytes'); AT_CheckWarning(rc)
%                         
%             [rc, obj.height] = AT_GetInt(obj.hndl, 'AOIHeight'); AT_CheckWarning(rc)
%                         
%             [rc, obj.width] = AT_GetInt(obj.hndl, 'AOIWidth'); AT_CheckWarning(rc)
%             
%             [rc, obj.stride] = AT_GetInt(obj.hndl, 'AOIStride'); AT_CheckWarning(rc)
%             
%             [rc, obj.clockFreq] = AT_GetInt(obj.hndl, 'TimestampClockFrequency'); AT_CheckWarning(rc)
%             
% %             obj.timeout = 300;
%             
%             for X = 1:10
%                 [rc] = AT_QueueBuffer(obj.hndl, obj.imageSize); AT_CheckWarning(rc);
%             end
%             
%             % reset the time stamp
%             [rc] = AT_Command(obj.hndl, 'TimestampClockReset'); AT_CheckWarning(rc);
%             obj.matlab_time = clock;
%             
%             [rc] = AT_Command(obj.hndl, 'AcquisitionStart'); AT_CheckWarning(rc);
%             
%             obj.num_restarts = 0;
            
        end
        
        function captureMex(~, owner, num_batches, batch_size, buffer_wheel)
                        
            if ~isa(buffer_wheel, 'file.BufferWheel')
                error('Must give a file.BufferWheel object as input 4. Was given a %s...', class(buffer_wheel));
            end
            
            if owner.use_error_log
                owner.error_log = zeros(num_batches,1); 
            else
                owner.error_log = [];
            end
            
            obs.cam.mex.capture(owner, owner.mex_flag, buffer_wheel.buf, buffer_wheel.index_rec_vec, num_batches, batch_size, owner.error_log);
            
        end
        
        function finishup(obj, owner)
            
%             [rc] = AT_Command(obj.hndl,'AcquisitionStop'); AT_CheckWarning(rc);
            
%             [rc] = AT_Flush(obj.hndl); AT_CheckWarning(rc);
            
        end
        
        function restart(obj)
            
            [rc] = AT_Command(obj.hndl,'AcquisitionStop'); AT_CheckWarning(rc);
            
            [rc] = AT_Flush(obj.hndl); AT_CheckWarning(rc);
            
            for X = 1:10
                [rc] = AT_QueueBuffer(obj.hndl, obj.imageSize); AT_CheckWarning(rc);
            end
                        
            [rc] = AT_Command(obj.hndl, 'AcquisitionStart'); AT_CheckWarning(rc);
            
            if strcmp(obj.getTriggerMode, 'Software')
                [rc] = AT_Command(obj.hndl,'SoftwareTrigger'); AT_CheckWarning(rc);
            end
        end
        
        function batch(obj, buf, owner)
            
            % you can add some Zyla specific processing but it will all
            % happen inside the mex file anyway. 
            
        end
        
        function batchX(obj, Nframes, top_level_obj) % old code, no longer in use since it all happens in the mex file
                       
            if nargin<3 || isempty(top_level_obj)
                top_level_obj = [];
            end
            
%             obj.images = zeros(obj.width, obj.height, Nframes, 'uint16'); % notice the axes are inverted because this is C vs. Matlab! 
            temp_images = zeros(obj.width, obj.height, Nframes, 'uint16'); % notice the axes are inverted because this is C vs. Matlab! 
            obj.timestamps = zeros(Nframes,1);
            obj.timestamps_matlab = zeros(Nframes,1);
            
            if strcmp(obj.getCycleMode, 'Fixed')
                [rc] = AT_SetInt(obj.hndl,'FrameCount',Nframes); AT_CheckWarning(rc);
            end
            
            counter = 0;
            
            [rc] = AT_QueueBuffer(obj.hndl,obj.imageSize); AT_CheckWarning(rc);

            if strcmp(obj.getTriggerMode, 'Software')
                [rc] = AT_Command(obj.hndl,'SoftwareTrigger'); AT_CheckWarning(rc);
            end
            
            for ii = 1:Nframes
                
%                 buf = zeros(obj.width, obj.height);
                buf = [];
                
                [rc, buf] = AT_WaitBuffer(obj.hndl,obj.timeout); % AT_CheckWarning(rc);
                
                if rc==0 % if we did not timeout or other error

                    if strcmp(obj.getTriggerMode, 'Software')
                        [rc] = AT_Command(obj.hndl,'SoftwareTrigger'); AT_CheckWarning(rc);
                    end
                    
                    counter = counter+1;
                    
                    [rc,buf2] = AT_ConvertMono16ToMatrix(buf, obj.height, obj.width, obj.stride); AT_CheckWarning(rc);
                    
%                     if ii<Nframes, 
                        [rc] = AT_QueueBuffer(obj.hndl,obj.imageSize); AT_CheckWarning(rc); 
%                     end
                    
                    temp_images(:,:,counter) = buf2;
                    
                    %Get timestamp and convert it into seconds
                    [rc,ticks] = AT_GetTimeStamp(buf, obj.imageSize); AT_CheckWarning(rc);
                    obj.timestamps(counter) = double(ticks)./double(obj.clockFreq);
                    obj.timestamps_matlab(counter) = etime(clock, obj.matlab_time);
%                     disp(['t= ' num2str(obj.timestamps(counter))]);

                else % if the acquisition got stuck... 
                    if rc==13
                        disp('ERROR 13: timeout...');
                    end
                    obj.num_restarts = obj.num_restarts + 1;
                    disp(['Restarting acquisition... num_restart= ' num2str(obj.num_restarts)]);
                    obj.restart;
                end
                                
                if ~isempty(top_level_obj)
%                     top_level_obj.updateGUI;
                    top_level_obj.live_counter = top_level_obj.live_counter + 1;
                    top_level_obj.gui.updateLive;
                    if rc==0
%                         disp('calling show');
                        top_level_obj.show(buf2); % also calls draw now after putting the image in place...
                    end
                    if top_level_obj.brake_bit || ~top_level_obj.gui.check
                        break;
                    end
                    
                end
                                
            end
            
            temp_images(:,:,counter+1:end) = [];
            obj.timestamps(counter+1:end) = [];
            obj.timestamps(counter+1:end) = [];
            
            obj.images = temp_images;
            
            if 0
                
                for ii = 1:length(obj.timestamps)
                    fprintf('t= %10.3f | tm= %10.3f\n', obj.timestamps(ii), obj.timestamps_matlab(ii));
                end
                
            end
              
        end
        
    end
    
end