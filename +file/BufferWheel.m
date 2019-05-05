classdef BufferWheel < file.AstroData
% Usage: obj = BufferWheel(num_buffers=5, pars=[])
% Holds several (e.g., 5) Buffer structs and rotates them as they get filled
% and emptied to file. 
%
% Uses index (and index_rec) to keep track of which buffer to use when
% recording or reading and writing to file. Each buffer struct can have all 
% the inputs and outputs (e.g. images, fluxes, psfs, etc.) as defined in 
% file.AstroData class.
%
% Model 1: input data from camera/sim/pipeline using input method, then
% immediately use save method which also brings up new buffer. If you try
% to input to a buffer that is still saving, you get an error. 
% (this is used by Acquisition and Deflator classes to save cutouts etc.).
% 
% Model 2: Give the camera all the buffers and have it track which one it
% is writing to using "index_rec". If buffer is still saving or analysing, 
% camera is delayed. When camera finishes, it gets the next buffer in line. 
% Run analysis on any buffers released from camera (if analysis slows down,
% you may need to catch up or skip buffers). Keep track of which buffer is
% under analysis at any time using "index". Use "save" to save current 
% buffer, e.g., when there is a trigger in the buffer you are now analyzing. 
%
% OBJECTS:
%   -buf: vector of Buffer structs. 
%   -pars: link back to parent object.
%   -gui: Optional GUI object.
%   -this_buf, next_buf, prev_buf: getters from "buf", based on index.
%    (NOTE: these shortcuts can ONLY be used for READING data from buffers). 
%
% SWITCHES & CONTROLS:
%   -product_type: string telling which kind of output (raw? cutouts?
%    fluxes? dark/flat? etc.). Will print this in the filename. 
%   -base_dir: set to [] to use getenv('DATA_TEMP') or getenv('DATA') 
%   -date_dir: set to [] to use this night's date dir. 
%   -target_dir: set to [] to use "pars.target_name" as target_dir. 
%     *** "full_path" is: base_dir/date_dir/target_dir. ***
%   -dir_extension: overrides the automatic extension to the dir name. 
%   -use_dir_type: automatically adds "product_type" to the dir name. 
%   -use_overwrite: silently delete existing files (default is true).
%   -use_write_pars: try to save the "pars" object into file.
%   -use_deflate: if non-zero, use deflation (can be anything from 1 to 10)
%   -use_async: write in separate thread (mex only!). 
%   -use_mex: use the mex function to write to file (HDF5 only, for now). 
%   -use_disable_mex_async_deflate: when this is enabled, doesn't let you
%    use deflate when using mex and async (avoids crashes). 
%   -chunk: used for deflation. Default is [64 64].
%
% NOTE: product types:
%
% There are five (main) file types: 
% 1) Full: full-frame data. Usually Raw data is saved straight from camera 
%          in case of trigger, or when observing slow mode or burst mode. 
% 3) Dark: Same as "Full", only for dark frames. These are always full raw. 
% 4) Flat: Same as "Full", only for flat frames. These are always full raw. 
% 5) Cutouts: including cutouts, positions, stacks, lightcurves, etc. This 
%             is the main product type for normal W-FAST operations. 
% 
%   *Additional qualifiers may include "Sim" for simulated images, "cal" or 
%    "proc" for calibrated/processed full or cutout frames. 
%
% NOTE: mex_flags are used to tell matlab that a certain buffer is still 
%       working on something and shouldn't be altered. Three things can 
%       occupy a buffer: 
%       * recording (i.e., by the camera in a separate thread) 
%       * reading by the analysis program 
%       * writing to disk (i.e., by a separate thread using mex-write). 
%       Each of these flags is a 2-element vector that is passed to mex
%       functions that can read and change the value under the hood. That's
%       why you must not change these vectors directly (then there will be
%       a different vector for matlab and for the mex thread). 
%       The first element of each vector is the lock (0- open, 1-locked).
%       The second is the number of seconds spent waiting for the lock. 

    properties(Transient=true)
        
        gui@file.gui.BufferGUI;
        
    end

    properties % objects
        
        pars@head.Parameters; % link back to owner object's pars
        buf@struct; % struct array with the actual data saved in each element's fields
        
    end
    
    properties % switches
        
        index = 1; % which "buf" is now active for reading/saving 
        product_type = 'Raw'; % can be Raw, ROI, Cutouts (including stacks), Dark (always raw), Dlat (always raw). Can append additional strings like "Sim" or "Cal" or "Proc"
        dir_extension = ''; % overrides any automatic directory extensions... 
        use_dir_types = 0; % automatically add a dir_extension equal to "type"
        
        % write parameters
        use_save_raw_images = 1; % if set to 0, will not save to disk the full frame raw images! 
        use_overwrite = 1; % delete existing files (this should not happen in normal operations because of unique filenames) 
        use_write_pars = 1; % write pars object to each file
        use_deflate = 0; % also can specify how much deflation you want, between 1 to 10 (more than 1 is usually slow and doesn't change much)
        use_async = 1; % write files in a different thread (mex only)
        use_mex = 1; % use mex interface to write files (use 0 as backup if mex fails or if you want other file types)
        chunk = 64; % chunk size (for deflate only)
        
        file_type = 'hdf5'; % support is added for "fits", "mat", and "tiff" in use_mex==0
                
        serial = 1; % keeps track of the number of files saved since object is initialized (should this go to 1 in "reset"?)
        
        save_time = 0; % time to save images, since run started (seconds)
        num_saved_batches = 0; % number of batches saved since run started. 
        
        timeout = 5; % defualt timeout, in seconds
        
        debug_bit = 1;
                
    end
    
    properties(Dependent=true)
        
        index_rec; % depends on hidden value that is updated internally by camera mex file
            
        this_buf; % link to current buffer (this is a copy of a struct, so do not change its properties)
        prev_buf; % link to previous buffer (this is a copy of a struct, so do not change its properties)
        next_buf; % link to next buffer (this is a copy of a struct, so do not change its properties)
               
        im_size; % depends on size of images that were given
        
        base_dir; % where is the data stored (default is getenv('DATA_TEMP') or getenev('DATA'))
        date_dir; % YYYY-MM-DD sub-directory
        target_dir; % name of target/object/run as subfolder below date_dir
        directory; % full folder name is base_dir/date_dir/target_dir
        filename; % filename composed of naming convention, or overriden by user
        prev_name; % full file name (including path) to last file that was saved (use h5disp on this to inspect file)
                
    end
    
    properties(Hidden=true)
        
        index_rec_vec = [0 0]; % do not change this vector (interacts under the hood with camera mex file)
        pars_struct_cell = {}; % transform the "pars" object into a cell of structs before saving it. 
        
        % naming convention override (user can override any of these)
        base_dir_override; 
        date_dir_override;
        target_dir_override;
        directory_override;
        filename_override;
        
        % what to call the directory for darks/flats (used when product_type is dark/flat)
        dark_dir_name = 'dark';
        flat_dir_name = 'flat'; 
        
        default_product_type; 
        default_dir_extension;
        default_use_dir_types;
        default_project = 'WFAST';
        default_camera = 'Zyla';
        
        default_use_overwrite;
        default_use_write_pars;
        default_use_deflate;
        default_use_async;
        default_use_mex;
        default_chunk;
        default_file_type;
        
        use_disable_mex_async_deflate = 1; % if you try to defalte using mex file and async write, it will silently cancel deflation... if you override this, Matlab may crash...
        
        camera_mex_flag; % used for interaction with camera mex. Do not change this vector!
        % vector elements: 
        %      (1) is the camera running 
        %      (2) error flag from camera 
        %      (3) how many milliseconds is camera waiting for write
        
        version = 1.02;
        
    end
    
    methods % constructor
        
        function obj = BufferWheel(num_buffers, pars)
        % Constructor for BufferWheel. Default number of buffers is 5. 
        % Second argument links the Parameters object with owner. Default
        % is to create a new pars object. 
        
            if nargin<1 || isempty(num_buffers)
                num_buffers = 5;
            end
            
            if nargin<2 || isempty(pars) || ~isa(pars, 'head.Parameters')
                pars = head.Parameters.empty;
            end
            
            if obj.debug_bit
                fprintf('BufferWheel constructor v%4.2f with %d buffers. ', obj.version, num_buffers);
                if isempty(pars) 
                    fprintf('Creating new pars object\n');
                else
                    fprintf('Using input pars.\n');
                end
            end
            
            if isempty(pars)
                obj.pars = head.Parameters;
            else
                obj.pars = pars; % link to pars handle as given
            end
            
            for ii = 1:num_buffers
                obj.addBuffer; % this generates each buffer with all the data it needs (based on file.AstroData class)
            end
            
            util.oop.save_defaults(obj); % save default value of property "X" to "default_X" if it exists
            
        end
        
    end
    
    methods % resetters
        
        function initialize(obj) % return to class defaults
                      
            util.oop.load_defaults(obj);
            obj.reset;
                        
        end
        
        function reset(obj) % prepare object for new run
                 
            reset@file.AstroData(obj);
            
            for ii = 1:length(obj.buf)
                obj.clearBuf(ii);
            end
            
            obj.index = 1;
            obj.index_rec = 1;
            
            obj.save_time = 0; 
            obj.num_saved_batches = 0; 
        
            obj.unlockMexFlag;
            
            obj.pars_struct_cell = {};
            
        end
                
        function clear(obj) % prepare object for new batch
        
            clear@file.AstroData(obj);
%             obj.clearBuf(obj.index);
            
        end
        
        function clearBuf(obj, idx) % clear data inside a single buffer

%             obj.buf(idx).images = [];
% 
%             obj.buf(idx).cutouts = [];
%             obj.buf(idx).positions = [];
% 
%             obj.buf(idx).cutouts_bg = [];
%             obj.buf(idx).positions_bg = [];
% 
%             obj.buf(idx).stack = [];
%             obj.buf(idx).num_sum = [];
% 
%             obj.buf(idx).timestamps = [];
%             obj.buf(idx).t_start = [];
%             obj.buf(idx).t_end = [];
%             obj.buf(idx).t_vec = [];
% 
%             obj.buf(idx).psfs = [];
%             obj.buf(idx).psf_sampling = [];

            list = properties(file.AstroData);
            for ii = 1:length(list)
                obj.buf(idx).(list{ii}) = [];
            end
            
            obj.buf(idx).latest_filename = [];

            util.vec.mex_change(obj.buf(idx).mex_flag_record, 2, 0);
            util.vec.mex_change(obj.buf(idx).mex_flag_read, 2, 0);
            util.vec.mex_change(obj.buf(idx).mex_flag_write, 2, 0);
            
        end
        
        function clearImages(obj) % get rid of full frame images only (to be depricated!)
            
            obj.images = [];
            for ii = 1:length(obj.buf)
                obj.buf(ii).images = [];
            end
            
        end
        
    end 
    
    methods % getters

        function val = isWriting(obj, buf) % check if mex process for saving to disk is still running on "buf" (default is "this_buf")
            
            if nargin<2 || isempty(buf)
                buf = obj.this_buf;
            end
            
            val = buf.mex_flag_write(1)==1;
            
%             flag = buf.mex_flag_write;
%             val = flag(1)==1; % && flag(2)~=1;
            
        end
        
        function val = isRecording(obj, buf) % check if mex process for recording in camera is still running on "buf" (default is "this_buf")
            
            if nargin<2 || isempty(buf)
                buf = obj.this_buf;
            end
            
            val = buf.mex_flag_record(1)==1;
            
%             flag = buf.mex_flag_record;
%             val = flag(1)==1; %&& flag(2)~=1;
%             val = val || ~buf.mex_flag_read(2); % add second check to see if there is any data in the buffer at all (if not, consider it still under recording...)
            
        end
        
        function val = isReading(obj, buf) % check if current analysis is still running on "buf" (default is "this_buf")
            
            if nargin<2 || isempty(buf)
                buf = obj.this_buf;
            end
            
            val = buf.mex_flag_read(1)==1;
            
%             flag = buf.mex_flag_read;
%             val = flag(1)==1; %&& flag(2)~=1;
            
        end
        
        function val = checkReadability(obj) % check if we can get information from buffers
            
            val = ~obj.isRecording;
            
        end
        
        function val = checkInputability(obj) % check if we can set information to buffers
            
            val = ~obj.isRecording && ~obj.isReading && ~obj.isWriting;
            
        end
        
        function val = checkWritablility(obj) % check if we can send buffer to write to disk
            
            val = ~obj.isRecording && ~obj.isWriting;
            
        end
          
        function val = get.index_rec(obj) % convert from mex-compatible vector to matlab index
            
            val = obj.index_rec_vec(1) + 1;
            
        end
        
        function val = mod(obj, index) % cycle index over the number of buffers in the wheel
            
            val = mod(index-1, length(obj.buf))+1;
            
        end
        
        function val = get.this_buf(obj) % buffer being read from right now (as set by "index")
            
            val = obj.buf(obj.mod(obj.index));
            
        end
        
        function val = get.next_buf(obj) % buffer for index+1
           
            val = obj.buf(obj.mod(obj.index+1));
            
        end
        
        function val = get.prev_buf(obj) % buffer for index-1
           
            val = obj.buf(obj.mod(obj.index-1));
            
        end
        
        function val = get.im_size(obj) % size of images in current buffer (if images is empty, use PSF size)
           
            if ~isempty(obj.this_buf.images)
                val = util.vec.imsize(obj.this_buf.images);
            elseif ~isempty(obj.this_buf.psfs)
                val = util.vec.imsize(obj.this_buf.psfs);
            else
                val = [];
            end
            
        end
        
        function val = get.use_deflate(obj) % can be overriden by "use_disable_mex_async_deflate" 
            
            val = obj.use_deflate;
            
            if obj.use_disable_mex_async_deflate && obj.use_mex && obj.use_async
                val = 0;
            end
            
        end                
        
        function val = get.base_dir(obj) % check for override. Default: use getenv('DATA_TEMP') on SSD or getenv('DATA') on storage or PWD if no environmentals are found 
            
            if ~isempty(obj.base_dir_override)
                val = obj.base_dir_override;
            elseif ~isempty(getenv('DATA_TEMP'))
                val = getenv('DATA_TEMP');
            elseif ~isempty(getenv('DATA'))
                val = getenv('DATA');
            else
                val = pwd;    
            end
            
        end
        
        function date = current_datetime(obj) % create a datetime object from camera's "t_end" datestr or from current time
            
            if ~isempty(obj.this_buf.t_end)
                date = util.text.str2time(obj.this_buf.t_end); % best time estimate is when the file has finished recording (matched from camera)
            else
                date = datetime('now', 'TimeZone', 'UTC'); % just use what ever time is right now...
            end
            
        end
        
        function val = get.date_dir(obj) % check for override. Default: get YYYY-MM-DD format dirname. Will change day at noon UTC! 
                       
            if ~isempty(obj.date_dir_override)
                val = obj.date_dir_override;                
            else
            
                val = util.sys.date_dir; 
                
            end
            
        end
        
        function val = get.target_dir(obj) % check for override. Default: get subfolder name from pars object (uses <target_name> or "dark" or "flat")
            
            if ~isempty(obj.target_dir_override)
                val = obj.target_dir_override;                
            elseif ~isempty(obj.pars) && ~isempty(obj.pars.folder_name)
                val = obj.pars.folder_name;
            else
                val = 'run1';
            end
            
            if isnumeric(val)
                val = num2str(val);
            end
            
            if ~isempty(obj.dir_extension)
                val = [val '_' obj.dir_extension];
            elseif obj.use_dir_types && ~isempty(obj.product_type)
                val = [val '_' obj.product_type];
            end
            
        end
        
        function val = get.directory(obj) % check for override. Default: directory is made from <base_dir>/<date_dir>/<target_dir>
            
            if ~isempty(obj.directory_override)
                val = obj.directory_override;
            else
                val = obj.makeFullpath;
            end
            
        end
        
        function val = get.filename(obj) %  check for override. Default: naming convention in "makeFilename"
            
            if ~isempty(obj.filename_override)
                val = obj.filename_override;
            else
                val = obj.makeFilename;
            end
            
        end
                
        function val = makeFullpath(obj) % directory is made from <base_dir>/<date_dir>/<target_dir>
            
            d = obj.target_dir;
            
            if util.text.cs(obj.product_type, 'dark')
                d = obj.dark_dir_name;
            elseif util.text.cs(obj.product_type, 'flat')
                d = obj.flat_dir_name;
            end
            
            val = fullfile(obj.base_dir, obj.date_dir, d);
            
        end
        
        function val = makeFilename(obj) % apply filename convention 
           
            import util.text.cs;
            
            date = obj.current_datetime; % datetime object of now (or buffer's t_end, if it exists)
            
            time_str = datestr(date, 'yyyymmdd-HHMMSS-FFF'); % convert using matlab's own datetime functions 
            
            if ~isempty(obj.pars)
                project = obj.pars.project;
                camera = obj.pars.cam_name;
            else
                project = obj.default_project;
                camera = obj.default_camera;
            end
            
            field_id = 0; % default field is zero, which means observing on non-defined field
            
            filter_name = 'unknown';
            if ~isempty(obj.pars) && ~isempty(obj.pars.filter_name)
                filter_name = obj.pars.filter_name;
            end
            
            if cs(obj.file_type, 'hdf5', 'h5')
                if obj.use_deflate
                    ext = 'h5z';
                else
                    ext = 'h5';
                end
            elseif cs(obj.file_type, 'fits')
                ext = 'fit';
            elseif cs(obj.file_type, 'mat')
                ext = 'mat';
            end
            
            % EXAMPLE FILENAME: WFAST_ZYLA_20190410-155223-035_clear_0_Raw.h5
            val = sprintf('%s_%s_%s_%s_%d_%s.%s', project, camera, time_str, filter_name, field_id, obj.product_type, ext);
            
        end
        
        function text_file_name = getReadmeFilename(obj, prefix) % generate filename for readme file appended to beginning and end of each folder
  
            if nargin<2 || isempty(prefix)
                prefix = 'A';
            end
            
            text_file_name = fullfile(obj.directory, [prefix '_README.txt']);
            
            for ii = 1:1000
                if exist(text_file_name, 'file')
                    text_file_name = fullfile(obj.directory, [prefix '_README' num2str(ii) '.txt']);
                else 
                    break;
                end
            end
            
        end
        
        function val = getMeanSaveTime(obj) % how much time is spent writing each file to disk, on average (in seconds)
            
            if obj.save_time>0
                val = obj.save_time./obj.num_saved_batches;
            else
                val = 0;
            end
            
        end
                
        function val = get.prev_name(obj) % shortcut to get the full filename of previous buf. Use this to inspect what is saved. 
            
            val = obj.prev_buf.latest_filename;
            
        end
        
    end
    
    methods % setters
        
        function set.this_buf(obj, val) % protect the current buffer in case it is reading/writing/recording
            
            if ~obj.checkInputability
                return;
            end
            
            obj.buf(obj.mod(obj.index)) = val;
            
        end
        
        function set.index_rec(obj, val) % protects index_rec from being changed by matlab (that could replace the vector and ruin everything)
            
            util.vec.mex_change(obj.index_rec_vec, 1, val-1);
            
        end
        
        function set.base_dir(obj, val) % setting this property just sets the override (hidden) property 
            
            obj.base_dir_override = val;
            
        end
        
        function set.date_dir(obj, val) % setting this property just sets the override (hidden) property 
            
            obj.date_dir_override = val;
            
        end
        
        function set.target_dir(obj, val) % setting this property just sets the override (hidden) property 
            
            obj.target_dir_override = val;
            
        end
        
        function set.directory(obj, val) % setting this property just sets the override (hidden) property 
            
            obj.directory_override = val;
            
        end
        
        function set.filename(obj, val) % setting this property just sets the override (hidden) property 
                        
            obj.filename_override = val;
            
        end
        
        function vec2times(obj, buf) % convert internal posix-time vector given by camera into date-strings

            if nargin<2 || isempty(buf)
                buf = obj.this_buf;
            end
            
            if ~isfield(buf, 't_vec')
                return;
            end
            
            if ~isempty(buf.t_vec) && buf.t_vec(1)>0
                obj.buf(buf.buf_number).t_start = util.text.time2str(datetime(buf.t_vec(1), 'ConvertFrom', 'posixtime', 'TimeZone', 'utc'));
            end
            
            if ~isempty(buf.t_vec) && buf.t_vec(2)>0
                obj.buf(buf.buf_number).t_end = util.text.time2str(datetime(buf.t_vec(2), 'ConvertFrom', 'posixtime', 'TimeZone', 'utc'));
            end
            
            if ~isempty(buf.t_vec) && buf.t_vec(3)>0
                obj.buf(buf.buf_number).t_end_stamp = buf.t_vec(3);
            end
            
        end
        
    end 
    
    methods % internal functionality
        
        function addBuffer(obj) % this is where we define the content of individual buffers (based on file.AstroData)
            
            idx = length(obj.buf)+1;

            list = properties(file.AstroData);
            temp_buf = struct;
            for ii = 1:length(list)
                temp_buf.(list{ii}) = [];
            end
            
            temp_buf.latest_filename = '';    
            temp_buf.mex_flag_record = [0 0];
            temp_buf.mex_flag_write = [0 0];
            temp_buf.mex_flag_read = [0 0];
            temp_buf.buf_number = idx;
            temp_buf.t_vec = [];
            
            if idx==1
                obj.buf = temp_buf;
            else
                obj.buf(idx) = temp_buf;
            end
            
        end
        
        function nextBuffer(obj, serial) % cycle to the next buffer (only call this when you are done processing the data!) 
            
            if obj.debug_bit>2
                disp(['moving on to buffer ' num2str(obj.next_buf.buf_number)]);
            end
            
%             obj.clear;
            util.vec.mex_change(obj.this_buf.mex_flag_read, 1, 0); % unlock the current buffer
            
            obj.index = obj.index + 1;
            if obj.index > length(obj.buf)
                obj.index = 1;
            end
            
            if nargin>1 && ~isempty(serial) % make sure the next buffer continues with the same serial number as the last (or override it)
                obj.this_buf.serial = serial;
            end
            
        end
        
        function loadDataFromBuffer(obj, idx) % copies pointers from the buffer struct into the object (remember copy on write)
            
            if nargin<2 || isempty(idx)
                idx = obj.index;
            end
            
%             util.vec.mex_change(obj.buf(idx).mex_flag_read, 1, 1); % lock the buffer for reading/processing?
            
            % load the pointers to the data stored in the buffers
            list = properties(file.AstroData);
            for ii = 1:length(list)
                obj.(list{ii}) = obj.buf(idx).(list{ii});
            end
            
        end
        
        function waitForRecording(obj, buf, timeout) % check if the current buffer has finished recording (or never started recording)
            
            if nargin<2 || isempty(buf)
                buf = obj.this_buf;
            end
            
            if nargin<3 || isempty(timeout)
                timeout = max([10, obj.pars.expT.*5*size(buf.images,3)]); % seconds
            end
            
            if obj.debug_bit>1, fprintf('waitForRecording. record_flag= %d %d\n', buf.mex_flag_record(1), buf.mex_flag_record(2)); end
            
            res = 0.01; % time resolution

            for ii = 1:timeout/res
                
                if ~obj.isRecording(buf)
                    return;
                end
                
                util.vec.mex_change(buf.mex_flag_record, 2); % how many seconds (total) have we waited for this flag...
                
                if ii==1
                    
                    if obj.debug_bit>2
                        disp(['BufferWheel: waiting for buffer ' num2str(buf.buf_number) ' to finish recording... record flag: ' util.text.print_vec(buf.mex_flag_record)]);                    
                    end
                    
                end
                
                if  ii>2 && ~isempty(obj.camera_mex_flag) && obj.camera_mex_flag(1)==0
                    timeout = ii*res;
                    if obj.debug_bit, disp('Camera stopped!, skipping wait time...'); end
                    break;
                elseif ~isempty(obj.camera_mex_flag) && obj.camera_mex_flag(2)~=0
                    timeout = ii*res;
                    if obj.debug_bit, disp(['Camera error ' num2str(obj.camera_mex_flag(2)) ', skipping wait time...']); end
                    break;
                end
                
                pause(res);
                
            end
            
            error('timeout (%4.2f sec) when waiting for buffer %d to clear from recording... ', timeout, buf.buf_number);
                        
        end
        
        function waitForReading(obj, buf, timeout) % check if the current buffer has finished reading (or never started reading)
            
            if nargin<2 || isempty(buf)
                buf = obj.this_buf;
            end
            
            if nargin<3 || isempty(timeout)
                timeout = 10; % seconds
            end
            
            res = 0.001; % time resolution
            
            for ii = 1:timeout/res
                
                if ~obj.isReading(buf)
                    return;
                end
                
                util.vec.mex_change(buf.mex_flag_read, 2); % how many seconds (total) have we waited for this flag...

                if ii==1
                    if obj.debug_bit>2
                        disp(['BufferWheel: waiting for buffer ' num2str(buf.buf_number) ' to finish reading... read flag: ' util.text.print_vec(buf.mex_flag_read)]);                    
                    end
                end
                
                pause(res);
                
            end
            
            error('timeout (%4.2f sec) when waiting for buffer %d to clear from reading... ', timeout, buf.buf_number);
                        
        end
        
        function waitForWriting(obj, buf, timeout) % check if the current buffer has finished writing (or never started writing)
            
            if nargin<2 || isempty(buf)
                buf = obj.this_buf;
            end
            
            if nargin<3 || isempty(timeout)
                timeout = 10; % seconds
            end
            
            res = 0.001; % time resolution
            
            for ii = 1:timeout/res
                
                if ~obj.isWriting(buf)
                    return;
                end
                
                util.vec.mex_change(buf.mex_flag_write, 2); % how many seconds (total) have we waited for this flag...
                
                if ii==1
                    
                    if obj.debug_bit
                        disp(['BufferWheel: waiting for buffer ' num2str(buf.buf_number) ' to finish writing... write flag: ' util.text.print_vec(buf.mex_flag_write)]);                    
                    end
                end
                
                pause(res);
                
            end
            
            error('timeout (%4.2f sec) when waiting for buffer %d to clear from writing... ', timeout, buf.buf_number);
                        
        end
        
        function unlockMexFlag(obj) % use this after errors that leave some buffers locked (careful! leave enough time to clear recording/writing to disk)
            
            for ii = 1:length(obj.buf)
                
                for jj = 1:2
                    
                    obj.buf(ii).mex_flag_record(jj) = 0;
                    obj.buf(ii).mex_flag_read(jj) = 0;
                    obj.buf(ii).mex_flag_write(jj) = 0;
                                        
                end
                
            end
            
        end
        
    end
    
    methods % actions!
        
        function input(obj, varargin) % give data to the buffer (use "save" to dump data to disk)
            
            if isempty(varargin)
                return; % do not call clear (or anything) if you didn't get any data... 
            elseif isa(varargin{1}, 'file.AstroData') % just grab the data from the object given (including derived classes e.g., img.Acquisition)
                input = varargin{1};
            else
                input = util.text.InputVars;
                input.setupDataInput; % default used for scanning all sort of inputs
                input.scan_vars(varargin{:});
            end
            
            obj.clear;
            
            list = properties(file.AstroData);
            for ii = 1:length(list) % copy the data from the input object into the buffer... 
                
                name = list{ii};
                
                if isfield(obj.this_buf, name)
                    obj.buf(obj.mod(obj.index)).(name) = input.(name);
                else
                    warning(['Cannot find field "' name '" in buffer struct']);
                end
                
            end
            
            obj.loadDataFromBuffer; % make sure the data is reflected in the object (not only the buffer)
            
            obj.gui.update;
            
        end
        
        function save(obj, buf) % dump existing data to disk
        
            import util.text.cs;
            
            if nargin<2 || isempty(buf)
                buf = obj.this_buf;
            end
            
            obj.gui.update;
            
            obj.waitForRecording(buf); % make sure the buffer is not being recorded into (by camera)
            obj.waitForWriting(buf); % make sure this buffer isn't already writing to disk...
            
            if isempty(buf.images)  && isempty(buf.cutouts) && isempty(buf.stack) && ...
                    isempty(buf.psfs) && isempty(buf.fluxes) % no images of any type, PSFs or fluxes are available. 
                return;
            end
            
            util.vec.mex_change(buf.mex_flag_write, 1, 1); % started writing...
            
            tStart = tic;
            
            % make dir if it doesn't exist yet...
            this_dir = obj.directory;            
            if ~isempty(this_dir) && ~exist(this_dir, 'dir') % create directory if doesn't exist
                if obj.debug_bit, disp(['creating directory ' this_dir]); end                
                mkdir(this_dir);
            end
            
            this_filename = fullfile(this_dir, obj.filename);
            if exist(this_filename, 'file') % what to do with pre-existing file?
                if obj.use_overwrite
                    disp(['file ' this_filename ' already exists. rewriting it now!']);
                    delete(this_filename);
                else
                    error(['file ' this_filename ' already exists, cannot save again (rewrite disabled!)']);
                end
            end
            
            obj.vec2times(buf); % transform camera output posix time vectors into formatted time strings for saving
            
            if ~isempty(obj.pars) % update the pars object when batch started/ended
                obj.pars.t_start = buf.t_start;
                obj.pars.t_end = buf.t_end;
            end
            
            obj.buf(buf.buf_number).latest_filename = this_filename; % keep track of the latest filename writen to disk
            
            if obj.debug_bit
                fprintf('Buffer %d saving file % 4d  : %s (async: %d defalte: %d)\n', buf.buf_number, obj.serial, this_filename, obj.use_async, obj.use_deflate);
            end
            
            if obj.use_mex
                
                if obj.use_write_pars && ~isempty(obj.pars) % need to replace this with util.oop.save to structure and then get mexWrite to 
                    obj.pars_struct_cell = util.oop.save(obj.pars, 'w', 'name', 'pars', 'format', 'struct');
                else
                    obj.pars_struct_cell = {};
                end
                
%                 file.mex.write(this_filename, buf.mex_flag_write, 'images', buf.images, ...
%                     'cutouts', buf.cutouts,  'positions', buf.positions,...
%                     'stack', buf.stack, 'num_sum', buf.num_sum,...
%                     'timestamps', buf.timestamps, 't_end_stamp', buf.t_end_stamp, 't_end', buf.t_end, 't_start', buf.t_start,...
%                     'psfs', buf.psfs, 'sampling_psf', buf.sampling_psf, 'fluxes', buf.fluxes, 'parameters', obj.pars_struct_cell,...
%                     'chunk', obj.chunk, 'deflate', obj.use_deflate, 'async_write', obj.use_async, 'debug_bit', obj.debug_bit);
                
                % right now mex write is only for HDF5 files    
                if obj.use_save_raw_images
                    file.mex.write(this_filename, buf.mex_flag_write, buf, 'parameters', obj.pars_struct_cell,...
                    'chunk', obj.chunk, 'deflate', obj.use_deflate, 'async_write', obj.use_async, 'debug_bit', obj.debug_bit);
                else
                    file.mex.write(this_filename, buf.mex_flag_write, buf, 'parameters', obj.pars_struct_cell,...
                    'chunk', obj.chunk, 'deflate', obj.use_deflate, 'async_write', obj.use_async, 'debug_bit', obj.debug_bit,...
                    'images', []); % last line overwrites the "images" in the buf and sets it to empty... 
                end
                
            else % use the old, non-mex writing method

                if isempty(obj.file_type) || cs(obj.file_type, {'hdf5', 'h5'})
                    obj.saveHDF5(this_filename);
                elseif cs(obj.file_type, 'fits')
                    obj.saveFits(this_filename);
                elseif cs(obj.file_type, 'mat')
                    obj.saveMatFile(this_filename);
                else
                    error('unknown file type: %s', obj.file_type);
                end
                
                util.vec.mex_change(buf.mex_flag_write, 1, 0); % done writing
            
            end
            
            obj.serial = obj.serial + 1; % keep track of the number of files saved (for all runs?)
            
            % track the average save time for this run
            this_save_time = toc(tStart);
            obj.save_time = obj.save_time + this_save_time;
            obj.num_saved_batches = obj.num_saved_batches + 1;
            
            if obj.debug_bit>1, display(['write time= ' num2str(this_save_time)]); end
            
            obj.gui.update;
            
        end
        
        function saveHDF5(obj, filename)
            
            if nargin<2 || isempty(filename)
                error('cannot run saveHDF5 without a filename!');
            end
            
            if obj.debug_bit>2, disp(['this is saveHDF5. deflate: ' num2str(obj.deflate) ' class(images)= ' class(obj.images)]); end
            
            if ~isempty(obj.images)
                
                chunk = [obj.chunk obj.chunk 1 1];
                chunk = chunk(1:ndims(obj.images));
                chunk = min(chunk, size(obj.images));
                
                if obj.use_deflate>0
                    h5create(filename, '/images', size(obj.images), 'DataType', class(obj.images), 'ChunkSize', chunk, 'Deflate', obj.use_deflate);
                else
                    h5create(filename, '/images', size(obj.images), 'DataType', class(obj.images));
                end
                
                h5write(filename, '/images', obj.images);
                
            end
            
            if ~isempty(obj.cutouts)
                
                chunk = [obj.chunk obj.chunk 1 1];
                chunk = chunk(1:ndims(obj.cutouts));
                chunk = min(chunk, size(obj.cutouts));
                
                if obj.deflate>0
                    h5create(filename, '/cutouts', size(obj.cutouts), 'DataType', class(obj.cutouts), 'ChunkSize', chunk, 'Deflate', obj.deflate);
                else
                    h5create(filename, '/cutouts', size(obj.cutouts), 'DataType', class(obj.cutouts));
                end
                
                h5write(filename, '/cutouts', obj.cutouts);
                
            end
            
            if ~isempty(obj.positions)
                h5create(filename, '/positions', size(obj.positions), 'DataType', 'double');
                h5write(filename, '/positions', obj.positions);
            end
            
            if ~isempty(obj.cutouts_bg)
                
                chunk = [obj.chunk obj.chunk 1 1];
                chunk = chunk(1:ndims(obj.cutouts_bg));
                chunk = min(chunk, size(obj.cutouts_bg));
                
                if obj.deflate>0
                    h5create(filename, '/cutouts_bg', size(obj.cutouts_bg), 'DataType', class(obj.cutouts_bg), 'ChunkSize', chunk, 'Deflate', obj.deflate);
                else
                    h5create(filename, '/cutouts_bg', size(obj.cutouts_bg), 'DataType', class(obj.cutouts_bg));
                end
                
                h5write(filename, '/cutouts_bg', obj.cutouts_bg);
                
            end
            
            if ~isempty(obj.positions_bg)
                h5create(filename, '/positions_bg', size(obj.positions_bg), 'DataType', 'double');
                h5write(filename, '/positions_bg', obj.positions_bg);
            end
            
            if ~isempty(obj.stack)
                
                chunk = [obj.chunk obj.chunk 1 1];
                chunk = chunk(1:ndims(obj.stack));
                chunk = min(chunk, size(obj.stack));
                
                if obj.deflate>0
                    h5create(filename, '/stack', size(obj.stack), 'DataType', class(obj.stack), 'ChunkSize', chunk, 'Deflate', obj.deflate);
                else
                    h5create(filename, '/stack', size(obj.stack), 'DataType', class(obj.stack));
                end
                
                h5write(filename, '/stack', obj.stack);
                
                h5writeatt(filename, '/stack', 'num_sum', obj.num_sum);
                
            end
            
            if ~isempty(obj.timestamps)
                h5create(filename, '/timestamps', size(obj.timestamps));
                h5write(filename, '/timestamps', obj.timestamps);
                
                if ~isempty(obj.t_start)
                    h5writeatt(filename, '/timestamps', 't_start', obj.t_start);
                end
                
                if ~isempty(obj.t_end_stamp)
                    h5writeatt(filename, '/timestamps', 't_end_stamp', obj.t_end_stamp);
                end
                
                if ~isempty(obj.t_end)
                    h5writeatt(filename, '/timestamps', 't_end', obj.t_end);
                end
                
            end
            
            if ~isempty(obj.psfs)
                
                chunk = [obj.chunk obj.chunk 1 1];
                chunk = chunk(1:ndims(obj.psfs));
                chunk = min(chunk, size(obj.psfs));
                
                if deflate>0
                    h5create(filename, '/psfs', size(obj.psfs), 'ChunkSize', chunk, 'Deflate', deflate);
                else
                    h5create(filename, '/psfs', size(obj.psfs));
                end
                
                h5write(filename, '/psfs', obj.psfs);
                
                h5writeatt(filename, '/psfs', 'psf_sampling', obj.psf_sampling);
                
            end
            
            if ~isempty(obj.fluxes)
                h5create(filename, '/fluxes', size(obj.fluxes,2))
                h5write(filename, '/fluxes', obj.fluxes);
            end
                        
            if obj.use_write_pars && ~isempty(obj.pars)
                util.oop.save(obj.pars, filename, 'name', 'pars', 'append', 1);
            end
            
        end
        
        function saveFits(obj, filename) % need to finish this!!!
            
            if nargin<2 || isempty(filename)
                error('cannot run saveFits without a filename!');    
            end
            
            file_ptr = matlab.io.fits.createFile(filename);
            file_cleanup = onCleanup(@()matlab.io.fits.closeFile(file_ptr));
            
            if ~isempty(obj.images)

                chunk = [obj.chunk obj.chunk 1 1];
                chunk = chunk(1:ndims(obj.images));
                chunk = min(chunk, size(obj.images));

                S = size(obj.images);

                if obj.use_deflate>0
                    matlab.io.fits.setCompressionType(file_ptr,'GZIP2');
                    matlab.io.fits.setTileDim(file_ptr, chunk);
                end

                if isa(obj.images, 'uint16')

                    matlab.io.fits.createImg(file_ptr, 'int16', S);
                    matlab.io.fits.writeImg(file_ptr, int16(double(obj.images)-2^15));

                    matlab.io.fits.writeKey(file_ptr, 'BSCALE', 1);
                    matlab.io.fits.writeKey(file_ptr, 'BZERO', 2^15);

                elseif isa(obj.images, 'double')

                    matlab.io.fits.createImg(file_ptr, 'double', S);
                    matlab.io.fits.writeImg(file_ptr, obj.images);

                end

            end

            if ~isempty(obj.positions)

            end

            if ~isempty(obj.timestamps)

            end

            if ~isempty(obj.psfs)

                if chunk(1)>size(obj.psfs,1)
                    chunk(1) = size(obj.psfs,1);
                end

                if chunk(2)>size(obj.psfs,2)
                    chunk(2) = size(obj.psfs,2);
                end

                chunk = chunk(1:ndims(obj.psfs));

                if obj.use_deflate>0

                else

                end

                % can we write these attributes into the FITS header?
%                     if isempty(psf_sampling)
%                         psf_sampling = 2;
%                     end

%                     h5writeatt(filename, '/psfs', 'psf_sampling', psf_sampling);

            end

            if ~isempty(obj.fluxes)

            end

            if obj.use_write_pars && ~isempty(obj.pars)
                obj.pars.writeFitsHeader(file_ptr);
            end

        end
        
        function saveMatFile(obj, filename)
                       
            if nargin<2 || isempty(filename)
                error('cannot run saveMatFile without a filename!');    
            end
            
            data = file.AstroData(obj); % cast this object into an AstroData to be saved
            
            util.oop.save(data, filename, 'format', 'mat');
            
        end
        
    end
    
    methods % gui stuff
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = file.gui.BufferGUI(obj);                
            end
            
            obj.gui.make;
            
        end
        
    end
      
    methods (Static=true) % static save (old, non-mex methods) to be depricated!!!
        
        function ok = saveHDF5Static(filename, deflate, chunk, images, positions, timestamps, t_end_stamp, t_end, t_start, psfs, psf_sampling, fluxes, pars, use_write_pars, debug_bit)
            
            ok = 0;
            
            if nargin<1
                error('cannot run saveHDF5Static without a filename!');
            end
            
            if nargin<2 || isempty(deflate), deflate = 0; end
            if nargin<3 || isempty(chunk), chunk = 64; end
            if nargin<4, images = []; end
            if nargin<5, positions = []; end
            if nargin<6, timestamps = []; end
            if nargin<7, t_end_stamp = []; end
            if nargin<8, t_end = []; end
            if nargin<9, t_start = []; end
            if nargin<10, psfs = []; end
            if nargin<11, psf_sampling = []; end
            if nargin<12, fluxes = []; end
            if nargin<13, pars = []; end
            if nargin<14 || isempty(use_write_pars), use_write_pars = 0; end            
            if nargin<15 || isempty(debug_bit), debug_bit = 0; end
            
            chunk = [chunk chunk 1 1];
            
            if debug_bit>2, disp(['this is saveHDF5Static. deflate: ' num2str(deflate) ' class(images)= ' class(images)]); end
                                   
            if ~isempty(images)
                
                chunk = chunk(1:ndims(images));
                chunk = min(chunk, size(images));
                                
                if deflate>0
                    h5create(filename, '/images', size(images), 'DataType', class(images), 'ChunkSize', chunk, 'Deflate', deflate);
                else
                    h5create(filename, '/images', size(images), 'DataType', class(images));
                end
                
                h5write(filename, '/images', images);
                
            end
               
            if ~isempty(positions)
                h5create(filename, '/positions', size(positions), 'DataType', 'double');
                h5write(filename, '/positions', positions);                
            end
            
            if ~isempty(timestamps)                
                h5create(filename, '/time', size(timestamps));
                h5write(filename, '/time', timestamps);
                
                if ~isempty(t_end_stamp)
                    h5writeatt(filename, '/time', 't_end_stamp', t_end_stamp);
                end
                
                if ~isempty(t_end)
                    h5writeatt(filename, '/time', 'file_write_timestr', t_end);
                end
                
                if ~isempty(t_start)
                    h5writeatt(filename, '/time', 'file_start_timestr', t_start);
                end
                
            end
            
            
            if ~isempty(psfs)
                
                if chunk(1)>size(psfs,1)
                    chunk(1) = size(psfs,1);
                end
                
                if chunk(2)>size(psfs,2)
                    chunk(2) = size(psfs,2);
                end
                
                chunk = chunk(1:ndims(psfs));
                
                if deflate>0
                    h5create(filename, '/psfs', size(psfs), 'ChunkSize', chunk, 'Deflate', deflate);
                else
                    h5create(filename, '/psfs', size(psfs));
                end
                
                h5write(filename, '/psfs', psfs);
                
                if ~isempty(psf_sampling)
                    h5writeatt(filename, '/psfs', 'psf_sampling', psf_sampling);
                end
                                
            end
                        
            if ~isempty(fluxes)
                h5create(filename, '/fluxes', size(fluxes,2))
                h5write(filename, '/fluxes', fluxes);
            end
                        
            if use_write_pars && ~isempty(pars)
                util.oop.save(pars, filename, 'append', 1);
            end
                        
            ok = 1;
            
        end
        
        function ok = saveFitsStatic(filename, deflate, chunk, images, positions, timestamps, t_end_stamp, t_end, t_start, psfs, psf_sampling, fluxes, pars, use_write_pars, debug_bit)
                            
            if nargin<1
                error(['cannot run saveFitsStatic without a filename!']);    
            end
            
            if nargin<2 || isempty(deflate), deflate = 0; end
            if nargin<3 || isempty(chunk), chunk = 64; end
            if nargin<4, images = []; end
            if nargin<5, positions = []; end
            if nargin<6, timestamps = []; end            
            if nargin<7, t_end_stamp = []; end
            if nargin<8, t_end = []; end
            if nargin<9, t_start = []; end
            if nargin<10, psfs = []; end
            if nargin<11, psf_sampling = []; end
            if nargin<12, fluxes = []; end
            if nargin<13, pars = []; end
            if nargin<14, use_write_pars = 0; end            
            if nargin<15 || isempty(debug_bit), debug_bit = 0; end        
            
            if isempty(images), return; end
            
            if debug_bit>2, disp(['this is saveFitsStatic. deflate: ' num2str(deflate) ' class(images)= ' class(images)]); end
            
            file_ptr = matlab.io.fits.createFile(filename);
                     
            chunk = [chunk chunk 1 1];
            
            try
                                   
                if ~isempty(images)
                    
                    S = size(images);
                    
                    if length(chunk)<length(S)
                        chunk = [chunk ones(1, length(S)-length(chunk))];
                    else
                        chunk = chunk(1:length(S));
                    end
                    
                    if deflate>0
                        matlab.io.fits.setCompressionType(file_ptr,'GZIP2');
                        matlab.io.fits.setTileDim(file_ptr, chunk);
                    end
                    
                    if isa(images, 'uint16')
                    
                        matlab.io.fits.createImg(file_ptr, 'int16', S);
                        matlab.io.fits.writeImg(file_ptr, int16(double(images)-2^15));

                        matlab.io.fits.writeKey(file_ptr, 'BSCALE', 1);
                        matlab.io.fits.writeKey(file_ptr, 'BZERO', 2^15);

                    elseif isa(images, 'double')
                        
                        matlab.io.fits.createImg(file_ptr, 'double', S);
                        matlab.io.fits.writeImg(file_ptr, images);

                    end
                        
                end

                if ~isempty(positions)
                    
                end

                if ~isempty(timestamps)
                    
                end

                if ~isempty(psfs)

                    if chunk(1)>size(psfs,1)
                        chunk(1) = size(psfs,1);
                    end

                    if chunk(2)>size(psfs,2)
                        chunk(2) = size(psfs,2);
                    end

                    chunk = chunk(1:ndims(psfs));

                    if deflate>0
                        
                    else
                        
                    end

                    % can we write these attributes into the FITS header?
%                     if isempty(psf_sampling)
%                         psf_sampling = 2;
%                     end

%                     h5writeatt(filename, '/psfs', 'psf_sampling', psf_sampling);

                end

                if ~isempty(fluxes)
                    
                end

                if use_write_pars && ~isempty(pars)
                    pars.writeFitsHeader(file_ptr);
                end

                if ~isempty(backup_filename)
                    if debug_bit>1, disp(['backing up to file: ' backup_filename]); end
                    if strcmp(filename, backup_filename)
                        if debug_bit, disp('filename and backup_filename are the same!'); end
                    else
                        copyfile(filename, backup_filename);
                        if debug_bit>1, disp('backup finished!'); end
                    end
                end

                ok = 1;

            catch ME
                matlab.io.fits.closeFile(file_ptr);    
                rethrow(ME);
            end
            
            matlab.io.fits.closeFile(file_ptr);
            
        end
        
        function ok = saveMatFileStatic(filename, images, positions, timestamps, t_end_stamp, t_end, t_start, psfs, psf_sampling, fluxes, pars, use_save_pars, debug_bit)
            
            ok = 0;
            
            if nargin<2, images = []; end
            if nargin<3, positions = []; end            
            if nargin<4, timestamps = []; end
            if nargin<5, t_end_stamp = []; end
            if nargin<6, t_end = []; end
            if nargin<7, t_start = []; end            
            if nargin<8, psfs = []; end
            if nargin<9, psf_sampling = []; end
            if nargin<10, fluxes = []; end
            if nargin<11, pars = []; end
            if nargin<12, use_save_pars = ''; end
            if nargin<13 || isempty(debuf_bit), debug_bit = 0; end
            
            name = util.text.extension(filename, '.mat');
            
            if use_save_pars
                save(name, 'images', 'positions', 'timestamps', 'psfs', 'psf_sampling', 'fluxes', 'pars');
            else
                save(name, 'images', 'positions', 'timestamps', 'psfs', 'psf_sampling', 'fluxes');
            end
            
            ok = 1;
            
        end
                
    end
    
end