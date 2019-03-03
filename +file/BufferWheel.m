classdef BufferWheel < file.AstroData
% Usage: obj = BufferWheel(num_buffers=5, pars=[])
% Holds several (e.g., 5) Buffer structs and rotates them as they get filled
% and emptied to file. 
%
% Uses index (and index_rec) to keep track of which buffer to use when
% recording, reading and writing to file. Each buffer struct can have all 
% the inputs and outputs (e.g. images, lightcurves, psfs, etc.). 
%
% Model 1: input data from camera/sim/pipeline using input method, then
% immediately use save method which also brings up new buffer. If you try
% to input to a buffer that is still saving, you get an error. 
% (this is used for all but FullRaw data products).
% 
% Model 2: Give the camera all the buffers and have it track which one it
% is writing to using index_rec. If buffer is still saving or in analysis, 
% delay camera. When camera finishes, it gets the next buffer in line. 
% Run analysis on any buffers released from camera (if analysis slows down,
% you may need to catch up or skip buffers). Keep track of which buffer is
% under analysis at any time using index. Use "save" to save one or more
% buffers (defined by num2save) when e.g. there is a trigger in the buffer
% you are now analyzing. 
%
% OBJECTS:
%   -buf: vector of Buffer structs. 
%   -pars: link back to parent object.
%   -gui: Optional GUI for all contained Buffers. 
%   -this_buf, next_buf, prev_buf: handles into "buf", based on index. 
%
% SWITCHES & CONTROLS:
%   -product_type: string telling which kind of output (raw? cutouts?
%    lightcurves? dark/flat? etc.). Will print this in the filename. 
%   -base_dir: set to [] to use getenv('DATA') as the base dir. 
%   -date_dir: set to [] to use this night's date dir. 
%   -target_dir: set to [] to use "pars.target_name" as target_dir. 
%     *** full-path will be set to base_dir/date_dir/target_dir. ***
%   -dir_extension: overrides the automatic extension to the dir name. 
%   -use_dir_type: automatically adds "product_type" to the dir name. 
%   -use_overwrite: silently delete existing files (default is true).
%   -use_write_pars: try to save the "pars" object into file.
%   -use_deflate: if non-zero, use deflation (can be anything from 1 to 10)
%   -use_async: write in separate thread/job. 
%   -use_mex: use the mex function to write to file (HDF5 only, for now). 
%   -use_disable_mex_async_deflate: when this is enabled, doesn't let you
%    add deflate when mex and async are used (avoids crashes). 
%   -chunk: used for deflation. Default is [64 64].

    properties(Transient=true)
        
        gui@file.gui.BufferGUI;
        
    end

    properties % objects
        
        pars@head.Parameters;
        buf@struct;
        pars_struct_cell = {};
        
    end
    
    properties % switches
        
%         dark_applied;
%         flat_applied;
%         is_sim;
%         is_dark;
%         is_flat;
        
        index = 1;
        product_type = 'RawFull'; % can be RawFull, RawCut, dark, flat, CalCut, CalSum, LCs, PSFs, etc...         
        dir_extension = ''; % overrides any automatic directory extensions... 
        use_dir_types = 0; % automatically add a dir_extension equal to "type"
        
        % write parameters
        use_overwrite = 1;
        use_write_pars = 1;        
        use_deflate = 0; % also can specify how much deflation you want (1 to 10)
        use_async = 0;
        use_mex = 1;
        use_datestrings = 1;
        chunk = 64;
        
        file_type = 'hdf5'; % support is added for "fits", "mat", and "tiff"
                
        serial = 1;
        
        save_time = 0;
        num_saved_batches = 0;
        
        timeout = 5; % defualt timeout, in seconds
        
        debug_bit = 1;
                
    end
    
    properties(Dependent=true)
        
        index_rec;
            
        this_buf;
        prev_buf;
        next_buf;
               
        im_size;
        
        base_dir;
        date_dir;
        target_dir;
        directory;
        filename;
        prev_name;
                
    end
    
    properties(Hidden=true)
        
        index_rec_vec = [0 0];
        
        % naming convention
        base_dir_override;
        date_dir_override;
        target_dir_override;
        directory_override;
        filename_override;
        dark_dir_name = 'dark'; % what to call the directory for darks
        flat_dir_name = 'flat'; % what to call the directory for flats
        
        default_product_type;
        default_dir_extension;
        default_use_dir_types;
        default_project = 'Kraar';
        default_camera = 'Zyla';
        
        default_use_overwrite;
        default_use_write_pars;
        default_use_deflate;
        default_use_asnc;
        default_use_mex;
        default_chunk;
        default_file_type;
        
        use_disable_mex_async_deflate = 1; % if you try to defalte using mex file and async write, it will silently cancel deflation... if you try to do this Matlab will crash...
        
%         pars_struct; % we must keep a copy of this object so it doesn't get deleted until it is written by mex file...
        
        camera_mex_flag;
        
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = BufferWheel(num_buffers, pars)
            
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
                pars = head.Parameters;
            end
            
            for ii = 1:num_buffers
                obj.addBuffer;
            end
            
            obj.pars = pars;
            
            util.oop.save_defaults(obj); % we don't have any defaults here
            
        end
        
    end
    
    methods % resetters
        
        function initialize(obj)
                      
            util.oop.load_defaults(obj);
            obj.reset;
                        
        end
        
        function reset(obj)
                        
            for ii = 1:length(obj.buf)
                
                obj.buf(ii).images = [];
                obj.buf(ii).images_raw = [];
                obj.buf(ii).images_cal = [];
                obj.buf(ii).cutouts_raw = [];
                obj.buf(ii).cutouts_cal = [];
                obj.buf(ii).full_sum = [];
                obj.buf(ii).positions = [];
                obj.buf(ii).num_sum = [];
                
                obj.buf(ii).timestamps = [];
                obj.buf(ii).t_start = [];
                obj.buf(ii).t_end = [];
                obj.buf(ii).t_vec = [];
                
                obj.buf(ii).psfs = [];
                obj.buf(ii).psf_sampling = [];
                obj.buf(ii).latest_filename = [];
                
                util.vec.mex_change(obj.buf(ii).mex_flag_record, 2, 0);
                util.vec.mex_change(obj.buf(ii).mex_flag_read, 2, 0);
                util.vec.mex_change(obj.buf(ii).mex_flag_write, 2, 0);
                
            end
            
            obj.index = 1;
            obj.index_rec = 1;
            
            obj.unlockMexFlag;
            
            obj.pars_struct_cell = {};
            
        end
                
        function clear(obj)
        
            clear@file.AstroData(obj);
            
            for ii = 1:length(obj.buf)
                obj.buf(ii).t_start = [];
                obj.buf(ii).t_end = [];
                obj.buf(ii).t_end_stamp = [];
            end
            
        end
        
    end 
    
    methods % getters

        function val = isWriting(obj, buf)
            
            if nargin<2 || isempty(buf)
                buf = obj.this_buf;
            end
            
            flag = buf.mex_flag_write;
            val = flag(1)==1; % && flag(2)~=1;
            
        end
        
        function val = isRecording(obj, buf)
            
            if nargin<2 || isempty(buf)
                buf = obj.this_buf;
            end
            
            flag = buf.mex_flag_record;
            val = flag(1)==1; %&& flag(2)~=1;
%             val = val || ~buf.mex_flag_read(2); % add second check to see if there is any data in the buffer at all (if not, consider it still under recording...)
            
        end
        
        function val = isReading(obj, buf)
            
            if nargin<2 || isempty(buf)
                buf = obj.this_buf;
            end
            
            flag = buf.mex_flag_read;
            val = flag(1)==1; %&& flag(2)~=1;
            
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
          
        function val = get.index_rec(obj)
            
            val = obj.index_rec_vec(1) + 1;
            
        end
        
        function val = mod(obj, index)
            
            val = mod(index-1, length(obj.buf))+1;
            
        end
        
        function val = get.this_buf(obj)
            
            val = obj.buf(obj.mod(obj.index));
            
        end
        
        function val = get.next_buf(obj)
           
            val = obj.buf(obj.mod(obj.index+1));
            
        end
        
        function val = get.prev_buf(obj)
           
            val = obj.buf(obj.mod(obj.index-1));
            
        end
        
        function val = get.im_size(obj)
           
            if ~isempty(obj.this_buf.images)
                val = util.vec.imsize(obj.this_buf.images);
            elseif ~isempty(obj.this_buf.psfs)
                val = util.vec.imsize(obj.this_buf.psfs);
            else
                val = [];
            end
            
        end
        
        function val = get.use_deflate(obj)
            
            val = obj.use_deflate;
            
            if obj.use_disable_mex_async_deflate && obj.use_mex && obj.use_async
                val = 0;
            end
            
        end                
        
        function val = get.base_dir(obj)
            
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
        
        function date = current_datetime(obj)
            
            if ~isempty(obj.this_buf.t_end)
                date = util.text.str2time(obj.this_buf.t_end); % best time estimate is when the file has finished recording (matched from camera)
            else
                date = datetime('now', 'TimeZone', 'UTC'); % just use what ever time is right now...
            end
            
        end
        
        function val = get.date_dir(obj) 
                       
            if ~isempty(obj.date_dir_override)
                val = obj.date_dir_override;                
            else
                
                date = obj.current_datetime;
                
                % make sure each night is in the same folder, even after midnight.
                if date.Hour<12 % this is 14:00 local time, pretty safe to say it is after midnight of the next night...
                    date.Day = date.Day - 1;
                end
                
                val = datestr(date, 'yyyy-mm-dd/');
                
            end
            
        end
        
        function val = get.target_dir(obj)
            
            if ~isempty(obj.target_dir_override)
                val = obj.target_dir_override;                
            elseif ~isempty(obj.pars) && ~isempty(obj.pars.folder_name)
                val = obj.pars.folder_name;
            else
                val = 'run1';
            end
            
            if ~isempty(obj.dir_extension)
                val = util.text.sub_append(val, obj.dir_extension);
            elseif obj.use_dir_types && ~isempty(obj.product_type)
                val = util.text.sub_append(val, obj.product_type);
            end
            
        end
        
        function val = get.directory(obj)
            
            if ~isempty(obj.directory_override)
                val = obj.directory_override;
            else
                val = obj.makeFullpath;
            end
            
        end
        
        function val = get.filename(obj)
            
            if ~isempty(obj.filename_override)
                val = obj.filename_override;
            else
                val = obj.makeFilename;
            end
            
        end
                
        function val = makeFullpath(obj)
            
            d = obj.target_dir;
            
            if util.text.cs(obj.product_type, 'dark')
                d = obj.dark_dir_name;
            elseif util.text.cs(obj.product_type, 'flat')
                d = obj.flat_dir_name;
            end
            
            val = fullfile(obj.base_dir, obj.date_dir, d);
            
        end
        
        function val = makeFilename(obj)
           
            import util.text.cs;
            
            date = obj.current_datetime;
            
            str = datestr(date, 'yyyymmdd-HHMMSS-FFF');
            
            if isempty(obj.pars)
                project = obj.default_project;
                camera = obj.default_camera;
            else
                project = obj.pars.project;
                camera = obj.pars.cam_name;
            end
                
            ccd_id = 0;
            amp_id = 0;
            
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
            
            val = sprintf('%s_%s_%s_%d_%d_%s_%06d.%s', project, camera, str, ccd_id, amp_id, obj.product_type, obj.serial, ext);
%             val = sprintf('test_%d.%s', obj.serial, ext); % debug
            
        end
        
        function text_file_name = getReadmeFilename(obj, prefix)
  
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
        
        function val = getMeanSaveTime(obj)
           
            if obj.save_time>0
                val = obj.save_time./obj.num_saved_batches;
            else
                val = 0;
            end
            
        end
                
        function val = get.prev_name(obj)
            
            val = obj.prev_buf.latest_filename;
            
        end
        
    end
    
    methods % setters
        
        function set.this_buf(obj, val)
            
            if ~obj.checkInputability
                return;
            end
            
            obj.buf(obj.mod(obj.index)) = val;
            
        end
        
        function set.index_rec(obj, val)
            
            util.vec.mex_change(obj.index_rec_vec, 1, val-1);
            
        end
        
        function set.base_dir(obj, val)
            
            obj.base_dir_override = val;
            
        end
        
        function set.date_dir(obj, val)
            
            obj.date_dir_override = val;
            
        end
        
        function set.target_dir(obj, val)
            
            obj.target_dir_override = val;
            
        end
        
        function set.directory(obj, val)
            
            obj.directory_override = val;
            
        end
        
        function set.filename(obj, val)
                        
            obj.filename_override = val;
            
        end
        
        function vec2times(obj, buf)

            if nargin<2 || isempty(buf)
                buf = obj.this_buf;
            end
            
            if ~isfield(buf, 't_vec')
                return;
            end
            
            if ~isempty(buf.t_vec) && buf.t_vec(1)>0 && isempty(buf.t_start)
                obj.buf(buf.buf_number).t_start = util.text.time2str(datetime(buf.t_vec(1), 'ConvertFrom', 'posixtime', 'TimeZone', 'utc'));
            end
            
            if ~isempty(buf.t_vec) && buf.t_vec(2)>0 && isempty(buf.t_end)
                obj.buf(buf.buf_number).t_end = util.text.time2str(datetime(buf.t_vec(2), 'ConvertFrom', 'posixtime', 'TimeZone', 'utc'));
            end
            
            if ~isempty(buf.t_vec) && buf.t_vec(3)>0 && isempty(buf.t_end_stamp)
                obj.buf(buf.buf_number).t_end_stamp = buf.t_vec(3);
            end
            
        end
        
    end 
    
    methods % actions
        
        function addBuffer(obj) % this is where we define the content of individual buffers
            
            idx = length(obj.buf)+1;

            list = properties(file.AstroData);
            temp_buf = struct;
            for ii = 1:length(list)
                temp_buf.(list{ii}) = [];
            end
            
            temp_buf.latest_filename = '';
%             temp_buf.mex_flag_record = [0 0 0 0];
%             temp_buf.mex_flag_write = [0 0 0 0];
%             temp_buf.mex_flag_read = [0 0 0 0];            
            temp_buf.mex_flag_record = [0 0];
            temp_buf.mex_flag_write = [0 0];
            temp_buf.mex_flag_read = [0 0];
            temp_buf.buf_number = idx;
            temp_buf.t_vec = [];
            
%             temp_buf = struct('images', [], 'images_raw', [], 'images_cal', [],...
%                     'cutouts_raw', [], 'cutouts_cal', [], 'positions', [], ...
%                     'full_sum', [], 'num_sum', [],...
%                     'timestamps', [], 't_start', [], 't_end', [], 't_end_stamp', [], 't_vec', [], ...
%                     'psfs', [], 'psf_sampling', [],...
%                     'lightcurves', [], 'latest_filename', [], ...
%                     'mex_flag_record', [0 0 0 0], 'mex_flag_write', [0 0 0 0], 'mex_flag_read', [0 0 0 0],...
%                     'buf_number', idx);
            
            if idx==1
                obj.buf = temp_buf;
            else
                obj.buf(idx) = temp_buf;
            end
            
        end
        
        function nextBuffer(obj, serial) % only call this when you are done processing the data! 
            
            if obj.debug_bit>2
                disp(['moving on to buffer ' num2str(obj.next_buf.buf_number)]);
            end
            
%             obj.loadDataFromBuffer(obj.next_buf);
            obj.clear;
%             obj.markAsRead; % the data in obj is now overwritten by the data from next_buf
            util.vec.mex_change(obj.this_buf.mex_flag_read, 1, 0)            
            obj.index = obj.index + 1;
            if obj.index > length(obj.buf)
                obj.index = 1;
            end
            
            if nargin>1 && ~isempty(serial)
                obj.this_buf.serial = serial;
            end
                        
        end
        
        function loadDataFromBuffer(obj, idx) % copies pointers from the buffer into obj
            
%             disp('loading data from buffer');
            
            if nargin<2 || isempty(idx)
                idx = obj.index;
            end
            
%             util.vec.mex_change(obj.buf(idx).mex_flag_read, 1, 1);
            
            % load the pointers to the data stored in the buffers
            list = properties(file.AstroData);
            for ii = 1:length(list)
                obj.(list{ii}) = obj.buf(idx).(list{ii});
            end
            
%             util.vec.mex_change(obj.buf(idx).mex_flag_read, 1, 0);
%             util.vec.mex_change(obj.buf(idx).mex_flag_read, 3, 0);

%             fprintf('Loading data finished, idx= %d, flag(1)= %d, flag(2)= %d\n', idx, obj.buf(idx).mex_flag_read(1), obj.buf(idx).mex_flag_read(2));
            
        end
        
%         function markAsRead(obj, buf) % release the buffer to rewrite (only do this when done processing!)
%             
%             if nargin<2 || isempty(buf)
%                 buf = obj.this_buf;
%             end
%             
% %             util.vec.mex_change(obj.buf(index).mex_flag_read, 1, 0);
%             util.vec.mex_change(buf.mex_flag_read, 2, 1);
%             util.vec.mex_change(buf.mex_flag_read, 3, 0);
%             
%         end
        
        function waitForRecording(obj, buf, timeout) % check if the current buffer has finished recording (or never started recording)
            
            if nargin<2 || isempty(buf)
                buf = obj.this_buf;
            end
            
            if nargin<3 || isempty(timeout)
                timeout = max(10, obj.pars.T.*5*size(buf.images,3)); % seconds
            end
            
            if obj.debug_bit>1, fprintf('waitForRecording. record_flag= %d %d\n', buf.mex_flag_record(1), buf.mex_flag_record(2)); end
            
            res = 0.01; % time resolution
%             buf.mex_flag_read(3)
            for ii = 1:timeout/res
                
                if ~obj.isRecording(buf)
                    return;
                end
                
                util.vec.mex_change(buf.mex_flag_record, 2); % how many milliseconds (total) have we waited for this flag...
                
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
                
                util.vec.mex_change(buf.mex_flag_record, 2); % how many milliseconds (total) have we waited for this flag...

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
                
                util.vec.mex_change(buf.mex_flag_write, 2); % how many milliseconds (total) have we waited for this flag...
                
                if ii==1
                    
                    if obj.debug_bit
                        disp(['BufferWheel: waiting for buffer ' num2str(buf.buf_number) ' to finish writing... write flag: ' util.text.print_vec(buf.mex_flag_write)]);                    
                    end
                end
                
                pause(res);
                
            end
            
            error('timeout (%4.2f sec) when waiting for buffer %d to clear from writing... ', timeout, buf.buf_number);
                        
        end
        
        function unlockMexFlag(obj)
            
            for ii = 1:length(obj.buf)
                
                for jj = 1:2
                    
                    obj.buf(ii).mex_flag_record(jj) = 0;
                    obj.buf(ii).mex_flag_read(jj) = 0;
                    obj.buf(ii).mex_flag_write(jj) = 0;
                                        
                end
                
            end
            
        end
        
        function input(obj, varargin)
            
            if isempty(varargin)
                % pass
            elseif isa(varargin{1}, 'file.AstroData')
                input = varargin{1};
                list = properties(file.AstroData);
            else
                input = util.text.InputVars;
                input.setupDataInput; % default used for scanning all sort of inputs
                input.scan_vars(varargin{:});
                list = input.list_added_properties;
            end
            
            obj.clear;
            
%             obj.setHasData;
            
            for ii = 1:length(list)
                
                name = list{ii};
                
                if isfield(obj.this_buf, name)
                    obj.buf(obj.mod(obj.index)).(name) = input.(name);
                else
                    warning(['Cannot find field "' name '" in buffer struct']);
                end
                
            end
            
            obj.gui.update;
            
        end
        
%         function setHasData(obj)
%             
%             util.vec.mex_change(obj.this_buf.mex_flag_read, 3, 1); % tells the buffer there is something to read
%             
%         end
        
        function startWriting(obj, buf)
            
            if nargin<2 || isempty(buf)
                buf = obj.this_buf;
            end
            
            util.vec.mex_change(buf.mex_flag_write, 1, 1);
            
        end
        
        function save(obj, buf)
        
            if nargin<2 || isempty(buf)
                buf = obj.this_buf;
            end
            
            import util.text.cs;
            
            obj.gui.update;
            
            obj.waitForRecording(buf);
            obj.waitForWriting(buf);
            
            if isempty(buf.images) && isempty(buf.images) && isempty(buf.images_proc) && ...
                    isempty(buf.cutouts) && isempty(buf.cutouts_proc) && isempty(buf.stack) && ...
                    isempty(buf.psfs) && isempty(buf.lightcurves) % no images of any type, PSFs or lightcurves are available. 
                return;
            end
            
%             util.vec.mex_change(buf.mex_flag_write, 2, 0); % not done writing...
            util.vec.mex_change(buf.mex_flag_write, 1, 1); % started writing...
            
            tStart = tic;
            
            % make dir if it doesn't exist yet...
            this_dir = obj.directory;            
            if ~isempty(this_dir) && ~exist(this_dir, 'dir') % create directory if doesn't exist
                if obj.debug_bit, disp(['creating directory ' this_dir]); end                
                mkdir(this_dir);
            end
            
            filename = fullfile(obj.directory, obj.filename);
            if exist(filename, 'file') % what to do with pre-existing file?
                if obj.use_overwrite
                    disp(['file ' filename ' already exists. rewriting it now!']);
                    delete(filename);
                else
                    warning(['file ' filename ' already exists, cannot save again (rewrite disabled!)']);
                end
            end
            
            obj.vec2times(buf);
            
            if ~isempty(obj.pars) % update the pars object when batch started/ended
                obj.pars.t_start = buf.t_start;
                obj.pars.t_end = buf.t_end;
            end
            
            obj.buf(buf.buf_number).latest_filename = filename;
            
            if obj.debug_bit
                fprintf('Buffer %d saving file % 4d  : %s (async: %d defalte: %d)\n', buf.buf_number, obj.serial, filename, obj.use_async, obj.use_deflate);
            end
            
            if obj.use_mex
                
                if obj.use_write_pars && ~isempty(obj.pars) % need to replace this with util.oop.save to structure and then get mexWrite to 
                    obj.pars_struct_cell = util.oop.save(obj.pars, 'w', 'name', 'pars', 'format', 'struct');
                else
                    obj.pars_struct_cell = {};
                end
                
                % right now mex write is only for HDF5 files
                file.mex.write(filename, buf.mex_flag_write, 'images', buf.images, 'images_proc', buf.cutouts_proc, ...
                    'cutouts', buf.cutouts, 'cutouts_proc', buf.cutouts_proc, 'positions', buf.positions,...
                    'stack', buf.stack, 'num_sum', buf.num_sum,...
                    'timestamps', buf.timestamps, 't_end_stamp', buf.t_end_stamp, 't_end', buf.t_end, 't_start', buf.t_start,...
                    'psfs', buf.psfs, 'sampling_psf', buf.sampling_psf, 'lightcurves', buf.lightcurves, 'parameters', obj.pars_struct_cell,...
                    'chunk', obj.chunk, 'deflate', obj.use_deflate, 'async_write', obj.use_async, 'debug_bit', obj.debug_bit);
                    % should also mark the mex_flag_write as finished writing... 
                    
            else % use the old saveHDF5Static+parpool saving method

                if isempty(obj.file_type) || cs(obj.file_type, {'hdf5', 'h5'})
                    obj.saveHDF5(filename);
                elseif cs(obj.file_type, 'fits')
                    obj.saveFits(filename);
                elseif cs(obj.file_type, 'mat')
                    obj.saveMatFile(filename);
                else
                    error('unknown file type: %s', obj.file_type);
                end
                
                util.vec.mex_change(buf.mex_flag_write, 1, 0); % done writing
            
            end
            
            obj.serial = obj.serial + 1;
            
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
            
            if ~isempty(obj.images_raw)
                
                chunk = [obj.chunk obj.chunk 1 1];
                chunk = chunk(1:ndims(obj.images_raw));
                chunk = min(chunk, size(obj.images_raw));
                                
                if obj.use_deflate>0
                    h5create(filename, '/images_raw', size(obj.images_raw), 'DataType', class(obj.images_raw), 'ChunkSize', chunk, 'Deflate', obj.use_deflate);
                else
                    h5create(filename, '/images_raw', size(obj.images_raw), 'DataType', class(obj.images_raw));
                end
                
                h5write(filename, '/images_raw', obj.images_raw);
                
            end
            
            if ~isempty(obj.cutouts_raw)
                
                chunk = [obj.chunk obj.chunk 1 1];
                chunk = chunk(1:ndims(obj.cutouts_raw));
                chunk = min(chunk, size(obj.cutouts_raw));
                
                if obj.deflate>0
                    h5create(filename, '/cutouts_raw', size(obj.cutouts_raw), 'DataType', class(obj.cutouts_raw), 'ChunkSize', chunk, 'Deflate', obj.deflate);
                else
                    h5create(filename, '/cutouts_raw', size(obj.cutouts_raw), 'DataType', class(obj.cutouts_raw));
                end
                
                h5write(filename, '/cutouts_raw', obj.cutouts_raw);
                
            end
            
            if ~isempty(obj.positions)
                h5create(filename, '/positions', size(obj.positions), 'DataType', 'double');
                h5write(filename, '/positions', obj.positions);
            end
            
            if ~isempty(obj.full_sum)
                
                chunk = [obj.chunk obj.chunk 1 1];
                chunk = chunk(1:ndims(obj.full_sum));
                chunk = min(chunk, size(obj.full_sum));
                
                if obj.deflate>0
                    h5create(filename, '/full_sum', size(obj.full_sum), 'DataType', class(obj.full_sum), 'ChunkSize', chunk, 'Deflate', obj.deflate);
                else
                    h5create(filename, '/full_sum', size(obj.full_sum), 'DataType', class(obj.full_sum));
                end
                
                h5write(filename, '/full_sum', obj.full_sum);
                
                h5writeatt(filename, '/full_sum', obj.num_sum);
                
            end
            
            if ~isempty(obj.timestamps)
                h5create(filename, '/time', size(obj.timestamps));
                h5write(filename, '/time', obj.timestamps);
                
                if ~isempty(obj.t_start)
                    h5writeatt(filename, '/time', 't_start', obj.t_start);
                end
                
                if ~isempty(obj.t_end_stamp)
                    h5writeatt(filename, '/time', 't_end_stamp', obj.t_end_stamp);
                end
                
                if ~isempty(obj.t_end)
                    h5writeatt(filename, '/time', 't_end', obj.t_end);
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
                        
            if ~isempty(obj.lightcurves)
                h5create(filename, '/lightcurves', size(obj.lightcurves,2))
                h5write(filename, '/lightcurves', obj.lightcurves);
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
                                   
            if ~isempty(obj.images_raw)

                chunk = [obj.chunk obj.chunk 1 1];
                chunk = chunk(1:ndims(obj.images_raw));
                chunk = min(chunk, size(obj.images_raw));

                S = size(obj.images_raw);

                if deflate>0
                    matlab.io.fits.setCompressionType(file_ptr,'GZIP2');
                    matlab.io.fits.setTileDim(file_ptr, chunk);
                end

                if isa(obj.images_raw, 'uint16')

                    matlab.io.fits.createImg(file_ptr, 'int16', S);
                    matlab.io.fits.writeImg(file_ptr, int16(double(obj.images_raw)-2^15));

                    matlab.io.fits.writeKey(file_ptr, 'BSCALE', 1);
                    matlab.io.fits.writeKey(file_ptr, 'BZERO', 2^15);

                elseif isa(obj.images_raw, 'double')

                    matlab.io.fits.createImg(file_ptr, 'double', S);
                    matlab.io.fits.writeImg(file_ptr, obj.images_raw);

                end

            end

            if ~isempty(positions)

            end

            if ~isempty(timestamps)

            end

            if ~isempty(obj.psfs)

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

            if ~isempty(lightcurves)

            end

            if obj.use_write_pars && ~isempty(obj.pars)
                obj.pars.writeFitsHeader(file_ptr);
            end

        end
        
        function saveMatFile(obj, filename)
                       
            if nargin<2 || isempty(filename)
                error('cannot run saveFits without a filename!');    
            end
            
            filename = util.text.extension(filename, '.mat');
            
            list_pars = {};
            
            if ~isempty(obj.images_raw)
                images_raw = obj.images_raw;
                list_pars{end+1} = 'images_raw';
            end
            
            if ~isempty(obj.cutouts_raw)
                cutouts_raw = obj.cutouts_raw;
                list_pars{end+1} = 'cutouts_raw';
            end
            
            if ~isempty(obj.positions)
                positions = obj.positions;
                list_pars{end+1} = 'positions';
            end
            
            if ~isempty(obj.full_sum)
                full_sum = obj.full_sum;
                list_pars{end+1} = 'full_sum';
            end
            
            if ~isempty(obj.num_sum)
                num_sum = obj.num_sum;
                list_pars{end+1} = 'num_sum';
            end
            
            if ~isempty(obj.timestamps)
                timestamps = obj.timestamps;
                list_pars{end+1} = 'timestamps';
            end
            
            if ~isempty(obj.t_start)
                t_start = obj.t_start;
                list_pars{end+1} = 't_start';
            end
            
            if ~isempty(obj.t_end)
                t_end = obj.t_end;
                list_pars{end+1} = 't_end';
            end
            
            if ~isempty(obj.t_end_stamp)
                t_end_stamp = obj.t_end_stamp;
                list_pars{end+1} = 't_end_stamp';
            end
            
            if ~isempty(obj.psfs)
                psfs = obj.psfs;
                list_pars{end+1} = 'psfs';
            end
            
            if ~isempty(obj.psf_sampling)
                psf_sampling = obj.psf_sampling;
                list_pars{end+1} = 'psf_sampling';
            end
            
            if ~isempty(obj.lightcurves)
                lightcurves = obj.lightcurves;
                list_pars{end+1} = 'lightcurves';
            end
            
            if obj.use_write_pars && ~isempty(obj.pars)
                pars = obj.pars;
                list_pars{end+1} = 'pars';
            end
            
            save(filename, list_pars{:});
            
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
        
        function ok = saveHDF5Static(filename, deflate, chunk, images, positions, timestamps, t_end_stamp, t_end, t_start, psfs, psf_sampling, lightcurves, pars, use_write_pars, debug_bit)
            
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
            if nargin<12, lightcurves = []; end
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
                        
            if ~isempty(lightcurves)
                h5create(filename, '/lightcurves', size(lightcurves,2))
                h5write(filename, '/lightcurves', lightcurves);
            end
                        
            if use_write_pars && ~isempty(pars)
                util.oop.save(pars, filename, 'append', 1);
            end
                        
            ok = 1;
            
        end
        
        function ok = saveFitsStatic(filename, deflate, chunk, images, positions, timestamps, t_end_stamp, t_end, t_start, psfs, psf_sampling, lightcurves, pars, use_write_pars, debug_bit)
                            
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
            if nargin<12, lightcurves = []; end
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

                if ~isempty(lightcurves)
                    
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
        
        function ok = saveMatFileStatic(filename, images, positions, timestamps, t_end_stamp, t_end, t_start, psfs, psf_sampling, lightcurves, pars, use_save_pars, debug_bit)
            
            ok = 0;
            
            if nargin<2, images = []; end
            if nargin<3, positions = []; end            
            if nargin<4, timestamps = []; end
            if nargin<5, t_end_stamp = []; end
            if nargin<6, t_end = []; end
            if nargin<7, t_start = []; end            
            if nargin<8, psfs = []; end
            if nargin<9, psf_sampling = []; end
            if nargin<10, lightcurves = []; end
            if nargin<11, pars = []; end
            if nargin<12, use_save_pars = ''; end
            if nargin<13 || isempty(debuf_bit), debug_bit = 0; end
            
            name = util.text.extension(filename, '.mat');
            
            if use_save_pars
                save(name, 'images', 'positions', 'timestamps', 'psfs', 'psf_sampling', 'lightcurves', 'pars');
            else
                save(name, 'images', 'positions', 'timestamps', 'psfs', 'psf_sampling', 'lightcurves');
            end
            
            ok = 1;
            
        end
                
    end
    
end