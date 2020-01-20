classdef Reader < file.AstroData
% Reads images (and other data) from HDF5 files. Usually used when there is
% more data than fits into memory and the processing is to be done in a
% loop. There are two main work modes:
% (1) Single file (default): when "num_files_per_batch==1" only one file is
%     loaded for every call to "batch". User can limit the number of frames
%     (and effectively load one file in multiple batches) by using
%     "num_frames_per_batch". By default it is empty, meaning it loads all
%     frames from the file. Use "frame_index_start" and "_finish" to read
%     only some of the frames in each file. Use "file_index_start" and
%     "_finish" to read only some files on the list. 
%     Use "num_batches" to limit the total number of batches to be read. 
%     Use "wrap_around" to keep reading the file list indefinitely. 
% (2) If "num_files_per_batch>1" then multiple files are loaded in each
%     batch. In this mode files cannot be split into multiple batches, but
%     limits such as "frame_index_start" and "_finish" are still respected.
%     Note that in both modes, the "frame_index" limits mean some images
%     are not read at all. 
% 
% Interface (API):
% Use "browseDir" or choose the directory yourself using the
% WorkingDirectory object. Then use "loadFiles". 
% You should use the "startup" and "finishup" commands at the beginning and
% end of each run (which also calls "loadFiles"). This conforms to the
% generic API (like the camera and simulator).
% After that (if there are any files that match the "glob_string") use
% "batch" to load the images into memory. This advances the files
% internally and automatically. 
% The "loop" command simply runs "batch" until the object is done reading
% all files (by setting "is_finished=1"), then calls "finishup". The
% command "run" also uses "startup" before the loop. 
%
% OUTPUTS: 
% -images: for older type files, but will also be a copy of what is saved
%  in any other images type (see next items) by the order: raw, cut, sum. 
%  If saving calibrated, full-frame images they could be saved here. 
% -images_raw: images straight from camera. This includes AOI images. 
% -images_cal: images after dark and flat calibration. 
% -cutouts_raw: cutouts in a 4D matrix. Without calibration.
% -cutouts_cal: calibrated cutouts. 
% -full_sum: sum of images in a batch. This should be calibrated. 
% -num_sum: how many images were put into the sum (for housekeeping). 
% -positions: a 2xN matrix of x and y positions of each cutout. 
% -timestamps: the starting time of each image in the batch (in seconds).
% -t_start: time of beginning of batch (can be Julian date or time-string).
% -t_end: time of end of batch (can be Julian date or time-string).
% -t_end_stamp: timestamp of the time of end of batch (to calibrate
%  timestamps to absolute time). 
% -psf: image of the PSF if it was measured / simulated. 
% -psf_sampling: in units of lambda/D per pixel (default is 2-> Nyquist).
% -fluxes: values of the amount of light extracted per star/cutout. 
%
% Other useful parameters incude:
% -glob_string (default= "*.h5*"): This string is used in glob-expansion to
% select only the correct files for reading. Other useful strings could be:
%  "Kraar*.h5*" | "WFAST*.h5*" | *.h5z* 
% -use_transpose: for images accidentally saved with row-columns flipped.
% -use_reshape: for images accidentally saved with wrong number of rows and
%  columns (bugs in going from 1D array in C to 3D matrix in matlab). 
% -custom_code: if user requires some additional pre-processing of images
%  from file before they are used, add matlab code in this cell array. 
% -custom_text: if custom_code needs embedded strings, pass them through
%  this cell array. 
% 
% A Graphic User Interface (GUI) for this class is also included in the 
% sub-package +gui. It is invoked using obj.makeGUI. 

    properties (Transient=true)
        
        gui@file.gui.ReaderGUI;
    
        latest_input@util.text.InputVars;
        
    end
    
    properties % object
        
        pars@head.Parameters;
        dir@util.sys.WorkingDirectory;
        
    end
    
    properties % file outputs
        
%         % make sure anything you add is added to clear...
%         images; % older files will just have generic "images" dataset
%         
%         images_raw;
%         images_cal;
%         
%         cutouts_raw;
%         cutouts_cal;
%         positions; % only for cutouts. a 2xN matrix (X then Y, N is the number of cutouts). 
%         
%         full_sum; % sum of the full frame image
%         num_sum; % if the images are summed, how many frames were added
%         
%         timestamps; % output timestamps (if available)
%         t_start; % absolute date and time (UTC) when first image is taken
%         t_end; % absolute date and time (UTC) when batch is finished
%         t_end_stamp; % timestamp when batch is finished
%                 
%         psfs; % output PSFs from file (if available)        
%         psf_sampling; % if loaded PSFs need to be binned (this is the binning factor)
%         
%         fluxes;
        
        filenames; % cell array of file names that match
        
        input_filename = '';
        
    end
    
    properties % internal counters and indices (for stopping and continuing)
        
        latest_batch_size = 100; % number of frames in each batch (default is 100, changes with latest file or by set batch limits)
        this_file_index = 1; % the next file to be read
        this_frame_index = 1; % the next frame inside the file (num_files_per_batch==1)
        counter_files = 0; % number of files that have been read completely. 
        counter_frames = 0; % number of frames read
        counter_batches = 0; % number of batches loaded
        
    end
    
    properties % switches
       
        glob_string = '*.h5*'; % match files based on this glob-expression
        
        % batch limits
        num_batches = []; % if you only want to read a few batches
        
        % file limits
        num_files_per_batch = 1; % number of files to load each batch
        file_index_start = 1; % start from this file on the list
        file_index_finish = []; % stop reading after you get to this file on the list (including)
        use_wrap_around = 0;        
        
        % frame limits
        num_frames_per_batch = []; % how many frames to read in each batch (works only when num_files_per_batch==1)
        frame_index_start = 1; % in each file, start reading from this index
        frame_index_finish = []; % in each file, finish reading on this index (or to the end of file)
        
        % Area Of Interest (to be depricated)
        AOI_left;
        AOI_width;
        AOI_top;
        AOI_height;
        
        ROI = [1 1 2560 2160];
        
        use_transpose = 0; % switches rows for columns
        use_reshape = 0; % for files saved with the wrong row/column sizes (correct data, wrong division into rows/columns)
        
        use_images_generic = 1; % if true, will put in "images" a copy of "images_cal" or "images_raw" or "images_sum" or "cutouts_cal" or "cutouts_raw" (in that order). 
        
        custom_code = {}; % after each batch, run eval on each cell. 
        custom_text = {}; % if you need to call some text in the eval, reference these cells.
        
        brake_bit = 1; % when set to 1 the reader stops the run at the end of the batch...
        
        debug_bit = 1;
                
    end
    
    properties(Dependent=true)
        
        num_files; % length of filenames
        current_dir; % string from WorkingDirectory.pwd
        this_filename; % file that is next in line to be read
        prev_filename; % file that was just now read
        
    end
    
    properties(Hidden=true)
        
        info; % lazy load this from h5info
        
        dataset_names; % a struct with fields for each possible output, containing cell arrays of strings where each string is a possible name for the dataset inside the h5 file
        attribute_names; % a struct with the names (and possible aliases) of important attributes we need to load 
        
        % defaults
        default_glob_str;
        default_num_batches; % if you only want to read a few batches
        
        % file limits
        default_num_files_per_batch; % number of files to load each batch
        default_file_index_start; % start from this file on the list
        default_use_wrap_around;        
        
        % frame limits
        default_num_frames_per_batch; % how many frames to read in each batch (works only when num_files_per_batch==1)
        default_frame_index_start; % in each file, start reading from this index
        default_frame_index_finish; % in each file, finish reading on this index (or to the end of file)
        
        version = 1.01;
        
    end
    
    properties (Hidden=true, Transient=true) % temporary controls (for custom runs like quick scan) we should move all of these to "latest_input"
        
        % temporary parameters (used for e.g. quick-scan)
        use_temp_limits = 0;
        temp_num_batches; % if you only want to read a few batches
        
        % file limits
        temp_num_files_per_batch; % number of files to load each batch
        temp_file_index_start; % start from this file on the list
        temp_file_index_finish; % stop reading after you get to this file on the list (including)
        temp_use_wrap_around;
        
        % frame limits
        temp_num_frames_per_batch; % how many frames to read in each batch (works only when num_files_per_batch==1)
        temp_frame_index_start; % in each file, start reading from this index
        temp_frame_index_finish; % in each file, finish reading on this index (or to the end of file)
        
        % temporary AOI parameters
        use_temp_AOI = 0;
        temp_AOI_left;
        temp_AOI_width;
        temp_AOI_top;
        temp_AOI_height;
                
        % temporary custom code
        use_temp_custom_code = 0;
        temp_custom_code = {};
        temp_custom_text = {};
                
    end
    
    methods % constructor
       
        function obj = Reader(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'file.Reader') % copy-constructor
                obj = util.oop.full_copy(varargin{1});
                if obj.debug_bit, fprintf('file.Reader copy-constructor v%4.2f\n', obj.version); end
            else
                
                util.oop.save_defaults(obj);
                                                
                if isempty(varargin) % default constructor
                    if obj.debug_bit, fprintf('file.Reader constructor v%4.2f\n', obj.version); end
                    obj.dir = util.sys.WorkingDirectory(getenv('DATA')); % try to start at the DATA directory (if fails, silently uses CWD). 
                elseif isa(varargin{1}, 'util.sys.WorkingDirectory') % directory constructor
                    obj.dir = util.oop.full_copy(varargin{1});
                    if obj.debug_bit, fprintf('file.Reader directory constructor v%4.2f\n', obj.version); end
                elseif ischar(varargin{1})
                    obj.dir = util.sys.WorkingDirectory(varargin{1});
                    if obj.debug_bit, fprintf('file.Reader directory constructor v%4.2f\n', obj.version); end
                else
                    error(['input to file.Reader constructor is of class "' class(varargin{1}) '". Try giving it a Reader or WorkingDirectory...']);
                end
                
                obj.initialize;
                
            end
            
        end
       
        function setupDataNames(obj)
            
            % aliases for the different data products (for backward compatibility)
            obj.dataset_names.images = {'cube','images','images_raw'};
            obj.dataset_names.images_proc = {'images_proc', 'images_cal'};
            obj.dataset_names.images = [obj.dataset_names.images, obj.dataset_names.images_proc];
            
            obj.dataset_names.timestamps = {'timestamps', 'times'};
            obj.attribute_names.t_start = {'t_start', 'file_start_datetime'};
            obj.attribute_names.t_end = {'t_end', 'file_write_datetime'};
            obj.attribute_names.t_end_stamp = {'t_end_stamp', 'file_write_timestamp'};
            
            obj.dataset_names.cutouts = {'cutouts', 'cutouts_raw', 'images_cut', 'images_stamps', 'image_cutouts'};
            obj.dataset_names.cutouts_proc = {'cutouts_proc', 'cutouts_cal'};
            obj.dataset_names.cutouts = [obj.dataset_names.cutouts, obj.dataset_names.cutouts_proc];
            
            obj.dataset_names.positions = {'positions'};
            obj.dataset_names.coordinates = {'coordinates'};
            obj.dataset_names.magnitudes = {'magnitudes'};
            obj.dataset_names.temperatures = {'temperatures'};
            
            obj.dataset_names.cutouts_bg = {'cutouts_bg'};
            obj.dataset_names.positions_bg = {'positions_bg'};
            
            obj.dataset_names.stack = {'stack', 'full_sum', 'images_sum', 'image_sum', 'summed_image', 'images_stack', 'image_stack'};
            obj.attribute_names.num_sum = {'num_sum'};
            
            obj.dataset_names.psfs = {'psfs'};
            obj.attribute_names.psf_sampling = {'psf_sampling', 'psf_binning'};
            
            obj.dataset_names.fluxes = {'fluxes', 'lightcurves'};
            
            obj.dataset_names.pars = {'pars'};
            
        end
        
    end
    
    methods % initialize/reset/clear
        
        function initialize(obj) % puts back all default values
            
            util.oop.load_defaults(obj);
            
            obj.filenames = {}; % cell array of file names that match
            
            % batch limits
            obj.num_batches = Inf; % if you only want to read a few batches
        
            % file limits
            obj.num_files_per_batch = 1; % number of files to load each batch
            obj.file_index_start = 1; % start from this file on the list
            obj.use_wrap_around = 0;        
        
            % frame limits
            obj.num_frames_per_batch = []; % how many frames to read in each batch (works only when num_files_per_batch==1)
            obj.frame_index_start = 1; % in each file, start reading from this index
            obj.frame_index_finish = Inf; % in each file, finish reading on this index (or to the end of file)
            
            obj.resetAOI;
            obj.use_transpose = 0; % switches rows for columns
            obj.use_reshape = 0;
        
            obj.custom_code = {}; % after each batch, run eval on each cell. 
            obj.custom_text = {}; % if you need to call some text in the eval, reference these cells.
        
            obj.brake_bit = 1; % when set to 1 the reader stops the run at the end of the batch...
            
%             obj.dir = util.sys.WorkingDirectory;
            obj.pars = head.Parameters;
            
            obj.setupDataNames
            
            obj.reset;
        
        end
        
        function reset(obj) % get ready for another run
            
            obj.this_file_index = obj.file_index_start; % the current file that is about to be read
            obj.this_frame_index = 1; % the next frame inside the file (for num_files_per_batch==1 only)
            obj.counter_files = 0;
            obj.counter_frames = 0;
            obj.counter_batches = 0;
            obj.latest_batch_size = 100; % default value
            
            obj.matchTempPars;
            
            obj.clear;
            
            obj.loadFiles;
            
        end
        
        function clear(obj) % get ready for another batch
            
%             obj.images = [];
%             obj.images_raw = [];
%             obj.images_cal = [];
%             obj.cutouts_raw = [];
%             obj.cutouts_cal = [];
%             obj.full_sum = [];
%             obj.positions = [];
%             
%             obj.timestamps = [];
%             obj.t_start = '';
%             obj.t_end = '';
%             obj.t_end_stamp = [];
%             
%             obj.psfs = [];
%             obj.psf_sampling = [];
%             
%             obj.fluxes = [];
            
            clear@file.AstroData(obj);

            obj.input_filename = '';
            
        end
                
        function resetAOI(obj)
           
            obj.AOI_left = [];
            obj.AOI_width = [];
            obj.AOI_top = [];
            obj.AOI_height = [];
            
        end
        
    end
    
    methods % getters
        
%         function val = get.images(obj)
%             
%             if ~isempty(obj.images)
%                 val = obj.images;
%             elseif obj.use_images_generic==0
%                 val = [];
%             elseif ~isempty(obj.images_cal)
%                 val = obj.images_cal;
%             elseif ~isempty(obj.images_raw)
%                 val = obj.images_raw;
%             elseif ~isempty(obj.full_sum)
%                 val = obj.full_sum;
%             elseif ~isempty(obj.cutouts_cal)
%                 val = obj.cutouts_cal;
%             elseif ~isempty(obj.cutouts_raw)
%                 val = obj.cutouts_raw;
%             else
%                 val = [];
%             end
%             
%         end
        
        function val = getNumBatches(obj) % an estimate of the number of batches to be loaded in this run
            
            Nfiles1 = length(obj.filenames);
            Nfiles2 = obj.num_batches*obj.temp_num_files_per_batch;
            Nfiles3 = obj.temp_file_index_finish-obj.temp_file_index_start;
            
            Nfiles = min([Nfiles1,Nfiles2,Nfiles3]); % which ever finishes first...
            
            val = ceil(Nfiles/obj.temp_num_files_per_batch);
            
        end
        
        function val = get.num_files(obj)
            
            val = length(obj.filenames);
            
        end
        
        function val = fileIndex(obj, idx)
            
            val = mod(idx-1, length(obj.filenames))+1;
            
        end
        
        function val = get.current_dir(obj)
            
            val = obj.dir.pwd;
            
        end
        
        function val = get.this_filename(obj)
            
            if isempty(obj.filenames)
                val = '';
            elseif obj.temp_use_wrap_around
                val = obj.filenames{obj.fileIndex(obj.this_file_index)};
            elseif obj.this_file_index>length(obj.filenames)
                val = '';
            elseif obj.this_file_index<1
                val = '';  
            else
                val = obj.filenames{obj.this_file_index};
            end
            
        end
        
        function val = get.prev_filename(obj)
            
            if isempty(obj.filenames)
                val = '';
            elseif obj.temp_use_wrap_around
                val = obj.filenames{obj.fileIndex(obj.this_file_index-1)};
            elseif obj.this_file_index-1>length(obj.filenames)
                val = '';
            elseif obj.this_file_index-1<1
                val = ''; 
            else
                val = obj.filenames{obj.this_file_index-1};
            end
            
        end
        
        function val = shortname(obj)
           
            [~, filename ext] = fileparts(obj.this_filename);
            
            val = [filename ext];
            
        end
        
    end
    
    methods % setters
       
        function set.dir(obj, val)
            
            obj.dir = val;
            obj.loadFiles;
            
        end
        
        function set.glob_string(obj, val)
            
            obj.glob_string = val;
            obj.loadFiles;
            
        end
        
        function matchTempPars(obj) % matches temporary parameters to user defined parameters
            
            obj.temp_num_batches = obj.num_batches;
            obj.temp_num_files_per_batch = obj.num_files_per_batch;
            obj.temp_file_index_start = obj.file_index_start;
            obj.temp_file_index_finish = obj.file_index_finish;
            obj.temp_use_wrap_around = obj.use_wrap_around;
            
            obj.temp_num_frames_per_batch = obj.num_frames_per_batch;
            obj.temp_frame_index_start = obj.frame_index_start;
            obj.temp_frame_index_finish = obj.frame_index_finish;
            
            obj.temp_AOI_left = obj.AOI_left;
            obj.temp_AOI_width = obj.AOI_width;
            obj.temp_AOI_top = obj.AOI_top;
            obj.temp_AOI_height = obj.AOI_height;
            
            obj.temp_custom_code = obj.custom_code;
            obj.temp_custom_text = obj.custom_text;
            
        end
        
        function set.num_batches(obj, val)
            
            obj.num_batches = val;
            obj.temp_num_batches = val;
            
        end
        
        function set.num_files_per_batch(obj, val)
            
            obj.num_files_per_batch = val;
            obj.temp_num_files_per_batch = val;
            
        end
        
        function set.file_index_start(obj, val)
            
            obj.file_index_start = val;
            obj.temp_file_index_start = val;
            
        end
        
        function set.file_index_finish(obj, val)
            
            obj.file_index_finish = val;
            obj.temp_file_index_finish = val;
            
        end
        
        function set.use_wrap_around(obj, val)
            
            obj.use_wrap_around = val;
            obj.temp_use_wrap_around = val;
            
        end
        
        function set.num_frames_per_batch(obj, val)
            
            obj.num_frames_per_batch = val;
            obj.temp_num_frames_per_batch = val;
            
        end
        
        function set.frame_index_start(obj, val)
            
            obj.frame_index_start = val;
            obj.temp_frame_index_start = val;
            
        end
        
        function set.frame_index_finish(obj, val)
            
            obj.frame_index_finish = val;
            obj.temp_frame_index_finish = val;
            
        end
        
        function set.AOI_left(obj, val)
            
            obj.AOI_left = val;
            obj.temp_AOI_left = val;
            
        end
        
        function set.AOI_width(obj, val)
            
            obj.AOI_width = val;
            obj.temp_AOI_width = val;
            
        end
        
        function set.AOI_top(obj, val)
            
            obj.AOI_top = val;
            obj.temp_AOI_top = val;
            
        end
        
        function set.AOI_height(obj, val)
            
            obj.AOI_height = val;
            obj.temp_AOI_height = val;
            
        end
                
        function set.custom_code(obj, val)
            
            obj.custom_code = val;
            obj.temp_custom_code = val;
            
        end
        
        function set.custom_text(obj, val)
            
            obj.custom_text = val;
            obj.temp_custom_text = val;
            
        end
        
    end
    
    methods % actions
        
        function dir = browseDir(obj)
           
            dir = obj.dir.browse;
            
            if ~isempty(dir) && ~isnumeric(dir)           
                obj.loadFiles;
            end
            
        end
        
        function loadFiles(obj)
            
            obj.filenames = obj.dir.match(obj.glob_string);

            if ~isempty(obj.filenames) && ~isempty(obj.filenames{1})
                obj.readPars; % load the parameters only for first file
            end
            
        end
        
        function advanceFile(obj)
            
            obj.info = []; % this is lazy loaded for each new file
            
            obj.this_file_index = obj.this_file_index + 1;
            
            % if we want to wrap around to the beginning
            if obj.temp_use_wrap_around
                if ~isempty(obj.temp_file_index_finish) && obj.this_file_index > obj.temp_file_index_finish
                    obj.this_file_index = obj.temp_file_index_start;
                elseif obj.this_file_index > length(obj.filenames)
                    obj.this_file_index = obj.temp_file_index_start;
                end
            end
            
            obj.this_frame_index = obj.temp_frame_index_start; % rewind the frame index, too
            
            obj.counter_files = obj.counter_files + 1;
            
        end
        
        function getAttribute(obj, filename, data_name, att_list, att_name)
            
            for jj = 1:length(att_list) % search for important attributes
            
                if any(strcmp(att_list{jj}, obj.attribute_names.(att_name)))
                    
                    obj.(att_name) = h5readatt(filename, util.text.sa('/', data_name), att_list{jj});
                    if iscell(obj.(att_name)) && ~isempty(obj.(att_name))
                        obj.(att_name) = obj.(att_name){end};
                    end
                    
                    break;
                
                end
            end
            
        end
        
        function readFile(obj)
            
            import util.text.cs;
            import util.text.sa;
            
            filename = obj.this_filename;
            num_images_on_file = 0; % by default, unless we later see there are any images
            
            if isempty(filename)
                return;
            end
            
            if ~exist(obj.this_filename,'file')
                error(['Filename: "' filename '" does not exist...']);
            end
            
            if obj.debug_bit>1
                disp(['Loading file: ' filename]);
            end
            
            obj.input_filename = filename;
            
            % check how many frames we need to read
            if obj.temp_num_files_per_batch>1 || isempty(obj.temp_num_frames_per_batch) || isinf(obj.temp_num_frames_per_batch)
                num_frames = Inf; % read the whole file
            else
                num_frames = obj.temp_num_frames_per_batch; % read only a part of the file
            end
            
            if ~isempty(obj.temp_frame_index_finish)
                num_frames = min(num_frames, obj.temp_frame_index_finish-obj.temp_frame_index_start+1);                
            end
            
            frame_start = obj.this_frame_index;
            
            % check the AOI we need to read
            if isempty(obj.temp_AOI_left), left = 1; else, left = obj.temp_AOI_left; end
            if isempty(obj.temp_AOI_width), width = Inf; else, width = obj.temp_AOI_width; end
            if isempty(obj.temp_AOI_top), top = 1; else, top = obj.temp_AOI_top; end
            if isempty(obj.temp_AOI_height), height = Inf; else, height = obj.temp_AOI_height; end
            
            [~,~,ext] = fileparts(obj.this_filename);
            
            loaded_pars = head.Parameters.empty;
            
            if cs(ext, '.h5', '.h5z', '.hdf5')
                
                if isempty(obj.info) % lazy load the info
                    try
                        obj.info = h5info(obj.this_filename);
                    catch ME
                        disp(['filename: ' obj.this_filename]);
                        rethrow(ME);
                    end
                end
                
                % maybe add a quick check that the number of frames in "images", "psfs", and "timestamps" (and "fluxes") all match?                                
                
                num_images_loaded = 0; % can be number of images or number of PSFs (if no images exist)                  
                
                for ii = 1:length(obj.info.Datasets) % go over all datasets in file, load each to the right matrix
                    
                    data_name = obj.info.Datasets(ii).Name; % the specific field we are now reading
                    data_size = obj.info.Datasets(ii).Dataspace.Size;
                    
                    if isempty(obj.info.Datasets(ii).Attributes)
                        att_names = {};
                    else
                        att_names = {obj.info.Datasets(ii).Attributes.Name};
                    end
                    
                    if any(strcmp(data_name, obj.dataset_names.images)) && all(data_size) % dataset_names.images may be a cell array of different optional names
                        
                        if length(data_size)>=3
                            num_images_on_file = data_size(3);
                        else
                            num_images_on_file = 1;
                        end
                        
                        if frame_start>num_images_on_file || (~isempty(obj.temp_frame_index_finish) && frame_start>obj.temp_frame_index_finish) % if the index of the first frame is already out of this file
                            obj.advanceFile; % this is not supposed to happen because of the "advanceFile" at the end of this function
                            return;
                        end
                        
                        num_frames = min(num_frames, num_images_on_file-obj.temp_frame_index_start+1); % make sure we don't ask for more frames than we have on file
                        start_vec = [top, left, frame_start];
                        step_vec = [height, width, num_frames];
                        if length(data_size)==2
                            start_vec = start_vec(1:2);
                            step_vec = step_vec(1:2);
                        elseif length(data_size)==4
                            start_vec = [start_vec 1];
                            step_vec = [step_vec data_size(4)];
                        end
                        
                        loaded_images = h5read(filename, sa('/', data_name), start_vec, step_vec);
                        
                        if cs(data_name, 'images_proc', 'images_cal', 8) % in the very odd case that we saved the processed images... 
                            obj.images_proc = cat(3, obj.images_proc, loaded_images); % append to the images from previous files
                        else
                            obj.images = cat(3, obj.images, loaded_images); % append to the images from previous files
                        end
                        
                        num_images_loaded = size(loaded_images,3);
                        
                    elseif any(strcmp(data_name, obj.dataset_names.cutouts_bg)) && all(data_size) % dataset_names.cutouts is a cell array of different optional names for this dataset
                        
                        num_images_on_file = data_size(3);
                        
                        if frame_start>num_images_on_file || (~isempty(obj.temp_frame_index_finish) && frame_start>obj.temp_frame_index_finish) % if the index of the first frame is already out of this file
                            obj.advanceFile; % this is not supposed to happen because of the "advanceFile" at the end of this function
                            return;
                        end
                        
                        num_frames = min(num_frames, num_images_on_file-obj.temp_frame_index_start+1); % make sure we don't ask for more frames than we have on file
                        start_vec = [top, left, frame_start];
                        step_vec = [height, width, num_frames];
                        
                        if length(data_size)==4
                            start_vec = [start_vec 1];
                            step_vec = [step_vec data_size(4)];
                        end
                        
                        loaded_cutouts = h5read(filename, sa('/', data_name), start_vec, step_vec);
                        
                        obj.cutouts_bg = cat(3, obj.cutouts_bg, loaded_cutouts); % append to the images from previous files
                        
                        if ~num_images_loaded % if the file doesn't have images but does have cutouts, count the cutout frames
                            num_images_loaded = size(loaded_cutouts,3);
                        end
                                                
                    elseif any(strcmp(data_name, obj.dataset_names.positions_bg)) && all(data_size) % data_names.positions may be a cell array of different optional names
                        
                        loaded_positions = h5read(filename, sa('/', data_name)); % read the entire "positions_bg" dataset, if it exists
                        obj.positions_bg = cat(3, obj.positions_bg, loaded_positions);
                        
                    elseif any(strcmp(data_name, obj.dataset_names.cutouts)) && all(data_size) % dataset_names.cutouts is a cell array of different optional names for this dataset
                        
                        num_images_on_file = data_size(3);
                        
                        if frame_start>num_images_on_file || (~isempty(obj.temp_frame_index_finish) && frame_start>obj.temp_frame_index_finish) % if the index of the first frame is already out of this file
                            obj.advanceFile; % this is not supposed to happen because of the "advanceFile" at the end of this function
                            return;
                        end
                        
                        num_frames = min(num_frames, num_images_on_file-obj.temp_frame_index_start+1); % make sure we don't ask for more frames than we have on file
                        start_vec = [top, left, frame_start];
                        step_vec = [height, width, num_frames];
                        
                        if length(data_size)==4
                            start_vec = [start_vec 1];
                            step_vec = [step_vec data_size(4)];
                        end
                        
                        loaded_cutouts = h5read(filename, sa('/', data_name), start_vec, step_vec);
                        
                        obj.cutouts = cat(3, obj.cutouts, loaded_cutouts); % append to the images from previous files
                        
                        if ~num_images_loaded % if the file doesn't have images but does have cutouts, count the cutout frames
                            num_images_loaded = size(loaded_cutouts,3);
                        end
                                                
                    elseif any(strcmp(data_name, obj.dataset_names.positions)) && all(data_size) % data_names.positions may be a cell array of different optional names
                        
                        loaded_positions = h5read(filename, sa('/', data_name)); % read the entire "positions" dataset, if it exists
                        obj.positions = loaded_positions(:,:,end); % to handle cases where we accidentally saved coordinates and positions together... 
                        
                    elseif any(strcmp(data_name, obj.dataset_names.coordinates)) && all(data_size) % data_names.coordinates may be a cell array of different optional names
                        
                        obj.coordinates = h5read(filename, sa('/', data_name)); % read the entire "cordinates" dataset, if it exists
                        
                    elseif any(strcmp(data_name, obj.dataset_names.magnitudes)) && all(data_size) % data_names.magnitudes may be a cell array of different optional names
                        
                        obj.magnitudes = h5read(filename, sa('/', data_name)); % read the entire "magnitudes" dataset, if it exists
                        
                    elseif any(strcmp(data_name, obj.dataset_names.temperatures)) && all(data_size) % data_names.temperatures may be a cell array of different optional names
                        
                        obj.temperatures = h5read(filename, sa('/', data_name)); % read the entire "temperatures" dataset, if it exists
                        
                    elseif any(strcmp(data_name, obj.dataset_names.stack)) && all(data_size)
                        
                        loaded_stack = h5read(filename, sa('/', data_name)); % if there are any stack images (or "full_sum" images) just get them from the file. 
                        obj.stack = cat(3, obj.stack, loaded_stack); % append to previous file's results. 
                        
                        obj.getAttribute(filename, data_name, att_names, 'num_sum'); % need to add this to the "stack" property
                        
                    elseif any(strcmp(data_name, obj.dataset_names.timestamps)) && all(data_size) % data_names.timestamps may be a cell array of different optional names
                        
                        num_images_on_file = max(data_size); % timestamps should be 1D vector or 2D 1xN vector
                        
                        num_frames = min(num_frames, num_images_on_file); % make sure we don't ask for more frames than we have on file
                        
                        if length(data_size)>1
                            frame_start = [frame_start 1];
                            num_frames = [num_frames 1];
                        end
                        
                        loaded_timestamps = h5read(filename, sa('/', data_name), frame_start, num_frames);
                        obj.timestamps = cat(1, obj.timestamps, loaded_timestamps); % append to the existing images
                        
                        obj.getAttribute(filename, data_name, att_names, 't_start');
                        obj.getAttribute(filename, data_name, att_names, 't_end');
                        obj.getAttribute(filename, data_name, att_names, 't_end_stamp');
                        
                    elseif any(strcmp(data_name, obj.dataset_names.psfs)) && all(data_size) % data_names.psfs may be a cell array of different optional names
                        
                        num_images_on_file = data_size(3);
                        % do we need a check for frame_start>num_images_on_file? maybe need to verify no images exist??
                        
                        num_frames = min(num_frames, num_images_on_file); % make sure we don't ask for more frames than we have on file
                        start_vec = [top, left, frame_start];
                        step_vec = [height, width, num_frames];
                        
                        loaded_psfs = h5read(filename, sa('/', data_name), start_vec, step_vec);
                        obj.psfs = cat(3, obj.psfs, loaded_psfs); % append to the existing images
                        
                        obj.getAttribute(filename, data_name, att_names, 'psf_sampling');
                                                
                        if ~num_images_loaded % if the file doesn't have images but does have PSFs, count the PSF frames
                            num_images_loaded = size(loaded_psfs,3);
                        end
                        
                    elseif any(strcmp(data_name, obj.dataset_names.fluxes)) && data_size(1) % data_names.fluxes may be a cell array of different optional names
                        
                        num_images_on_file = data_size(1); % fluxes should be 2D with number of frames in the 1st dim
                        
                        num_frames = min(num_frames, num_images_on_file); % make sure we don't ask for more frames than we have on file
                        
                        loaded_fluxes = h5read(filename, sa('/', data_name), [frame_start 1], [num_frames Inf]); % must check the dimensions on file fit what I think I am saving...
                        obj.fluxes = cat(1, obj.fluxes, loaded_fluxes); % append to the existing images
                        
                    elseif any(strcmp(data_name, obj.dataset_names.pars)) 
                        obj.loadParsHDF5(filename, loaded_pars, data_name, att_names);
                    end
                    
                end % scan all datasets
                
                for ii = 1:length(obj.info.Groups) % go over all groups in file, looking for a "pars" object
                    
                    group_name = obj.info.Groups(ii).Name; % the specific field we are now reading
                    
                    if isempty(obj.info.Groups(ii).Attributes)
                        att_names = {};
                    else
                        att_names = {obj.info.Groups(ii).Attributes.Name};
                    end
                    
                    if any(strcmp(group_name, strcat('/', obj.dataset_names.pars)))
                        obj.loadParsHDF5(filename, loaded_pars, group_name, att_names);
                    end
                    
                end
                
                obj.this_frame_index = obj.this_frame_index + num_images_loaded; % update the frame index with the actual numbe of frames read
                obj.counter_frames = obj.counter_frames + num_images_loaded;
                
                % check if file is done
                if obj.this_frame_index>=num_images_on_file || (~isempty(obj.temp_frame_index_finish) && frame_start(1)>obj.temp_frame_index_finish)
                    obj.advanceFile;
                end
                
            elseif cs(ext, '.mat')
                error('Reading MAT files is still not implemented!');
            elseif cs(ext, '.fits')
                error('Reading FITS files is still not implemented!');
            else
                error(['Filename ' obj.filenames{obj.this_file_index} ' has unknown format... use HDF5 or mat']);
            end
        end
        
        function readPars(obj)
            
            loaded_pars = util.oop.load(obj.filenames{1}, 'location', '/pars', 'class', class(obj.pars));
            util.oop.copy_props(obj.pars, loaded_pars);  % make a shallow copy (sub-objects are referenced, then loaded_pars is destroyed)
            
        end
        
        function loadParsHDF5(obj, filename, loaded_pars, data_name, att_names)
            
            import util.text.cs;
            import util.text.sa;
            
            location = sa('/', data_name);
            
            try
                loaded_pars = util.oop.load(filename, 'location', location, 'class', class(obj.pars));

                loaded_pars.ephem.time = util.text.str2time(loaded_pars.STARTTIME); % fix the bug in read/write of datetime objects we used to have (only rely on times stored as strings)
                loaded_pars.ephem.updateSecondaryCoords;
            catch ME
                % if we can't read this it is OK, there is still the README file. Can't give warnings on every file now...
                disp('Could not read parameter object...');
                warning(getReport(ME));
            end
            
            if ~isempty(loaded_pars)
                util.oop.copy_props(obj.pars, loaded_pars); % make a shallow copy (sub-objects are referenced, then loaded_pars is destroyed)
                
                % additional loading of attributes (mostly Dependent) for backward compatibility with older files.
                props = {'RA_DEG', 'DEC_DEG'};
                
                for ii = 1:length(att_names)                    
                    for jj = 1:length(props)                        
                        if cs(att_names{ii}, props{jj})
                            try 
                                loaded_pars.(props{jj}) = h5readatt(filename, location, att_names{ii});
                                if iscell(loaded_pars.(props{jj}))
                                    loaded_pars.(props{jj}) = loaded_pars.(props{jj}){1};
                                end
                            end
                        end
                    end
                end
                
            end
            
            loaded_pars.ephem.time = util.text.str2time(loaded_pars.STARTTIME); % fix the bug in read/write of datetime objects we used to have (only rely on times stored as strings)
            loaded_pars.ephem.updateSecondaryCoords;
            
        end
        
        function batch(obj)
            
            if obj.is_finished
                obj.brake_bit = 1;
                return;
            end
            
            obj.clear;
            
            frames_now = obj.counter_frames; % how many frames were there at the start of the batch
            
            for ii = 1:length(obj.num_files_per_batch)
                obj.readFile;
                drawnow;
            end
            
            obj.counter_batches = obj.counter_batches + 1;
            
            num_frames_in_batch = obj.counter_frames-frames_now;
            
            if num_frames_in_batch
                obj.latest_batch_size = num_frames_in_batch;
            end
            
            if obj.gui.check
                obj.show;
            end
            
        end
        
        function loop(obj) % continue to loop from where you stopped
        
            while ~obj.is_finished
                
                obj.batch;
                
                drawnow;
                if obj.brake_bit
                    break;
                end
                
            end
            
            obj.finishup;
            
        end
        
        function startup(obj, varargin) % begin a new run (or continue the same run)
            
            obj.latest_input = obj.makeInputVars(varargin{:});
            
%             obj.reset;
            obj.brake_bit = 0;
            
            % read pars from text file? 
            
        end
        
        function finishup(obj) % finish a run
            
            obj.brake_bit = 1;
            % anything else??
            
        end
            
        function run(obj) % calls startup, loop and finishup
            
            if obj.num_files==0
                disp(['No files found in dir: "' obj.current_dir '" with glob= "' obj.glob_string '"']);
                return;
            end
            
            obj.startup;
                        
            obj.loop;
                                    
        end
        
        function quickScan(obj) % reads one frame from each file
            
            if obj.num_files==0
                disp(['No files found in dir: "' obj.current_dir '" with glob= "' obj.glob_string '"']);
                return;
            end
            
            obj.startup;
            
            obj.temp_num_files_per_batch = 1;
            obj.temp_frame_index_start = 1;
            obj.temp_frame_index_finish = 1;
            
            
            obj.loop;
            
        end
        
        function convert2fits(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('RA', [], 'right ascention');
            input.input_var('DE', [], 'declination');
            input.input_var('name', '', 'target_name');
            input.input_var('filter', '');
            input.input_var('num_batches', []);
            input.input_var('folder', [], 'directory');
            input.scan_vars(varargin{:});
            
            if isempty(input.folder)
                input.folder = pwd;
            elseif isa(input.folder, 'util.sys.WorkingDirectory')
                input.folder = input.folder.pwd;
            end
            
            obj.startup;
            
            N = obj.num_batches;
            if isinf(N), N = 1e6; end
            
            for ii = 1:N
                
                if obj.brake_bit || obj.is_finished
                    break;
                end
                
                obj.batch;
                
                if ~isempty(input.RA)
                    obj.pars.RA = input.RA;
                end

                if ~isempty(input.DE)
                    obj.pars.DE = input.DE;
                end

                if ~isempty(input.name)
                    obj.pars.target_name = input.name;
                end
                
                if ~isempty(input.filter)
                    obj.pars.filter = input.filter;
                end
                
                [~, name, ~] = fileparts(obj.prev_filename);
                
                ext = '.fit';
                
                for jj = 1:size(obj.images,3)
                    
                    I = obj.images(:,:,jj);
                    I = int32(I);
%                     I = I - 2^15;
%                     I = int16(I);
                    
                    if size(obj.images,3)==1
                        filename = fullfile(input.folder, sprintf('%s%s', name, ext));
                    elseif size(obj.images,3)<=10
                        filename = fullfile(input.folder, sprintf('%s_%d%s', name, jj-1, ext));
                    elseif size(obj.images,3)<=100
                        filename = fullfile(input.folder, sprintf('%s_%02d%s', name, jj-1, ext));
                    elseif size(obj.images,3)<=1000
                        filename = fullfile(input.folder, sprintf('%s_%03d%s', name, jj-1, ext));
                    end
                    
                    if obj.debug_bit, fprintf('Writing FITS file: %s\n', filename); end
                    
                    fitswrite(I, filename, 'WriteMode', 'overwrite', 'Compression', 'gzip');
                    
                    obj.pars.writeFITS(filename, obj.timestamps(jj)-obj.timestamps(1));
                    
                    drawnow;
                    
                end
                
            end
            
            obj.finishup;
            
        end
        
        function val = is_finished(obj)
                        
            if obj.counter_batches >= obj.num_batches
                val = 1;
            elseif obj.temp_use_wrap_around==0
                if obj.this_file_index>length(obj.filenames) || (~isempty(obj.temp_file_index_finish) && obj.this_file_index>obj.temp_file_index_finish)
                    val = 1;
                else
                    val = 0;
                end                
            else
                val = 0;
            end
            
        end
        
        function input = makeInputVars(obj, varargin)
            
            idx = util.text.InputVars.isInputVars(varargin);
            
            if ~isempty(varargin) && any(idx)
                idx = find(idx, 1, 'first');
                input = varargin{idx}; 
                varargin(find(idx, 1, 'first')) = []; % remove the one cell with InputVars in it
            else
                input = util.text.InputVars;
                input.use_ordered_numeric = 1;
                input.input_var('num_batches', obj.num_batches, 'Nbatches');
                % add other options??
            end
            
            input.scan_vars(varargin{:});
            
        end
        
    end
    
    methods % GUI and plotting tools
       
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = file.gui.ReaderGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
        function show(obj, varargin)
            
            if ~isempty(obj.images) || ~isempty(obj.stack)
                if obj.gui.check
                    
                    if ~isempty(obj.images)
                        I = obj.images(:,:,1);
                    elseif ~isempty(obj.stack)
                        I = obj.stack;
                    end
                    
                    im = findobj(obj.gui.axes_image, 'type', 'Image');
                    if isempty(im)
                        imagesc(obj.gui.axes_image, I);
                        axis(obj.gui.axes_image, 'image');
                        colorbar(obj.gui.axes_image, 'on');
                    else
                        im.CData = I;
                    end
                    
                    obj.gui.update;
                    
                end
            end
            
        end
        
    end
    
end