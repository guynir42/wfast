classdef Deflator < file.AstroData
% This class brings together a file.Reader and file.BufferWheel, and loops 
% over a list of files and folders, reading them from one location, and 
% writing them to another location, after deflating them. 
% This is done because we currently cannot deflate and save the files in 
% asynchronuous mode. The idea is to run deflation during the day. 
%
% NOTE: deflator skips HDF5 files that already exist in the destination/backup
%       folder, so if you accidentally deflate the same folder twice, it only
%       re-copies the catalog, readme files and calibration files. 
% 
% Using the Deflator
%--------------------
% Manual control: Use the GUI or select the working directories to the places
% where you are storing the source and output (and backup) files. 
% You can use the WorkingDirectory's own cd() or browse() methods. 
% Once these are selected, simply use run() or run_async() to loop over all 
% the folders in "src_dir". 
% 
% Automated start: When the Acquisition class starts a run it sets up the 
% deflator to a delayed start using setup_timer(hour=5), where hour is the 
% time of day to start deflation (in UTC, whereas 5 is 7am or 8am). 
% In this case you should make sure that "auto_deflate_destination" and
% "auto_deflate_backup_dir" are all set to proper folders (these are strings
% and not WorkingDirectory objects!). 
% Also choose "use_parallel_auto_deflate" for doing the deflation using 
% run() or run_async() (default is true --> use async). 
% NOTE: this goes over all data in DATA_TEMP, but the deflator skips files 
% that are already present in the destination. 
% 
% OTHER OPTIONS:
%----------------
% -use_auto_backup: make a third copy of the data (using matlab's copyfile).
%                   This is stored to the "out_dir_backup" or subfolders of
%                   the "auto_deflate_backup_dir". Default is true. 
% 
% -use_auto_delete: delete the original files after deflation. Default false. 
%
% -use_copy_text_files or _catalogs or _calibration: These are used to copy
%                 over the auxiliary files generated with the data. 
%                 The default is true for all three.
% 
% -auto_deflate_destination and _backup_dir: these are strings specifying the
%              destination for autodeflate calls. This can be the same as 
%              the main data folder but that would mean the Dropbox would
%              need to be constantly uploading files (we can't do that). 
%              Instead these are set, by default, to environmental variables
%              DATA_EXTRAS and DATA_BACKUP, that should be on two separate 
%              drives, outside the Dropbox folder. The backup should be to 
%              a non-removable HDD. 
%              Deflator will automatically make date-folders at those destinations
%              for any existing date folders it finds inside DATA_TEMP. 
%
% -use_parallel_auto_deflate: when the automatic deflate timer goes off in 
%                             the morning, use the async mode (Default true). 
%
% NOTE: the HDF5 files are deflated and then their extension changes from 
%       .h5 to .h5z to show they have been compressed. 

    
    properties (Transient=true)
        
        gui@file.gui.DeflatorGUI;
        
        prog@util.sys.ProgressBar;
        
        futures = {}; % for using parfeval
        
        timer;
        
    end
    
    properties % objects and resources
        
        head@head.Header; % link back to owner's header (also shared with reader and buffers)
        
        src_dir@util.sys.WorkingDirectory; % point this to a date-folder with the data you wish to deflate
        out_dir@util.sys.WorkingDirectory; % point this to a new date-folder you generated and wish to deflate into
        out_dir_backup@util.sys.WorkingDirectory; % point this to a new date-folder you generated and wish to copy the deflated files into
        
        reader@file.Reader; % this object just reads the HDF5 files
        buffers@file.BufferWheel; % this is used to deflate and write the data into new files
        
    end
    
    properties % switches and controls
        
        auto_deflate_destination = getenv('DATA_EXTRAS'); % automatically create date-folders inside this destination (should be on removable HDD)
        auto_deflate_backup_dir = getenv('DATA_BACKUP'); % automatically create date-folders inside this destination (should be on permanent HDD, outside Dropbox)
        
        use_parallel_auto_deflate = 1; % when setting up a delayed, automatic defaltion, use run_asyn() instead of run()
        
        use_auto_delete = 0; % delete files from source after they are successfully deflated
        
        use_auto_backup = 1; % use matlab's copyfile to make a second copy of newly deflated files
        
        use_copy_text_files = 1; % make sure to copy the readme text files from each folder
        use_copy_catalogs = 1; % make sure to copy the catalog MAT-files from each folder
        use_copy_calibration = 1; % make sure to copy the calibration MAT-files from each folder
        use_copy_lightcurves = 1; % make sure to copy the stack lightcurves MAT-files from each folder
        use_copy_microflares = 1; % make sure to copy the micro flares summary MAT-files from each folder
        debug_bit = 1;
        
    end
    
    properties(Hidden=true)
        
        src_subdir@util.sys.WorkingDirectory; % this is automatically pointed at subfolder of src_dir
        out_subdir@util.sys.WorkingDirectory; % this is automatically pointed at subfolder of out_dir
        out_subdir_backup@util.sys.WorkingDirectory; % this is automatically pointed at subfolder of out_dir_backup
        
        brake_bit = 1;
        
        version = 1.02;
                
    end
    
    methods % constructor
        
        function obj = Deflator(other)
            
            if nargin==0 || isempty(other)
                if obj.debug_bit>1, fprintf('Deflator constructor v%4.2f\n', obj.version); end
                
                obj.setupDirs;
                obj.setupBuffers;
                obj.setupHeader;
                obj.setupReader;
                obj.setupProgressBar;
                
            end
            
        end
        
        function setupHeader(obj)
           
            obj.head = head.Header;
            
        end
        
        function setupDirs(obj, source, output, backup)
            
            if nargin<2 || isempty(source)
                source = getenv('DATA_TEMP');
                if isempty(source)
                    source = getenv('DATA');
                end
            end
            
            if nargin<3 || isempty(output)
                output = getenv('DATA_EXTRAS');
                if isempty(output)
                    output = getenv('DATA');
                end
            end
            
            if nargin<4 || isempty(backup)
                backup = getenv('DATA_BACKUP'); 
                if isempty(backup)
                    backup = getenv('DATA');
                end
            end
            
            obj.src_dir = util.sys.WorkingDirectory(source);
            obj.src_subdir = util.sys.WorkingDirectory(source);
            
            obj.out_dir = util.sys.WorkingDirectory(output);            
            obj.out_subdir = util.sys.WorkingDirectory(output);
            
            obj.out_dir_backup = util.sys.WorkingDirectory(backup);
            obj.out_subdir_backup = util.sys.WorkingDirectory(backup);
            
        end
        
        function setupBuffers(obj, num_buffers, deflate_level)
            
            if nargin<2 || isempty(num_buffers)
                num_buffers = 5;
            end
            
            if nargin<3 || isempty(deflate_level)
                deflate_level = 1;
            end
            
            obj.buffers = file.BufferWheel(num_buffers, obj.head);
            
            obj.buffers.use_deflate = deflate_level;
            
            obj.buffers.use_async = 0;
            obj.buffers.debug_bit = 1;
                        
        end
        
        function setupReader(obj)
            
            obj.reader = file.Reader;
            
        end
        
        function setupProgressBar(obj)
            
            obj.prog = util.sys.ProgressBar;
            
        end
        
    end
    
    methods % reset
        
        function clear(obj)
           
            clear@file.AstroData(obj);
            obj.clearBuffers;
            
        end
        
        function clearBuffers(obj)
        
            obj.buffers.clear;         
            
        end
        
        function resetFilenames(obj)
            
            obj.reader.dir = util.sys.WorkingDirectory(obj.src_subdir);
            obj.reader.reset;
            obj.reader.num_files_per_batch = 1;
            
        end
        
        function reset_futures(obj)
            
            obj.futures = {};
            
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
        function set.head(obj, header) % make sure all subobjects (reader/buffers) share this header
            
            obj.head= header;

            obj.buffers.head= header;
                        
        end
        
    end
    
    methods % API
        
        function run(obj, varargin) % loop over all date-folders inside src_dir and make deflated copies in out_dir
        % Usage: run(obj, varargin)
        % Find all folders inside the src_dir and make copies of them into 
        % the out_dir (possibly making a second copy to out_dir_backup). 
        % 
        % Will copy over text and MAT-files but will read HDF5 files and 
        % write them back to disk after deflation is done. 
        % Note that this is done using BufferWheel's use_async turned off. 
        % 
        % Ideally this should be called after setting the src_dir and out_dir
        % to a date folder in the DATA_TEMP and in the DATA_EXTRAS folders. 
        % 
        % Optional Arguments:
        %---------------------
        % out_dir: supply a full path to the output directory (making it if
        %          it does not exist yet. 
        % src_dir: supply a full path to the source directory. 
        %
            
            import util.text.sa;
            
            input = util.text.InputVars;
            input.input_var('out_dir', '', 'output_directory');
            input.input_var('src_dir', '', 'source_directory')
            input.scan_vars(varargin{:});
            
            if ~isempty(input.src_dir)
                
                if ~exist(input.src_dir, 'dir')
                    error('Source directory %s does not exist!', input.src_dir);
                end
                
                obj.src_dir.cd(pwd);
                obj.src_dir.cd(input.src_dir);
                
            end
            
            if ~isempty(input.out_dir)
                
                if ~exist(input.out_dir, 'dir')
                    mkdir(input.out_dir);
                end
                
                obj.out_dir.cd(pwd);
                obj.out_dir.cd(input.out_dir);
                
            end
            
            obj.brake_bit = 0;
            
            try
                
                d = obj.src_dir.dir;
                d{end+1} = ''; % add the base dir itself to the list... 
                
                N = 0;
                
                for jj = 1:length(d) % estimate number of files
                    obj.src_subdir.cd(sa(obj.src_dir.pwd, d{jj}));
                    N = N + length(obj.src_subdir.match('*.h5'));
                end
                
                if isempty(obj.prog)
                    obj.prog = util.sys.ProgressBar;
                end
                
                obj.prog.start(N);
                
                for jj = 1:length(d)

                    obj.src_subdir.cd(sa(obj.src_dir.pwd, d{jj}));
                    
                    temp_dir = sa(obj.out_dir.pwd, d{jj});

                    if ~exist(temp_dir, 'dir')
                        mkdir(temp_dir);
                    end
                    
                    obj.out_subdir.cd(temp_dir);

                    if obj.use_auto_backup

                        temp_dir = sa(obj.out_dir_backup.pwd, d{jj});

                        if ~exist(temp_dir, 'dir')
                            mkdir(temp_dir);
                        end
                        
                        obj.out_subdir_backup.cd(temp_dir);

                    end
                    
                    util.text.date_printf('Going over directory: %s.', obj.src_subdir.pwd);
                    
                    if obj.brake_bit
                        return;
                    end

                    obj.resetFilenames;
                    
                    if obj.use_copy_text_files
                        f = obj.src_subdir.match('*.txt');
                        for ii = 1:length(f)
                            
                            [~, name, ext] = fileparts(f{ii});
                            util.text.date_printf('Copying file %s%s from %s to %s.', name, ext, obj.src_subdir.pwd, obj.out_subdir.pwd);
                            copyfile(sa(obj.src_subdir.pwd,name,ext), sa(obj.out_subdir.pwd, name, ext));
                            
                            if obj.use_auto_backup
                                util.text.date_printf('Copying file %s%s from %s to %s (backup).', name, ext, obj.src_subdir.pwd, obj.out_subdir_backup.pwd);
                                copyfile(sa(obj.src_subdir.pwd,name,ext), sa(obj.out_subdir_backup.pwd, name, ext));
                            
                            end
                            
                        end
                    end

                    if obj.use_copy_catalogs
                        f = obj.src_subdir.match('catalog.mat');
                        if ~isempty(f)
                            [~, name, ext] = fileparts(f{1});
                            util.text.date_printf('Copying file %s from %s to %s.', name, ext, obj.src_subdir.pwd, obj.out_subdir.pwd);
                            copyfile(sa(obj.src_subdir.pwd,name,ext), sa(obj.out_subdir.pwd, name, ext));
                        end
                    end
                    
                    if obj.use_copy_calibration
                        f = obj.src_subdir.match('calibration*.mat');
                        if ~isempty(f)
                            [~, name, ext] = fileparts(f{1});
                            util.text.date_printf('Copying file %s%s from %s to %s.', name, ext, obj.src_subdir.pwd, obj.out_subdir.pwd);
                            copyfile(sa(obj.src_subdir.pwd,name,ext), sa(obj.out_subdir.pwd, name, ext));
                        end
                    end
                    
                    if obj.use_copy_lightcurves
                        f = obj.src_subdir.match('lightcurves*.mat');
                        if ~isempty(f)
                            [~, name, ext] = fileparts(f{1});
                            util.text.date_printf('Copying file %s%s from %s to %s.', name, ext, obj.src_subdir.pwd, obj.out_subdir.pwd);
                            copyfile(sa(obj.src_subdir.pwd,name,ext), sa(obj.out_subdir.pwd, name, ext));
                        end
                    end
                    
                    if obj.use_copy_microflares
                        f = obj.src_subdir.match('micro_flares*.mat');
                        if ~isempty(f)
                            [~, name, ext] = fileparts(f{1});
                            util.text.date_printf('Copying file %s%s from %s to %s.', name, ext, obj.src_subdir.pwd, obj.out_subdir.pwd);
                            copyfile(sa(obj.src_subdir.pwd,name,ext), sa(obj.out_subdir.pwd, name, ext));
                        end
                    end
                    
                    obj.buffers.serial = 1;
                    
                    for ii = 1:100000 % go over image files

                        if obj.brake_bit
                            return;
                        end

                        if obj.reader.is_finished
                            break;
                        end
                        
                        obj.clear;
                        
                        origin_filename = obj.reader.this_filename;
                        [output_exists, output_filename] = obj.checkFileExists(origin_filename, obj.out_subdir.pwd);
                        
                        if ~output_exists % load and write to output
                            obj.batch;
                        else % file exists, check if we need to back it up too
                            obj.reader.advanceFile; % if this file exists we don't need to re-read it
                        end
                        
                        if obj.use_auto_backup 
                            
                            [backup_exists, backup_filename] = obj.checkFileExists(origin_filename, obj.out_subdir_backup.pwd);
                        
                            if ~backup_exists % must make a backup of this file
                                
                                if ~exist(obj.out_subdir_backup.pwd, 'dir')
                                    mkdir(obj.out_subdir_backup.pwd);
                                end
                                
                                copyfile(output_filename, backup_filename);
                                
                            end
                            
                        end

                        obj.prog.advance;
                        obj.prog.showif;
                        
                        if ~isempty(obj.gui) && obj.gui.check
                            obj.gui.updateGUI;
                        end
                        
                    end % for ii (files)

                end % for jj (directories)
            
            catch ME
                obj.brake_bit = 1;
                rethrow(ME);
            end
                
            obj.brake_bit = 1;
            
            obj.prog.finish;
            
            util.text.date_printf('done deflating files...');
            
        end
        
        function run_async(obj, varargin) % same as run(), only it works on a parallel worker
        
            obj.futures{end+1} = parfeval(@obj.run, 0, varargin{:}); 
            
        end
        
        function deleteTempFiles(obj, folder, backup_folder)
            
            d_temp = getenv('DATA_TEMP'); 
            l_temp = length(d_temp);
            
            if isempty(d_temp)
                error('Cannot find the DATA_TEMP folder! Set the environmental variable...'); 
            end
            
            if nargin<2 || isempty(folder)
                folder = d_temp; 
            end
            
            if nargin<3 || isempty(backup_folder)
                backup_folder = obj.auto_deflate_destination;
            end
            
            if length(folder)<l_temp || ~strcmp(d_temp, folder(1:l_temp)) % folder must start with DATA_TEMP address
                folder = fullfile(d_temp, folder); 
            end
            
            if ~exist(folder, 'dir')
                error('Could not find the folder "%s". ', folder); 
            end
            
            if ~exist(backup_folder, 'dir')
                error('Could not find the backup folder "%s". ', backup_folder); 
            end
            
            d = util.sys.WorkingDirectory(folder); 
            if isempty(d.dir)
                return;
            end
            
            % check that deflator futures are all done
            for ii = 1:length(obj.futures)
                if strcmp(obj.futures{ii}.State, 'running')
                    return;
                end
            end

            % prompt the user for confirmation (at least for now)
%             rep = questdlg('Ready to delete temp data files', 'Delete files?', 'Yes', 'No', 'Yes');
%             
%             if isempty(rep) || strcmp(rep, 'No')
%                 return;
%             end
            
            if obj.debug_bit, util.text.date_printf('Deleting old temporary files.'); end
            
            % run auto-delete!
            list = flip(d.walk);
            
            for ii = 1:length(list)
                
                if strcmp(d_temp, list{ii}) || strcmp([d_temp '\'], list{ii}) || strcmp([d_temp '/'], list{ii})
                    continue; % skip the root DATA_TEMP
                end 
                
                d.cd(list{ii}); 
                
                files = d.files('', 1);
                
                d_bkp = list{ii}(l_temp+1:end);
                d_bkp = fullfile(backup_folder, d_bkp);
                
                for jj = 1:length(files)
                    
                    if obj.checkFileExists(files{jj}, d_bkp)
                        
                        if obj.debug_bit>1, util.text.date_printf('Found copy of file "%s" in "%s". Deleting it now!', files{jj}, d_bkp); end
                        
                        delete(files{jj}); 
                        
                    end
                    
                end
                
                if d.is_empty
                    if obj.debug_bit, util.text.date_printf('Deleting empty folder "%s".', d.pwd); end
                    rmdir(d.pwd); 
                end
                
            end
            
        end
        
        function setup_timer(obj, hour) % sets up the deflator to automatically deflate latest files in DATA_TEMP at the given hour today (UTC!)
            
            if nargin<2 || isempty(hour)
                hour = 5; % equal to 7 or 8am local time
            end
            
            if isempty(getenv('DATA_TEMP'))
                warning('Cannot auto-deflate without environmental variable DATA_TEMP!');
                return;
            end
            
            if ~isempty(obj.timer) && isvalid(obj.timer)
                stop(obj.timer); 
                delete(obj.timer);
            end
            
            obj.timer = timer('TimerFcn', @obj.callback_timer); 
            obj.timer.Name = 'auto-deflate-timer';
            obj.timer.BusyMode = 'queue';
            
            t = datetime(util.sys.date_dir('now'), 'TimeZone', 'UTC');
            
            t = t + days(1); % deflate on the next morning!
            t.Hour = hour;
            t.TimeZone = 'local'; % must convert to local time before giving it to startat()
            
            startat(obj.timer, [t.Year, t.Month, t.Day, t.Hour, 0, 0]);
            
        end
        
        function callback_timer(obj, ~, ~) % this is called when the timer triggers
            
            util.text.date_printf('Running auto-deflate!');
            
            obj.makeGUI;
            
            d = util.sys.WorkingDirectory(getenv('DATA_TEMP')); 
%             list = d.dir('',1);
%             d.cd(list{end});
            list = d.dir('',1);
            
            for ii = 1:length(list)
                
                obj.src_dir.cd(list{ii});
                
                dest = fullfile(obj.auto_deflate_destination, obj.src_dir.tail);
                bkp = fullfile(obj.auto_deflate_backup_dir, obj.src_dir.tail); 
                
                if ~exist(dest, 'dir')
                    mkdir(dest);
                end
                
                obj.out_dir.cd(dest);
                
                fprintf('\nSRC DIR: %s\n\nDEST DIR: %s', obj.src_dir.pwd, obj.out_dir.pwd);
                
                if obj.use_auto_backup && ~isempty(obj.auto_deflate_backup_dir)
                    
                    bkp = fullfile(obj.auto_deflate_backup_dir, obj.src_dir.tail); 
                
                    if ~exist(bkp, 'dir')
                        mkdir(bkp);
                    end

                    obj.out_dir_backup.cd(bkp); 
                    
                    fprintf(' (backup: %s)', bkp); 
                    
                end
                
                fprintf('\n\n'); 
                
                if obj.use_parallel_auto_deflate
                    obj.run_async;
                else
                    obj.run; 
                end
                
            end
            
            
        end
        
        function makeGUI(obj)
           
            if isempty(obj.gui)
                obj.gui = file.gui.DeflatorGUI(obj);
            end
            
            obj.gui.makeGUI;
            
        end
        
    end
    
    methods % internal calculations
        
        function batch(obj) % read one file and write it to disk after deflating it
            
%             t1 = tic;
            obj.reader.batch;
            obj.copyFrom(obj.reader);
            obj.head = obj.reader.head; % is this necessary?
%             read_time = toc(t1); 
            
%             t2 = tic;
            obj.save;
%             save_time = toc(t2); 
            
%             fprintf('read time: %f | save time: %f\n', read_time, save_time); 
            
            if ~isempty(obj.gui)
                obj.gui.updateGUI;
            end
            
            obj.reader.clear;
            
        end
                
        function save(obj) % order the buffer object to save the file
                       
            if ~isempty(obj.gui) && ~isempty(obj.buffers.gui)
                obj.buffers.gui.update;
            end
            
%             obj.head = head.Header(obj.reader.head);
                        
            obj.buffers.directory = obj.out_subdir.pwd;
            
            [~, filename, ext] = fileparts(obj.reader.prev_filename);
            
%             obj.buffers.input('images', obj.images, 'images_raw', obj.images_raw, 'imagesd_cut', obj.images_cut, ...
%                 'images_sum', obj.images_sum, 'num_sum', obj.num_sum, 'cut_pos', obj.cut_pos,... 
%                 'timestamps', obj.timestamps, 't_end_stamp', obj.t_end_stamp, 't_end', obj.t_end,...
%                 'psfs', obj.psfs, 'psf_sampling', obj.psf_sampling, 'lightcurves', obj.lightcurves);
            
            if obj.buffers.use_deflate && strcmpi(ext, '.h5')
                ext = [ext 'z'];
            end

            obj.buffers.input(obj);

            obj.buffers.directory = obj.out_subdir.pwd; % do we need this??
            obj.buffers.filename = [filename ext];
            
            obj.buffers.save;
            obj.buffers.nextBuffer;
            
        end
                
        function [c, new_file] = checkFileExists(obj, file_src, directory) % check if file exists in the destination directory
            
            [~, name, ext] = fileparts(file_src);
            
            if strcmpi(ext, '.h5')
                new_file = [fullfile(directory, name) ext 'z'];
            else
                new_file = [fullfile(directory, name) ext];
            end
            
            if exist(new_file, 'file')
                c = 1;
                f1 = dir(file_src);
                f2 = dir(new_file);
                if f1.bytes>f2.bytes*50 % if source file (f1) is much bigger than new file (f2) that means f2 is not properly copied 
                    c = 0;
                end
            else
                c = 0;
            end
            
        end
                
    end
    
end