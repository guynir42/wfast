classdef Deflator < file.AstroData
   
    properties (Transient=true)
        
        gui@file.gui.DeflatorGUI;
        
        prog@util.sys.ProgressBar;
        
        futures = {}; % for using parfeval
        
        timer;
        
    end
    
    properties % objects and resources
        
        head@head.Header;
        
        src_dir@util.sys.WorkingDirectory;
        out_dir@util.sys.WorkingDirectory;
        out_dir_backup@util.sys.WorkingDirectory;
        
        src_subdir@util.sys.WorkingDirectory;
        out_subdir@util.sys.WorkingDirectory;
        out_subdir_backup@util.sys.WorkingDirectory;
        
        reader@file.Reader;
        buffers@file.BufferWheel;
        
    end
    
    properties(Dependent=true)
        
    end
    
    properties % switches and controls
        
        auto_deflate_destination = getenv('DATA_EXTRAS');
        auto_deflate_backup_dir = getenv('DATA_BACKUP'); 
        
        use_auto_delete = 0;
        
        use_parallel_auto_deflate = 1;
        
        use_auto_backup = 1;
        
        use_copy_text_files = 1;
        use_copy_catalogs = 1;
        use_copy_calibration = 1;
        
        debug_bit = 1;
        
    end
    
    properties(Hidden=true)
        
        brake_bit = 1;
        
        version = 1.02;
                
    end
    
    methods % constructor
        
        function obj = Deflator(other)
            
            if nargin==0 || isempty(other)
                if obj.debug_bit, fprintf('Deflator constructor v%4.2f\n', obj.version); end
                
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
                output = getenv('DATA');
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
        
        function set.head(obj, header)
            
            obj.head= header;

            obj.buffers.head= header;
                        
        end
        
    end
    
    methods % API
        
        function run(obj, varargin)
            
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
                    
                    disp(['going over directory: ' obj.src_subdir.pwd]);
                    
                    if obj.brake_bit
                        return;
                    end

                    obj.resetFilenames;
                    
                    if obj.use_copy_text_files
                        f = obj.src_subdir.match('*.txt');
                        for ii = 1:length(f)
                            
                            [~, name, ext] = fileparts(f{ii});
                            disp(['copying file ' name ext ' from ' obj.src_subdir.pwd ' to ' obj.out_subdir.pwd])
                            copyfile(sa(obj.src_subdir.pwd,name,ext), sa(obj.out_subdir.pwd, name, ext));
                            
                            if obj.use_auto_backup
                                disp(['copying file ' name ext ' from ' obj.src_subdir.pwd ' to ' obj.out_subdir_backup.pwd ' (backup)'])
                                copyfile(sa(obj.src_subdir.pwd,name,ext), sa(obj.out_subdir_backup.pwd, name, ext));
                            
                            end
                            
                        end
                    end

                    if obj.use_copy_catalogs
                        f = obj.src_subdir.match('catalog.mat');
                        if ~isempty(f)
                            [~, name, ext] = fileparts(f{1});
                            disp(['copying file ' name ext ' from ' obj.src_subdir.pwd ' to ' obj.out_subdir.pwd])
                            copyfile(sa(obj.src_subdir.pwd,name,ext), sa(obj.out_subdir.pwd, name, ext));
                        end
                    end
                    
                    if obj.use_copy_calibration
                        f = obj.src_subdir.match('calibration*.mat');
                        if ~isempty(f)
                            [~, name, ext] = fileparts(f{1});
                            disp(['copying file ' name ext ' from ' obj.src_subdir.pwd ' to ' obj.out_subdir.pwd])
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
            
            disp('done deflating files...');
            
        end
        
        function run_async(obj)
            
            obj.futures{end+1} = parfeval(@obj.run, 0); 
            
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
            
            t = datetime(util.sys.date_dir('now'));
            
            t = t + days(1); % deflate on the next morning!
            
            startat(obj.timer, [t.Year, t.Month, t.Day, hour, 0, 0]);
            
        end
        
        function callback_timer(obj, ~, ~)
            
            disp('Running auto-deflate!');
            
            obj.makeGUI;
            
            d = util.sys.WorkingDirectory(getenv('DATA_TEMP')); 
            list = d.dir('',1);
            d.cd(list{end});
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
        
        function batch(obj)
            
            t1 = tic;
            obj.reader.batch;
            obj.copyFrom(obj.reader);
            obj.head = obj.reader.head;
            read_time = toc(t1); 
            
            t2 = tic;
            obj.save;
            save_time = toc(t2); 
            
%             fprintf('read time: %f | save time: %f\n', read_time, save_time); 
            
            if ~isempty(obj.gui)
                obj.gui.updateGUI;
            end
            
            obj.reader.clear;
            
        end
                
        function save(obj)
                       
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
                
        function [c, new_file] = checkFileExists(obj, file_src, directory)
            
            [~, name, ext] = fileparts(file_src);
            
            new_file = [fullfile(directory, name) ext 'z'];
            if exist(new_file, 'file')
                c = 1;
                f1 = dir(file_src);
                f2 = dir(new_file);
                if f1.bytes>f2.bytes*4
                    c = 0;
                end
            else
                c = 0;
            end
            
        end
                
    end
    
end