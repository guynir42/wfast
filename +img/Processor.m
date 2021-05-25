classdef Processor < dynamicprops
    
    properties(Transient=true)
        
        gui;
        
        prog@util.sys.ProgressBar;
        
    end
    
    properties % objects
        
        head@head.Header;
        cat@head.Catalog;
        
        reader@file.Reader;
        cal@img.Calibration;
        
        phot@img.Photometry;
        light@img.Lightcurves;
        
        flux_buf@util.vec.CircularBuffer;
        
        pars@util.text.InputVars; % store all controls in this
        last_pars@util.text.InputVars; % parameters used in latest run
        
        data@struct; % file.AstroData; % keep track of all the regular data objects, and can dynamically add more (clear all every file)
        buffer@struct; % file.AstroData; % keep a few sets of data from previous files to calculate deeper stats
        data_logs@struct; % struct with data collected for all images across entire run 
        
        timings = struct; % collect how much time (seconds) each part of the analysis takes
        
        futures = {}; % used for running asynchronuous tasks
        
        func_handle; % generic function that takes as input the "data" object and modifies it, and also gets the "pars" object (can also give a cell of function handles)
        
        % ... add streak detection maybe?
        
    end
    
    properties % inputs/outputs
        
        file_index = 1;
        
        to_process_index = 1; 
                
        % references for quick_align 
        ref_image = [];
        ref_positions = [];
        
        forced_coordinates = {};
        forced_indices = []; 
        
        % table with find_stars results
        stars; 
        
    end
    
    properties % switches/controls
        
        % the rest are stored in the pars object        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        brake_bit = 1; % when running, this is 0, when it is turned to 1 (e.g., by GUI) the run is paused
        
        buf_idx = 1; % which buffer are we looking at right now
        prev_buf_idx = 0; % previously processed buffer
        
        phot_ap_idx = []; % index of photometric types ("aperture index")
        
        % these are interpreted from e.g., coordinates 
        fits_roi_position = []; % center of the ROI x,y position in pixels
        fits_roi_size = []; % height and width of the ROI in pixels
        
        forced_cutout_motion_per_file = []; 
        
        failed_file_counter = 0;
        
        backup_pars@util.text.InputVars;
        
        output_folder_name = ''; % folder containing the processed data and logs (e.g., process_2021-05-21)
        
        version = 1.02;
        
    end
    
    methods % constructor and make pars/data objects
        
        function obj = Processor(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'img.Processor')
                if obj.debug_bit>1, fprintf('Processor copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('Processor constructor v%4.2f\n', obj.version); end
                
                obj.head = head.Header;
                obj.cat = head.Catalog;
                obj.cat.head = obj.head;

                obj.reader = file.Reader;
                obj.reader.head = obj.head;
                
                obj.cal = img.Calibration;

                obj.phot = img.Photometry;
                
                obj.light = img.Lightcurves;
                obj.light.head = obj.head;
                obj.light.cat = obj.cat;

                obj.flux_buf = util.vec.CircularBuffer;

                obj.prog = util.sys.ProgressBar;
                
                obj.makeParsObject;
                obj.reset;
                
            end
            
        end
             
        function makeParsObject(obj)
            
            obj.pars = util.text.InputVars;
            
            obj.pars.input_var('reset', true); obj.pars.add_comment('reset the object data at the beginning of the run'); 
            
            obj.pars.input_var('num_files', [], 6); obj.pars.add_comment('maximum number of files to run. Run may end before that if files run out or stars are lost'); 
            obj.pars.input_var('cut_size', 15); obj.pars.add_comment('number of pixels in a cutout'); 
            obj.pars.input_var('num_stars', 2000, 6); obj.pars.add_comment('maximum number of stars to extract from the images');
            obj.pars.input_var('threshold', 5); obj.pars.add_comment('how many "sigma" above the noise'); 
            obj.pars.input_var('saturation', 5e4); obj.pars.add_comment('for a single image, if stack is read, this is multiplied by num_sum'); 
            obj.pars.input_var('use_bg_map', false, 6); obj.pars.add_comment('if true, will adjust each point in the stack by the interpolated b/g values, instead of the median value'); 
            
            obj.pars.input_var('use_coadds', true); obj.pars.add_comment('produce deeper coadds and use them for finding stars, CR mask, etc.'); 
            obj.pars.input_var('coadd_size', 10); obj.pars.add_comment('how many files should be loaded into one coadd'); 
            obj.pars.input_var('coadd_method', 'sum'); obj.pars.add_comment('how to coadd images, use "sum", "median", etc...'); 
            
            obj.pars.input_var('redo_photometry_distance', 2); obj.pars.add_comment('maximum pixel shift above which we redo photometry for that file'); 
            obj.pars.input_var('realign_attempts', 2); obj.pars.add_comment('how many times to try to realign the positions before failing that file'); 
            
            obj.pars.input_var('lock_adjust', 'stars'); 
            obj.pars.add_comment(sprintf('which cutouts should be moved together (locked), which can move separately (unlocked).\n Choose "all", "none", or "stars" that only locks targets with GAIA magnitudes')); 
            
            obj.pars.input_var('use_astrometry', true, 6); obj.pars.add_comment('make sure to solve the astrometry before continuing'); 
            obj.pars.input_var('use_require_astrometry', true, 6); obj.pars.add_comment('will throw an error if no astrometric solution is found'); 
            obj.pars.input_var('use_remove_bad_matches', false, 6); obj.pars.add_comment('remove stars that have no GAIA match if their S/N is too low'); 
            obj.pars.input_var('bad_match_min_snr', 30); obj.pars.add_comment('if an object has no GAIA match, it can be removed if its S/N is below this threshold'); 
            obj.pars.input_var('use_forced_photometry', true); obj.pars.add_comment('add a cutout at the given coordinates, even if there is no source found there'); 
            obj.pars.input_var('forced_extra_positions', {}); obj.pars.add_comment('to give additional positions (besides header RA/Dec)\n give coordinates as [RA,Dec] pairs (in deg) inside a cell array'); 
            obj.pars.input_var('use_interp_forced', true, 6, 'use_interpolate_forced', 'interpolate_forced', 'interp_forced'); obj.pars.add_comment('check the position of the forced cutout at time of end of run, if it moves, move the cutout linearly with time'); 
            
            obj.pars.input_var('use_auto_load_cal', true, 6); obj.pars.add_comment('load calibration based on the date of the files'); 
            obj.pars.input_var('use_auto_aperture', true, 6); obj.pars.add_comment('lset the photometric aperture based on the measured psd width X3'); 
            
            obj.pars.input_var('num_stars_filter_kernel', 200, 6); obj.pars.add_comment('use a rough star finding on the brightest stars, just to estimate the PSF width'); 
            obj.pars.input_var('initial_guess_psf_width', 1.0); obj.pars.add_comment('typical value just for finding those stars on the 1st attempt'); 
            obj.pars.input_var('fwhm_indices', 100:200); obj.pars.add_comment('which stars to use for calculating FWHM (other stars will have NaN results)'); 
            
            obj.pars.input_var('use_check_flux', true, 6); obj.pars.add_comment('check if the star flux is gone for a few files and cut the run short'); 
            obj.pars.input_var('max_failed_files', 3, 6); obj.pars.add_comment('if star flux is lost for more than this number of files, quit the run'); 
            
            obj.pars.input_var('display_num_rect_stars', 30, 10); obj.pars.add_comment('how many rectangles to display'); 
            obj.pars.input_var('use_display_rect_text', false, 6); obj.pars.add_comment('show the clip number for each square'); 
            
            obj.pars.input_var('use_save_results', false, 6); obj.pars.add_comment('save the lightcurves etc. to file at end of run'); 
            obj.pars.input_var('use_overwrite', false, 6); obj.pars.add_comment('if true, will happily clear existing "processing" folders with the same date'); 
            obj.pars.input_var('output_folder', ''); obj.pars.add_comment('folder inside original data, typically "processor_YYYY-MM-DD" for storing the results and log files'); 
            
            obj.pars.input_var('fits', util.text.InputVars); obj.pars.add_comment('nested parameters related to making FITS files'); 
            obj.pars.fits.input_var('use_save', false, 10); obj.pars.fits.add_comment('save FITS files in a folder below the raw images'); 
            obj.pars.fits.input_var('use_roi',false, 6); obj.pars.fits.add_comment('use a region of interest (ROI) to cut before saving the FITS'); 
            obj.pars.fits.input_var('roi_size', 512, 6); obj.pars.fits.add_comment('define the size of the ROI as scalar (square region) or [height,width] in pixels');
            obj.pars.fits.input_var('roi_position', [], 6); obj.pars.fits.add_comment('define the center point of the ROI as [x,y] in pixels');             
            obj.pars.fits.input_var('roi_coordinates', [], 6, 'roi_coords'); obj.pars.fits.add_comment('define the center of the ROI as [RA, Dec] in degrees. Or "header" to copy coordinates from the header.'); 
            obj.pars.fits.input_var('use_flip', false, 6); obj.pars.fits.add_comment('flip the image by 180 degrees before saving to FITS (this is done after cutting out the ROI'); 
            obj.pars.fits.input_var('use_coadds', false, 6); obj.pars.fits.add_comment('save FITS images of coadded images'); 
            obj.pars.fits.input_var('directory', 'FITS'); obj.pars.fits.add_comment('what to call the folder where FITS files are saved (relative to original file directory or absolute path'); 
            obj.pars.fits.input_var('rename', '', 'name', 'filename'); obj.pars.fits.add_comment('rename each FITS file to this string, followed by a zero padded serial number'); 
            obj.pars.fits.input_var('use_finding', true); obj.pars.fits.add_comment('save a finding chart PNG of the first image'); 
            
            obj.pars.input_var('header', {}, 'header corrections'); obj.pars.add_comment('cell array with keyword-value pairs to modify the header info from file'); 
            
            obj.pars.input_var('func_pars', [], 'function parameters', 'func parameters', 'function pars'); obj.pars.add_comment('additional parameters to "func handle". Could be a struct.'); 
            
            
            
            obj.backup_pars = obj.pars; % a copy of the handle, pointing to the same object, used to restore the parameters at the end of each run
            
        end
        
        function data_out = makeDataStruct(obj)
            
            data_out = struct; % make a new struct
            
            % these fields basically follow what is in file.AstroData, but it is easier to work with a struct
            
            data_out.images = []; % this is raw images and it is usually what we save on file
        
            data_out.timestamps = []; % output timestamps (if available)
            data_out.t_start = []; % absolute date and time (UTC) when first image is taken
            data_out.t_end = []; % absolute date and time (UTC) when batch is finished
            data_out.t_end_stamp = []; % timestamp when batch is finished (hopefully, this is the same time that t_end is recorded). 
         
            data_out.juldates = []; % translation of timestamps to julian date using t_end_stamp
        
            data_out.stack  = [];% sum of the full frame image
            data_out.num_sum = []; % if the images are summed, how many frames were added. If equal to 1, the sum is the same as the images. 
        
            data_out.cutouts = []; % this is raw cutouts and it is usually what we save on file
            data_out.positions = []; % only for cutouts. a 2xN matrix (X then Y, N is the number of cutouts). 
            data_out.object_idx = []; % what is the index of the object closest to the coordinates given
            
            data_out.image_proc = []; % image after calibration
            data_out.image_cut = []; % image after calibration and removing cutouts 
            data_out.cutouts_sub = []; % cutouts after background subtraction
            data_out.positions_new = []; % positions after adjustments based on star positions
            data_out.fwhm = []; % FWHM for each cutout, in arcsec
            data_out.fwhm_interp = []; % FWHM intepolated using a 2D polynomial fit to the non-NaN FWHM values
            data_out.psf_width = []; % best estimate for the average PSF width (gaussian sigma)
            
            data_out.exposure_time = []; 
            data_out.file_index = []; 
            
            data_out.coadd_image = []; 
            data_out.coadd_number = [];
            data_out.coadd_exposure = []; 
            
        end
        
        function data_out = makeDataObject(obj) % to be deprecated! 
            
            data_out = file.AstroData; 
            data_out.addProp('image_proc'); % image after calibration
            data_out.addProp('image_cut'); % image after calibration and removing cutouts 
            data_out.addProp('cutouts_sub'); % cutouts after background subtraction
            data_out.addProp('positions_new'); % positions after adjustments based on star positions
            data_out.addProp('fwhm'); % FWHM for each cutout, in arcsec
            data_out.addProp('fwhm_interp'); % FWHM intepolated using a 2D polynomial fit to the non-NaN FWHM values
            data_out.addProp('psf_width'); % best estimate for the average PSF width (gaussian sigma)
            
            data_out.addProp('exposure_time'); 
            data_out.addProp('file_index'); 
            
            data_out.addProp('coadd_image'); 
            data_out.addProp('coadd_number'); 
            data_out.addProp('coadd_exposure'); 
            
        end
           
    end
    
    methods % reset/clear
        
        function restoreBackupPars(obj)
            
            obj.pars = obj.backup_pars; 
            
        end
        
        function resetToDefaultPars(obj)
            
            obj.pars = obj.makeParsObject;
            
        end
        
        function reset(obj) % call this at the start of a new run
            
            % objects
            obj.reader.reset;
            obj.reader.loadFiles; 
            
            obj.cat.reset;
            obj.phot.reset;
            obj.light.reset;
            obj.flux_buf.reset;
            
            % the index of the type of photometry
            obj.phot_ap_idx = [];
            
            % buffer
            obj.buffer = obj.makeDataStruct; % making a single (empty content) struct means we can later assign full data structures into this
            
            % info about the run
            obj.data_logs = struct; 
            obj.data_logs.run_date = ''; % date when images were taken (date of beginning of night) in YYYY-MM-DD format
            obj.data_logs.camera = ''; % name of camera (can be e.g., Zyla or Balor)
            obj.data_logs.project = ''; % name of telescope/project (can be e.g., Kraar, LAST or WFAST)
            obj.data_logs.fwhm = []; % for each file get the median FWHM (in arcsec)
            obj.data_logs.background = []; % for each file get the median background value from all stars
            obj.data_logs.variance = []; % for each file get the median variance value from all stars
            
            % counters
            obj.file_index = 1;
            obj.to_process_index = 1;
            obj.failed_file_counter = 0;
            
            % parameters
            obj.last_pars = util.text.InputVars.empty; 
            
            % timing data
            obj.timings = struct; 
            
            % clear data logs
            obj.data_logs = struct; 
            
            % other stuff
            obj.forced_coordinates = {}; 
            obj.forced_indices = []; 
            
            obj.clear;
            
        end
        
        function clear(obj) % call this each batch
            
            obj.phot.clear;
            obj.cat.clear;
            obj.reader.clear;
            
            obj.data = struct.empty; % remove the handle
            
%             obj.data.timestamps = []; % the time of the beginning of the acquisition of the image/stack 
%             obj.data.juldates = []; % matched to global time
%             obj.data.image = []; % raw images or stack from file
%             obj.data.image_proc = []; % after calibration, background subtraction, etc...
%             obj.data.image_cut = []; % image after removing all the stars using mexCutout
%             obj.data.num_sum = 1; % if we read something else, it will be altered
%             obj.data.cutouts = []; % cutouts are only made from the stack, not from the full frame rate data
            
        end
        
    end
    
    methods % getters
        
        function val = getNumFiles(obj) % get the number of files we can read, given file list and the limit given in "pars"
            
            if isempty(obj.pars.num_files)
                val = length(obj.reader.filenames);
            else
                val = min(obj.pars.num_files, length(obj.reader.filenames)); 
            end
            
        end
        
        function val = this_filename(obj, full_path)
            
            if nargin<2 || isempty(full_path)
                full_path = 1;
            end
            
            val = obj.reader.filenames{obj.file_index}; 
            
            if full_path==0
                [~,val, ext] = fileparts(val);
                val = [val, ext]; 
            end
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % commands
        
        function obj = run(obj, varargin)
            
            obj.pars = util.oop.full_copy(obj.pars); % generate a new (temporary) copy of the parameters
            obj.pars.scan_vars(varargin{:}); % scan any user-defined parameter changes
            
            on_cleanup = onCleanup(@obj.finishup); % on finishup() the original "pars" is restored
            
            try
            
                obj.startup; 

                for ii = obj.file_index:obj.getNumFiles % main loop

                    if obj.brake_bit, break; end

                    obj.load_data;
                    
                    drawnow; 
                    
                    if mod(obj.file_index - 1, obj.pars.coadd_size)==0 % finished collecting N files
                        obj.calculateCoadds; 
                        obj.process; % go over unprocessed files
                    end
                    
                end
            
            catch ME
                
                if obj.pars.use_save_results 
                    obj.write_log(ME.getReport); % make sure to also save the error in the log file 
                end
                
                if obj.file_index >= obj.getNumFiles % run has finished all useful files
                    warning(ME.getReport); % this is some after-processing error, we can live with it
                else % run has caught an error mid-way
                    rethrow(ME); 
                end
            end
            
        end
        
        function cont(obj)
            
            if isempty(obj.last_pars)
                obj.run('reset', 0); % use the current parameters 
            else
                obj.run(obj.last_pars.output_vars{:}, 'reset', 0); % continue with the same parameters as we used in the last run
            end
            
        end
        
        function run_async(obj)
            
        end
    
        function stop(obj)
            
            obj.brake_bit = 1;
            
        end
        
    end
    
    methods % utilities
        
        function browse(obj)
            
            obj.reader.browseDir; 
            
            obj.reader.loadFiles; 
            
        end
        
        function saveResults(obj)
            
            % save the lightcurves, catalog, summary text file, what else?
            
        end
        
        function str = printout(obj)
            
            str = sprintf('file: %d/%d', obj.file_index, obj.getNumFiles); 
            
            str = sprintf('%s | folder: %s', str, strrep(obj.reader.dir.two_tail, '\', ' \ ')); 
            
            str = sprintf('%s | coords: %s%s', str, obj.head.RA, obj.head.Dec); 
            
            if ~isempty(obj.head.AIRMASS)
                str = sprintf('%s | a.m.= %4.2f', str, obj.head.AIRMASS); 
            end
            
        end
        
    end
    
    methods(Hidden=false) % internal processing
        
        function startup(obj, varargin)
            
            if isempty(obj.prog)
                obj.prog = util.sys.ProgressBar;
            end
            
            if obj.pars.reset % all the things that must happen when a run begins (and not happen when a run continues)
                
                obj.displayInfo(sprintf('Starting new run with %d files', obj.getNumFiles));
                
                obj.reset;
                obj.prog.reset(obj.getNumFiles);
            
                obj.data_logs.run_date = obj.find_folder_date;
                obj.data_logs.camera = obj.find_camera;
                obj.data_logs.project = obj.find_project;

                if obj.pars.use_auto_load_cal
                    obj.displayInfo(sprintf('Loading calibration file for date %s', obj.data_logs.run_date)); 
                    obj.cal.loadByDate(obj.data_logs.run_date, obj.data_logs.camera, obj.data_logs.project); 
                end

                if ~obj.cal.checkDark
                    error('Cannot start a new run without loading darks into calibration object!');
                end

                if obj.pars.use_save_results

                    d = obj.pars.output_folder; 

                    if isempty(d)
                        t = datetime('now', 'TimeZone', 'UTC'); 
                        d = fprintf('processor_%4d-%2d-%2d', t.Year, t.Month, t.Day); 
                    end

                    obj.output_folder_name = fullfile(obj.reader.dir.cwd, d); 

                    if exist(obj.output_folder_name, 'dir') && obj.pars.use_overwrite==0
                        [~, d] = fileparts(obj.output_folder_name); 
                        error('Cannot save results, folder %s already exists!', d); 
                    end

                end

            end
            
            obj.prog.unpause;
            
            obj.brake_bit = 0;
            
        end
        
        function finishup(obj)
            
            obj.last_pars = obj.pars; % store the parameters used in last run (used to continue run with same parameters)
            
            obj.pars = obj.backup_pars; % restore the parameters of the original object before parsing inputs
            
            obj.prog.finish;
            
            obj.light.finishup;
            
            if obj.pars.use_save_results
                obj.saveResults;
            end
            
            if obj.debug_bit, fprintf('Finished run with %d files.\n', obj.file_index); end
            
        end
        
        function load_data(obj)
            
            if obj.debug_bit>1
                util.text.date_printf('Loading file number %d.', obj.file_index); 
            end
            
            try % get the data
                
                obj.clear; % get rid of existing datasets
                
                t0 = tic;
                
                obj.reader.batch; % load the images from file

                idx = mod(obj.file_index-1,obj.pars.coadd_size) + 1; % which buffer to fill                
                idx_prev = mod(obj.file_index-2,obj.pars.coadd_size) + 1; % which buffer to fill
                obj.data = obj.makeDataStruct; % dump the old data and construct a new object
                
                obj.data.images = obj.reader.images;
                obj.data.stack = obj.reader.stack;
                obj.data.num_sum = obj.reader.num_sum;
                obj.data.timestamps = obj.reader.timestamps;
                obj.data.juldates = obj.reader.juldates;
                obj.data.t_start = obj.reader.t_start;
                obj.data.t_end = obj.reader.t_end;
                obj.data.t_end_stamp = obj.reader.t_end_stamp;
                obj.data.exposure_time = obj.head.EXPTIME; 
                
                obj.data.file_index = obj.file_index; 
                
                if isempty(obj.data.juldates)
                    % recaclculate it from timestamps
                end
                
                if size(obj.data.images,3)>1
                    error('I don''t know what to do with a data cube!');
                end
                
                obj.addTimingData('read', toc(t0)); 
                
            catch ME
                warning(ME.getReport);
                obj.failed_file_counter = obj.failed_file_counter + 1;
                obj.file_index = obj.file_index + 1;
                obj.reader.advanceFile;
                return; % skip this batch, report it, and continue! 
            end
            
            % apply head corrections, if any
            for ii = 1:2:length(obj.pars.header)

                if length(obj.pars.header)>ii
                    if isprop(obj.head, obj.pars.header{ii})
                        obj.head.(obj.pars.header{ii}) = obj.pars.header{ii+1}; 
                    else
                        error('Header does not have a property named "%s".', obj.pars.header{ii}); 
                    end
                else
                    error('Broken key-value pair for header correction'); 
                end

            end
            
            % keep track of some header data for the whole run
            list = {'airmass', 'alt', 'jd', 'sensor_temp', 'light', 'pressure', ...
                'humid_out', 'wind_speed', 'wind_dir', 'temp_out', 'focus_pos'};
            
            for jj = 1:length(list)
                s = substruct('.', list{jj}, '()', {obj.file_index,1}); 
                obj.data_logs = subsasgn(obj.data_logs, s, obj.head.(list{jj})); 
            end
            
            t0 = tic; %%%%%% run calibration on full-frame images %%%%%%%%

            if ~isempty(obj.data.images)
                obj.data.image_proc = obj.cal.input(obj.data.images); 
                obj.data.num_sum = 1; 
            elseif ~isempty(obj.data.stack)
                obj.data.image_proc = obj.cal.input(obj.data.stack, 'sum', obj.data.num_sum); 
            else
                error('Could not find images or stack to process...'); 
            end

            obj.addTimingData('calibration', toc(t0));
    
            if idx_prev>1 && idx_prev<=length(obj.buffer)
                obj.data.positions = obj.buffer(idx_prev).positions_new;
            end
            
            obj.buffer(idx) = obj.data; % note that the stuct remaining in "data" is not to be used unless re-loaded from the buffer
                
            obj.file_index = obj.file_index + 1;
            
        end
        
        function setBufferIndices(obj)
            
            obj.prev_buf_idx = find([obj.buffer.file_index]==obj.to_process_index - 1); % last buffer index
            obj.buf_idx = find([obj.buffer.file_index]==obj.to_process_index); % current buffer index we are processing

        end
        
        function process(obj)
            
            if length(obj.buffer)<obj.pars.coadd_size
                return; % skip processing until buffer is full
            end
            
            for ii = 1:obj.pars.coadd_size
                
                if obj.to_process_index>obj.file_index
                    break; % finished processing backlog in buffer
                end
                
                if obj.debug_bit>1
                    util.text.date_printf('Processing file number %d.', obj.to_process_index); 
                end
                
                obj.setBufferIndices;
                
                obj.data = obj.buffer(obj.buf_idx); % create a separate copy of the data loading it into "data" struct
                
                if ~isempty(obj.prev_buf_idx) && obj.prev_buf_idx>0 && obj.prev_buf_idx<=length(obj.buffer) % if a previous buffer exists...
                    obj.data.positions = obj.buffer(obj.prev_buf_idx).positions_new; % make sure to propagate the positions to the new buffer
                end
                
                %%%%%%%%%%%%%%%%%%%%% PROCESSING %%%%%%%%%%%%%%%%%%%%%

                if obj.pars.fits.use_save
                    obj.saveFITS; 
                end

                obj.calculateCutouts; 
                
                obj.calculatePhotometry; 
                
                obj.calculatePositions; 
                
                obj.calculateFWHM; 

                obj.calculateLightcurves; 
                
                if obj.pars.use_save_results
                    obj.write_log; % by default it should save the file_summary() result
                end

                if obj.debug_bit>1, disp(obj.file_summary); end

                if ~isempty(obj.gui)
                    obj.gui.update;
                end
            
                obj.prog.showif(obj.to_process_index);
                
                obj.buffer(obj.buf_idx) = obj.data; % store the results of all calculations back in the buffer
                
                obj.to_process_index = obj.to_process_index + 1; 
                
                drawnow;
                
            end % loop over buffer data
            
        end
        
        function calculateCoadds(obj)
            
            % make a deep coadd
            if obj.debug_bit>1
                util.text.date_printf('Making a deep coadd...'); 
            end
            
            t0 = tic; % make the deep coadd
            
            obj.setBufferIndices;
            
            obj.data = obj.buffer(obj.buf_idx); 
            
            I = zeros([size(obj.data.image_proc), length(obj.buffer)], 'like', obj.data.image_proc); 
            T = 0; % number of seconds total integration
            N = 0; % number of images coadded
            
            for ii = 1:length(obj.buffer)
                I(:,:,ii) = obj.buffer(ii).image_proc; 
                T = T + obj.buffer(ii).exposure_time; 
                N = N + 1; 
            end
            
            if util.text.cs(obj.pars.coadd_method, 'sum')
                obj.data.coadd_image = nansum(I,3); 
            elseif util.text.cs(obj.pars.coadd_method, 'median')
                obj.data.coadd_image = nanmedian(I,3); 
            else
                error('Unknown coadd method "%s". Use "sum" or "median", etc.', obj.pars.coadd_method); 
            end
            
            % calculate cosmic ray mask here... 
            
            obj.data.coadd_number = N;
            obj.data.coadd_exposure = T; 
            
            obj.addTimingData('coaddition', toc(t0));
            
            if isempty(obj.data.positions)
                
                obj.findStars;

                if obj.pars.use_astrometry

                    obj.solveAstrometry;

                    if obj.pars.use_forced_photometry

                        obj.forced_coordinates = horzcat({[obj.head.RA_DEG, obj.head.DEC_DEG]}, util.vec.torow(obj.pars.forced_extra_positions));

                        for ii = 1:length(obj.forced_coordinates)

                            ra_dec = obj.forced_coordinates{ii};
                            xy = obj.cat.coo2xy(ra_dec(1), ra_dec(2)); 

                            obj.data.positions(end+1,:) = xy; % add another row to the coordinate list
                            obj.forced_indices(end+1) = size(obj.data.positions,1); % keep track of indices of forced cutouts

                        end

                        obj.forced_cutout_motion_per_file = []; 

                        if obj.pars.use_interp_forced

                            h = util.oop.full_copy(obj.head); 
                            h.ephem.timeTravelHours(obj.getNumFiles.*obj.head.EXPTIME/3600); % move this Ephemeris object forward to the end of the run
                            h.ephem.debug_bit = 0; 
                            h.ephem.resolve; 

                            pos1 = obj.data.positions(obj.forced_indices(1),:); 
                            pos2 = obj.cat.coo2xy(h.RA_DEG, h.DEC_DEG);

                            if sqrt(sum((pos1-pos2).^2))<min(size(obj.data.image_proc))
                                obj.forced_cutout_motion_per_file = (pos2-pos1)/obj.getNumFiles; % should be a 2-vector
                            end

                        end

                    end

                end

                % store the first image to keep a record of the positions of all stars in case we need to re-align
                obj.ref_image = obj.data.image_proc;
                obj.ref_positions = obj.data.positions;

            end

            obj.buffer(obj.buf_idx) = obj.data; 
                            
        end
        
        function saveFITS(obj)
            
            t0 = tic; 
            
            %%%%% preprocessing of the image %%%%%%
            
            if obj.pars.fits.use_coadds
                
                if isempty(obj.data.coadd_image)
                    return; % not all buffers would have a coadd
                end
                
                I = double(obj.data.coadd_image); 
                N = obj.coadd_number;
                
            else
                I = double(obj.data.image_proc); 
                N = 1; 
            end
            
            if obj.pars.fits.use_roi

                if ~isempty(obj.pars.fits.roi_position)
                    obj.fits_roi_position = obj.pars.fits.roi_position; % if position is given explicitely
                elseif ~isempty(obj.pars.fits.roi_coordinates)

                    if obj.cat.success==0
                        error('Can not translate coordinates to pixel position without solving astrometry...');
                    end

                    coords = obj.pars.fits.roi_coordinates; % shorthand

                    if ischar(coords) % either a string with the coordinates (separated by +/-) or a command like "header"

                        if util.text.cs(coords, 'header', 'auto')
                            obj.fits_roi_position = obj.cat.coo2xy; % the default is to use the header info 
                        elseif util.text.cs(coords, 'center')
                            obj.fits_roi_position = floor(size(I)/2)+1; % image center
                        elseif (ismember('+', coords) || ismember('+', coords))
                            idx = regexp(coords, '[-+]');
                            if ~isempty(idx)
                                obj.fits_roi_position = obj.cat.coo2xy(coords(1:idx-1), coords(idx+1:end)); % use the given coordinates
                            else
                                error('Could not find a "+" or "-" in coordinate string.'); 
                            end
                        else
                            error('Must give coordinates as "RA+Dec" or "RA-Dec" or use "header" or "center"'); 
                        end

                    elseif isnumeric(coords) % assume we got [RA,Dec] in degrees
                        obj.fits_roi_position = obj.cat.coo2xy(coords); 
                    elseif iscell(coords) % a cell with 2 elements for RA and Dec
                        if length(coords)==2
                            obj.fits_roi_position = obj.cat.coo2xy(coords{:}); % use the given coordinates
                        else
                            error('Must input coordinates as as cell with two sexagesimal strings'); 
                        end
                    else
                        error('Coordinates must be given as a string or 2-element cell, or use "header" or "center" instead. '); 
                    end

                else % didn't get position or coordinates for the ROI
                    obj.fits_roi_position = floor(size(I)/2)+1; % just center the ROI by default
                end
                
            end
            
            obj.fits_roi_position = round(obj.fits_roi_position); 
            
            obj.fits_roi_size = obj.pars.fits.roi_size;
            if isscalar(obj.fits_roi_size)
                obj.fits_roi_size = [1 1].*obj.fits_roi_size;
            end

            x1 = obj.fits_roi_position(1) - floor(obj.fits_roi_size(2)/2);
            x2 = obj.fits_roi_position(1) + floor((obj.fits_roi_size(2)-1)/2);
            y1 = obj.fits_roi_position(2) - floor(obj.fits_roi_size(1)/2);
            y2 = obj.fits_roi_position(2) + floor((obj.fits_roi_size(1)-1)/2);
            
            I = I(y1:y2, x1:x2); 
            
            if obj.pars.fits.use_flip
                I = rot90(I,2); 
            end

            %%%%% figure out the path %%%%%%
            
            [base_dir, filename] = fileparts(obj.this_filename);
            
            d = obj.pars.fits.directory; 
            
            if ispc % on Windows
                if isempty(regexp(d, '^[A-Z]:\', 'once'))
                    d = fullfile(base_dir, d); % assume a relative path, so add the original filename's path 
                end
            else % on linux 
                if isempty(regexp(d, '^/', 'once'))
                    d = fullfile(base_dir, d); % assume a relative path, so add the original filename's path 
                end
            end
            
            d = strrep(d, ' (Weizmann Institute)', ''); % remove the dropbox's annoying spaces and parenthesis that FITS cannot handle (must have a link set up "Dropbox" --> "Dropbox (Weizmann Institute)"
            
            if ~isempty(obj.pars.fits.rename)
                filename = sprintf('%s_%04d', obj.pars.fits.rename, obj.to_process_index); 
            end
            
            filename = [filename, '.fits'];
            
            fullname = fullfile(d, filename); 
            
            %%%%% now handle the actual saving %%%%%%%
            
            if ~exist(d, 'dir')
                mkdir(d);
            end
            
            if obj.debug_bit>1
                util.text.date_printf('Saving FITS file to %s\n', fullname); 
            end
            
            fitswrite(I, fullname); 
            obj.head.writeFITS(fullname, obj.data.timestamps(1), obj.data.num_sum*N);
            
            %%%%%%%% additional housekeeping %%%%%%%%%%%
            
            try % save a finding chart
                                
                if obj.pars.fits.use_finding && obj.to_process_index==1
                    
                    f = figure('Name', 'finding chart', 'Position', [100 100 1000 1000]);
                    
                    ax = axes('Parent', f); 
                    
                    util.plot.show(I, 'auto', 1, 'monochrome', 1, 'fancy', 0); 
                                       
                    hold(ax, 'on'); 
                    
                    plot(ax, size(I,2)/2, size(I,1)/2, 'og', 'MarkerSize', 25, 'LineWidth', 2); 
                    
                    side = obj.head.ephem.side;
                    if obj.pars.fits.use_flip
                        if strcmpi(side, 'east')
                            side = 'west'; 
                        elseif strcmpi(side, 'west')
                            side = 'east'; 
                        end
                    end
                    
                    util.plot.compass(side, 'margin', [0.2 0.2], 'figure', f); 
                    
                    x = obj.fits_roi_size(2).*0.15;
                    y = obj.fits_roi_size(1).*0.95; 
                    w = 100; % width of the bar, in arcsec
                    
                    errorbar(ax, [1 1].*x, [1 1].*y, [1 1].*w/obj.head.SCALE, [0 0], '.m', 'Horizontal', 'LineWidth', 3);
                    text(ax, x, y, sprintf('  %d"', w), 'Color', 'm', 'FontSize', 20); 
                    
                    saveas(f, fullfile(d, 'finding_chart.png'));
                    
                    delete(f); 
                    
                end
                
            catch ME
                warning(ME.getReport); 
            end
            
            try % add entry in the text file

                if obj.file_index==1
                    fid = fopen(fullfile(d, '00README.txt'), 'wt'); % overwrite existing log file
                    on_cleanup = onCleanup(@() fclose(fid));
                    
                    fprintf(fid, '%s/%s observations taken on %s. Reduced on %s. Aperture: %dcm, filter: %d-%dnm, EXPTIME: %4.3fs, pix-scale: %4.2f"\n', ...
                        obj.head.PROJECT, obj.head.INST, obj.head.STARTTIME, util.text.time2str('now'), round(obj.head.TEL_APER), ...
                        round(obj.head.filter_range), obj.head.EXPTIME, obj.head.SCALE); 
                    
                    fprintf(fid, 'finding_chart.png\n'); 
                    
                else
                    fid = fopen(fullfile(d, '00README.txt'), 'at'); % append to log file 
                    on_cleanup = onCleanup(@() fclose(fid));
                end
                
                fprintf(fid, '%s\n', filename); 
                
            catch ME
                warning(ME.getReport); 
            end
            
            obj.addTimingData('fits', toc(t0));
            
        end
        
        function val = getAverageWidth(obj) % flux weighted average PSF width
        
            if isempty(obj.phot.widths)
                val = [];
            else
                val = util.vec.weighted_average(obj.phot.widths(:,:,obj.phot_ap_idx), obj.phot.fluxes(:,:,obj.phot_ap_idx), 2);
            end
            
        end
        
        function val = getAverageFWHM(obj) % flux weighted average PSF width
        
            if isempty(obj.data.fwhm)
                val = [];
            else
                val = util.vec.weighted_average(obj.data.fwhm, obj.phot.fluxes(:,:,obj.phot_ap_idx), 2);
            end
            
        end
        
        function val = getAverageOffsets(obj) % flux weighted average offsets [dx,dy]
        
            if isempty(obj.phot.offsets_x) || isempty(obj.phot.offsets_y)
                val = [];
            else
                dx = util.vec.weighted_average(obj.phot.offsets_x(:,:,obj.phot_ap_idx), obj.phot.fluxes(:,:,obj.phot_ap_idx), 2);
                dy = util.vec.weighted_average(obj.phot.offsets_y(:,:,obj.phot_ap_idx), obj.phot.fluxes(:,:,obj.phot_ap_idx), 2);
                val = [dx, dy];
            end
            
        end
        
        function findStars(obj)
            
            t0 = tic;
            
            obj.displayInfo('Finding stars.'); 
            
            if obj.pars.use_coadds
                I = obj.data.coadd_image;
                N = obj.data.coadd_number.*obj.data.num_sum; 
            else
                I = obj.data.image_proc; 
                N = obj.data.num_sum; 
            end
            
            T = util.img.quick_find_stars(I, 'psf', obj.pars.initial_guess_psf_width, 'number', obj.pars.num_stars_filter_kernel,...
               'dilate', obj.pars.cut_size-5, 'saturation', obj.pars.saturation.*N, 'unflagged', 1); 
            
            % estimate the real PSF width from found stars
            C = util.img.mexCutout(I, T.pos, obj.pars.cut_size, NaN, NaN); 
            
            s = util.img.photometry2(C, 'aperture', obj.pars.initial_guess_psf_width.*3, 'use_forced', 1); 

            W = util.img.fwhm(C-permute(s.forced_photometry.background, [3,4,1,2]), ...
                'method', 'filters', 'step_size', 0.01, 'min_size', 0.5, 'max_size', 4); % arbitrary range for getting an estimate of the PSF width
            
%             obj.data.psf_width = util.vec.weighted_average(s.forced_photometry.width, s.forced_photometry.flux, 2);
            obj.data.psf_width = util.vec.weighted_average(W/2.355, s.forced_photometry.flux, 2);
            
            if obj.pars.use_auto_aperture
                obj.phot.aperture = floor(obj.data.psf_width.*3);
            end
            
            T = util.img.quick_find_stars(I, 'psf', obj.data.psf_width, 'number', obj.pars.num_stars,...
               'dilate', obj.pars.cut_size-5, 'saturation', obj.pars.saturation.*N, 'sigma', obj.pars.threshold, 'unflagged', 0); 
            
            if isempty(T)
                error('Could not find any stars using quick_find_stars!');
            end

            obj.data.positions = T.pos;
            
            obj.stars = T; % store a copy of this table

            obj.addTimingData('find stars', toc(t0)); 
            
        end
        
        function solveAstrometry(obj)
            
            t0 = tic;
            
            obj.displayInfo('Running astrometry'); 
            
            if util.text.cs(obj.head.cam_name, 'balor')
                obj.cat.input_rotation = -60; 
            elseif util.text.cs(obj.head.cam_name, 'zyla')
                obj.cat.input_rotation = 15; 
            end

            obj.cat.input(obj.stars); 

            if obj.cat.success==0
                error('Could not find astrometric solution!'); 
            end

            if obj.pars.use_remove_bad_matches

                bad_idx = isnan(obj.cat.magnitudes) & obj.stars.snr<obj.pars.bad_match_min_snr;

                obj.cat.data = obj.cat.data(~bad_idx, :); 
                obj.cat.magnitudes = obj.cat.magnitudes(~bad_idx); 
                obj.cat.coordinates = obj.cat.coordinates(~bad_idx, :); 
                obj.cat.temperatures = obj.cat.temperatures(~bad_idx); 
                obj.cat.positions = obj.cat.positions(~bad_idx, :); 
                obj.data.positions = obj.data.positions(~bad_idx, :); 
                obj.stars = obj.stars(~bad_idx, :); 

            end
            
            if obj.cat.success
                str = sprintf('Found astrometric solution: %s%s', head.Ephemeris.deg2hour(obj.cat.central_RA), head.Ephemeris.deg2sex(obj.cat.central_Dec)); 
            else
                str = 'Cannot find an astrometric solution!'; 
            end

            obj.displayInfo(str);
            
            obj.addTimingData('astrometry', toc(t0));
            
        end
        
        function calculateCutouts(obj)

            t0 = tic; 
            
            [obj.data.cutouts, obj.data.image_cut] = util.img.mexCutout(obj.data.image_proc, ...
                obj.data.positions, obj.pars.cut_size, NaN, NaN); % replace and fill up using NaNs (it is safe, the processed image is single precision

            obj.addTimingData('cutouts', toc(t0)); 

        end
        
        function calculatePhotometry(obj)
            
            t0 = tic; 

            obj.phot.input(obj.data.cutouts, 'positions', obj.data.positions, 'timestamps', obj.data.timestamps, 'juldates', obj.data.juldates); 

            if isempty(obj.phot_ap_idx)
                obj.phot_ap_idx = obj.find_phot_ap;
            end

            obj.data.cutouts_sub = obj.data.cutouts - permute(obj.phot.backgrounds(:,:,obj.phot_ap_idx), [3,4,1,2]); % subtract the background from the cutouts

            obj.data_logs.background(obj.to_process_index,1) = util.stat.median2(obj.phot.backgrounds(:,:,obj.phot_ap_idx)); 

            obj.flux_buf.input(obj.phot.fluxes(:,:,obj.phot_ap_idx)); % update the flux buffer (used to check the fluxes didn't suddenly disappear)

            if ~isempty(obj.phot.gui) && obj.phot.gui.check, obj.phot.gui.update; end
            
            obj.addTimingData('photometry', toc(t0)); 

        end
        
        function calculatePositions(obj)
            
            t0 = tic; 

            for ii = obj.pars.realign_attempts

                if obj.checkFluxes
                    obj.adjustPositions; % use the centroids to push the cutouts a little
                else
                    if obj.debug_bit, disp('Lost star positions, using quick_align...'); end
                    obj.realignPositions; 
                end

                if any(median(abs(obj.data.positions - obj.data.positions_new))>=obj.pars.redo_photometry_distance) % new positions are too far (in x or y) from old positions

                    obj.flux_buf.back; % roll back the index to put in new flux measurements
                    
                    obj.data.positions = obj.data.positions_new; % try again with updated positions
                    
                    obj.calculateCutouts;
                    obj.calculatePhotometry;
                    
                else
                    break; % if the positions did not change dramatically, no need for more iterations
                end

            end

            obj.addTimingData('adjust_positions', toc(t0)); 
            
        end
        
        function adjustPositions(obj) % minor position adjustment based on photometric centroids
            
            import util.text.cs;

            mean_offsets = obj.getAverageOffsets;
            if any(isnan(mean_offsets))
                mean_offsets = 0;
            end
            
            dx = obj.phot.offsets_x(:,:,obj.phot_ap_idx)'; 
            dx(isnan(dx)) = 0; 
            dy = obj.phot.offsets_y(:,:,obj.phot_ap_idx)'; 
            dy(isnan(dy)) = 0; 
            
            if cs(obj.pars.lock_adjust, 'all') || obj.pars.use_astrometry==0 % by default, "stars" is replaced by "all" when not running astrometry
                obj.data.positions_new = double(obj.data.positions + mean_offsets);
            elseif cs(obj.pars.lock_adjust, 'none')
                obj.data.positions_new = double(obj.data.positions + [dx, dy]); % check dimensionality in case we have multiple apertures! 
            elseif cs(obj.pars.lock_adjust, 'stars') % assume astrometry has succeeded
                
                star_idx = ~isnan(obj.cat.magnitudes); % all stars that have a non NaN magnitude in GAIA
                
                star_idx = vertcat(star_idx, ones(size(obj.data.positions,1)-length(star_idx),1)); % if positions is longer than the catalog (e.g., forced photometry)
                
                new_pos_free = double(obj.data.positions + [dx, dy]); % check dimensionality in case we have multiple apertures! 
                
                new_pos_lock = double(obj.data.positions + mean_offsets); 
                
                obj.data.positions_new = star_idx.*new_pos_lock + (~star_idx).*new_pos_free; % logical indexing to add the correct adjusted position to stars and non-stars
                
            else
                error('Unknown "lock_adjust" option: "%s". Choose "all", "none" or "stars". ', obj.pars.lock_adjust); 
            end
            
            if size(obj.data.positions,3)>1
                error('dimensionality issue...'); 
            end
            
            if ~isempty(obj.forced_cutout_motion_per_file)
                obj.data.positions_new(obj.forced_indices(1), :) = obj.data.positions_new(obj.forced_indices(1), :) + obj.forced_cutout_motion_per_file;
            end
            
        end
        
        function realignPositions(obj) % use the reference image and positions to realign
            
            [~,shift] = util.img.quick_align(obj.data.image_proc, obj.ref_image);
            obj.data.positions_new = double(obj.ref_positions + flip(shift));
            
        end
        
        function val = checkFluxes(obj) % check that we can still see most of the stars, compared to the flux buffer
            
            if size(obj.flux_buf.data, 2)~=size(obj.phot.fluxes,2)
                obj.flux_buf.reset;
            end
            
            if is_empty(obj.flux_buf)
                val = 1;
            else
                
                mean_fluxes = obj.flux_buf.mean;
                mean_fluxes(mean_fluxes<=0) = NaN;

                new_fluxes = obj.phot.fluxes(:,:,obj.phot_ap_idx);
                
                % remove all points where the reference fluxes are NaN
                new_fluxes(isnan(mean_fluxes)) = [];
                mean_fluxes(isnan(mean_fluxes)) = [];

                flux_lost = sum(new_fluxes<0.5*mean_fluxes)>0.5*numel(mean_fluxes); % lost half the flux in more than half the stars...
                % want to add more tests...?

                val = ~flux_lost;
                
            end
            
        end
        
        function calculateFWHM(obj)

            t0 = tic;
            
            N_stars = size(obj.data.cutouts,4); 

            idx = obj.pars.fwhm_indices; 

            idx(idx>N_stars) = []; % remove indices outside the bounds

            if length(idx)<20 % arbitratry cutoff
                idx = 1:N_stars; % replace with just all stars instead
            end

            C = obj.data.cutouts_sub(:,:,:,idx); % cutouts chosen for FWHM calculation

            pix = obj.head.SCALE; 
            
            w = NaN(1,N_stars); % all stars that were not chosen for calculation are set to NaN
            w(idx) = util.img.fwhm(C, 'method', 'filters', 'defocus', obj.data.psf_width, ...                 
                'step_size', 0.1./pix, 'min_size', 0.5.*pix, 'max_size', obj.pars.cut_size).*obj.head.SCALE; 

            obj.data.fwhm = w;

            fr = util.fit.surf_poly(obj.data.positions(idx,1), obj.data.positions(idx,2), w(idx)', 'sigma', 3, 'order', 2, 'double', 1); 
            
            obj.data.fwhm_interp = fr.func(obj.data.positions(:,1), obj.data.positions(:,2));
            
            obj.data_logs.fwhm(obj.to_process_index,1) = nanmedian(w); 

            obj.addTimingData('fwhm', toc(t0)); 

        end
        
        function calculateLightcurves(obj)
            
            t0 = tic;

            obj.light.index_flux = obj.phot_ap_idx; % tell this object which aperture index we want to use
            obj.light.getData(obj.phot); 
            if obj.light.gui.check, obj.light.gui.update; end
    
            obj.addTimingData('lightcurves', toc(t0)); 

        end
        
        function calculateFunction(obj)
            
            t0 = tic; 

            if ~isempty(obj.func_handle)

                if isa(obj.func_handle, 'function_handle')
                    obj.func_handle(obj.data, obj.pars)
                elseif iscell(obj.func_handle)
                    for jj = 1:length(obj.func_handle)
                        obj.func_handle{ii}(obj.data, obj.pars)
                    end
                end
            end

            obj.addTimingData('custom_func', toc(t0)); 
            
        end
        
    end
    
    methods(Hidden=true) % internal utility functions
        
        function val = find_folder_date(obj)
            
            base_dir = obj.reader.current_dir;

            for ii = 1:3 % try to figure out this run's own date

                [base_dir, end_dir] = fileparts(base_dir);

                if isempty(end_dir), break; end

                [idx1,idx2] = regexp(end_dir, '\d{4}-\d{2}-\d{2}');
                if isempty(idx1), continue; end

                val = end_dir(idx1:idx2);

                if isempty(base_dir), break; end

            end
            
        end
        
        function val = find_camera(obj)
            
            val = '';
            
            if ~isempty(obj.reader.filenames)
                
                f = lower(obj.reader.filenames{1});

                if contains(f, {'balor'})
                    val = 'Balor';
                elseif contains(f, {'zyla'})
                    val = 'Zyla';
                end

            end
            
        end
        
        function val = find_project(obj)
            
            val = '';
            
            if ~isempty(obj.reader.filenames)
                
                f = lower(obj.reader.filenames{1});

                if contains(f, {'wfast'})
                    val = 'WFAST';
                elseif contains(f, {'kraar'})
                    val = 'Kraar';
                elseif contains(f, {'last'})
                    val = 'LAST';
                end

            end
            
        end
        
        function val = find_filter(obj)
            
        end
        
        function val = find_phot_ap(obj)
            
            if isempty(obj.data.psf_width)
                val = [];
            else
                a = floor(obj.data.psf_width.*3); 
                val = find(strcmp(obj.phot.pars_struct.types, sprintf('forced %4.2f', a))); 
            end
            
        end
        
        function addTimingData(obj, name, seconds)
            
            name = strtrim(name); 
            name = strrep(name, ' ', '_'); 
            
            if ~isfield(obj.timings, name)
                obj.timings.(name) = 0; 
            end
            
            obj.timings.(name) = obj.timings.(name) + seconds; 
            
        end
        
        function displayInfo(obj, str) % show the info on screen, write to log (if saving data) and show it on the GUI
            
            obj.write_log(util.text.date_printf(str)); 
            
            try 
                
                if ~isempty(obj.gui) && obj.gui.check
                    obj.gui.panel_info.button_info.String = str;
                end
            
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function write_log(obj, str)
            
            if obj.pars.use_save_results==0
                return; 
            end
            
            if nargin<2 || isempty(str)
                str = obj.file_summary;
            end
            
            filename = fullfile(obj.output_folder_name, 'process_log.txt'); 
            
            fid = fopen(filename, 'at');
            on_cleanup = onCleanup(@() fclose(fid));
            
            str = [util.text.time2str('now') ' : ' str]; 
            
            if isempty(str) || str(end)~=newline
                str(end+1) = '\n';
            end
            
            fprintf(fid, str); 
            
        end
        
        function str = file_summary(obj)
           
            str = sprintf('file: %d/%d', obj.to_process_index, obj.getNumFiles); 
            
            offsets = obj.getAverageOffsets; 
            fwhm = obj.getAverageFWHM; 
            
            if ~isempty(offsets)
                str = sprintf('%s | dx= %4.2f | dy= %4.2f | FWHM= %4.2f"', str, offsets(1), offsets(2), fwhm); 
            end
            
            if ~isempty(obj.head.AIRMASS)
                str = sprintf('%s | a.m.= %4.2f', str, obj.head.AIRMASS); 
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('autodyn', false); 
            input.scan_vars(varargin{:});
            
            if isempty(obj.data)
                return;
            end
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            I = obj.data.image_proc;
            
            % any preprocessing you want to do? 
            
            util.plot.setImage(I, input.ax);
            
            if input.autodyn
                input.ax.CLim = util.img.autodyn(I); 
            end
            
            delete(findobj(input.ax, 'type', 'rectangle'));
            delete(findobj(input.ax, 'tag', 'clip number'));
            
            if ~isempty(obj.data.positions)

                P = obj.data.positions;

                for ii = 1:min(size(P,1), obj.pars.display_num_rect_stars)
                    
                    if ii<=size(obj.stars,1) && obj.stars.flag(ii)==2
%                         rectangle('Position', [P(ii,:)-0.5-obj.pars.cut_size/2 obj.pars.cut_size obj.pars.cut_size], 'Parent', input.ax, 'EdgeColor', 'yellow'); 
                        c = 'yellow'; % yellow is for point sources
                    elseif ii<=size(obj.stars,1) && obj.stars.flag(ii)==3
%                         rectangle('Position', [P(ii,:)-0.5-obj.pars.cut_size/2 obj.pars.cut_size obj.pars.cut_size], 'Parent', input.ax, 'EdgeColor', 'green'); 
                        c = 'green'; % green is for extended sources
                    elseif ii<=size(obj.stars,1) && obj.stars.flag(ii)==1
%                         rectangle('Position', [P(ii,:)-0.5-obj.pars.cut_size/2 obj.pars.cut_size obj.pars.cut_size], 'Parent', input.ax, 'EdgeColor', 'red'); 
                        c = 'red'; % red is for saturated stars
                    elseif ii<=size(obj.stars,1) && ~isempty(obj.cat.magnitudes) && isnan(obj.cat.magnitudes(ii))
%                         rectangle('Position', [P(ii,:)-0.5-obj.pars.cut_size/2 obj.pars.cut_size obj.pars.cut_size], 'Parent', input.ax, 'EdgeColor', 'white'); 
                        c = 'white'; % white is for no match to GAIA
                    else    
%                         rectangle('Position', [P(ii,:)-0.5-obj.pars.cut_size/2 obj.pars.cut_size obj.pars.cut_size], 'Parent', input.ax, 'EdgeColor', 'black'); 
                        c = 'black'; % black is good stars
                    end
                    
                    if obj.pars.use_display_rect_text, text(P(ii,1), P(ii,2), sprintf('clip %d', ii),'FontSize', 16, 'Parent', input.ax, 'Color', c, 'Tag', 'clip number'); end
                    rectangle('Position', [P(ii,:)-0.5-obj.pars.cut_size/2 obj.pars.cut_size obj.pars.cut_size], 'Parent', input.ax, 'EdgeColor', c)
                    
                end
                
                for ii = 1:length(obj.forced_indices)
                    c = 'magenta'; 
                    idx = obj.forced_indices(ii); 
                    if obj.pars.use_display_rect_text, text(P(idx,1), P(idx,2), sprintf('forced %d', ii),'FontSize', 16, 'Parent', input.ax, 'Color', c, 'Tag', 'clip number'); end
                    rectangle('Position', [P(idx,:)-0.5-obj.pars.cut_size/2 obj.pars.cut_size obj.pars.cut_size], 'Parent', input.ax, 'EdgeColor', c); 
                    rectangle('Position', [P(idx,:)-0.5-obj.pars.cut_size obj.pars.cut_size*2 obj.pars.cut_size*2], 'Parent', input.ax, 'EdgeColor', c, 'LineWidth', 3); 
                    
                end
                
            end
            
        end
        
        function showFindingChart(obj, varargin)
            
        end
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = img.gui.ProcGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
    end
    
end

 