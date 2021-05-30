classdef Analysis < file.AstroData
% Analyze existing HDF5 files, searching for KBOs, and creating lightcurves. 
% 
% Choose a folder with some raw cutouts/stacks using chooseDir() function 
% or using the "browse" button on the GUI. 
% 
% Start a new analysis run using run(...), optionally with a call to reset()
% before you start. This will clear all data in case the object was already
% used in a different analysis run. Without this, it will continue the same
% analysis run from where it left off (this is a little buggy still). 
% 
% The run() command accepts some arguments:
%   -reset: equivalent to calling reset() at the start of the run. 
%   -logging: create an analysis folder with a text file loggint the progress
%             of the analysis run. 
%   -save: create an analysis folder and save the results of the run. 
%   -overwrite: if false, will raise an error if trying to create an 
%               analysis folder with the same date as an existing one. 
%               This prevents running the same analysis in the same day. 
%               If a previous run was cut short and you want to restart it, 
%               use this option. Make sure the previous run was stopped!
%
% NOTE: all these arguments are false by default. 
%
% There are actually a lot more parameters that can be changed in this 
% object. They are specified in the "switches/controls" block. 
% Many other parameters for specific searches are saved inside the relevant
% objects, e.g., the trig.EventFinder object "finder" defines many parameters
% for the KBO/Oort cloud search. 
% 
% The run starts with a call to startup() and ends with finishup(). 
% Each batch in the run is called by batch() which loads a file or other
% dataset (see the "reader" object) and then runs all sorts of analysis 
% functions (all named analysisXXX). This includes making stacks and cutouts, 
% calibrating them, running photometry on the cutouts, finding KBOs, etc. 
%
% Another option is to use async_run(). This is the same as run() only it
% creates a copy of the Analysis object and sends it to a parallel worker. 
% This function accepts the same arguments as run() but also lets you choose
% the "worker" option, for specifying which worker to use out of the pool. 
% The GUI button for "run async" also pops up a window for choosing these
% parameters, and it automatically finds an available worker. 
%
% This class inherits from file.AstroData which defines all the basic data
% properties such as "images", "cutouts", etc. It also defines the copyFrom()
% that copies all these data from the reader. 
%
    
    properties(Transient=true)
        
        gui;
        audio@util.sys.AudioControl;
        
        pool; 
        futures; 
        futures_dir;
        futures_analysis_folder;
        futures_batches;
        
        monitor@img.gui.FutureMonitor; 
        
        aux_figure;
        
    end
    
    properties % objects
        
        head@head.Header;
        cat@head.Catalog;
        reader@file.Reader;
        cal@img.Calibration;
        clip@img.Clipper;
        clip_bg@img.Clipper;
        back@img.Background;
        phot@img.Photometry;
        phot_stack@img.Photometry;
        flux_buf@util.vec.CircularBuffer;
        mean_buf@util.vec.CircularBuffer;
        var_buf@util.vec.CircularBuffer;
        back_buf@util.vec.CircularBuffer;
        width_buf@util.vec.CircularBuffer;
        
        lightcurves@img.Lightcurves;
        light_stack@img.Lightcurves; 
        buf@file.BufferWheel; 
        
        model_psf@img.ModelPSF;
        
        finder@trig.EventFinder;
        
        sky_pars;
        
        cutout_store@learn.CutoutStorage; 
        
        prog@util.sys.ProgressBar;
        
        image_mextractor;
        matched_gaia;
        catalog;
        
        func; % any function that takes first argument this object and runs custom analysis
        
    end
    
    properties % inputs/outputs
        
        cutouts_proc;
        cutouts_sub;
        cutouts_bg_proc;
        
        stack_cutouts; 
        stack_cutouts_sub;
        stack_cutouts_bg;
        stack_proc;
        deep_stack;
        deep_stack_aligned;
        subtract_stack;
        
        prev_stack;
        
        FWHM; % latest measured full width half maximum (in arcsec)
        
        FWHM_log; % keep track of the FWHM over entire run
        juldates_full; % keep a record of the juldates across entire run
        
        batch_counter = 0;
        
    end
    
    properties % switches/controls
        
        num_stars = []; % if this is empty, use all the stars
        cut_size = 21;
        saturation_value = 50000; % consider any pixels above this to be saturated
        
        use_background_stack = 1; % subtract b/g from the full-frame stack
        use_background_cutouts = 0; % subtract b/g from the cutouts (and stack cutouts!)
        use_refine_bg = 0; % need to figure out exactly how to do this
        
        use_check_flux = 1;
        max_failed_batches = 3; % if star flux is lost for more than this number of batches, quit the run
        
        use_astrometry = 1;
        use_save_astrometry = 1; % save the new astrometry, updating the catalog file
        use_require_astrometry = 1; % when true, will error if there is no astrometric solution
        
        use_cutouts = 1;
        use_photometry = 1;
        
        use_analysis_dir_save = 0; % do we want to save analysis results
        use_analysis_dir_log = 0; % do we want analysis log file
        
        use_full_lightcurves = 0; % use the Lightcurves object to hold all the fluxes and other measurements for the entire run
        use_save_full_lightcurves = 1; % save these full lightcurves for the entire run as a single MAT file (ONLY when use_full_lightcurves=1 and use_analysis_dir_save=1)
        use_save_batched_lightcurves = 1; % save each batch's photometric result in a separate file (ONLY when use_analysis_dir_save=1)
        use_stack_lightcurves = 1; % keep a lightcurves object for the stack images
        
        use_psf_model = 1;
        
        use_event_finding = 1;
        
        use_auto_load_cal = 1;
        
        use_fits_save = 0;
        use_fits_flip = 0;
        use_fits_roi = 0;
        fits_roi = []; % [top, left, height, width]
        
        use_max_bad_batches = 1; 
        max_num_bad_batches = 10; % if more than this number of batches have high background or bad focus, stop the analysis
        
        use_audio = 0;
        
        use_display_flip = 0;
        display_num_rect_stars = 30;
        
        use_fwhm_stop = false;
        max_fwhm_stop = 10; % if seeing is worse than this number (in arcsec) quit the run
        
        brake_bit = 1;
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        filename;
        directory;
        num_batches;
        average_width;
        average_offsets;
        
        num_stars_used; % min of num_stars and number of cutouts
        
    end
    
    properties(Hidden=true)
       
        max_num_workers = 10; 
        
        ref_stack; % for quick align
        ref_positions; 
        
        failed_batch_counter = 0;
        
        use_cutouts_all = 0;
        use_cutouts_all_proc = 1;
        use_cutouts_store = 0; 
        
        use_duplicate_filter = 1;
        duplicate_batches = [];
        
        cutouts_all;
        positions_x_all;
        positions_y_all;
        
        use_stack_all = 0;
        use_stack_all_proc = 1;
        stack_all;
        
        num_bad_batches = 0; 
        
        % these are used for backward compatibility with older versions of 
        % img.Clipper that made a slightly different cutout based on the 
        % same positions. This affects ONLY the cutting of CALIBRATION! 
        use_cutout_adjustment = 0; % turn adjustments on/off
        cutout_adjustment_pixels = 0; % how many pixels to push (back or forward) relative to today's positions
        use_cutout_adjustment_floor = 0; % use floor of positions before (instead of) using round(). 
        
        log_dir; % where to save the analysis results and log file
        log_name; % name of log file
        log_obj; % name of object text dump
        
        prev_average_width;
        
        num_batches_limit;
        
        version = 1.06;
        
    end
    
    methods % constructor
        
        function obj = Analysis(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'img.Analysis')
                if obj.debug_bit>1, fprintf('Analysis copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('Analysis constructor v%4.2f\n', obj.version); end
            
                obj.reader = file.Reader;
                obj.cal = img.Calibration;
%                 obj.cal.load;
                obj.clip = img.Clipper;
                obj.clip.use_adjust = 0; % this should be disabled and depricated!
                obj.clip_bg = img.Clipper;
                obj.clip_bg.use_adjust = 0; % this should be disabled and depricated!
                obj.back = img.Background;
                obj.phot = img.Photometry;
                obj.phot.aperture = [3 5 7]; 
                obj.phot.index = 1;
                
                obj.phot_stack = img.Photometry;
                obj.phot_stack.aperture = [3 5 7]; 
                obj.phot_stack.saturation_value = 5e6; 
                obj.phot_stack.index = 2;
                
                obj.flux_buf = util.vec.CircularBuffer;
                obj.mean_buf = util.vec.CircularBuffer;
                obj.var_buf = util.vec.CircularBuffer;
                obj.back_buf = util.vec.CircularBuffer;
                obj.width_buf = util.vec.CircularBuffer;
                
                obj.lightcurves = img.Lightcurves; 
                obj.light_stack = img.Lightcurves;
                obj.buf = file.BufferWheel;
                obj.buf.use_async = 0;
                obj.buf.use_deflate = 1;
                
                obj.model_psf = img.ModelPSF;
                
                obj.finder = trig.EventFinder;
%                 obj.finder.loadFilterBank;
                
                obj.cutout_store = learn.CutoutStorage; 
                
                obj.prog = util.sys.ProgressBar;
                obj.audio = util.sys.AudioControl;
                
                obj.head = head.Header; % this also gives "head" to all sub-objects
                obj.cat = head.Catalog(obj.head);
                obj.finder.cat = obj.cat;
                obj.finder.head = obj.head;
                obj.finder.store.checker.setupSensor; 
%                 obj.lightcurves.head = obj.head;
                obj.lightcurves.cat = obj.cat;
                obj.light_stack.cat = obj.cat;
                
                util.oop.save_defaults(obj); % make sure each default_XXX property is updated with the current XXX property value. 
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            list = properties(obj);
            
            for ii = 1:length(list)
                
                if isobject(obj.(list{ii})) && ~isempty(obj.(list{ii})) && ismethod(obj.(list{ii}), 'reset') 
                    obj.(list{ii}).reset;
                end
                
            end
            
            obj.batch_counter = 0;

            obj.clear;
            
            obj.prev_stack = [];
            obj.ref_stack = [];
            obj.ref_positions = [];
            
            obj.use_analysis_dir_save = 0;
            obj.use_analysis_dir_log = 0;
            
            obj.failed_batch_counter = 0;
            
            obj.cutouts_all = [];
            obj.positions_x_all = [];
            obj.positions_y_all = [];
            obj.stack_all = [];
            
            obj.duplicate_batches = []; 
            
            obj.finder.store.checker.setupSensor; % make sure the quality checker has the right bad rows/columns
            
            obj.FWHM_log = []; 
            obj.juldates_full = [];
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            list = properties(obj);
            
            for ii = 1:length(list)
                
                if isobject(obj.(list{ii})) && ~isempty(obj.(list{ii})) && ismethod(obj.(list{ii}), 'clear') 
                    obj.(list{ii}).clear;
                end
                
            end
            
            clear@file.AstroData(obj);
            
            obj.cutouts_proc = [];
            obj.cutouts_sub = [];
            obj.cutouts_bg_proc = [];
            obj.stack_cutouts = [];
            obj.stack_cutouts_sub = [];
            obj.stack_cutouts_bg = [];
            obj.stack_proc = [];
            
            obj.sky_pars = [];
            
            obj.FWHM = []; 
            
        end
        
    end
    
    methods % getters
        
        function val = get.filename(obj)
            
            if isempty(obj.reader) || isempty(obj.reader.prev_filename)
                val = '';
            else
                [~,file, ext] = fileparts(obj.reader.prev_filename);
                val = [file, ext];
            end
            
        end
        
        function val = get.directory(obj)
           
            if isempty(obj.reader) || isempty(obj.reader.prev_filename)
                val = '';
            else
                val = fileparts(obj.reader.prev_filename);
            end
            
        end
        
        function val = get_obj_name(obj)
            
            if ~isempty(obj.head) && ~isempty(obj.head.OBJECT)
                val = obj.head.OBJECT;
            else
                [~, val] = fileparts(obj.directory); 
                
                idx = regexp(val, '_run\d+$');
                if ~isempty(idx) && idx>1
                    val = val(1:idx-1);
                end
                
            end
            
        end
        
        function val = get.num_batches(obj)
            
            if isempty(obj.reader) && isempty(obj.num_batches_limit)
                val = [];
            elseif ~isempty(obj.num_batches_limit)
                val = obj.num_batches_limit;
            elseif ~isempty(obj.reader)
                val = obj.reader.getNumBatches;
            else
                val = min([obj.num_batches_limit, obj.reader.getNumBatches]);
            end
            
        end
        
        function val = seeing(obj)
            
            if isempty(obj.head)
                val = obj.FWHM.*1.24;
            else
                val = obj.FWHM.*obj.head.SCALE;
            end
            
        end
        
        function val = thisFilename(obj)
            
            if isempty(obj.reader)
                val = '';
            else
                val = obj.reader.prev_filename;
            end
            
        end
        
        function val = getWidthEstimate(obj)
            
            if isempty(obj.prev_average_width) || ~isreal(obj.prev_average_width) || obj.prev_average_width<=0
                val = 1.6;
            else
                val = obj.prev_average_width;
            end
            
        end
        
        function val = get.average_width(obj)
            
            if isempty(obj.model_psf)
                val = [];
            else
                val = (obj.model_psf.maj_axis+obj.model_psf.min_axis)/2;
            end
            
        end
        
        function val = get.average_offsets(obj)
            
            if isempty(obj.phot_stack)
                val = [];
            else
                val = [obj.phot_stack.average_offset_x obj.phot_stack.average_offset_y];
            end
            
        end
        
        function val = get.num_stars_used(obj)
            
            if isempty(obj.num_stars)
                val = size(obj.cutouts, 4);
            else
                val = min(obj.num_stars, size(obj.cutouts,4));
            end
            
        end
        
    end
    
    methods % setters
        
        function set.head(obj,val)
            
            obj.head = val;
            
            list = properties(obj);
            
            for ii = 1:length(list)
                
                try
                    if isobject(obj.(list{ii})) && ~isempty(obj.(list{ii})) && ~istable(obj.(list{ii})) && isprop(obj.(list{ii}), 'head') 
                        obj.(list{ii}).head = val;
                    end
                catch ME
                    
                    disp(['trouble setting "head" into variable: ' list{ii}]);
                    warning(ME.getReport);
                    
                end
                
            end
            
        end
        
        function set.brake_bit(obj, val)
            
            obj.brake_bit = val;
            obj.reader.brake_bit = val;

        end
    
        function set.num_batches(obj, val)
            
            obj.num_batches_limit = val;
            
        end
        
    end
    
    methods % utilities
        
        function chooseDir(obj, dirname)
            
            if nargin<2 || isempty(dirname)
                dirname = '';
            end
            
            if isempty(dirname)
                obj.reader.browseDir;
            else
                if ~obj.reader.dir.cd(dirname)
                    error('cannot find directory %s', dirname);
                end
            end
            
%             obj.cal.reader_dark.dir.cd(r.dir.pwd);
%             obj.cal.reader_dark.dir.cd('..');
%             if obj.cal.reader_dark.dir.smart_cd('dark')
%                 obj.cal.reader_dark.loadFiles;
%             end
%             
%             obj.cal.reader_flat.dir.cd(r.dir.pwd);
%             obj.cal.reader_flat.dir.cd('..');
%             if obj.cal.reader_flat.dir.smart_cd('flat')
%                 obj.cal.reader_flat.loadFiles;
%             end
%             
%             obj.cal.load; 
            
        end
        
        function write_log(obj, filename, str)
            
            if nargin<3 || isempty(str)
                str = '';
            end
            
            fid = fopen(filename, 'at');
            on_cleanup = onCleanup(@() fclose(fid));
            
            obs_date = obj.head.STARTTIME;
            read_date = util.text.time2str(datetime('now', 'TimeZone', 'UTC'));
            
            f = [0 0 0];
            if ~isempty(obj.fluxes), f(1) = obj.fluxes(1,1); end
            if size(obj.fluxes,2)>=10, f(2) = obj.fluxes(1,10); end
            if size(obj.fluxes,2)>=100, f(3) = obj.fluxes(1,100); end

            fprintf(fid, '%s: ', read_date);
            fprintf(fid, 'Batch: %04d, ObsDate: %s, Flux: [% 9.1f % 8.1f % 7.1f]', obj.batch_counter+1, obs_date, f(1), f(2), f(3));
            fprintf(fid, '%s\n', str); 
            
%             if isempty(str)
%                 ev_str = '';
% 
%                 for ii = 1:length(obj.finder.last_events)
% 
%                     if obj.finder.last_events(ii).keep
%                         star_str = '*';
%                     else
%                         star_str = '';
%                     end
% 
%                     ev_str = sprintf('%s%4.2f%s ', ev_str, obj.finder.last_events(ii).snr, star_str);
% 
%                 end
% 
%                 f = [0 0 0];
%                 if ~isempty(obj.fluxes), f(1) = obj.fluxes(1,1); end
%                 if size(obj.fluxes,2)>=10, f(2) = obj.fluxes(1,10); end
%                 if size(obj.fluxes,2)>=100, f(3) = obj.fluxes(1,100); end
% 
%                 fprintf(fid, 'Batch: %04d, ObsDate: %s, Flux: [% 9.1f % 8.1f % 7.1f]', obj.batch_counter+1, obs_date, f(1), f(2), f(3));
%                 
%                 if ~isempty(obj.sky_pars)
%                     if isfield(obj.sky_pars, 'zero_point'), zp = obj.sky_pars.zero_point; else, zp = NaN; end
%     %                 if isfield(obj.sky_pars, 'noise_level'), nl = obj.sky_pars.noise_level; else, nl = NaN; end
%                     if isfield(obj.sky_pars, 'limiting_mag'), lm = obj.sky_pars.limiting_mag; else, lm = NaN; end
% 
%                     fprintf(fid, ' | seeing: %4.2f" | back: %5.3f | area: %4.2f | zp: %6.4g | lim. mag: %4.2f', ...
%                         obj.sky_pars.seeing, obj.sky_pars.background, obj.sky_pars.area, zp, lm);
%                 end
%                 
%                 fprintf(fid, ' | Events S/N: [%s], ReadDate: %s\n', ev_str, read_date);
%             
%             else
%                 fprintf(fid, '%s: %s\n', read_date, str);
%             end
            
            % if there is no object-dump file, create one now! 
            if ~exist(obj.log_obj, 'file')
                util.oop.save(obj, obj.log_obj, 'hidden', 1); 
            end
            
        end
        
        function saveResults(obj)
            
            if ~exist(obj.log_dir, 'dir')
                mkdir(obj.log_dir); % if we call this function we are ignoring "overwrite analysis folder" mechanism (from start of run() function) and creating a folder if needed!~
            end
            
            name = obj.get_obj_name;
            
            try % save the full lightcurves
                if obj.use_full_lightcurves && obj.use_save_full_lightcurves
                    obj.lightcurves.saveAsMAT(fullfile(obj.log_dir, ['lightcurves_' name]));
                end
            catch ME
                warning(ME.getReport);
            end
            
            try % save the stack lightcurves
                if obj.use_stack_lightcurves 
                    obj.light_stack.saveAsMAT(fullfile(obj.log_dir, ['light_stack_' name]));
                end
            catch ME
                warning(ME.getReport);
            end
            
            try % save the event finder
                
                obj.finder.finishup; 

                summary = obj.finder.produceSummary;
                
                save(fullfile(obj.log_dir, 'summary.mat'), 'summary', '-v7.3'); 
                
                try 
                    util.oop.save(summary, fullfile(obj.log_dir, 'summary.txt')); 
                catch ME
                    warning(ME.getReport); 
                end

                cand = obj.finder.cand; 
                
                save(fullfile(obj.log_dir, 'candidates.mat'), 'cand', '-v7.3'); 
                
                finder = obj.finder;

%                 save(fullfile(obj.log_dir, ['finder_' name]), 'finder', '-v7.3');
                save(fullfile(obj.log_dir, 'finder.mat'), 'finder', '-v7.3');
                
            catch ME
                warning(ME.getReport);
            end
            
            try % save a text file with any duplicates
                if ~isempty(obj.duplicate_batches)
                    fid = fopen(fullfile(obj.log_dir, 'duplicate_batches.txt'), 'wt');
                    on_cleanup = onCleanup(@() fclose(fid));
                    fprintf(fid, '%d\n', obj.duplicate_batches);
                end
            catch ME
                warning(ME.getReport);
            end
            
            % if there is no object-dump file, create one now! 
            if ~exist(obj.log_obj, 'file')
                util.oop.save(obj, obj.log_obj, 'hidden', 1); 
            end
            
        end
        
        function saveSummary(obj)
            
            if ~exist(obj.log_dir, 'dir')
                mkdir(obj.log_dir); % if we call this function we are ignoring "overwrite analysis folder" mechanism (from start of run() function) and creating a folder if needed!~
            end
            
%             filename = ['summary_' obj.get_obj_name '.txt'];
%             
%             fid = fopen(fullfile(obj.log_dir, filename), 'wt');
%             
%             if fid<0
%                 warning('Cannot open file %s', fullfile(obj.log_dir, filename));
%             else
% 
%                 onc = onCleanup(@() fclose(fid));
% 
%                 fprintf(fid, 'Summary for run %s, with %d batches.\n', obj.head.OBJECT, obj.batch_counter);
% 
%                 v = abs(obj.finder.snr_values);
% 
%                 fprintf(fid, 'S/N for the last %d batches is distrubuted: min= %f median= %f max= %f\n', numel(v), nanmin(v), nanmedian(v), nanmax(v));
% 
%                 v = abs([obj.finder.cand.snr]);
% 
%                 fprintf(fid, 'S/N for %d triggered events is distrubuted: min= %f median= %f max= %f\n', numel(v), nanmin(v), nanmedian(v), nanmax(v));
% 
%                 v = abs([obj.finder.kept.snr]);
% 
%                 fprintf(fid, 'S/N for %d kept events is distrubuted: min= %f median= %f max= %f\n', numel(v), nanmin(v), nanmedian(v), nanmax(v));
% 
%                 fprintf(fid, 'Number of events: total= %d | kept= %d\n', length(obj.finder.cand), length(obj.finder.kept));
%             
% %                 fprintf(fid, 'Star hours (above stellar S/N of %4.2f): %4.2f \n', obj.finder.min_star_snr, obj.finder.star_hours_total);
% %                 
% %                 fprintf(fid, 'Star hours (above stellar S/N of %4.2f): %4.2f \n', obj.finder.min_star_snr*2, obj.finder.star_hours_total_better);
% %                 
% %                 fprintf(fid, 'Star hours (above stellar S/N of %4.2f): %4.2f \n', obj.finder.min_star_snr*4, obj.finder.star_hours_total_best);
%                 
%             end
            
            % if there is no object-dump file, create one now! 
            if ~exist(obj.log_obj, 'file')
                util.oop.save(obj, obj.log_obj, 'hidden', 1); 
            end
            
        end
        
        function print_futures(obj)
            
            for ii = 1:length(obj.futures)
                
                if ~isempty(obj.futures{ii}) && isa(obj.futures{ii}, 'parallel.Future') && isvalid(obj.futures{ii})
                
                    asterisk = ' ';
                    if obj.futures{ii}.Read==0
                        asterisk = '*';
                    end
                    
                    finish_datetime = obj.futures{ii}.FinishDateTime;
                    if isempty(finish_datetime)
                        finish_datetime = datetime('now', 'TimeZone', 'Local');
                        asterisk = ' ';
                    end
                    
                    fprintf('Future{%2d}: State= %14s%s | Error= %d | runtime= %9s', ii, obj.futures{ii}.State, asterisk, ~isempty(obj.futures{ii}.Error),...
                        char(finish_datetime-obj.futures{ii}.StartDateTime));
                
                    if length(obj.futures_batches)>=ii && ~isempty(obj.futures_batches{ii})
                        fprintf(' | N= %4d', obj.futures_batches{ii});
                    end
                    
                    if length(obj.futures_dir)>=ii && ~isempty(obj.futures_dir{ii})
                        fprintf(' | dir= %s', obj.futures_dir{ii});
                    end
                    
                    if strcmp(obj.futures{ii}.State, 'running')
                        fprintf('  <-----(running)---- \n');
                    elseif ~isempty(obj.futures{ii}.Error)
                        fprintf('      !!! Error !!! \n');
                    else
                        fprintf('\n');
                    end
                end
                
                
                
            end
            
        end
        
        function list = list_analysis_dirs(obj, directory)
            
            if nargin<2 || isempty(directory)
                directory = fullfile(getenv('DATA'), 'WFAST'); 
            end
            
            d = util.sys.WorkingDirectory(directory);
            
            full_list = d.walk; 
            
            idx = ~cellfun(@isempty, regexp(full_list, 'analysis_\d{4}-\d{2}-\d{2}$'));
            
            full_list = full_list(idx);
            
            list = {};
            
            prev_path = '';
            prev_date = NaT;
            
            for ii = 1:length(full_list)
                
                [leading_path,analysis_dir] = fileparts(full_list{ii});
                
                date_new = datetime(analysis_dir(10:19), 'Format', 'yyyy-MM-dd');
                
                if strcmp(leading_path, prev_path) % if this analysis folder is for the same dataset as previous one
                    if date_new>prev_date % update newer analysis folder
                        list{end} = full_list{ii};
                    end
                else % new dataset, just add this analysis folder to end of list
                    list{end+1} = full_list{ii};
                end

                prev_path = leading_path; 
                prev_date = date_new;

            end
            
            list = list';
            
        end
        
        function [star_hours, hour_dist] = summarize_star_hours(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('threshold', 5); 
            input.input_var('directory', '');
            input.scan_vars(varargin{:});
            
            list = obj.list_analysis_dirs(input.directory);
            
            hour_dist = [];
            
            d = util.sys.WorkingDirectory;
            
            for ii = 1:length(list)
                
                d.cd(list{ii});
                filename = d.match('summary*');
                
                if ~isempty(filename)
                    
                    f = fopen(filename{1}, 'r');
%                     oc = onCleanup(@() fclose(f));
                    
                    if f<0
                        fclose(f);
                        continue;
                    end
                    
                    for jj = 1:100
                        
                        line = fgetl(f);
                        
                        if isnumeric(line) && line<0
                            fclose(f);
                            break;
                        end
                        
                        if ~isempty(regexpi(line, 'star hours'))
                            
                            c = util.text.extract_numbers(line);
                            
                            thresh = c{1}(1);
                            hours = c{1}(end); 
                            
                            if isempty(hour_dist)
                                hour_dist = [thresh hours];
                                continue;
                            end
                            
                            idx = hour_dist(:,1)==thresh;
                            
                            if any(idx)
                                hour_dist(idx,2) = hour_dist(idx,2) + hours;
                            else
                                hour_dist(end+1,:) = [thresh hours];
                            end
                            
                        end
                        
                    end
                    
                end
                                
            end
            
            [mn,idx] = min(hour_dist(:,1)); 
            
            star_hours = hour_dist(idx,2); 
            
        end
        
        function auto_fits_roi(obj, roi_size)
            
            if nargin<2 || isempty(roi_size)
                roi_size = 500; 
            end
            
            if isscalar(roi_size)
                roi_size = roi_size.*[1 1];
            end
            
            xy = obj.head.wcs.coo2xy(obj.head.RA, obj.head.DEC);

            x1 = round(xy(1)-roi_size(2)/2); 
            y1 = round(xy(2)-roi_size(1)/2); 

            if x1<1
                x1 = 1;
            end

            if y1<1
                y1 = 1;
            end

            obj.fits_roi = [y1 x1 roi_size(1), roi_size(2)]; 
            
        end
        
    end
    
    methods % calculations
        
        function startup(obj)
            
            if isempty(obj.cat.data) && (isempty(obj.use_astrometry) || obj.use_astrometry==0) % in case use_astrometry==1 we will redo the astrometry anyway
                % try to get the catalog file
                filename = fullfile(obj.reader.dir.pwd, 'catalog.mat');
                if exist(filename, 'file')
                    obj.cat.loadMAT(filename)
                end
            end
            
            obj.brake_bit = 0;
            
            if obj.use_audio
                try obj.audio.playTakeForever; catch ME, warning(ME.getReport); end
            end
            
            if obj.batch_counter==0
                obj.prog.reset(obj.num_batches);
            end
            
            obj.prog.unpause;
            
        end
        
        function finishup(obj)
            
            obj.prog.finish;
            
            obj.finder.finishup;
            
            obj.brake_bit = 1;
            
            if obj.use_audio
                try
                    obj.audio.playShowsOver;
                catch ME
                    warning(ME.getReport);
                end
            end
            
            if obj.use_full_lightcurves
                obj.lightcurves.finishup;
            end
            if obj.debug_bit, disp(['Finished run with ' num2str(obj.batch_counter) ' batches.']); end
            
%             diary('off'); 
            
        end
        
        function async_run(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('reset', 0); % reset the object before running (i.e., start a new run)
            input.input_var('logging', []); % create log files in the analysis folder
            input.input_var('save', []); % save the events and lightcurves from this run
            input.input_var('overwrite', 0); % delete the existing analysis folder without warning (make sure there is no ongoing analysis!)
            input.input_var('worker', []); % index of worker/future you want to use
            input.input_var('output', true); % choose if you want to get back an Analysis object as output from the calculation
            input.scan_vars(varargin{:});
            
            if isempty(input.worker)
                input.worker = obj.findWorker;
            end
            
            obj.futures{input.worker} = parfeval(obj.pool, @obj.run, double(input.output), 'reset', input.reset, 'logging', input.logging, 'save', input.save, 'overwrite', input.overwrite); 
            obj.futures_dir{input.worker} = obj.reader.dir.two_tail;
            obj.futures_analysis_folder{input.worker} = fullfile(obj.reader.dir.pwd, ['analysis_' char(datetime('now', 'TimeZone', 'UTC'), 'yyyy-MM-dd')]);
            obj.futures_batches{input.worker} = obj.num_batches; 
            
        end
        
        function idx = findWorker(obj, varargin)
            
            if isempty(gcp('nocreate'))
                obj.pool = parpool;
                obj.pool.IdleTimeout = 360;
            end
            
            N = obj.pool.NumWorkers; 
            
            if ~isempty(obj.max_num_workers) && obj.max_num_workers<N
                N = obj.max_num_workers; 
            end
            
            idx = [];
            
            for ii = 1:N
                
                if length(obj.futures)<ii || ~isa(obj.futures{ii}, 'parallel.Future') || ~isvalid(obj.futures{ii})...
                        || ( strcmp(obj.futures{ii}.State, 'finished') && obj.futures{ii}.Read==1)
                    idx = ii;
                    return;
                end
                
            end
            
            if ii==N
                error('Cannot find a free worker to run analysis...');
            end
            
            
        end
        
        function idx = findWorkerUnread(obj)
            
            if isempty(gcp('nocreate'))
                obj.pool = parpool;
                obj.pool.IdleTimeout = 360;
            end
            
            N = obj.pool.NumWorkers; 
            
            if ~isempty(obj.max_num_workers) && obj.max_num_workers<N
                N = obj.max_num_workers; 
            end
            
            idx = [];
            
            for ii = 1:N
                
                if length(obj.futures)<ii || ~isa(obj.futures{ii}, 'parallel.Future') || ~isvalid(obj.futures{ii})...
                        || (strcmp(obj.futures{ii}.State, 'finished') && isempty(obj.futures{ii}.Error))
                    
%                     if strcmp(obj.futures{ii}.State, 'finished')
%                         
%                         if ~isempty(obj.futures_dir{ii}) % we know the analysis folder, we can add a diary file to it
%                             
%                         end
%                         
%                     end
                    
                    idx = ii;
                    return;
                end
                
            end
            
        end
        
        function clearWorkerErrors(obj, number)
            
            if nargin<2
                number = 1:length(obj.futures);
            end
            
            for ii = number
                
                if ii<=length(obj.futures)
                    
                    if ~isa(obj.futures{ii}, 'parallel.Future') || ~isvalid(obj.futures{ii})...
                        || (strcmp(obj.futures{ii}.State, 'finished') && ~isempty(obj.futures{ii}.Error))
                    
                        obj.futures{ii} = []; 
                    
                    end
                    
                end
                
            end
            
        end
        
        function obj_out = run(obj, varargin)
            
            try
            
                input = util.text.InputVars;
                input.input_var('reset', false); % reset the object before running (i.e., start a new run)
                input.input_var('logging', []); % create log files in the analysis folder
                input.input_var('save', []); % save the events and lightcurves from this run
                input.input_var('overwrite', false); % if true, will quietly delete previous analysis folder from current date (otherwise, throws an error)
%                 input.input_var('diary', false); % save standard output to file
                input.scan_vars(varargin{:});

                if input.reset
                    obj.reset;
                end

                if obj.use_auto_load_cal

                    date = datetime.empty;

                    base_dir = obj.reader.current_dir;

                    for ii = 1:3 % try to figure out this run's own date

                        [base_dir, end_dir] = fileparts(base_dir);

                        if isempty(end_dir), break; end

                        [idx1,idx2] = regexp(end_dir, '\d{4}-\d{2}-\d{2}');
                        if isempty(idx1), continue; end

                        date = datetime(end_dir(idx1:idx2));

                        if isempty(base_dir), break; end

                    end

                    camera = 'Zyla';
                    project = 'WFAST'; 
                    
                    filenames = obj.reader.dir.match('*.h5*'); 
                    
                    if ~isempty(filenames)
                        
                        f = lower(filenames{1}); 
                        
                        if contains(f, {'balor'})
                            camera = 'Balor';
                        elseif contains(f, {'zyla'})
                            camera = 'Zyla';
                        end
                        
                        if contains(f, {'wfast', 'w-fast', 'w_fast'})
                            project = 'WFAST';
                        elseif contains(f, {'kraar'})
                            project = 'Kraar';
                        end
                        
                    end
                    
                    if ~isempty(date)
                        obj.cal.loadByDate(datestr(date, 'yyyy-mm-dd'), camera, project, 0); % last argument is to NOT reload if date is consistent
                    end

                end

                if ~obj.cal.checkDark
                    error('Cannot start a new run without loading darks into calibration object!');
                end

                % update hidden variables in case we use GUI to stop then continue this run
                if ~isempty(input.logging), obj.use_analysis_dir_log = input.logging; end            
                if ~isempty(input.save), obj.use_analysis_dir_save = input.save; end

                % must be ready to save the analysis results event if the
                % run started without this mode tunred on... 
                log_time = datetime('now', 'TimeZone', 'UTC');
                obj.log_dir = fullfile(obj.reader.dir.pwd, ['analysis_' char(log_time, 'yyyy-MM-dd')]);
                obj.log_name = fullfile(obj.log_dir, 'analysis_log.txt');
                obj.log_obj = fullfile(obj.log_dir, 'analysis_parameters.txt');

                if obj.use_analysis_dir_log || obj.use_analysis_dir_save
                    
                    % must check pre existing analysis folder for the same
                    % date, and either overwrite or throw an error
                    if obj.batch_counter==0
                        if exist(obj.log_dir, 'dir') && input.overwrite==0
                            error('Folder %s already exists! Is this run already being processed?', obj.log_dir);
                        elseif exist(obj.log_dir, 'dir') && input.overwrite
                            rmdir(obj.log_dir, 's'); 
                            mkdir(obj.log_dir);
                        else
                            mkdir(obj.log_dir);
                        end
                    end

                end

%                 if obj.use_analysis_dir_log && input.diary
%                     diary(fullfile(obj.log_dir, 'diary.txt')); 
%                 end
                
                cleanup = onCleanup(@() obj.finishup);
                obj.startup;

                for ii = obj.batch_counter+1:obj.num_batches

                    if obj.brake_bit
                        break;
                    end

                    obj.batch;

                    if obj.use_analysis_dir_log
                        obj.write_log(obj.log_name);
                    end

                    if ~isempty(obj.func)

                        if isa(obj.func, 'function_handle')
                            feval(obj.func, obj);
                        elseif iscell(obj.func)
                            for jj = 1:length(obj.func)
                                feval(obj.func{jj}, obj);
                            end
                        end
                    end
                    
                    obj.prog.showif(ii);

                    drawnow;

                    obj.batch_counter = obj.batch_counter + 1;

                end

                % Skip this part in case the run is stopped mid way (e.g., by user input, but not by detecting flux is lost)
                if obj.batch_counter>=obj.num_batches || obj.failed_batch_counter>obj.max_failed_batches || obj.num_bad_batches>obj.max_num_bad_batches

                    if obj.use_analysis_dir_save
                        obj.saveResults;
                    end

                    if obj.use_analysis_dir_log
                        obj.saveSummary;
                    end

                end

            catch ME
                
                if exist(obj.log_name, 'file') && obj.use_analysis_dir_log
                    obj.write_log(obj.log_name, ME.getReport);
                end
                
                if obj.batch_counter>100 % if we managed to go through a big part of the run, might as well save the results
                
                    if obj.use_analysis_dir_save
                        obj.saveResults;
                    end

                    if obj.use_analysis_dir_log
                        obj.saveSummary;
                    end
                    
                end
                
                if obj.batch_counter<obj.num_batches && obj.failed_batch_counter<obj.max_failed_batches
                    rethrow(ME); % critical error in pipeline
                else
                    warning(ME.getReport); % minor error occured after pipeline is done
                end
                
            end
            
            if nargout>0
                obj_out = obj;
            end
            
        end
        
        function batch(obj)
            
            try
                obj.getData;
            catch ME
                warning(ME.getReport);
                obj.reader.advanceFile;
                return; % skip this batch, report it, and continue! 
            end
            
            obj.analysisStack;
            
            if obj.batch_counter==0 && (isempty(obj.use_astrometry) || obj.use_astrometry)
                try
                    obj.analysisAstrometry;
                catch ME
                    if obj.use_require_astrometry
                        rethrow(ME); % need to have astrometry to continue
                    else
                        warning(ME.getReport); % non essential to  successfully running the data analysis
                    end
                end
            end
           
            if ~isempty(obj.cat) % other tests??
                obj.magnitudes = obj.cat.magnitudes;
                obj.temperatures = obj.cat.temperatures;
                obj.coordinates = obj.cat.coordinates;
            end

            if obj.use_cutouts
               
                obj.analysisCutouts;
            
                if obj.use_photometry

                    obj.analysisPhotometry;

                    if obj.use_psf_model
                        obj.analysisModelPSF;
                    end

                    if obj.use_event_finding
                        try
                            obj.analysisEventFinding;
                            
                            if obj.use_max_bad_batches
                                
                                if obj.finder.checkBatchGood
                                    obj.num_bad_batches = 0;
                                else
                                    obj.num_bad_batches = obj.num_bad_batches + 1;
                                    if obj.num_bad_batches>obj.max_num_bad_batches
                                        if obj.debug_bit, fprintf('Found %d bad batches in a row. Quitting run...\n', obj.num_bad_batches); end
                                        obj.brake_bit = 1; % quit this run after the end of this batch 
                                    end
                                end
                                
                            end
                            
                        catch ME
%                             warning(ME.getReport); 
                            rethrow(ME); 
                        end
                    end
                    
                    if obj.use_cutouts_store
                        obj.analysisCutoutsStore; % maybe add a try-catch later? 
                    end

                end

            end
            
            if obj.use_analysis_dir_save && obj.use_save_batched_lightcurves && obj.use_photometry
                obj.phot.cutouts = [];
                obj.buf.input(obj.phot);
                obj.buf.directory = obj.log_dir; 
                obj.buf.save;
            end
            
            if obj.use_fits_save
                obj.analysisSaveFITS;
            end
            
            obj.updateHeader;
            
            obj.analysisDisplayGUI;
            
        end
        
        function getData(obj)
            
            %%%%%%%%%%%%%%%%%%%%% GET DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            t = tic;
            
            obj.clear;
            
            obj.reader.batch;
            obj.copyFrom(obj.reader); 
            
            obj.head.run_identifier = util.text.run_id(obj.thisFilename); 
            
            if isempty(obj.images) && isempty(obj.stack)
                disp(['empty batch in filename: ' obj.thisFilename]);
                return;
            elseif ~isempty(obj.images) && isempty(obj.stack) % we got images, need to produce the stack ourselves
                
                obj.num_sum = size(obj.images,3);
%                 obj.stack = util.stat.sum_single(obj.images); % sum along the 3rd dimension directly into single precision
                obj.stack = single(sum(obj.images,3));
%                 obj.positions = obj.clip.positions;
%                 obj.positions_bg = obj.clip_bg.positions;
            end
            
            obj.clip.positions = obj.positions;
            obj.clip.cut_size = size(obj.cutouts,1);

            if ~isempty(obj.positions_bg)
                obj.clip_bg.positions = obj.positions_bg;
                obj.clip_bg.cut_size = size(obj.cutouts_bg,1);
            end

            if isempty(obj.num_sum) 
                if ~isempty(obj.cutouts)
                    obj.num_sum = size(obj.cutouts,3);
                else
                    error('Unknown num_sum, and no cutouts to figure it out!');
                end
            end
               
            
%             if obj.use_cutout_adjustment
%                 
%                 if obj.use_cutout_adjustment_floor
%                     obj.clip.positions = floor(obj.positions) + obj.cutout_adjustment_pixels;
%                 else
%                     obj.clip.positions = obj.positions + obj.cutout_adjustment_pixels;
%                 end
%                 
%             end
            
            if isempty(obj.clip_bg.positions)
                obj.clip_bg.num_stars = 50;
                obj.clip_bg.cut_size = 20;
                obj.clip_bg.arbitraryPositions('imsize', size(obj.stack));
                obj.positions_bg = obj.clip_bg.positions;
            end
            
            if ~isempty(obj.num_stars) && ~isinf(obj.num_stars) && obj.num_stars<size(obj.cutouts,4) % use a smaller number of stars than we got
                
                if ~isempty(obj.positions)
                    obj.positions = obj.positions(1:obj.num_stars,:);
                    obj.clip.positions = obj.positions;
                end
                
                if ~isempty(obj.cutouts) 
                    obj.cutouts = obj.cutouts(:,:,:,1:obj.num_stars);
                end
                
                if ~isempty(obj.stack_cutouts)
                    obj.stack_cutouts = obj.stack_cutouts(:,:,:,1:obj.num_stars);
                end
                
            end
            
            if obj.debug_bit>1, fprintf('Time to get data: %f seconds\n', toc(t)); end
            
        end
        
        function analysisStack(obj)
            
            %%%%%%%%%%%%%%%%%%%%% STACK ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            t = tic;
            
            if isempty(obj.stack)
                return;
            end
            
            if isempty(obj.head.ROI)
                obj.cal.use_roi = 0;
            else
                obj.cal.use_roi = 1;
                obj.cal.ROI = obj.head.ROI; 
            end
            
            % calibrate the stack if needed
            if nnz(isnan(obj.stack)) % stack is already calibrated (has NaN values...)
                obj.stack_proc = obj.stack;
            else
                obj.stack_proc = obj.cal.input(obj.stack, 'sum', obj.num_sum);
            end
            
            if isempty(obj.positions) % if we didn't get star positions from file
                obj.findStars;
            end

            % cutouts of the stack
            obj.stack_cutouts = obj.clip.input(obj.stack_proc); 
            obj.stack_cutouts_bg = obj.clip_bg.input(obj.stack_proc); 
            
            % make background model based on stack cutouts
            if obj.use_background_cutouts || obj.use_background_stack
                obj.back.input(obj.stack_cutouts_bg, obj.clip_bg.positions);
            end
            
            % subtract background from stack (full frame)
            if obj.use_background_stack
                B = obj.back.getImage(size(obj.stack));
                obj.stack_proc = obj.stack_proc - B;
            end
            
            % subtract background from cutouts stack
            if obj.use_background_cutouts % if we already used background subtraction on all pixels, why not just cutout from stack_proc??
                BC = obj.back.getPoints(obj.clip.positions);
                BC = permute(BC, [4,3,2,1]); % turn the column vector into a 4D vector
                obj.stack_cutouts_sub = obj.stack_cutouts - BC; 
            else
                obj.stack_cutouts_sub = obj.stack_cutouts;
            end
            
            if obj.batch_counter==0
                obj.ref_stack = obj.stack_proc;
                obj.ref_positions = obj.positions;
            end
            
            obj.phot_stack.input(obj.stack_cutouts, 'positions', obj.positions, ...
                'timestamps', mean(obj.timestamps), 'filename', obj.reader.this_filename, ...
                't_start', obj.t_start, 't_end', obj.t_end, 't_end_stamp', obj.t_end_stamp, ...
                'juldates', mean(obj.juldates), 'variance', single(2.5)); % run photometry on the stack to verify flux and adjust positions

            if isempty(obj.cutouts) && ~isempty(obj.images) % only if we have images and not pre-cut cutouts
                obj.adjustPositions; % got a chance to realign the positions before making cutouts from raw images
            end
            
            if obj.use_check_flux % check if fluxes are still visible
                if obj.checkFluxes % if yes, reset the counter
                    obj.failed_batch_counter = 0;
                else % if not, check if this is happening many times in a row
                    obj.failed_batch_counter = obj.failed_batch_counter + 1;
                    if obj.failed_batch_counter>obj.max_failed_batches
                        if obj.debug_bit, fprintf('Cannot find stars %d times in a row. Quiting run...\n', obj.failed_batch_counter ); end
                        obj.brake_bit = 1; % finish this batch and then quit the run
                    end
                end
            end
            
            obj.flux_buf.input(obj.phot_stack.fluxes(:,:,end));% store the latest fluxes from the stack cutouts (to verify we didn't lose the stars)
            
            if obj.use_stack_all
                if obj.use_stack_all_proc
                    obj.stack_all = cat(3, obj.stack_all, obj.stack_proc);
                else
                    obj.stack_all = cat(3, obj.stack_all, obj.stack);
                end
            end
            
            if obj.debug_bit>1, fprintf('Time for stack analysis: %f seconds\n', toc(t)); end
            
        end
        
        function analysisAstrometry(obj)
            
            %%%%%%%%%%%%%%%%%%%%% ASTROMETRY ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%
            
            t = tic;
            
                ast = obj.use_astrometry;

                if isempty(obj.use_astrometry)  % automatically determine if we need to run astrometry
                    if isempty(obj.cat.data) && isempty(obj.cat.success) % no astrometry was attempted
                        ast = 1; % run astrometry 
                    else
                        ast = 0; % no need because either we have the data or the previous attempts have failed
                    end
                else
                    ast = obj.use_astrometry;
                end

                if ast 

                    obj.getDetectionParameters;
                    
%                     obj.cat.input(obj.stack_proc);
                    obj.cat.use_matched_only = 0; % cannot throw away stars with no match, we have cutouts already!! 
                    
                    if util.text.cs(obj.head.cam_name, 'balor')
                        obj.cat.input_rotation = -60; 
                    elseif util.text.cs(obj.head.cam_name, 'zyla')
                        obj.cat.input_rotation = -15; 
                    end
                    
                    obj.cat.inputPositions(obj.positions); 

                    if ~isempty(obj.cat.data) && obj.cat.success % successfully filled the catalog

                        obj.cat.num_stars = obj.num_stars;
%                         obj.cat.findStars(obj.positions); 

                        obj.positions = obj.cat.positions; % usually we will already have positions so this should do nothing (unless this analysis is on full frame rate images)

                        filename = fullfile(obj.reader.dir.pwd, 'catalog.mat');

                        if obj.use_save_astrometry
                            if isempty(obj.use_astrometry)
                                if ~exist(filename, 'file') % in auto-mode, only save if there was no catalog file
                                    obj.cat.saveMAT(filename);
                                end
                            elseif obj.use_astrometry % in force-astrometry mode must update the catalog file
                                obj.cat.saveMAT(filename);
                            end
                        end
                        
                        obj.head.THRESH_DETECTION = obj.cat.detection_threshold;
                        obj.head.LIMMAG_DETECTION = obj.cat.detection_limit; 

                        obj.magnitudes = obj.cat.magnitudes;
                        obj.coordinates = obj.cat.coordinates;
                        obj.temperatures = obj.cat.temperatures;
                       
                    elseif obj.cat.success==0
                        if obj.use_require_astrometry
                            error('Could not find an astrometric solution!'); 
                        else
                            fprintf('Warning: could not find an astrometric solution!\n'); 
                        end
                    end

                end

                        
            if obj.debug_bit>1, fprintf('Time for astrometry: %f seconds\n', toc(t)); end
            
        end
        
        function getDetectionParameters(obj)
            
            import util.text.cs;
            
            filename = fullfile(obj.reader.dir.pwd, 'A_README.txt'); 
            
            if ~exist(filename, 'file')
                return;
            end
            
            f = fopen(filename, 'r'); 
            
            on_cleanup = onCleanup( @()fclose(f)); 
            
            for ii = 1:1e6
                
                line = fgetl(f);
                
                if line==-1, break; end
                
                line = strip(line);
                
                if ~isempty(line)
                    
                    c = strsplit(line, ':'); 
                    
                    if length(c)>1 && ~isempty(c{2}) && ~cs(c{2}, '[]') 
                        
                        if cs(c{1}, 'detect_thresh', 'detection_threshold')
                            obj.cat.detection_threshold = str2double(c{2}); % consider a test if this field is already filled??
                        elseif cs(c{1}, 'NAXIS3', 6)
                            obj.cat.detection_stack_number = str2double(c{2}); % consider a test if this field is already filled?? 
                        elseif cs(c{1}, 'EXPTIME')
                            obj.cat.detection_exposure_time = str2double(c{2}); % consider a test if this field is already filled?? 
                        end
                    
                    end
                    
                end
                
            end
            
        end
        
        function analysisCutouts(obj)
        
            %%%%%%%%%%%%%%%%%%%%% CUTOUT ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%
            
            t = tic;
            
            if isempty(obj.cutouts)
                if ~isempty(obj.images)
                    obj.cutouts = obj.clip.input(obj.images);
                    obj.cutouts_bg = obj.clip_bg.input(obj.images);
                else
%                     error('Cannot produce cutouts without images!');
                    return;
                end
            end
            
            obj.cutouts_proc = obj.cal.input(obj.cutouts, 'clip', obj.clip);
            
            if ~isempty(obj.cutouts_bg)
                obj.cutouts_bg_proc = obj.cal.input(obj.cutouts_bg, 'clip', obj.clip_bg);
            end
            
            B = obj.back.getPoints(obj.clip.positions); 
            B = permute(B, [4,3,2,1]); % turn the column vector into a 4D vector
            % can also get variance from background object...
            
            if obj.use_refine_bg
                % use bg_cutouts to calculate overall differences between
                % frames to correct for regional results from background
                % object (based on the stack). 
            end
            
            if obj.use_background_cutouts
                obj.cutouts_sub = obj.cutouts_proc - B./obj.num_sum;
            else
                obj.cutouts_sub = obj.cutouts_proc;
            end
            
            if obj.use_cutouts_all
                
                if obj.use_cutouts_all_proc
                    obj.cutouts_all = cat(3, obj.cutouts_all, obj.cutouts_proc);
                else
                    obj.cutouts_all = cat(3, obj.cutouts_all, obj.cutouts);
                end
                
                obj.positions_x_all = cat(1, obj.positions_x_all, obj.positions(:,1)');
                obj.positions_y_all = cat(1, obj.positions_y_all, obj.positions(:,2)');
                
            end
            
            if obj.debug_bit>1, fprintf('Time for cutouts: %f seconds\n', toc(t)); end
            
        end
           
        function analysisPhotometry(obj)
        
            %%%%%%%%%%%%%%%%%%%%% PHOTOMETRY ANALYSIS %%%%%%%%%%%%%%%%%%%%%

            t = tic;

            obj.phot.input('images', obj.cutouts_sub, 'timestamps', obj.timestamps, 'filename', obj.reader.this_filename, ...
                't_start', obj.t_start, 't_end', obj.t_end, 't_end_stamp', obj.t_end_stamp, ...
                'juldates', obj.juldates, 'positions', obj.positions, 'variance', single(2.5)); % need to add the sky background too

            if obj.use_duplicate_filter
                if obj.findDuplicateFrames
                    obj.duplicate_batches = vertcat(obj.duplicate_batches, obj.batch_counter + 1); 
                end
            end
            
            if obj.use_full_lightcurves
                
                obj.lightcurves.getData(obj.phot);

                if obj.lightcurves.gui.check, obj.lightcurves.gui.update; end
                
            end
            
            if obj.use_stack_lightcurves
                
                obj.light_stack.getData(obj.phot_stack); 
                
                if obj.light_stack.gui.check, obj.light_stack.gui.update; obj.light_stack.gui.adjustAxes; end
                
            end
            
            if obj.debug_bit>1, fprintf('Time for photometry: %f seconds\n', toc(t)); end
            
            t = tic;
            
            obj.mean_buf.input(nanmean(obj.phot.fluxes));
            obj.var_buf.input(nanvar(obj.phot.fluxes));
            obj.back_buf.input(nanmean(obj.phot.backgrounds));
            obj.width_buf.input(nanmean(obj.phot.widths));
            
%             try
%                 obj.calcSkyParameters; % get an estimate of the zero point, seeing, background, limiting magnitude, etc. 
%             catch ME
%                 warning(ME.getReport); 
%             end
            
            if obj.debug_bit>1, fprintf('Time to calculate sky parameters: %f seconds\n', toc(t)); end
            
        end
        
        function [val, counter] = findDuplicateFrames(obj, star_index, aperture_index)
            
            if nargin<2 || isempty(star_index)
                star_index = 1;
            end
            
            if nargin<3 || isempty(aperture_index)
                aperture_index = 1;
            end
            
            f = obj.phot.fluxes(:,star_index,aperture_index); % pick any flux you like
            
            counter = 0;
            
            for ii = 1:size(f,1)-10
                
                if f(ii)==f(ii+10)
                    counter = counter + 1;
                end
                
            end
            
            if counter>3
                val = 1;
            else
                val = 0;
            end
            
        end
        
        function calcSkyParameters(obj) % take the photometery (and possible the catalog) and calculate seeing, background and zeropoint
            
            if ~isempty(obj.mean_buf) && ~isempty(obj.cat) && ~isempty(obj.cat.magnitudes) && ~isempty(obj.cat.success) && obj.cat.success==1

                S = nanmean(obj.phot.fluxes(:,:,end)-obj.phot.areas(:,:,end).*obj.phot.backgrounds(:,:,end))'; % signal (instrumental)
                N = nanstd(obj.phot.fluxes(:,:,end))'; % noise (instrumental)
                M = obj.cat.magnitudes; % magnitude from catalog

%                 obj.head.ZEROPOINT = util.stat.median2(S*10.^(0.4.*M'));
                instr_mag = -2.5*log10(S); 
                instr_mag(S<0) = NaN;
                instr_mag = real(instr_mag); 
                
                obj.head.ZEROPOINT = nanmedian(M-instr_mag); 
                
                idx = ~isnan(S) & ~isnan(N) & S./N>1.5; % choose stars that are actually measureable 
                S2 = S(idx);
                N2 = N(idx);
                M2 = M(idx);

                if isempty(S2)
                    return;
                end

%                 if isempty(obj.head.THRESH_INDIVIDUAL)
%                     obj.head.THRESH_INDIVIDUAL = obj.finder.min_star_snr;
%                 end
                
                if ~isempty(obj.aux_figure) && isvalid(obj.aux_figure)
                    delete(obj.aux_figure.Children);
                    ax = axes('Parent', obj.aux_figure);
                    obj.head.LIMMAG_INDIVIDUAL = head.limiting_magnitude(M2, S2./N2, obj.head.THRESH_INDIVIDUAL, 'maximum', 10, 'plot', 1, 'axes', ax, 'marker', 'pm'); 
                else
                    obj.head.LIMMAG_INDIVIDUAL = head.limiting_magnitude(M2, S2./N2, obj.head.THRESH_INDIVIDUAL, 'maximum', 10); 
                end
                
                % consider allowing the user to choose the type of photometry to take (the 3rd index)
                f = obj.phot_stack.fluxes(:,:,end)';
                
                b = nanmedian(squeeze(util.stat.median2(obj.stack_cutouts_bg))); % the background is attained by the median of the background cutouts
                obj.head.BACKGROUND = b./obj.head.NAXIS3; % the background is given per frame, not for the stack! 
                
                v = obj.phot_stack.variances(:,:,end)'; % the variance per star
                v(v==0) = NaN; % sometimes we get zero variance instead of NaN value
                
                a = obj.phot_stack.areas(:,:,end)'; % the area per star
                
                if ~obj.use_background_cutouts % if we didn't subtract background for the stack cutouts
                    f = f-a.*b; % correct the fluxes for background values.
                end
                
                S_stack = f; % signal (instrumental)
                N_stack = sqrt(a.*v); % noise (instrumental)
                
                idx = ~isnan(S_stack./N_stack) & S_stack./N_stack>1.5;
                S_stack = S_stack(idx);
                N_stack = N_stack(idx);
                M_stack = M(idx);
                
                if ~isempty(S_stack)

                    if isempty(obj.head.THRESH_STACK)
                        obj.head.THRESH_STACK = obj.head.THRESH_INDIVIDUAL;
                    end

                    if ~isempty(obj.aux_figure) && isvalid(obj.aux_figure)
    %                     delete(obj.aux_figure.Children);
    %                     ax = axes('Parent', obj.aux_figure);
                        ax.NextPlot = 'add';
                        obj.head.LIMMAG_STACK = head.limiting_magnitude(M_stack, S_stack./N_stack, obj.head.THRESH_STACK, 'maximum', Inf, 'var', 'snr', 'plot', 1, 'axes', ax); 
                        ax.NextPlot = 'replace';
                    else
                        obj.head.LIMMAG_STACK = head.limiting_magnitude(M_stack, S_stack./N_stack, obj.head.THRESH_STACK, 'maximum', Inf, 'var', 'snr'); 
                    end

                end
                
                drawnow;
                
%                 fprintf('the flux ratio of LIMMAG for individual and stack is %4.2f\n', 10.^(0.4*(obj.head.LIMMAG_STACK-obj.head.LIMMAG_INDIVIDUAL))); 
                
%                 T = table(-2.5*log10(S2)+obj.head.ZEROPOINT, S2./N2, M2, 'VariableNames', {'Measured_mag', 'SNR', 'GAIA_mag'});
%                 T2 = util.vec.bin_table_stats(T, 30);
%                 T2 = T2(T2.N>=5,:);
%                 fr = util.fit.polyfit(T2.SNR_nanmedian, T2.Measured_mag_nanmedian, 'order', 2);
%                 SNR = S2./N2;
%                 M3 = M2(SNR<10);
%                 SNR = SNR(SNR<10); 
%                 
%                 fr = util.fit.polyfit(SNR, M3, 'order', 2, 'plot', 1); 
%                 
%                 obj.head.MAG_LIMIT = fr.coeffs(1)+fr.coeffs(2).*thresh+fr.coeffs(3)*thresh.^2;
%                 
%                 if ~isempty(obj.aux_figure) && isvalid(obj.aux_figure)
%                     delete(obj.aux_figure.Children);
%                     ax = axes('Parent', obj.aux_figure);
%                     x_extrap = min(fr.x):-0.1:0;
%                     y_extrap = fr.coeffs(1)+fr.coeffs(2).*x_extrap + fr.coeffs(3).*x_extrap.^2;
%                     x_th = thresh;
%                     y_th = fr.coeffs(1)+fr.coeffs(2).*x_th + fr.coeffs(3).*x_th.^2;
% %                     plot(ax, T2.SNR_nanmedian, T2.Measured_mag_nanmedian, 'p', fr.x, fr.ym, 'r-', ...
% %                         x_extrap, y_extrap, 'r:', x_th, y_th, 'k+');
%                     plot(ax, SNR, M3, 'p', fr.x, fr.ym, 'r-', ...
%                         x_extrap, y_extrap, 'r:', x_th, y_th, 'k+');
%                     xlabel(ax, 'measured flux/rms, binned'); 
%                     ylabel(ax, 'measured magnitude, binned'); 
%                     drawnow;
%                 end

            end
            
        end
        
        function analysisModelPSF(obj)
            
            %%%%%%%%%%%%%%%%%%%%% PSF modeling %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            t = tic;

            obj.model_psf.input(obj.cutouts_sub, obj.phot.offsets_x(:,:,end), ...
                obj.phot.offsets_y(:,:,end), obj.phot.fluxes(:,:,end), obj.positions);

            obj.FWHM = obj.model_psf.fwhm;
            obj.FWHM_log = [obj.FWHM_log; obj.FWHM]; 
            
            T_mid = mean(obj.timestamps);
            N_frames = length(obj.timestamps); % number of frames per batch
            
            obj.juldates_full = [obj.juldates_full; obj.head.get_juldates(T_mid)];
            
            if obj.use_full_lightcurves
                obj.lightcurves.batch_timestamps(end+1,1) = T_mid; 
                obj.lightcurves.batch_midframe(end+1,1) = length(obj.lightcurves.timestamps) + N_frames/2; % add half the frame number to the number of existing frames
                obj.lightcurves.fwhm_coeffs_pix(end+1,:) = obj.model_psf.surf_coeffs;
                obj.lightcurves.fwhm_x_center = obj.model_psf.surf_fit.xc;
                obj.lightcurves.fwhm_y_center = obj.model_psf.surf_fit.yc;
            end
            
            if obj.use_stack_lightcurves
                obj.light_stack.batch_timestamps(end+1,1) = T_mid; 
                obj.light_stack.batch_midframe(end+1,1) = length(obj.light_stack.timestamps); % add half the frame number to the number of existing frames
                obj.light_stack.fwhm_coeffs_pix(end+1,:) = obj.model_psf.surf_coeffs';
                obj.light_stack.fwhm_x_center = obj.model_psf.surf_fit.xc;
                obj.light_stack.fwhm_y_center = obj.model_psf.surf_fit.yc;
            end
            
            if obj.debug_bit>1, fprintf('Time for PSF model: %f seconds\n', toc(t)); end
            
        end
           
        function analysisEventFinding(obj)
        
            %%%%%%%%%%%%%%%%%%%%% Event finding %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            t = tic;
            
            obj.head.PHOT_PARS = obj.phot.pars_struct; 

            obj.finder.input(obj.phot, obj.model_psf); 
            
            obj.finder.total_batches = obj.batch_counter; % keep track of all batches, not just those that were processed
            
            if obj.debug_bit>1, fprintf('Time to find events: %f seconds\n', toc(t)); end

            if ~isempty(obj.finder.gui) && obj.finder.gui.check
                obj.finder.gui.update;
            end
            
        end
        
        function analysisCutoutsStore(obj)
        
            obj.cutout_store.input(obj.phot); 
            
        end
        
        function analysisSaveFITS(obj)
            
            %%%%%%%%%%%%%%%%%%%% save FITS files of stacks %%%%%%%%%%%%%%%%
            
            [d, f] = fileparts(obj.thisFilename);

            d = strrep(d, ' (Weizmann Institute)', '');

            d = fullfile(d, 'FITS/');

            if ~exist(d, 'dir')
                mkdir(d);
            end

            fullname = fullfile(d,[f,'.fits']);
            fprintf('Saving "stack_proc" in FITS file: %s\n', fullname);
            
            I = double(obj.stack_proc);
            
            if obj.use_fits_roi && ~isempty(obj.fits_roi)
                I = I(obj.fits_roi(1):obj.fits_roi(1)+obj.fits_roi(3)-1,obj.fits_roi(2):obj.fits_roi(2)+obj.fits_roi(4)-1);
            end

            if obj.use_fits_flip
                I = rot90(I,2);
            end

            fitswrite(I, fullname); 
            obj.head.writeFITS(fullname, obj.timestamps(1), obj.num_sum);

            
        end
        
        function analysisDisplayGUI(obj)
            
            %%%%%%%%%%%%%%%%%%%% Update GUI and show stuff %%%%%%%%%%%%%%%%
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.update; 
            end

            if ~isempty(obj.finder.gui) && obj.finder.gui.check
                obj.finder.showLatest;
%             elseif ~isempty(obj.aux_figure) && isvalid(obj.aux_figure)
%                 obj.finder.showLatest(obj.aux_figure);
            end
            
            drawnow;
            
        end
        
        function findStars(obj)
            
            T = util.img.quick_find_stars(obj.stack_proc, 'psf', obj.getWidthEstimate, 'number', obj.num_stars,...
               'dilate', obj.cut_size-5, 'saturation', obj.saturation_value.*obj.num_sum, 'unflagged', 1); 
            
            if isempty(T)
                error('Could not find any stars using quick_find_stars!');
            end

            obj.clip.positions = T.pos;
            obj.positions = T.pos;

        end
        
        function adjustPositions(obj)
            
            if ~isempty(obj.phot_stack.gui) && obj.phot_stack.gui.check, obj.phot_stack.gui.update; end

            obj.checkRealign;

            obj.clip.positions = double(obj.clip.positions + obj.average_offsets);
            
            obj.model_psf.input(obj.stack_cutouts, obj.phot_stack.offsets_x, obj.phot_stack.offsets_y);
            
            obj.prev_average_width = obj.average_width;
            
        end
        
        function val = checkFluxes(obj)
            
            if size(obj.flux_buf.data, 2)~=size(obj.phot_stack.fluxes,2)
                obj.flux_buf.reset;
            end
            
            if is_empty(obj.flux_buf)
                val = 1;
            else
                
                mean_fluxes = obj.flux_buf.mean;
                mean_fluxes(mean_fluxes<=0) = NaN;

                new_fluxes = obj.phot_stack.fluxes(:,:,end);
                new_fluxes(isnan(mean_fluxes)) = [];
                mean_fluxes(isnan(mean_fluxes)) = [];

                flux_lost = sum(new_fluxes<0.5*mean_fluxes)>0.5*numel(mean_fluxes); % lost half the flux in more than half the stars...
                % want to add more tests...?

                val = ~flux_lost;
                
            end
            
        end
        
        function checkRealign(obj)
            
            if ~obj.checkFluxes
                
                disp('Lost star positions, using quick_align');

                [~,shift] = util.img.quick_align(obj.stack_proc, obj.ref_stack);
                obj.clip.positions = double(obj.ref_positions + flip(shift));

                obj.stack_cutouts = obj.clip.input(obj.stack_proc);

                obj.phot_stack.input(obj.stack_cutouts, 'positions', obj.positions, ...
                    'timestamps', mean(obj.timestamps), 'filename', obj.reader.this_filename, ...
                    't_start', obj.t_start, 't_end', obj.t_end, 't_end_stamp', obj.t_end_stamp, ...
                    'juldates', mean(obj.juldates), 'variance', single(2.5)); % run photometry on the stack to verify flux and adjust positions
                
                if ~isempty(obj.phot_stack.gui) && obj.phot_stack.gui.check, obj.phot_stack.gui.update; end

            end
            
        end
        
        function updateHeader(obj)
           
            % to be updated...
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.scan_vars(varargin{:});
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            I = obj.stack_proc;
            
            if obj.use_display_flip
                I = rot90(I,2);
            end
            
            util.plot.setImage(I, input.ax);

            obj.clip.showRectangles('color', 'black', 'ax', input.ax, 'delete', 1, 'text', 1, 'num', obj.display_num_rect_stars, 'flip', obj.use_display_flip);
            obj.clip_bg.showRectangles('color', 'red', 'ax', input.ax, 'delete', 0, 'text', 0, 'flip', obj.use_display_flip);
            
        end
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = img.gui.AnalysisGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
        function makeMonitor(obj)
            
            if isempty(obj.monitor)
                obj.monitor = img.gui.FutureMonitor(obj);
            end
            
            obj.monitor.show;
            
        end
        
    end    
    
end

