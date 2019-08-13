classdef Analysis < file.AstroData

    properties(Transient=true)
        
        gui;
        audio@util.sys.AudioControl;
        
        pool; 
        futures; 
        futures_dir;
        
        aux_figure;
        
    end
    
    properties % objects
        
        pars@head.Parameters;
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
        
%         light_original@img.Lightcurves;
%         light_basic@img.Lightcurves;
        lightcurves@img.Lightcurves;
%         light_gauss@img.Lightcurves;
%         light_cosqrt@img.Lightcurves;
        % light_fit@img.Lightcurves;
        
        model_psf@img.ModelPSF;
        
        finder@trig.Finder;
        
        sky_pars;
        
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
        
        FWHM; % latest measured full width half maximum
        
        batch_counter = 0;
        
    end
    
    properties % switches/controls
        
        num_stars = 500;
        cut_size = 21;
        saturation_value = 50000; % consider any pixels above this to be saturated
        
        use_background_stack = 1; % subtract b/g from the full-frame stack
        use_background_cutouts = 1; % subtract b/g from the cutouts (and stack cutouts!)
        use_refine_bg = 0; % need to figure out exactly how to do this
        
        use_check_flux = 1;
        max_failed_batches = 3; % if star flux is lost for more than this number of batches, quit the run
        
        use_astrometry = 0;
        use_cutouts = 1;
        use_photometry = 1;
        use_psf_model = 1;
        use_event_finding = 1;
        
        use_auto_load_cal = 1;
        
        use_fits_save = 0;
        use_fits_flip = 0;
        use_fits_roi = 0;
        fits_roi = [];
        
        use_audio = 0;
        
        use_display_flip = 0;
        display_num_rect_stars = 30;
        
        brake_bit = 1;
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        filename;
        directory;
        num_batches;
        average_width;
        average_offsets;
        
    end
    
    properties(Hidden=true)
       
        ref_stack; % for quick align
        ref_positions; 
        
        failed_batch_counter = 0;
        
        % these are used for backward compatibility with older versions of 
        % img.Clipper that made a slightly different cutout based on the 
        % same positions. This affects ONLY the cutting of CALIBRATION! 
        use_cutout_adjustment = 0; % turn adjustments on/off
        cutout_adjustment_pixels = 0; % how many pixels to push (back or forward) relative to today's positions
        use_cutout_adjustment_floor = 0; % use floor of positions before (instead of) using round(). 
        
        analysis_dir_save = 0;
        analysis_dir_log = 0;
        
        prev_average_width;
        
        num_batches_limit;
        version = 1.02;
        
    end
    
    methods % constructor
        
        function obj = Analysis(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'img.Analysis')
                if obj.debug_bit, fprintf('Analysis copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Analysis constructor v%4.2f\n', obj.version); end
            
                obj.reader = file.Reader;
                obj.cal = img.Calibration;
                obj.cal.load;
                obj.clip = img.Clipper;
                obj.clip.use_adjust = 0; % this should be disabled and depricated!
                obj.clip_bg = img.Clipper;
                obj.clip_bg.use_adjust = 0; % this should be disabled and depricated!
                obj.back = img.Background;
                obj.phot = img.Photometry;
                obj.phot.use_basic = 0;
                obj.phot.use_aperture = 1;
                obj.phot.use_gaussian = 0;
                
                obj.phot_stack = img.Photometry;
                obj.flux_buf = util.vec.CircularBuffer;
                obj.mean_buf = util.vec.CircularBuffer;
                obj.var_buf = util.vec.CircularBuffer;
                obj.back_buf = util.vec.CircularBuffer;
                obj.width_buf = util.vec.CircularBuffer;
                
%                 obj.light_original = img.Lightcurves; 
%                 obj.light_basic = img.Lightcurves; obj.light_basic.signal_method = 'square'; obj.light_basic.background_method = 'corners';
                obj.lightcurves = img.Lightcurves; obj.lightcurves.signal_method = 'aperture'; obj.lightcurves.background_method = 'annulus';
%                 obj.light_gauss = img.Lightcurves; obj.light_gauss.signal_method = 'gauss'; obj.light_gauss.background_method = 'annulus';
                
                obj.model_psf = img.ModelPSF;
                
                obj.finder = trig.Finder;
                obj.finder.loadFilterBank;
                
                obj.prog = util.sys.ProgressBar;
                obj.audio = util.sys.AudioControl;
                
                obj.pars = head.Parameters; % this also gives "pars" to all sub-objects
                obj.cat = head.Catalog(obj.pars);
                obj.finder.cat = obj.cat;
                
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
            
            obj.analysis_dir_save = 0;
            obj.analysis_dir_log = 0;
            
            obj.failed_batch_counter = 0;
            
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
            
            if isempty(obj.pars)
                val = obj.FWHM.*1.24;
            else
                val = obj.FWHM.*obj.pars.SCALE;
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
        
    end
    
    methods % setters
        
        function set.pars(obj,val)
            
            obj.pars = val;
            
            list = properties(obj);
            
            for ii = 1:length(list)
                
                try
                    if isobject(obj.(list{ii})) && ~isempty(obj.(list{ii})) && ~istable(obj.(list{ii})) && isprop(obj.(list{ii}), 'pars') 
                        obj.(list{ii}).pars = val;
                    end
                catch ME
                    
                    disp(['trouble setting "pars" into variable: ' list{ii}]);
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
            
            obs_date = obj.pars.STARTTIME;
            read_date = util.text.time2str(datetime('now', 'TimeZone', 'UTC'));
            
            if isempty(str)
                ev_str = '';

                for ii = 1:length(obj.finder.last_events)

                    if obj.finder.last_events(ii).keep
                        star_str = '*';
                    else
                        star_str = '';
                    end

                    ev_str = sprintf('%s%4.2f%s ', ev_str, obj.finder.last_events(ii).snr, star_str);

                end

                f = [0 0 0];
                if ~isempty(obj.fluxes), f(1) = obj.fluxes(1,1); end
                if size(obj.fluxes,2)>=10, f(2) = obj.fluxes(1,10); end
                if size(obj.fluxes,2)>=100, f(3) = obj.fluxes(1,100); end

                fprintf(fid, 'Batch: %04d, ObsDate: %s, Flux: [% 9.1f % 8.1f % 7.1f]', obj.batch_counter+1, obs_date, f(1), f(2), f(3));
                
                fprintf(fid, ' | seeing: %4.2f" | back: %5.3f | zp: %6.4g | noise: %5.2f | lim. mag: %4.2f', ...
                    obj.sky_pars.seeing, obj.sky_pars.background, obj.sky_pars.zero_point, obj.sky_pars.noise_level, obj.sky_pars.limiting_mag);
                
                fprintf(fid, ' | Events S/N: [%s], ReadDate: %s\n', ev_str, read_date);
            
            else
                fprintf(fid, '%s: %s\n', read_date, str);
            end
                
        end
        
        function saveResults(obj, dirname, varargin)
            
            % add parsing of varargin later
            
            
            name = obj.pars.OBJECT;
            
            obj.lightcurves.saveAsMAT(fullfile(dirname, ['lightcurves_' name]));
            
            finder = obj.finder;
            save(fullfile(dirname, ['finder_' name]), 'finder', '-v7.3');
            
        end
        
        function saveSummary(obj, dirname)
            
            filename = ['summary_' obj.pars.OBJECT '.txt'];
            
            fid = fopen(fullfile(dirname, filename), 'wt');
            
            if fid<0
                warning('Cannot open file %s', fullfile(dirname, filename));
            else

                onc = onCleanup(@() fclose(fid));

                fprintf(fid, 'Summary for run %s, with %d batches.\n', obj.pars.OBJECT, obj.batch_counter);

                v = abs(obj.finder.snr_values);

                fprintf(fid, 'S/N for the last %d batches is distrubuted: min= %f median= %f max= %f\n',...
                    numel(v), min(v, [], 'omitnan'), median(v, 'omitnan'), max(v, [], 'omitnan'));

                v = abs([obj.finder.all_events.snr]);

                fprintf(fid, 'S/N for %d triggered events is distrubuted: min= %f median= %f max= %f\n',...
                    numel(v), min(v, [], 'omitnan'), median(v, 'omitnan'), max(v, [], 'omitnan'));

                v = abs([obj.finder.kept_events.snr]);

                fprintf(fid, 'S/N for %d kept events is distrubuted: min= %f median= %f max= %f\n',...
                    numel(v), min(v, [], 'omitnan'), median(v, 'omitnan'), max(v, [], 'omitnan'));

                fprintf(fid, 'Number of events: total= %d | kept= %d\n', length(obj.finder.all_events), length(obj.finder.kept_events));
            
                fprintf(fid, 'Star hours above stellar S/N of %4.2f: total= %4.2f | lost= %4.2f | kept= %4.2f\n', obj.finder.min_star_snr, ...
                    obj.finder.star_hours_total, obj.finder.star_hours_lost, obj.finder.star_hours_total-obj.finder.star_hours_lost);
                
                fprintf(fid, 'Star hours above stellar S/N of %4.2f: total= %4.2f | lost= %4.2f | kept= %4.2f\n', obj.finder.min_star_snr*2, ...
                    obj.finder.star_hours_total_better, obj.finder.star_hours_lost_better, obj.finder.star_hours_total_better-obj.finder.star_hours_lost_better);
                
                fprintf(fid, 'Star hours above stellar S/N of %4.2f: total= %4.2f | lost= %4.2f | kept= %4.2f\n', obj.finder.min_star_snr*4, ...
                    obj.finder.star_hours_total_best, obj.finder.star_hours_lost_best, obj.finder.star_hours_total_best-obj.finder.star_hours_lost_best);
                
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
            
            if obj.debug_bit, disp(['Finished run with ' num2str(obj.batch_counter) ' batches.']); end
            
        end
        
        function async_run(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('reset', 0); % reset the object before running (i.e., start a new run)
            input.input_var('logging', []); % create log files in the analysis folder
            input.input_var('save', []); % save the events and lightcurves from this run
            input.input_var('worker', []); % index of worker/future you want to use
            input.scan_vars(varargin{:});
            
            if isempty(input.worker)
                input.worker = obj.findWorker;
            end
            
            obj.futures{input.worker} = parfeval(obj.pool, @obj.run, 1, 'reset', input.reset, 'logging', input.logging, 'save', input.save); 
            obj.futures_dir{input.worker} = obj.reader.dir.two_tail;
            
        end
        
        function idx = findWorker(obj, varargin)
            
            if isempty(gcp('nocreate'))
                obj.pool = parpool;
                obj.pool.IdleTimeout = 360;
            end
            
            N = obj.pool.NumWorkers; 
            
            idx = [];
            
            for ii = 1:N
                
                if length(obj.futures)<ii || ~isvalid(obj.futures{ii}) || strcmp(obj.futures{ii}.State, 'finished')
                    idx = ii;
                    return;
                end
                
            end
            
            if ii==N
                error('Cannot find a free worker to run analysis...');
            end
            
            
        end
        
        function obj = run(obj, varargin)
            
            try
            
            input = util.text.InputVars;
            input.input_var('reset', 0); % reset the object before running (i.e., start a new run)
            input.input_var('logging', []); % create log files in the analysis folder
            input.input_var('save', []); % save the events and lightcurves from this run
            input.scan_vars(varargin{:});
            
            if input.reset
                obj.reset;
            end
            
            % update hidden variables in case we use GUI to stop then continue this run
            if ~isempty(input.logging), obj.analysis_dir_log = input.logging; end            
            if ~isempty(input.save), obj.analysis_dir_save = input.save; end
            
            if obj.analysis_dir_log || obj.analysis_dir_save
                
                log_time = datetime('now', 'TimeZone', 'UTC');
                log_dir = fullfile(obj.reader.dir.pwd, ['analysis_' char(log_time, 'yyyy-MM-dd')]);
                log_name = fullfile(log_dir, 'analysis_log.txt');
                log_obj = fullfile(log_dir, 'analysis_parameters.txt');
                
                if obj.batch_counter==0
                    if exist(log_dir, 'dir')
                        error('Folder %s already exists! Is this run already being processed?', log_dir);
                    else
                        mkdir(log_dir);
                    end
                end
                
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
                
                if ~isempty(date)
                    obj.cal.loadByDate(datestr(date, 'yyyy-mm-dd'), 0); % last argument is to NOT reload if date is consistent
                end
                
            end
            
            if ~obj.cal.checkDark
                error('Cannot start a new run without loading darks into calibration object!');
            end
            
            cleanup = onCleanup(@() obj.finishup);
            obj.startup;
            
            for ii = obj.batch_counter+1:obj.num_batches
                
                if obj.brake_bit
                    break;
                end
                
                obj.batch;
                
                if obj.analysis_dir_log
                    obj.write_log(log_name);
                    if ~exist(log_obj, 'file')
                        util.oop.save(obj, log_obj, 'hidden', 1); 
                    end
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
            if obj.batch_counter>=obj.num_batches || obj.failed_batch_counter>obj.max_failed_batches
                
                if obj.analysis_dir_save
                    obj.saveResults(log_dir);
                end

                if obj.analysis_dir_log
                    obj.saveSummary(log_dir);
                end
                
            end
            
            catch ME
                
                if obj.analysis_dir_log
                    obj.write_log(log_name, ME.getReport);
                end
                
                rethrow(ME);
            end
            
        end
        
        function batch(obj)
            
            obj.getData;
            
            obj.analysisStack;
            
            if obj.use_cutouts
               
                obj.analysisCutouts;
            
                if obj.use_photometry

                    obj.analysisPhotometry;

                    if obj.use_psf_model
                        obj.analysisModelPSF;
                    end

                    if obj.use_event_finding
                        obj.analysisEventFinding;
                    end

                end

            end
            
            if obj.use_fits_save
                obj.analysisSaveFITS;
            end
            
            obj.analysisDisplayGUI;
            
        end
        
        function getData(obj)
            
            %%%%%%%%%%%%%%%%%%%%% GET DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            t = tic;
            
            obj.clear;
            
            obj.reader.batch;
            obj.copyFrom(obj.reader); 
            
            if isempty(obj.images) && isempty(obj.stack)
                disp(['empty batch in filename: ' obj.thisFilename]);
                return;
            elseif ~isempty(obj.images) % we got images, need to produce cutouts, stack and positions outselves
                
                obj.num_sum = size(obj.images,3);
                obj.stack = util.stat.sum_single(obj.images); % sum along the 3rd dimension directly into single precision
                obj.positions = obj.clip.positions;
                obj.positions_bg = obj.clip_bg.positions;
                
            elseif ~isempty(obj.stack) % got stack (and assume we got cutouts and positions, too)
                
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
                obj.clip_bg.arbitraryPositions;
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
            
            % calibrate the stack if needed
            if nnz(isnan(obj.stack)) % stack is already calibrated (has NaN values...)
                obj.stack_proc = obj.stack;
            else
                obj.stack_proc = obj.cal.input(obj.stack, 'sum', obj.num_sum);
            end
            
            if isempty(obj.use_astrometry) || obj.use_astrometry
                try
                    obj.analysisAstrometry;
                catch ME
                    warning(ME.getReport); % non essential to  successfully running the data analysis
                end
            end
           
            if ~isempty(obj.cat) % other tests??
                obj.magnitudes = obj.cat.magnitudes;
                obj.temperatures = obj.cat.temperatures;
                obj.coordinates = obj.cat.coordinates;
            end

            if isempty(obj.positions) % if astrometry failed or was not used, we need to get positions somehow
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
            
            obj.phot_stack.input(obj.stack_cutouts_sub, 'positions', obj.clip.positions); % run photometry on the stack to verify flux and adjust positions
            
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
            
            obj.flux_buf.input(obj.phot_stack.fluxes);% store the latest fluxes from the stack cutouts (to verify we didn't lose the stars)
            
            if obj.debug_bit>1, fprintf('Time for stack analysis: %f seconds\n', toc(t)); end
            
        end
        
        function analysisAstrometry(obj)
            
            %%%%%%%%%%%%%%%%%%%%% ASTROMETRY ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%
            
            t = tic;
            
            if obj.batch_counter==0
            
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

                    obj.cat.input(obj.stack_proc);

                    if ~isempty(obj.cat.data) % successfully filled the catalog

                        obj.cat.num_stars = obj.num_stars;
                        obj.cat.findStars(obj.positions); 

                        obj.positions = obj.cat.positions; % usually we will already have positions so this should do nothing (unless this analysis is on full frame rate images)

                    end

                    filename = fullfile(obj.reader.dir.pwd, 'catalog.mat');

                    if ~isempty(obj.cat.data)
                        if isempty(obj.use_astrometry)
                            if ~exist(filename, 'file') % in auto-mode, only save if there was no catalog file
                                obj.cat.saveMAT(filename);
                            end
                        elseif obj.use_astrometry % in force-astrometry mode must update the catalog file
                            obj.cat.saveMAT(filename);
                        end
                    end

                end

            end
            
            if ~isempty(obj.cat.data)
                obj.magnitudes = obj.cat.magnitudes;
                obj.coordinates = obj.cat.coordinates;
                obj.temperatures = obj.cat.temperatures; 
            end
            
            if obj.debug_bit>1, fprintf('Time for astrometry: %f seconds\n', toc(t)); end
            
        end
        
        function analysisCutouts(obj)
        
            %%%%%%%%%%%%%%%%%%%%% CUTOUT ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%
            
            t = tic;
            
            if isempty(obj.cutouts)
                if ~isempty(obj.images)
                    obj.cutouts = obj.clip.input(obj.images);
                    obj.cutouts_bg = obj.clip_bg.input(obj.images);
                else
                    error('Cannot produce cutouts without images!');
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
            
            
            if obj.debug_bit>1, fprintf('Time for cutouts: %f seconds\n', toc(t)); end
            
        end
           
        function analysisPhotometry(obj)
        
            %%%%%%%%%%%%%%%%%%%%% PHOTOMETRY ANALYSIS %%%%%%%%%%%%%%%%%%%%%

            t = tic;

            obj.phot.input('images', obj.cutouts_sub, 'timestamps', obj.timestamps, ...
                'positions', obj.positions, 'variance', single(2.5)); % need to add the sky background too

            obj.lightcurves.getData(obj.phot);
            if obj.lightcurves.gui.check, obj.lightcurves.gui.update; end

            if obj.debug_bit>1, fprintf('Time for photometry: %f seconds\n', toc(t)); end
            
            t = tic;
            
            obj.mean_buf.input(mean(obj.phot.fluxes,1,'omitnan'));
            obj.var_buf.input(var(obj.phot.fluxes,[],1,'omitnan'));
            obj.back_buf.input(mean(obj.phot.backgrounds,1,'omitnan'));
            obj.width_buf.input(mean(obj.phot.widths,1,'omitnan'));
            
            obj.calcSkyParameters;
            
            if obj.debug_bit>1, fprintf('Time to calculate sky parameters: %f seconds\n', toc(t)); end
            
        end
        
        function calcSkyParameters(obj) % take the stack photometery (and possible the catalog) and calculate seeing, background and zeropoint
            
            import util.stat.median2;
            
            if isempty(obj.mean_buf) % any other tests??
                obj.sky_pars = [];
            else
                
                obj.sky_pars = struct;
%                 pars.seeing = median(obj.width_buf.median,2,'omitnan').*obj.pars.SCALE.*2.355;
                obj.sky_pars.seeing = median2(obj.phot.widths).*obj.pars.SCALE.*2.355;
%                 pars.background = median(obj.back_buf.median, 2,'omitnan');
                obj.sky_pars.background = median2(obj.phot.backgrounds);
                
                if ~isempty(obj.cat) && ~isempty(obj.cat.magnitudes) % this is a fairly good indicator that mextractor/astrometry worked
                    
%                     S = double(obj.mean_buf.median)';
%                     N = double(sqrt(obj.var_buf.median))';

                    S = double(mean(obj.phot.fluxes, 1, 'omitnan'))';
                    N = double(std(obj.phot.fluxes, [], 1, 'omitnan'))';
                    
%                     pars.zero_point = median(obj.mean_buf.median.*10.^(0.4.*obj.cat.magnitudes') ,2,'omitnan');
                    obj.sky_pars.zero_point = median2(S*10.^(0.4.*obj.cat.magnitudes'));
                    
                    idx = ~isnan(S) & ~isnan(N) & S./N>3 & S./N<8; 
                    S2 = S(idx);
                    N2 = N(idx);

                    if ~isempty(S)
                        
                        thresh = 5; % minimal S/N 
                        
                        fr = util.fit.polyfit(S2, N2.^2, 'order', 2, 'sigma', 2, 'iterations', 10); 
                        obj.sky_pars.noise_level = sqrt(fr.coeffs(1));
                        
                        % solve the equation N(S) = c1 + c2*S + c3*S^2 == S/thresh
                        a = fr.coeffs(3)-1./thresh.^2;
                        b = fr.coeffs(2);
                        c = fr.coeffs(1);
                        
                        s1 = (-b-sqrt(b.^2-4*a*c))./(2*a);
                        s2 = (-b+sqrt(b.^2-4*a*c))./(2*a);
                        
                        obj.sky_pars.limiting_signal = min(s1,s2);
                        if obj.sky_pars.limiting_signal<0
                            obj.sky_pars.limiting_signal = max(s1,s2);
                        end
                        
                        
                        if ~isreal(obj.sky_pars.limiting_signal) || isnan(obj.sky_pars.limiting_signal) || obj.sky_pars.limiting_signal<0
                            obj.sky_pars.limiting_mag = 2.5*log10(obj.sky_pars.zero_point./(obj.sky_pars.noise_level*thresh));
                        else
                            obj.sky_pars.limiting_mag = 2.5*log10(obj.sky_pars.zero_point./obj.sky_pars.limiting_signal);
                        end
                    end

                end
                
            end
            
        end
        
        function analysisModelPSF(obj)
            
            %%%%%%%%%%%%%%%%%%%%% PSF modeling %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            t = tic;

            obj.model_psf.input(obj.cutouts_sub, obj.phot.offsets_x, obj.phot.offsets_y);

            obj.FWHM = util.img.fwhm(obj.model_psf.stack);

            if obj.debug_bit>1, fprintf('Time for PSF model: %f seconds\n', toc(t)); end
            
        end
           
        function analysisEventFinding(obj)
        
            %%%%%%%%%%%%%%%%%%%%% Event finding %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            t = tic;
            
            f = obj.phot.fluxes;
            e = obj.phot.errors;
            a = obj.phot.areas;
            b = obj.phot.backgrounds;
            v = obj.phot.variances;
            x = obj.phot.offsets_x;
            y = obj.phot.offsets_y;
            w = obj.phot.widths;
            p = obj.phot.bad_pixels;
            phot_pars = obj.phot.pars_struct; % maybe also give this to model_psf??

            r = [];
            g = [];

            if obj.phot.use_gaussian
                g = obj.phot.gauss_sigma;
            elseif obj.phot.use_aperture

                r = obj.phot.aperture;
            end

            obj.finder.input(f, e, a, b, v, x, y, w, p, r, g, ...
                obj.timestamps, obj.cutouts_proc, obj.positions, obj.stack_proc, ...
                obj.batch_counter+1, 'filename', obj.thisFilename, ...
                't_end', obj.t_end, 't_end_stamp', obj.t_end_stamp,...
                'used_background', obj.phot.use_backgrounds, 'pars', phot_pars);

            if obj.debug_bit>1, fprintf('Time to find events: %f seconds\n', toc(t)); end

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
            if obj.use_fits_flip
                I = rot90(I,2);
            end

            if obj.use_fits_roi && ~isempty(obj.fits_roi)
                I = I(obj.fits_roi(1):obj.fits_roi(1)+obj.fits_roi(3)-1,obj.fits_roi(2):obj.fits_roi(2)+obj.fits_roi(4)-1);
            end

            fitswrite(I, fullname); 
            obj.pars.writeFITS(fullname, [], obj.num_sum);

            
        end
        
        function analysisDisplayGUI(obj)
            
            %%%%%%%%%%%%%%%%%%%% Update GUI and show stuff %%%%%%%%%%%%%%%%
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.update; 
            end

            if ~isempty(obj.finder.gui) && obj.finder.gui.check
                obj.finder.showLatest(obj.finder.gui.panel_image);
            elseif ~isempty(obj.aux_figure) && isvalid(obj.aux_figure)
                obj.finder.showLatest(obj.aux_figure);
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

                new_fluxes = obj.phot_stack.fluxes;
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

                obj.phot_stack.input(obj.stack_cutouts, 'positions', obj.positions); % run photometry on the stack to verify flux and adjust positions
                if ~isempty(obj.phot_stack.gui) && obj.phot_stack.gui.check, obj.phot_stack.gui.update; end

            end
            
        end
        
        function S = runMextractor(obj, I)
            
            if isempty(which('mextractor'))
                error('Cannot load the MAAT package. Make sure it is on the path...');
            end
            
            if nargin<2 || isempty(I)
                I = obj.stack_proc;
            end
            
            if isempty(I)
                error('Must supply an image to run mextractor (or fill stack_proc).');
            end
            
            I = regionfill(I, isnan(I));
            
            S = SIM;
            S.Im = I;
            evalc('S = mextractor(S);');
            
            SN = S.Cat(:,find(strcmp(S.ColCell, 'SN')));
            SN2 = S.Cat(:,find(strcmp(S.ColCell, 'SN_UNF')));
            S.Cat = S.Cat(SN>SN2-2,:);
            
            obj.image_mextractor = S;
            
        end
        
        function SS = runAstrometry(obj, S)
            
            if isempty(which('astrometry'))
                error('Cannot load the MAAT package. Make sure it is on the path...');
            end
            
            if nargin<2 || isempty(S)
                if ~isempty(obj.image_mextractor)
                    S = obj.image_mextractor;
                else
                    error('Must supply a SIM object to run astrometry (or fill image_mextractor).');
                end
            end
            
%             addpath(fullfile(getenv('DATA'), 'GAIA\DR2'));
            
            [~,S]=astrometry(S, 'RA', obj.pars.RA_DEG/180*pi, 'Dec', obj.pars.DEC_DEG/180*pi, 'Scale', obj.pars.SCALE,...
                'Flip',[1 1;1 -1;-1 1;-1 -1], 'RefCatMagRange', [7 17], 'BlockSize', [3000 3000], 'ApplyPM', false, ...
                'MinRot', -25, 'MaxRot', 25);
            
            % update RA/Dec in catalog according to WCS
            obj.image_mextractor = update_coordinates(S);
            
            %  Match sources with GAIA
            SS = catsHTM.sources_match('GAIADR2',obj.image_mextractor);
            
            obj.matched_gaia = SS;
            
            
            
        end
        
        function T = makeCatalog(obj)
            
            S = obj.image_mextractor;
            SS = obj.matched_gaia;
            
            T = array2table([SS.Cat, S.Cat], 'VariableNames', [SS.ColCell, S.ColCell]);
            
            T = T(~isnan(T{:,1}),:);
             
            [~, idx] = unique(T{:,1:2}, 'rows');
            T = T(idx,:);

%             T.Properties.VariableNames; % change variable names??

            T.RA = T.RA.*180/pi;
            T.Dec = T.Dec.*180/pi;
            T.Dist = T.Dist.*180/pi*3600;
            
            T.ALPHAWIN_J2000 = T.ALPHAWIN_J2000.*180/pi;
            T.DELTAWIN_J2000 = T.DELTAWIN_J2000.*180/pi;
            
            T.Properties.VariableUnits = {'deg', 'deg', 'year', '"', '"', '"', '"', '"', '"', '"', '"', '', '', '', ...
                'mag', 'mag', 'mag', 'mag', 'mag', 'mag', 'km/s', 'km/s', '', 'K', 'K', 'K', '', '', '"', '',  ...
                'pix', 'pix', 'pix', 'pix', 'pix', 'pix', 'pix', 'deg', '', ...
                'deg', 'deg', 'counts', 'counts', '', '', '', '','', 'counts', 'counts', 'mag', 'mag', '', '', '', ...
                'counts', 'counts', 'counts', 'counts', 'counts', 'counts', 'counts', ...
                'counts', 'counts', 'counts', 'counts', 'counts', 'counts', 'counts', ...
                '', '', 'arcsec'}; % input units for all variables
            
            obj.catalog = T;

        end
        
        function findStarsMAAT(obj) % to be depricated! 
            
            if isempty(which('mextractor'))
                error('Cannot load the MAAT package. Make sure it is on the path...');
            end
             
            % add additional tests to remove irrelvant stars
            
            if obj.min_star_temp
                T = T(T{:,'Teff'}>=obj.min_star_temp,:); % select only stars with temperature above minimal level (hotter stars have smaller angular scale)
            end
            
            T = sortrows(T, 'Mag_G'); % sort stars from brightest to faintest
            
            obj.positions = [T.XPEAK_IMAGE T.YPEAK_IMAGE];
            obj.clip.positions = obj.positions;
            
            obj.magnitudes = T{:,'Mag_G'};
            obj.coordinates = [T.RA T.Dec];
            
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
        
    end    
    
end

