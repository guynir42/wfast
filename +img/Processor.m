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
        lightcurves@img.Lightcurves;
        
        flux_buf@util.vec.CircularBuffer;
        
        pars@util.text.InputVars; % store all controls in this
        last_pars@util.text.InputVars; % parameters used in latest run
        
        func; 
        
        % ... add streak detection maybe?
        
    end
    
    properties % inputs/outputs
        
        data; % struct with all input/output data 
        futures = {}; % used for running asynchronuous tasks
        
    end
    
    properties % switches/controls
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        brake_bit = 1; % when running, this is 0, when it is turned to 1 (e.g., by GUI) the run is paused
        
        current_batch = 1;
        
        backup_pars@util.text.InputVars;
        
        failed_batch_counter = 0;
        
        successful_batches = 0;
        
        output_folder_name = ''; 
        
        version = 1.00;
        
    end
    
    methods % constructor
        
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
                
                obj.lightcurves = img.Lightcurves;
                obj.lightcurves.head = obj.head;
                obj.lightcurves.cat = obj.cat;

                obj.flux_buf = util.vec.CircularBuffer;

                obj.prog = util.sys.ProgressBar;
                
                obj.makeParsObject;
                obj.reset;
                
            end
            
        end
             
        function makeParsObject(obj)
            
            obj.pars = util.text.InputVars;
            
            obj.pars.input_var('reset', true); obj.pars.add_comment('reset the object data at the beginning of the run'); 
            
            obj.pars.input_var('num_batches', [], 6); obj.pars.add_comment('maximum number of batches to run. Run may end before that if files run out or stars are lost'); 
            obj.pars.input_var('cut_size', 15); obj.pars.add_comment('number of pixels in a cutout'); 
            obj.pars.input_var('num_stars', 2000, 6); obj.pars.add_comment('maximum number of stars to extract from the images');
            obj.pars.input_var('threshold', 5); obj.pars.add_comment('how many "sigma" above the noise'); 
            obj.pars.input_var('saturation', 5e4); obj.pars.add_comment('for a single image, if stack is read, this is multiplied by num_sum'); 
            obj.pars.input_var('use_bg_map', false, 6); obj.pars.add_comment('if true, will adjust each point in the stack by the interpolated b/g values, instead of the median value'); 
            
            obj.pars.input_var('lock_adjust', 'all'); 
            obj.pars.add_comment(sprintf('which cutouts should be moved together (locked), which can move separately (unlock).\n Choose "all", "none", or "stars" that only locks targets with GAIA matches')); 
            
            obj.pars.input_var('use_astrometry', true, 6); obj.pars.add_comment('make sure to solve the astrometry before continuing'); 
            obj.pars.input_var('use_remove_bad_matches', true, 6); obj.pars.add_comment('remove stars that have no GAIA match if their S/N is too low'); 
            obj.pars.input_var('bad_match_min_snr', 30); obj.pars.add_comment('if an object has no GAIA match, it can be removed if its S/N is below this threshold'); 
            
            obj.pars.input_var('use_auto_load_cal', true, 6); obj.pars.add_comment('load calibration based on the date of the files'); 
            
            obj.pars.input_var('num_stars_filter_kernel', 200, 6); obj.pars.add_comment('use a rough star finding on the brightest stars, just to estimate the PSF width'); 
            obj.pars.input_var('initial_guess_psf_width', 1.0); obj.pars.add_comment('typical value just for finding those stars on the 1st attempt'); 
            
            obj.pars.input_var('use_check_flux', true, 6); obj.pars.add_comment('check if the star flux is gone for a few batches and cut the run short'); 
            obj.pars.input_var('max_failed_batches', 3, 6); obj.pars.add_comment('if star flux is lost for more than this number of batches, quit the run'); 
            
            obj.pars.input_var('display_num_rect_stars', 30, 10); obj.pars.add_comment('how many rectangles to display'); 
            obj.pars.input_var('use_display_rect_text', false, 6); obj.pars.add_comment('show the clip number for each square'); 
            
            obj.pars.input_var('use_save_results', false, 6); obj.pars.add_comment('save the lightcurves etc. to file at end of run'); 
            obj.pars.input_var('use_overwrite', false, 6); obj.pars.add_comment('if true, will happily clear existing "processing" folders with the same date'); 
            obj.pars.input_var('output_folder', ''); obj.pars.add_comment('folder inside original data, typically "processor_YYYY-MM-DD" for storing the results and log files'); 
            
            obj.pars.input_var('use_fits_save', false, 10); obj.pars.add_comment('save FITS files in a folder below the raw images'); 
            obj.pars.input_var('use_fits_roi',false, 10); obj.pars.add_comment('use a region of interest (ROI) to cut before saving the FITS'); 
            obj.pars.input_var('fits_roi', []); obj.pars.add_comment('define the region of interest (ROI) in [left, top, width, height] in the matlab image'); 
            obj.pars.input_var('use_fits_flip', false, 10); obj.pars.add_comment('flip the image by 180 degrees before saving to FITS (this is done after cutting out the FITS_ROI'); 
            
            obj.backup_pars = obj.pars;
            
        end
           
    end
    
    methods % reset/clear
        
        function restorePars(obj)
            
            obj.pars = obj.backup_pars; 
            
        end
        
        function reset(obj) % call this at the start of a new run
            
            % objects
            obj.reader.reset;
            obj.cat.reset;
            obj.phot.reset;
            obj.lightcurves.reset;
            obj.flux_buf.reset;
            
            % info about the run
            obj.data.date = ''; % date when images were taken (date of beginning of night) in YYYY-MM-DD format
            obj.data.camera = ''; % name of camera (can be e.g., Zyla or Balor)
            obj.data.project = ''; % name of telescope/project (can be e.g., Kraar, LAST or WFAST)
            
            % related to find stars 
            obj.data.bg_mean = []; % background mean using im_stats (scalar or map)
            obj.data.bg_var = []; % background variance using im_stats (scalar or map)
            obj.data.psf_width = []; % calculate this from a sample of bright stars then use it in quick_find_stars
            obj.data.positions = []; % positions are found from the stack not from what is given in the file
            obj.data.found_stars = []; % a table with the results from quick_find_stars
            
            % references for quick_align 
            obj.data.ref_image = []; % save a copy of the first image in the run (used to quick_align)
            obj.data.ref_positions = []; % save a copy of the first positions in the run (used to quick_align)
            
            % counters
            obj.current_batch = 1;
            obj.failed_batch_counter = 0;
            obj.successful_batches = 0;
            
        end
        
        function clear(obj) % call this each batch
            
            obj.phot.clear;
            obj.cat.clear;
            obj.reader.clear;
            
            obj.data.timestamps = []; % the time of the beginning of the acquisition of the image/stack 
            obj.data.juldates = []; % matched to global time
            obj.data.image = []; % raw images or stack from file
            obj.data.image_proc = []; % after calibration, background subtraction, etc...
            obj.data.image_cut = []; % image after removing all the stars using mexCutout
            obj.data.num_sum = 1; % if we read something else, it will be altered
            obj.data.cutouts = []; % cutouts are only made from the stack, not from the full frame rate data
            
        end
        
    end
    
    methods % getters
        
        function val = getNumBatches(obj) % get the number of batches we can read, given file list and the limit given in "pars"
            
            if isempty(obj.pars.num_batches)
                val = length(obj.reader.filenames);
            else
                val = min(obj.pars.num_batches, length(obj.reader.filenames)); 
            end
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % commands
        
        function obj = run(obj, varargin)
            
            obj.pars = util.oop.full_copy(obj.pars); 
            obj.pars.scan_vars(varargin{:}); 
            
            on_cleanup = onCleanup(@obj.finishup);
            
            try
            
            obj.startup; 
            
            for ii = obj.current_batch:obj.getNumBatches
                
                if obj.brake_bit, break; end
                
                obj.batch;
                
                if ii==1, obj.first_batch; end
                
                obj.successful_batches = obj.successful_batches + 1;
                
                obj.prog.showif(ii); 
                
            end
            
            catch ME
                
                if obj.pars.use_save_results 
                    obj.write_log(ME.getReport); % make sure to also save the error in the log file 
                end
                
                if obj.successful_batches>=obj.getNumBatches || obj.failed_batch_counter>=obj.pars.max_failed_batches % run has finished all useful batches
                    warning(ME.getReport); % this is some after-processing error, we can live with it
                else % run has caught an error mid-way
                    rethrow(ME); 
                end
            end
            
        end
        
        function cont(obj)
            
            obj.run(obj.last_pars.output_vars{:}, 'reset', 0); 
            
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
            
        end
        
        function saveResults(obj)
            
            % save the lightcurves, catalog, summary text file, what else?
            
        end
        
        function str = printout(obj)
            
            str = sprintf('batch: %d/%d', obj.successful_batches, obj.getNumBatches); 
            
            str = sprintf('%s | folder: %s', str, strrep(obj.reader.dir.two_tail, '\', ' \ ')); 
            
            str = sprintf('%s | coords: %s%s', str, obj.head.RA, obj.head.Dec); 
            
            if ~isempty(obj.head.AIRMASS)
                str = sprintf('%s | a.m.= %4.2f', str, obj.head.AIRMASS); 
            end
            
        end
        
    end
    
    methods(Hidden=true) % internal calculations
        
        function startup(obj, varargin)
            
            if obj.pars.reset % all the things that must happen when a run begins (and not happen when a run continues)
                
                obj.displayInfo(sprintf('starting new run with %d batches', obj.getNumBatches));
                
                obj.reset;
                obj.prog.reset(obj.getNumBatches);
            
                obj.data.date = obj.find_folder_date;
                obj.data.camera = obj.find_camera;
                obj.data.project = obj.find_project;

                if obj.pars.use_auto_load_cal
                    obj.displayInfo(sprintf('Loading calibration file for date %s', obj.data.date)); 
                    obj.cal.loadByDate(obj.data.date, obj.data.camera, obj.data.project); 
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
            
            obj.last_pars = obj.pars; % store the parameters used in last run
            
            obj.pars = obj.backup_pars; % restore the parameters of the original object before parsing inputs
            
            obj.prog.finish;
            
            obj.lightcurves.finishup;
            
            if obj.pars.use_save_results
                obj.saveResults;
            end
            
            if obj.debug_bit, fprintf('Finished run with %d batches.\n', obj.successful_batches); end
            
        end
        
        function batch(obj)
            
            try % get the data
                
                obj.clear; % get rid of existing datasets
                
                obj.reader.batch; % load the images from file
                
                % copy the relevant data (images only)
                if ~isempty(obj.reader.stack)
                    obj.data.image = obj.reader.stack;
                    obj.data.num_sum = obj.reader.num_sum;
                elseif ~isempty(obj.reader.images)
                    if size(obj.reader.images,3)==1
                        obj.data.image = obj.reader.images;
                    else
                        error('I don''t know what to do with a data cube!'); 
                    end
                else
                    error('Could not find any images or stack in the Reader...'); 
                end
                
                obj.data.timestamps = obj.reader.timestamps(1); 
                
                if ~isempty(obj.reader.juldates)
                    obj.data.juldates = obj.reader.juldates(1); 
                else
                    % TODO: recalculate it from the timestamps...
                end
                
            catch ME
                warning(ME.getReport);
                obj.reader.advanceFile;
                return; % skip this batch, report it, and continue! 
            end
            
            % processing
            
            obj.data.image_proc = obj.cal.input(obj.data.image, 'sum', obj.data.num_sum); 
            
            if isempty(obj.data.positions)
                
                obj.findStars;
                
                if obj.pars.use_astrometry
                    obj.solveAstrometry;
                end
                
            end
            
            obj.data.ref_image = obj.data.image_proc;
            obj.data.ref_positions = obj.data.positions;
            
            [obj.data.cutouts, obj.data.image_cut] = util.img.mexCutout(obj.data.image_proc, obj.data.positions, obj.pars.cut_size, NaN, NaN); % replace and fill up using NaNs (it is safe, the processed image is single precision
            
            obj.phot.input(obj.data.cutouts, 'positions', obj.data.positions, 'timestamps', obj.data.timestamps, 'juldates', obj.data.juldates); 
            if ~isempty(obj.phot.gui) && obj.phot.gui.check, obj.phot.gui.update; end

            obj.adjustPositions; % use the centroids to push the cutouts a little, and the fluxes to check if we still see the stars (if not, call re-align)
            
            obj.flux_buf.input(obj.phot.fluxes);
            
            obj.lightcurves.getData(obj.phot); 
            if obj.lightcurves.gui.check, obj.lightcurves.gui.update; end
            
            if ~isempty(obj.func)

                if isa(obj.func, 'function_handle')
                    feval(obj.func, obj);
                elseif iscell(obj.func)
                    for jj = 1:length(obj.func)
                        feval(obj.func{jj}, obj);
                    end
                end
            end

            if obj.pars.use_save_results
                obj.write_log(obj.log_name);
            end
            
            if obj.pars.use_fits_save
                % TODO: need to add this...
            end
            
            if ~isempty(obj.gui)
                obj.gui.update;
            end
            
        end
        
        function first_batch(obj)
            
        end
        
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
        
        function write_log(obj)
            
        end
        
        function val = getAverageWidth(obj) % flux weighted average PSF width
        
            if isempty(obj.phot.widths)
                val = [];
            else
                val = util.vec.weighted_average(obj.phot.widths, obj.phot.fluxes, 2);
            end
            
        end
        
        function val = getAverageOffsets(obj) % flux weighted average offsets [dx,dy]
        
            if isempty(obj.phot.offsets_x) || isempty(obj.phot.offsets_y)
                val = [];
            else
                dx = util.vec.weighted_average(obj.phot.offsets_x, obj.phot.fluxes, 2);
                dy = util.vec.weighted_average(obj.phot.offsets_y, obj.phot.fluxes, 2);
                val = [dx, dy];
            end
            
        end
        
        function findStars(obj)
            
            T = util.img.quick_find_stars(obj.data.image_proc, 'psf', obj.pars.initial_guess_psf_width, 'number', obj.pars.num_stars_filter_kernel,...
               'dilate', obj.pars.cut_size-5, 'saturation', obj.pars.saturation.*obj.data.num_sum, 'unflagged', 1); 
            
            % estimate the real PSF width from found stars
            C = util.img.mexCutout(obj.data.image_proc, T.pos, obj.pars.cut_size, NaN, NaN); 
            
            s = util.img.photometry2(C, 'aperture', obj.pars.initial_guess_psf_width.*3); 
            
            obj.data.psf_width = util.vec.weighted_average(s.forced_photometry.width, s.forced_photometry.flux, 2);
            
            T = util.img.quick_find_stars(obj.data.image_proc, 'psf', obj.data.psf_width, 'number', obj.pars.num_stars,...
               'dilate', obj.pars.cut_size-5, 'saturation', obj.pars.saturation.*obj.data.num_sum, 'unflagged', 0); 
            
            if isempty(T)
                error('Could not find any stars using quick_find_stars!');
            end

            obj.data.positions = T.pos;
            
            obj.data.found_stars = T; % store a copy of this table

        end
        
        function adjustPositions(obj)
            
            import util.text.cs;
            
            obj.checkRealign;

            if cs(obj.pars.lock_adjust, 'all') || obj.pars.use_astrometry==0 % by default, "stars" is replaced by "all" when not running astrometry
                obj.data.positions = double(obj.data.positions + obj.getAverageOffsets);
            elseif cs(obj.pars.lock_adjust, 'none')
                obj.data.positions = double(obj.data.positions + [obj.phot.offsets_x, obj.phot.offsets_y]); % check dimensionality in case we have multiple apertures! 
            elseif cs(obj.pars.lock_adjust, 'stars') % assume astrometry has succeeded
                
                star_idx = ~isnan(obj.cat.magnitudes); % all stars have a non NaN magnitude in GAIA
                
                new_pos_free = double(obj.data.positions + [obj.phot.offsets_x, obj.phot.offsets_y]); % check dimensionality in case we have multiple apertures! 
                
                new_pos_lock = double(obj.data.positions + obj.getAverageOffsets); 
                
                obj.data.positions = star_idx.*new_pos_lock + (~star_idx).*new_pos_free; % logical indexing to add the correct adjusted position to stars and non-stars
                
            else
                error('Unknown "lock_adjust" option: "%s". Choose "all", "none" or "stars". ', obj.pars.lock_adjust); 
            end
            
            if size(obj.data.positions,3)>1
                error('dimensionality issue...'); 
            end
            
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

                new_fluxes = obj.phot.fluxes;
                new_fluxes(isnan(mean_fluxes)) = [];
                mean_fluxes(isnan(mean_fluxes)) = [];

                flux_lost = sum(new_fluxes<0.5*mean_fluxes)>0.5*numel(mean_fluxes); % lost half the flux in more than half the stars...
                % want to add more tests...?

                val = ~flux_lost;
                
            end
            
        end
        
        function checkRealign(obj) % check the stars are visible, if not, realign and re-do photometry
            
            if ~obj.checkFluxes
                
                if obj.debug_bit, disp('Lost star positions, using quick_align and doing photometry again...'); end

                [~,shift] = util.img.quick_align(obj.data.image_proc, obj.data.ref_image);
                obj.data.positions = double(obj.data.ref_positions + flip(shift));

                obj.data.cutouts = util.img.mexCutout(obj.data.image_proc, obj.data.positions, obj.pars.cut_size, NaN, NaN); % replace and fill up using NaNs (it is safe, the processed image is single precision

                obj.phot.input(obj.data.cutouts, 'positions', obj.data.positions, 'timestamps', obj.data.timestamps, 'juldates', obj.data.juldates); % run photometry again
                if ~isempty(obj.phot.gui) && obj.phot.gui.check, obj.phot.gui.update; end

            end
            
        end
        
        function solveAstrometry(obj)
            
            obj.displayInfo('Running astrometry'); 
            
            if util.text.cs(obj.head.cam_name, 'balor')
                obj.cat.input_rotation = -60; 
            elseif util.text.cs(obj.head.cam_name, 'zyla')
                obj.cat.input_rotation = -15; 
            end

            obj.cat.input(obj.data.positions); 

            if obj.cat.success==0
                error('Could not find astrometric solution!'); 
            end

            if obj.pars.use_remove_bad_matches

                bad_idx = isnan(obj.cat.magnitudes) & obj.data.found_stars.snr<obj.pars.bad_match_min_snr;

                obj.cat.data = obj.cat.data(~bad_idx, :); 
                obj.cat.magnitudes = obj.cat.magnitudes(~bad_idx); 
                obj.cat.coordinates = obj.cat.coordinates(~bad_idx, :); 
                obj.cat.temperatures = obj.cat.temperatures(~bad_idx); 
                obj.cat.positions = obj.cat.positions(~bad_idx, :); 
                obj.data.positions = obj.data.positions(~bad_idx, :); 
                obj.data.found_stars = obj.data.found_stars(~bad_idx, :); 

            end
            
            if obj.cat.success
                str = sprintf('Found astrometric solution: %s%s', head.Ephemeris.deg2hour(obj.cat.central_RA), head.Ephemeris.deg2sex(obj.cat.central_Dec)); 
            else
                str = 'Cannot find an astrometric solution!'; 
            end

            obj.displayInfo(str);
            disp(str); 

        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('autodyn', false); 
            input.scan_vars(varargin{:});
            
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
            if obj.pars.use_display_rect_text, delete(findobj(input.ax, 'type', 'text')); end
            
            if ~isempty(obj.data.positions)

                P = obj.data.positions;

                for ii = 1:min(size(P,1), obj.pars.display_num_rect_stars)
                    if obj.pars.use_display_rect_text, text(P(ii,1), P(ii,2), sprintf('clip %d', ii),'FontSize', 16, 'Parent', input.ax); end
%                     if obj.pars.use_display_rect_text, text(P(ii,1), P(ii,2), sprintf('S/N %4.2f', obj.data.found_stars.snr(ii)),'FontSize', 16, 'Parent', input.ax); end
                    if isnan(obj.cat.magnitudes(ii))
                        rectangle('Position', [P(ii,:)-0.5-obj.pars.cut_size/2 obj.pars.cut_size obj.pars.cut_size], 'Parent', input.ax, 'EdgeColor', 'white'); % white is for no match to GAIA
                    elseif obj.data.found_stars.flag(ii)==1
                        rectangle('Position', [P(ii,:)-0.5-obj.pars.cut_size/2 obj.pars.cut_size obj.pars.cut_size], 'Parent', input.ax, 'EdgeColor', 'red'); % red is for saturated stars
                    elseif obj.data.found_stars.flag(ii)==2
                        rectangle('Position', [P(ii,:)-0.5-obj.pars.cut_size/2 obj.pars.cut_size obj.pars.cut_size], 'Parent', input.ax, 'EdgeColor', 'yellow'); % yellow is for point sources
                    elseif obj.data.found_stars.flag(ii)==3
                        rectangle('Position', [P(ii,:)-0.5-obj.pars.cut_size/2 obj.pars.cut_size obj.pars.cut_size], 'Parent', input.ax, 'EdgeColor', 'green'); % green is for extended sources
                    else    
                        rectangle('Position', [P(ii,:)-0.5-obj.pars.cut_size/2 obj.pars.cut_size obj.pars.cut_size], 'Parent', input.ax, 'EdgeColor', 'black'); % black is good stars
                    end
                end
                
            end
            
        end
        
        function showFindingChart(obj, varargin)
            
        end
        
        function displayInfo(obj, str)
            
            try 
                
                if ~isempty(obj.gui) && obj.gui.check
                    obj.gui.panel_info.button_info.String = str;
                end
            
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = img.gui.ProcGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
    end
    
end

 