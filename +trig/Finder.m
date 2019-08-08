classdef Finder < handle

    properties(Transient=true)
        
        gui@trig.gui.FinderGUI;
        hist_fig;
        bank@occult.ShuffleBank;
        
    end
    
    properties % objects
        
        pars@head.Parameters;
        cat@head.Catalog;
        
        cal@trig.Calibrator;
        
        all_events@trig.Event;
        new_events@trig.Event;
        last_events@trig.Event;
        
        var_buf@util.vec.CircularBuffer;
        
        phot_pars; % a struct with some housekeeping about how the photometry was done
        
        sim_bank@occult.FilterBank;
        sim_events = {};
        
    end
    
    properties % inputs/outputs
        
        black_list_stars;
        black_list_batches;
         
        timestamps;
        cutouts;
        positions;
        stack;
        batch_index;
        filename;
        
        fluxes;
        errors;
        areas;
        backgrounds;
        variances;
        offsets_x;
        offsets_y;
        widths;
        bad_pixels;
        
        prev_timestamps;
        prev_cutouts;
        prev_positions;
        prev_stack;
        prev_batch_index;
        prev_filename;
        
        prev_fluxes;
        prev_errors
        prev_areas;
        prev_backgrounds;
        prev_variances;
        prev_offsets_x;
        prev_offsets_y;
        prev_widths;
        prev_bad_pixels;
        
        dt; 
        coverage_total;
        coverage_lost;
        star_hours_total;
        star_hours_lost;
        snr_values;
        
    end
    
    properties % switches/controls
        
        min_star_snr = 5; % stars with lower S/N are not even tested for events
        threshold = 7.5; % threshold (in units of S/N) for peak of event 
        time_range_thresh = -2.5; % threshold for including area around peak (in continuous time)
        kern_range_thresh = -1; % area threshold (in kernels, discontinuous) NOTE: if negative this will be relative to "threshold"
        star_range_thresh = -1; % area threshold (in stars, discontinuous) NOTE: if this is higher than "threshold" there will be no area around peak
        
        % additional cuts on events
        max_events = 5; % how many events can we have triggered on the same 2-batch window?
        max_stars = 5; % how many stars can we afford to have triggered at the same time? 
        max_frames = 50; % maximum length of trigger area (very long events are disqualified)
        max_num_nans = 1;
        max_corr = 0.75;
        
        use_var_buf = 1; % use PSD tracking! 
        
        num_hits_black_list = 4;
        
        use_conserve_memory = 1;
        
        use_sim = 0; % run simulation on incoming data... 
        
        frame_rate = 25; % if timestamps are not given explicitely
        
        display_event_idx = [];
        use_display_kept_events = 0;
        
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        total_batches;
        kept_events;
        
    end
    
    properties(Hidden=true)
       
        t_end;
        t_end_stamp;
        used_background_sub;
        aperture;
        gauss_sigma;
        
        version = 1.02;
        
    end
    
    methods % constructor
        
        function obj = Finder(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.Finder')
                if obj.debug_bit, fprintf('Finder copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Finder constructor v%4.2f\n', obj.version); end
            
                obj.cal = trig.Calibrator;
                
                obj.var_buf = util.vec.CircularBuffer;
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.all_events = trig.Event.empty;
            obj.new_events = trig.Event.empty;
            obj.last_events = trig.Event.empty;
            
            obj.black_list_stars = [];
            obj.black_list_batches = [];
            
            obj.cal.reset;
            obj.var_buf.reset;
            
            obj.prev_fluxes = [];
            obj.prev_errors = [];
            obj.prev_areas = [];
            obj.prev_backgrounds = [];
            obj.prev_variances = [];
            obj.prev_offsets_x = [];
            obj.prev_offsets_y = [];
            obj.prev_widths = [];
            obj.prev_bad_pixels = [];
            obj.prev_timestamps = [];
            obj.prev_cutouts = [];
            obj.prev_positions = [];
            obj.prev_stack = [];
            obj.prev_batch_index = [];
            obj.prev_filename = [];
            
            obj.dt = []; 
            obj.coverage_total = 0;
            obj.coverage_lost = 0;
            obj.star_hours_total = 0;
            obj.star_hours_lost = 0;
            obj.snr_values = [];
            
            obj.display_event_idx = [];
            
            if ~isempty(obj.gui) && obj.gui.check
                delete(obj.gui.panel_image.Children);
                obj.gui.update;
            end
            
            obj.sim_events = {};
            
            if ~isempty(obj.sim_bank)
                obj.sim_bank.reset;
            end
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            obj.fluxes = [];
            obj.errors = [];
            obj.areas = [];
            obj.backgrounds = [];
            obj.variances = [];
            obj.offsets_x = [];
            obj.offsets_y = [];
            obj.widths = [];
            obj.bad_pixels = [];
            obj.timestamps = [];
            obj.cutouts = [];
            obj.positions = [];
            obj.stack = [];
            obj.batch_index = [];
            obj.filename = [];
            
            obj.cal.clear;
            if ~isempty(obj.bank)
                obj.bank.clear;
            end
            
        end
        
    end
    
    methods % getters
        
        function val = getTimeThresh(obj)
            
            if obj.time_range_thresh>0
                val = obj.time_range_thresh;
            else
                val = obj.threshold + obj.time_range_thresh;
            end
            
        end
        
        function val = getKernThresh(obj)
            
            if obj.time_range_thresh>0
                val = obj.kern_range_thresh;
            else
                val = obj.threshold + obj.kern_range_thresh;
            end
            
        end
        
        function val = getStarThresh(obj)
            
            if obj.star_range_thresh>0
                val = obj.star_range_thresh;
            else
                val = obj.threshold + obj.star_range_thresh;
            end
            
        end
        
        function val = get.total_batches(obj)
            
            val = length(obj.snr_values);
            
        end
        
        function val = get.kept_events(obj)
            
            val = obj.all_events([obj.all_events.keep]==1);
            
        end
        
        function val = num_events(obj)
            
            val = length(obj.all_events);
            
        end
        
        function val = num_kept(obj)
            
            val = length(obj.kept_events);
            
        end
        
        function val = obs_pars_str(obj)
            
            if isempty(obj.pars)
                val = '';
            else
                if isempty(obj.pars.STARTTIME)
                    date_str = '';
                else
                    date_str = obj.pars.STARTTIME(1:10);
                end
                val = sprintf('date: %s', date_str); % need to expand this...
            end
            
        end
        
        function val = phot_pars_str(obj)
            
            if isempty(obj.phot_pars)
                val = '';
            else
                
                val = '';
                
                val = [val 'b/g sub: ' num2str(obj.phot_pars.used_bg_sub)];
                
                if strcmpi(obj.phot_pars.signal_method, 'aperture')
                    val = [val sprintf(' | aperture %4.2f pix radius', obj.phot_pars.radius)];
                end
                
                if strcmpi(obj.phot_pars.signal_method, 'aperture')
                    val = [val sprintf(' | annulus %4.2f-%4.2f pixels', obj.phot_pars.annulus, obj.phot_pars.annulus_outer)];
                end
                
                val = [val ' | iter= ' num2str(obj.phot_pars.iterations)];
                
            end
            
        end
        
        function val = event_pars_str(obj)
            
            val = '';
            
        end
        
        function val = getSimulatedSNR(obj)
           
        end
        
    end
    
    methods % setters
        
        function set.display_event_idx(obj, val)
            
            obj.display_event_idx = val;
            if ~isempty(obj.display_event_idx) && ~isempty(obj.gui) && obj.gui.check
                obj.all_events(obj.display_event_idx).show('parent', obj.gui.panel_image);
            end

        end
        
    end
    
    methods % calculations
        
        function loadFilterBank(obj)

            f = fullfile(getenv('DATA'), '/WFAST/saved/FilterBankShuffle.mat');
            if exist(f, 'file')
                load(f, 'bank');
                obj.bank = bank;
            else
                error('Cannot load kernels from ShuffleBank object'); 
            end

        end
        
        function input(obj, varargin)
                
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('fluxes', []);
            input.input_var('errors', []);
            input.input_var('areas', []);
            input.input_var('backgrounds', []);
            input.input_var('variances', []);
            input.input_var('offsets_x', [], 'offset_x', 'dx');
            input.input_var('offsets_y', [], 'offset_y', 'dy');
            input.input_var('widths', []);
            input.input_var('bad_pixels', []);
            input.input_var('aperture', [], 'radius');
            input.input_var('gauss_sigma', [], 'sigma', 'gaussian_sigma');
            input.input_var('timestamps', []); 
            input.input_var('cutouts', []);
            input.input_var('positions', []);
            input.input_var('stack', []);
            input.input_var('batch_index', [], 'batch_idx', 'batch_number');
            input.input_var('filename', '');
            input.input_var('t_end', [], 8);
            input.input_var('t_end_stamp', [], 8);
            input.input_var('used_background_sub', []); 
            input.input_var('phot_pars', [], 'pars_struct');
            input.input_var('lightcurves', []);
            input.input_var('index', []);
            input.input_var('length', 100);
            input.scan_vars(varargin{:});
            
            if isempty(obj.bank)
                obj.loadFilterBank;
            end
            
            if obj.use_sim
                
                if isempty(obj.sim_bank)
                    obj.sim_bank = occult.FilterBank;
                end
                
                if isempty(obj.sim_bank.bank)
                    obj.sim_bank.makeBank;
                end
                
                if isempty(obj.sim_bank.filtered_bank) || numel(obj.sim_bank.filtered_bank)~=numel(obj.bank.kernels)*obj.sim_bank.num_pars
                    
                    if obj.debug_bit, fprintf('Cross filtering all kernels in ShuffleBank (%d) with all templates in FilterBank (%d)...', size(obj.bank.kernels,2), obj.sim_bank.num_pars); end
                    
                    obj.sim_bank.filtered_bank = util.vec.convolution(obj.bank.kernels, permute(obj.sim_bank.bank-1, [1,6,2,3,4,5])); 
                    
                end
                
            end
            
            if ~isempty(input.lightcurves) && ~isempty(input.index) && ~isempty(input.length) % parse a "lightcurve" struct/object
                
                input.fluxes = input.lightcurves.fluxes(input.index:input.index+input.length-1, :);
                input.errors= input.lightcurves.errors(input.index:input.index+input.length-1, :);
                input.areas = input.lightcurves.areas(input.index:input.index+input.length-1, :);
                input.backgrounds = input.lightcurves.backgrounds(input.index:input.index+input.length-1, :);
                input.variances = input.lightcurves.variances(input.index:input.index+input.length-1, :);
                
                if isfield(input.lightcurves, 'offsets_x') && ~isempty(input.lightcurves.offsets_x)
                    input.offsets_x = input.lightcurves.offsets_x(input.index:input.index+input.length-1, :);
                else % reconstruct the offsets from centroids by assuming a constant base-pixel for the whole batch... 
                    input.offsets_x = input.lightcurves.centroids_x(input.index:input.index+input.length-1, :)...
                        - round(median(input.lightcurves.centroids_x(input.index:input.index+input.length-1, :), 'omitnan'));
                end
                
                if isfield(input.lightcurves, 'offsets_y') && ~isempty(input.lightcurves.offsets_y)
                    input.offsets_y = input.lightcurves.offsets_y(input.index:input.index+input.length-1, :);
                else % reconstruct the offsets from centroids by assuming a constant base-pixel for the whole batch... 
                    input.offsets_y = input.lightcurves.centroids_y(input.index:input.index+input.length-1, :)...
                        - round(median(input.lightcurves.centroids_y(input.index:input.index+input.length-1, :), 'omitnan')); 
                end
                    
                input.widths = input.lightcurves.widths(input.index:input.index+input.length-1, :);
                input.bad_pixels = input.lightcurves.bad_pixels(input.index:input.index+input.length-1, :);
                
            end
            
            obj.clear;
            
            obj.fluxes = input.fluxes;
            obj.errors = input.errors;
            obj.areas = input.areas;
            obj.backgrounds = input.backgrounds;
            obj.variances = input.variances;
            obj.offsets_x = input.offsets_x;
            obj.offsets_y = input.offsets_y;
            obj.widths = input.widths;
            obj.bad_pixels = input.bad_pixels;
            obj.timestamps = input.timestamps;
            obj.cutouts = input.cutouts;
            obj.positions = input.positions;
            obj.stack = input.stack;
            obj.batch_index = input.batch_index;
            obj.filename = input.filename;
            obj.phot_pars = input.phot_pars;
            obj.t_end = input.t_end;
            obj.t_end_stamp = input.t_end_stamp;
            obj.aperture = input.aperture;
            obj.gauss_sigma = input.gauss_sigma;
            obj.used_background_sub = input.used_background_sub;
            
            if all(obj.timestamps==0)
%                 obj.timestamps = []; % in case we didn't properly store the times
                if isempty(obj.prev_timestamps)
                    obj.timestamps = (1:size(input.fluxes,1))'./obj.frame_rate;
                else
                    obj.timestamps = obj.prev_timestamps(end)+(1:size(input.fluxes,1))'./obj.frame_rate;
                end
            end
            
            %%%%%%%%%%%%%%% done parsing inputs %%%%%%%%%%%%%%%%%%%%%%%
            
            if ~isempty(obj.prev_fluxes) % skip first batch! 
            
                t = tic;
                obj.cal.input(vertcat(obj.prev_fluxes, obj.fluxes), vertcat(obj.prev_errors, obj.errors), vertcat(obj.prev_timestamps, obj.timestamps)); 
                if obj.debug_bit>2, fprintf('Calibration time: %f seconds.\n', toc(t)); end
                
                t = tic;
                obj.bank.input(obj.cal.fluxes_detrended, obj.cal.stds_detrended, obj.cal.timestamps); % use the filter bank on the fluxes

                if nnz(~isnan(obj.bank.fluxes_filtered))==0
                    error('Filtered fluxes in ShuffleBank are all NaN!');
                end
                
                if obj.debug_bit>2, fprintf('Filtering time: %f seconds.\n', toc(t)); end
                
                if obj.use_sim
                    
                    obj.sim_bank.filtered_index = 0; % start by running events without any added occultations
                    
                    snr_vec = zeros(1,obj.sim_bank.num_pars);
                    
                    for ii = 1:1+obj.sim_bank.num_pars
                        
                        obj.findEvents;
                        
                        if obj.sim_bank.filtered_index==0
                            obj.last_events = obj.new_events;
                            obj.all_events = [obj.all_events obj.new_events];
                        else % store new events in different bins for different simulation parameters
%                             if length(obj.sim_events)<obj.sim_bank.filtered_index || isempty(obj.sim_events{obj.sim_bank.filtered_index})
%                                 obj.sim_events{obj.sim_bank.filtered_index} = obj.new_events;
%                             else
%                                 obj.sim_events{obj.sim_bank.filtered_index} = [obj.sim_events{obj.sim_bank.filtered_index} obj.new_events];
%                             end
                        end
                        
                        if obj.sim_bank.filtered_index>0
                            if ~isempty(obj.new_events)
                                [mx1, idx] = max([obj.new_events.snr]); % find the best event (if there are more than one)
                                mx2 = obj.new_events(idx).star_snr;
                                snr_vec(obj.sim_bank.filtered_index) = mx1/mx2; % save the normalized S/N
                            end
                        end
                        
                        obj.sim_bank.filtered_index = obj.sim_bank.filtered_index + 1; % each time the index points at a different template to add
                    
                    end
                    
                    obj.sim_bank.snr_sim_full(end+1,:) = snr_vec;
                    
                else
                    
                    obj.findEvents; % just find real events! 
                    obj.last_events = obj.new_events;
                    obj.all_events = [obj.all_events obj.new_events];
                    
                end
                
                
            end
            
            % store these for next time
            obj.prev_fluxes = obj.fluxes;
            obj.prev_errors = obj.errors;
            obj.prev_areas = obj.areas;
            obj.prev_backgrounds = obj.backgrounds;
            obj.prev_variances = obj.variances;
            obj.prev_offsets_x = obj.offsets_x;
            obj.prev_offsets_y = obj.offsets_y;
            obj.prev_widths = obj.widths;
            obj.prev_bad_pixels = obj.bad_pixels;
            obj.prev_timestamps = obj.timestamps;
            obj.prev_cutouts = obj.cutouts;
            obj.prev_positions = obj.positions;
            obj.prev_stack = obj.stack;
            obj.prev_batch_index = obj.batch_index;
            obj.prev_filename = obj.filename;
            
            if ~isempty(obj.bank.fluxes_filtered)
                obj.var_buf.input(var(obj.bank.fluxes_filtered, [], 1, 'omitnan')); % PSD tracking: keep a running buffer of the variance of previous filter results
            end
            
        end
        
        function findEvents(obj)
            
            obj.new_events = trig.Event.empty;
           
            ff = obj.bank.fluxes_filtered; % dim 1 is time, dim 2 is kernels, dim 3 is stars
            
            good_stars = obj.cal.star_snrs>obj.min_star_snr;
            if all(~good_stars)
                return;
            end
            
            t = tic; 
            
            if obj.use_sim && obj.sim_bank.filtered_index>0
                
                t_sim = tic;
                
                for ii = 1:1e6
                    star_index_sim = randi(size(ff,3)); % random choice of star
                    if good_stars(star_index_sim), break; end % repick a new star if this one is low S/N
                end
                
                if obj.debug_bit>1 
                    s = obj.sim_bank.sim_pars;
                    fprintf('Adding flux with R= %4.2f, r= %4.2f, b= %4.2f, v= %4.2f to flux of star %d!\n', s.R, s.r, s.b, s.v, star_index_sim);
                end
                
                % add the occultation on top of the data!
                ff_sim = util.img.pad2size(obj.sim_bank.filtered_bank(:,:,obj.sim_bank.filtered_index), size(ff)); 
                ff_sim = ff_sim.*obj.cal.star_snrs(star_index_sim); 
                max_shift = size(ff,1)-5;
                shift_frames = randi(max_shift)-floor(max_shift/2); % shift by a random number of frames
                ff_sim = util.img.shift(ff_sim, 0, shift_frames, 0);
                ff = ff(:,:,star_index_sim) + ff_sim;
                
                
                if obj.debug_bit>3, fprintf('Adding simulated template time: %f seconds.\n', toc(t_sim)); end

            else
                star_index_sim = []; % this indicates we are not running sim! 
            end
            
            if obj.use_var_buf && ~obj.var_buf.is_empty
     
                t_psd = tic;

                if isempty(star_index_sim)
                    real_rms = sqrt(obj.var_buf.mean);
                else
                    real_rms = sqrt(obj.var_buf.mean(:,:,star_index_sim));
                end
                
                ff = ff./real_rms;
                
                if obj.debug_bit>3, fprintf('PSD drift correction time: %f seconds.\n', toc(t_psd)); end

            end
            
            obj.dt = median(diff(obj.timestamps));
            
            for ii = 1:obj.max_events
                
                t_max = tic;
                
                if isempty(star_index_sim)
                    [mx, idx] = util.stat.maxnd(abs(ff.*util.vec.topages(good_stars))); % note we are triggering on negative and positive events
                else
                    [mx, idx] = util.stat.max2(abs(ff)); % scan only the star we added simulation to
                    idx(3) = 1; % in sim-mode we only scan a 2D ff matrix, but we need to index into it as a 3D matrix later on...  
                end
                
                if obj.debug_bit>3, fprintf('Find max time: %f seconds.\n', toc(t_max)); end

                if ii==1
                    if isempty(star_index_sim) % only store S/N values when not in sim-mode
                        obj.snr_values(end+1) = ff(idx(1),idx(2),idx(3));
                    end
                end
                
                if mx>=obj.threshold
                
                    t_ev = tic;
                    
                    ev = trig.Event;
                    ev.snr = mx; % note this is positive even for negative filter responses! 

                    ev.time_index = idx(1); 
                    ev.kern_index = idx(2);
                    
                    if isempty(star_index_sim) % (note: in sim mode, idx(3)==1 by definition, but star_index is not!)
                        ev.star_index = idx(3);
                    else
                        ev.star_index = star_index_sim;
                    end

                    ev.time_indices = obj.findTimeRange(ff, idx(1), idx(3)); % find continuous area that is above time_range_thresh

                    ev.kern_indices = find(max(abs(ff(ev.time_indices, :, idx(3))))>obj.getKernThresh);

                    ev.star_indices = find(max(abs(ff(ev.time_indices, idx(2), :)))>obj.getStarThresh);

                    ev.timestamps = obj.bank.timestamps;
                    ev.time_step = obj.dt;

                    if isempty(ev.time_indices)
                        warning('time_indices for latest event are empty!')
                        ev.duration = 0;
                    else
                        ev.duration =  obj.dt + obj.bank.timestamps(ev.time_indices(end))-obj.bank.timestamps(ev.time_indices(1));
                    end

                    ev.flux_filtered = ff(:,idx(2),idx(3));
                    if obj.use_var_buf && ~is_empty(obj.var_buf)
                        ev.previous_std = sqrt(obj.var_buf.mean(1,ev.kern_index, ev.star_index));
                    end

                    ev.flux_detrended = obj.cal.fluxes_detrended(:,ev.star_index); 
                    ev.std_flux = std(ev.flux_detrended, [], 'omitnan');
                    ev.flux_raw_all = obj.cal.fluxes;
                    % somewhere around here we MUST make use of the flux errors

                    obj.new_events(end+1) = ev; % add this event to the list

                    if ~obj.use_sim || obj.sim_bank.filtered_index==0 % only store star hours when not in sim-mode
                        obj.coverage_lost = obj.coverage_lost + ev.duration; % how much time is "zeroed out"
                        obj.star_hours_lost = obj.star_hours_lost + (ev.duration)*sum(good_stars)/3600;
                    end
                    
                    if obj.debug_bit>3, fprintf('New event time: %f seconds.\n', toc(t_ev)); end
            
                    ff(ev.time_indices, :, :) = 0; % don't look at the same region twice
            
                else
                    break;
                end
                
            end
            
            if isempty(star_index_sim) % only store star hours when not in sim-mode
                obj.coverage_total = obj.coverage_total + obj.timestamps(end) - obj.timestamps(1) + obj.dt;
                obj.star_hours_total = obj.star_hours_total + (obj.timestamps(end) - obj.timestamps(1) + obj.dt)*sum(good_stars)/3600;
                sim_pars = [];
            else
                sim_pars = obj.sim_bank.sim_pars;
            end
            
            if obj.debug_bit>2, fprintf('Triggering time: %f seconds.\n', toc(t)); end

            obj.storeEventHousekeeping(sim_pars); % add some data from this batch to the triggered events
            
            if ~obj.use_sim || obj.sim_bank.filtered_index==0 % don't check for duplicates in sim-mode
                obj.checkEvents; % check each event by itself + make sure we aren't taking events that already exist
            else
                for ii = 1:length(obj.new_events), obj.new_events(ii).self_check; end % only check the events with internal checks
            end

        end
        
        function time_range = findTimeRange(obj, ff, time_index, star_index)
            
            N = size(ff,1); % time length
            
            thresh = obj.getTimeThresh;
            time_range = [];

            for jj = 0:N % go backward in time

                idx = time_index - jj;

                if idx<1, break; end

                if any(abs(ff(idx, :, star_index))>=thresh)
                    time_range = [time_range, idx];
                else 
                    break;
                end

            end

            time_range = flip(time_range);

            for jj = 1:N % go forward in time

                idx = time_index + jj;

                if idx>N, break; end

                if any(abs(ff(idx, :, star_index))>=thresh)
                    time_range = [time_range, idx];
                else 
                    break;
                end

            end

        end
        
        function storeEventHousekeeping(obj, sim_pars) % store all the additional metadata about this event
            
            if nargin<2 || isempty(sim_pars)
                sim_pars = [];
            end
            
            t = tic;
            
            for ii = 1:length(obj.new_events)
                
                ev = obj.new_events(ii);
                
                ev.serial = length(obj.all_events) + ii; % just keep a running index
                
                % fluxes and timestamps
                ev.is_positive = obj.bank.fluxes_filtered(ev.time_index, ev.kern_index, ev.star_index)>0;
                
                ev.peak_timestamp = ev.timestamps(ev.time_index);
                ev.best_kernel = obj.bank.kernels(:,ev.kern_index);
                
                ev.threshold = obj.threshold;
                ev.used_background_sub = obj.used_background_sub;
                ev.time_range_thresh = obj.time_range_thresh;
                ev.kern_range_thresh = obj.kern_range_thresh;
                ev.star_range_thresh = obj.star_range_thresh;
                ev.star_snr = obj.cal.star_snrs(ev.star_index);
                
                % housekeeping data from photometry 
                a = vertcat(obj.prev_areas, obj.areas);
                ev.areas_at_peak = a(ev.frame_index, :);
                ev.areas_at_star = a(:, ev.star_index);
                ev.areas_time_average = mean(a, 1, 'omitnan');
                ev.areas_star_average = mean(a, 2, 'omitnan');
                
                b = vertcat(obj.prev_backgrounds, obj.backgrounds);
                ev.backgrounds_at_peak = b(ev.frame_index, :);
                ev.backgrounds_at_star = b(:, ev.star_index);
                ev.backgrounds_time_average = mean(b, 1, 'omitnan');
                ev.backgrounds_star_average = mean(b, 2, 'omitnan');
                
                v = vertcat(obj.prev_variances, obj.variances);
                ev.variances_at_peak = v(ev.frame_index, :);
                ev.variances_at_star = v(:, ev.star_index);
                ev.variances_time_average = mean(v, 1, 'omitnan');
                ev.variances_star_average = mean(v, 2, 'omitnan');
                
                dx = vertcat(obj.prev_offsets_x, obj.offsets_x);
                ev.offsets_x_at_peak = dx(ev.frame_index, :);
                ev.offsets_x_at_star = dx(:, ev.star_index);
                ev.offsets_x_time_average = mean(dx, 1, 'omitnan');
                ev.offsets_x_star_average = mean(dx, 2, 'omitnan');
                
                dy = vertcat(obj.prev_offsets_y, obj.offsets_y);
                ev.offsets_y_at_peak = dy(ev.frame_index, :);
                ev.offsets_y_at_star = dy(:, ev.star_index);
                ev.offsets_y_time_average = mean(dy, 1, 'omitnan');
                ev.offsets_y_star_average = mean(dy, 2, 'omitnan');
                
                w = vertcat(obj.prev_widths, obj.widths);
                ev.widths_at_peak = w(ev.frame_index, :);
                ev.widths_at_star = w(:, ev.star_index);
                ev.widths_time_average = mean(w, 1, 'omitnan');
                ev.widths_star_average = mean(w, 2, 'omitnan');
                
                p = vertcat(obj.prev_bad_pixels, obj.bad_pixels);
                ev.bad_pixels_at_peak = p(ev.frame_index, :);
                ev.bad_pixels_at_star = p(:, ev.star_index);
                ev.bad_pixels_time_average = mean(p, 1, 'omitnan');
                ev.bad_pixels_star_average = mean(p, 2, 'omitnan');
                
                ev.aperture = obj.aperture;
                ev.gauss_sigma = obj.gauss_sigma;
                
                if isempty(sim_pars)
                    if ~isempty(obj.prev_cutouts), ev.cutouts_first = obj.prev_cutouts(:,:,:,ev.star_index); end
                    if ~isempty(obj.cutouts), ev.cutouts_second = obj.cutouts(:,:,:,ev.star_index); end
                    if ~isempty(obj.prev_positions), ev.positions_first = obj.prev_positions; end
                    if ~isempty(obj.positions), ev.positions_second = obj.positions; end
                    if ~isempty(obj.prev_stack), ev.stack_first = obj.prev_stack; end
                    if ~isempty(obj.stack), ev.stack_second = obj.stack; end
                    
                    ev.batch_index_first = obj.prev_batch_index;
                    ev.batch_index_second = obj.batch_index;
                    ev.filename_first = obj.prev_filename;
                    ev.filename_second = obj.filename;
                    
                end
                
                ev.sim_pars = sim_pars;
                
                if ischar(obj.t_end) && isnumeric(obj.t_end_stamp)
                    ev.t_end = obj.t_end;
                    ev.t_end_stamp = obj.t_end_stamp;
                end
                
                ev.phot_pars = obj.phot_pars;
                
                ev.max_num_nans = obj.max_num_nans;
                ev.max_corr = obj.max_corr;
                
                if ~isempty(obj.cat) && ~isempty(obj.cat.data)
                    
                    if ~isempty(obj.cat.magnitudes) 
                        ev.magnitude = obj.cat.magnitudes(ev.star_index); 
                    end
                    
                    if ~isempty(obj.cat.coordinates)
                        ev.RA = head.Ephemeris.deg2hour(obj.cat.coordinates(ev.star_index, 1));
                        ev.Dec = head.Ephemeris.deg2sex(obj.cat.coordinates(ev.star_index, 2));
                    end
                    
                end
                
                % any other info that needs to be saved along with the event object?
                % ...
                
            end
            
            if obj.debug_bit>2, fprintf('Housekeeping time: %f seconds.\n', toc(t)); end

        end
        
        function checkEvents(obj)
            
            t = tic;
            
            % go over new events and check if they are duplicates
            for ii = 1:length(obj.new_events)
                
                for jj = 1:length(obj.last_events)
                    
                    if obj.new_events(ii).is_same(obj.last_events(jj))
                        
                        if obj.new_events(ii).snr<obj.last_events(jj).snr % old event is better, mark off the new one
                            obj.new_events(ii).keep = 0;
                            obj.new_events(ii).is_duplicate = 1;
                            obj.new_events(ii).addNote(['duplicate of ' num2str(obj.last_events(jj).serial)]);
                        else % new event is better, mark off the old one
                            obj.last_events(jj).keep = 0;
                            obj.last_events(jj).is_duplicate = 1;
                            obj.last_events(jj).addNote(['duplicate of ' num2str(obj.new_events(ii).serial)]);
                        end
                        
                        break;
                        
                    end
                    
                end
                
                obj.new_events(ii).self_check; % should this move to before duplicate check?
                
            end
            
            if obj.debug_bit>2, fprintf('runtime "check_new_events": %f seconds\n', toc(t)); end
            
        end
        
        function finishup(obj) % remove black list events, clear images for unkept events
            
            t = tic;
            
            stars = [obj.all_events.star_index];
            if ~isempty(stars)
                [N,E] = histcounts(stars, 'BinWidth', 1, 'BinLimits', [1 max(stars)]);
                obj.black_list_stars = [obj.black_list_stars E(N>=obj.num_hits_black_list)];
            end
            
            batches = [obj.all_events.batch_index];
            if ~isempty(batches)
                [N,E] = histcounts(batches, 'BinWidth', 1, 'BinLimits', [1 max(batches)]);
                obj.black_list_batches = [obj.black_list_batches E(N>=obj.num_hits_black_list)];
            end
            
            for ii = 1:length(obj.all_events)
                
                if ismember(obj.all_events(ii).star_index, obj.black_list_stars)
                    obj.all_events(ii).keep = 0;
                    obj.all_events(ii).is_bad_star = 1;
                    obj.all_events(ii).addNote(sprintf('star %d is on black list', obj.all_events(ii).star_index));
                end
                
                if ismember(obj.all_events(ii).batch_index, obj.black_list_stars)
                    obj.all_events(ii).keep = 0;
                    obj.all_events(ii).is_bad_batch = 1;
                    obj.all_events(ii).addNote(sprintf('batch %d is on black list', obj.all_events(ii).batch_index));
                end
                
            end
            
            if obj.use_conserve_memory
                for ii = 1:length(obj.all_events)
                    if obj.all_events(ii).keep==0
                        obj.all_events(ii).clearImages;
                    end
                end
            end
            
            if obj.debug_bit>2, fprintf('runtime "finishup": %f seconds\n', toc(t)); end
            
        end
        
        function saveEvents(obj, filename)
            
            for ii = 1:length(obj.all_events)
                if obj.all_events(ii).keep
                    events(ii) = obj.all_events(ii).reduce_memory;
                else
                    events(ii) = obj.all_events(ii);
                end
            end
            
            save(filename, 'events');
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = trig.gui.FinderGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
        function show(obj, parent)
            
            if isempty(obj.this_event)
                return;
            end
            
            if nargin<2 || isempty(parent)
                if ~isempty(obj.gui) && obj.gui.check
                    parent = obj.gui.panel_image;
                else
                    parent = gcf;
                end
            end
            
            obj.this_event.show('Parent', parent);
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.update;
            end
            
            drawnow;
            
        end
        
        function val = this_event(obj)
            
            if isempty(obj.display_event_idx) || obj.display_event_idx>obj.num_events
                val = trig.Event.empty;
            else
                val = obj.all_events(obj.display_event_idx);
            end
            
        end
        
        function display_prev_event(obj)
           
            if obj.use_display_kept_events
                if isempty(obj.kept_events)
                    obj.display_event_idx = [];
                    return;
                end
            else
                if isempty(obj.all_events)
                    obj.display_event_idx = [];
                    return;
                end
            end
            
            if isempty(obj.display_event_idx)
                idx = 1;
            else
                idx = obj.display_event_idx;
            end
            
            idx = idx - 1;
            
            if idx<1
                idx = obj.num_events;
            end
            
            if obj.use_display_kept_events
                
                for ii = idx:-1:1
                    if obj.all_events(ii).keep
                        idx = ii;
                        break;
                    end
                end
                
                if ii==1 % loop back around...
                
                    for ii = obj.num_events:-1:idx
                        if obj.all_events(ii).keep
                            idx = ii;
                            break;
                        end
                    end

                end
                
            end
            
            obj.display_event_idx = idx; % this also calls "show" from the setter
            
        end
        
        function display_next_event(obj)
           
            if obj.use_display_kept_events
                if isempty(obj.kept_events)
                    obj.display_event_idx = [];
                    return;
                end
            else
                if isempty(obj.all_events)
                    obj.display_event_idx = [];
                    return;
                end
            end
            
            if isempty(obj.display_event_idx)
                idx = 0;
            else
                idx = obj.display_event_idx;
            end
            
            idx = idx + 1;
            
            if idx>obj.num_events
                idx = 1;
            end
            
            if obj.use_display_kept_events
                
                for ii = idx:obj.num_events
                    if obj.all_events(ii).keep
                        idx = ii;
                        break;
                    end
                end
                
                if ii==obj.num_events % loop back around...
                
                    for ii = 1:idx
                        if obj.all_events(ii).keep
                            idx = ii;
                            break;
                        end
                    end

                end
                
            end

            obj.display_event_idx = idx; % this also calls "show" via the setter
            
        end
        
        function histogram(obj, parent)
            
            if nargin<2 || isempty(parent)
                
                if isempty(obj.hist_fig) || ~isa(obj.hist_fig, 'matlab.ui.Figure') || ~isvalid(obj.hist_fig)
                    obj.hist_fig = figure;
                end
                
                parent = obj.hist_fig;
                
            end
            
            delete(parent.Children);
            ax = axes('Parent', parent);
            
            histogram(ax, abs(obj.snr_values), 'BinWidth', 0.2);
            
        end
        
        function showLatest(obj, parent)
            
            if isempty(parent) || ~isvalid(parent)
                if ~isempty(obj.gui) && obj.gui.check
                    parent = obj.gui.panel_image;
                else
                    parent = gcf;
                end
            end
            
            for ii = 1:length(obj.last_events)
                
                if ii>1, pause(2); end
                
                obj.display_event_idx = obj.last_events(ii).serial;
                obj.show(parent);
            end
            
            
        end
        
        function showSimulatedSNR(obj, varargin)
            
        end
        
    end    
    
end

