classdef Finder < handle

    properties(Transient=true)
        
        gui@trig.gui.FinderGUI;
        hist_fig;
        bank@occult.ShuffleBank;
        
    end
    
    properties % objects
        
        pars@head.Parameters;
        
        cal@trig.Calibrator;
        all_events@trig.Event;
        new_events@trig.Event;
        last_events@trig.Event;
        
        phot_pars; % a struct with some housekeeping about how the photometry was done
        
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
        
        use_conserve_memory = 1;
        
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
       
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = Finder(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.Finder')
                if obj.debug_bit, fprintf('Finder copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Finder constructor v%4.2f\n', obj.version); end
            
                obj.cal = trig.Calibrator;
                                
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
            obj.bank.clear;
            
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
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function setupKernels(obj)
            
            if isempty(obj.bank)
                f = fullfile(getenv('DATA'), '/WFAST/saved/FilterBankShuffle.mat');
                if exist(f, 'file')
                    load(f, 'bank');
                    obj.bank = bank;
                else
                    error('Cannot load kernels from ShuffleBank object'); 
                end
            end
            
            % some easy parameters
%             obj.bank.R = 0;
%             obj.bank.r = 1;
%             obj.bank.b = 0;
%             obj.bank.v = 5:5:30;
%             
%             obj.bank.getLightCurves;
%             obj.filt.kernels = obj.bank.lc.flux - 1;
            
            % consider changing the default parameters of bank
%             obj.bank.makeBank;
%             obj.filt.kernels = single(reshape(obj.bank.bank-1, [size(obj.bank.bank,1), obj.bank.num_pars]));
            
            

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
            input.scan_vars(varargin{:});
            
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
            
            if all(obj.timestamps==0)
%                 obj.timestamps = []; % in case we didn't properly store the times
                if isempty(obj.prev_timestamps)
                    obj.timestamps = (1:size(input.fluxes,1))'./obj.frame_rate;
                else
                    obj.timestamps = obj.prev_timestamps(end)+(1:size(input.fluxes,1))'./obj.frame_rate;
                end
            end
            
            if ~isempty(obj.prev_fluxes) % skip first batch! 
            
                t = tic;
                obj.cal.input(vertcat(obj.prev_fluxes, obj.fluxes), vertcat(obj.prev_errors, obj.errors), vertcat(obj.prev_timestamps, obj.timestamps)); 
                if obj.debug_bit>1, fprintf('Calibration time: %f seconds.\n', toc(t)); end
                
                t = tic;
                obj.bank.input(obj.cal.fluxes_detrended, obj.cal.stds_detrended, obj.cal.timestamps); % use the filter bank on the fluxes
                if obj.debug_bit>1, fprintf('Filtering time: %f seconds.\n', toc(t)); end
                
                t = tic;
                obj.findEvents;
                if obj.debug_bit>1, fprintf('Triggering time: %f seconds.\n', toc(t)); end
                
                t = tic;
                obj.storeEventHousekeeping(input); % add some data from this batch to the triggered events
                
                obj.checkEvents; % make sure we aren't taking events that already exist
            
                obj.last_events = obj.new_events;
                obj.all_events = [obj.all_events obj.new_events];
                obj.new_events = trig.Event.empty;
            
                if obj.debug_bit>1, fprintf('Housekeeping time: %f seconds.\n', toc(t)); end
            
                t = tic;
                
                if obj.debug_bit>1, fprintf('Checking time: %f seconds.\n', toc(t)); end
                
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
            
        end
        
        function findEvents(obj)
            
            obj.new_events = trig.Event.empty;
            
            ff = obj.bank.fluxes_filtered; % dim 1 is time, dim 2 is kernels, dim 3 is stars
            
            good_stars = obj.cal.star_snrs>obj.min_star_snr;
            
            [~,idx] = util.stat.maxnd(abs(ff.*util.vec.topages(good_stars))); % only get the best S/N for first batch (for record keeping)
            
            if isempty(idx)
                return; % this can happen if all stars are below the S/N threshold
            end
            
            obj.snr_values(end+1) = ff(idx(1),idx(2),idx(3));
            
            obj.dt = median(diff(obj.timestamps));
            
            for ii = 1:obj.max_events
                
                [mx, idx] = util.stat.maxnd(abs(ff.*util.vec.topages(good_stars))); % note we are triggering on negative and positive events
                
                if mx<obj.threshold, break; end 
                
                ev = trig.Event;
                ev.snr = mx; % note this is positive even for negative filter responses! 
                
                ev.time_index = idx(1); 
                ev.kern_index = idx(2);
                ev.star_index = idx(3);
                
                ev.time_indices = obj.findTimeRange(ff, ev.time_index, ev.star_index); % find continuous area that is above time_range_thresh
                
                ev.kern_indices = find(max(abs(ff(ev.time_indices, :, ev.star_index)))>obj.getKernThresh);
                
                ev.star_indices = find(max(abs(ff(ev.time_indices, ev.kern_index, :)))>obj.getStarThresh);
                
                ev.timestamps = obj.bank.timestamps;
                ev.time_step = obj.dt;
                ev.duration =  obj.dt + obj.bank.timestamps(ev.time_indices(end))-obj.bank.timestamps(ev.time_indices(1));
                
                ev.flux_filtered = obj.bank.fluxes_filtered(:,ev.kern_index,ev.star_index);
%                 ev.flux_raw_all = permute(obj.bank.fluxes, [1,3,2]);
%                 ev.stds_raw_all = permute(obj.bank.stds, [1,3,2]);
                ev.flux_detrended = obj.cal.fluxes_detrended(:,ev.star_index); 
                ev.std_flux = std(ev.flux_detrended, [], 'omitnan');
                ev.flux_raw_all = obj.cal.fluxes;
                % somewhere around here we MUST make use of the flux errors
                
                obj.new_events(end+1) = ev; % add this event to the list
                
                obj.coverage_lost = obj.coverage_lost + ev.duration; % how much time is "zeroed out"
                obj.star_hours_lost = obj.star_hours_lost + (ev.duration)*sum(good_stars)/3600;
                
                ff(ev.time_indices, :, :) = 0; % don't look at the same region twice
                
            end
            
            obj.coverage_total = obj.coverage_total + obj.timestamps(end) - obj.timestamps(1) + obj.dt;
            obj.star_hours_total = obj.star_hours_total + (obj.timestamps(end) - obj.timestamps(1) + obj.dt)*sum(good_stars)/3600;
            
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
        
        function storeEventHousekeeping(obj, input) % store all the additional metadata about this event
            
            for ii = 1:length(obj.new_events)
                
                ev = obj.new_events(ii);
                
                ev.serial = length(obj.all_events) + ii; % just keep a running index
                
                % fluxes and timestamps
                ev.is_positive = obj.bank.fluxes_filtered(ev.time_index, ev.kern_index, ev.star_index)>0;
                
                ev.peak_timestamp = ev.timestamps(ev.time_index);
                ev.best_kernel = obj.bank.kernels(:,ev.kern_index);
                
                ev.threshold = obj.threshold;
                ev.used_background_sub = input.used_background_sub;
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
                
                b = vertcat(obj.prev_backgrounds, input.backgrounds);
                ev.backgrounds_at_peak = b(ev.frame_index, :);
                ev.backgrounds_at_star = b(:, ev.star_index);
                ev.backgrounds_time_average = mean(b, 1, 'omitnan');
                ev.backgrounds_star_average = mean(b, 2, 'omitnan');
                
                v = vertcat(obj.prev_variances, input.variances);
                ev.variances_at_peak = v(ev.frame_index, :);
                ev.variances_at_star = v(:, ev.star_index);
                ev.variances_time_average = mean(v, 1, 'omitnan');
                ev.variances_star_average = mean(v, 2, 'omitnan');
                
                dx = vertcat(obj.prev_offsets_x, input.offsets_x);
                ev.offsets_x_at_peak = dx(ev.frame_index, :);
                ev.offsets_x_at_star = dx(:, ev.star_index);
                ev.offsets_x_time_average = mean(dx, 1, 'omitnan');
                ev.offsets_x_star_average = mean(dx, 2, 'omitnan');
                
                dy = vertcat(obj.prev_offsets_y, input.offsets_y);
                ev.offsets_y_at_peak = dy(ev.frame_index, :);
                ev.offsets_y_at_star = dy(:, ev.star_index);
                ev.offsets_y_time_average = mean(dy, 1, 'omitnan');
                ev.offsets_y_star_average = mean(dy, 2, 'omitnan');
                
                w = vertcat(obj.prev_widths, input.widths);
                ev.widths_at_peak = w(ev.frame_index, :);
                ev.widths_at_star = w(:, ev.star_index);
                ev.widths_time_average = mean(w, 1, 'omitnan');
                ev.widths_star_average = mean(w, 2, 'omitnan');
                
                p = vertcat(obj.prev_bad_pixels, input.bad_pixels);
                ev.bad_pixels_at_peak = p(ev.frame_index, :);
                ev.bad_pixels_at_star = p(:, ev.star_index);
                ev.bad_pixels_time_average = mean(p, 1, 'omitnan');
                ev.bad_pixels_star_average = mean(p, 2, 'omitnan');
                
                ev.aperture = input.aperture;
                ev.gauss_sigma = input.gauss_sigma;
                
                ev.cutouts_first = obj.prev_cutouts(:,:,:,ev.star_index);
                ev.cutouts_second = input.cutouts(:,:,:,ev.star_index);
                ev.positions_first = obj.prev_positions;
                ev.positions_second = input.positions;
                ev.stack_first = obj.prev_stack;
                ev.stack_second = input.stack;
                ev.batch_index_first = obj.prev_batch_index;
                ev.batch_index_second = input.batch_index;
                ev.filename_first = obj.prev_filename;
                ev.filename_second = input.filename;
                if ischar(input.t_end) && isnumeric(input.t_end_stamp)
                    ev.t_end = input.t_end;
                    ev.t_end_stamp = input.t_end_stamp;
                end
                
                ev.phot_pars = obj.phot_pars;
                
                ev.max_num_nans = obj.max_num_nans;
                ev.max_corr = obj.max_corr;
                
                % any other info that needs to be saved along with the event object?
                % ...
                
            end
        
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
            
            if obj.debug_bit>1, fprintf('runtime "check_new_events": %f seconds\n', toc(t)); end
            
        end
        
        function finishup(obj) % remove black list events, clear images for unkept events
            
            t = tic;
            
            stars = [obj.all_events.star_index];
            if ~isempty(stars)
                [N,E] = histcounts(stars, 'BinWidth', 1, 'BinLimits', [1 max(stars)]);
                obj.black_list_stars = [obj.black_list_stars E(N>5)];
            end
            
            batches = [obj.all_events.batch_index];
            if ~isempty(batches)
                [N,E] = histcounts(batches, 'BinWidth', 1, 'BinLimits', [1 max(batches)]);
                obj.black_list_batches = [obj.black_list_batches E(N>5)];
            end
            
            for ii = 1:length(obj.all_events)
                
                if ismember(obj.all_events(ii).star_index, obj.black_list_stars)
                    obj.all_events(ii).keep = 0;
                    obj.all_events(ii).is_real = 0;
                    obj.all_events(ii).is_bad_star = 1;
                    obj.all_events(ii).addNote(sprintf('%s, star %d is on black list', obj.all_events(ii).notes, obj.all_events(ii).star_index));
                end
                
                if ismember(obj.all_events(ii).batch_index, obj.black_list_stars)
                    obj.all_events(ii).keep = 0;
                    obj.all_events(ii).is_real = 0;
                    obj.all_events(ii).is_bad_batch = 1;
                    obj.all_events(ii).addNote(sprintf('%s, batch %d is on black list', obj.all_events(ii).notes, obj.all_events(ii).batch_index));
                end
                
            end
            
            if obj.use_conserve_memory
                for ii = 1:length(obj.all_events)
                    if obj.all_events(ii).keep==0
                        obj.all_events(ii).clearImages;
                    end
                end
            end
            
            if obj.debug_bit>1, fprintf('runtime "finishup": %f seconds\n', toc(t)); end
            
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
                obj.display_event_idx = 1;
            end
            
            obj.display_event_idx = obj.display_event_idx - 1;
            
            if obj.display_event_idx<1
                obj.display_event_idx = obj.num_events;
            end
            
            if obj.use_display_kept_events
                
                for ii = obj.display_event_idx:-1:1
                    if obj.all_events(ii).keep
                        obj.display_event_idx = ii;
                        break;
                    end
                end
                
                % loop back around...
                for ii = obj.num_events:-1:obj.display_event_idx
                    if obj.all_events(ii).keep
                        obj.display_event_idx = ii;
                        break;
                    end
                end
                
            end
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.all_events(obj.display_event_idx).show('parent', obj.gui.panel_image);
            end
            
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
                obj.display_event_idx = 0;
            end
            
            obj.display_event_idx = obj.display_event_idx + 1;
            
            if obj.display_event_idx>obj.num_events
                obj.display_event_idx = 1;
            end
            
            if obj.use_display_kept_events
                
                for ii = obj.display_event_idx:obj.num_events
                    if obj.all_events(ii).keep
                        obj.display_event_idx = ii;
                        break;
                    end
                end
                
                % loop back around...
                for ii = 1:obj.display_event_idx
                    if obj.all_events(ii).keep
                        obj.display_event_idx = ii;
                        break;
                    end
                end
            
            end

            if ~isempty(obj.gui) && obj.gui.check
                obj.all_events(obj.display_event_idx).show('parent', obj.gui.panel_image);
            end
            
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
        
    end    
    
end

