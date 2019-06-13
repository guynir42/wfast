classdef Finder < handle

    properties(Transient=true)
        
        generator@occult.CurveGenerator;
        
    end
    
    properties % objects
        
        cal@trig.Calibrator;
        filt@trig.Filter;
        ev@trig.Event;
        new_events@trig.Event;
        last_events@trig.Event;
        
        
    end
    
    properties % inputs/outputs
        
        black_list_stars;
        black_list_batches;
        
        prev_timestamps;
        prev_cutouts;
        prev_positions;
        prev_stack;
        prev_batch_index;
        prev_filename;
        
        prev_fluxes;
        prev_backgrounds;
        prev_variances;
        prev_offsets_x;
        prev_offsets_y;
        prev_widths;
        prev_bad_pixels;
        
    end
    
    properties % switches/controls
        
        threshold = 5; % threshold (in units of S/N) for peak of event 
        time_range_thresh = -2; % threshold for including area around peak (in continuous time)
        kern_range_thresh = -1; % area threshold (in kernels, discontinuous) NOTE: if negative this will be relative to "threshold"
        star_range_thresh = -1; % area threshold (in stars, discontinuous) NOTE: if this is higher than "threshold" there will be no area around peak
        
        max_events = 5; % how many events can we have triggered on the same 2-batch window?
        max_stars = 5; % how many stars can we afford to have triggered at the same time? 
        max_frames = 50; % maximum length of trigger area (very long events are disqualified)
        
        use_conserve_memory = 0;
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        ev_kept;
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Finder(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.Finder')
                if obj.debug_bit, fprintf('Finder copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Finder constructor v%4.2f\n', obj.version); end
            
                obj.cal = trig.Calibrator;
                obj.filt = trig.Filter;
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.ev = trig.Event.empty;
            obj.new_events = trig.Event.empty;
            obj.last_events = trig.Event.empty;
            
            obj.black_list_stars = [];
            obj.black_list_batches = [];
            
            obj.cal.reset;
            obj.filt.reset;
            
            obj.prev_fluxes = [];
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
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            obj.cal.clear;
            obj.filt.clear;
            
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
        
        function val = get.ev_kept(obj)
            
            val = obj.ev([obj.ev.keep]==1);
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function setupKernels(obj)
            
            if isempty(obj.generator)
                obj.generator = occult.CurveGenerator;
            end
            
            % some easy parameters
            obj.generator.R = 0;
            obj.generator.r = 1;
            obj.generator.b = 0;
            obj.generator.v = 5:5:30;
            
            obj.generator.getLightCurves;
            obj.filt.kernels = obj.generator.lc.flux - 1;
            
        end
        
        function input(obj, varargin)
                
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('fluxes', []);
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
            input.scan_vars(varargin{:});
            
            obj.clear;
            
            if ~isempty(obj.prev_fluxes) % skip first batch! 
            
                t = tic;
                obj.cal.input(vertcat(obj.prev_fluxes, input.fluxes), vertcat(obj.prev_timestamps, input.timestamps)); 
                if obj.debug_bit>1, fprintf('Calibration time: %f seconds.\n', toc(t)); end
                
                t = tic;
                obj.filt.input(obj.cal.fluxes_subtracted, obj.cal.timestamps); 
                if obj.debug_bit>1, fprintf('Filtering time: %f seconds.\n', toc(t)); end
                
                t = tic;
                obj.findEvents;
                if obj.debug_bit>1, fprintf('Triggering time: %f seconds.\n', toc(t)); end
                
                t = tic;
                obj.storeEventHousekeeping(input); % add some data from this batch to the triggered events
                
                obj.last_events = obj.new_events;
                obj.ev = [obj.ev obj.new_events];
                obj.new_events = trig.Event.empty;
            
                if obj.debug_bit>1, fprintf('Housekeeping time: %f seconds.\n', toc(t)); end
            
                t = tic;
                
                obj.checkEvents; % make sure we aren't taking events that already exist
            
                if obj.debug_bit>1, fprintf('Checking time: %f seconds.\n', toc(t)); end
                
            end
            
            % store these for next time
            obj.prev_fluxes = input.fluxes;
            obj.prev_backgrounds = input.backgrounds;
            obj.prev_variances = input.variances;
            obj.prev_offsets_x = input.offsets_x;
            obj.prev_offsets_y = input.offsets_y;
            obj.prev_widths = input.widths;
            obj.prev_bad_pixels = input.bad_pixels;
            obj.prev_timestamps = input.timestamps;
            obj.prev_cutouts = input.cutouts;
            obj.prev_positions = input.positions;
            obj.prev_stack = input.stack;
            obj.prev_batch_index = input.batch_index;
            obj.prev_filename = input.filename;
            
        end
        
        function findEvents(obj)
            
            obj.new_events = trig.Event.empty;
            
            ff = obj.filt.fluxes_filtered; % dim 1 is time, dim 2 is kernels, dim 3 is stars
            
            for ii = 1:obj.max_events
                
                [mx, idx] = util.stat.maxnd(abs(ff)); % note we are triggering on negative and positive events
                
                if mx<obj.threshold, break; end 
                
                ev = trig.Event;
                ev.snr = mx; % note this is positive even for negative filter responses! 
                
                ev.time_index = idx(1); 
                ev.kern_index = idx(2);
                ev.star_index = idx(3);
                
                ev.time_indices = obj.findTimeRange(ff, ev.time_index, ev.kern_index, ev.star_index); % find continuous area that is above time_range_thresh
                
                ev.kern_indices = find(max(abs(ff(ev.time_indices, :, ev.star_index)))>obj.getKernThresh);
                
                ev.star_indices = find(max(abs(ff(ev.time_indices, ev.kern_index, :)))>obj.getStarThresh);
                
                ev.timestamps = obj.filt.timestamps;
                ev.flux_filtered = obj.filt.fluxes_filtered(:,ev.kern_index,ev.star_index);
                ev.flux_raw_all = permute(obj.filt.fluxes, [1,3,2]);
                ev.stds_raw_all = permute(obj.filt.stds, [1,3,2]);
                
                ff(ev.time_indices, :, :) = 0; % don't look at the same region twice
                
                obj.new_events(end+1) = ev; % add this event to the list
                
            end
            
            
            
        end
        
        function time_range = findTimeRange(obj, ff, time_index, kern_index, star_index)

            N = size(ff,1); % time length
            
            thresh = obj.getTimeThresh;
            time_range = [];

            for jj = 0:N % go backward in time

                idx = time_index - jj;

                if idx<1, break; end

                if abs(ff(idx, kern_index, star_index))>=thresh
                    time_range = [time_range, idx];
                else 
                    break;
                end

            end

            time_range = flip(time_range);

            for jj = 1:N % go forward in time

                idx = time_index + jj;

                if idx>N, break; end

                if abs(ff(idx, kern_index, star_index))>=thresh
                    time_range = [time_range, idx];
                else 
                    break;
                end

            end

        end
        
        function storeEventHousekeeping(obj, input) % store all the additional metadata about this event
            
            for ii = 1:length(obj.new_events)
                
                ev = obj.new_events(ii);
                
                ev.serial = length(obj.ev) + ii; % just keep a running index
                
                % fluxes and timestamps
                ev.is_positive = obj.filt.fluxes_filtered(ev.time_index, ev.kern_index, ev.star_index)>0;
                
                ev.peak_timestamp = ev.timestamps(ev.time_index);
                ev.best_kernel = obj.filt.kernels(:,ev.kern_index);
                
                ev.threshold = obj.threshold;
                ev.time_range_thresh = obj.time_range_thresh;
                ev.kern_range_thresh = obj.kern_range_thresh;
                ev.star_range_thresh = obj.star_range_thresh;
                
                % housekeeping data from photometry
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
                
                ev.cutouts_first = obj.prev_cutouts;
                ev.cutouts_second = input.cutouts;
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
                            obj.new_events(ii).notes = [obj.new_events(ii).notes ', duplicate of ' num2str(obj.last_events(jj).serial)];
                        else % new event is better, mark off the old one
                            obj.last_events(jj).keep = 0;
                            obj.last_events(jj).is_duplicate = 1;
                            obj.last_events(jj).notes = [obj.last_events(jj).notes ', duplicate of ' num2str(obj.new_events(ii).serial)];
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
            
            stars = [obj.ev.star_index];
            if ~isempty(stars)
                [N,E] = histcounts(stars, 'BinWidth', 1, 'BinLimits', [1 max(stars)]);
                obj.black_list_stars = [obj.black_list_stars E(N>5)];
            end
            
            batches = [obj.ev.batch_index];
            if ~isempty(batches)
                [N,E] = histcounts(batches, 'BinWidth', 1, 'BinLimits', [1 max(batches)]);
                obj.black_list_batches = [obj.black_list_batches E(N>5)];
            end
            
            for ii = 1:length(obj.ev)
                
                if ismember(obj.ev(ii).star_index, obj.black_list_stars)
                    obj.ev(ii).keep = 0;
                    obj.ev(ii).is_real = 0;
                    obj.ev(ii).is_bad_star = 1;
                    obj.ev(ii).notes = sprintf('%s, star %d is on black list', obj.ev(ii).notes, obj.ev(ii).star_index);
                end
                
                if ismember(obj.ev(ii).batch_index, obj.black_list_stars)
                    obj.ev(ii).keep = 0;
                    obj.ev(ii).is_real = 0;
                    obj.ev(ii).is_bad_batch = 1;
                    obj.ev(ii).notes = sprintf('%s, batch %d is on black list', obj.ev(ii).notes, obj.ev(ii).batch_index);
                end
                
            end
            
            if obj.use_conserve_memory
                for ii = 1:length(obj.ev)
                    if obj.ev(ii).keep==0
                        obj.ev(ii).clearImages;
                    end
                end
            end
            
            if obj.debug_bit>1, fprintf('runtime "finishup": %f seconds\n', toc(t)); end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

