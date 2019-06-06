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
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
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
            obj.cal.reset;
            obj.filt.reset;
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            obj.cal.clear;
            obj.filt.clear;
            
        end
        
    end
    
    methods % getters
        
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
            
                obj.cal.input(vertcat(obj.prev_fluxes, input.fluxes), vertcat(obj.prev_timestamps, input.timestamps)); 
                obj.filt.input(obj.cal.fluxes_subtracted, obj.cal.timestamps); 
                
                obj.new_events = obj.new_events;
                
                for ii = 1:length(obj.new_events)
                    
                    this_ev = obj.new_events(ii);
                    this_ev.serial = length(obj.ev) + ii;
                    this_ev.flux_raw_all = fluxes; 
                    this_ev.cutouts_first = obj.prev_cutouts;
                    this_ev.cutouts_second = input.cutouts;
                    this_ev.positions_first = obj.prev_positions;
                    this_ev.positions_second = input.positions;
                    this_ev.stack_first = obj.prev_stack;
                    this_ev.stack_second = input.stack;
                    this_ev.batch_index_first = obj.prev_batch_index;
                    this_ev.batch_index_second = input.batch_index;
                    this_ev.filename_first = obj.prev_filename;
                    this_ev.filename_second = input.filename;
                    if ischar(input.t_end) && isnumeric(input.t_end_stamp)
                        this_ev.t_end = input.t_end;
                        this_ev.t_end_stamp = input.t_end_stamp;
                    end

                    b = vertcat(obj.prev_backgrounds, input.backgrounds);
                    this_ev.background_at_peak = b(this_ev.which_frame, this_ev.which_star);
                    this_ev.background_time_average = mean(b, 1, 'omitnan'); 
                    this_ev.background_space_average = mean(b, 2, 'omitnan'); 
                    
                    v = vertcat(obj.prev_variances, input.variances);
                    this_ev.variance_at_peak = v(this_ev.which_frame, this_ev.which_star);
                    this_ev.variance_time_average = mean(v, 1, 'omitnan'); 
                    this_ev.variance_space_average = mean(v, 2, 'omitnan'); 
                    
                    dx = vertcat(obj.prev_offsets_x, input.offsets_x);
                    this_ev.offset_x.at_peak = dx(this_ev.which_frame, this_ev.which_star);
                    this_ev.offset_x_time_average = mean(dx, 1, 'omitnan'); 
                    this_ev.offset_x_space_average = mean(dx, 2, 'omitnan'); 
                    
                    dy = vertcat(obj.prev_offsets_y, input.offsets_y);
                    this_ev.offset_y_at_peak = dy(obj.which_frame, obj.which_star);
                    this_ev.offset_y_time_average = mean(dy, 1, 'omitnan'); 
                    this_ev.offset_y_space_average = mean(dy, 2, 'omitnan'); 
                    
                    w = vertcat(obj.prev_widths, input.widths);
                    this_ev.width_at_peak = w(this_ev.which_frame, this_ev.which_star);
                    this_ev.width_time_average = mean(w, 1, 'omitnan'); 
                    this_ev.width_space_average = mean(w, 2, 'omitnan'); 
                    
                    p = vertcat(obj.prev_bad_pixels, input.bad_pixels);
                    this_ev.bad_pixels_at_peak = p(this_ev.which_frame, this_ev.which_star); 
                    this_ev.bad_pixels_time_average = mean(p, 1, 'omitnan'); 
                    this_ev.bad_pixels_space_average = mean(p, 2, 'omitnan'); 
                    
                    % any other info that needs to be saved along with the event object? 
                    % ...

                end
                
                obj.check_new_events; % make sure we aren't taking events that already exist
                
                obj.last_events = obj.new_events;
                obj.ev = [obj.ev obj.new_events];
            
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
        
        function check_new_events(obj)
            
            % go over new events and check if they are duplicates
            for ii = 1:length(obj.new_events)
                
                for jj = 1:length(obj.last_events)
                    
                    if obj.new_events(ii).is_same(obj.last_events(jj))
                        
                        if obj.new_events(ii).snr<obj.last_events(jj).snr % old event is better, mark off the new one
                            obj.new_events(ii).keep = 0;
                            obj.new_events(ii).notes = [obj.new_events(ii).notes ', duplicated of ' num2str(obj.last_events(jj).serial)];
                        else % new event is better, mark off the old one
                            obj.last_events(ii).keep = 0;
                            obj.last_events(ii).notes = [obj.last_events(ii).notes ', duplicated of ' num2str(obj.new_events(jj).serial)];
                        end
                        
                        break;
                        
                    end
                    
                end
                
                obj.new_events(ii).self_check;
                
            end
            
            obj.last_events = obj.new_events; % this would be added to the list of found events
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

