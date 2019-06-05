classdef Finder < handle

    properties(Transient=true)
        
        generator@occult.CurveGenerator;
        
    end
    
    properties % objects
        
        cal@trig.Calibrator;
        filt@trig.Filter;
        ev@trig.Event;
        last_events@trig.Event;
        
    end
    
    properties % inputs/outputs
        
        prev_cutouts;
        prev_positions;
        prev_stack;
        prev_batch_index;
        prev_filename;
        
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
        
        function input(obj, fluxes, varargin)
            
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
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
            
            obj.cal.input(fluxes, input.timestamps); 
            obj.filt.input(obj.cal.fluxes_subtracted, input.timestamps); 
            
            obj.check_last_events; % make sure we aren't taking events that already exist
            
            for ii = 1:length(obj.last_events)
                
                obj.ev(end+1) = obj.last_events(ii);
                obj.ev(end).flux_raw_all = fluxes; 
                obj.ev(end).cutouts_first = obj.prev_cutouts;
                obj.ev(end).cutouts_second = input.cutouts;
                obj.ev(end).positions_first = obj.prev_positions;
                obj.ev(end).positions_second = input.positions;
                obj.ev(end).stack_first = obj.prev_stack;
                obj.ev(end).stack_second = input.stack;
                obj.ev(end).batch_index_first = obj.prev_batch_index;
                obj.ev(end).batch_index_second = input.batch_index;
                obj.ev(end).filename_first = obj.prev_filename;
                obj.ev(end).filename_second = input.filename;
                if ischar(input.t_end) && isnumeric(input.t_end_stamp)
                    obj.ev(end).t_end = input.t_end;
                    obj.ev(end).t_end_stamp = input.t_end_stamp;
                end
                
                % any other info that needs to be saved along with the event object? 
                % ...
                
            end
            
            % store these for next time
            obj.prev_cutouts = input.cutouts;
            obj.prev_positions = input.positions;
            obj.prev_stack = input.stack;
            obj.prev_batch_index = input.batch_index;
            obj.prev_filename = input.filename;
            
        end
        
        function check_last_events(obj)
            
            temp_events = obj.filt.found_events;
            
            idx_keep = true(length(temp_events),1); 
            
            for ii = 1:length(temp_events)
                if temp_events(ii).is_same(obj.last_events)
                    idx_keep(ii) = false;
                end
            end
            
            obj.last_events = temp_events(idx_keep);
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

