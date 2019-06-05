classdef Finder < handle

    properties(Transient=true)
        
        generator@occult.CurveGenerator;
        
    end
    
    properties % objects
        
        cal@trig.Calibrator;
        filt@trig.Filter;
        ev@trig.Events;
        
    end
    
    properties % inputs/outputs
        
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
            
            obj.filt.kernels = obj.generator.getLightCurves - 1;
            
        end
        
        function input(obj, fluxes, timestamps)
            
            obj.cal.input(fluxes, timestamps); 
            obj.filt.input(obj.cal.fluxes_subtracted, timestamps); 
            
            for ii = 1:length(obj.filt.found_events)
                obj.ev(end+1) = obj.filt.found_events(ii);
                obj.ev(end).flux_raw_full = fluxes; 
                % any other info that needs to be saved along with the event object? 
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

