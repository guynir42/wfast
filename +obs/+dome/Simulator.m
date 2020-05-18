classdef Simulator < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        status = 1;
        id = 'dome(sim)';
        
        shutter1_deg = 90;
        shutter2_deg = 90;
        
        is_closed = 1;
        
        
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
        
        function obj = Simulator(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.dome.Simulator')
                if obj.debug_bit>1, fprintf('Simulator (dome) copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('Simulator (dome) constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

