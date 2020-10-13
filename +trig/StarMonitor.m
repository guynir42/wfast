classdef StarMonitor < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        cand@trig.Candidate;
        
    end
    
    properties % inputs/outputs
        
        flux; 
        aux; 
        
    end
    
    properties % switches/controls
        
        star_indices = []; 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = StarMonitor(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.StarMonitor')
                if obj.debug_bit>1, fprintf('StarMonitor copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('StarMonitor constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.cand = trig.Candidate.empty;
            
            obj.flux = []; 
            obj.aux = []; 
            
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function input(obj, flux, aux, cutouts, candidates)
            
            if isempty(obj.star_indices)
                return;
            end
            
            obj.flux = vertcat(obj.flux, flux(:,obj.star_indices)); 
            obj.aux = vertcat(obj.aux, aux(:,obj.star_indices,:)); 
            
            for ii = 1:length(obj.star_indices)
                for jj = 1:length(candidates)
                    if candidates(jj).star_index==obj.star_indices(ii)
                        c = util.oop.full_copy(candidates(jj));
                        c.cutouts_all = cutouts;
                        obj.cand = vertcat(obj.cand, c);
                    end
                end
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

