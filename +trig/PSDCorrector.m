classdef PSDCorrector < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        fluxes_input;
        fluxes_raw; 
        psd; 
        
        fluxes_deredened;
        std_deredened;
        
        fluxes_blued;
        stds_blued;
        
        
    end
    
    properties % switches/controls
        
        N = 1000; % number of samples in the buffer
        
                
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = PSDCorrector(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.PSDCorrector')
                if obj.debug_bit, fprintf('PSDCorrector copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('PSDCorrector constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.fluxes_raw = [];
            obj.psd = [];
            
            obj.clear;
            
        end
        
        function clear(obj)
           
            obj.fluxes_deredened = [];
            obj.std_deredened = [];
            
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function input(obj, fluxes)
            
            if ndims(fluxes)>3
                error('Must add some more : marks to allow fluxes with over 3 dimensions!'); 
            end
            
            obj.clear;
            
            fluxes = fillmissing(fluxes, 'spline'); 
            
            obj.fluxes_input = fluxes;
            
            obj.fluxes_raw = vertcat(obj.fluxes_raw, fluxes);
            
            if size(obj.fluxes_raw,1)>obj.N
                obj.fluxes_raw = obj.fluxes_raw(end-obj.N+1:end,:,:); % be careful with fluxes with more than 3 dimensions! 
            end
            
            obj.calculate;
            
        end
        
        function calculate(obj)
           
            obj.psd = pwelch(obj.fluxes_raw, 128, [], size(obj.fluxes_input,1), 'twosided')*2; % this factor of 2 needs to be explained at some point! 
            
            obj.fluxes_deredened = ifft(fft(obj.fluxes_input./sqrt(obj.psd))); 
            
            obj.std_deredened = std(obj.fluxes_deredened); 
            
            % these are divided twice by the sqrt(PSD) to account for the filter being deredened as well.
            obj.fluxes_blued = ifft(fft(obj.fluxes_input./(obj.psd))); 
            
            obj.std_blued = std(obj.fluxes_deredened); 
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

