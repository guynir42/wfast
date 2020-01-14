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
        stds_deredened;
        
        fluxes_blued;
        stds_blued;
        
        
    end
    
    properties % switches/controls
        
        num_frames_to_add = 100; % how many frames from the input flux should be added to the PSD (in case we get multiple, overlapping batches)
        
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
            obj.stds_deredened = [];
            
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
            
            obj.fluxes_raw = vertcat(obj.fluxes_raw, fluxes(1:obj.num_frames_to_add,:));
            
            if size(obj.fluxes_raw,1)>obj.N
                obj.fluxes_raw = obj.fluxes_raw(end-obj.N+1:end,:,:); % be careful with fluxes with more than 3 dimensions! 
            end
            
            obj.calculate;
            
        end
        
        function calculate(obj)
           
            f = obj.fluxes_input;
            f = util.img.pad2size(f, [2 1].*size(f)); % zero pad the fluxes! 
            
            obj.psd = pwelch(obj.fluxes_raw, 200, [], size(obj.fluxes_input,1)*2, 'twosided');
                        
            obj.fluxes_deredened = util.img.crop2size(ifft(fft(f./sqrt(obj.psd))), size(obj.fluxes_input)); 
            obj.stds_deredened = std(obj.fluxes_deredened); 
            
            % these are divided twice by the sqrt(PSD) to account for the filter being deredened as well.
            obj.fluxes_blued = util.img.crop2size(ifft(fft(f./(obj.psd))), size(obj.fluxes_input)); 
            obj.stds_blued = std(obj.fluxes_blued); 
            
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

