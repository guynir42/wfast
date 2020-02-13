classdef PSDCorrector < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        fluxes_input;
        fluxes_calibrated; % removing a linear fit to the data
        stds_calibrated; % from the fit, not including outliers! 
        
        flux_buffer; 
        psd; 
        
        fluxes_deredened;
        stds_deredened;
        
        fluxes_blued;
        stds_blued;
        
        
    end
    
    properties % switches/controls
        
        num_frames_to_add = 100; % how many frames from the input flux should be added to the PSD (in case we get multiple, overlapping batches)
        
        N_buf = 20000; % number of samples in the buffer
        N_min = 5000; % number of samples needed to apply PSD corrections
        
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
            
            obj.flux_buffer = [];
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
            
            obj.fluxes_input = fillmissing(fluxes, 'spline'); 
            
            fr = util.fit.polyfit(1:size(obj.fluxes_input), obj.fluxes_input, 'order', 1); 

            obj.fluxes_calibrated = obj.fluxes_input-[fr.ym]; % subtract the model 
            obj.stds_calibrated = sqrt([fr.variance]); % the model variance is without outliers...  
            
            obj.calculate;
            
            obj.flux_buffer = vertcat(obj.flux_buffer, obj.fluxes_calibrated(1:obj.num_frames_to_add,:));
            
            if size(obj.flux_buffer,1)>obj.N_buf
                obj.flux_buffer = obj.flux_buffer(end-obj.N_buf+1:end,:,:); % be careful with fluxes with more than 3 dimensions! 
            end
            
        end
        
        function calculate(obj)
            
            if isempty(obj.N_min) || size(obj.flux_buffer,1)>=obj.N_min

                f = obj.fluxes_calibrated;

                S = size(f);
                S(1) = S(1)*2; % add padding twice the size of the data
                f = util.img.pad2size(f, S); % zero pad the fluxes! 

                obj.psd = pwelch(obj.flux_buffer, size(obj.fluxes_input, 1), [], S(1), 'twosided');

                ff = fft(f);

                ff(1,:) = 0; % zero frequency contains just noise, and should be zero after detrending with linear fitter
                
                obj.fluxes_deredened = real(util.img.crop2size(ifft(ff./sqrt(obj.psd)), size(obj.fluxes_input))); 
                obj.stds_deredened = std(obj.fluxes_deredened); 

                % these are divided twice by the sqrt(PSD) to account for the filter being deredened as well.
                obj.fluxes_blued = real(util.img.crop2size(ifft(ff./obj.psd), size(obj.fluxes_input))); 
                obj.stds_blued = std(obj.fluxes_blued); 

            else % the flux buffer is not large enough to get a reasonable estimate for the PSD, just detrend them using polyfit
                
                obj.fluxes_deredened = obj.fluxes_calibrated; % subtract the model 
                obj.stds_deredened = obj.stds_calibrated;
                
                obj.fluxes_blued = obj.fluxes_calibrated;
                obj.stds_blued = obj.stds_calibrated; 
                
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

