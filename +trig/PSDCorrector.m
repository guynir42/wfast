classdef PSDCorrector < handle
% Use Welch's method on a buffer of previous fluxes to get an estimate for 
% the Power Spectral Desnity (PSD) of the flux of each star. 
% You will need to keep quite a few samples of previous fluxes to get a good
% estimate of the PSD. This is controlled by "N_buf" and "B_min". 
% The former determines the maximum buffer size, and the latter determines
% the minimal buffer size that allows using the PSD. 
% In the first few batches the PSD is not used, until reaching "N_min". 
% Then after we reach "N_buf" samples in the buffer, the PSD correction is
% at full power, and remains so until the end of the run. Any new flux data
% that comes into this object is stored in a FIFO buffer. 
%
% Use input(flux) to give the newest batch of fluxes to this object. Make 
% sure to define "num_frames_to_add" as the number of flux points in a single
% batch, so it only inputs one batch at a time to the PSD calculation, not
% the overlapping two-batch we usually work with. The two-batch is indeed 
% given to input() but only so it can be calibrated (de-reddened) using the
% PSD we have already calculated on previous batches. 
%
% Use the "fluxes_deredened" for fluxes corrected by the sqrt of the PSD. 
% Use the "fluxes_blued" for fluxes with double correction (extra red removal) 
% which is useful for accounting also for the filter that should be divided 
% by sqrt of the PSD. Thus the "blued" fluxes can be match-filtered directly. 
% See Barak's paper: https://ui.adsabs.harvard.edu/abs/2019arXiv190805644Z/abstract
%
% NOTE: we use a linear fit (excluding outliers) to calibrate the fluxes
% before putting them into the PSD. If the flux buffer is too short for 
% getting a good PSD estimate (at the beginning of the run) then we return
% the fluxes with only this simple calibration. If the PSD is applied, it is
% applied to the calibrated fluxes. 
% 
% The PSD itself is stored in "psd" and can be plotted for inspection. 
% 
    
    
    properties(Transient=true)
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        fluxes_input;
        fluxes_calibrated; % removing a linear fit to the data
        stds_calibrated; % from the fit, not including outliers! 
        
        flux_buffer; % keep a FIFO buffer of fluxes from previos batches to run welch on 
        psd; % Power Spectra Density for each star (1st dim is freq, 2nd dim is star index)
        freq; % the frequency axis for the PSD
        
        fluxes_deredened;
        stds_deredened;
        
        fluxes_blued;
        stds_blued;
        
        
    end
    
    properties % switches/controls
        
        frame_rate = 25; 
        
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
                if obj.debug_bit>1, fprintf('PSDCorrector copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('PSDCorrector constructor v%4.2f\n', obj.version); end
            
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

                [obj.psd, obj.freq] = pwelch(obj.flux_buffer, size(obj.fluxes_input, 1), [], S(1), 'twosided', obj.frame_rate);

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

