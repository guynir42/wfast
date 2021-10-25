classdef CorrectPSD < handle
% Calculate the Power Spectral Density of lightcurves, and apply a correction
% to new lightcurves that are given to it. 
% 
% Use addToBuffer() with the detrended fluxes (and buffer length) to update
% the flux buffer that is used to calculate the PSD. 
% Use calcPSD() to produce an updated PSD estimate. 
% Use input() with the new, short flux interval to get corrected fluxes. 
% All fluxes are expected to be detrended (i.e., removed a linear fit). 
% Handles 2D matrices where each column is a lightcurve for a different star. 
% 
% The object parameters, "window_size", "overlap" and "num_points" are passed
% to the pwelch() function used to calculate the PSD. 
% If the size of the buffer is smaller than the window_size, we use regular
% FFT to get the power spectrum. If it is longer, we apply Welch's method. 
    
    
    properties(Transient=true)
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        flux_buffer;
        frame_rate;
        power_spectrum; 
        freq; 
                
    end
    
    properties % switches/controls

        window_size = 100;
        overlap = []; % take the default value
        num_points = 200; 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.03;
        
    end
    
    methods % constructor
        
        function obj = CorrectPSD(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'tno.CorrectPSD')
                if obj.debug_bit>1, fprintf('CorrectPSD copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('CorrectPSD constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.flux_buffer = [];
            obj.frame_rate = [];
            obj.power_spectrum = [];
            obj.freq = [];
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            
            
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function flux_corrected = input(obj, flux, psd_power, star_indices)
        % Give a short interval of fluxes to be corrected. 
        % INPUTS: 
        %   -flux: must be detrended. Can be 2D. 
        %   -psd_power: what power of the PSD to use when dividing in 
        %               Fourier space (default is 1, but also 0.5 is
        %               sometimes useful). Use 1 if you plan to match
        %               filter the data. 
        %   -star_indices: run the correction only on a subset of the stars. 
        
            if nargin<3 || isempty(psd_power)
                psd_power = 1; 
            end
            
            if nargin<4 || isempty(star_indices)
                star_indices = 1:size(flux,2);
            end
            
            flux = flux - nanmean(flux); 
            
            S = size(flux); 
            ff = fft(util.img.pad2size(flux, [size(obj.power_spectrum,1) S(2:end)])); % zero pad the flux, then FFT it
            
            ff_corrected = ff./(obj.power_spectrum(:,star_indices)).^psd_power;
            
            flux_corrected = util.img.crop2size(real(ifft(ff_corrected)), S); 
            
        end
        
        function addToBuffer(obj, flux, buffer_size)
        % Add new fluxes (detrended!) to the buffer used to produce the PSD. 
        % The buffer_size determines the maximum size of the buffer and
        % when new data increase its size too much, it tosses out the
        % oldest data (FIFO). 
        % Any NaN values in the fluxes are removed at this stage. 
        
            M = nanmedian(flux); 
            S = mad(flux); 
            
            outlier_idx = abs(flux-M)./S > 5;
            flux(outlier_idx) = NaN; % remove outliers
            
            flux = flux - nanmean(flux); % remove the mean value after removing outliers
            
            if isempty(obj.power_spectrum)
                flux(isnan(flux)) = 0; 
            else
                nan_indices = find(sum(isnan(obj.flux_buffer),1));
                flux = obj.inpaintNaNs(flux, nan_indices);
                flux(isnan(flux)) = 0; % just to verify all NaNs were eliminated
            end
            
            obj.flux_buffer = vertcat(obj.flux_buffer, flux); % add the new fluxes
            
            if size(obj.flux_buffer,1)>buffer_size
                obj.flux_buffer = obj.flux_buffer(end-buffer_size+1:end,:); % remove old fluxes from start of buffer
            end
            
        end
        
        function flux_painted = inpaintNaNs(obj, flux, star_indices)
        % generate a simulated flux from a randomized instance of the PSD 
        % FFTd back to time-domain, and use that to in-paint any NaN values
        % in the flux buffer.
            
            if nargin<2 || isempty(star_indices)
                star_indices = 1:size(flux,2);
            end
            
            f_real = flux(:,star_indices);
            s_real = nanstd(obj.flux_buffer(:,star_indices),[],1); % the RMS of the real fluxes 
            
            P = obj.power_spectrum(:,star_indices); 
            f_sim = real(fft(sqrt(P).*exp(1i*2*pi*rand(size(P))))); % draw random phases from the power spectrum to make simulated flux
            s_sim = nanstd(f_sim,[],1); 
            
            f_sim = f_sim.*s_real./s_sim; % normalize by the overall noise level
            
            start_index = size(f_sim,1) - size(f_real,1) + 1;
            
            f_sim = f_sim(start_index:end,:); 
            
            frame_indices = isnan(f_real);
            f_real(frame_indices) = f_sim(frame_indices); % replace all NaNs inside the selected star LCs
            
            flux_painted = flux; 
            flux_painted(:,star_indices) = f_real; 
            
        end
        
        function calcPSD(obj)
            
            if size(obj.flux_buffer,1)<obj.window_size
                obj.power_spectrum = abs(fft(obj.flux_buffer).^2); 
            else
                [obj.power_spectrum, obj.freq] = pwelch(obj.flux_buffer, obj.window_size, obj.overlap, obj.num_points, obj.frame_rate, 'twosided');
            end
            
        end
        
        function calcPSD_old(obj, flux_buffer, timestamps, window, number)
            
            M = nanmedian(flux_buffer); 
            S = mad(flux_buffer); 
            
            outlier_idx = abs(flux_buffer-M)./S > 5;
            flux_buffer(outlier_idx) = NaN; % remove outliers
            
            obj.flux_buffer = flux_buffer - nanmean(flux_buffer);
%             obj.flux_buffer(outlier_idx) = 0; % remove these NaNs after subtracting the mean

            if isempty(obj.power_spectrum)
                obj.flux_buffer = fillmissing(obj.flux_buffer, 'linear'); % this in unecessary if we already remove NaNs at input time
            else
                nan_indices = find(sum(isnan(obj.flux_buffer),1)); 
                obj.inpaintNaNs(nan_indices); 
            end
            
            obj.timestamps = timestamps;
            
            if nargin>3 && ~isempty(window)
                obj.window_size = window;
            end
            
            if nargin>4 && ~isempty(number)
                obj.num_points = number;
            end
            
            dt = median(diff(obj.timestamps)); 
            frame_rate = 1./dt; 
            
            [obj.power_spectrum, obj.freq] = pwelch(obj.flux_buffer, obj.window_size, obj.overlap, obj.num_points, frame_rate, 'twosided');
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

