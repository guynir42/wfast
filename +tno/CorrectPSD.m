classdef CorrectPSD < handle
% Calculate the Power Spectral Density of lightcurves, and apply a correction
% to new lightcurves that are given to it. 
% 
% Use calcPSD() by giving it a long baseline flux matrix and timestamps. 
% Use input() with the new, short flux interval to get corrected fluxes. 
% Handles 2D matrices where each column is a lightcurve for a different star. 
% 
% The object parameters, "window_size", "overlap" and "num_points" are passed
% to the pwelch() function used to calculate the PSD. 
% 
    
    
    properties(Transient=true)
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        timestamps; 
        flux_buffer;
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
            
            obj.timestamps = [];
            obj.flux_buffer = [];
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
            
            if nargin<3 || isempty(psd_power)
                psd_power = 1; 
            end
            
            if nargin<4 || isempty(star_indices)
                star_indices = 1:size(flux,2);
            end
            
            flux = flux - nanmean(flux); 
            
            S = size(flux); 
            ff = fft(util.img.pad2size(flux, [obj.num_points S(2:end)])); % zero pad the flux, then FFT it
            
            ff_corrected = ff./(obj.power_spectrum(:,star_indices)).^psd_power;
            
            flux_corrected = util.img.crop2size(real(ifft(ff_corrected)), S); 
            
        end
        
        function calcPSD(obj, flux_buffer, timestamps, window, number)
            
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
        
        function f_sim = inpaintNaNs(obj, star_indices)
        % generate a simulated flux from a randomized instance of the PSD 
        % FFTd back to time-domain, and use that to in-paint any NaN values
        % in the flux buffer.
            
            if nargin<2 || isempty(star_indices)
                star_indices = 1:size(obj.flux_buffer,2);
            end
            
            f_real = obj.flux_buffer(:,star_indices); 
            s_real = nanstd(f_real,[],1); 
            
            P = obj.power_spectrum(:,star_indices); 
            f_sim = real(fft(sqrt(P).*exp(1i*pi*rand(size(P))))); 
            s_sim = nanstd(f_sim,[],1); 
            
            f_sim = f_sim.*s_real./s_sim; % normalize by the overall noise level
            
            f_inpaint = zeros(size(f_real), 'like', f_real); 
            
            start_index = size(f_real,1) - size(f_sim,1) + 1;
            
            if start_index<1 % the simulated data is bigger than the flux buffer! 
                f_sim = f_sim(2-start_index:end,:); 
                start_index = 1; 
            end
            
            f_inpaint(start_index:end,:) = f_sim; 
            
            frame_indices = isnan(f_real);
            f_real(frame_indices) = f_inpaint(frame_indices); % replace all NaNs inside the selected star LCs
            
            obj.flux_buffer(:,star_indices) = f_real; % return the adjusted LCs into the buffer
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

