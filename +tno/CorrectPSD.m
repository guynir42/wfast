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
       
        version = 1.00;
        
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
            
            flux_buffer(abs(flux_buffer-M)./S > 5) = NaN; % remove outliers
            
            obj.flux_buffer = flux_buffer - nanmean(flux_buffer);
            obj.flux_buffer = fillmissing(obj.flux_buffer, 'linear'); 
            
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

