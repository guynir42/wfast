classdef MicroFlare < handle

    properties(Transient=true)
        gui;
    end
    
    properties % inputs/outputs
        
        filename; % where the flare was detected
        file_index; % index in the run
        frame_index; % what frame did the peak appear in
        pos; % x and y position of the peak pixel in the full-frame image
        peak; % value of the peak pixel
        pixel_var; % variance of underlying pixel (from calibration file)
        cutouts; % cutout images around the peak and all frames before/after it in that batch
        
        % from photometry2 or from rough estimates on the whole cutout
        flux;
        error;
        background;
        area;
        offset;        
        width;
        bad_pixels;
        
        mean; % mean flux calculated outside the peak region
        std; % RMS flux calculated outside the peak region
        
        image; % single cutout at the frame of the peak
        num_peaks; % how many distinct peaks can we find in the one image
        num_frames; % how many continuous frames we have around the peak
        num_pixels; % how many pixels above the threshold do we have 
        
        fwhm; % calculted from the image, using dedicated FWHM function (e.g., using filters)
        
        type = ''; % can be cosmic ray, flare, satellite, pixel... 
        
    end
    
    properties % switches/controls
        
        pixel_thresh = 256; % threshold of the brightest pixel, used for triggering 
        flux_thresh = 7.5; % threshold for peak intensity (of the flux, not of single pixel)
        side_thresh = 3; % threshold for flux in frames adjacent to the peak
        interval_peak = 2; % how many frames in either direction from peak should be excluded for calculating the mean/std of the flux
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = MicroFlare(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'img.MicroFlare')
                if obj.debug_bit>1, fprintf('MicroFlare copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('MicroFlare constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods(Static=true) % utilities
        
        function obj = readStruct( st)
            
            if ~isstruct(st)
                error('Must supply a struct array to this function! Instead got a "%s" type object.', class(st)); 
            end
            
            list = properties(img.MicroFlare); 
            
            for ii = 1:length(st)
                obj(ii) = img.MicroFlare; 
                for jj = 1:length(list)
                    if isfield(st(ii), list{jj})
                        obj(ii).(list{jj}) = st(ii).(list{jj}); 
                    end
                end
            end
            
        end
        
    end
    
    methods % calculations
        
        function calculate(obj, varargin)
            
            if numel(obj)>1
                
                for ii = 1:numel(obj)
                    obj(ii).calculate(varargin{:}); % recursively call this for individual objects in the array
                end
                
                return;
                
            end
            
            input = util.text.InputVars;
            input.input_var('pixel_thresh', [], 'pixel_threshold'); 
            input.input_var('flux_thresh', [], 'flux_threshold'); 
            input.input_var('side_thresh', [], 'side_threshold'); 
            input.input_var('interval_peak', []); 
            input.input_var('debug_bit', []); 
            input.scan_vars(varargin{:}); 
            
            for ii = 1:length(input.list_added_properties)
                if ~isempty(input.(input.list_added_properties{ii})) && isprop(obj, input.list_added_properties{ii})
                    obj.(input.list_added_properties{ii}) = input.(input.list_added_properties{ii}); 
                end
            end
            
            obj.image = obj.cutouts(:,:,obj.frame_index); 
            
            indices = obj.frame_index;
            for ii = 1:obj.interval_peak
                
                idx = obj.frame_index - ii;
                if idx<1, break; end
                
                indices = [indices idx]; 
                    
            end
            
            for ii = 1:obj.interval_peak
                
                idx = obj.frame_index + ii;
                if idx>length(obj.flux)
                    break;
                end
                
                indices = [indices idx]; 
                    
            end
            
            indices = sort(indices);
            flux_bg = obj.flux;
            flux_bg(indices) = NaN;
            
            obj.mean = nanmean(flux_bg);
            obj.std = nanstd(flux_bg); 
            
            flux_norm = (obj.flux-obj.mean)./obj.std;
            
            % how many frames above the threshold
            counter = 1;
            for ii = 1:length(flux_norm)
                idx = obj.frame_index - ii;
                
                if idx<1, break; end
                if flux_norm(idx)<obj.side_thresh, break; end
                
                counter = counter + 1;
                
            end
            
            for ii = 1:length(flux_norm)
                idx = obj.frame_index + ii;
                
                if idx>length(flux_norm), break; end
                if flux_norm(idx)<obj.side_thresh, break; end
                
                counter = counter + 1;
                
            end
            
            obj.num_frames = counter;
            
            % how many pixels in the image are above the threshold
            BW = obj.image>obj.pixel_thresh; % black/white image of pixels above the threshold
            N = nnz(BW); 
            
            if N<=1
                obj.num_pixels = 1;
            else
                
                C = bwconncomp(BW); 
                [X,Y] = meshgrid((1:size(obj.image,2))-floor(size(obj.image,2)/2)-1, (1:size(obj.image,1))-floor(size(obj.image,1)/2)-1);
                
                for ii = 1:C.NumObjects
                    
                    r(ii) = min(sqrt(X(C.PixelIdxList{ii}).^2 + Y(C.PixelIdxList{ii}).^2)); 
                    
                end
                
                [~,idx] = min(r); % find the connected component that is nearest to the center of the cutout
                
                N = numel(C.PixelIdxList{idx}); 
                
            end
            
            obj.num_pixels = N;
            
            obj.fwhm = util.img.fwhm(obj.image, 'method', 'filters'); 
            
            if obj.num_pixels==1
                obj.type = 'bad pixel';
            elseif obj.num_frames==1
                obj.type = 'cosmic ray';
            else
                obj.type = 'flare';
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

