classdef AstroData < dynamicprops
% Contains all the outputs from image analysis (images, cutouts, etc.)
% Since the long and ever expanding list of data products is shared by
% several objects (file.Reader, file.BufferWheel, file.Deflator,
% img.PipeLine, etc.) it is helpful to keep these outputs all in one place.
% 

    properties % these include all the images, cutouts, lightcurves, etc that we use in image analysis
        
        % make sure anything you add is added to clear...
        images; % this is raw images and it is usually what we save on file
        
        timestamps; % output timestamps (if available)
        
        t_start; % absolute date and time (UTC) when first image is taken
        t_end; % absolute date and time (UTC) when batch is finished
        t_end_stamp; % timestamp when batch is finished (hopefully, this is the same time that t_end is recorded). 
        
        juldates; % translation of timestamps to julian date using t_end_stamp
        
        cutouts; % this is raw cutouts and it is usually what we save on file
        positions; % only for cutouts. a 2xN matrix (X then Y, N is the number of cutouts). 
        obj_idx; % what is the index of the star closest to the run's intended RA/DEC 
        coordinates; % match each star/cutout position with RA/DEC (in degrees)
        magnitudes; % each star's magnitude (from catalog)
        temperatures; % each star's temperature in K (from catalog)
        object_idx; % what is the index of the object closest to the coordinates given
        
        % these are photometric products for each batch
        fluxes; % amount of light (total) in each cutout
        errors; % error estimate on the flux measurements
        areas; % aperture area 
        backgrounds; % value of the average number of photons per pixel (of the background)
        variances; % value of the average variance per pixel (of the background)
        offsets_x; % measured inside the cutout (in pixels)
        offsets_y; % measured inside the cutout (in pixels)
        centroids_x; % measured relative to the full image (in pixels)
        centroids_y; % measured relative to the full image (in pixels)
        widths; % the second moment of the PSF
        bad_pixels; % how many bad pixels are included in the aperture
        flags; % if any of the measurements (1st or 2nd order) don't make sense, it is flagged by the photometry pipeline.  
        
        cutouts_bg; % samples of the raw images at random locations to calculate the backgrounds
        positions_bg; % only for bg_cutouts. a 2xN matrix (X then Y, N is the number of cutouts). 
        
        stack % sum of the full frame image
        num_sum; % if the images are summed, how many frames were added. If equal to 1, the sum is the same as the images. 
         
        psfs; % output PSFs from file (if available)
        sampling_psf; % if loaded PSFs need to be binned (this is the binning factor)
        
    end
    
%     properties(Access=private)
%         
%         positions_;
%         positions_bg_;
%         
%     end
    
    methods % constructor
        
        function obj = AstroData(varargin)
        
            if ~isempty(varargin) && isa(varargin{1}, 'file.AstroData')
%                 obj = util.oop.full_copy(varargin{1});
                
                % only save the properties that exist in AstroData
                % this means derived objects can be given to this
                % constructor and we get back an AstroData object.
                list = properties(file.AstroData);
                
                for ii = 1:length(list)
                    obj.(list{ii}) = varargin{1}.(list{ii});
                end

            else
                
            end
            
        end
        
    end
    
    methods % getters
%         
%         function val = get.positions(obj)
%             
%             val = obj.getPositions;
%             
%         end
%         
%         function val = get.positions_bg(obj)
%             
%             val = obj.getPositionsBG;
%             
%         end
        
    end
    
    methods % setters
        
%         function set.positions(obj, val)
%             
%             obj.setPositions(val);
%             
%         end
%         
%         function set.positions_bg(obj, val)
%             
%             obj.setPositionsBG(val);
%             
%         end
        
    end 
    
    methods % other utilities
        
        function takeFrom(obj, other)
            error('Do not use this!');
            obj.copyFrom(other);
            other.clear; % this is very dangerous, especially when taking from BufferWheel, as it resets the buffers as well... 
            
        end
        
        function copyTo(obj, other)
            
            copyFrom(other, obj);
            
        end
        
        function copyFrom(obj, other)
            
            if nargin<2 || isempty(other)
                error('must give a non-empty argument to "copyFrom"');
            end
            
            if ~isa(other, 'file.AstroData')
                error('Other object must be a subclass of "file.AstroData" but is a "%s" instead...', class(other));
            end
            
            list = properties(file.AstroData);
            
            for ii = 1:length(list)
                if ~isempty(other.(list{ii}))
                    obj.(list{ii}) = other.(list{ii});
                end
            end
            
        end
        
        function reset(obj)
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            obj.images = [];
            
            obj.timestamps = [];
            obj.t_start = [];
            obj.t_end = [];
            obj.t_end_stamp = [];
            obj.juldates = [];
            
            obj.cutouts = [];
            
            obj.positions = [];
            obj.coordinates = [];
            obj.magnitudes = [];
            obj.temperatures = [];
        
            obj.cutouts_bg = [];
            obj.positions_bg = [];
            
            obj.stack = [];
            obj.num_sum = [];
            
            obj.psfs = [];
            obj.sampling_psf = [];
            
            obj.fluxes = [];
            obj.errors = [];
            obj.areas = [];
            obj.backgrounds = [];
            obj.variances = [];
            obj.offsets_x = [];
            obj.offsets_y = [];
            obj.centroids_x = [];
            obj.centroids_y = [];
            obj.widths = [];
            obj.bad_pixels = [];
            obj.flags = [];
            
        end
                
    end
    
    methods(Access=protected)
%         
%         function val = getPositions(obj)
%             
%             val = obj.positions_;
%             
%         end
%         
%         function setPositions(obj, val)
%             
%             obj.positions_ = val;
%             
%         end
%         
%         function val = getPositionsBG(obj)
%             
%             val = obj.positions_bg_;
%             
%         end
%         
%         function setPositionsBG(obj, val)
%             
%             obj.positions_bg_ = val;
%             
%         end
%         
    end
    
end