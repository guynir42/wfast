classdef AstroData < dynamicprops
% Contains all the outputs from image analysis (images, cutouts, etc.)
% Since the long and ever expanding list of data products is shared by
% several objects (file.Reader, file.BufferWheel, file.Deflator,
% img.PipeLine, etc.) it is helpful to keep these outputs all in one place.
% 

    properties % these include all the images, cutouts, lightcurves, etc that we use in image analysis
        
        % make sure anything you add is added to clear...
        images; % this is raw images and it is usually what we save on file
        
        cutouts; % this is raw cutouts and it is usually what we save on file
        positions; % only for cutouts. a 2xN matrix (X then Y, N is the number of cutouts). 
        coordinates; % match each star/cutout position with RA/DEC (in degrees)
        magnitudes; % each star's magnitude (from catalog)
        temperature; % each star's temperature in K (from catalog)
        fluxes; % amount of light (total) in each cutout
        
        cutouts_bg; % samples of the raw images at random locations to calculate the backgrounds
        positions_bg; % only for bg_cutouts. a 2xN matrix (X then Y, N is the number of cutouts). 
        backgrounds; % value of the average number of photons per pixel (of the background)
        
        stack % sum of the full frame image
        num_sum; % if the images are summed, how many frames were added. If equal to 1, the sum is the same as the images. 
        
        timestamps; % output timestamps (if available)
        t_start; % absolute date and time (UTC) when first image is taken
        t_end; % absolute date and time (UTC) when batch is finished
        t_end_stamp; % timestamp when batch is finished (hopefully, this is the same time that t_end is recorded). 
                
        psfs; % output PSFs from file (if available)
        sampling_psf; % if loaded PSFs need to be binned (this is the binning factor)
        
        
    end
    
    properties(Access=private)
        
        positions_;
        positions_bg_;
        
    end
    
    methods % constructor
        
        function obj = AstroData(varargin)
        
            if ~isempty(varargin) && isa(varargin{1}, 'file.AstroData')
                obj = util.oop.full_copy(varargin{1});
            else
                
            end
            
        end
        
    end
    
    methods % getters
        
        function val = get.positions(obj)
            
            val = obj.getPositions;
            
        end
        
        function val = get.positions_bg(obj)
            
            val = obj.getPositionsBG;
            
        end
        
    end
    
    methods % setters
        
        function set.positions(obj, val)
            
            obj.setPositions(val);
            
        end
        
        function set.positions_bg(obj, val)
            
            obj.setPositionsBG(val);
            
        end
        
    end 
    
    methods % other utilities
        
        function takeFrom(obj, other)
            
            obj.copyFrom(other);
            other.clear;
            
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
                obj.(list{ii}) = other.(list{ii});
            end
            
        end
        
        function clear(obj)
            
            obj.images = [];
%             obj.images_proc = [];
            
            obj.cutouts = [];
%             obj.cutouts_proc = [];
            obj.positions = [];

            obj.cutouts_bg = [];
%             obj.bg_cutouts_proc = [];
            obj.positions_bg = [];
%             obj.backgrounds = [];
            
            obj.stack = [];
            obj.num_sum = [];

            obj.timestamps = [];
            obj.t_start = [];
            obj.t_end = [];
            obj.t_end_stamp = [];

            obj.psfs = [];
            obj.sampling_psf = [];

            obj.fluxes = [];
            
        end
                
    end
    
    methods(Access=protected)
        
        function val = getPositions(obj)
            
            val = obj.positions_;
            
        end
        
        function setPositions(obj, val)
            
            obj.positions_ = val;
            
        end
        
        function val = getPositionsBG(obj)
            
            val = obj.positions_bg_;
            
        end
        
        function setPositionsBG(obj, val)
            
            obj.positions_bg_ = val;
            
        end
        
    end
    
end