classdef AstroData < util.oop.DynamicObject
% Contains all the outputs from image analysis (images, cutouts, etc.)
% Since the long and ever expanding list of data products is shared by
% several objects (file.Reader, file.BufferWheel, file.Deflator,
% img.PipeLine, etc.) it is helpful to keep these outputs all in one place.
% 

    properties % these include all the images, cutouts, lightcurves, etc that we use in image analysis
        
        % make sure anything you add is added to clear...
        images; % this is raw images and it is usually what we save on file
        images_proc; % this is after calibration, subtraction, etc. If no processing is done, put a pointer to the raw data. 
        
        cutouts; % this is raw cutouts and it is usually what we save on file
        cutouts_proc; % this is after calibration, subtraction, etc. If no processing is done, put a pointer to the raw data. 
        positions; % only for cutouts. a 2xN matrix (X then Y, N is the number of cutouts). 
        
        bg_cutouts; % samples of the raw images at random locations to calculate the backgrounds
        bg_cutouts_proc; % this is after calibration, subtraction, etc. If no processing is done, put a pointer to the raw data. 
        bg_positions; % only for bg_cutouts. a 2xN matrix (X then Y, N is the number of cutouts). 
        backgrounds; % value of the average number of photons per pixel (of the background)
        
        stack % sum of the full frame image
        stack_proc; % this is after calibration, subtraction, etc. If no processing is done, put a pointer to the raw data. 
        num_sum; % if the images are summed, how many frames were added. If equal to 1, the sum is the same as the images. 
        
        timestamps; % output timestamps (if available)
        t_start; % absolute date and time (UTC) when first image is taken
        t_end; % absolute date and time (UTC) when batch is finished
        t_end_stamp; % timestamp when batch is finished (hopefully, this is the same time that t_end is recorded). 
                
        psfs; % output PSFs from file (if available)
        sampling_psf; % if loaded PSFs need to be binned (this is the binning factor)
        
        lightcurves; % amount of light (total) in each cutout. 
        
    end
    
    methods % constructor
        
        function obj = AstroData(varargin)
        
            if ~isempty(varargin) && isa(varargin{1}, 'file.AstroData')
                obj = util.oop.full_copy(varargin{1});
            else
                
            end
            
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
            obj.images_proc = [];
            
            obj.cutouts = [];
            obj.cutouts_proc = [];
            obj.positions = [];

            obj.bg_cutouts = [];
            obj.bg_cutouts_proc = [];
            obj.bg_positions = [];
            obj.backgrounds = [];
            
            obj.stack = [];
            obj.num_sum = [];

            obj.timestamps = [];
            obj.t_start = [];
            obj.t_end = [];
            obj.t_end_stamp = [];

            obj.psfs = [];
            obj.sampling_psf = [];

            obj.lightcurves = [];
            
        end
                
    end
    
end