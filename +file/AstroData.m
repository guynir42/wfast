classdef AstroData < dynamicprops
% Contains all the outputs from image analysis (images, cutouts, etc.)
% Since the long and ever expanding list of data products is shared by
% several objects (file.Reader, file.BufferWheel, file.Deflator,
% img.PipeLine, etc.) it is helpful to keep these outputs all in one place.
% 
% These properties are all meant to contain data for a single batch, 
% so it would have one image in the "stack", 100 images and timestamps, etc. 
%
% NOTE: this object and all its sub-clasess (inheritors) get the dynamicprops
%       ability to add properties when needed. This is useful when doing 
%       custom analysis and need temporary variables to be stored with the 
%       regular data products. 

    properties % these include all the images, cutouts, lightcurves, etc that we use in image analysis
        
        % make sure anything you add is added to clear...
        images; % this is raw images and it is usually what we save on file
        
        timestamps; % output timestamps (if available)
        t_start; % absolute date and time (UTC) when first image is taken
        t_end; % absolute date and time (UTC) when batch is finished
        t_end_stamp; % timestamp when batch is finished (hopefully, this is the same time that t_end is recorded). 
         
        juldates; % translation of timestamps to julian date using t_end_stamp
        
        stack % sum of the full frame image
        num_sum; % if the images are summed, how many frames were added. If equal to 1, the sum is the same as the images. 
        
        cutouts; % this is raw cutouts and it is usually what we save on file
        positions; % only for cutouts. a 2xN matrix (X then Y, N is the number of cutouts). 
        object_idx; % what is the index of the object closest to the coordinates given
        forced_indices; % if any forced targets are added to the positions field, list their indices in a vector
        unlocked_indices; % if any targets are set to adjust cutout position independently (not locked to the mean drift), list their indices in a vector
        dynamic_indices; % these positions can change with each batch, e.g., for tracking cosmic rays 
        
        coordinates; % match each star/cutout position with RA/DEC (in degrees)
        magnitudes; % each star's magnitude (from catalog)
        temperatures; % each star's temperature in K (from catalog)
        
        cutouts_bg; % samples of the raw images at random locations to calculate the backgrounds
        positions_bg; % only for bg_cutouts. a 2xN matrix (X then Y, N is the number of cutouts). 
        
        psfs; % output PSFs from file (if available)
        sampling_psf; % if loaded PSFs need to be binned (this is the binning factor)
        
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
        
        
    end
    
    properties(Hidden=true)
    
        list_dynamic_props = {}; 
        list_permanent_props = {'list_permanent_props', 'list_dynamic_props', ...
            'positions', 'positions_bg', 'object_idx', ... 
            'forced_indices', 'unlocked_indices', 'dynamic_indices', ...
            'coordinates', 'magnitudes', 'temperatures'}; 
        
    end
    
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

    end
    
    methods % setters
        
    end 
    
    methods % other utilities
        
        function addProp(obj, name, value, make_permanent)
            
            if nargin<3
                value = [];
            end
            
            if nargin<4 || isempty(make_permanent)
                make_permanent = 0;
            end
            
            if ~isprop(obj, name)
                addprop(obj, name); 
                obj.list_dynamic_props{end+1} = name; 
                
                if make_permanent
                    obj.list_permanent_props{end+1} = name; % add to the list of things that don't get cleared each batch
                end
                
            end
            
            obj.(name) = value; 
            
        end
        
%         function takeFrom(obj, other)
%             error('Do not use this!');
%             obj.copyFrom(other);
%             other.clear; % this is very dangerous, especially when taking from BufferWheel, as it resets the buffers as well... 
%             
%         end
        
        function copyTo(obj, other) % moves all the data from this object to another object
            
            copyFrom(other, obj);
            
        end
        
        function copyFrom(obj, other) % moves all the data from another object to this one
        % Usage: copyFrom(obj, other)
        % Takes all the data products (defined as the list of properties of 
        % the AstroData class) and copies them from one object to another. 
        % The two objects must be of AstroData type or any class that 
        % inherits from it. 
        % For example, the Analysis object can simply copyFrom the Reader 
        % object, and it is guaranteed to get all the data products. 
        
            if nargin==1, help('file.AstroData.copyFrom'); return; end
        
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
        
        function clear(obj) % set to empty all the data products, and get ready for a new batch
            
%             list = properties(obj);
            list = properties(file.AstroData);
            
            for ii = 1:length(list)
                
%                 P = findprop(obj, list{ii});
                
%                 if ~isempty(P.DefiningClass) && isequal(P.DefiningClass.Name, 'file.AstroData') && ...
                  if ~contains(list{ii}, obj.list_permanent_props) % also exclude some other properties
                    obj.(list{ii}) = []; % clear everything, not including dynamic properties
                end
                
            end
            
            % clear the dynamically allocated properties as well... 
            for ii = 1:length(obj.list_dynamic_props)
                obj.(obj.list_dynamic_props{ii}) = []; 
            end
            
            
%             obj.images = [];
%             
%             obj.timestamps = [];
%             obj.t_start = [];
%             obj.t_end = [];
%             obj.t_end_stamp = [];
%             obj.juldates = [];
%             
%             obj.cutouts = [];
%             
%             obj.positions = [];
%             obj.coordinates = [];
%             obj.magnitudes = [];
%             obj.temperatures = [];
%         
%             obj.cutouts_bg = [];
%             obj.positions_bg = [];
%             
%             obj.stack = [];
%             obj.num_sum = [];
%             
%             obj.psfs = [];
%             obj.sampling_psf = [];
%             
%             obj.fluxes = [];
%             obj.errors = [];
%             obj.areas = [];
%             obj.backgrounds = [];
%             obj.variances = [];
%             obj.offsets_x = [];
%             obj.offsets_y = [];
%             obj.centroids_x = [];
%             obj.centroids_y = [];
%             obj.widths = [];
%             obj.bad_pixels = [];
%             obj.flags = [];
            
        end
                
    end
    
end