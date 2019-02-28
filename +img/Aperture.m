classdef Aperture < handle
   
    properties
       
        mask;
        weight;        
        
        plateau_size=50;
        gaussian_size=70;
        annulus_size;
        tile_size=128;
        
        use_moments = 1;
        
        debug_bit = 1;
        
    end
    
    properties(Hidden=true)
         
        mask_f;
        
        version = 1.00;
        
    end
    
    properties(Dependent=true)
        
    end
    
    methods % constructor
        
        function obj = Aperture(other)
            
            if nargin>0 && isa(other, 'img.Aperture')
               
                if obj.debug_bit, fprintf('Aperture copy-constructor v%4.2f\n', obj.version); end
                
                obj = util.oop.full_copy(obj, other);
                
            else
        
                if obj.debug_bit, fprintf('Aperture constructor v%4.2f\n', obj.version); end
                                
            end
            
        end
        
    end
    
    methods % getters
        
        function M = get.mask(obj)
            
            if isempty(obj.tile_size)
                M = [];
                return;
            end
            
            if isempty(obj.mask)
                obj.mask = obj.makeMask;
            end
            
            M = obj.mask;
            
        end
        
        function M = get.mask_f(obj)
           
            if isempty(obj.mask_f)
                obj.mask_f = util.fft.fftshift2(fft2(obj.mask));
            end
            
            M = obj.mask_f;
            
        end
        
        function val = get.weight(obj)
           
            if isempty(obj.weight)
                obj.weight = obj.findWeight;
            end
            
            val = obj.weight;                
            
        end
        
        function val = pixels(obj, val)
            
            if ~isempty(val) && val>0 && val<=1 && ~isempty(obj.tile_size)                
                val = val.*obj.tile_size(1);
            end
            
        end
        
    end
    
    methods % setters
    
        function reset(obj)
           
            obj.mask = [];
            obj.mask_f = [];
            obj.weight = [];
            
        end
        
        function set.plateau_size(obj, aperture)
            
            aperture = obj.pixels(aperture);
            
            if obj.plateau_size~=aperture
                obj.reset;
            end
            
            obj.plateau_size = aperture;
            
        end
        
        function set.gaussian_size(obj, aperture)
            
            aperture = obj.pixels(aperture);
            
            if obj.gaussian_size~=aperture
                obj.reset;
            end
            
            obj.gaussian_size = aperture;
            
        end
        
        function set.annulus_size(obj, aperture)
            
            aperture = obj.pixels(aperture);
            
            if obj.annulus_size~=aperture
                obj.reset;
            end
            
            obj.annulus_size = aperture;
            
        end
        
        function set.tile_size(obj, s)
            
            if any(obj.tile_size~=s)
                obj.tile_size = s;
                obj.reset;
            end
                        
        end
        
    end
    
    methods % mask generation / calculations
    
        function f = makeMask(obj)
            
%             if isempty(obj.gaussian_size) % || obj.gaussian_size<obj.plateau_size
%                 gaus_r = obj.pixels(obj.plateau_size)/2; % no tapering
%             else
%                 gaus_r = obj.pixels(obj.gaussian_size)/2; % tapering (gaussian radius)
%             end
            
            if isempty(obj.plateau_size)
                plat_r = Inf;
            else
                plat_r = obj.pixels(obj.plateau_size)/2; % flat area
            end

            [x,y] = meshgrid(-floor(obj.tile_size/2):ceil(obj.tile_size/2)-1);
            
            dist = sqrt(x.^2 + y.^2);
            
            if isempty(obj.pixels(obj.gaussian_size))
                f = zeros(obj.tile_size);
            else
                gaus_r = obj.pixels(obj.gaussian_size)/2;
                f = exp(-0.5*(dist-plat_r).^2./gaus_r.^2);
            end
            
            f(dist<plat_r) = 1;
            
            if ~isempty(obj.annulus_size)
                f(dist<obj.pixels(obj.annulus_size)/2) = 0; % inner radius that is zero...
            end
            
        end
        
        function val = findWeight(obj)
           
            val = util.stat.sum2(obj.mask);
            
        end
        
    end
    
    methods % photometry
               
        function [LC, timestamps, images_adjusted] = input(obj, varargin) % to be depricated! (this is now done in img.Photometry)
            
            error('to be depricated! (this is now done in img.Photometry)');
            
            images = [];
            timestamps = [];            
            
            for ii = 1:length(varargin)
                
                if isa(varargin{ii}, 'img.DataSet') && util.check_str(varargin{ii}.type, 'images')
                    images = varargin{ii}.data;
                elseif isnumeric(varargin{ii}) && ~isvector(varargin{ii})
                    images = varargin{ii};
                elseif isa(varargin{ii}, 'img.DataSet') && util.check_str(varargin{ii}.type, 'timestamps')
                    timestamps = varargin{ii}.data;
                elseif isnumeric(varargin{ii}) && isvector(varargin{ii})
                    timestamps = varargin{ii};
                end
                
            end
           
            if isempty(images)
                error('must give some image inputs to Aperture!');
            end
            
            lightcurves = util.sum2(bsxfun(@times, obj.mask, images)); % assuming no shifts in positions are required... 
            images_adjusted = images;
            
            if obj.use_moments % shift the images after applying mask then sum
            
                csize = size(images,1);
                [x,y] = meshgrid((1:csize)-floor(csize/2)-1, (1:csize)-floor(csize/2)-1);
                
                for k = 1:obj.use_moments % this switch also tells number of iterations... 
                    
                    masked_stars = bsxfun(@times, obj.mask, images_adjusted); % dim 1,2 are y,x. dim 3 is frames, dim 4 is cutouts
                                        
                    moment_x = round(util.sum2(bsxfun(@times, x, masked_stars))./lightcurves);
                    moment_y = round(util.sum2(bsxfun(@times, y, masked_stars))./lightcurves);
                    
                    % readjust positions...
                    for ii = 1:size(images_adjusted, 3)
                        for jj = 1:size(images_adjusted,4)
                            images_adjusted(:,:,ii,jj) = circshift(images_adjusted(:,:,ii,jj), [-moment_y(1,1,ii,jj), -moment_x(1,1,ii,jj)]);
                        end
                    end
                    masked_stars = bsxfun(@times, obj.mask, images_adjusted); % dim 1,2 are y,x. dim 3 is frames, dim 4 is cutouts
                    
                    lightcurves = util.sum2(masked_stars); % updated to new position... 
                    
                end
                
            end
            
%             lightcurves = util.sum2(bsxfun(@times, obj.mask, images))./obj.weight;            
            
            lightcurves = permute(lightcurves, [3, 4, 1, 2]);
            
            if nargout==1
                LC = img.DataSet(lightcurves);
                LC.N = size(lightcurves,1);
                LC.axes_labels = {'timestamps', 'cutouts'};
                LC.axes_values = {timestamps, []};
                LC.axes_units = {'seconds', ''};
            else
                LC = lightcurves;                
            end
            
        end
        
    end
    
    methods % plotting tools (and radial profiles...)
       
        function [out_s, out_m] = showRadialProfile(obj, I)
 
            r = 2:2:floor(obj.tile_size/2);
            
            m = zeros(length(r),1);
            s = zeros(length(r),1);
            
            for ii = 1:length(r)
               
                p = util.img.annulusPixels(I, r(ii)-2, r(ii));
                
                m(ii) = mean(p);
                s(ii) = std(p);
                
            end
            
            errorbar(r,m,s);
            title('Radial profile');
            xlabel('distance from center of image (pixels)');
            ylabel('mean and std of pixels in annulus');
            
            if nargout>0
                out_s = s;
            end
            
            if nargout>1
                out_m = m;
            end
            
        end
        
        function showMask(obj)
           
            p = obj.mask(end/2+1:end, end/2+1);
            
            plot(p, '+');
            
        end
        
    end
    
end