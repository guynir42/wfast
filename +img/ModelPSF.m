classdef ModelPSF < handle
% Calculate some properties of the PSF using the cutouts and the photometry. 
% Use input(cutouts, dx, dy, flux, positions) to provide all necessary data
% for calculating the stacked PSF (which is only good if the PSF is fairly 
% uniform across the field) and the FWHM as a function of cutout position. 
% 
% Use num_stars to determine how many cutouts (from the brightest first)
% to use in the FWHM calculation. Too many cutouts take a long time to 
% process and will contain faint stars where the width is un-measurable. 
%
% The surf_fit and surf_coeffs contain data on the FWHM (in arcsec) as a 
% function of star position. surf_fit.xc and .yc are the x and y centers. 
% Use x-xc and y-yc when inputting the star position into the polynomial 
% fit result to get the correct FWHM. 
% Take the 1st coefficient to get the seeing at the center of the field. 
%
% The calculation uses util.img.align based on the photometric centroids, 
% then sums all frames for each star to increase S/N and reduce runtime, 
% then calculates the width using util.img.fwhm with three different filter
% types: gaussian, generalized gaussian, and defocus annulus. 
% 
% 
    
    properties(Transient=true)
        
    end
    
    properties % objects
        
        head@head.Header; 
        
    end
    
    properties % inputs/outputs
        
        offsets_x;
        offsets_y;
        fluxes;
        positions;
        cutouts;
               
        cutouts_shifted;
        stack; 
        mask;
        
        fwhm; % arcsec
        fwhm_pix; % pixels
        
        m2x;
        m2y;
        mxy;
        maj_axis;
        min_axis;
        angle;
        
        surf_fit; 
        surf_coeffs;
        
    end
    
    properties % switches/controls
        
        use_mex = 1;
        radius = 5;
        
        use_gaussian = 0;
        gauss_sigma = 5;
        
        num_stars = 100; 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = ModelPSF(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'img.ModelPSF')
                if obj.debug_bit>1, fprintf('ModelPSF copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('ModelPSF constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            obj.offsets_x = [];
            obj.offsets_y = [];
            obj.cutouts = [];
            obj.cutouts_shifted = [];
            obj.stack = []; 
            obj.mask = [];

%             obj.fwhm = [];

%             obj.m2x = [];
%             obj.m2y = [];
%             obj.mxy = [];
%             obj.maj_axis = [];
%             obj.min_axis = [];
%             obj.angle = [];
    
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function input(obj, varargin)
            
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('cutouts', [], 'images');
            input.input_var('offsets_x', [], 'dx', 9);
            input.input_var('offsets_y', [], 'dx', 9);
            input.input_var('fluxes', []); 
            input.input_var('positions', []); 
            input.input_var('radius', [], 'aperture'); 
            input.scan_vars(varargin{:});
            
            if isempty(input.cutouts)
                return;
            end
            
            N = min(obj.num_stars, size(input.cutouts,4)); 
            
            obj.cutouts = input.cutouts(:,:,:,1:N);
            obj.offsets_x = input.offsets_x(:,1:N);
            obj.offsets_y = input.offsets_y(:,1:N);
            obj.fluxes = input.fluxes(:,1:N); 
            if ~isempty(input.positions)
                obj.positions = input.positions(1:N,:); 
            end
            
            S = util.vec.imsize(obj.cutouts);
            
            if ~isempty(input.radius)
                obj.radius = input.radius;
            end
            
            dx = obj.offsets_x;
            dy = obj.offsets_y;
            
            % average offsets (I would prefer the flux-weighted average)
            Adx = repmat(nanmean(dx,2), [1,size(dx,2)]);
            Ady = repmat(nanmean(dy,2), [1,size(dy,2)]); 
            
            dx = (abs(dx)<=S(2)/2).*dx + (isnan(dx) | abs(dx)>S(2)/2).*Adx;
            dy = (abs(dy)<=S(1)/2).*dy + (isnan(dy) | abs(dy)>S(1)/2).*Ady;
            
            if isempty(dx) || isempty(dy)
                obj.cutouts_shifted = obj.cutouts;
            else
                obj.cutouts_shifted = util.img.align(obj.cutouts, dx, dy); 
            end
            
            obj.calcStack;
            
            if ~isempty(obj.positions)
                obj.calcSurfaceFit; 
            end
            
        end
        
        function calcStack(obj)
            
            import util.stat.sum2;
            
            obj.stack = nansum(nansum(obj.cutouts_shifted,3),4);
            
            if nnz(obj.stack)==0 || nnz(~isnan(obj.stack))==0
                obj.angle = NaN;
                obj.maj_axis = NaN;
                obj.min_axis = NaN;
                obj.fwhm = NaN;
                obj.fwhm_pix = NaN;
                return; 
            end
            
            if obj.use_gaussian
                obj.stack = obj.stack.*util.shapes.gaussian(obj.gauss_sigma, 'size', size(obj.cutouts)); 
            end
            
            obj.mask = util.shapes.ellipse('radius', obj.radius, 'size', size(obj.cutouts));
            
            c = util.vec.imsize(obj.cutouts); 
            [X, Y] = meshgrid((1:c(2))-floor(c(2)/2)-1, (1:c(1))-floor(c(1)/2)-1); 
            
            I = obj.stack.*obj.mask;
            
            S = sum2(I);
            
            obj.m2x = sum2(X.^2.*I)./S;
            obj.m2y = sum2(Y.^2.*I)./S;
            obj.mxy = sum2(X.*Y.*I)./S;
            
            M = [obj.m2x obj.mxy; obj.mxy obj.m2y];
            
            [R,E] = eig(M);
            
            obj.angle = asind(R(4));
            obj.maj_axis = max(sqrt(diag(E)));
            obj.min_axis = min(sqrt(diag(E)));
            
%             obj.fwhm = sqrt(mean([obj.m2x,obj.m2y])).*2.355;
            obj.fwhm_pix = util.img.fwhm(I, 'method', 'filters');
            
            if ~isempty(obj.head) && ~isempty(obj.head.SCALE)
                obj.fwhm = obj.fwhm_pix.*obj.head.SCALE;
            end
            
        end
        
        function calcSurfaceFit(obj)
            
            C = nansum(obj.cutouts_shifted,3); % stack the individual frames
            
            C = C - util.stat.corner_median(C); % can we figure out a better way to remove the background? 
            
            w = util.img.fwhm(C, 'method', 'filters', 'defocus', 1, 'generalized', 5, 'step', 0.25, 'min_size', 1);
            w = util.vec.tocolumn(w); 
            
            x = obj.positions(:,1); 
            y = obj.positions(:,2); 
            
            if ~isempty(obj.head)
                xc = floor(obj.head.NAXIS2/2)+1; 
                yc = floor(obj.head.NAXIS1/2)+1;                 
            else
                xc = nanmean(x);
                yc = nanmean(y); 
            end
            
            F = util.vec.tocolumn(nansum(obj.fluxes,1)); % average flux is used as weight
            
            obj.surf_fit = util.fit.surf_poly(x-xc, y-yc, w, 'weights', sqrt(abs(F)), 'order', 2, 'sigma', 3, 'iterations', 2, 'plot', 0); 
            obj.surf_fit.xc = xc;
            obj.surf_fit.yc = yc;
            
            obj.surf_coeffs = obj.surf_fit.coeffs;
            
            obj.fwhm_pix = obj.surf_coeffs(1);
            
            if ~isempty(obj.head) && ~isempty(obj.head.SCALE)
                obj.fwhm = obj.fwhm_pix.*obj.head.SCALE;
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

