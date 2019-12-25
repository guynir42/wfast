classdef Photometry < handle

    properties(Transient=true)
        
        gui;
        
    end
    
    properties % objects
        
        pars_struct; % a struct with some housekeeping about how the photometry was done
        
    end
    
    properties % inputs/outputs
        
        % inputs
        positions; % cutout center positions, given from analysis/acquisition
        
        cutouts; % cutouts as given from analysis/acquisition (after calibration)
        cutouts_proc; % subtract backround, kill spurius negative pixels, etc
        
        var_map; % can be scalar or a map of the same size as cutouts
        gain; % can be scalar or a map of the same size as cutouts
        
        timestamps;
        
        psf; % either given or found from the image (to be depricated)
        
        % outputs
        fluxes;
        errors;
        areas;
        backgrounds;
        variances;
        centroids_x;
        centroids_y;
        offsets_x;
        offsets_y;
        widths;
        bad_pixels;
        flags; 
        
        average_flux;
        average_background;
        average_variance;
        average_offset_x;
        average_offset_y;
        average_width;
        
    end
    
    properties % switches/controls
        
        use_mex = 1; % use the (now old) mex function for faster processing 
        use_new_method = 1; % use the new mex function that is much faster and more accurate
        num_threads = 2; % for multithreaded mex photometry
        use_backgrounds = 1; % remove background from individual cutout
        use_self_psf = 0; % use image as a proxy for its own PSF
        
        iterations = 2; % repeat the photometry and relocate the aperutre to the centroid
        use_basic = 1;
        use_centering = 1;
        use_gaussian = 1;
        use_aperture = 1;
        use_forced = 1;
        
        corner_size = 0.15; % fraction of the cut_size or pixel value (must be smaller than cut_size!)
        aperture = 5;
        annulus = 8;
        annulus_outer = 11; % empty means take the rest of the cutout
        gauss_sigma = 2;
        gauss_thresh = 1e-6;
        
        use_fitter = 0;
%         use_gaussian_psf = 1; % just use a simple PSF 
        
        % fitting parameters
        bg_noise = 1;
        saturation_value = 50000;
        
        % not really using these...
        use_pixel_calibration = 0;
        pixel_cal_iterations = 3;
        
        percentile = 0.5; % what fraction of the best flux are we taking for averages
        
        show_num_stars = 10;
        show_num_frames = 3;
        
        show_bit = 0; % not really using this one...
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        cut_size;
        
    end
    
    properties(Hidden=true)
        
        default_iterations;
        
        default_corner_size;
        default_aperture;
        default_annulus;
        default_gauss_sigma;
        default_precentile;
        
        default_seeing = 1;
        default_ap_multiplier = 3.5;
        
        X; % output from meshgrid
        Y; % output from meshgrid
        
        cut_size_latest;
        
        fluxes_basic;
        errors_basic;
        areas_basic;
        backgrounds_basic;
        variances_basic;
        offsets_x_basic;
        offsets_y_basic;
        centroids_x_basic;
        centroids_y_basic;
        widths_basic;
        bad_pixels_basic;
        flags_basic;
        
        fluxes_ap;
        errors_ap;
        areas_ap;
        backgrounds_ap;
        variances_ap;
        offsets_x_ap;
        offsets_y_ap;
        centroids_x_ap;
        centroids_y_ap;
        widths_ap;
        bad_pixels_ap;
        flags_ap;
        
        fluxes_gauss;
        errors_gauss;
        areas_gauss;
        backgrounds_gauss;
        variances_gauss;
        offsets_x_gauss;
        offsets_y_gauss;
        centroids_x_gauss;
        centroids_y_gauss;
        widths_gauss;
        bad_pixels_gauss;
        flags_gauss;

        fluxes_forced;
        errors_forced;
        areas_forced;
        backgrounds_forced;
        variances_forced;
        offsets_x_forced;
        offsets_y_forced;
        centroids_x_forced;
        centroids_y_forced;
        widths_forced;
        bad_pixels_forced;
        flags_forced;
        
        fluxes_fit;
        errors_fit;
        areas_fit;
        backgrounds_fit;
        variances_fit;
        offsets_x_fit;
        offsets_y_fit;
        centroids_x_fit;
        centroids_y_fit;
        widths_fit;
        bad_pixels_fit;
        flags_fit;
        
        version = 1.06;
        
    end
    
    methods % constructor
        
        function obj = Photometry(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'img.Photometry')
                if obj.debug_bit, fprintf('Photometry copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Photometry constructor v%4.2f\n', obj.version); end
                util.oop.save_defaults(obj);
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.var_map = [];
            obj.gain = [];

            obj.average_flux = [];
            obj.average_background = [];
            obj.average_variance = [];
            obj.average_offset_x = [];
            obj.average_offset_y = [];
            obj.average_width = [];
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            obj.cutouts = [];
            obj.cutouts_proc = [];
            obj.timestamps = [];
            obj.psf = [];
            
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
            
            obj.fluxes_basic = [];
            obj.errors_basic = [];
            obj.areas_basic = [];
            obj.backgrounds_basic = [];
            obj.variances_basic = [];
            obj.offsets_x_basic = [];
            obj.offsets_y_basic = [];
            obj.centroids_x_basic = [];
            obj.centroids_y_basic = [];
            obj.widths_basic = [];
            obj.bad_pixels_basic = [];
            obj.flags_basic = [];
            
            obj.fluxes_ap = [];
            obj.errors_ap = [];
            obj.areas_ap = [];
            obj.backgrounds_ap = [];
            obj.variances_ap = [];
            obj.offsets_x_ap = [];
            obj.offsets_y_ap = [];
            obj.centroids_x_ap = [];
            obj.centroids_y_ap = [];
            obj.widths_ap = [];
            obj.bad_pixels_ap = [];
            obj.flags_ap = [];
            
            obj.fluxes_gauss = [];
            obj.errors_gauss = [];
            obj.areas_gauss = [];
            obj.backgrounds_gauss = [];
            obj.variances_gauss = [];
            obj.offsets_x_gauss = [];
            obj.offsets_y_gauss = [];
            obj.centroids_x_gauss = [];
            obj.centroids_y_gauss = [];
            obj.widths_gauss = [];
            obj.bad_pixels_gauss = [];
            obj.flags_gauss = [];
            
            obj.fluxes_fit = [];
            obj.errors_fit = [];
            obj.areas_fit = [];
            obj.backgrounds_fit = [];
            obj.variances_fit = [];
            obj.offsets_x_fit = [];
            obj.offsets_y_fit = [];
            obj.centroids_x_fit = [];
            obj.centroids_y_fit = [];
            obj.widths_fit = [];
            obj.bad_pixels_fit = [];
            obj.flags_fit = [];
            
        end
        
    end
    
    methods % getters
        
        function val = get.cut_size(obj)
            
            val = obj.cut_size_latest;
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function input(obj, varargin) 
            
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('cutouts', [], 'images');
            input.input_var('variance', [], 'var_map');
            input.input_var('gain', [], 'gain_map', 'gain_scalar');
            input.input_var('positions', []);
            input.input_var('timestamps', [], 'times');
            input.scan_vars(varargin{:});
            
            if isa(input.cutouts, 'single')
                obj.cutouts = input.cutouts;
            else
                obj.cutouts = single(input.cutouts);
            end
            
            obj.positions = input.positions;
            obj.var_map = input.variance;
            obj.gain = input.gain;
            obj.timestamps = input.timestamps;
            
            obj.cut_size_latest = size(input.cutouts);
            obj.cut_size_latest = obj.cut_size(1:2);
            
            if obj.use_new_method 
                
                s = util.img.photometry2(single(obj.cutouts), 'iterations', obj.iterations, ...
                    'radii', obj.aperture, 'annulus', [obj.annulus, obj.annulus_outer], 'sigma', obj.gauss_sigma, ...
                    'use_gaussian', obj.use_gaussian, 'use_centering', obj.use_centering, ...
                    'use_apertures', obj.use_aperture, 'use_forced', obj.use_forced, ...
                    'threads', obj.num_threads, 'debug_bit', obj.debug_bit); 
                
                obj.pars_struct = s.parameters;
                
                obj.fluxes_basic = s.raw_photometry.flux;
                obj.areas_basic = s.raw_photometry.area;
                obj.errors_basic = s.raw_photometry.error;
                obj.backgrounds_basic = s.raw_photometry.background;
                obj.variances_basic = s.raw_photometry.variance;
                obj.offsets_x_basic = s.raw_photometry.offset_x;
                obj.offsets_y_basic = s.raw_photometry.offset_y;
                obj.widths_basic = s.raw_photometry.width;
                obj.bad_pixels_basic = s.raw_photometry.bad_pixels;
                obj.flags_basic = s.raw_photometry.flag;
                
                if ~isempty(obj.positions)
                    obj.centroids_x_basic = obj.offsets_x_basic + obj.positions(:,1)';
                    obj.centroids_y_basic = obj.offsets_y_basic + obj.positions(:,2)';
                end
                
                % update the newest values
                obj.fluxes = obj.fluxes_basic;
                obj.areas = obj.areas_basic;
                obj.errors = obj.errors_basic;
                obj.backgrounds = obj.backgrounds_basic;
                obj.variances = obj.variances_basic;
                obj.offsets_x = obj.offsets_x_basic;
                obj.offsets_y = obj.offsets_y_basic;
                obj.widths = obj.widths_basic;
                obj.bad_pixels = obj.bad_pixels_basic;
                obj.flags = obj.flags_basic;
                                
                if obj.use_gaussian

                    obj.fluxes_gauss = s.gaussian_photometry.flux;
                    obj.areas_gauss = s.gaussian_photometry.area;
                    obj.errors_gauss = s.gaussian_photometry.error;
                    obj.backgrounds_gauss = s.gaussian_photometry.background;
                    obj.variances_gauss = s.gaussian_photometry.variance;
                    obj.offsets_x_gauss = s.gaussian_photometry.offset_x;
                    obj.offsets_y_gauss = s.gaussian_photometry.offset_y;
                    obj.widths_gauss = s.gaussian_photometry.width;
                    obj.bad_pixels_gauss = s.gaussian_photometry.bad_pixels;
                    obj.flags_gauss = s.gaussian_photometry.flag;

                    if ~isempty(obj.positions)
                        obj.centroids_x_gauss = obj.offsets_x_gauss + obj.positions(:,1)';
                        obj.centroids_y_gauss = obj.offsets_y_gauss + obj.positions(:,2)';
                    end
                    
                    % update the newest values
                    obj.fluxes = obj.fluxes_gauss;
                    obj.areas = obj.areas_gauss;
                    obj.errors = obj.errors_gauss;
                    obj.backgrounds = obj.backgrounds_gauss;
                    obj.variances = obj.variances_gauss;
                    obj.offsets_x = obj.offsets_x_gauss;
                    obj.offsets_y = obj.offsets_y_gauss;
                    obj.widths = obj.widths_gauss;
                    obj.bad_pixels = obj.bad_pixels_gauss;
                    obj.flags = obj.flags_gauss;
                    
                end
                
                if obj.use_aperture

                    obj.fluxes_ap = s.apertures_photometry.flux;
                    obj.areas_ap = s.apertures_photometry.area;
                    obj.errors_ap = s.apertures_photometry.error;
                    obj.backgrounds_ap = s.apertures_photometry.background;
                    obj.variances_ap = s.apertures_photometry.variance;
                    obj.offsets_x_ap = s.apertures_photometry.offset_x;
                    obj.offsets_y_ap = s.apertures_photometry.offset_y;
                    obj.widths_ap = s.apertures_photometry.width;
                    obj.bad_pixels_ap = s.apertures_photometry.bad_pixels;
                    obj.flags_ap = s.apertures_photometry.flag;

                    if ~isempty(obj.positions)
                        obj.centroids_x_ap = obj.offsets_x_ap + obj.positions(:,1)';
                        obj.centroids_y_ap = obj.offsets_y_ap + obj.positions(:,2)';
                    end 
                    
                    % update the newest values
                    obj.fluxes = obj.fluxes_ap;
                    obj.areas = obj.areas_ap;
                    obj.errors = obj.errors_ap;
                    obj.backgrounds = obj.backgrounds_ap;
                    obj.variances = obj.variances_ap;
                    obj.offsets_x = obj.offsets_x_ap;
                    obj.offsets_y = obj.offsets_y_ap;
                    obj.widths = obj.widths_ap;
                    obj.bad_pixels = obj.bad_pixels_ap;
                    obj.flags = obj.flags_ap;
                    
                end
                
                if obj.use_forced
                    
                    obj.fluxes_forced = s.forced_photometry.flux;
                    obj.areas_forced = s.forced_photometry.area;
                    obj.errors_forced = s.forced_photometry.error;
                    obj.backgrounds_forced = s.forced_photometry.background;
                    obj.variances_forced = s.forced_photometry.variance;
                    obj.offsets_x_forced = s.forced_photometry.offset_x;
                    obj.offsets_y_forced = s.forced_photometry.offset_y;
                    obj.widths_forced = s.forced_photometry.width;
                    obj.bad_pixels_forced = s.forced_photometry.bad_pixels;
                    obj.flags_forced = s.forced_photometry.flag;

                    if ~isempty(obj.positions)
                        obj.centroids_x_forced = obj.offsets_x_forced + obj.positions(:,1)';
                        obj.centroids_y_forced = obj.offsets_y_forced + obj.positions(:,2)';
                    end 
                    
                end
                
                if ~isempty(obj.positions)
                    obj.centroids_x = obj.offsets_x + obj.positions(:,1)';
                    obj.centroids_y = obj.offsets_y + obj.positions(:,2)';
                end
                
            else
                
                if obj.use_basic
                    obj.calcBasic;
                end

                if obj.use_aperture
                    obj.calcAperture;
                end

                if obj.use_gaussian
                    obj.calcGaussian;
                end

                if obj.use_fitter
                    obj.calcFit;
                end
                
                obj.updatePars;
                obj.flags = zeros(size(obj.fluxes)); 
                
            end
            
            obj.updateAverages;
            
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.update;
            end
            
        end
        
        function preprocess(obj) % we don't even use this... 
            
            I = obj.cutouts;
            
            I(I<-3.*util.stat.std2(I)) = NaN;
            
            obj.cutouts_proc = I;
            
        end
        
        function [fluxes, errors, areas, backgrounds, variances, offsets_x, offsets_y, widths, bad_pixels] = calculate(obj, shape, bg_shape, iterations)

            import util.stat.sum2;
            import util.stat.median2;
            import util.stat.var2;
            import util.text.cs;
            
            if nargin<2 || isempty(shape)
                shape = 'circle';
            end
            
            if nargin<3 || isempty(bg_shape)
                bg_shape = 'corner';
            end
            
            if cs(shape, 'none')
                ap_shape = @(x,y) ones(obj.cut_size, 'like', obj.cutouts);
            elseif cs(shape, 'circle', 'aperture')
                ap_shape = @(x,y) obj.makeCircle(x,y);
            elseif cs(shape, 'gaussian')
                ap_shape = @(x,y) obj.makeGaussian(x,y);
            else
                error('Unknown shape "%s". Use "none", "circle" or "gaussian"', shape);
            end
            
            if cs(bg_shape, 'corner')
                make_bg_shape = @(x,y) obj.makeCorners(x,y);
            elseif cs(bg_shape, 'annulus')
                make_bg_shape = @(x,y) obj.makeAnnulus(x,y);
            else
                error('Unknown background shape "%s". Use "corner" or "annulus"', bg_shape);
            end
            
            if ~isempty(obj.fluxes)
                fluxes = obj.fluxes;
            else
                fluxes = NaN(size(obj.cutouts,3), size(obj.cutouts,4));    
            end
            
            if ~isempty(obj.errors)
                errors = obj.errors;
            else
                errors = NaN(size(obj.cutouts,3), size(obj.cutouts,4));    
            end
            
            if ~isempty(obj.areas)
                areas = obj.areas;
            else
                areas = NaN(size(obj.cutouts,3), size(obj.cutouts,4));    
            end
            
            if ~isempty(obj.backgrounds)
                backgrounds = obj.backgrounds;
            else
                backgrounds = NaN(size(obj.cutouts,3), size(obj.cutouts,4));    
            end
            
            if ~isempty(obj.variances)
                variances = obj.variances;
            else
                variances = NaN(size(obj.cutouts,3), size(obj.cutouts,4));    
            end
            
            if ~isempty(obj.offsets_x)
                offsets_x = obj.offsets_x;
            else
                offsets_x = NaN(size(obj.cutouts,3), size(obj.cutouts,4));    
            end
            
            if ~isempty(obj.offsets_y)
                offsets_y = obj.offsets_y;
            else
                offsets_y = NaN(size(obj.cutouts,3), size(obj.cutouts,4));    
            end
            
            if ~isempty(obj.widths)
                widths = obj.widths;
            else
                widths = NaN(size(obj.cutouts,3), size(obj.cutouts,4));    
            end
            
            if ~isempty(obj.bad_pixels)
                bad_pixels = obj.bad_pixels;
            else
                bad_pixels = NaN(size(obj.cutouts,3), size(obj.cutouts,4));    
            end
            
            % go over each cutout and shift it to the center
            for ii = 1:size(obj.cutouts,4) % go over different stars

                for jj = 1:size(obj.cutouts,3) % go over frames of the same star
                    
                    dx = offsets_x(jj,ii);
                    if isnan(dx), dx = 0; end
                    dy = offsets_y(jj,ii);
                    if isnan(dy), dy = 0; end
                    
                    for kk = 1:iterations
                        
                        I = obj.cutouts(:,:,jj,ii);
                        
                        
                        ap = ap_shape(dx,dy); % should we also let radius/psf_sigma change with the width?? 
                        
                        bad_mask = ap.*isnan(I); % how many (fractional) bad pixels are in the aperture
                        
                        ap(isnan(I)) = nan; % aperture must have NaN values in the same places (for area calculation)
                        
                        ann = make_bg_shape(dx,dy);
                        ann(isnan(I)) = nan; 
                        
                        B = median2(I.*ann); % background per pixel...
                        V = var2(I.*ann); 
                        
                        if isempty(obj.var_map)
                            v = ones(size(obj.cutouts,1), size(obj.cutouts,2), 'single').*V;
                        elseif isscalar(obj.var_map)
                            v = ones(size(obj.cutouts,1), size(obj.cutouts,2), 'single').*obj.var_map;
                        else
                            v = obj.var_map(:,:,jj,ii);
                        end
                        
                        if isempty(obj.gain)
                            g = 1;
                        elseif isscalar(obj.gain)
                            g = obj.gain;
                        else
                            g = obj.gain(:,:,jj,ii);
                        end
                        
                        if obj.use_backgrounds
                            I = I - B;
                        end
                        
                        if obj.use_self_psf
                            e = (g.*I+v);
                            weight = ap.*g.*I./e;
                        else
                            e = v;
                            weight = ap./e;
                        end
                        
                        % new method
                        E = sum2(e);
                        I = I.*weight./sum2(weight); % reweigh the whole image
                        
                        % old method:
                        S = sum2(ap);
%                         ap = ap./S;
%                         I = I.*ap./sum2(ap.^2);

                        m0 = sum2(I);
                        m1x = sum2(I.*obj.X)./m0;
                        m1y = sum2(I.*obj.Y)./m0;
                        m2x = sum2(I.*(obj.X-m1x).^2)./m0;
                        m2y = sum2(I.*(obj.Y-m1y).^2)./m0;
                        mxy = sum2(I.*(obj.X-m1x).*(obj.Y-m1y))./m0;
                        
                        % quality checks:
                        if m0==0
                            if obj.debug_bit>3, disp('m0 is zero!'); end
                            m1x = 0;
                            m1y = 0;
                            m2x = NaN;
                            m2y = NaN;
                            mxy = NaN;
                        end
                        
                        if m2x<0, m2x = NaN; if obj.debug_bit>3, disp('m2x is negative!'); end, end
                        if m2y<0, m2y = NaN; if obj.debug_bit>3, disp('m2y is negative!'); end, end
                        
                        dx = m1x;
                        dy = m1y;
                        
%                         W = sqrt(mean([m2x,m2y], 'omitnan')); % should we add mxy??
                        
                        M = [m2x mxy; mxy m2y];
                        
                        if any(isnan(M(:)))
                            W = NaN;
                        else
                            [X, D] = eig(M); % use X to calculate rotation angle? 
                            W = mean(sqrt(diag(D)), 'omitnan');
                        end
                        
                    end

                    fluxes(jj,ii) = m0;
                    areas(jj,ii) = S;
                    errors(jj,ii) = E; 
                    backgrounds(jj,ii) = B;
                    variances(jj,ii) = V;
                    offsets_x(jj,ii) = dx;
                    offsets_y(jj,ii) = dy;
                    widths(jj,ii) = W;                    
                    bad_pixels(jj,ii) = sum2(bad_mask);
                        
                    % add break point if dx and dy don't change much...
                    
                end

            end
            
        end
        
        function makeAxes(obj)
            
            c = size(obj.cutouts); c = c(1:2);
            [obj.X,obj.Y] = meshgrid((1:c(2))-floor(c(2)/2)-1, (1:c(1))-floor(c(1)/2)-1); 
            
        end
        
        function val = makeCorners(obj, ~, ~)
            
            if isempty(obj.X) || isempty(obj.Y)
                obj.makeAxes;
            end
            
            val = zeros(obj.cut_size);
            
            if obj.corner_size<1
                pixels = obj.corner_size.*obj.cut_size;
            elseif obj.corner_size<obj.cut_size
                pixels = obj.corner_size;
            else
                error('Wrong corner_size= %d. Input a fraction or a pixel number smaller than cut_size= %d', obj.corner_size, obj.cut_size);
            end
            
            pixels = round(pixels);
            
            val(1:pixels, 1:pixels) = 1;
            val(end-pixels+1:end, 1:pixels) = 1;            
            val(1:pixels, end-pixels+1:end) = 1;
            val(end-pixels+1:end, end-pixels+1:end) = 1;
            
            if isa(obj.cutouts, 'single')
                val = single(val);
            else
                val = double(val);
            end
            
            val(val==0) = NaN;
            
        end
        
        function val = makeCircle(obj, dx, dy)
            
            if nargin<2 || isempty(dx)
                dx = 0;
            end
            
            if nargin<3 || isempty(dy)
                dy = 0;
            end
            
            if isempty(obj.X) || isempty(obj.Y)
                obj.makeAxes;
            end
            
            R = obj.aperture;
            
            if isnan(R) % use NaN to indicate automatically choose radius
                if ~isempty(obj.widths)
                    R = obj.default_ap_multiplier*util.stat.median2(obj.widths);
                elseif ~isempty(obj.default_aperture) && ~isnan(obj.default_aperture)
                    R = obj.default_aperture;
                else
                    R = obj.default_seeing.*obj.default_ap_multiplier;
                end
            end
            
            r = sqrt((obj.X-dx).^2+(obj.Y-dy).^2);
            val = R+0.5-r;
            val(val>1) = 1;
%             val(val<0) = 0;
            val(val<1) = 0;
            
            if isa(obj.cutouts, 'single')
                val = single(val);
            end
            
            val(val==0) = NaN;
            
        end
        
        function val = makeAnnulus(obj, dx, dy)
           
            if nargin<2 || isempty(dx)
                dx = 0;
            end
            
            if nargin<3 || isempty(dy)
                dy = 0;
            end
            
            if isempty(obj.X) || isempty(obj.Y)
                obj.makeAxes;
            end
            
            R1 = obj.annulus;
            
            if isempty(obj.annulus_outer)
                R2 = Inf;
            else
                R2 = obj.annulus_outer;
            end
            
            r = sqrt((obj.X-dx).^2+(obj.Y-dy).^2);
            
            val = r>R1 & r<R2;
            
            if isa(obj.cutouts, 'single')
                val = single(val);
            else
                val = double(val);
            end
            
            val(val==0) = NaN;
            
        end
        
        function val = makeGaussian(obj, dx, dy)
            
            if nargin<2 || isempty(dx)
                dx = 0;
            end
            
            if nargin<3 || isempty(dy)
                dy = 0;
            end
            
            if isempty(obj.X) || isempty(obj.Y)
                obj.makeAxes;
            end
            
            sig = obj.gauss_sigma;
            
            if isnan(sig) %use NaN to indicate automatically choose sigma
                if ~isempty(obj.widths)
                    sig = mean(obj.widths, 1, 'omitnan');
                elseif ~isempty(obj.default_gauss_sigma) && ~isnan(obj.default_gauss_sigma)
                    sig = obj.default_gauss_sigma;
                else
                    sig = obj.default_seeing*2.355;
                end
            end
            
            r = sqrt((obj.X-dx).^2+(obj.Y-dy).^2);
            
            val = exp(-0.5.*(r./sig).^2);
            
            if isa(obj.cutouts, 'single')
                val = single(val);
            end
            
            val(val<obj.gauss_thresh) = NaN;
            
        end
        
        function calcBasic(obj)
            
            import util.stat.sum2;
            
            if obj.use_mex
                [f,e,a,b,v,x,y,wd,p] = util.img.photometry(single(obj.cutouts), 'square', [], 'corners', obj.corner_size,...
                    'widths', obj.widths, 'subtract', obj.use_backgrounds, ...
                    'var_map', obj.var_map, 'use_self', obj.use_self_psf, 'threads', obj.num_threads, 'debug_bit', obj.debug_bit);
            else
                [f,e,a,b,v,x,y,wd,p] = obj.calculate('none', 'corner', 1);
            end
            
            x(abs(x)>size(obj.cutouts,2)) = NaN;
            y(abs(y)>size(obj.cutouts,1)) = NaN;
            wd(wd<0 | wd>size(obj.cutouts,1) | wd>size(obj.cutouts,2)) = NaN;
            
            obj.fluxes_basic = f;
            obj.errors_basic = e;
            obj.areas_basic = a;
            obj.backgrounds_basic = b;
            obj.variances_basic = v;
            
            obj.offsets_x_basic = x;
            obj.offsets_y_basic = y;
            
            if ~isempty(obj.positions)
                obj.centroids_x_basic = x + obj.positions(:,1)';
                obj.centroids_y_basic = y + obj.positions(:,2)';
            end
            
            obj.widths_basic = wd;
            obj.bad_pixels_basic = p;
            
            % update the newest values
            obj.fluxes = obj.fluxes_basic;
            obj.errors = obj.errors_basic;
            obj.areas = obj.areas_basic;
            obj.backgrounds = obj.backgrounds_basic;
            obj.variances = obj.variances_basic;
            
            obj.offsets_x = obj.offsets_x_basic;
            obj.offsets_y = obj.offsets_y_basic;
            
            if ~isempty(obj.positions)
                obj.centroids_x = obj.centroids_x_basic;
                obj.centroids_y = obj.centroids_y_basic;
            end
            
            obj.widths = obj.widths_basic;
            obj.bad_pixels = obj.bad_pixels_basic;
            
        end
        
        function calcAperture(obj)
            
            if obj.use_mex
                [f,e,a,b,v,x,y,wd,p] = util.img.photometry(single(obj.cutouts), 'circle', obj.aperture, 'annulus', obj.annulus,...
                    'widths', obj.widths, 'iterations', obj.iterations, 'subtract', obj.use_backgrounds, ...
                    'var_map', obj.var_map, 'use_self', obj.use_self_psf, 'threads', obj.num_threads, 'debug_bit', obj.debug_bit);
            else
                [f,e,a,b,v,x,y,wd,p] = obj.calculate('circle', 'annulus', obj.iterations);
            end
            
            x(abs(x)>size(obj.cutouts,2)/2+obj.aperture) = NaN;
            y(abs(y)>size(obj.cutouts,1)/2+obj.aperture) = NaN;
            wd(wd<0 | wd>size(obj.cutouts,1) | wd>size(obj.cutouts,2)) = NaN;
            
            obj.fluxes_ap = f;
            obj.errors_ap = e;
            obj.areas_ap = a;
            obj.backgrounds_ap = b;
            obj.variances_ap = v;
            
            obj.offsets_x_ap = x;
            obj.offsets_y_ap = y;
            
            if ~isempty(obj.positions)
                obj.centroids_x_ap = x + obj.positions(:,1)';
                obj.centroids_y_ap = y + obj.positions(:,2)';
            end
            
            obj.widths_ap = wd;
            obj.bad_pixels_ap = p;
            
            % update the newest values
            obj.fluxes = obj.fluxes_ap;
            obj.errors = obj.errors_ap;
            obj.areas = obj.areas_ap;
            obj.backgrounds = obj.backgrounds_ap;
            obj.variances = obj.variances_ap;
            
            obj.offsets_x = obj.offsets_x_ap;
            obj.offsets_y = obj.offsets_y_ap;
            
            if ~isempty(obj.positions)
                obj.centroids_x = obj.centroids_x_ap;
                obj.centroids_y = obj.centroids_y_ap;
            end
            
            obj.widths = obj.widths_ap;
            obj.bad_pixels = obj.bad_pixels_ap;
            
        end
        
        function calcGaussian(obj)
            
            if obj.use_mex
                [f,e,a,b,v,x,y,wd,p] = util.img.photometry(single(obj.cutouts), 'gauss', obj.gauss_sigma, 'annulus', obj.annulus,...
                    'iterations', obj.iterations, 'widths', obj.widths, 'subtract', obj.use_backgrounds, ...
                    'var_map', obj.var_map, 'use_self', obj.use_self_psf, 'threads', obj.num_threads, 'debug_bit', obj.debug_bit);
            else
                [f,e,a,b,v,x,y,wd,p] = obj.calculate('gaussian', 'annulus', obj.iterations);
            end
            
            x(abs(x)>size(obj.cutouts,2)/2+obj.gauss_sigma*3) = NaN;
            y(abs(y)>size(obj.cutouts,1)/2+obj.gauss_sigma*3) = NaN;
            wd(wd<0 | wd>size(obj.cutouts,1) | wd>size(obj.cutouts,2)) = NaN;
            
            obj.fluxes_gauss = f;
            obj.errors_gauss = e;
            obj.areas_gauss = a;
            obj.backgrounds_gauss = b;
            obj.variances_gauss = v;
            
            obj.offsets_x_gauss = x;
            obj.offsets_y_gauss = y;
            
            if ~isempty(obj.positions)
                obj.centroids_x_gauss = x + obj.positions(:,1)';
                obj.centroids_y_gauss = y + obj.positions(:,2)';
            end
            
            obj.widths_gauss = wd*sqrt(2); % correction given because weighing by PSF makes the second moment smaller
            obj.bad_pixels_gauss = p;
            
            % update the newest values
            obj.fluxes = obj.fluxes_gauss;
            obj.errors = obj.errors_gauss;
            obj.areas = obj.areas_gauss;
            obj.backgrounds = obj.backgrounds_gauss;
            obj.variances = obj.variances_gauss;
            
            obj.offsets_x = obj.offsets_x_gauss;
            obj.offsets_y = obj.offsets_y_gauss;
            
            if ~isempty(obj.positions)
                obj.centroids_x = obj.centroids_x_gauss;
                obj.centroids_y = obj.centroids_y_gauss;
            end
            
            obj.widths = obj.widths_gauss;
            obj.bad_pixels = obj.bad_pixels_gauss;
            
        end
        
        function updateAverages(obj)
            
            F = obj.fluxes;
            
            idx = F>util.stat.max2(F).*obj.percentile & ~isnan(F) & obj.flags; % choose only good flux values
            
            % 1D vectors containing the good values only...
            F = obj.fluxes(idx);
            B = obj.backgrounds(idx);
            V = obj.variances(idx);
            W = obj.widths(idx);
            DX = obj.offsets_x(idx);
            DY = obj.offsets_y(idx);
            
            M = mean(F, 'omitnan'); % why not median too?
            
            obj.average_flux = median(F, 'omitnan'); 
            obj.average_background = median(B, 'omitnan');
            obj.average_variance = median(V, 'omitnan');
            obj.average_offset_x = median(F./M.*DX, 'omitnan');
            obj.average_offset_y = median(F./M.*DY, 'omitnan');
            obj.average_width = median(W.*F./M, 'omitnan');
            
        end
        
        function updatePars(obj)
            
            obj.pars_struct = struct;
            obj.pars_struct.used_bg_sub = obj.use_backgrounds;
            obj.pars_struct.use_self_psf = obj.use_self_psf;
            
            if ~isempty(obj.var_map)
                if isscalar(obj.var_map)
                    obj.pars_struct.variance = obj.var_map;
                else
                    obj.pars_struct.variance = 'map';
                end
            end
            
            if ~isempty(obj.gain)
                if isscalar(obj.gain)
                    obj.pars_struct.gain = obj.gain;
                else
                    obj.pars_struct.gain = 'map';
                end
            end
            
            if obj.use_fitter
                obj.pars_struct.signal_method = 'fitter';
                % any other parameters?
            elseif obj.use_gaussian
                obj.pars_struct.signal_method = 'gaussian';
                obj.pars_struct.radius = obj.gauss_sigma;
                obj.pars_struct.radius = obj.gauss_thresh;
                obj.pars_struct.background_method = 'annulus';
                obj.pars_struct.annulus = obj.annulus;
                obj.pars_struct.annulus_outer = obj.annulus_outer;
                obj.pars_struct.iterations = obj.iterations;
            elseif obj.use_aperture
                obj.pars_struct.signal_method = 'aperture';
                obj.pars_struct.radius = obj.aperture;
                obj.pars_struct.background_method = 'annulus';
                obj.pars_struct.annulus = obj.annulus;
                obj.pars_struct.annulus_outer = obj.annulus_outer;
                obj.pars_struct.iterations = obj.iterations;
            elseif obj.use_basic
                obj.pars_struct.signal_method = 'square';
                obj.pars_struct.background_method = 'corners';
                obj.pars_struct.corner_size = obj.corner_size;
            end
            
        end
        
        % everything below this are functions we don't use 
        function calcFit(obj) % need to finish this or remove it
            
            if obj.use_gaussian_psf
                
                for ii = 1:size(obj.cutouts,4) % go over different stars
                
                    disp(['ii= ' num2str(ii)]);
                    
                    for jj = 1:size(obj.cutouts,3) % go over frames of the same star

                        I = obj.cutouts(:,:, jj, ii); % image
                        N = obj.fluxes(jj, ii); % normalization
                        V = abs(I).*obj.gain + obj.bg_noise; % variance estimate
                        B = obj.backgrounds(jj,ii);
                        
                        func = @(p) obj.compareModelGaussian(I, N, V, B, jj,ii, p);

                        p_initial = [obj.offsets_x(jj,ii), obj.offsets_y(jj,ii), 1, obj.widths(jj,ii), 1];

                        p_found = fminsearch(func, p_initial);

                        obj.offsets_x_fit(jj,ii) = p_found(1);
                        obj.offsets_y_fit(jj,ii) = p_found(2);
                        obj.fluxes_fit(jj,ii) = p_found(3).*obj.fluxes(jj,ii);
                        obj.widths_fit(jj,ii) = p_found(4);
                        obj.backgrounds_fit(jj,ii) = p_found(5).*B;
                        
                        % test if fit was good??
                        obj.offsets_x(jj,ii) = obj.offsets_x_fit(jj,ii);
                        obj.offsets_y(jj,ii) = obj.offsets_y_fit(jj,ii);
                        obj.fluxes(jj,ii) = obj.fluxes_fit(jj,ii);
                        obj.widths(jj,ii) = obj.widths_fit(jj,ii);
                        obj.backgrounds(jj,ii) = obj.backgrounds_fit(jj,ii);
                        
                    end

                end
                
            else
                obj.calcPSF;
            
                for ii = 1:size(obj.cutouts,4) % go over different stars
                
                    for jj = 1:size(obj.cutouts,3) % go over frames of the same star

                        func = @(p) obj.compareModelPSF(jj,ii, p);

                        p_initial = [obj.offsets_x(jj,ii), obj.offsets_y(jj,ii), obj.fluxes(jj,ii), 0];

                        p_found = fminsearch(func, p_initial);

                        obj.offsets_x(jj,ii) = p_found(1);
                        obj.offsets_y(jj,ii) = p_found(2);
                        obj.fluxes(jj,ii) = p_found(3);

                    end

                end
            
            end
            
            
        end
        
        function chi2 = compareModelGaussian(obj, I, N, V, B, frame_index, star_index, par_vec)
            
            M = obj.modelGaussian(par_vec); % model

            I(I>obj.saturation_value) = NaN;
            
            D = M - (I-B.*par_vec(5))./N;
            
            chi2 = util.stat.sum2(D.^2./V.*N^2); 
            
            if obj.show_bit
                util.plot.show(D, 'autodyn');
                title(num2str(par_vec));
                xlabel(['star= ' num2str(star_index) ' | frame= ' num2str(frame_index) ' | chi2= ' num2str(chi2)]);
                drawnow;
            end
            
        end
        
        function M = modelGaussian(obj, par_vec)
            
            x_shift = par_vec(1);
            y_shift = par_vec(2);
            flux = par_vec(3);
            sigma = par_vec(4);
            
            if sigma<0
                M = realmax;
            else
                M = util.img.gaussian2('sigma', sigma, 'x_shift', x_shift, 'y_shift', y_shift, 'size', obj.cut_size, 'norm', 1)*flux;
            end
            
        end
        
        function chi2 = compareModelPSF(obj, frame_index, star_index, par_vec)
            
            if size(obj.psf,3)>1
                psf = obj.psf(:,:,frame_index);
            else
                psf = obj.psf;
            end
            
            M = obj.modelPSF(psf, par_vec);
            
            D = M - obj.cutouts(:,:, frame_index, star_index);
            
            chi2 = util.stat.sum2(D.^2);
            
            if obj.show_bit
                util.plot.show(D, 'autodyn');
                title(num2str(par_vec));
                xlabel(['star= ' num2str(star_index) ' | frame= ' num2str(frame_index) ' | chi2= ' num2str(chi2)]);
                drawnow;
            end
            
        end
        
        function M = modelPSF(obj, psf, par_vec)
            
            x_shift = par_vec(1);
            y_shift = par_vec(2);
            flux = par_vec(3);
            bg = par_vec(4);
            
            P = util.img.FourierShift2D(psf, -[x_shift, y_shift]);
            
            M = bg + flux.*P;
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function plot(obj, parent, varargin)
            
            if isempty(obj.fluxes)
                return;
            end
            
            if nargin<2 || isempty(parent) || ~isgraphics(parent)
                if ~isempty(obj.gui) && obj.gui.check
                    parent = obj.gui.panel_image;
                else
                    error('Cannot plot photometry without a parent figure/panel');
                end
            end
            
%             input = util.text.InputVars;
%             input.input_var('parent', []);
%             input.scan_vars(varargin{:});
%             
%             if isempty(input.parent)
%                 input.parent = gcf;
%             end
            
            delete(parent.Children);
            
            Ns = obj.show_num_stars;
            if Ns>size(obj.fluxes,2), Ns = size(obj.fluxes,2); end
                
            Nf = obj.show_num_frames;
            if Nf>size(obj.fluxes,1), Nf = size(obj.fluxes,1); end
            
            ax1 = axes('Parent', parent, 'Position', [0.15 0.52 0.83 0.4]);
            delete(ax1.Children);
            ax1.NextPlot = 'add';
            
            for ii = 1:Nf
                h(ii) = errorbar(ax1, 1:Ns, obj.fluxes(ii, 1:Ns), obj.errors(ii, 1:Ns), '*');
                h(ii).UserData = ['flux= %4.2f, frame= ' num2str(ii)];
                h(ii).ButtonDownFcn = @obj.callback_touch_point;
            end
            
            ax1.XTick = [];
            ax1.XLim = [0.5,Ns+0.5];
            
%             util.plot.inner_title(ax1, 'Flux', 'position', 'right');
            
            ylabel(ax1, 'flux (counts)');

            ax2 = axes('Parent', parent, 'Position', [0.15 0.27 0.83 0.22]);
%             plot(ax2, 1:Ns, obj.backgrounds(1:Nf, 1:Ns), 'o', 1:Ns, obj.variances(1:Nf, 1:Ns), 'x');
            
            h = plot(ax2, 1:Ns, obj.backgrounds(1:Nf, 1:Ns), 'o');
            for ii = 1:length(h)
                h(ii).UserData = ['background= %4.2f, frame= ' num2str(ii)];
                h(ii).ButtonDownFcn = @obj.callback_touch_point;
            end
            
            ax2.XTick = [];
            ax2.XLim = [0.5,Ns+0.5];
            
%             util.plot.inner_title(ax2, 'b/g', 'position', 'right');
            
            ylabel(ax2, 'b/g (counts)');

            ax3 = axes('Parent', parent, 'Position', [0.15 0.02 0.83 0.22]);
            ax3.NextPlot = 'add';
            
            h = plot(ax3, 1:Ns, obj.offsets_x(1:Nf, 1:Ns), 'x'); 
            
            for ii = 1:length(h)
                h(ii).UserData = ['offset_x= %4.2f, frame= ' num2str(ii)];
                h(ii).ButtonDownFcn = @obj.callback_touch_point;
            end
            
            ax3.ColorOrderIndex = 1;
            h = plot(ax3, 1:Ns, obj.offsets_y(1:Nf, 1:Ns), '+');
            
            for ii = 1:length(h)
                h(ii).UserData = ['offset_y= %4.2f, frame= ' num2str(ii)];
                h(ii).ButtonDownFcn = @obj.callback_touch_point;
            end
            
            ax3.ColorOrderIndex = 1;
            h = plot(ax3, 1:Ns, obj.widths(1:Nf, 1:Ns), 'o');
            
            for ii = 1:length(h)
                h(ii).UserData = ['width= %4.2f, frame= ' num2str(ii)];
                h(ii).ButtonDownFcn = @obj.callback_touch_point;
            end
            
            ax3.XTick = [];
            ax3.XLim = [0.5,Ns+0.5];
            
%             util.plot.inner_title(ax3, 'Offsets (x,+)', 'position', 'right');
%             util.plot.inner_title(ax3, 'Widths (o)', 'position', 'left');
            
            ylabel(ax3, 'offset/width');

        end
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = img.gui.PhotGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
        function callback_touch_point(obj, hndl, event)
            
%             str = [hndl.UserData, ', star= ' num2str(event.IntersectionPoint(1)) ', value= ' num2str(event.IntersectionPoint(2))];
            
            str = sprintf([hndl.UserData, ', star= %d'], event.IntersectionPoint(2),event.IntersectionPoint(1));
            
            if obj.debug_bit
                disp('callback: touch_point');
                disp(str);
            end
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.panel_show.button_picker.String = str;
            end
            
        end
        
    end
    
end

 