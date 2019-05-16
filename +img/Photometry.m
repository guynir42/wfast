classdef Photometry < handle

    properties(Transient=true)
        
        gui;
        
    end
    
    properties % objects
        
        
        
    end
    
    properties % inputs/outputs
        
        % inputs
        positions; % cutout center positions, given from analysis/acquisition
        
        cutouts; % cutouts as given from analysis/acquisition (after calibration)
        cutouts_proc; % subtract backround, kill spurius negative pixels, etc
        
        cutouts_ap;
        cutouts_psf;
        
        timestamps;
        
        psf; % either given or found from the image
        
        % outputs
        fluxes;
        weights;
        centroids_x;
        centroids_y;
        offsets_x;
        offsets_y;
        widths;
        backgrounds;
        
    end
    
    properties % switches/controls
        
        use_backgrounds = 0; % remove background from individual cutout
        corner_size = 0.15; % fraction of the cut_size or pixel value (must be smaller than cut_size!)
        
        use_aperture = 1;
        aperture = 5;
        annulus = 8;
        annulus_outer = []; % empty means take the rest of the cutout
        iterations = 3; % 
        
        use_gaussian = 1;
        gauss_sigma = 2;
        gauss_thresh = 1e-6;
        
        use_fitter = 0;
%         use_gaussian_psf = 1; % just use a simple PSF 
        
        bg_noise = 1;
        gain = 1;
        saturation_value = 50000;
        
        use_pixel_calibration = 0;
        pixel_cal_iterations = 3;
        
        show_bit = 0;
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        cut_size;
        
    end
    
    properties(Hidden=true)
        
        X; % output from meshgrid
        Y; % output from meshgrid
        
        cut_size_latest;
        
        fluxes_basic;
        weights_basic;
        offsets_x_basic;
        offsets_y_basic;
        centroids_x_basic;
        centroids_y_basic;
        widths_basic;
        backgrounds_basic;
        
        fluxes_ap;
        weights_ap;
        offsets_x_ap;
        offsets_y_ap;
        centroids_x_ap;
        centroids_y_ap;
        widths_ap;
        backgrounds_ap;
        
        fluxes_psf;
        weights_psf;
        offsets_x_psf;
        offsets_y_psf;
        centroids_x_psf;
        centroids_y_psf;
        widths_psf;
        backgrounds_psf;
        
        fluxes_fit;
        weights_fit;
        offsets_x_fit;
        offsets_y_fit;
        centroids_x_fit;
        centroids_y_fit;
        widths_fit;
        backgrounds_fit;
        
        default_ap_size = 3;
        
        version = 1.02;
        
    end
    
    methods % constructor
        
        function obj = Photometry(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'img.Photometry')
                if obj.debug_bit, fprintf('Photometry copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Photometry constructor v%4.2f\n', obj.version); end
                                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            obj.cutouts = [];
            obj.cutouts_proc = [];
            obj.cutouts_ap = [];
            obj.cutouts_psf = [];
            obj.timestamps = [];
            obj.psf = [];
            
            obj.fluxes = [];
            obj.weights = [];
            obj.offsets_x = [];
            obj.offsets_y = [];
            obj.widths = [];
            obj.backgrounds = [];
            
            obj.fluxes_basic = [];
            obj.weights_basic = [];
            obj.offsets_x_basic = [];
            obj.offsets_y_basic = [];
            obj.widths_basic = [];
            obj.backgrounds_basic = [];
            
            obj.fluxes_ap = [];
            obj.weights_ap = [];
            obj.offsets_x_ap = [];
            obj.offsets_y_ap = [];
            obj.widths_ap = [];
            obj.backgrounds_ap = [];
            
            obj.fluxes_psf = [];
            obj.weights_psf = [];
            obj.offsets_x_psf = [];
            obj.offsets_y_psf = [];
            obj.widths_psf = [];
            obj.backgrounds_psf = [];
            
            obj.fluxes_fit = [];
            obj.weights_fit = [];
            obj.offsets_x_fit = [];
            obj.offsets_y_fit = [];
            obj.widths_fit = [];
            obj.backgrounds_fit = [];
            
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
            input.input_var('positions', []);
            input.input_var('timestamps', [], 'times');
            input.scan_vars(varargin{:});
            
            if isa(input.cutouts, 'single')
                obj.cutouts = input.cutouts;
            else
                obj.cutouts = double(input.cutouts);
            end
            
            obj.positions = input.positions;
            obj.timestamps = input.timestamps;
            
            obj.cut_size_latest = size(input.cutouts);
            obj.cut_size_latest = obj.cut_size(1:2);
            
            obj.calcBasic;
            
            if obj.use_aperture
                obj.calcAperture;
            end
            
            if obj.use_gaussian
                obj.calcGaussian;
            end
            
            if obj.use_fitter
                obj.calcFit;
            end
            
        end
        
        function preprocess(obj)
            
            I = obj.cutouts;
            
            I(I<-3.*util.stat.std2(I)) = NaN;
            
            obj.cutouts_proc = I;
            
        end
        
        function [flux, weight, offset_x, offset_y, width, background] = calculate(obj, shape, bg_shape, iterations)

            import util.stat.sum2;
            import util.stat.median2;
            import util.text.cs;
            
            if nargin<2 || isempty(shape)
                shape = 'circle';
            end
            
            if nargin<3 || isempty(bg_shape)
                bg_shape = 'corner';
            end
            
            if cs(shape, 'none')
                make_shape = @(x,y) ones(obj.cut_size, 'like', obj.cutouts);
            elseif cs(shape, 'circle', 'aperture')
                make_shape = @(x,y) obj.makeCircle(x,y);
            elseif cs(shape, 'gaussian')
                make_shape = @(x,y) obj.makeGaussian(x,y);
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
            
            c = size(obj.cutouts); c = c(1:2);
            [obj.X,obj.Y] = meshgrid((1:c(2))-floor(c(2)/2)-1, (1:c(1))-floor(c(1)/2)-1); 
            
            if ~isempty(obj.fluxes)
                flux = obj.fluxes;
            else
                flux = NaN(size(obj.cutouts,3), size(obj.cutouts,4));    
            end
            
            if ~isempty(obj.weights)
                weight = obj.weights;
            else
                weight = NaN(size(obj.cutouts,3), size(obj.cutouts,4));    
            end
            
            if ~isempty(obj.offsets_x)
                offset_x = obj.offsets_x;
            else
                offset_x = NaN(size(obj.cutouts,3), size(obj.cutouts,4));    
            end
            
            if ~isempty(obj.offsets_y)
                offset_y = obj.offsets_y;
            else
                offset_y = NaN(size(obj.cutouts,3), size(obj.cutouts,4));    
            end
            
            if ~isempty(obj.widths)
                width = obj.widths;
            else
                width = NaN(size(obj.cutouts,3), size(obj.cutouts,4));    
            end
            
            if ~isempty(obj.backgrounds)
                background = obj.backgrounds;
            else
                background = NaN(size(obj.cutouts,3), size(obj.cutouts,4));    
            end
            
            % go over each cutout and shift it to the center
            for ii = 1:size(obj.cutouts,4) % go over different stars

                for jj = 1:size(obj.cutouts,3) % go over frames of the same star
                    
                    dx = offset_x(jj,ii);
                    if isnan(dx), dx = 0; end
                    dy = offset_y(jj,ii);
                    if isnan(dy), dy = 0; end
                    
                    for kk = 1:iterations
                        
                        I = obj.cutouts(:,:,jj,ii);
                        
                        ap = make_shape(dx,dy); % should we also let radius/psf_sigma change with the width?? 
                        
                        ap(isnan(I)) = nan; % aperture must have NaN values in the same places (for weight calculation)
                        
                        ann = make_bg_shape(dx,dy);
                        ann(isnan(I)) = nan; 
                        
                        B = median2(I.*ann); % background per pixel...
                        
                        if obj.use_backgrounds
                            I = I - B;
                        end
                        
                        S = sum2(ap);
                        ap = ap./S;
                        I = I.*ap./sum2(ap.^2);
                        m0 = sum2(I);
                        m1x = sum2(I.*obj.X)./m0;
                        m1y = sum2(I.*obj.Y)./m0;
                        m2x = sum2(I.*(obj.X-m1x).^2)./m0;
                        m2y = sum2(I.*(obj.Y-m1y).^2)./m0;
%                         mxy = sum2(I.*(obj.X-m1x).*(obj.Y-m1y))./m0;
                        
                        % quality checks:
                        if m0==0
                            m1x = 0;
                            m1y = 0;
                            m2x = NaN;
                            m2y = NaN;
%                             mxy = NaN;
                        end
                        
                        if m2x<0, m2x = NaN; end
                        if m2y<0, m2y = NaN; end
                        
                        dx = m1x;
                        dy = m1y;
                        
                        W = sqrt(mean([m2x,m2y], 'omitnan')); % should we add mxy??
                        
                    end

                    flux(jj,ii) = m0;
                    weight(jj,ii) = S;
                    offset_x(jj,ii) = dx;
                    offset_y(jj,ii) = dy;
                    width(jj,ii) = W;
                    background(jj,ii) = B;
                    
                    % add break point if dx and dy don't change much...
                    
                end

            end
            
        end
        
        function val = makeCorners(obj, ~, ~)
            
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
            
            R = obj.aperture;
            
            r = sqrt((obj.X-dx).^2+(obj.Y-dy).^2);
            val = R+0.5-r;
            val(val>1) = 1;
            val(val<0) = 0;
            
            if isa(obj.cutouts, 'single')
                val = single(val);
            end
            
            val(val==0) = NaN;
            
        end
        
        function val = makeAnnulus(obj, dx, dy)
           
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
        
            sig = obj.gauss_sigma;
            
            r = sqrt((obj.X-dx).^2+(obj.Y-dy).^2);
            
            val = exp(-0.5.*(r./sig).^2);
            
            if isa(obj.cutouts, 'single')
                val = single(val);
            end
            
            val(val<obj.gauss_thresh) = NaN;
            
        end
        
        function calcBasic(obj)
            
            import util.stat.sum2;
            
            [f,w,x,y,W,b] = obj.calculate('none', 'corner', 1);
            
            obj.fluxes_basic = f;
            obj.weights_basic = w;
            obj.offsets_x_basic = x;
            obj.offsets_y_basic = y;
            if ~isempty(obj.positions)
                obj.centroids_x_basic = x + obj.positions(:,1)';
                obj.centroids_y_basic = y + obj.positions(:,2)';
            end
            obj.widths_basic = W; 
            obj.backgrounds_basic = b;
            
            % update the newest values
            obj.fluxes = obj.fluxes_basic;
            obj.weights = obj.weights_basic;
            obj.offsets_x = obj.offsets_x_basic;
            obj.offsets_y = obj.offsets_y_basic;
            if ~isempty(obj.positions)
                obj.centroids_x = obj.centroids_x_basic;
                obj.centroids_y = obj.centroids_y_basic;
            end
            obj.widths = obj.widths_basic;
            obj.backgrounds = obj.backgrounds_basic;
            
        end
        
        function calcAperture(obj)
            
            [f,w,x,y,W,b] = obj.calculate('circle', 'annulus', obj.iterations);
            
            obj.fluxes_ap = f;
            obj.weights_ap = w;
            obj.offsets_x_ap = x;
            obj.offsets_y_ap = y;
            if ~isempty(obj.positions)
                obj.centroids_x_ap = x + obj.positions(:,1)';
                obj.centroids_y_ap = y + obj.positions(:,2)';
            end
            obj.widths_ap = W; 
            obj.backgrounds_ap = b;
            
            % update the newest values
            obj.fluxes = obj.fluxes_ap;
            obj.weights = obj.weights_ap;
            obj.offsets_x = obj.offsets_x_ap;
            obj.offsets_y = obj.offsets_y_ap;
            if ~isempty(obj.positions)
                obj.centroids_x = obj.centroids_x_ap;
                obj.centroids_y = obj.centroids_y_ap;
            end
            obj.widths = obj.widths_ap;
            obj.backgrounds = obj.backgrounds_ap;

        end
        
        function calcGaussian(obj)
           
            [f,w,x,y,W,b] = obj.calculate('gaussian', 'annulus', obj.iterations);
            
            obj.fluxes_psf = f;
            obj.weights_psf = w;
            obj.offsets_x_psf = x;
            obj.offsets_y_psf = y;
            if ~isempty(obj.positions)
                obj.centroids_x_psf = x + obj.positions(:,1)';
                obj.centroids_y_psf = y + obj.positions(:,2)';
            end
            obj.widths_psf = W;
            obj.backgrounds_psf =b;
            
            % update the newest values
            obj.fluxes = obj.fluxes_psf;
            obj.weights = obj.weights_psf;
            obj.offsets_x = obj.offsets_x_psf;
            obj.offsets_y = obj.offsets_y_psf;
            if ~isempty(obj.positions)
                obj.centroids_x = obj.centroids_x_psf;
                obj.centroids_y = obj.centroids_y_psf;
            end
            obj.widths = obj.widths_psf;
            obj.backgrounds = obj.backgrounds_psf;
            
        end
        
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
        
        function makePSF(obj) % make an image based PSF
            
            obj.psf = sum(obj.cutouts_aperture, 4, 'omitnan');
            obj.psf = obj.psf./util.stat.sum2(obj.psf);
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end
    
end

