classdef Photometry < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        aperture@img.Aperture;
        
    end
    
    properties % inputs/outputs
        
        % inputs
        cutouts;
        cutouts_proc;
        cutouts_shifted;
        cutouts_aperture;
        cut_size;
        psf; % either given or found from the image
        
        % outputs
        fluxes;
        weights;
        offsets_x;
        offsets_y;
        widths;
        backgrounds;
        
        % save the data from the entire run
%         fluxes_all;
%         weights_all;
%         offsets_x_all;
%         offsets_y_all;
%         widths_all;
%         backgrounds_all;
        
    end
    
    properties % switches/controls
        
        use_subtract_corner = 0;
        
        use_aperture = 1;
        ap_iterations = 3;
        
        use_fitter = 0;
        use_gaussian_psf = 1; % just use a simple PSF 
        
        bg_noise = 1;
        gain = 1;
        saturation_value = 50000;
        
        use_pixel_calibration = 0;
        pixel_cal_iterations = 3;
        
        show_bit = 0;
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        ap_size; % how many pixels is the aperture?
        
    end
    
    properties(Hidden=true)
        
        fluxes_raw;
        weights_raw;
        offsets_x_raw;
        offsets_y_raw;
        widths_raw;
        backgrounds_raw;
        
        fluxes_ap;
        weights_ap;
        offsets_x_ap;
        offsets_y_ap;
        widths_ap;
        backgrounds_ap;
        
        fluxes_fit;
        weights_fit;
        offsets_x_fit;
        offsets_y_fit;
        widths_fit;
        backgrounds_fit;
        
        default_ap_size = 3;
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Photometry(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'img.Photometry')
                if obj.debug_bit, fprintf('Photometry copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Photometry constructor v%4.2f\n', obj.version); end
                
                obj.aperture = img.Aperture;
                obj.aperture.plateau_size = obj.default_ap_size;
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            obj.cutouts = [];
            obj.cutouts_shifted = [];
            obj.cutouts_aperture = [];
            obj.cut_size = [];
            
            obj.fluxes = [];
            obj.weights = [];
            obj.offsets_x = [];
            obj.offsets_y = [];
            obj.widths = [];
            obj.backgrounds = [];
            
            obj.fluxes_raw = [];
            obj.weights_raw = [];
            obj.offsets_x_raw = [];
            obj.offsets_y_raw = [];
            obj.widths_raw = [];
            obj.backgrounds_raw = [];
            
            obj.fluxes_ap = [];
            obj.weights_ap = [];
            obj.offsets_x_ap = [];
            obj.offsets_y_ap = [];
            obj.widths_ap = [];
            obj.backgrounds_ap = [];
            
            obj.fluxes_fit = [];
            obj.weights_fit = [];
            obj.offsets_x_fit = [];
            obj.offsets_y_fit = [];
            obj.widths_fit = [];
            obj.backgrounds_fit = [];
            
        end
        
    end
    
    methods % getters
        
        function val = get.ap_size(obj)
            
            val = obj.aperture.plateau_size;
            
        end
        
    end
    
    methods % setters
        
        function set.ap_size(obj, val)
            
            obj.aperture.plateau_size = val;
            
        end
        
    end
    
    methods % calculations
        
        function input(obj, varargin) 
            
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('cutouts', [], 'images');
            input.scan_vars(varargin{:});
            
            obj.cutouts = double(input.cutouts);
            
            obj.cut_size = size(input.cutouts);
            obj.cut_size = obj.cut_size(1:2);
            
            obj.calcRaw;
            
            if obj.use_aperture
                obj.calcAperture;
            end
            
            if obj.use_fitter
                obj.calcFit;
            end
            
        end
        
        function calcRaw(obj)
            
            import util.stat.sum2;

            I = obj.cutouts;
            B = util.stat.corner_median(obj.cutouts);
            
            if obj.use_subtract_corner
                I = I - B;
            end
            
            I(I<-3.*util.stat.std2(I)) = NaN;
            
            obj.cutouts_proc = I;
            
            [m1x, m1y, m2x, m2y, ~] = util.img.moments(I);
            
            obj.fluxes_raw = permute(sum2(I), [3,4,2,1]);
            obj.weights_raw = permute(sum2(~isnan(I)), [3,4,2,1]);
            
            obj.offsets_x_raw = permute(m1x, [3,4,2,1]);
            obj.offsets_y_raw = permute(m1y, [3,4,2,1]);

            obj.widths_raw = permute(sqrt(m2x+m2y), [3,4,2,1]); % should we add mxy too?

            obj.backgrounds_raw = permute(B, [3,4,2,1]);
            
            obj.fluxes = obj.fluxes_raw;
            obj.weights = obj.weights_raw;
            obj.offsets_x = obj.offsets_x_raw;
            obj.offsets_y = obj.offsets_y_raw;
            obj.widths = obj.widths_raw;
            obj.backgrounds = obj.backgrounds_raw;
            
        end
        
        function calcAperture(obj)
            
            % use the first moment to adjust positions inside the aperture?
            
            import util.stat.sum2;

            obj.aperture.tile_size = obj.cut_size;
            
            obj.cutouts_shifted = zeros(size(obj.cutouts));
            
            for ii = 1:obj.ap_iterations
            
                % go over each cutout and shift it to the center
                for kk = 1:size(obj.cutouts,4) % go over different stars

                    for jj = 1:size(obj.cutouts,3) % go over frames of the same star

                        if ~isnan(obj.offsets_x(jj,kk)) && ~isnan(obj.offsets_y(jj,kk))
                            obj.cutouts_shifted(:,:,jj,kk) = util.img.imshift(obj.cutouts_proc(:,:,jj,kk), -obj.offsets_y(jj,kk), -obj.offsets_x(jj,kk));
                        else
                            obj.cutouts_shifted(:,:,jj,kk) = obj.cutouts_proc(:,:,jj,kk);
                        end
                        
                    end

                end

                % calculate new moments
                obj.cutouts_aperture = obj.cutouts_shifted.*obj.aperture.mask;
                
                I = obj.cutouts_aperture;
                
                [m1x, m1y, m2x, m2y, ~] = util.img.moments(I);
                
                % calculate new fluxes
                obj.fluxes_ap = permute(sum2(I), [3,4,2,1]);
                obj.weights_ap = repmat(obj.aperture.weight, size(obj.fluxes_ap));
                % need to think of a way to also count the NaN pixels in the weight...

                obj.offsets_x_ap = permute(m1x, [3,4,2,1]);
                obj.offsets_y_ap = permute(m1y, [3,4,2,1]);

                obj.widths_ap = permute(sqrt(m2x+m2y), [3,4,2,1]); % should we add mxy too?

                % update the newest values
                obj.fluxes = obj.fluxes_ap;
                obj.weights = obj.weights_ap;
                obj.offsets_x = obj.offsets_x_ap;
                obj.offsets_y = obj.offsets_y_ap;
                obj.widths = obj.widths_ap;
                % add annulus if you want a background estimate here... 
                
            end % go over iterations
            
        end
        
        function calcPSF(obj) % make a Gaussian PSF or an image based PSF (both normalized to unity)
            
            obj.psf = sum(obj.cutouts_aperture, 4, 'omitnan');
            obj.psf = obj.psf./util.stat.sum2(obj.psf);
            
        end
                
        function calcFit(obj)
            
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
        
    end
    
end

