classdef LimitingMagnitude < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        dir@util.sys.WorkingDirectory;
        head@head.Header;
        cat@head.Catalog;
        cal@img.Calibration; 
        
        fig@util.plot.FigHandler;
        
    end
    
    properties % inputs/outputs
        
        files = {}; % list of filenames from the dir
        file_idx = 1; % which file was loaded
        date = ''; 
        
        I_final; % the sum/image used in the final calculation
        I_removed; % image after removing all stars
        
        I_raw; % raw images 
        I_cal; % calibrated images
        im_idx = 1; % which image from the data cube was used
                
        S_raw; % raw sum of images
        S_cal; % calibrated sum of images
        num_sum; % how many images were loaded into the sum
        
        cutouts; % cutouts from the images or sum
        
        bg_map; % map of the mean background of each point in the image
        bg_var_map; % map of the average variance of each point in the image
        bg_mean; % estimate for the background mean
        bg_std; % estimate for the background std
        fwhm; % estimate the average PSF FWHM in arcsec
        
        stars; % table output of quick_find_stars()
        snr; % signal to noise ratio for each star
        mag; % magnitude from Gaia for each star
        
        fit_results; % from fitting mag log10(abs(snr))
        snr_func; % from inverting the fit coefficients to get a function from mag to snr
        
        zero_point; % the magnitude at which the star's flux is 1 count / image
        sky_mag; % magnitude of the sky background per arcsec^2
        limmag; 
        
        stack_str = ''; % a note telling if I_final was from a stack of images
        
    end
    
    properties % switches/controls
        
        use_sum = false; % if we want to do all calculations on the sum instead of individual images
        
        % parameters for quick_find_stars()
        threshold = 5; % detection threshold
        psf_sigma = 1; % gaussian width of filter
        saturation = 5e4; % for a single image
        max_mag = 18; % do not use stars fainter than this for any calculations (they are likely bad matches)
        use_find_twice = false; % run quick_find_stars again after adjusting the PSF sigma from the calculated width
        
        camera = 'Balor'; % can choose Zyla if you really wanted to... 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        expT; 
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = LimitingMagnitude(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'img.LimitingMagnitude')
                if obj.debug_bit>1, fprintf('LimitingMagnitude copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('LimitingMagnitude constructor v%4.2f\n', obj.version); end
                
                obj.dir = util.sys.WorkingDirectory;
                
                obj.cat = head.Catalog;
                obj.head = head.Header;

                obj.cal = img.Calibration; 

            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            list = properties(obj); 
            
            for ii = 1:length(list)
                
                prop = obj.(list{ii});
                
                if ~isempty(prop) && isobject(prop) && ~istable(prop) && ismethod(prop, 'reset')
                    obj.(list{ii}).reset;
                end
                
            end
            
        end
        
        function clear(obj)
            
            obj.cat.reset; 
            
            obj.date = ''; 
            obj.files = {};
            obj.file_idx = 1;

            obj.I_final = []; 
            obj.I_removed = [];
            
            obj.I_raw = [];
            obj.I_cal = [];
            obj.im_idx = 1; 

            obj.S_raw = [];
            obj.S_cal = [];
            obj.num_sum = [];

            obj.cutouts = [];

            obj.bg_map = [];
            obj.bg_var_map = []; 
            obj.bg_mean = [];
            obj.bg_std = [];
            obj.fwhm = [];

            obj.stars = [];
            obj.snr = [];
            obj.mag = [];

            obj.fit_results = []; 
            
            obj.zero_point = [];
            obj.sky_mag = []; 
            obj.limmag = [];
            
            obj.stack_str = ''; 
            
        end
        
    end
    
    methods % getters
        
        function val = get.expT(obj)
            
            if isempty(obj.head)
                val = [];
            else
                val = obj.head.EXPTIME;
            end
            
        end
        
        function val = getDescription(obj)
            
            val = sprintf('Lim-mag is %4.2f for %s%4.2fs images', obj.limmag, obj.stack_str, obj.head.EXPTIME); 
            
        end
        
    end
    
    methods % setters
        
        function set.head(obj, val)
            
            if isa(val, 'head.Header')
                
                obj.head = val;
                if ~isempty(obj.cat)
                    obj.cat.head = val; 
                end
                
            end
            
        end
        
        function set.debug_bit(obj, val)
            
            obj.debug_bit = val;
            
            list = properties(obj); 
            
            for ii = 1:length(list)
                
                prop = obj.(list{ii});
                
                if ~isempty(prop) && isobject(prop) && ~istable(prop) && isprop(prop, 'debug_bit')
                    obj.(list{ii}).debug_bit = val;
                end
                
            end
            
        end
        
    end
    
    methods % calculations
        
        function run(obj, varargin)
        
            input = util.text.InputVars;
            input.input_var('folder', obj.dir.pwd, 'directory'); 
            input.input_var('file_idx', obj.file_idx, 'file_number');
            input.input_var('im_idx', obj.im_idx, 'image_index', 'im_index', 'im_numbder', 'image_number'); 
            input.input_var('use_sum', obj.use_sum, 'sum'); 
            input.input_var('threshold', obj.threshold);
            input.input_var('psf_sigma', obj.psf_sigma, 'sigma_psf');
            input.input_var('saturation', obj.saturation);
            input.input_var('max_mag', obj.max_mag); 
            input.input_var('use_find_twice', obj.use_find_twice, 'find_twice', 'twice'); 
            input.input_var('camera', obj.camera); % can choose Zyla if you really wanted to... 
            input.input_var('debug_bit', obj.debug_bit); 
            input.scan_vars(varargin);
            
            % update the object with the input parameters
            obj.dir.cd(input.folder); 
            obj.file_idx = input.file_idx;
            obj.im_idx = input.im_idx; 
            obj.use_sum = input.use_sum; 
            obj.threshold = input.threshold; 
            obj.psf_sigma = input.psf_sigma;
            obj.saturation = input.saturation;
            obj.max_mag = input.max_mag; 
            obj.use_find_twice = input.use_find_twice;
            obj.camera = input.camera; 
            obj.debug_bit = input.debug_bit;
            
            obj.clear; 
            
            obj.files = obj.dir.match('*.h5*'); 
            
            obj.loadHeader; 
            
            obj.loadImages; 
            
            obj.findStars;
            
            obj.runAstrometry;
            
            obj.fitMagnitudesSNR;
            
            obj.calcZeroPointSkyMag;
            
            
        end
        
        function loadHeader(obj)
            
            obj.head = util.oop.load(obj.files{obj.file_idx}, 'location', '/header'); 
            
            obj.date = obj.head.RUNSTART(1:10); 
            
            if obj.debug_bit, fprintf('Loaded header for "%s" taken at %s\n', obj.head.OBJECT, obj.date); end
            
            obj.cal.loadByDate(obj.date, obj.camera); 
            
        end
        
        function loadImages(obj)
            
            f = obj.files{obj.file_idx}; % shorthand for the filename
            
            if obj.use_sum % prefer to load the summed images
                
                try 
                    obj.S_raw = h5read(f, '/stack'); 
                    obj.num_sum = h5readatt(f, '/stack', 'num_sum'); 
                catch
                    
                    if obj.debug_bit, fprintf('Could not find stack images, generating sum from image data...\n'); end
                    
                    obj.I_raw = h5read(f, '/images'); 
                    obj.S_raw = single(sum(obj.I_raw,3)); 
                    obj.num_sum = size(obj.I_raw,3); 
                    
                end
                
                obj.S_cal = obj.cal.input(obj.S_raw, 'sum', obj.num_sum); 
                
                obj.I_final = obj.S_cal; 
                obj.stack_str = sprintf('a %d stack of ', obj.num_sum); 
                
            else
                
                try
                    
                    obj.I_raw = h5read(f, '/images'); 
                    obj.I_raw = obj.I_raw(:,:,obj.im_idx); 
                    obj.I_cal = obj.cal.input(obj.I_raw); 
                    obj.num_sum = 1; 
                    
                    obj.I_final = obj.I_cal; 
                    
                catch
                    
                    if obj.debug_bit, fprintf('Could not find raw images, using stack instead...\n'); end
                                        
                    obj.S_raw = h5read(f, '/stack'); 
                    obj.num_sum = h5readatt(f, '/stack', 'num_sum'); 
                    
                    obj.S_cal = obj.cal.input(obj.S_raw, 'sum', obj.num_sum); 
                    
                    obj.I_final = obj.S_cal; 
                    obj.stack_str = sprintf('a %d stack of ', obj.num_sum); 
                    
                end
                
            end
            
            if obj.debug_bit, fprintf('Loaded %d images\n', obj.num_sum); end
            
        end
        
        function findStars(obj)
            
            [obj.bg_map, obj.bg_var_map] = util.img.im_stats(obj.I_final, 'tile', 250, 'overlap', 100, 'method', 'median', 'output', 'map');

            obj.bg_mean = util.stat.median2(obj.bg_map);
            obj.bg_std = sqrt(util.stat.median2(obj.bg_var_map)); 
            
            Nstars = 100; 
            
            psf_width = obj.psf_sigma; 
            I = obj.I_final - obj.bg_map; 
            
            if obj.use_find_twice
    
                obj.stars = util.img.quick_find_stars(I, 'number', Nstars, ... % only look for the brightest stars at this point
                    'thresh', obj.threshold, 'saturation', obj.saturation.*obj.num_sum, ...
                    'psf', psf_width, 'unflagged', 1, ...
                    'mean', 0, 'std', sqrt(obj.bg_var_map)); % the background is estimated for all locations
                
                Nstars = min(Nstars, height(obj.stars)); 
                
                [obj.cutouts, obj.I_removed] = util.img.mexCutout(I, obj.stars.pos(1:Nstars,:), 15, 0, NaN); % should replace 100 with a parameter (and test it isn't bigger than the number of stars!)

                w = obj.head.SCALE.*squeeze(util.img.fwhm(obj.cutouts, 'method', 'filters',...
                    'defocus', 1, 'generalized', 5, 'min', 0.25, 'step', 0.25, 'max', 10)); 

                obj.fwhm = util.vec.weighted_average(w', obj.stars.snr(1:length(w)), 1); 
                
                psf_width = obj.fwhm/2.355/obj.head.SCALE; 
                
                if obj.debug_bit, fprintf('Found %d stars, FWHM= %4.2f", background mean/std= %4.2f/%4.2f\n', ...
                        height(obj.stars), obj.fwhm, obj.bg_mean, obj.bg_std); end 
            end
                            
            obj.stars = util.img.quick_find_stars(I, ...
                'thresh', obj.threshold, 'saturation', obj.saturation.*obj.num_sum, ...
                'psf', psf_width, 'unflagged', 1, ...
                'mean', 0, 'std', sqrt(obj.bg_var_map)); % the background is estimated for all locations

            [obj.cutouts, obj.I_removed] = util.img.mexCutout(I, obj.stars.pos(1:100,:), 15); % should replace 100 with a parameter (and test it isn't bigger than the number of stars!)

            w = obj.head.SCALE.*squeeze(util.img.fwhm(obj.cutouts, 'method', 'filters',...
                'defocus', 1, 'generalized', 5, 'min', 0.25, 'step', 0.25, 'max', 10)); 

            obj.fwhm = util.vec.weighted_average(w', obj.stars.snr(1:length(w)), 1); 

            if obj.debug_bit, fprintf('Found %d stars, FWHM= %4.2f", background mean/std= %4.2f/%4.2f\n', ...
                    height(obj.stars), obj.fwhm, obj.bg_mean, obj.bg_std); end 

            
            
        end
        
        function runAstrometry(obj)
            
            obj.cat.input(obj.stars); 
            
            obj.mag = obj.cat.magnitudes;
            obj.snr = obj.stars.flux; 
            
        end
        
        function calcZeroPointSkyMag(obj)
            
            % get the zero point
            f = obj.stars.flux;
            dm = obj.mag + 2.5.*log10(f); % diffence between instrumental magnitude and Gaia magnitude
            obj.zero_point = util.vec.weighted_average(dm, sqrt(abs(f))); 

            % get the sky magnitude
            B = obj.bg_mean./obj.head.SCALE.^2; % convert to background flux per arcsecond^2
            
            obj.sky_mag = obj.zero_point - 2.5*log10(B); 
           
            if obj.debug_bit, fprintf('Zero point= %4.2f, sky mag= %4.2f\n', obj.zero_point, obj.sky_mag); end
            
        end
        
        function fitMagnitudesSNR(obj)
            
            m = obj.mag;
            s = log10(abs(obj.snr)); 
            
            % remove spurious detections
            m(obj.mag>obj.max_mag | obj.snr<obj.threshold) = []; 
            s(obj.mag>obj.max_mag | obj.snr<obj.threshold) = []; 
            
            obj.fit_results = util.fit.polyfit(s, m, 'order', 2, 'double', 1, 'sigma', 2.5, 'iterations', 3, 'var', nanmax(s)./s); 
            
            obj.limmag = obj.fit_results.func(log10(obj.threshold)); 
            
            if obj.debug_bit, fprintf('Limiting magnitude is %4.2f\n', obj.limmag); end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('parent', [], 'figure'); 
            input.input_var('font_size', 20); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.parent)
                input.parent = gcf;
            end
            
            delete(input.parent.Children); 
            
            ax1 = axes('Parent', input.parent, 'Position', [0.05, 0.15, 0.85, 0.8]); 
            
            obj.showMagnitudes('ax', ax1, 'font_size', input.font_size); 
            
            ax1.YAxisLocation = 'right'; 
            
            util.plot.inner_title(obj.getDescription, 'Position', 'NorthWest', ...
                'ax', ax1, 'FontSize', input.font_size); 
            
            ax2 = axes('Parent', input.parent, 'Position', [0.15 0.45 0.4 0.4]); 
            
            obj.showFit('ax', ax2, 'font_size', input.font_size); 
            
            
            
        end
        
        function showFit(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis'); 
            input.input_var('font_size', 18); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            m_min = nanmin(obj.mag);
            s_max = nanmax(obj.snr);
            
            s = obj.threshold:0.01:s_max; 
            m = obj.fit_results.func(log10(s)); 
            
            s(m<m_min | m>obj.max_mag) = [];
            m(m<m_min | m>obj.max_mag) = [];
            
            semilogy(input.ax, obj.mag(1:1:end), obj.snr(1:1:end), '.', ...
                m, s, '-', 'LineWidth', 2, 'MarkerSize', 10); 
            
            hold(input.ax, 'on'); 
            
            plot(input.ax, [m_min, obj.max_mag], obj.threshold.*[1 1], 'g--', 'HandleVisibility', 'off'); 
            text(input.ax, double(m_min)+0.2, double(obj.threshold)*2, ...
                sprintf('threshold= %4.1f \\sigma', obj.threshold), ...
                'Color', 'g', 'FontSize', input.font_size-4); 
            
            plot(input.ax, obj.limmag*[1 1], [1, 1e10], 'g--', 'HandleVisibility', 'off'); 
            text(input.ax, double(obj.limmag+0.5), double(s_max), ...
                sprintf('lim-mag= %4.1f', obj.limmag), ...
                'Color', 'g', 'FontSize', input.font_size-4, 'Rotation', 270); 
            
            hold(input.ax, 'off'); 
            
            xlabel(input.ax, 'Gaia Mag Bp'); 
            ylabel(input.ax, 'detection S/N'); 
            
            input.ax.XLim = [m_min, obj.max_mag]; 
            input.ax.YLim = [1, s_max*2]; 
            input.ax.YTick = 10.^(1:2:log10(s_max*2)); 
            input.ax.FontSize = input.font_size; 
    
        end
        
        function showMagnitudes(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis'); 
            input.input_var('font_size', 18); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            h = histogram(input.ax, obj.mag, 'BinWidth', 0.1); 
            
            hold(input.ax, 'on'); 
            
            plot(input.ax, obj.limmag.*[1 1], [0 1e10], 'g--'); 
            text(input.ax, double(obj.limmag+0.4), max(h.Values), ...
                sprintf('limiting mag= %4.1f', obj.limmag), ...
                'Color', 'g', 'FontSize', input.font_size, 'Rotation', 270); 
            
            input.ax.XLim = [4 obj.max_mag]; 
            input.ax.YLim = [0 max(h.Values)+10]; 
            
            xlabel(input.ax, 'Gaia Mag Bp'); 
            ylabel(input.ax, 'Number of stars'); 

            input.ax.FontSize = input.font_size; 
            
            hold(input.ax, 'off'); 
            
        end
        
    end    
    
end

