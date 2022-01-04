classdef Lightcurves < handle
% Stores and displays photometric results over the entire run. 
%
% The first job of this class is to store fluxes, backgrounds, etc, from 
% photometry on multiple batches. Each time the batch is done, use the input
% metohd to add the new data to this object. 
%
% Example (assume "phot" is an img.Photometry object with new data):
% >> lightcurves.input(phot.fluxes, phot.erros, timestamps, phot.areas, ...
% You can input the data by order, or by keyword-value pairs. 
% Finally, you can just use getData to extract everything from "phot":
% >> lightcurves.getData(phot);
%
% The new data must expand the existing matrices, that can become quite big. 
% To save some time on data reallocation, you can use "use_double_up", which 
% makes the hidden, underlying matrices double their allocated size using 
% NaN padding, and only displaying the part of the matrix that is really 
% filled up. To this end the class keeps "frame_index" and "num_frames". 
% The storage is doubled up when needed, and in between doubling no allocation 
% is required. The data is kept in "fluxes_full" and other such matrices. 
% 
% NOTE: the photometry object usually also provides a "phot_pars" object 
%       that contains some of the crucial information on how the photometry
%       was actually performed: aperture sizes, gaussian sigma, etc. 
%       Make sure you use this info when documenting the results!
%
% The data types that are collected are:
%  timestamps, fluxes, errors, areas, backgrounds, variances, offsets_x, 
%  offsets_y, centroids_x, centroids_y, widths, bad_pixels, flags. 
%  (see description of util.img.photometry2 for more info). 
% The timestamps are one vector, but the others are 2D or 3D matrices. 
% The 2nd dim is for stars, the 3rd dim is for different apertures. 
% Choose the aperture you want using "index_flux". The default is "end", 
% which specifies using the forced-aperture with the largest radius. 
% 
% The second job of this class is to clean up the fluxes. 
% The raw fluxes are transformed using a few cleanup methods (all optional) 
% in a few stages: 
%  -fluxes_sub: remove the background (see "use_subtract_backgrounds"). 
%  -fluxes_rem: remove outliers (see "use_skip_flagged", "use_remove_outliers", 
%               and "use_skip_bad_times"). 
%  -fluxes_cal: apply "zero point correction" type self-calibration. 
%
% NOTE: the derived fluxes are all calculated once and lazy loaded. This is
%       important to keep in mind also when adding new processing methods, 
%       that also you must add a setter for the new control parameters, so
%       a change in the controls clears the saved fluxes. 
%
% The third job of this class is to calculate the power spectrum and a poly 
% fit to the data. These results can later be used for more advanced 
% statistical analysis. These are all lazy loaded as well. 
%
% The fourth job of this class is to calculate statistics on each star. 
% <This part needs to be expanded once we write the code!> 
%
% The last job of this class is to display the data in a useful way. 
% *You can specify how many stars to show, or show a list of stars, using 
%  "show_num_stars" or "show_indices" (which overrules the latter). 
% *Use "show_what" to change the display from fluxes to other parameters. 
%  You can choose the "power spectrum" to show the flux PS in frequency space. 
% *Use "show_flux_type" to choose which processing of the flux should be used. 
%  You can choose "all" to show raw, rem and cal fluxes superimposed. 
% *Use "show_for" to switch from "time" (default) which shows data over time, 
%  to "stars" that shows data for different stars. 
% *Use "use_show_log" to display the data in log scale. For most 
%  results this only turns the Y axis to log. Power spectra is always shown 
%  in log Y scale, but then this switch controls the X axis scale. 
% *Use "use_smooth" and "smooth_interval" to show the results binned over 
%  some number of samples. This declutters the data and shows long term trends, 
%  at the expense of diluting large transient effects. This is display only! 

    properties(Transient=true)
        
        gui@img.gui.LightGUI; % graphic user interface
        
        figures_spawned = {};
        
    end
    
    properties % objects
        
        head@head.Header; % keep a copy of the header metadata
        cat@head.Catalog; % keep a copy of the star catalog
        
        phot_pars; % a struct with some housekeeping about how the photometry was done
        
    end
    
    properties % inputs/outputs
        
        frame_index = 1;
        
        batch_timestamps; % timestamp for the middle of each batch
        batch_midframe; % frame number for the central frame of each batch (for interpolating back onto the full frame rate)
        fwhm_coeffs_pix; % coefficients needed to calculate the FWHM (in pixels!)
        fwhm_x_center; % must subtract this from the x positions when using the polynomial coeffs
        fwhm_y_center; % must subtract this from the y positions when using the polynomial coeffs
        
        magnitudes; % keep a copy of the star magnitudes for doing statistics
        colors; % the Bp-Rp color of each star
        
        zp_fit_results; % the fit results for the spatial fit to the zero point
        mean_zp; % the average ZP for each frame
        
    end
    
    properties % switches/controls
        
        % these controls are used for inputting the data:
        use_double_up = 0; % choose if you want to expand the data storage by factor of 2 each time when space runs out... 
        use_preallocate = 1; % if this is used, we preallocate the entire storage at the beginning (must call startup(N) with the number of frames)
        
        % processing steps for fluxes_sub:
        use_subtract_backgrounds = 1; 
        use_background_median = 1;
        use_background_fit = 0; % fit the background to a 2nd order polynomial for each frame
        
        % processing steps for fluxes_rem:
        minimal_length = 10; % minimal number of frames to run any analysis more comlicated than b/g removal
        use_skip_flagged = 0; % skip samples that have a bad photometry flag (for forced, this may be overkill)
        use_remove_outliers = 1; % remove outliers with "N sigmas" above the noise, calculated by polyfit
        outlier_sigma = 5; % this is NOT used by polyfit, only for making fluxes_sub
        use_skip_bad_times = 0; % if a large fraction of stars have NaN in this epoch, set all to NaN in this epoch
        bad_times_fraction = 0.1; % what fraction of stars need to be NaN to be considered bad times
        
        % processing steps for fluxes_cal:
        use_airmass_correction = 1;
        use_color_correction  = 1; 
        use_psf_correction = 0;
        use_zero_point = 1;
        zero_point_spatial_order = 2; % fit a polynomial to the ZP
        use_offset_fit = 1; 
        use_sysrem = 0;
        sysrem_iterations = 1;
        use_self_exclude = 0;
        
        missing_method = 'linear'; % can also choose 'previous', 'next', 'nearest', 'linear', 'spline', 'pchip', 'makima' (see the help section of fillmissing())
%         use_savitzky_golay = 1;
%         sg_order = 3;
%         sg_length = 125; % take 5 seconds as long enough to smooth over
        
        index_flux = 'end'; % which aperture to use, default "end" is for forced photometery with the largest aperture
        
        use_welch = 0;
        
        sampling_jump = 2; % how many data points to bin together when calculating RE on larger and larger time-intervals
        num_points_rms = 20; % how many data points of binned data to use for calculating local relative error
        
        show_for = 'time'; % can also choose "stars" to get star by star stats
        show_what = 'flux'; % can choose "flux", "areas", "backgrounds", "variances", "centroids", "offsets", "widths", "power spectra"
        show_flux_type = 'raw'; % can choose "raw" or "sub" or "rem" or "cal" or "all".
        
        show_num_stars = 10; % up to this number of stars are shown.
        show_indices = []; % if this is not empty, it is used as a list of stars to show
        
        use_show_datetime = 1; % instead of timestamps, if there are juldates available
        
        use_smooth = 1; % smoothing is applied to the plotted data only! 
        smooth_interval = 10; % how many samples to average over when smoothing
        
        use_show_log = 1; % show the plot in log scale. Power spectra always shows log Y, but can toggle X. For others, toggles log Y only (X is linear). 
        
%         use_show_re_bins = 0;
        show_mag_limit = 15; 
        
        use_show_bin_fits = 1;
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true, Transient=true)
        
        timestamps;
        juldates;
        
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
        
        fluxes_sub; % remove the measured background
        fluxes_rem; % remove outliers using a few methods
        fluxes_cal; % full calibration using airmass, PSF width, Sys-Rem, etc
        
        fit_results;
        
        power_spectra; 
        psd_freq;
        
        total_RE; % relative error calculated on the entire length of the data
        local_RE; % relative error calculated on a short interval and averaged over all intervals
        bin_widths_seconds; % the time-interval used for each calculation of RE, starting from the native time sampling
        
        RE_fit_results;
        bad_star_indices;
        
    end
    
    properties(Hidden=true)
        
        width_vector; % save a copy of the flux-weighted mean width of all stars. 
        airmass_vector; % save a copy of the airmass 
        
        % These full datasets can include a lot of extra allocated memory, 
        % filled with NaN, to allow less allocations of the data. 
        % To use the actual data call the same names without "_full". 
        timestamps_full;
        timestamps_full_original; % if we use repair 
        
        juldates_full;
        juldates_full_original;
        
        fluxes_full;
        errors_full;
        areas_full;
        backgrounds_full;
        variances_full;
        centroids_x_full;
        centroids_y_full;
        offsets_x_full;
        offsets_y_full;
        widths_full;
        bad_pixels_full;
        flags_full;
        
        show_what_list = {'fluxes', 'areas', 'backgrounds', 'variances', 'centroids', 'offsets', 'widths', 'bad_pixels', 'power spectra', 'statistics', 'binning'};
        show_flux_type_list = {'raw', 'sub', 'rem', 'cal'};
        
        sysrem_a_vec;
        sysrem_c_vec; 
        
        last_shown = '';
        last_xlim = [];
        last_ylim = [];
        
        version = 1.10;
        
    end
    
    properties(Hidden=true, Transient=true)
        
        % These are lazy loaded from XXX_full and are cleared whenever the 
        % underlying matrices are changed. This saves time as just slicing
        % into the XXX_full matrices is slow when they are large. 
        % These are not saved if the object is dumped to memory. 
        timestamps_;
        juldates_;
        
        fluxes_;
        errors_;
        areas_;
        backgrounds_;
        variances_;
        centroids_x_;
        centroids_y_;
        offsets_x_;
        offsets_y_;
        widths_;
        bad_pixels_;
        flags_;
        
        fluxes_sub_;
        fluxes_rem_;
        fluxes_cal_;
        
        power_spectra_; 
        psd_freq_;
        
        fit_results_; % save this from the polyfit (lazy load)
        
        total_RE_;
        local_RE_;
        bin_widths_seconds_;
        
        RE_fit_results_;
        bad_star_indices_;
        
        
        
    end
    
    methods % constructor
        
        function obj = Lightcurves(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'Lightcurves')
                if obj.debug_bit>1, fprintf('Lightcurves copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('Lightcurves constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj) % remove the long term storage for the entire run
            
            obj.timestamps_full = [];
            obj.timestamps_full_original = [];
            obj.juldates_full = [];
            obj.juldates_full_original = [];
            obj.magnitudes = [];
            
            obj.fluxes_full = [];
            obj.errors_full = [];
            obj.areas_full = [];
            obj.backgrounds_full = [];
            obj.variances_full = [];
            obj.offsets_x_full = [];
            obj.offsets_y_full = [];
            obj.centroids_x_full = [];
            obj.centroids_y_full = [];
            obj.widths_full = [];
            obj.bad_pixels_full = [];
            obj.flags_full = [];
            
            obj.frame_index = 1;
            
            obj.batch_timestamps = [];
            obj.batch_midframe = [];
            obj.fwhm_coeffs_pix = [];
            obj.fwhm_x_center = [];
            obj.fwhm_y_center = [];
        
            obj.clearIntermidiate;
            obj.clearFits;
            obj.clearFluxes;
            obj.clearPSD;
            obj.clearRE;
            obj.clearREfit;
            
        end
        
        function clearIntermidiate(obj) % clear all the data sliced from the full matrices
            
            obj.timestamps_ = [];
            obj.juldates_ = [];
            
            obj.fluxes_ = [];
            obj.fluxes_sub_ = [];
            obj.fluxes_cal_ = [];
            obj.errors_ = [];
            obj.areas_ = [];
            obj.backgrounds_ = [];
            obj.variances_ = [];
            obj.centroids_x_ = [];
            obj.centroids_y_ = [];
            obj.offsets_x_ = [];
            obj.offsets_y_ = [];
            obj.widths_ = [];
            obj.bad_pixels_ = [];
            obj.flags_ = [];
            
        end
        
        function clearFits(obj) % clear the polyfit results
            
            obj.fit_results_ = [];
            
        end
        
        function clearFluxes(obj) % clear the process fluxes (sub, rem, cal)

            obj.fluxes_sub_ = [];
            obj.fluxes_rem_ = [];
            obj.fluxes_cal_ = [];
            
            obj.clearREfit;
            
        end
        
        function clearPSD(obj) % clear the power spectrum and frequency axis

            obj.power_spectra_ = []; 
            obj.psd_freq_ = [];
            
        end
        
        function clearRE(obj)
           
            obj.total_RE_ = [];
            obj.local_RE_ = [];
            obj.bin_widths_seconds_ = [];
            
        end
        
        function clearREfit(obj)
            
            obj.RE_fit_results = [];
            obj.bad_star_indices = [];
            
        end
        
    end
    
    methods % getters
        
        function val = index_flux_number(obj)
            
            if isempty(obj.index_flux)
                val = size(obj.fluxes_full,3);
            elseif ischar(obj.index_flux) && strcmpi(obj.index_flux, 'end')
                val = size(obj.fluxes_full,3);
            elseif ischar(obj.index_flux) && strcmpi(obj.index_flux, 'first')
                val = 1;
            elseif isnumeric(obj.index_flux)
                val = obj.index_flux;
            else
                val = 1;
            end
            
        end
        
        function val = num_frames(obj)
            
            val = obj.frame_index - 1;
            
        end
        
        function val = get.timestamps(obj)
            
            if isempty(obj.timestamps_full)
                val = [];
                return;
            end
            
            if isempty(obj.timestamps_)
                obj.timestamps_ = obj.timestamps_full(1:obj.frame_index-1);
            end
            
            val = obj.timestamps_;
            
            if all(isnan(val))
                val = 1:obj.num_frames;
                return;
            end
            
        end
        
        function val = get.juldates(obj)
            
            if isempty(obj.juldates_full)
                val = [];
                return;
            end
            
            if isempty(obj.juldates_)
                obj.juldates_ = obj.juldates_full(1:obj.frame_index-1);
            end
            
            val = obj.juldates_;
                        
        end
        
        function val = get.fluxes(obj)
            
            if isempty(obj.fluxes_full)
                val = [];
                return;
            end
            
            if isempty(obj.fluxes_)
%                 disp('slicing fluxes'); 
                obj.fluxes_ = obj.fluxes_full(1:obj.frame_index-1,:,:);
            end
            
            val = obj.fluxes_;
            
            if all(isnan(val))
                val = [];
                return;
            end
            
        end
        
        function val = get.fluxes_sub(obj)
            
            if isempty(obj.fluxes)
                val = [];
            else
                
                if isempty(obj.fluxes_sub_)
                    
                    t_sub = tic;
                    
                    obj.fluxes_sub_ = obj.fluxes(:,:,obj.index_flux_number); 

                    if obj.use_subtract_backgrounds
                        
                        a = obj.areas(:,:,obj.index_flux_number); 
                        b = NaN(size(a), 'like', a); 
                        
                        if obj.use_background_fit
                            for ii = 1:size(a,1) % number of frames
                                x = obj.centroids_x(ii,:,obj.index_flux_number); 
                                y = obj.centroids_y(ii,:,obj.index_flux_number); 
                                
                                x = (x-nanmean(x))./nanstd(x); 
                                y = (y-nanmean(y))./nanstd(y); 
                                fr = util.fit.surf_poly(x,y,obj.backgrounds(ii,:,obj.index_flux_number), ...
                                    'order', 2, 'iterations', 2, 'double', 1);
                                b(ii,:) = fr.vm'; % the model values are used for each star. 
                            end
                        elseif obj.use_background_median
                            b = nanmedian(obj.backgrounds(:,:,obj.index_flux_number),2); 
                        else
                            b = obj.backgrounds(:,:,obj.index_flux_number); 
                        end
                        
                        obj.fluxes_sub_ = obj.fluxes_sub_ - a.*b;
                        
                        obj.fluxes_sub_(a==0) = 0; 
                        
                    end
                    
                    if obj.debug_bit, fprintf('Time to get fluxes_sub was %f seconds\n', toc(t_sub)); end
                
                end
                
                val = obj.fluxes_sub_;
                
            end
             
            if all(isnan(val))
                val = [];
                return;
            end
            
        end
        
        function val = get.fluxes_rem(obj)
           
            if isempty(obj.fluxes)
                val = [];
            elseif size(obj.fluxes_sub, 1)<obj.minimal_length
                val = obj.fluxes_sub;
            else
                
                if isempty(obj.fluxes_rem_)
                    
                    t_rem = tic;
                    
                    obj.fluxes_rem_ = obj.fluxes_sub;
                    
                    if obj.use_skip_flagged
                        
                        obj.fluxes_rem_(logical(obj.flags(:,:,obj.index_flux_number))) = NaN;
                        
                    end
                    
                    if obj.use_remove_outliers
                        
                        r = obj.fit_results; % polyfit to a 5th order polynomial with 3 sigma outlier rejection
                        
                        model = [r.ym];
                        
                        residuals = obj.fluxes_sub-model; 
                        
                        sigma = sqrt([r.variance]); 
                        
                        idx = abs(residuals)>obj.outlier_sigma.*sigma; 
                        
                        obj.fluxes_rem_(idx) = NaN;
                        
                    end
                    
                    if obj.use_skip_bad_times
                        
                        count = sum(isnan(obj.fluxes_rem_),2); % sum the number of NaNs in each epoch
                        
                        idx = (count./size(obj.fluxes_rem_,2))>obj.bad_times_fraction; 
                        
                        obj.fluxes_rem_(idx,:) = NaN;
                        
                    end
                    
                    if obj.debug_bit, fprintf('Time to get fluxes_rem was %f seconds\n', toc(t_rem)); end
                
                end
                
                val = obj.fluxes_rem_;
                
            end
            
            if all(isnan(val))
                val = [];
                return;
            end
            
            
        end
        
        function val = get.fluxes_cal(obj)

            if isempty(obj.fluxes)
                val = [];
            elseif size(obj.fluxes_sub, 1)<obj.minimal_length
                val = obj.fluxes_sub;
            else
            
                if isempty(obj.fluxes_cal_)
                    
                    t_cal = tic;
                    
                    obj.fluxes_cal_ = obj.calculateCalibration(obj.fluxes_rem); 
                    
                    if obj.debug_bit, fprintf('Time to get fluxes_cal was %f seconds\n', toc(t_cal)); end
                    
                end
                    
                val = obj.fluxes_cal_;
                
            end
            
        end
        
        function val = get.errors(obj)
            
            if isempty(obj.errors_full)
                val = [];
                return;
            end
            
            if isempty(obj.errors_)
                obj.errors_ = obj.errors_full(1:obj.frame_index-1,:,:);
            end
            
            val = obj.errors_; 
            
            if all(isnan(val))
                val = [];
                return;
            end
            
        end
        
        function val = get.areas(obj)
            
            if isempty(obj.areas_full)
                val = [];
                return;
            end
            
            if isempty(obj.areas_)
                obj.areas_ = obj.areas_full(1:obj.frame_index-1,:,:);
            end
            
            val = obj.areas_;
            
            if all(isnan(val))
                val = [];
            end
            
        end
        
        function val = get.backgrounds(obj)
            
            if isempty(obj.backgrounds_full)
                val = [];
                return;
            end
            
            if isempty(obj.backgrounds_)
                obj.backgrounds_ = obj.backgrounds_full(1:obj.frame_index-1,:,:);
            end
            
            val = obj.backgrounds_;
            
            if all(isnan(val))
                val = [];
                return;
            end
            
        end
        
        function val = get.variances(obj)
            
            if isempty(obj.variances_full)
                val = [];
                return;
            end
            
            if isempty(obj.variances_)
                obj.variances_ = obj.variances_full(1:obj.frame_index-1,:,:);
            end
            
            val = obj.variances_;
            
            if all(isnan(val))
                val = [];
                return;
            end
            
        end
        
        function val = get.offsets_x(obj)
                        
            if isempty(obj.offsets_x_full)
                val = [];
                return;
            end
            
            if isempty(obj.offsets_x_)
                obj.offsets_x_ = obj.offsets_x_full(1:obj.frame_index-1,:,:);
            end
            
            val = obj.offsets_x_;
            
            if all(isnan(val))
                val = [];
                return;
            end
            
        end
        
        function val = get.offsets_y(obj)
            
            if isempty(obj.offsets_y_full)
                val = [];
                return;
            end
            
            if isempty(obj.offsets_y_)
                obj.offsets_y_ = obj.offsets_y_full(1:obj.frame_index-1,:,:);
            end
            
            val = obj.offsets_y_;
            
            if all(isnan(val))
                val = [];
                return;
            end
            
        end
        
        function val = get.centroids_x(obj)
            
            if isempty(obj.centroids_x_full)
                val = [];
                return;
            end
            
            if isempty(obj.centroids_x_)
                obj.centroids_x_ = obj.centroids_x_full(1:obj.frame_index-1,:,:);
            end
            
            val = obj.centroids_x_;
            
            if all(isnan(val))
                val = [];
            end
            
        end
        
        function val = get.centroids_y(obj)
            
            if isempty(obj.centroids_y_full)
                val = [];
                return;
            end
            
            if isempty(obj.centroids_y_)
                obj.centroids_y_ = obj.centroids_y_full(1:obj.frame_index-1,:,:);
            end
            
            val = obj.centroids_y_;
            
            if all(isnan(val))
                val = [];
                return;
            end
            
        end
        
        function val = get.widths(obj)
                        
            if isempty(obj.widths_full)
                val = [];
                return;
            end
            
            if isempty(obj.widths_)
                obj.widths_ = obj.widths_full(1:obj.frame_index-1,:,:);
            end
            
            val = obj.widths_;
            
            if all(isnan(val))
                val = [];
                return;
            end
            
        end
        
        function val = get.bad_pixels(obj)
                        
            if isempty(obj.bad_pixels_full)
                val = [];
                return;
            end
            
            if isempty(obj.bad_pixels_)
                obj.bad_pixels_ = obj.bad_pixels_full(1:obj.frame_index-1,:,:);
            end
            
            val = obj.bad_pixels_;
            
            if all(isnan(val))
                val = [];
            end
            
        end
        
        function val = get.flags(obj)
            
            if isempty(obj.flags_full)
                val = [];
                return;
            end
            
            if isempty(obj.flags_)
                obj.flags_ = obj.flags_full(1:obj.frame_index-1,:,:);
            end
            
            val = obj.flags_;
            
            if all(isnan(val))
                val = [];
                return;
            end
            
        end
        
        function val = get.fit_results(obj)
            
            if size(obj.fluxes_sub, 1)<obj.minimal_length
                val = [];
            else
                if isempty(obj.fit_results_)
                    obj.fit_results_ = util.fit.polyfit(obj.timestamps, abs(obj.fluxes_sub), 'order', 5, 'sigma', 3, 'double', 1); % must use double precision for the sizes and length of these fluxes!
                end
            
                val = obj.fit_results_;
            end
            
        end
        
        function val = get.power_spectra(obj)
            
            if isempty(obj.fluxes_sub) || size(obj.fluxes_sub, 1)<obj.minimal_length
                val = [];
                return; 
            end
            
            if isempty(obj.power_spectra_) || isempty(obj.psd_freq_)
%                 disp('calculating power spectrum (from get.power_spectra)'); 
                [obj.power_spectra_, obj.psd_freq_] = obj.calculatePSD;
            end
            
            val = obj.power_spectra_; 
            
        end
        
        function val = get.psd_freq(obj)
            
            if isempty(obj.fluxes_sub) || size(obj.fluxes_sub, 1)<obj.minimal_length
                val = [];
                return; 
            end
            
            if isempty(obj.power_spectra_) || isempty(obj.psd_freq_)
                disp('calculating power spectrum (from get.power_spectra)'); 
                [obj.power_spectra_, obj.psd_freq_] = obj.calculatePSD; 
            end
            
            val = obj.psd_freq_; 
            
        end
        
        function val = get.local_RE(obj)
            
            if isempty(obj.local_RE_) || size(obj.fluxes_sub, 1)<obj.minimal_length
                obj.calculateRelativeError;
            end
            
            val = obj.local_RE_;
            
        end
        
        function val = get.total_RE(obj)
            
            if isempty(obj.total_RE_) || size(obj.fluxes_sub, 1)<obj.minimal_length
                obj.calculateRelativeError;
            end
            
            val = obj.total_RE_;
            
        end
        
        function val = get.bin_widths_seconds(obj)
            
            if isempty(obj.bin_widths_seconds_) || size(obj.fluxes_sub, 1)<obj.minimal_length
                obj.calculateRelativeError;
            end
            
            val = obj.bin_widths_seconds_;
            
        end
        
        function str_out = print_pars(obj)

            import util.text.print_vec;

            if isempty(obj.phot_pars)
                str = sprintf('Nframes: %d ', obj.num_frames); 
            else

%                 str = sprintf('Nframes: %d | Nstars: %d | Cutous: %dx%d ', obj.num_frames, obj.phot_pars.cutout_size(4), obj.phot_pars.cutout_size(1), obj.phot_pars.cutout_size(2));
                str = sprintf('Nframes: %d | Nstars: %d | Cutous: %dx%d ', obj.num_frames, size(obj.fluxes_full,2), obj.phot_pars.cutout_size(1), obj.phot_pars.cutout_size(2));

                if isfield(obj.phot_pars, 'aperture_radius') && ~isempty(obj.phot_pars.aperture_radius)
                    str = sprintf('%s | ap: %s', str, print_vec(obj.phot_pars.aperture_radius));
                end

                if isfield(obj.phot_pars, 'forced_radius') && ~isempty(obj.phot_pars.forced_radius)
                    str = sprintf('%s | f.ap: %s', str, print_vec(obj.phot_pars.forced_radius)); 
                end

                if isfield(obj.phot_pars, 'gauss_sigma') && ~isempty(obj.phot_pars.gauss_sigma)
                    str = sprintf('%s | gauss: %4.2f', str, obj.phot_pars.gauss_sigma);
                end

                if isfield(obj.phot_pars, 'annulus_radii') && ~isempty(obj.phot_pars.annulus_radii)
                    str = sprintf('%s | ann: %s', str, print_vec(obj.phot_pars.annulus_radii)); 
                end

                if isfield(obj.phot_pars, 'iterations') && ~isempty(obj.phot_pars.iterations)
                    str = sprintf('%s | iter: %d', str, obj.phot_pars.iterations);
                end

                if isfield(obj.phot_pars, 'shift_resolution') && ~isempty(obj.phot_pars.shift_resolution)
                    str = sprintf('%s | shift res: %4.2f', str, obj.phot_pars.shift_resolution); 
                end

            end
            
            if nargout>0
                str_out = str;
            else
                disp(str); 
            end

        end

        function val = get.RE_fit_results(obj)
            
            if size(obj.fluxes_sub,1)<obj.minimal_length
                val = [];
            else

                if isempty(obj.RE_fit_results_)
                    obj.findOutliers;
                end

                val = obj.RE_fit_results_;

            end
            
        end
        
        function val = get.bad_star_indices(obj)

            if size(obj.fluxes_sub,1)<obj.minimal_length
                val = [];
            else

                if isempty(obj.bad_star_indices_)
                    obj.findOutliers;
                end

                val = obj.bad_star_indices_;

            end
            
        end
        
        function val = getFluxTypeIndex(obj, type)
            
            import util.text.cs;
            
            if cs(type, 'raw')
                val = 1;
            elseif cs(type, 'subtracted')
                val = 2;
            elseif cs(type, 'removed outliers')
                val = 3;
            elseif cs(type, 'calibrated')
                val = 4;
            else
                error('Unknown flux type "%s"... Use "raw" or "sub" or "rem" or "cal".', type); 
            end
            
        end
        
        function val = getFluxType(obj, index)
            
            if index==1
                val = obj.fluxes;
            elseif index==2
                val = obj.fluxes_sub;
            elseif index==3
                val = obj.fluxes_rem;
            elseif index==4
                val = obj.fluxes_cal;
            else
                error('Unknown flux index %d, use a number in the range 1-4.', index); 
            end
            
        end
        
        function val = getAirmass(obj)
            
            if isempty(obj.juldates)
                obj.calcJuldates; 
            end
            
            if isempty(obj.juldates) || isempty(obj.head) || isempty(obj.head.RA) || isempty(obj.head.DEC)
                val = [];
            else
                val = celestial.coo.airmass(obj.juldates, obj.head.RA, obj.head.DEC, [obj.head.longitude, obj.head.latitude]./180.*pi);
            end
            
        end
        
    end
    
    methods % setters
        
        function set.cat(obj, val)
            
            obj.cat = val;
            
            if ~isempty(obj.cat.data)
                obj.magnitudes = obj.cat.magnitudes; 
                obj.colors = obj.cat.data.Mag_BP - obj.cat.data.Mag_RP; 
            end
            
        end
        
        function set.frame_index(obj, val)
            
            if ~isequal(val, obj.frame_index)
                
                obj.frame_index = val;
                
                obj.clearIntermidiate;
                obj.clearPSD;
                obj.clearRE; 
                
            end
            
        end
        
        function set.fluxes_full(obj, val)
            
            obj.fluxes_full = val;
            obj.clearFits;
            obj.clearFluxes;
            obj.clearPSD;
            obj.clearRE;
            
        end
        
        function set.timestamps_full(obj, val)
            
            obj.timestamps_full = val;
            obj.clearFits;
            obj.clearFluxes;
            obj.clearPSD;
            obj.clearRE;
            
        end
        
        function set.backgrounds_full(obj, val)
            
            obj.backgrounds_full = val;
            obj.clearFluxes;
            obj.clearPSD;
            obj.clearRE;
            
        end
        
        function set.areas_full(obj, val)
            
            obj.areas_full = val;
            obj.clearFluxes;
            obj.clearPSD;
            obj.clearRE;
            
        end
        
        function set.use_subtract_backgrounds(obj, val)
            
            if ~isequal(obj.use_subtract_backgrounds, val)
                obj.use_subtract_backgrounds = val;
                obj.clearFluxes;
                obj.clearPSD;
                obj.clearRE;
            end
            
        end
        
        function set.use_background_median(obj, val)
            
            if ~isequal(obj.use_background_median, val)
                obj.use_background_median = val;
                obj.clearFluxes;
                obj.clearPSD;
                obj.clearRE;
            end
            
        end
        
        function set.use_background_fit(obj, val)
            
            if ~isequal(obj.use_background_fit, val)
                obj.use_background_fit = val;
                obj.clearFluxes;
                obj.clearPSD;
                obj.clearRE;
            end
            
        end
        
        function set.use_skip_flagged(obj, val)
            
            if ~isequal(obj.use_skip_flagged, val)
                obj.use_skip_flagged = val;
                obj.fluxes_rem_ = [];
                obj.fluxes_cal_ = [];
                obj.clearREfit;
                obj.clearPSD;
                obj.clearRE;
            end
            
        end
        
        function set.use_remove_outliers(obj, val)
            
            if ~isequal(obj.use_remove_outliers, val)
                obj.use_remove_outliers = val;
                obj.fluxes_rem_ = [];
                obj.fluxes_cal_ = [];
                obj.clearREfit;
                obj.clearPSD;
                obj.clearRE;
            end
            
        end
        
        function set.outlier_sigma(obj, val)
            
            if ~isequal(obj.outlier_sigma, val)
                obj.outlier_sigma = val;
                obj.fluxes_rem_ = [];
                obj.fluxes_cal_ = [];
                obj.clearREfit;
                obj.clearPSD;
                obj.clearRE;
            end
            
        end
        
        function set.use_skip_bad_times(obj, val)
            
            if ~isequal(obj.use_skip_bad_times, val)
                obj.use_skip_bad_times = val;
                obj.fluxes_rem_ = [];
                obj.fluxes_cal_ = [];
                obj.clearREfit;
                obj.clearPSD;
                obj.clearRE;
            end
            
        end
        
        function set.bad_times_fraction(obj, val)
            
            if ~isequal(obj.bad_times_fraction, val)
                obj.bad_times_fraction = val;
                obj.fluxes_rem_ = [];
                obj.fluxes_cal_ = [];
                obj.clearREfit;
                obj.clearPSD;
                obj.clearRE;
            end
            
        end
        
        function set.use_welch(obj, val)
            
            if ~isequal(obj.use_welch, val)
                
                obj.use_welch = val;
                obj.clearPSD;
                
            end
            
        end
        
        function set.use_airmass_correction(obj, val)
            
            if ~isequal(obj.use_airmass_correction, val)
                obj.use_airmass_correction = val;
                obj.fluxes_cal_ = [];
                obj.clearPSD;
                obj.clearRE;
            end
            
        end
        
        function set.use_color_correction(obj, val)
            
            if ~isequal(obj.use_color_correction, val)
                obj.use_color_correction = val;
                obj.fluxes_cal_ = [];
                obj.clearPSD;
                obj.clearRE;
            end
            
        end
        
        function set.use_psf_correction(obj, val)
            
            if ~isequal(obj.use_psf_correction, val)
                obj.use_psf_correction = val;
                obj.fluxes_cal_ = [];
                obj.clearPSD;
                obj.clearRE;
            end
            
        end
        
        function set.use_offset_fit(obj, val)
            
            if ~isequal(obj.use_offset_fit, val)
                obj.use_offset_fit = val;
                obj.fluxes_cal_ = [];
                obj.clearPSD;
                obj.clearRE;
            end
            
        end
        
        function set.use_zero_point(obj, val)
            
            if ~isequal(obj.use_zero_point, val)
                obj.use_zero_point = val;
                obj.fluxes_cal_ = [];
                obj.clearPSD;
                obj.clearRE;
            end
            
        end
        
        function set.zero_point_spatial_order(obj, val)
            
            if ~isequal(obj.zero_point_spatial_order, val)
                obj.zero_point_spatial_order = val;
                obj.fluxes_cal_ = [];
                obj.clearPSD;
                obj.clearRE;
            end
            
        end
        
        function set.use_sysrem(obj, val)
            
            if ~isequal(obj.use_sysrem, val)
                obj.use_sysrem = val;
                obj.fluxes_cal_ = [];
                obj.clearPSD;
                obj.clearRE;
            end
            
        end
        
        function set.sysrem_iterations(obj, val)
            
            if ~isequal(obj.sysrem_iterations, val)
                obj.sysrem_iterations = val;
                obj.fluxes_cal_ = [];
                obj.clearPSD;
                obj.clearRE;
            end
            
        end
        
        function set.use_self_exclude(obj, val)
            
            if ~isequal(obj.use_self_exclude, val)
                obj.use_self_exclude = val;
                obj.fluxes_cal_ = [];
                obj.clearPSD;
                obj.clearRE;
            end
            
        end
        
%         function set.use_savitzky_golay(obj, val)
%             
%             if ~isequal(obj.use_savitzky_golay, val)
%                 obj.use_savitzky_golay = val;
%                 if obj.use_zero_point
%                     obj.fluxes_cal_ = [];
%                     obj.clearPSD;
%                     obj.clearRE;
%                 end
%             end
%             
%         end
%         
%         function set.sg_order(obj, val)
%             
%             if ~isequal(obj.sg_order, val)
%                 obj.sg_order = val;
%                 if obj.use_zero_point
%                     obj.fluxes_cal_ = [];
%                     obj.clearPSD;
%                 end
%             end
%             
%         end
%         
%         function set.sg_length(obj, val)
%             
%             if ~isequal(obj.sg_length, val)
%                 obj.sg_length = val;
%                 if obj.use_zero_point
%                     obj.fluxes_cal_ = [];
%                     obj.clearPSD;
%                     obj.clearRE;
%                 end
%             end
%             
%         end
%         
        function set.index_flux(obj, val)
            
            if ~isequal(obj.index_flux, val)
                obj.index_flux = val;
                obj.clearFluxes;
                obj.clearPSD;
                obj.clearRE;
            end
            
        end
        
        function set.sampling_jump(obj, val)
            
            if ~isequal(obj.sampling_jump, val)
               
                obj.sampling_jump = val;
                
                obj.clearRE;
                
            end
            
        end
        
        function set.num_points_rms(obj, val)
            
            if ~isequal(obj.num_points_rms, val)
               
                obj.num_points_rms = val;
                
                obj.clearRE;
                
            end
            
        end
        
        function set.RE_fit_results(obj, val)
            
            obj.RE_fit_results_ = val;
            
        end
        
        function set.bad_star_indices(obj, val)
            
            obj.bad_star_indices_ = val;
            
        end
        
    end
    
    methods % startup/finishup/get data functions
        
        function startup(obj, num_frames, num_stars, num_apertures) % preallocate the required arrays
            
            if nargin<2 || isempty(num_frames)
                num_frames = [];
            end
            
            if nargin<3 || isempty(num_stars)
                num_stars = [];
            end
            
            if nargin<4 || isempty(num_apertures)
                num_apertures = 1; 
            end
            
            if obj.use_preallocate && obj.frame_index==1 && ~isempty(num_frames) && ~isempty(num_stars) && ~isempty(num_apertures)
                obj.timestamps_full = NaN(num_frames, 1, 'double'); 
                obj.juldates_full = NaN(num_frames, 1, 'double'); 
                obj.fluxes_full = NaN(num_frames, num_stars, num_apertures, 'single'); 
                obj.errors_full = NaN(num_frames, num_stars, num_apertures, 'single');
                obj.areas_full = NaN(num_frames, num_stars, num_apertures, 'single');
                obj.backgrounds_full = NaN(num_frames, num_stars, num_apertures, 'single');
                obj.variances_full = NaN(num_frames, num_stars, num_apertures, 'single');
                obj.centroids_x_full = NaN(num_frames, num_stars, num_apertures, 'single');
                obj.centroids_y_full = NaN(num_frames, num_stars, num_apertures, 'single');
                obj.offsets_x_full = NaN(num_frames, num_stars, num_apertures, 'single');
                obj.offsets_y_full = NaN(num_frames, num_stars, num_apertures, 'single');
                obj.widths_full = NaN(num_frames, num_stars, num_apertures, 'single');
                obj.bad_pixels_full = NaN(num_frames, num_stars, num_apertures, 'single');
                obj.flags_full = NaN(num_frames, num_stars, num_apertures, 'single');
            end
            
        end
        
        function finishup(obj) % truncate trailing (unused) NaN values in all arrays
            
            obj.timestamps_full = obj.timestamps;
            obj.juldates_full = obj.juldates;
            
            obj.fluxes_full = obj.fluxes; 
            obj.errors_full = obj.errors;
            obj.areas_full = obj.areas;
            obj.backgrounds_full = obj.backgrounds;
            obj.variances_full = obj.variances;
            obj.centroids_x_full = obj.centroids_x;
            obj.centroids_y_full = obj.centroids_y;
            obj.offsets_x_full = obj.offsets_x;
            obj.offsets_y_full = obj.offsets_y;
            obj.widths_full = obj.widths;
            obj.bad_pixels_full = obj.bad_pixels;
            obj.flags_full = obj.flags; 
            
            if ~isempty(obj.cat) && ~isempty(obj.cat.magnitudes) && size(obj.cat.magnitudes,1)>=size(obj.fluxes,2)
                obj.magnitudes = obj.cat.magnitudes;
            elseif ~isempty(obj.cat) && ~isempty(obj.cat.data) && any(contains(obj.cat.data.Properties.VariableNames, 'Mag_BP')) && size(obj.cat.data.Mag_BP,1)==size(obj.fluxes,2)
                obj.magnitudes = obj.cat.data.Mag_BP;
            end
            
        end
        
        function input(obj, varargin) % give the photometric results explicitely, they are added to the storage
            
            import util.vec.insert_matrix;
            
            if isscalar(varargin) && isa(varargin{1}, 'util.text.InputVars')
                input = varargin{1};
            else
                input = util.text.InputVars;
                input.use_ordered_numeric = 1;
                input.input_var('timestamps', [], 'times');
                input.input_var('juldates', [], 'juliandates'); 
                input.input_var('fluxes', [], 'fluxes_raw');
                input.input_var('errors', []);
                input.input_var('areas', []);
                input.input_var('backgrounds', []);
                input.input_var('variances', []);
                input.input_var('centroids_x', []);
                input.input_var('centroids_y', []);
                input.input_var('offsets_x', []);
                input.input_var('offsets_y', []);
                input.input_var('widths', []);
                input.input_var('bad_pixels', []);
                input.input_var('flags', []); 
                input.input_var('pars_struct', [], 'phot_pars');
                input.scan_vars(varargin{:});
            end
            
            N = size(input.fluxes,1);
            
            use_copy = 0; % the last input is zero telling insert_matrix NOT to make a copy of the big matrices and modify them in place
            
            obj.timestamps_full = insert_matrix(obj.timestamps_full, input.timestamps, [obj.frame_index 1], NaN, obj.use_double_up, use_copy);
            obj.juldates_full = insert_matrix(obj.juldates_full, input.juldates, [obj.frame_index, 1], NaN, obj.use_double_up, use_copy); 

            obj.fluxes_full = insert_matrix(obj.fluxes_full, input.fluxes, [obj.frame_index,1,1], NaN, obj.use_double_up, use_copy);
            obj.errors_full = insert_matrix(obj.errors_full, input.errors, [obj.frame_index,1,1], NaN, obj.use_double_up, use_copy);
            obj.areas_full = insert_matrix(obj.areas_full, input.areas, [obj.frame_index,1,1], NaN, obj.use_double_up, use_copy);
            obj.backgrounds_full = insert_matrix(obj.backgrounds_full, input.backgrounds, [obj.frame_index,1,1], NaN, obj.use_double_up, use_copy);
            obj.variances_full = insert_matrix(obj.variances_full, input.variances, [obj.frame_index,1,1], NaN, obj.use_double_up, use_copy);
            obj.centroids_x_full = insert_matrix(obj.centroids_x_full, input.centroids_x, [obj.frame_index,1,1], NaN, obj.use_double_up, use_copy);
            obj.centroids_y_full = insert_matrix(obj.centroids_y_full, input.centroids_y, [obj.frame_index,1,1], NaN, obj.use_double_up, use_copy);
            obj.offsets_x_full = insert_matrix(obj.offsets_x_full, input.offsets_x, [obj.frame_index,1,1], NaN, obj.use_double_up, use_copy);
            obj.offsets_y_full = insert_matrix(obj.offsets_y_full, input.offsets_y, [obj.frame_index,1,1], NaN, obj.use_double_up, use_copy);
            obj.widths_full = insert_matrix(obj.widths_full, input.widths, [obj.frame_index,1,1], NaN, obj.use_double_up, use_copy);
            obj.bad_pixels_full = insert_matrix(obj.bad_pixels_full, input.bad_pixels, [obj.frame_index,1,1], NaN, obj.use_double_up, use_copy);
            obj.flags_full = insert_matrix(obj.flags_full, input.flags, [obj.frame_index,1,1], NaN, obj.use_double_up, use_copy);
            
            obj.frame_index = obj.frame_index + N;
            
            if ~isempty(obj.cat) && ~isempty(obj.cat.success)
                obj.magnitudes = obj.cat.magnitudes;
            end
            
            if ~isempty(input.pars_struct)
                obj.phot_pars = input.pars_struct;
            end
            
        end
        
        function getData(obj, photometry, type) % copy the data from an img.Photometry object
            
            if nargin<3 || isempty(type)
                type = '';
            end
            
            list = {'fluxes', 'errors', 'areas', 'backgrounds', 'variances', 'offsets_x', 'offsets_y', 'centroids_x', 'centroids_y', 'widths', 'bad_pixels', 'flags'};
            
            if ~isempty(type)
                list2 = strcat(list, ['_' type]);
            else
                list2 = list;
            end
            
            input = util.text.InputVars;
            
            for ii = 1:length(list)
                input.input_var(list{ii}, photometry.(list2{ii}));
            end
            
            input.input_var('timestamps', photometry.timestamps);
            input.input_var('juldates', photometry.juldates);
            input.input_var('pars_struct', photometry.pars_struct);
            
            obj.input(input);
            
        end
        
        function getAperturesAndForced(obj, photometry) % get the data from an img.Photometry object, only focusing on aperture and forced photometry (this is already done in getData)
            
            input = util.text.InputVars;           
            
            list = {'fluxes', 'errors', 'areas', 'backgrounds', 'variances', 'offsets_x', 'offsets_y', 'centroids_x', 'centroids_y', 'widths', 'bad_pixels', 'flags'};
            list2 = strcat(list, '_ap');
            list3 = strcat(list, '_forced');
            
            for ii = 1:length(list)
                input.input_var(list{ii}, cat(3, photometry.(list2{ii}), photometry.(list3{ii})));
            end
            
            input.input_var('timestamps', photometry.timestamps);
            input.input_var('pars_struct', photometry.pars_struct);
            
            obj.input(input);
            
        end
        
        function keepStarIndices(obj, indices)
            
            obj.fluxes_full = obj.fluxes_full(:,indices,:);
            obj.errors_full = obj.errors_full(:,indices,:);
            obj.areas_full = obj.areas_full(:,indices,:);
            obj.backgrounds_full = obj.backgrounds_full(:,indices,:);
            obj.variances_full = obj.variances_full(:,indices,:);
            obj.centroids_x_full = obj.centroids_x_full(:,indices,:);
            obj.centroids_y_full = obj.centroids_y_full(:,indices,:);
            obj.offsets_x_full = obj.offsets_x_full(:,indices,:);
            obj.offsets_y_full = obj.offsets_y_full(:,indices,:);
            obj.widths_full = obj.widths_full(:,indices,:);
            obj.bad_pixels_full = obj.bad_pixels_full(:,indices,:);
            obj.flags_full = obj.flags_full(:,indices,:);
            
            if ~isempty(obj.magnitudes)
                obj.magnitudes = obj.magnitudes(indices); 
            end
            
            if ~isempty(obj.cat)
                
                if ~isempty(obj.cat.data)
                    obj.cat.data = obj.cat.data(indices,:); 
                end
                
                if ~isempty(obj.cat.magnitudes)
                    obj.cat.magnitudes = obj.cat.magnitudes(indices);
                end
                
                if ~isempty(obj.cat.positions)
                    obj.cat.positions = obj.cat.positions(indices, :);
                end
                
                if ~isempty(obj.cat.coordinates)
                    obj.cat.coordinates = obj.cat.coordinates(indices, :);
                end
                
                if ~isempty(obj.cat.temperatures)
                    obj.cat.temperatures = obj.cat.temperatures(indices); 
                end
                
            end
            
            obj.clearIntermidiate;
            obj.clearFits;
            obj.clearFluxes;
            obj.clearPSD;
            obj.clearRE;
            obj.clearREfit;
            
        end
            
        function keepFrameIndices(obj, indices)
            
            obj.timestamps_full = obj.timestamps_full(indices);
            obj.juldates_full = obj.juldates(indices);
            
            obj.fluxes_full = obj.fluxes_full(indices,:,:);
            obj.errors_full = obj.errors_full(indices,:,:);
            obj.areas_full = obj.areas_full(indices,:,:);
            obj.backgrounds_full = obj.backgrounds_full(indices,:,:);
            obj.variances_full = obj.variances_full(indices,:,:);
            obj.centroids_x_full = obj.centroids_x_full(indices,:,:);
            obj.centroids_y_full = obj.centroids_y_full(indices,:,:);
            obj.offsets_x_full = obj.offsets_x_full(indices,:,:);
            obj.offsets_y_full = obj.offsets_y_full(indices,:,:);
            obj.widths_full = obj.widths_full(indices,:,:);
            obj.bad_pixels_full = obj.bad_pixels_full(indices,:,:);
            obj.flags_full = obj.flags_full(indices,:,:);
            
            obj.frame_index = size(obj.fluxes_full,1)+1; 
            
            obj.clearIntermidiate;
            obj.clearFits;
            obj.clearFluxes;
            obj.clearPSD;
            obj.clearRE;
            obj.clearREfit;
            
        end
        
        % add a method to remove 3rd dimension?
        
    end
    
    methods % calculations
        
        function [idx, fr] = findOutliers(obj, mean_flux, mean_error, iterations, sigma)
            
            if nargin<2 || isempty(mean_flux)
                mean_flux = nanmean(obj.fluxes_rem);
            end
            
            if nargin<3 || isempty(mean_error)
                mean_error = nanmedian(util.series.binning(obj.fluxes_rem, 100, 'func', 'std'),1);
            end
            
            if nargin<4 || isempty(iterations)
                iterations = 2;
            end
            
            if nargin<5 || isempty(sigma)
                sigma = 3;
            end
            
            F = mean_flux;
            E = mean_error;
            
            idx = F<=0 | isnan(F) | isinf(F) | E<=0 | isnan(E) | isinf(E);
            
            % find outliers that have too many measurements equal zero!
            frac_zero = nansum(obj.fluxes_rem==0, 1)./size(obj.fluxes_rem, 1); 
            idx = idx | frac_zero>0.1;
            
            % find outliers with many bad pixels a serious fraction of the time
            frac_bad_px = nansum(obj.bad_pixels(:,:,obj.index_flux_number)>20, 1)./size(obj.bad_pixels(:,:,obj.index_flux_number), 1);
            idx = idx | frac_bad_px>0.1;
            
            X = F;
            Y = E./F;
            err = 1./sqrt(F); 
            
            X(idx) = NaN;
            Y(idx) = NaN;
            err(idx) = NaN;
            
            for ii = 1:iterations

                fr = util.fit.power_law(X, Y, err, 'slopes', [-1, 0], 'breaks', 1000);

                relative_residuals = (fr.y-fr.func(fr.x))./fr.func(fr.x); 

                rms_rr = nanstd(relative_residuals);

                norm_residuals = abs(relative_residuals./rms_rr); 

                idx = idx | norm_residuals>sigma; % outliers are more than 3 sigma from the curve

            end
            
            obj.RE_fit_results = fr;
            obj.bad_star_indices = idx;
            
        end
        
        function flux = calculateCalibration(obj, flux)
            
            if obj.use_airmass_correction
                
                obj.airmass_vector = obj.getAirmass; 
                
                [L, b] = util.units.flux2lup(flux); % this also finds a reasonable estimate for the softening parameter b
                
%                 if obj.use_color_correction
                    
%                 else
                    fr = util.fit.polyfit(obj.airmass_vector, L, 'order', 2, 'double', 1, 'sigma', 3); 
%                 end
                
                model = [fr.ym]; 
                
                L_det = L - model + nanmean(L, 1); 
                
                flux = util.units.lup2flux(L_det, b); 
                
            end
            
            if obj.use_psf_correction
                
                F = nanmean(flux, 1); % weight for each star
                
                W = util.vec.weighted_average(obj.widths(:,:,obj.index_flux_number), F, 2); 
                
                % apply the same correction we applied to the fluxes before
                if obj.use_airmass_correction 
                    
                    fr_w = util.fit.polyfit(obj.airmass_vector, W, 'order', 2, 'double', 1, 'sigma', 3); 
                    
                    W = W - fr_w.ym + nanmedian(W); 
                
                end
                
                obj.width_vector = W;
                
                N = 20; 
                window_size = size(flux,1)/N; 
                window_size = max(window_size, 64); 
                
                flux = util.series.remove_aux(flux, obj.width_vector, 'window', window_size); 
                
%                 [L, b] = util.units.flux2lup(flux); % this also finds a reasonable estimate for the softening parameter b
%                 
%                 fr = util.fit.polyfit(W, L, 'order', 2, 'double', 1, 'sigma', 3); 
%                 
%                 model = [fr.ym]; 
%                 
%                 L_det = L - model + nanmean(L, 1); 
%                 
%                 flux = util.units.lup2flux(L_det, b);
                
            end
            
            if obj.use_zero_point
                
                if obj.zero_point_spatial_order==0

                    w = nanmean(flux); % weights are given by the average flux of each star
                    w = sqrt(abs(w)); 

                    f_average = nanmean(flux.*w,2); % weighted average of each frame
                    f_average_norm = f_average./nanmean(f_average); % normalize by the average of averages

                    flux = flux./f_average_norm;

                    % need to keep a record of the ZP for each frame here as well
                    
                else
                    
                    if isempty(obj.cat) || isempty(obj.cat.success) || obj.cat.success==0
                        error('Must have a catalog with magnitudes to calculate the ZP fit!'); 
                    end
                    
                    M = obj.cat.magnitudes; 
                    f = flux; 
                    M = vertcat(M, nan(size(f,2) - size(M,1),1)); % add NaN magnitudes to forced targets
                    f(f<=0) = NaN; % get rid of all negative flux values (complex magnitudes are no good)
                    zp = M + 2.5*log10(f)'; 
                    
                    x = obj.centroids_x(:,:,obj.index_flux_number)';                    
                    y = obj.centroids_y(:,:,obj.index_flux_number)';
                
                    if ~isempty(obj.head)
                        xc = floor(obj.head.NAXIS2/2)+1;
                        yc = floor(obj.head.NAXIS1/2)+1;
                    else
                        xc = nanmean(x(:));
                        yc = nanmean(y(:));
                    end
                    
                    e = ones(size(flux,1),1).*(nanstd(flux)./nanmean(flux)); 
                    
                    obj.zp_fit_results = util.fit.surf_poly(x-xc,y-xc,zp,'errors', e', 'order', obj.zero_point_spatial_order, ...
                        'double', 1, 'warning', 0); 
                    
                    ZP = [obj.zp_fit_results.vm]'; % the zero point we found for each star and each frame
                    C = [obj.zp_fit_results.coeffs]'; 
                    
                    obj.mean_zp = C(:,1); % the average zero point per image
                    
%                     M_inst = (nanmean(C(:,1),1) - 2.5.*log10(abs(nanmean(flux,1))));
                    
                    flux = flux.*10.^(0.4.*(nanmean(ZP(:))-ZP)); 
                    
                end
                
            end
            
            if obj.use_offset_fit
                
                x = obj.offsets_x(:,:,obj.index_flux_number); 
                y = obj.offsets_y(:,:,obj.index_flux_number); 
                
                fr = util.fit.surf_poly(x, y, flux, 'order', 2, 'warning', false, 'iter', 3, 'sigma', 3, 'double', 1); 
                
                c = [fr.coeffs]; % coefficients: the first is a constant term that we don't want to subtract
                flux = flux - [fr.vm] + c(1,:); 
                
            end
            
            if obj.use_sysrem
                                                
                [flux, a_vec, c_vec] = obj.calculateSysrem(flux); 

                obj.sysrem_a_vec = a_vec;
                obj.sysrem_c_vec = c_vec;
            
            end
            
        end
        
        function [flux, a_vec, c_vec] = calculateSysrem(obj, flux)
            
            idx = obj.bad_star_indices; 
            % maybe add a search for bad batches?

%            flux(:,idx) = NaN;

            a_vec = [];
            c_vec = [];

            for ii = 1:obj.sysrem_iterations

                e = nanmedian(util.series.binning(flux, 100, 'func', 'std'),1);
                e(e==0 | isnan(e)) = Inf; 
%                             e = sqrt(1./nanmean(f));
                [flux,a,c] = util.series.sysrem(flux, e, 'iterations', 3, 'stars', idx, 'self', obj.use_self_exclude, 'fft', 0); 

                a_vec(:,ii) = nanmean(a,2);
                c_vec(ii,:) = c;

            end
            
        end
        
        function Im = getZeroPointImages(obj, im_size)
            
            if nargin<2 || isempty(im_size)
                im_size = 400;
            end
            
            im_size = util.vec.imsize(im_size); 
            
            x = linspace(-obj.head.NAXIS2/2, obj.head.NAXIS2, im_size(2)); 
            y = linspace(-obj.head.NAXIS1/2, obj.head.NAXIS1, im_size(1)); 
            [X,Y] = meshgrid(x,y); 
            
            C = [obj.zp_fit_results.coeffs]; 
            px = obj.zp_fit_results(1).coeff_powers_x;
            py = obj.zp_fit_results(1).coeff_powers_y;
            
            N = size(C,2); 
            
            Im = zeros([im_size, N]); 
            
            for ii = 1:N
                
                new_image = zeros(im_size); 
                
                for jj = 1:size(C,1)
                     new_image = new_image + C(jj,ii).*X.^px(jj).*Y.^py(jj);
                end
                
                Im(:,:,ii) = new_image;
                
            end
            
            
        end
        
        function [f, f_transits] = injectionTest(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('number', 10); 
            input.input_var('depth', 0.05); 
            input.input_var('width', 100); 
            input.input_var('surround', 5); 
            input.input_var('stars', 1, 'star_index', 'star_indices');
            input.input_var('frames', [], 'frame_index', 'frame_indices'); 
            input.input_var('font_size', 26); 
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('show_indices', 1:10); 
            input.scan_vars(varargin{:});
            
            if isempty(input.stars)
                input.stars = 1:size(obj.fluxes,2); 
            end
            
            if isempty(input.frames)
                input.frames = 1:size(obj.fluxes,1); 
            end
                        
            if ~isequal(input.frames(1):input.frames(end), input.frames)
                error('The input "frames" must be a continuous range of frame indices with step=1'); 
            end
            
            if input.frames(end)>size(obj.fluxes,1)
                error('Last frame index, %d, is out of bounds for the number of available frames, %d.', input.frames(end), size(obj.fluxes,1)); 
            end
            
            if input.surround<3
                error('Input "surround" must be at least 3, to allow space around the transit to get the linear fit'); 
            end
            
            if length(input.frames)<input.width*input.surround*input.number
                error('frame indices cover only %d data points, which is too short to use with width %d, surround %d and number %d (total length needed is %d). ', ...
                    length(input.frames), input.width, input.surround, input.number, input.width.*input.surround.*input.number);
            end
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            time_indices = 1:input.width*input.surround:size(obj.fluxes,1); % all possible intervals that can fit a transit
            
            frames_picked = time_indices(randperm(length(time_indices)-1, input.number)); % transit frames are chosen without repitition
            
            stars_picked = input.stars(randi(length(input.stars), [1, input.number])); % star frames can be chosen with repitition
            
            f = obj.fluxes_rem; % get a copy of the data before calibration
            
            for ii = 1:input.number
                
                idx_region = ( frames_picked(ii):(frames_picked(ii)+ceil(input.width.*input.surround)-1) )';
                f_region = f(idx_region, stars_picked(ii)); 
                
                fr = util.fit.polyfit(idx_region, f_region, 'order', 1);
                                
                start_frame = frames_picked(ii)+randi(ceil(input.width.*(input.surround-1)));
                
                idx_transit(:,ii) = (start_frame:(start_frame+ceil(input.width)-1)); 
                
                missing_light(:,ii) = input.depth.*fr.func(idx_transit(:,ii)); 
                
                f(idx_transit(:,ii), stars_picked(ii)) = f(idx_transit(:,ii), stars_picked(ii))-missing_light(:,ii);
                
            end
            
            t = tic; 
            
            f = obj.calculateCalibration(f); 
            
            fprintf('Time to recalibrate fluxes was %f seconds...\n', toc(t)); 
            
            f_transits = ones(size(f,1),1).*nanmedian(f,1); 
            
            for ii = 1:input.number
                
                f_transits(idx_transit(:,ii),stars_picked(ii)) = f_transits(idx_transit(:,ii),stars_picked(ii)).*(1-input.depth); 
                
            end
            
            hold(input.ax, 'off');
            plot(input.ax, 1:size(f,1), f(:,input.show_indices), '.');
            
            input.ax.ColorOrderIndex = 1;
            hold(input.ax, 'on'); 
            plot(input.ax, 1:size(f,1), f_transits(:, input.show_indices), '-');
            hold(input.ax, 'off'); 
            
            if nargout==0
                clear f;
                clear f_transits;
            end
            
        end
        
        function [pxx, freq] = calculatePSD(obj)
            
            f = obj.fluxes_cal;
            f(isinf(f)) = NaN;
            f = fillmissing(f, 'linear'); 
            
            if obj.use_welch
            
%                 [pxx, freq] = pwelch(f, [], [], size(f,1), 1./median(diff(obj.timestamps)), 'power'); 
                [pxx, freq] = pwelch(f, [], [], size(f,1), 1./median(diff(obj.timestamps))); 
                
            else
                pxx = abs(fft(f-nanmean(f,1)));
                pxx = pxx(1:floor(size(f,1)/2)+1,:); 
                T = (obj.timestamps(end)-obj.timestamps(1)); 
                freq = (0:1/T:(size(pxx,1)-1)/T)';
            end
            
            pxx = real(pxx); 
                        
        end
        
        function calculateRelativeError(obj)
            
            import util.series.binning;
            
            obj.clearRE;

            for jj = 1:4 % go over flux types

%                 f = obj.fluxes_cal;
                f = obj.getFluxType(jj); 
                if size(f,3)>1
                    f = f(:,:,obj.index_flux_number); 
                end
                t = obj.timestamps;

                if isempty(f) || isempty(t)
                    return;
                end

                F = nanmedian(f); % the average flux is used to divide the rms and get the relative error (RE)

                for ii = 1:1e3 % arbitrary loop length with a break condition

                    obj.bin_widths_seconds_(ii) = nanmedian(diff(t),1); % get the average time sampling

                    obj.local_RE_(ii,:,jj) = nanmedian(binning(f, obj.num_points_rms, 'function', 'std'),1)./F; % the local relative error

                    obj.total_RE_(ii,:,jj) = nanstd(f)./abs(F); % the relative error on the entire dataset

                    f = binning(f, obj.sampling_jump); % downsample the lightcurves by a set factor 
                    t = binning(t, obj.sampling_jump); % make sure the timesteps are binned like the fluxes

                    if size(f,1)<obj.num_points_rms
                        break; % if the newly binned flux is too short to calculate local rms, there's no point to keep going! 
                    end

                end

            end
            
        end
        
        function rho = getCorrelationCoeff(obj, type) % calculate Pearson's coefficient with respect to flux and any other quantity
            
            f1 = obj.fluxes;
            f1(obj.flags>0) = NaN;
            f2 = obj.(type); 
            f2(obj.flags>0) = NaN;
            
            M1 = nanmean(f1);
            M2 = nanmean(f2);
            s1 = nanstd(f1);
            s2 = nanstd(f2); 
            C = nanmean((f1-M1).*(f2-M2)); 
            
            rho = C./s1./s2;
            
        end
        
        function flux_corr = calibrateFlux(obj, type) % not working yet
            
            rho = obj.getCorrelationCoeff(type); 
            
            Av = nanmean(obj.(type) .* ~obj.flags, 2); % maybe replace this with the flux weighted sum?
            
            flux_corr = obj.fluxes.*rho./(Av./nanmean(Av));
            
        end
        
        function calcJuldates(obj)
            
            if isempty(obj.head) || isempty(obj.head.ENDTIME) || isempty(obj.head.END_STAMP)
                error('Must have header with ENDTIME and END_STAMP to calculate julian dates'); 
            end
            
            t0 = util.text.str2time(obj.head.ENDTIME) - seconds(obj.head.END_STAMP); % datetime at timestamps=0
            
            t = t0 + seconds(obj.timestamps_full); 
            
            obj.juldates_full = juliandate(t); 
            
        end
        
        function repairTimestamps(obj)
            
            t = obj.timestamps_full;
            
            dt = prctile(diff(t), 70); % find the common time step (should be 0.04 seconds)
            
            t0 = t(1); 
            
            t = t0 + dt.*(0:length(t)-1)'; % adjust the times to be interpolated from the point we found
            
            if isempty(obj.timestamps_full_original)
                obj.timestamps_full_original = obj.timestamps_full;
            end
            
            obj.timestamps_full = t;
            
            if isempty(obj.juldates_full_original)
                obj.juldates_full_original = obj.juldates_full;
            end
            
            obj.juldates_full = obj.juldates_full(1) + (obj.timestamps_full - obj.timestamps_full(1))./24./3600;
            
            obj.timestamps_ = [];
            obj.juldates_ = [];
                        
        end
        
    end
    
    methods % load/save
        
        function saveDialog(obj, filename)
            
            if nargin<2 || isempty(filename)
                
                filename = 'Lightcurves'; 

                if isempty(obj.head) || ~ischar(obj.head.STARTTIME) || length(obj.head.STARTTIME)<10
                    run_date = ''; 
                else
                    run_date = obj.head.STARTTIME(1:10);
                end
                
                if isempty(obj.head) || isempty(obj.head.OBJECT)
                    run_name = ''; 
                else
                    run_name = obj.head.OBJECT;
                end
                
                if ~isempty(run_date)
                    filename = [filename '_' run_date];
                end
                
                if ~isempty(run_name)
                    filename = [filename '_' run_name];
                end

            end
            
            [filepath,filename,ext] = fileparts(filename);
            
            if ~isempty(ext)
                filename = [filename ext];
            else
                filename = [filename '.mat']; 
            end
            
            if isempty(filepath)
                filepath = fullfile(getenv('DATA'), 'WFAST/saved/lightcurves/'); 
            end
            
            [filename, filepath] = uiputfile(fullfile(filepath, filename)); 
            
            if ~isequal(filename, 0)
                obj.saveAsMAT(fullfile(filepath, filename)); 
            end
            
        end
        
        function saveAsMAT(obj, filename, star_indices, frame_indices, aperture_indices) % save only the data, in raw matrices, to a MAT file, along with phot_pars and a short readme. 
            
            if nargin<2 || isempty(filename)
                disp('Must supply a filename! Usage: obj.saveAsMAT(filename, star_indices, frame_indices, aperture_indices)'); 
            end
            
            if nargin<3 || isempty(star_indices) || (ischar(star_indices) && strcmpi(star_indices, 'all'))
                star_indices = [];
            end
            
            if nargin<4 || isempty(frame_indices) || (ischar(frame_indices) && strcmpi(frame_indices, 'all'))
                frame_indices = [];
            end
            
            if nargin<5 || isempty(aperture_indices) || (ischar(aperture_indices) && strcmpi(aperture_indices, 'all'))
                aperture_indices = [];
            end
            
            timestamps = obj.timestamps;
            
            if isempty(star_indices)
                fluxes = obj.fluxes;
                errors = obj.errors;
                areas = obj.areas;
                backgrounds = obj.backgrounds;
                variances = obj.variances;
                offsets_x = obj.offsets_x;
                offsets_y = obj.offsets_y;
                centroids_x = obj.centroids_x;
                centroids_y = obj.centroids_y;
                widths = obj.widths;
                bad_pixels = obj.bad_pixels;
                flags = obj.flags;
            else
                fluxes = obj.fluxes(:,star_indices,:);
                errors = obj.errors(:,star_indices,:);
                areas = obj.areas(:,star_indices,:);
                backgrounds = obj.backgrounds(:,star_indices,:);
                variances = obj.variances(:,star_indices,:);
                offsets_x = obj.offsets_x(:,star_indices,:);
                offsets_y = obj.offsets_y(:,star_indices,:);
                centroids_x = obj.centroids_x(:,star_indices,:);
                centroids_y = obj.centroids_y(:,star_indices,:);
                widths = obj.widths(:,star_indices,:);
                bad_pixels = obj.bad_pixels(:,star_indices,:);
                flags = obj.flags(:,star_indices,:);
            end
            
            if ~isempty(frame_indices)
                
                fluxes = fluxes(frame_indices,:,:);
                errors = errors(frame_indices,:,:);
                areas = areas(frame_indices,:,:);
                backgrounds = backgrounds(frame_indices,:,:);
                variances = variances(frame_indices,:,:);
                offsets_x = offsets_x(frame_indices,:,:);
                offsets_y = offsets_y(frame_indices,:,:);
                centroids_x = centroids_x(frame_indices,:,:);
                centroids_y = centroids_y(frame_indices,:,:);
                widths = widths(frame_indices,:,:);
                bad_pixels = bad_pixels(frame_indices,:,:);
                flags = obj.flags(frame_indices,:,:);
                
                timestamps = timestamps(frame_indices);
            
            end
            
            if ~isempty(aperture_indices)
                
                fluxes = fluxes(:,:,aperture_indices);
                errors = errors(:,:,aperture_indices);
                areas = areas(:,:,aperture_indices);
                backgrounds = backgrounds(:,:,aperture_indices);
                variances = variances(:,:,aperture_indices);
                offsets_x = offsets_x(:,:,aperture_indices);
                offsets_y = offsets_y(:,:,aperture_indices);
                centroids_x = centroids_x(:,:,aperture_indices);
                centroids_y = centroids_y(:,:,aperture_indices);
                widths = widths(:,:,aperture_indices);
                bad_pixels = bad_pixels(:,:,aperture_indices);
                flags = obj.flags(:,:,aperture_indices);
                
            end
            
            readme = 'Some info about the data stored in this file:';
            readme = sprintf('%s\n *timestamps: a column vector with %d elements, with the relative time of each frame', readme, size(timestamps, 1));
            readme = sprintf('%s\n *fluxes: the actual lightcurves, dim 1 is %d frames dim 2 is %d stars, dim 3 is for the different apertures+forced aperture',readme, size(fluxes,1), size(fluxes,2));
            readme = sprintf('%s\n *errors: calculated from the variance map and source intensity. ', readme); 
            readme = sprintf('%s\n *areas: size of aperture, removing bad pixels and so on.', readme);
            readme = sprintf('%s\n *backgrounds: measured intensity, per pixel, in the annulus.', readme);
            readme = sprintf('%s\n *variances: measured noise variance, per pixel, in the annulus.', readme);
            readme = sprintf('%s\n *offsets_x/y: distance from center of cutout (in pixels)', readme); 
            readme = sprintf('%s\n *centroids_x/y: position within image (in pixels)', readme);
            readme = sprintf('%s\n *widths: calculated using the 2nd moments, averaging the minor and major axis', readme);
            readme = sprintf('%s\n *bad_pixels: number of bad pixels in the aperture area.', readme);
            readme = sprintf('%s\n *flags: this is marked as non-zero if the photometry did not go well (e.g., large offsets, negative widths)', readme); 
            
            phot_pars = obj.phot_pars;
            
            header = util.oop.full_copy(obj.head);
            header.phot_pars = phot_pars; 
            
            if ~isempty(obj.cat)
                cat = obj.cat.data;
            end
            
            save(filename, 'timestamps', 'fluxes', 'errors', 'areas', 'backgrounds', 'variances',...
                'offsets_x', 'offsets_y', 'centroids_x', 'centroids_y', 'widths', ...
                'bad_pixels', 'flags', 'readme', 'phot_pars', 'header', 'cat', '-v7.3');
            
        end
        
        function loadFromMAT(obj, filename)
            
            if nargin<2 || isempty(filename)
                disp('Must supply a filename!'); 
            end
                        
            loaded = load(filename); 
            
            obj.reset; % reset after making sure we managed to load this file! 
            
            if isfield(loaded, 'timestamps'), obj.timestamps_full = loaded.timestamps; end
            if isfield(loaded, 'fluxes'), obj.fluxes_full = loaded.fluxes; end
            if isfield(loaded, 'errors'), obj.errors_full = loaded.errors; end
            if isfield(loaded, 'areas'), obj.areas_full = loaded.areas; end
            if isfield(loaded, 'backgrounds'), obj.backgrounds_full = loaded.backgrounds; end
            if isfield(loaded, 'variances'), obj.variances_full = loaded.variances; end
            if isfield(loaded, 'offsets_x'), obj.offsets_x_full = loaded.offsets_x; end
            if isfield(loaded, 'offsets_y'), obj.offsets_y_full = loaded.offsets_y; end
            if isfield(loaded, 'centroids_x'), obj.centroids_x_full = loaded.centroids_x; end
            if isfield(loaded, 'centroids_y'), obj.centroids_y_full = loaded.centroids_y; end
            if isfield(loaded, 'widths'), obj.widths_full = loaded.widths; end
            if isfield(loaded, 'bad_pixels'), obj.bad_pixels_full = loaded.bad_pixels; end
            if isfield(loaded, 'flags'), obj.flags_full = loaded.flags; end
            if isfield(loaded, 'phot_pars'), obj.phot_pars = loaded.phot_pars; end
            
            obj.head = head.Header.empty;
            if isfield(loaded, 'pars') && isa(loaded.pars, 'head.Parameters'), obj.head = loaded.pars; end
            if isfield(loaded, 'pars') && isa(loaded.pars, 'struct'), obj.head.loadFromStruct(loaded.pars); end 
            if isfield(loaded, 'header') && isa(loaded.header, 'head.Header'), obj.head = loaded.header; end
            if isfield(loaded, 'header') && isa(loaded.header, 'struct'), obj.head.loadFromStruct(loaded.header); end
            
            obj.cat = head.Catalog.empty;
            if isfield(loaded, 'cat') && isa(loaded.cat, 'head.Catalog')
                obj.cat = loaded.cat; 
                obj.cat.head = obj.head; 
            end
            
            if isfield(loaded, 'cat') && isa(loaded.cat, 'table')
                obj.cat = head.Catalog;
                obj.cat.head = obj.head; 
                obj.cat.loadFromTable(loaded.cat); 
            end
            
            obj.frame_index = size(obj.fluxes_full,1) + 1; 
            
        end
        
        function loadHDF5(obj, folder, use_reset)
            
            if isa(folder, 'util.sys.WorkingDirectory')
                % do nothing! 
            elseif ischar(folder) && exist(folder, 'dir')
                folder = util.sys.WorkingDirectory(folder); 
            else
                error('Must input a valid folder name or a WorkingDirectory object'); 
            end
            
            input = util.text.InputVars;
            % ...
            
            if nargin<3 || isempty(use_reset)
                use_reset = 1; % by default must reset this object before adding another folder!
            end
            
            if use_reset, obj.reset; end
            
            files = folder.match('*.h5*'); 
            list = {'timestamps', 'fluxes', 'errors', 'areas', 'backgrounds', 'variances', ...
                'centroids_x', 'centroids_y', 'offsets_x', 'offsets_y', 'widths', ...
                'bad_pixels', 'flags'}; 
            
            prog = util.sys.ProgressBar;
            prog.start(length(files)); 
            
            for ii = 1:length(files)
                
                if ii==1 % preallocate on the first batch...
                    f = h5read(files{1}, '/fluxes'); 
                    obj.startup(length(files)*size(f,1), size(f,2), size(f,3)); % preallocate! 
                end
                
                for jj = 1:length(list)
                    
                    name = list{jj}; 
                    name_slash = ['/' name]; 
                    name_full = [name '_full']; 
                    
                    use_copy = 0; % the last input is zero telling insert_matrix NOT to make a copy of the big matrices and modify them in place
                    
                    if strcmp(name, 'timestamps')
                        start_vec = [obj.frame_index 1];
                    else
                        start_vec = [obj.frame_index 1 1];
                    end
                    
                    data = h5read(files{ii}, name_slash);
                    if ~isempty(data)
                        N = size(data,1);
                        obj.(name_full) = util.vec.insert_matrix(obj.(name_full), data, start_vec, NaN, obj.use_double_up, use_copy);
                    end
                    
                end
                
                obj.frame_index = obj.frame_index + N;
                
                prog.showif(ii); 
                
            end
            
            prog.finish; 
            
        end
        
        function saveFigures(obj)

            idx = cellfun(@isvalid, obj.figures_spawned); 
            obj.figures_spawned = obj.figures_spawned(idx); 
            
            d = fullfile(getenv('DATA'), '/WFAST/figures');
            
            for ii = 1:length(obj.figures_spawned)
                savefig(fullfile(d, [strrep(obj.figures_spawned{ii}.UserData, {':','.'}, '_') '.fig']));
            end
            
        end
        
    end
    
    methods (Static=true)
        
        function obj = load(filename_or_folder, varargin)
            
            folder = [];
            filename = []; 
            
            if ischar(filename_or_folder) && exist(filename_or_folder, 'file')
                filename = filename_or_folder; 
            elseif isa(filename_or_folder, 'util.sys.WorkingDirectory')
                folder = filename_or_folder;
            elseif ischar(filename_or_folder) && exist(filename_or_folder, 'dir')
                folder = util.sys.WorkingDirectory(filename_or_folder); 
            else
                error('Must input a valid file or folder name or a WorkingDirectory object'); 
            end
            
            if isempty(filename) && isempty(folder)
                error('Must input a valid file or folder name or a WorkingDirectory object'); 
            end
            
            obj = img.Lightcurves; 
            
            if isempty(filename)
                obj.loadFromMAT(filename); 
            elseif isempty(folder)
                obj.loadHDF5(folder); 
            end
            
            if nargout==0
                assignin('base', 'lightcurves',obj)
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin) % plot the data, according to the different object parameters. Varargin can include "font_size" and "ax" (defaults to the GUI axis if it exists)
            
            import util.text.cs;
            
            input = util.text.InputVars;
            input.input_var('font_size', 26); 
            input.input_var('ax', [], 'axes', 'axis');
            input.scan_vars(varargin{:});
            
            if isempty(input.ax)
                if ~isempty(obj.gui) && obj.gui.check
                    input.ax = obj.gui.axes_image;
                else
                    input.ax = gca;
                end
            end
            
            if ~isempty(obj.last_shown)
                obj.last_xlim = input.ax.XLim;
                obj.last_ylim = input.ax.YLim;
            end
            
            % this is from: http://undocumentedmatlab.com/articles/determining-axes-zoom-state
            app_info = getappdata(input.ax, 'matlab_graphics_resetplotview'); % information about zoom state (empty if unzoomed!)
            
            cla(input.ax, 'reset'); 
            
%            delete(allchild(input.ax));
%            delete(findobj(input.ax.Parent, 'type', 'legend'));
            
            if obj.use_show_datetime && ~isempty(obj.juldates)
                xlabel(input.ax, ''); 
            else
                xlabel(input.ax, 'time (seconds)');
            end
            
            if cs(obj.show_what, 'fluxes')
                
                if cs(obj.show_flux_type, 'raw')
                    obj.addPlots(input.ax, obj.fluxes(:,:,obj.index_flux_number));
                elseif cs(obj.show_flux_type, 'sub')
                    obj.addPlots(input.ax, obj.fluxes_sub);
                elseif cs(obj.show_flux_type, 'rem')
                    obj.addPlots(input.ax, obj.fluxes_rem);
                elseif cs(obj.show_flux_type, 'cal')
                    obj.addPlots(input.ax, obj.fluxes_cal);
                elseif cs(obj.show_flux_type, 'all')
                    obj.addPlots(input.ax, obj.fluxes_sub, [], '.');
                    obj.addPlots(input.ax, obj.fluxes_rem, [], '-');
                    obj.addPlots(input.ax, obj.fluxes_cal, [], '--');
                    legend(input.ax, {'flux raw', 'flux rem', 'flux cal'}, 'location', 'NorthEast');                    
                end
                
                ylabel(input.ax, 'flux (counts)');
            
            elseif cs(obj.show_what, 'areas')
                obj.addPlots(input.ax, obj.areas);
                ylabel(input.ax, 'areas (pixels in aperture)');
            elseif cs(obj.show_what, 'offsets')
                obj.addPlots(input.ax, obj.offsets_x, [], '-');
                obj.addPlots(input.ax, obj.offsets_y, [], ':');
                ylabel(input.ax, 'offset inside cutouts (pixels)');
            elseif cs(obj.show_what, 'centroids')
                obj.addPlots(input.ax, obj.centroids_x, [], '-');
                obj.addPlots(input.ax, obj.centroids_y, [], ':');
                ylabel(input.ax, 'centroid positions in image (pixels)');
            elseif cs(obj.show_what, 'widths')
                obj.addPlots(input.ax, obj.widths);
                ylabel(input.ax, 'widths using 2nd moment (pixels)');
            elseif cs(obj.show_what, 'backgrounds')
                obj.addPlots(input.ax, obj.backgrounds);
                ylabel(input.ax, 'background (counts/pixel)');
            elseif cs(obj.show_what, 'variances')
                obj.addPlots(input.ax, obj.variances);
                ylabel(input.ax, 'variance (counts^2/pixel)'); 
            elseif cs(obj.show_what, 'bad_pixels')
                obj.addPlots(input.ax, obj.bad_pixels);
                ylabel(input.ax, 'number of bad pixels in aperture'); 
            elseif cs(obj.show_what, 'power spectra')
                
                obj.addPlots(input.ax, sqrt(obj.power_spectra), obj.psd_freq);
                if obj.use_welch
                    ylabel(input.ax, 'Power spectrum^{(1/2)} (Welch)'); 
                else
                    ylabel(input.ax, 'Power spectrum^{(1/2)} (FFT)'); 
                end
                xlabel(input.ax, 'Frequency [Hz]'); 
                
            elseif cs(obj.show_what, 'stats', 'statistics')
                obj.showStats('ax', input.ax, 'font_size', input.font_size);
            elseif cs(obj.show_what, 'binning')
                obj.showBinning('ax', input.ax, 'font_size', input.font_size); 
            else
                error('Unknown option to to show_what: "%s", use "fluxes" or "offset" etc...', obj.show_what);
            end
            
            if cs(obj.show_what, 'power spectra')
                
                input.ax.YScale = 'log';
                
                if obj.use_show_log
                    input.ax.XScale = 'log';
                else
                    input.ax.XScale = 'linear';
                end
                
            elseif cs(obj.show_what, 'stats', 'statistics', 'binning')
                % pass
            else
                
                input.ax.XScale = 'linear';
                
                if obj.use_show_log
                    input.ax.YScale = 'log';
                else
                    input.ax.YScale = 'linear';
                end
                
            end
            
            if cs(obj.show_what, 'stats', 'statistics', 'binning')
                use_restore_x = 0; % for the statistics plots we don't need to keep the zoom position
            elseif isequal(obj.show_what, obj.last_shown) % same type of data, must restore x
                use_restore_x = 1;
            elseif ~cs(obj.show_what, 'power spectra', 'stats', 'statistics') && ~cs(obj.last_shown, 'power spectra', 'stats', 'statistics', 'binning') % data types like flux/background/bad_pixels are all time based
                use_restore_x = 1;
            else
                use_restore_x = 0;
            end

            if use_restore_x && ~isempty(obj.last_xlim) && isequal(class(input.ax.XLim), class(obj.last_xlim))
                
                max_xlim = input.ax.XLim; 
                input.ax.XLim = obj.last_xlim;
                
                if ~isempty(app_info) && ~isfield(app_info, 'XLim') % if previous plot held a manual zoom state, make sure to update it! 
                    app_info.XLim = max_xlim;
                end
                
            end
            
            if ~isempty(obj.last_ylim) && isequal(obj.show_what, obj.last_shown) && ~cs(obj.show_what, 'stats', 'statistics', 'binning') % only restore the same YLim if viewing the same data type (and not for stats!)
                
                max_ylim = input.ax.YLim; 
                input.ax.YLim = obj.last_ylim;
                
                if ~isempty(app_info) && ~isfield(app_info, 'YLim') % if previous plot held a manual zoom state, make sure to update it! 
                    app_info.YLim = max_ylim;
                end
                
            end
            
            setappdata(input.ax, 'matlab_graphics_resetplotview', app_info);
            
            hold(input.ax, 'off');
            
            input.ax.FontSize = input.font_size;
            
            box(input.ax, 'on'); 
            
            obj.last_shown = obj.show_what;
            
        end
        
        function h = addPlots(obj, ax, data, x_axis, line_str) % internal function for adding another dataset to the plot
            
            if nargin<4 || isempty(x_axis)
                if obj.use_show_datetime && ~isempty(obj.juldates)
                    x_axis = datetime(obj.juldates, 'convertFrom', 'juliandate'); 
                else
                    x_axis = obj.timestamps;
                end
            end
                        
            if nargin<5 || isempty(line_str)
                line_str = '-';
            end
            
            if isempty(data)
                return;
            end
            
            ax.NextPlot = 'add';
            ax.ColorOrderIndex = 1;
            
            if isempty(obj.show_indices)
                idx_vec = 1:obj.show_num_stars;
            else
                idx_vec = obj.show_indices;
            end
            
            for ii = 1:length(idx_vec)
                
                data_star = data(:,idx_vec(ii));
                
%                 if obj.use_skip_flagged
%                     data_star(logical(obj.flags(:,ii,obj.index_flux_number))) = NaN;
%                 end
                
                if obj.use_smooth && size(data,1)>obj.smooth_interval.*4 % at least 5 data points for smoothing! 
                    
                    data_smoothed = util.series.binning(data_star, obj.smooth_interval); 
                    time_smoothed = util.series.binning(x_axis, obj.smooth_interval); 
                    
                    if obj.use_show_log
                        data_smoothed(data_smoothed<0) = 0;
                    end
                    
                    h(ii) = plot(ax, time_smoothed, data_smoothed, line_str, 'LineWidth', 1);
                    
                else
                    
                    if obj.use_show_log
                        data_star(data_star<0) = 0;
                    end
                    
                    h(ii) = plot(ax, x_axis, data_star, line_str, 'LineWidth', 1);
                
                end
                
                if idx_vec(ii)>1 && ~isempty(h(ii)), h(ii).HandleVisibility = 'off'; end

                if ~isempty(h(ii))
                    h(ii).UserData = idx_vec(ii);
                    h(ii).ButtonDownFcn = @obj.callback_cutout_number;
                end
                
                
            end
            
        end
        
        function showStats(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('font_size', 26); 
            input.input_var('ax', [], 'axes', 'axis');
            input.scan_vars(varargin{:});
            
            if isempty(input.ax)
                if ~isempty(obj.gui) && obj.gui.check
                    input.ax = obj.gui.axes_image;
                else
                    input.ax = gca;
                end
            end
            
            flux_index = obj.getFluxTypeIndex(obj.show_flux_type); % which flux to use: raw, sub, rem, or cal
            
            bin_indices = [1 ceil(length(obj.bin_widths_seconds)/2) length(obj.bin_widths_seconds)]; % which bins to use (first, middle and last)
            
            if ~isempty(obj.magnitudes)
                x = obj.magnitudes'; 
                x(x>obj.show_mag_limit) = NaN;
                used_mag = 1;
                x = horzcat(x, nan(1, size(obj.local_RE,2) - size(x,2)));
            else
                x = nanmean(obj.fluxes_cal_); % no magnitudes are given, use the average fluxes instead
                used_mag = 0;
            end

            h1 = plot(input.ax, x, obj.local_RE(bin_indices,:,flux_index), '.', 'HandleVisibility', 'off');

            for ii = 1:length(h1)
                h1(ii).ButtonDownFcn = @obj.callback_local_RE;
                h1(ii).UserData = ii;
            end

            input.ax.ColorOrderIndex = 1;

            input.ax.NextPlot = 'add';

            h2 = plot(input.ax, x, obj.total_RE(bin_indices,:,flux_index), 'p'); 

            for ii = 1:length(h2)
                h2(ii).ButtonDownFcn = @obj.callback_total_RE;
                h2(ii).UserData = ii;
                h2(ii).DisplayName = sprintf('Bin width= %4.2f seconds', obj.bin_widths_seconds(bin_indices(ii)));
            end

            input.ax.NextPlot = 'replace';

            if used_mag
                xlabel(input.ax, 'Magnitude'); 
                input.ax.XScale = 'linear';
                legend(input.ax, 'Location', 'SouthEast'); 
            else
                xlabel(input.ax, 'Flux [counts]'); 
                input.ax.XScale = 'log';
                legend(input.ax, 'Location', 'NorthEast');
            end
            
            ylabel(input.ax, 'Relative Error'); 
            input.ax.YScale = 'log';

        end
        
        function showBinning(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('font_size', 26); 
            input.input_var('ax', [], 'axes', 'axis');
            input.scan_vars(varargin{:});
            
            if isempty(input.ax)
                if ~isempty(obj.gui) && obj.gui.check
                    input.ax = obj.gui.axes_image;
                else
                    input.ax = gca;
                end
            end
            
            flux_index = obj.getFluxTypeIndex(obj.show_flux_type); 
            
            h1 = obj.addPlots(input.ax, obj.total_RE(:,:,flux_index), obj.bin_widths_seconds, 'v');
            xlabel(input.ax, 'Bin size [seconds]'); 
            input.ax.XScale = 'log';

            new_data = obj.total_RE(1,:,flux_index)./sqrt(obj.bin_widths_seconds'./obj.bin_widths_seconds(1)); 

            h2 = obj.addPlots(input.ax, new_data, obj.bin_widths_seconds, ':'); 

            ylabel(input.ax, 'Relative Error'); 
            input.ax.YScale = 'log';

            if obj.use_show_bin_fits
                
                x = log10(obj.bin_widths_seconds');
                fr = util.fit.polyfit(x, log10(obj.total_RE(:,:,flux_index)), 'order', 3);
                
                h3 = obj.addPlots(input.ax, 10.^[fr.ym], obj.bin_widths_seconds, '--');
                
                str = {};
                
                for ii = 1:length(obj.show_indices)
                    h1(ii).HandleVisibility = 'off';
                    h2(ii).HandleVisibility = 'off';
                    h3(ii).HandleVisibility = 'on';
                    
                    slope = fr(obj.show_indices(ii)).coeffs(2); % the linear part of the fit is used at proxy for the slope
                    
                    C = fr(obj.show_indices(ii)).coeffs;
                    N = length(C); 
                    
                    derivative = sum((1:N-1).*C(2:N)'.*(x.^(0:N-2)),2);
                    
                    slope = derivative(1); 
                    
                    str{ii} = sprintf('slope= %4.2f', slope); 
                    
                end
                
                hl = legend(input.ax, str);
                hl.FontSize = 12;
                hl.Location = 'SouthWest';
                
                grid(input.ax, 'on'); 
                
%                 for ii = 1:length(obj.show_indices)
%                     
%                     slope = fr(obj.show_indices(ii)).coeffs(2); % the linear part of the fit is used at proxy for the slope
%                     h = text(obj.bin_widths_seconds(1), obj.total_RE(1,obj.show_indices(ii),flux_index), sprintf('slope= %4.2f', slope));
%                     aspect = log10(input.ax.DataAspectRatio(1)./input.ax.DataAspectRatio(2));
%                     limits = diff(log10(input.ax.YLim))./diff(log10(input.ax.XLim)); 
%                     h.VerticalAlignment = 'top';
%                     h.Rotation = -atand(limits./aspect);
%                     h.FontSize = 14;
%                     
%                 end
                
            end
            
        end
        
        function callback_cutout_number(obj, hndl, ~) % add this callback to each line handle to make it display which star/cutout it represents
            
            if obj.gui.check
                obj.gui.button_cut_number.String = ['star: ' num2str(hndl.UserData)];
            else
                disp(['Star number is ' num2str(hndl.UserData)]);
            end
            
        end
        
        function callback_total_RE(obj, hndl, event) % add this callback to each line handle to make it display the star number of the specific point
            
%             M = event.IntersectionPoint(1);
            RE = event.IntersectionPoint(2);

            flux_index = obj.getFluxTypeIndex(obj.show_flux_type); 

            idx = find(abs((RE-obj.total_RE(hndl.UserData,:,flux_index))./RE)<1e-5);
            
            if length(idx)>1 % two or more stars with the same RE?!
                % check magnitudes or something 
            end
            
            if obj.gui.check
                obj.gui.button_cut_number.String = ['star: ' num2str(idx)];
            else
                disp(['Star number is ' num2str(idx)]);
            end
            
        end
        
        function callback_local_RE(obj, hndl, event) % add this callback to each line handle to make it display the star number of the specific point
            
%             M = event.IntersectionPoint(1);
            RE = event.IntersectionPoint(2);

            flux_index = obj.getFluxTypeIndex(obj.show_flux_type); 

            idx = find(abs((RE-obj.local_RE(hndl.UserData,:,flux_index))./RE)<1e-5);
            
            if length(idx)>1 % two or more stars with the same RE?!
                % check magnitudes or something 
                idx = [];
            end
            
            if obj.gui.check
                obj.gui.button_cut_number.String = ['star: ' num2str(idx)];
            else
                disp(['Star number is ' num2str(idx)]);
            end
            
        end
        
        function plotRelativeError(obj, varargin)
            
            % add varargin parsing
            
            fr = obj.RE_fit_results;
            idx = obj.bad_star_indices; 
            
            loglog(fr.x, fr.y, 'p', sort(fr.x), fr.func(sort(fr.x)), '-r', fr.x(idx), fr.y(idx), 'go'); 
            
        end
        
        function h_out = plotAirmass(obj, varargin)
            
            if isempty(obj.head) || isempty(obj.juldates)
                error('Cannot show the airmass without a header and the julian dates of all observations!'); 
            end
            
            input = util.text.InputVars;
            input.input_var('font_size', 26); 
            input.input_var('ax', [], 'axes', 'axis');
            input.scan_vars(varargin{:});
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            input.ax.NextPlot = 'replace';
            
            t = datetime(obj.juldates, 'convertFrom', 'juliandate');
            
            A = obj.getAirmass;
            
            a = 10.^(0.4.*obj.sysrem_a_vec(:,1));
%             a = obj.sysrem_a_vec(:,1);
            
            if ~isempty(obj.sysrem_a_vec)
                h = plot(input.ax, t, (a-nanmean(a))./nanstd(a).*nanstd(A)+nanmean(A), '.', t, A); 
                hl = legend(input.ax, {'sysrem 1st a vector'}, 'Location', 'NorthWest', 'Airmass'); 
            else
                h = plot(input.ax, t, A); 
            end
            
            input.ax.FontSize = input.font_size;
            
            ylabel(input.ax, 'Airmass'); 
            
            if nargout>0
                h_out = h;
            end
                        
        end
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = img.gui.LightGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
        function spawnFigure(obj)
            
            if obj.gui.check
                
                f = figure;
                f.Name = sprintf('spawned on %s', datetime('now'));
                f.Units = obj.gui.fig.fig.Units;
                f.Position = obj.gui.fig.fig.Position;

                p2 = copyobj(obj.gui.panel_image, f); 
                delete(findobj(p2.Children, 'type', 'UIControl'));
                
                p3 = copyobj(obj.gui.panel_info.panel, f); 
                h = findobj(p3.Children, 'type', 'UIControl'); 
                for ii = 1:length(h)
                    h(ii).Enable = 'inactive';
                end
                
                p4 = copyobj(obj.gui.panel_process.panel, f);
                h = findobj(p4.Children, 'type', 'UIControl'); 
                for ii = 1:length(h)
                    h(ii).Enable = 'inactive';
                end
                
                p5 = copyobj(obj.gui.panel_display.panel, f);
                h = findobj(p5.Children, 'type', 'UIControl'); 
                for ii = 1:length(h)
                    h(ii).Enable = 'inactive';
                end
                
                p6 = copyobj(obj.gui.panel_psd.panel, f);
                h = findobj(p6.Children, 'type', 'UIControl'); 
                for ii = 1:length(h)
                    h(ii).Enable = 'inactive';
                end
                
                p7 = copyobj(obj.gui.panel_stats.panel, f);
                h = findobj(p7.Children, 'type', 'UIControl'); 
                for ii = 1:length(h)
                    h(ii).Enable = 'inactive';
                end
                
                p_tag = uipanel(f, 'Position', obj.gui.panel_control.panel.Position);
                c = uicontrol(p_tag, 'Units', 'Normalized', 'Position', [0 0 1 1], ...
                    'Style', 'text', 'String', sprintf('Plot copy made on:\n %s', datetime('now')),...
                    'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times'); 
                
                if isempty(obj.head)
                    identifier = 'unknown_run';
                else
                    identifier = obj.head.run_identifier;
                end
                                
                f.UserData = sprintf('lightcurves_%s_spawned_%s', identifier, util.text.time2str(datetime('now'))); 
                
                % pop the GUI back on screen! 
                figure(obj.gui.fig.fig); 
                
                % add to list of figures
                idx = cellfun(@isvalid, obj.figures_spawned); 
                obj.figures_spawned = obj.figures_spawned(idx); 
                
                obj.figures_spawned{end+1} = f; 
                
            end
            
        end
        
    end
    
end
