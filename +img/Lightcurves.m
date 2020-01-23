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
% This part needs to be expanded once we write the code! 
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
        
    end
    
    properties % objects
        
        phot_pars; % a struct with some housekeeping about how the photometry was done
        
    end
    
    properties % inputs/outputs
        
        frame_index = 1;
        num_frames = 0;
        
    end
    
    properties % switches/controls
        
        % these controls are used for inputting the data:
        use_double_up = 1; % choose if you want to expand the data storage by factor of 2 each time when space runs out... 
        
        % processing steps for fluxes_sub:
        use_subtract_backgrounds = 1; % 
        
        % processing steps for fluxes_rem:
        use_skip_flagged = 0; % skip samples that have a bad photometry flag (for forced, this may be overkill)
        use_remove_outliers = 1; % remove outliers with "N sigmas" above the noise, calculated by polyfit
        outlier_sigma = 5; % this is NOT used by polyfit, only for making fluxes_sub
        use_skip_bad_times = 0; % if a large fraction of stars have NaN in this epoch, set all to NaN in this epoch
        bad_times_fraction = 0.1; % what fraction of stars need to be NaN to be considered bad times
        
        % processing steps for fluxes_cal:
        use_psf_correction = 1;
        use_polynomial = 1;
        use_zero_point = 1;
        
        index_flux = 'end'; % which aperture to use, default "end" is for forced photometery with the largest aperture
        
        show_for = 'time'; % can also choose "stars" to get star by star stats
        show_what = 'flux'; % can choose "flux", "areas", "backgrounds", "variances", "centroids", "offsets", "widths", "power spectra"
        show_flux_type = 'raw'; % can choose "raw" or "sub" or "rem" or "cal" or "all".
        
        show_num_stars = 10; % up to this number of stars are shown.
        show_indices = []; % if this is not empty, it is used as a list of stars to show
        
        use_smooth = 1; % smoothing is applied to the plotted data only! 
        smooth_interval = 10; % how many samples to average over when smoothing
        
        use_show_log = 0; % show the plot in log scale. Power spectra always shows log Y, but can toggle X. For others, toggles log Y only (X is linear). 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        timestamps;
        
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
        fluxes_cal; % 
        
        fit_results;
        
        power_spectra; 
        psd_freq;
        
    end
    
    properties(Hidden=true)
        
        % These full datasets can include a lot of extra allocated memory, 
        % filled with NaN, to allow less allocations of the data. 
        % To use the actual data call the same names without "_full". 
        timestamps_full;
        
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
        
        fit_results_; % save this from the polyfit (lazy load)
        
        show_what_list = {'fluxes', 'areas', 'backgrounds', 'variances', 'centroids', 'offsets', 'widths', 'bad_pixels', 'power spectra'};
        show_flux_type_list = {'raw', 'sub', 'rem', 'cal', 'all'};
        
        version = 1.07;
        
    end
    
    properties(Hidden=true, Transient=true)
        
        % These are lazy loaded from XXX_full and are cleared whenever the 
        % underlying matrices are changed. This saves time as just slicing
        % into the XXX_full matrices is slow when they are large. 
        % These are not saved if the object is dumped to memory. 
        timestamps_;
        
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
        
    end
    
    methods % constructor
        
        function obj = Lightcurves(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'Lightcurves')
                if obj.debug_bit, fprintf('Lightcurves copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Lightcurves constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj) % remove the long term storage for the entire run
            
            obj.fluxes_full = [];
            obj.errors_full = [];
            obj.timestamps_full = [];
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
            
            obj.clearIntermidiate;
            obj.clearFits;
            obj.clearFluxes;
            obj.clearPSD
            
        end
        
        function clearIntermidiate(obj) % clear all the data sliced from the full matrices
            
            obj.timestamps_ = [];

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
            
            
        end
        
        function clearPSD(obj) % clear the power spectrum and frequency axis

            obj.power_spectra_ = []; 
            obj.psd_freq_ = [];
            
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
                        b = obj.backgrounds(:,:,obj.index_flux_number); 
                        obj.fluxes_sub_ = obj.fluxes_sub_ - a.*b;
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
            else
                
                if isempty(obj.fluxes_rem_)
                    
                    t_rem = tic;
                    
                    obj.fluxes_rem_ = obj.fluxes_sub;
                    
                    if obj.use_skip_flagged
                        
                        obj.fluxes_rem_(logical(obj.flags(:,:,obj.index_flux_number))) = NaN;
                        
                    end
                    
                    if obj.use_remove_outliers
                        
                        r = obj.fit_results;
                        
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
            else
            
                if isempty(obj.fluxes_cal_)
                    
                    t_cal = tic;
                    
                    f = obj.fluxes_rem;

                    if obj.use_psf_correction
                        
                        F = nanmean(obj.fluxes_rem); % average flux is used as weight
                        W = nansum(F.*obj.widths(:,:,obj.index_flux_number),2)./nansum(F,2); % average PSF width, calculated per epoch over all stars
                        
                        Wclip = W;
                        Wclip(abs((W-nanmean(W))./nanstd(W))>3) = NaN; % flag all widths beyond 3sigma from the mean
                        
                        for ii = 1:size(f,2) % each star separately
                        
                            r = util.fit.polyfit(Wclip, f(:,ii), 'double', 1, 'order', 5); % find the coefficients that best correspond W to this star's flux
                        
                            fw = r.func(Wclip); % calculate the best estimate for the flux based on the measure width
                            
                            M = nanmean(f(:,ii)); 
                            
%                             f(:,ii) = f(:,ii) - fw + M; % subtract the fit to the best estimate 
                            
                            rw = util.fit.polyfit(obj.timestamps, fw, 'double', 1, 'order', 5); 
                            f(:,ii) = f(:,ii) - rw.ym + M; % subtract the fit to the best estimate 
                            
                        end
                        
                        
                    end
                    
                    if obj.use_polynomial
                        f = util.series.self_calibrate(f, obj.timestamps, 'zero point', false, 'welch', false, 'polynomial', true, 'order', 5); 
                    end
                    
                    if obj.use_zero_point
                        
                        f_frames_average = nansum(f,2); 
                        f_frames_average_norm = f_frames_average./nanmean(f_frames_average); 
                        f = f./f_frames_average_norm;

                    end
                    
                    obj.fluxes_cal_ = f;
                    
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
            
            val = obj.centroids_y_full;
            
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
            
            if isempty(obj.fit_results_)
                obj.fit_results_ = util.fit.polyfit(obj.timestamps, obj.fluxes_sub, 'order', 5, 'sigma', 3, 'double', 1); % must use double precision for the sizes and length of these fluxes!
            end
            
            val = obj.fit_results_;
            
        end
        
        function val = get.power_spectra(obj)
            
            if isempty(obj.fluxes_sub)
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
            
            if isempty(obj.fluxes_sub)
                val = [];
                return; 
            end
            
            if isempty(obj.power_spectra_) || isempty(obj.psd_freq_)
                disp('calculating power spectrum (from get.power_spectra)'); 
                [obj.power_spectra_, obj.psd_freq_] = obj.calculatePSD; 
            end
            
            val = obj.psd_freq_; 
            
        end
        
        function [pxx, freq] = calculatePSD(obj)
            
            f = obj.fluxes_sub;
            f = fillmissing(f, 'linear'); 
            
            [pxx, freq] = pwelch(f, [], [], [], 1./median(diff(obj.timestamps))); 
            
        end
        
        function str_out = print_pars(obj)

            import util.text.print_vec;

            str = sprintf('Nframes: %d | Nstars: %d | Cutous: %dx%d | ap: %s | ann: %s | f.ap: %s | gauss: %4.2f | iter: %d | shift_res: %4.2f', ...
                obj.num_frames, obj.phot_pars.cutout_size(4), obj.phot_pars.cutout_size(1), obj.phot_pars.cutout_size(2), ...
                print_vec(obj.phot_pars.aperture_radii), print_vec(obj.phot_pars.annulus_radii), print_vec(obj.phot_pars.forced_radius), ...
                obj.phot_pars.gauss_sigma, obj.phot_pars.iterations, obj.phot_pars.shift_resolution); 

            if nargout>0
                str_out = str;
            else
                disp(str); 
            end

        end

    end
    
    methods % setters
        
        function set.frame_index(obj, val)
            
            if ~isequal(val, obj.frame_index)
                
                obj.frame_index = val;
                
                obj.clearIntermidiate;
                obj.clearPSD;
                
            end
            
        end
        
        function set.fluxes_full(obj, val)
            
            obj.fluxes_full = val;
            obj.clearFits;
            obj.clearFluxes;
            obj.clearPSD;
            
        end
        
        function set.timestamps_full(obj, val)
            
            obj.timestamps_full = val;
            obj.clearFits;
            obj.clearFluxes;
            obj.clearPSD;
            
        end
        
        function set.backgrounds_full(obj, val)
            
            obj.backgrounds_full = val;
            obj.clearFluxes;
            obj.clearPSD;
            
        end
        
        function set.areas_full(obj, val)
            
            obj.areas_full = val;
            obj.clearFluxes;
            obj.clearPSD;
            
        end
        
        function set.use_subtract_backgrounds(obj, val)
            
            if ~isequal(obj.use_subtract_backgrounds, val)
                obj.use_subtract_backgrounds = val;
                obj.clearFluxes;
                obj.clearPSD;
            end
            
        end
        
        function set.use_skip_flagged(obj, val)
            
            if ~isequal(obj.use_skip_flagged, val)
                obj.use_skip_flagged = val;
                obj.fluxes_rem_ = [];
                obj.fluxes_cal_ = [];
                obj.clearPSD;
            end
            
        end
        
        function set.use_remove_outliers(obj, val)
            
            if ~isequal(obj.use_remove_outliers, val)
                obj.use_remove_outliers = val;
                obj.fluxes_rem_ = [];
                obj.fluxes_cal_ = [];
                obj.clearPSD;
            end
            
        end
        
        function set.outlier_sigma(obj, val)
            
            if ~isequal(obj.outlier_sigma, val)
                obj.outlier_sigma = val;
                obj.fluxes_rem_ = [];
                obj.fluxes_cal_ = [];
                obj.clearPSD;
            end
            
        end
        
        function set.use_skip_bad_times(obj, val)
            
            if ~isequal(obj.use_skip_bad_times, val)
                obj.use_skip_bad_times = val;
                obj.fluxes_rem_ = [];
                obj.fluxes_cal_ = [];
                obj.clearPSD;
            end
            
        end
        
        function set.bad_times_fraction(obj, val)
            
            if ~isequal(obj.bad_times_fraction, val)
                obj.bad_times_fraction = val;
                obj.fluxes_rem_ = [];
                obj.fluxes_cal_ = [];
                obj.clearPSD;
            end
            
        end
        
        function set.use_psf_correction(obj, val)
            
            if ~isequal(obj.use_psf_correction, val)
                obj.use_psf_correction = val;
                obj.fluxes_cal_ = [];
                obj.clearPSD;
            end
            
        end
        
        function set.use_polynomial(obj, val)
            
            if ~isequal(obj.use_polynomial, val)
                obj.use_polynomial = val;
                obj.fluxes_cal_ = [];
                obj.clearPSD;
            end
            
        end
        
        function set.use_zero_point(obj, val)
            
            if ~isequal(obj.use_zero_point, val)
                obj.use_zero_point = val;
                obj.fluxes_cal_ = [];
                obj.clearPSD;
            end
            
        end
        
        function set.index_flux(obj, val)
            
            if ~isequal(obj.index_flux, val)
                obj.index_flux = val;
                obj.clearFluxes;
                obj.clearPSD;
            end
            
        end
        
    end
    
    methods % calculations
        
        function input(obj, varargin) % give the photometric results explicitely, they are added to the storage
            
            import util.vec.insert_matrix;
            
            if isscalar(varargin) && isa(varargin{1}, 'util.text.InputVars')
                input = varargin{1};
            else
                input = util.text.InputVars;
                input.use_ordered_numeric = 1;
                input.input_var('timestamps', [], 'times');
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
            
            obj.fluxes_full = insert_matrix(obj.fluxes_full, input.fluxes, [obj.frame_index,1,1], NaN, obj.use_double_up);
            obj.errors_full = insert_matrix(obj.errors_full, input.errors, [obj.frame_index,1,1], NaN, obj.use_double_up);
            obj.timestamps_full = insert_matrix(obj.timestamps_full, input.timestamps, [obj.frame_index,1,1], NaN, obj.use_double_up);
            obj.areas_full = insert_matrix(obj.areas_full, input.areas, [obj.frame_index,1,1], NaN, obj.use_double_up);
            obj.backgrounds_full = insert_matrix(obj.backgrounds_full, input.backgrounds, [obj.frame_index,1,1], NaN, obj.use_double_up);
            obj.variances_full = insert_matrix(obj.variances_full, input.variances, [obj.frame_index,1,1], NaN, obj.use_double_up);
            obj.centroids_x_full = insert_matrix(obj.centroids_x_full, input.centroids_x, [obj.frame_index,1,1], NaN, obj.use_double_up);
            obj.centroids_y_full = insert_matrix(obj.centroids_y_full, input.centroids_y, [obj.frame_index,1,1], NaN, obj.use_double_up);
            obj.offsets_x_full = insert_matrix(obj.offsets_x_full, input.offsets_x, [obj.frame_index,1,1], NaN, obj.use_double_up);
            obj.offsets_y_full = insert_matrix(obj.offsets_y, input.offsets_y, [obj.frame_index,1,1], NaN, obj.use_double_up);
            obj.widths_full = insert_matrix(obj.widths_full, input.widths, [obj.frame_index,1,1], NaN, obj.use_double_up);
            obj.bad_pixels_full = insert_matrix(obj.bad_pixels_full, input.bad_pixels, [obj.frame_index,1,1], NaN, obj.use_double_up);
            obj.flags_full = insert_matrix(obj.flags_full, input.flags, [obj.frame_index,1,1], NaN, obj.use_double_up);
            obj.frame_index = obj.frame_index + N;
            
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
        
        function saveAsMAT(obj, filename, star_indices, frame_indices) % save only the data, in raw matrices, to a MAT file, along with phot_pars and a short readme. 
            
            if nargin<2 || isempty(filename)
                disp('Must supply a filename! Usage: obj.saveAsMAT(filename, star_indices, frame_indices'); 
            end
            
            if nargin<3 || isempty(star_indices)
                star_indices = [];
            end
            
            if nargin<4 || isempty(frame_indices)
                frame_indices = [];
            end
            
            timestamps = obj.timestamps;
            
            if isempty(star_indices)
                fluxes = obj.fluxes;
                errors = obj.errors;
                areas = obj.areas;
                backgrounds = obj.backgrounds;
                variances = obj.variances;
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
                centroids_x = centroids_x(frame_indices,:,:);
                centroids_y = centroids_y(frame_indices,:,:);
                widths = widths(frame_indices,:,:);
                bad_pixels = bad_pixels(frame_indices,:,:);
                flags = obj.flags(frame_indices,:,:);
                
                timestamps = timestamps(frame_indices,:);
            
            end
            
            readme = 'Some info about the data stored in this file:';
            readme = sprintf('%s\n *timestamps: a column vector with %d elements, with the relative time of each frame', readme, size(timestamps, 1));
            readme = sprintf('%s\n *fluxes: the actual lightcurves, dim 1 is %d frames dim 2 is %d stars, dim 3 is for the different apertures+forced aperture',readme, size(fluxes,1), size(fluxes,2));
            readme = sprintf('%s\n *errors: calculated from the variance map and source intensity. ', readme); 
            readme = sprintf('%s\n *areas: size of aperture, removing bad pixels and so on.', readme);
            readme = sprintf('%s\n *backgrounds: measured intensity, per pixel, in the annulus.', readme);
            readme = sprintf('%s\n *variances: measured noise variance, per pixel, in the annulus.', readme);
            readme = sprintf('%s\n *centroids_x/y: position within image (in pixels)', readme);
            readme = sprintf('%s\n *widths: calculated using the 2nd moment, averaging the minor and major axis', readme);
            readme = sprintf('%s\n *bad_pixels: number of bad pixels in the aperture area.', readme);
            readme = sprintf('%s\n *flags: this is marked as non-zero if the photometry did not go well (e.g., large offsets, negative widths)', readme); 
            
            save(filename, 'timestamps', 'fluxes', 'errors', 'areas', 'backgrounds', 'variances',...
                'centroids_x', 'centroids_y', 'widths', 'bad_pixels', 'readme', '-v7.3');
            
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
            
            delete(allchild(input.ax));
            delete(findobj(input.ax.Parent, 'type', 'legend'));
            
            if cs(obj.show_what, 'power spectra')
                
                input.ax.YScale = 'log';
                
                if obj.use_show_log
                    input.ax.XScale = 'log';
                else
                    input.ax.XScale = 'linear';
                end
                
            else
                
                input.ax.XScale = 'linear';
                
                if obj.use_show_log
                    input.ax.YScale = 'log';
                else
                    input.ax.YScale = 'linear';
                end
                
            end
            
            xlabel(input.ax, 'time (seconds)');
            
            if cs(obj.show_what, 'fluxes')
                
                if cs(obj.show_flux_type, 'raw')
                    obj.addPlots(input.ax, obj.fluxes);
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
                
                obj.addPlots(input.ax, obj.power_spectra, obj.psd_freq);
                ylabel(input.ax, 'power spectra (Welch)'); 
                xlabel(input.ax, 'Frequency [Hz]'); 
                
            else
                error('Unknown data to show "%s", use "fluxes" or "offset" etc...', obj.show_what);
            end
            
            hold(input.ax, 'off');
            
            input.ax.FontSize = input.font_size;
            
            
            
            box(input.ax, 'on'); 
            
        end
        
        function addPlots(obj, ax, data, x_axis, line_str) % internal function for adding another dataset to the plot
            
            if nargin<4 || isempty(x_axis)
                x_axis = obj.timestamps;
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
            
            for ii = idx_vec
                
                data_star = data(:,ii);
                
%                 if obj.use_skip_flagged
%                     data_star(logical(obj.flags(:,ii,obj.index_flux_number))) = NaN;
%                 end
                
                if obj.use_smooth
                    
                    data_smoothed = util.series.binning(data_star, obj.smooth_interval); 
                    time_smoothed = util.series.binning(x_axis, obj.smooth_interval); 
                    
                    if obj.use_show_log
                        data_smoothed(data_smoothed<0) = 0;
                    end
                    
                    h = plot(ax, time_smoothed, data_smoothed, line_str, 'LineWidth', 1);
                    if ii>min(idx_vec) && ~isempty(h), h.HandleVisibility = 'off'; end
                    
                else
                    
                    if obj.use_show_log
                        data_star(data_star<0) = 0;
                    end
                    
                    h = plot(ax, x_axis, data_star, line_str, 'LineWidth', 1);
                    if ii>min(idx_vec) && ~isempty(h), h.HandleVisibility = 'off'; end                    
                
                end
                
                if ~isempty(h)
                    h.UserData = ii;
                    h.ButtonDownFcn = @obj.callback_cutout_number;
                end
                
                
            end
            
        end
        
        function callback_cutout_number(obj, hndl, ~) % add this callback to each line handle to make it display which star/cutout it represents
            
            if obj.gui.check
                obj.gui.button_cut_number.String = ['star: ' num2str(hndl.UserData)];
            else
                disp(['Star number is ' num2str(hndl.UserData)]);
            end
            
        end
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = img.gui.LightGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
    end
    
end
