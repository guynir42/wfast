classdef QualityChecker < handle
% Data quality checker class. 
% This object is fed data from the DataStore and uses a bunch of cuts to 
% define what region of each lightcurve is useful or disqualified. 
% 
% The only way to input data to this object is by giving it a DataStore 
% object as the sole argument of the input() method. 
%
% This object has two purposes: 
% (a) To check the data and flag only the useful data that should be 
%     counted as good star-time. These values are reported to the StarHour
%     object and counted. This is important so we know how much time we 
%     really scanned for occultations. 
%     This can be used even when not doing full event-finding, only checking
%     the number of star hours. 
% (b) To disqualify event candidates that occur during bad times in the data. 
%     Each candidate that crosses the threshold in the EventFinder can be
%     checked agains the quality cuts in this object. If the peak of the 
%     event coincides with times/stars that are marked as bad, then the 
%     candidate is marked as un-kept. 
%
% The user-defined parameters of this class are saved as a struct "pars". 
% These parameters are defined in the "reset/clear" methods block, inside
% the resetPars() function. Refer to the inline docs for info on each parameter. 
%
% Some notable parameters are:
% -use_auto_count_hours: automatically add the new data into the StarHours 
%  object. This is useful if you are only scanning the data quality. When
%  using this object in event-finding, you may want to have control over if
%  a batch is added to the star hours (based on not having many candidates). 
%  The default is false (as most reasonable for event finding). 
% -use_dilate and dilate_region: any bad regions in the data are expanded
%  by this amount on either side, representing the fact that even if the 
%  event peaks close to a bad region we still disqualify it. 
%  Defaults are true and 5. 
% -use_subtract_mean_offsets: calculate the mean offsets and then subtract
%  them from each star's individual offset. This is useful as it reveals 
%  the real star's motion, and removes the tracking errors from this estimate. 
%  Default is true. 
% -use_XXX and thresh_XXX: these control the different cuts. A few exceptions
%  must be noted: (a) some cuts don't have a threshold. They are logical 
%  cuts, whenever the value is non-zero they are considered to be triggered. 
%  (b) The correlations cuts are all treated with a single use_correlations
%  and thresh_correlations, event though they include multiple cuts. 
% -corr_types and corr_timescales: the correlation cuts correlate the flux
%  to one of the types of auxiliary (default is b, x, y, r, and w). 
%  The timescales denote how many frames are used in calculating those 
%  correlations, defaulting to 10, 25, and 50. 
% -bad_rows and bad_columns: indices of these rows and columns mark bad 
%  regions in the sensor. They can be filled automatically using the 
%  setupSensor() function. Make sure the header is up to date and contains
%  the correct camera in the INST field. 
%  
% The cut results are stored in matrices where dim1 is time, dim2 is star
% index, and dim3 is the cut index. The index corresponding to each cut is
% defined in "cut_names" and "cut_indices". You can find the right cut:
% >> offset_size_cut = checker.cut_values_matrix(:,:,checker.cut_indices.offset_size); 
%
% The cut_flag_matrix is a logical matrix that shows where the cuts exceeded
% their respective threshold. This region is also expanded if using dilate. 
% The logical sum of the 3rd dimension of this matrix is accessible using
% the bad_times() function, which tells for each frame and star if that data
% is usable or disqualified. 
% 
% A note on histograms: this object also keeps a histogram of each cut, with
% the number of times each value range was measured, for each star separately. 
% This is used for debugging and optimizing the cuts. 
% To view these histograms use the plotting tool showHistogram(). 
% Some parameters for this function are:
% -type: which cut should be displayed. Use the index or name. Default 1. 
% -sum: if false, show the values for each star. If true, sum the data
%       into a 1D histogram, showing the values for all stars together. 
%       Default is false (show each star separately). 
% -axes: which axes to draw into. Default is gca(). 
% -font_size: to use on the axes. Default is 18. 
% -log: show the values on a logarithmic scale. Default is true. 
% 


    properties(Transient=true)
        
    end
    
    properties % objects
        
        head@head.Header; 
        hours@tno.StarHours;
        pars; % a struct with all the user-defined parameters
        
    end
    
    properties % inputs/outputs
        
        correlations; % 4D matrix: dim1 is frame index, dim2 is star index, dim3 is aux name, dim4 is timescale
        corr_indices; % generate a struct using setupCuts()
        
        delta_t; % difference between average time-step (dt) and average time-step (usually 0.04s). 
        shakes; % the offset size for each star
        defocus; % the size of the PSF
        slope; % the slope at different points in the flux
        offset_size; % the size of the offset in each star
        
        nan_flux; % places where the flux is NaN
        nan_offsets; % places where the offset_x or offset_y is nan
        photo_flag; % places flagged by photometry (mainly offsets and widths that are NaN)
        near_bad_rows_cols; % flag places where the centroids are too close to a bad row/column
        linear_motion; % are the x/y offsets moving in a linear line across the cutout (like a satellite?)
        background_intensity; % the background per pixel value (remove frames with star in annulus or when the sky is too bright)
        bad_pixels; % how many bad pixels are in the aperture
        repeating_columns; % how many columns are repeated for each frame/star
        
        flux_corr; % correlations of flux between stars
        
        mean_x; % the weighted average offset for all stars
        mean_y; % the weighted average offset for all stars
        
        % rename flags to cuts? make it 2D single matrix
        % dim1 of flag is frame index, dim2 is star index 
        % flags show what kind of bad time affects each frame
        % these are logicals and are summed up (ORed) to give "bad_times"
        cut_flag_matrix; % include all flags in this big matrix, dim1 frame number, dim2 star index, dim3 which cut is used
        cut_values_matrix; % these are the actual measurements for the different cuts, same dimensions as flags
        cut_thresholds; % set the threshold for each cut in one place
        cut_two_sided; % mark each cut if it is only for large positive or both directions
        cut_names = {}; % automatically assemble these at startup
        cut_indices;  % generate a struct using setupCuts()
        
        % these we get from the DataStore
        background_flux; % cut out the background flux for calculating the variance outside the search region
        background_detrend; % same flux, with a linear fit subtracted from each batch individually
        background_aux; % cut out the background aux for calculating the variance outside the search region
        background_timestamps; % timestamps for the background region
        
        extended_flux; % cutout of the flux from the flux_buffer extended around the search region
        extended_aux; % cutout of the aux from the aux_buffer extended around the search region 
        extended_detrend; % same flux, with a linear fit subtracted from each batch individually
        extended_timestamps; % timestamps for the extended batch region
        cutouts; % cutouts for the extended region
        
        search_start_idx; % starting index for the search region out of the EXTENDED BATCH!
        search_end_idx;  % end index for the search region out of the EXTENDED BATCH!
        
        search_flux; % the flux in the search region
        search_detrend; % same flux, with a linear fit subtracted from each batch individually
        search_aux; % the aux in the search region
        search_timestamps;  % timestamps for the search region
        search_juldates; % julian dates for each timestamp in the search region
        
        aux_names; % get the names (and indices below) from DataStore
        aux_indices; % struct with field=number for each of the above aux names
        
        histograms = []; % cut values accumulated over the run, with dim1 the edges defined for each cut type, dim2 the star number and dim3 is the cut type
        hist_edges = []; % edges for each cut (each row is edges for a different cut)
        
        juldate_log = []; % keep the julian date for each of the defocus log measurements
        defocus_log = []; % track the defocus result for each batch
        mean_width_values; % track the mean width for all batches in this run
        mean_background_values;  % track the mean background for all batches in this run
        
        fwhm; % from the ModelPSF (in arcsec)
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        num_frames;
        num_stars;
        num_cuts;
        
    end
    
    properties(Hidden=true)
        
        version = 1.02;
        
    end
    
    methods % constructor
        
        function obj = QualityChecker(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'tno.QualityChecker')
                if obj.debug_bit>1, fprintf('QualityChecker copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('QualityChecker constructor v%4.2f\n', obj.version); end

                obj.hours = tno.StarHours; 
        
                obj.resetPars;
                
                obj.setupSensor; 
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function resetPars(obj)
           
            obj.pars = struct;
            
            % if true, will automatically add new run time on each call to input() 
            % (default is false, we want to manually add hours after checking the batch is good!)
            obj.pars.use_auto_count_hours = false; 

            obj.pars.use_dilate = true; % do we want to spread out bad regions a few frames in either direction? 
            obj.pars.dilate_region = 10; % how many frames, on either side, do we flag next to each bad point? 
            
            obj.pars.use_subtract_mean_offsets = true; % before doing any calculations on the offsets, remove the mean offsets that also affect the forced photometry centroids

            % do we want to apply all these cuts? 
            obj.pars.use_delta_t = true;
            obj.pars.use_shakes = true;
            obj.pars.use_defocus = false;
            obj.pars.use_fwhm = true; 
            obj.pars.use_slope = true;
            obj.pars.use_near_bad_rows_cols = true;            
            obj.pars.use_offset_size = true;
            obj.pars.use_linear_motion = true;
            obj.pars.use_background_intensity = true;
            obj.pars.use_nan_flux = false;
            obj.pars.use_nan_offsets = false;
            obj.pars.use_photo_flag = false; 
            obj.pars.use_bad_pixels = true; 
            obj.pars.use_repeating_columns = true; 
            obj.pars.use_flux_corr = true;
            obj.pars.use_correlations = true;

            obj.pars.corr_types = {'b', 'x', 'y', 'r', 'w'}; % types of auxiliary we will use for correlations (r is derived from x and y)
            obj.pars.corr_timescales = [10, 25, 50]; % number of frames to run each correlation sum

            obj.pars.thresh_delta_t = 0.3; % events where the difference in timestamps, relative to the mean time-step are disqualified
            obj.pars.thresh_shakes = 5; % events where the mean offset r is larger than this are disqualified
            obj.pars.thresh_defocus = 3; % events with PSF width above this value are disqualified
            obj.pars.thresh_fwhm = 10; % batches with too big FWHM (from ModelPSF, in arcsec) are disqualified
            obj.pars.thresh_slope = 5; % events where the slope is larger than this value (in abs. value) are disqualified
            obj.pars.thresh_offset_size = 4; % events with offsets above this number are disqualified (after subtracting mean offsets)
            obj.pars.thresh_linear_motion = 2; % events showing linear motion of the centroids are disqualified
            obj.pars.thresh_background_intensity = 10; % events where the background per pixel is above this threshold are disqualified
            obj.pars.thresh_flux_corr = 3.0; % events where the flux of one star has 95% percentile higher than this are disqualified
            obj.pars.thresh_correlation = 4; % correlation max/min of flux (with e.g., background) with value above this disqualifies the region

            obj.pars.thresh_tracking_error = 4; % for each event, post-detection, check correlation with other stars
            
            obj.pars.smoothing_slope = 50; % number of frames to average over when calculating slope
            obj.pars.distance_bad_rows_cols = 5; % how many pixels away from a bad row/column would we still disqualify a star? (on either side, inclusive)
            obj.pars.linear_timescale = 25; % timescale for linear_motion cut

            obj.pars.num_hist_edges = 200; % including 100 negative and 100 positive values

            obj.pars.bad_columns = []; % which columns are considered bad
            obj.pars.bad_rows = []; % which rows are considered bad

            obj.pars.num_stars_defocus = 100; % how many stars (max) to use for calculating defocus
            
            obj.setupSensor;
            
            obj.reset; 
            
        end
        
        function reset(obj)
            
            obj.hours.reset;
            
            obj.setupCuts;
            
            obj.hist_edges = [];
            obj.histograms = [];
            
            obj.defocus_log = [];
            obj.mean_width_values = [];
            obj.mean_background_values = [];
            
            obj.clear;
            
        end
        
        function setupCuts(obj)
            
            obj.corr_indices = struct;
            for ii = 1:length(obj.pars.corr_types)
                obj.corr_indices.(obj.pars.corr_types{ii}) = ii; 
            end
            
            obj.cut_names = {}; 
            obj.cut_thresholds = [];
            obj.cut_two_sided = logical.empty;
            
            if obj.pars.use_delta_t
                obj.cut_names{end+1} = 'delta_t';
                obj.cut_thresholds(end+1) = obj.pars.thresh_delta_t;
                obj.cut_two_sided(end+1) = true;
            end
            
            if obj.pars.use_shakes 
                obj.cut_names{end+1} = 'shakes'; 
                obj.cut_thresholds(end+1) = obj.pars.thresh_shakes;
                obj.cut_two_sided(end+1) = false;
            end
            
            if obj.pars.use_defocus 
                obj.cut_names{end+1} = 'defocus'; 
                obj.cut_thresholds(end+1) = obj.pars.thresh_defocus;
                obj.cut_two_sided(end+1) = false;
            end
            
            if obj.pars.use_fwhm
                obj.cut_names{end+1} = 'fwhm'; 
                obj.cut_thresholds(end+1) = obj.pars.thresh_fwhm;
                obj.cut_two_sided(end+1) = false;
            end
            
            if obj.pars.use_slope
                obj.cut_names{end+1} = 'slope'; 
                obj.cut_thresholds(end+1) = obj.pars.thresh_slope;
                obj.cut_two_sided(end+1) = true;
            end
            
            if obj.pars.use_near_bad_rows_cols 
                obj.cut_names{end+1} = 'near_bad_rows_cols'; 
                obj.cut_thresholds(end+1) = 1; 
                obj.cut_two_sided(end+1) = false;
            end
            
            if obj.pars.use_offset_size
                obj.cut_names{end+1} = 'offset_size'; 
                obj.cut_thresholds(end+1) = obj.pars.thresh_offset_size;
                obj.cut_two_sided(end+1) = false;
            end
            
            if obj.pars.use_linear_motion
                obj.cut_names{end+1} = 'linear_motion'; 
                obj.cut_thresholds(end+1) = obj.pars.thresh_linear_motion; 
                obj.cut_two_sided(end+1) = false;
            end
            
            if obj.pars.use_background_intensity
                obj.cut_names{end+1} = 'background_intensity'; 
                obj.cut_thresholds(end+1) = obj.pars.thresh_background_intensity;
                obj.cut_two_sided(end+1) = false;
            end
            
            if obj.pars.use_nan_flux 
                obj.cut_names{end+1} = 'nan_flux'; 
                obj.cut_thresholds(end+1) = 1;
                obj.cut_two_sided(end+1) = false;
            end
            
            if obj.pars.use_nan_offsets 
                obj.cut_names{end+1} = 'nan_offsets'; 
                obj.cut_thresholds(end+1) = 1;
                obj.cut_two_sided(end+1) = false;
            end
            
            if obj.pars.use_photo_flag 
                obj.cut_names{end+1} = 'photo_flag'; 
                obj.cut_thresholds(end+1) = 1;
                obj.cut_two_sided(end+1) = false;
            end
            
            if obj.pars.use_bad_pixels
                obj.cut_names{end+1} = 'bad_pixels'; 
                obj.cut_thresholds(end+1) = 1;
                obj.cut_two_sided(end+1) = false;
            end
            
            if obj.pars.use_repeating_columns
                obj.cut_names{end+1} = 'repeating_columns'; 
                obj.cut_thresholds(end+1) = 1; % maybe at some point we may want to only trigger on multiple repeated rows per cutout? 
                obj.cut_two_sided(end+1) = false;
            end
            
            if obj.pars.use_flux_corr
                
                obj.cut_names{end+1} = 'flux_corr_25'; 
                obj.cut_thresholds(end+1) = obj.pars.thresh_flux_corr;
                obj.cut_two_sided(end+1) = false;
                
                obj.cut_names{end+1} = 'flux_corr_50'; 
                obj.cut_thresholds(end+1) = obj.pars.thresh_flux_corr;
                obj.cut_two_sided(end+1) = false;
                
            end
            
            for ii = 1:length(obj.pars.corr_types)
                
                for jj = 1:length(obj.pars.corr_timescales)
                    name = sprintf('corr_%s_%d', obj.pars.corr_types{ii}, obj.pars.corr_timescales(jj)); % e.g., corr_a_25 or corr_x_100 
                    obj.cut_names{end+1} = name;
                    obj.cut_thresholds(end+1) = obj.pars.thresh_correlation; 
                    obj.cut_two_sided(end+1) = true;
                end
                
            end
            
            obj.cut_indices = struct;
            for ii = 1:length(obj.cut_names)
                obj.cut_indices.(obj.cut_names{ii}) = ii; 
            end
            
        end
        
        function clear(obj)
            
            obj.background_flux = [];
            obj.background_detrend = [];
            obj.background_aux = [];
            obj.background_timestamps = [];

            obj.extended_flux = [];
            obj.extended_detrend = [];
            obj.extended_aux = [];
            obj.extended_timestamps = [];
            obj.cutouts = []; 
            
            obj.search_start_idx = [];
            obj.search_end_idx = [];

            obj.search_flux = [];
            obj.search_detrend = [];
            obj.search_aux = [];
            obj.search_timestamps = [];
            obj.search_juldates = [];
            
            obj.aux_names = [];
            obj.aux_indices = [];
            
            obj.correlations = [];
            obj.cut_flag_matrix = [];            
            obj.cut_values_matrix = [];
            
            obj.delta_t = []; 
            obj.shakes = [];
            obj.defocus = []; 
            obj.slope = [];
            obj.offset_size = [];
            obj.nan_flux = [];
            obj.nan_offsets = [];
            obj.near_bad_rows_cols = []; 
            
            obj.mean_x = [];
            obj.mean_y = [];
            
        end
        
    end
    
    methods % getters
        
        function val = get.num_frames(obj)
            
            if isempty(obj.extended_flux)
                val = 0;
            else
                val = size(obj.extended_flux,1); 
            end
            
        end
        
        function val = get.num_stars(obj)
            
            if isempty(obj.extended_flux)
                val = 0;
            else
                val = size(obj.search_flux,2); 
            end
            
        end
        
        function val = get.num_cuts(obj)
            
            if isempty(obj.cut_names)
                val = 0;
            else
                val = length(obj.cut_names); 
            end
            
        end
        
        function val = bad_times(obj)
            
            val = logical(sum(obj.cut_flag_matrix,3)); 
            
        end
       
        function val = getCorrNames(obj)
            
            val = {}; 
            
            for ii = 1:length(obj.pars.corr_types)
                
                for jj = 1:length(obj.pars.corr_timescales)
                    
                    val{end+1} = sprintf('corr_%s_%d', obj.pars.corr_types{ii}, obj.pars.corr_timescales(jj)); 
                    
                end
                
            end
            
        end
        
    end
    
    methods % setters
        
        function set.head(obj, val)
            
            obj.head = val;
            
            obj.hours.head = val;
            
            obj.setupSensor;
            
        end
        
        function setupSensor(obj, camera)
            
            if nargin<2 || isempty(camera)
                if ~isempty(obj.head) && ~isempty(obj.head.INST)
                    camera = obj.head.INST;
                else                
                    camera = 'Zyla'; 
                end
            end
            
            obj.pars.bad_columns = [];
            obj.pars.bad_rows = [];
            
            if util.text.cs('Balor', camera)
                obj.pars.bad_columns = [1:20, 1960, 4085:4104]; 
                obj.pars.bad_rows = [1 1161 1765 3716 4128]; 
            elseif util.text.cs('Zyla', camera)
                % I am not sure we have any bad rows/columns in the Zyla...
                % instead just put the first and last row/column. 
                obj.pars.bad_columns = [1 2156];
                obj.pars.bad_rows = [1 2560]; 
            else
                error('Unknown camera "%s". Use "Balor" or "Zyla...', camera);
            end
            
        end
        
    end
    
    methods % calculations
        
        function input(obj, store) % do all the calculations
            
            if nargin<2 || isempty(store) || ~isa(store, 'tno.DataStore')
                error('Must supply a valid DataStore!'); 
            end
            
            obj.clear;
            
            obj.ingestStore(store); 
            
            %%%%%%%%%%%%% load the data into shorthands %%%%%%%%%%%%%%%%%%%
            
            t = obj.extended_timestamps; 
            f = obj.extended_flux; 
            
            a = obj.extended_aux(:,:,obj.aux_indices.areas); 
            b = obj.extended_aux(:,:,obj.aux_indices.backgrounds); 
            x = obj.extended_aux(:,:,obj.aux_indices.offsets_x); 
            y = obj.extended_aux(:,:,obj.aux_indices.offsets_y); 
            X = obj.extended_aux(:,:,obj.aux_indices.centroids_x); 
            Y = obj.extended_aux(:,:,obj.aux_indices.centroids_y); 
            w = obj.extended_aux(:,:,obj.aux_indices.widths); 
            p = obj.extended_aux(:,:,obj.aux_indices.bad_pixels); 
            
            F = nanmean(f - a.*b,1); % mean flux, after subtracting the background
            
            %%%%%%%%%%%%% calculate / prepare the data %%%%%%%%%%%%%%%%%%%%
            
            idx = obj.search_start_idx:obj.search_end_idx; % indices inside the search region
            
            obj.mean_x = util.vec.weighted_average(x,F.^2,2); % maybe use square of flux instead??
            obj.mean_y = util.vec.weighted_average(y,F.^2,2); % last arg is for dimension 2 (average over stars)

            if obj.pars.use_defocus
                obj.defocus = obj.calculateDefocus; 
                obj.defocus_log = vertcat(obj.defocus_log, obj.defocus); 
            end
            
            obj.juldate_log = vertcat(obj.juldate_log, nanmean(obj.search_juldates)); % keep track of the julian date of each batch
            
            W = util.vec.weighted_average(w,F.^2,2);
            obj.mean_width_values = vertcat(obj.mean_width_values, W(idx)); % keep a log of the mean widths for the entire run
            obj.mean_background_values = vertcat(obj.mean_background_values, nanmedian(b(idx,:),2)); % keep a log of the sky background level for the entire run
            
            if obj.pars.use_subtract_mean_offsets
                x = x - obj.mean_x; 
                y = y - obj.mean_y;
            end

            w = fillmissing(w, 'spline'); 
            
            r = sqrt(x.^2 + y.^2); % size of offsets
            
            rr = cat(3,r,circshift(r,1),circshift(r,-1)); % pile up the r matrices with a small shift 
            obj.offset_size = min(rr,[],3); % find the minimal offset in unshifted, and slightly shifted matrices (e.g., the smallest value out of nearest neighbors)
            
            obj.shakes = util.vec.weighted_average(r,F.^2,2); % the average offset size (weighted by the flux squared, making it more reliable as it uses bright stars)

            mean_dt = median(diff(t)); % get the cadence (usually this is 0.04 seconds)
            dt_vector = [mean_dt; diff(t)]; % add one mean_dt in the first place to make dt_vector the same size as t
            
            obj.delta_t = (dt_vector - mean_dt)./mean_dt;
            
            k = (-obj.pars.smoothing_slope/2:obj.pars.smoothing_slope/2)'; 
            k = k./sqrt(sum(k.^2)); 
            
            FN = (nanmean(f,2)-nanmean(F))./nanstd(F); 
            
            obj.slope = filter2(k, FN); 

            obj.nan_flux = isnan(f); 
            obj.nan_offsets = isnan(x) | isnan(y);  
            obj.photo_flag = obj.extended_aux(:,:,obj.aux_indices.flags); 
            
            obj.near_bad_rows_cols = false(size(X)); % initialize to all zero
            
            for ii = 1:length(obj.pars.bad_rows)
                obj.near_bad_rows_cols = obj.near_bad_rows_cols | (abs(Y-obj.pars.bad_rows(ii))<obj.pars.distance_bad_rows_cols);
            end
            
            for ii = 1:length(obj.pars.bad_columns)
                obj.near_bad_rows_cols = obj.near_bad_rows_cols | (abs(X-obj.pars.bad_columns(ii))<obj.pars.distance_bad_rows_cols);
            end
            
            k = (-obj.pars.linear_timescale/2:obj.pars.linear_timescale/2)'; 
            k = k./sqrt(sum(k.^2)); 
            
            LX = filter2(k, x - nanmean(x)); 
            LY = filter2(k, y - nanmean(y)); 
            
            ff = filter2(ones(obj.pars.linear_timescale,1)./obj.pars.linear_timescale, (f-nanmean(f))./std(f)); % filter the normalized flux with a smoothing window
            
            obj.linear_motion = sqrt(LX.^2 + LY.^2).*ff; % linear motion is set to be proportional to the filtered flux
            
            obj.background_intensity = b; % just the background level
           
            obj.bad_pixels = p;
            
            if obj.pars.use_repeating_columns % we only trigger this if we have to, it takes a while to calculate
                margins = floor((size(obj.extended_flux,1) - size(obj.search_flux,1))/2); % the margins on either edge of the search region
                margins = margins- obj.pars.dilate_region; % make the margins smaller to allow dilation from outside the search region
                obj.repeating_columns = img.find_repeating_columns(obj.cutouts, 'margins', margins); % additional arguments may change the fraction of a column that is tested (for cutouts 50% is fine)
            end
            
            try 
                obj.flux_corr = obj.calculateFluxCorrSampling(obj.extended_detrend, [25 50], 1000); 
            catch ME
                 if strcmp(ME.identifier, 'MATLAB:nomem')

                    util.text.date_printf('Out of memory error... pasuing for 10s and trying again');
                    
                    pause(10); 
                    
                    obj.flux_corr = obj.calculateFluxCorrSampling(obj.extended_detrend, [25 50], 1000);
                    
                else
                    rethrow(ME);
                end
                    
            end
            
            aux = zeros(size(f,1),size(f,2),length(obj.pars.corr_types), 'like', f); 
%             aux(:,:,obj.corr_indices.a) = a;
            aux(:,:,obj.corr_indices.b) = b;
            aux(:,:,obj.corr_indices.x) = x;
            aux(:,:,obj.corr_indices.y) = y;
            aux(:,:,obj.corr_indices.r) = r;
            aux(:,:,obj.corr_indices.w) = w;
            
            aux = fillmissing(aux, 'spline'); 
            aux(isnan(aux)) = 0; % must make sure to get rid of all NaNs
            
            obj.correlations = zeros(size(aux,1), size(aux,2), size(aux,3), length(obj.pars.corr_timescales), 'like', aux); 
            
            for ii = 1:length(obj.pars.corr_timescales)
                % correlation of (mean reduced) x and y is x.*y./sqrt(sum(x.^2)*sum(y.^2))
                % the correction term sqrt(N) is to compensate because the 
                % shorter time scales always have higher correlations. 
                
                norm = sqrt(obj.pars.corr_timescales(ii));
                if isa(f, 'single') && isa(aux, 'single')
                    norm = single(norm);
                end
                
                try
%                     if nnz(isnan(f))
%                         disp('The flux has NaN values...'); 
%                     end
                    
%                     for jj = 1:size(aux,3)
%                         if nnz(isnan(aux(:,:,jj)))
%                             fprintf('The %s have NaN values...\n', obj.aux_names{jj});
%                         end
%                     end
                    
                    obj.correlations(:,:,:,ii) = util.series.correlation(obj.extended_detrend, aux, obj.pars.corr_timescales(ii)).*norm; 
                    
                catch ME
                    disp('here'); 
                    rethrow(ME); 
                end
                
            end
            
            obj.storeCutValues;
            
            obj.fillHistograms; 
            
            if obj.pars.use_auto_count_hours
                obj.hours.input(obj); % count the star hours 
            end
            
        end
        
        function storeCutValues(obj)
        
            obj.cut_flag_matrix = false(obj.num_frames, obj.num_stars, obj.num_cuts); % preallocate the flag matrix with zeros
            obj.cut_values_matrix = zeros(obj.num_frames, obj.num_stars, obj.num_cuts, 'single'); % preallocate the cut values matrix with zeros
            
            all_stars = true(obj.num_frames, obj.num_stars); % use this to expand vector results to apply to all stars
            
            if obj.pars.use_delta_t
                obj.cut_values_matrix(:,:,obj.cut_indices.delta_t) = obj.delta_t.* all_stars; % any large time delay/jump is considered a region with bad timestamps
            end
            
            if obj.pars.use_shakes
                obj.cut_values_matrix(:,:,obj.cut_indices.shakes) = obj.shakes.*all_stars; 
            end
            
            if obj.pars.use_defocus
                obj.cut_values_matrix(:,:,obj.cut_indices.defocus) = obj.defocus.*all_stars; 
            end
            
            if obj.pars.use_fwhm && ~isempty(obj.fwhm)
                obj.cut_values_matrix(:,:,obj.cut_indices.fwhm) = obj.fwhm.*all_stars; 
            end
            
            if obj.pars.use_slope
                obj.cut_values_matrix(:,:,obj.cut_indices.slope) = obj.slope.*all_stars;
            end
            
            if obj.pars.use_near_bad_rows_cols
                obj.cut_values_matrix(:,:,obj.cut_indices.near_bad_rows_cols) = obj.near_bad_rows_cols; 
            end
            
            if obj.pars.use_offset_size
                obj.cut_values_matrix(:,:,obj.cut_indices.offset_size) = obj.offset_size;
            end

            if obj.pars.use_linear_motion
                obj.cut_values_matrix(:,:,obj.cut_indices.linear_motion) = obj.linear_motion; 
            end
            
            if obj.pars.use_background_intensity
                obj.cut_values_matrix(:,:,obj.cut_indices.background_intensity) = obj.background_intensity; 
            end
            
            if obj.pars.use_nan_flux
                obj.cut_values_matrix(:,:,obj.cut_indices.nan_flux) = obj.nan_flux;
            end
            
            if obj.pars.use_nan_offsets
                obj.cut_values_matrix(:,:,obj.cut_indices.nan_offsets) = obj.nan_offsets;
            end
            
            if obj.pars.use_photo_flag
                obj.cut_values_matrix(:,:,obj.cut_indices.photo_flag) = obj.photo_flag; 
            end
            
            if obj.pars.use_bad_pixels
                obj.cut_values_matrix(:,:,obj.cut_indices.bad_pixels) = obj.bad_pixels;
            end
            
            if obj.pars.use_repeating_columns
                obj.cut_values_matrix(:,:,obj.cut_indices.repeating_columns) = obj.repeating_columns;
            end
            
            if obj.pars.use_flux_corr
                obj.cut_values_matrix(:,:,obj.cut_indices.flux_corr_25) = obj.flux_corr(:,:,1); 
                obj.cut_values_matrix(:,:,obj.cut_indices.flux_corr_50) = obj.flux_corr(:,:,2); 
            end
            
            if obj.pars.use_correlations
                
                for ii = 1:length(obj.pars.corr_types)

                    for jj = 1:length(obj.pars.corr_timescales)

                        name = sprintf('corr_%s_%d', obj.pars.corr_types{ii}, obj.pars.corr_timescales(jj)); % e.g., corr_a_25 or corr_x_100 

%                         obj.cut_flag_matrix(:,:,obj.cut_indices.(name)) = obj.correlations(:,:,ii,jj) > obj.thresh_correlation; 
                        obj.cut_values_matrix(:,:,obj.cut_indices.(name)) = obj.correlations(:,:,ii,jj);
                        
                    end

                end
                
            end
            
            for ii = 1:length(obj.cut_names)
                
                if obj.cut_two_sided(ii)
                    obj.cut_flag_matrix(:,:,ii) = abs(obj.cut_values_matrix(:,:,ii)) >= obj.cut_thresholds(ii); 
                else
                    obj.cut_flag_matrix(:,:,ii) = obj.cut_values_matrix(:,:,ii) >= obj.cut_thresholds(ii); 
                end
                
            end
            
            if obj.pars.use_dilate
                obj.cut_flag_matrix = imdilate(obj.cut_flag_matrix, ones(1+2.*obj.pars.dilate_region,1)); 
            end
            
        end
            
        function val = calculateDefocus(obj, cutouts, x, y, f)
            
            if nargin<2 || isempty(cutouts)
                cutouts = obj.cutouts;
            end
            
            if nargin<3 || isempty(x)
                x = obj.extended_aux(:,:, obj.aux_indices.offsets_x); 
            end
            
            if nargin<4 || isempty(y)
                y = obj.extended_aux(:,:, obj.aux_indices.offsets_y); 
            end
            
            if nargin<5 || isempty(f)
                f = obj.extended_flux; 
            end
            
            N = min(obj.pars.num_stars_defocus,size(cutouts,4)); % use only 100 stars, or less if there are not enough stars            
            cutouts = cutouts(:,:,:,1:N);
            x = x(:,1:N); 
            y = y(:,1:N); 
            f = f(:,1:N); 
            
            cutouts = util.img.align(cutouts, -x, -y); % re-center all cutouts according to the x,y offsets
            
            C = nansum(cutouts,3); % sum up all frames for each star
            
%             w = util.img.fwhm(C,'method', 'filters', 'generaized', 5, 'step_size', 0.25)./2.355; % use generalized gaussian to find the width
            w = util.img.fwhm(C,'method', 'filters', 'defocus', 1, 'step_size', 0.25)./2.355; % use defocus annulus to find the width
            
            F = nanmean(f,1); % the mean flux is used as a weight
            
            val = util.vec.weighted_average(w, sqrt(abs(F)), 2); 
            
        end
        
        function val = calculateFluxCorrFastMoreMemory(obj, flux, widths) 
            
            if nargin<3 || isempty(widths)
                widths = 50; 
            end
            
            if ndims(flux)>2
                error('Can only handle 2D fluxes'); 
            end
            
            flux = flux - nanmean(flux); % make sure the flux is normalized to zero
            
            f1 = flux; % 2D matrix 
            f2 = permute(flux, [1,3,2]); % turn the star index into 3rd dimension
            
            F = f1.*f2; % big, 3D matrix
            F1 = f1.^2; % this should still be 2D
            F2 = f2.^2; % this has only dim1 and dim3 (dim2 is singleton!)
            
            val = zeros([size(flux), length(widths)], 'like', flux); % preallocate enough room for multiple widths
            
            for ii = 1:length(widths)
                
                w = widths(ii); 

                numer = movmean(F, w, 'omitnan'); 
                sum1 = movmean(F1, w, 'omitnan'); 
                sum2 = movmean(F2, w, 'omitnan'); 
                
                bad_idx = sum1<0 | sum2<0;
                
                try 
                    FC = numer./sqrt(sum1.*sum2).*sqrt(w); % flux correlations
                catch ME
                    
                    if strcmp(ME.identifier, 'MATLAB:nomem')
                        
                        util.text.date_printf('Out of memory error... Using a loop instead!');
                        
                        FC = zeros(size(numer), 'like', numer); 
                        
                        for jj = 1:size(numer,3)
                            FC(:,:,jj) = numer(:,:,jj)./sqrt(sum1.*sum2(:,:,jj)); 
                        end
                        
                    else
                        rethrow(ME);
                    end
                    
                end
                
                
                FC(bad_idx) = 0; 
                
                for jj = 1:size(flux,2)
                    FC(:,jj,jj) = 0; % remove auto-correlations (should be all equal to 1 anyway
                end
                
                corr_value = zeros(size(FC,1), size(FC,2), 'like', FC); 
                
                for jj = 1:size(FC,3)
                    corr_value(:,jj) = prctile(FC(:,jj,:), 95, 3); % should we include anti-correlations (by setting abs() first)?
                end
                
                val(:,:,ii) = corr_value; 

            end
            
        end
        
        function val = calculateFluxCorr(obj, flux, widths)
            %
            % flux  - a 2D matrix of (time x star)
            % width - vector of time scales (e.g., [10 25 50])
            
            if nargin<3 || isempty(widths)
                widths = 50; 
            end
            
            if ndims(flux)>2
                error('Can only handle 2D fluxes'); 
            end
            
            flux = flux - nanmean(flux); % make sure the flux is normalized to zero
            
            F2 = flux.^2; % this should still be 2D
            
            val = zeros([size(flux), length(widths)], 'like', flux); % preallocate enough room for multiple widths
            
            for ii = 1:size(flux,2)
            
                f = flux(:,ii).*flux; % multiply one column in "flux" with all others

                for jj = 1:length(widths)

                    w = widths(jj); 

                    numer = movmean(f, w, 'omitnan'); 
                    S = movmean(F2, w, 'omitnan'); % sum of the flux^2 
                    
                    bad_idx = S<0;

                    FC = numer./sqrt(S(:,ii).*S).*sqrt(w); % flux correlations
                    
                    FC(bad_idx) = 0; 

                    FC(:,ii) = 0; % ignore self correlations
                    
                    % should we include anti-correlations (by setting abs() first)?
                    corr_value = prctile(FC, 95, 2); % get the points with the 95th percentile correlation 

                    val(:,ii,jj) = corr_value; 

                end
                
            end
            
        end
        
        function val = calculateFluxCorrSampling(obj, flux, widths, num_stars)
            % Calculate the cross correlation between fluxes of different
            % stars over intervals of size "widths" (can be a vector of 
            % width values). 
            % Inputs:
            %   *flux: a 2D matrix of (time x star)
            %   *widths: a scalar or vector of time scales (e.g., [25 50])
            %   *num_stars: maximum number of stars to test each flux
            %               against. If smaller than number of fluxes it is 
            %               ignored. Default is all fluxes. 
            %
            % Use formula for rho_ij from: https://www.sciencedirect.com/topics/earth-and-planetary-sciences/cross-correlation
            
            if nargin<3 || isempty(widths)
                widths = 50; 
            end
            
            if nargin<4 || isempty(num_stars)
                num_stars = size(flux,2); 
            end
            
            if num_stars>size(flux,2)
                num_stars = size(flux,2); 
            end
            
            if ndims(flux)>2
                error('Can only handle 2D fluxes'); 
            end
            
            val = zeros([size(flux), length(widths)], 'like', flux); % preallocate enough room for multiple widths
            
            for ii = 1:length(widths)

                w = widths(ii); 
                midpoints = floor(w/2):w:size(flux,1); % indices of mid points of intervals for xcorr calculation
                best_corr = zeros(length(midpoints), size(flux,2)); 
                
                for jj = 1:length(midpoints)
                    
                    idx = midpoints(jj)-floor(w/2)+1:midpoints(jj)+ceil(w/2); % which indices to pull out from the flux matrix
                    
                    f = flux(idx, :); 
                    
                    f = f - nanmean(f,1); 
                    fp = permute(f(:,1:num_stars), [1,3,2]); % turn 2nd dim into 3rd, and truncate the number of fluxes
                    
                    f2 = f.^2; 
                    f2p = permute(f2(:,1:num_stars), [1,3,2]); % turn 2nd dim into 3rd, and truncate the number of fluxes
                    
                    FC = nansum(f.*fp, 1)./sqrt(nansum(f2,1).*nansum(f2p,1)).*sqrt(w); % full 3D matrix (3rd dim is shorter than 2nd if using num_stars)
                                       
                    for kk = 1:size(flux,2)
                        FC(1,kk,kk) = 0; % remove auto-correlations (should be all equal to 1 anyway)
                    end
                    
                    % should we include anti-correlations (by setting abs() first)?
                    best_corr(jj,:) = prctile(FC, 95, 3); % get the points with the 95th percentile correlation from all stars in the 3rd dim, for each sampling point
                                        
                end % for jj (time jumps)
                
                InterpObj = griddedInterpolant({midpoints, 1:size(flux,2)}, best_corr);
                val(:,:,ii) = InterpObj({1:size(flux,1), 1:size(flux,2)});  
                
            end % for ii (over widths)
            
        end
        
        function val = calculateFluxCorrSamplingNotWorking(obj, flux, widths, num_stars)
            %
            % flux  - a 2D matrix of (time x star)
            % width - vector of time scales (e.g., [10 25 50])
            % num_stars - maximum number of stars to test each flux
            % against. If smaller than number of fluxes it is ignored. 
            % Default is all fluxes. 
            
            if nargin<3 || isempty(widths)
                widths = 50; 
            end
            
            if nargin<4 || isempty(num_stars)
                num_stars = size(flux,2); 
            end
            
            if num_stars>size(flux,2)
                num_stars = size(flux,2); 
            end
            
            if ndims(flux)>2
                error('Can only handle 2D fluxes'); 
            end
            
            flux = flux - nanmean(flux); % make sure the flux is normalized to zero
            
            f2 = flux.^2; % this should still be 2D
            
            val = zeros([size(flux), length(widths)], 'like', flux); % preallocate enough room for multiple widths
            
            for jj = 1:length(widths)

                w = widths(jj); 
  
%                 f2_bin = util.series.binning(f2(:,1:num_stars), w); 
                % replace the above line with quick binning:
                f2_bin = reshape(f2, [w, floor(size(f2,1)/w), size(f2,2)]);
                f2_bin = nanmean(f2_bin, 1); 
                f2_bin = permute(f2_bin, [2,3,1]); 
                
                for ii = 1:size(flux,2)
            
                    ff = flux(:,ii).*flux(:,1:num_stars); % multiply one column in "flux" with all others

                    idx = ceil(w/2):w:size(ff,1); % places in the ligtcurve where you want to calculate the correlations/moving mean 

%                     ff_bin = util.series.binning(ff, w); 
                    % replace the above line with quick binning:
                    ff_bin = reshape(ff, [w, floor(size(ff,1)/w), size(ff,2)]);
                    ff_bin = nanmean(ff_bin, 1); 
                    ff_bin = permute(ff_bin, [2,3,1]); 
                    
                    FC = ff_bin./sqrt(f2_bin(:,ii).*f2_bin(:,1:num_stars)).*sqrt(w); % flux correlations
                    
                    FC(f2_bin(:,1:num_stars)<0) = 0; % remove places where the denominator is negative

                    FC(:,ii) = 0; % ignore self correlations
                    
                    % should we include anti-correlations (by setting abs() first)?
                    corr_value = prctile(FC, 95, 2); % get the points with the 95th percentile correlation 
                    InterpObj = griddedInterpolant(idx, corr_value);
                    val(:,ii,jj) = InterpObj(1:size(ff,1));  

                end
                
            end
            
        end
        
        function ingestStore(obj, store) % parse the data from the store into similar containers in this object
            
            obj.background_flux = store.background_flux;
            obj.background_detrend = store.background_detrend; 
            obj.background_aux = store.background_aux;
            obj.background_timestamps = store.background_timestamps;

            obj.extended_flux = store.extended_flux;
            obj.extended_detrend = store.extended_detrend; 
            obj.extended_aux = store.extended_aux;
            obj.extended_timestamps = store.extended_timestamps;
            obj.cutouts = store.cutouts;
            
            obj.search_start_idx = store.search_start_idx;
            obj.search_end_idx = store.search_end_idx;
            obj.search_flux = store.search_flux;
            obj.search_detrend = store.search_detrend; 
            obj.search_aux = store.search_aux;
            obj.search_timestamps = store.search_timestamps;
            obj.search_juldates = store.search_juldates;
            
            obj.aux_names = store.aux_names;
            obj.aux_indices = store.aux_indices;
            
        end
        
        function makeHistograms(obj) %  make histograms for each cut type, keeping the number of times a cut had such-and-such value for each star
            
            obj.histograms = zeros(obj.pars.num_hist_edges, obj.num_stars, length(obj.cut_names), 'single'); 
            
            for ii = 1:length(obj.cut_names) 
                
                name = obj.cut_names{ii};
                
                thresh = obj.cut_thresholds(ii)*2; 
                
                if obj.cut_two_sided(ii)
                    obj.hist_edges(obj.cut_indices.(name),:) = linspace(-thresh, thresh, obj.pars.num_hist_edges+1); 
                else
                    obj.hist_edges(obj.cut_indices.(name),:) = linspace(0, thresh, obj.pars.num_hist_edges+1); 
                end
            end
                
        end
        
        function fillHistograms(obj) % put data into the cut value histograms
            
            if isempty(obj.histograms)
                obj.makeHistograms;
            end
            
            idx = obj.search_start_idx:obj.search_end_idx; % indices inside the search region
            
            star_edges = 1:obj.num_stars+1; 
            star_edges_rep = repmat(star_edges, [length(idx), 1]);
            star_edges_rep = star_edges_rep(:,1:end-1); 
            
            for ii = 1:length(obj.cut_names) 
                
                values = obj.cut_values_matrix(:,:,ii); 
                
                values = values(idx,:);
                obj.histograms(:,:,ii) = obj.histograms(:,:,ii) + histcounts2(values, star_edges_rep, obj.hist_edges(ii,:), star_edges-0.5); 
                
            end
            
        end
        
        function val = checkBatchGood(obj)
            
            val = 1; 
            
            if ~isempty(obj.extended_aux) 
                B = obj.extended_aux(obj.search_start_idx:obj.search_end_idx,:,obj.aux_indices.backgrounds);
                if B>5*obj.pars.thresh_background_intensity % this is 5 times higher than the regular threshold (wide margin)
                    val = 0;
                end
            end
            
            if obj.pars.use_defocus && ~isempty(obj.defocus) && obj.defocus>2*obj.pars.thresh_defocus
                val = 0; 
            end
            
            if obj.pars.use_fwhm && ~isempty(obj.fwhm) && obj.fwhm>2*obj.pars.thresh_fwhm
                val = 0; 
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function showHistogram(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('type', 1); 
            input.input_var('sum', false);
            input.input_var('axes', [], 'axis'); 
            input.input_var('font_size', 18); 
            input.input_var('log', true, 'use_log'); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            if ischar(input.type)
                input.type = obj.cut_indices.(input.type); 
            end
            
            name = obj.cut_names{input.type}; 
                        
            if strcmp(input.axes.NextPlot, 'replace')
                reset(input.axes); 
            end

            if nnz(obj.histograms(:,:,input.type))==0
                util.plot.inner_title(strrep(sprintf('Cut "%s" is empty!', name), '_', ' '), 'ax', input.axes, 'position', 'North');                                 
            else
            
                if util.text.cs(name(1:4), 'corr')
                    thresh = obj.pars.thresh_correlation;
                elseif length(name)>=9 && util.text.cs(name(1:9), 'flux_corr')
                    thresh = obj.pars.thresh_flux_corr; 
                else
                    if isfield(obj.pars, ['thresh_' name])
                        thresh = obj.pars.(['thresh_' name]); 
                    else
                        thresh = 0.5;
                    end
                end

                if input.sum

                    values = sum(obj.histograms(:,:,input.type),2); 
                    bar(input.axes, obj.hist_edges(input.type, 1:end-1), values); 
                    mx = max(values); 
                    mx = mx.*1.1;


                    colorbar(input.axes, 'off'); 
                    
                    if input.log                        
                        input.axes.YScale = 'log'; 
                        input.axes.YLim = [1 mx]; 
                    else
                        input.axes.YLim = [0 mx]; 
                    end

                else

                    values = obj.histograms(:,:,input.type);
                    imagesc(input.axes, 'XData', obj.hist_edges(input.type, 1:end-1), 'CData', values');
                    ylabel(input.axes, 'Star index'); 
                    input.axes.XLim = [obj.hist_edges(input.type,1), obj.hist_edges(input.type,end)]; 
                    input.axes.YLim = [1 obj.num_stars]; 
                    colorbar(input.axes); 

                    mx = obj.num_stars; 
                    
                    input.axes.YScale = 'linear';
                    
                    if input.log
                        input.axes.ColorScale = 'log'; 
                    end

                end

                yyaxis(input.axes, 'right'); 

                plot(input.axes, thresh.*[1 1], [0 1e20], 'r--'); 

                if obj.cut_two_sided(input.type)
                    hold(input.axes, 'on'); 
                    plot(input.axes,  -thresh.*[1 1], [0 1e20], 'r--');     
                end

                input.axes.YTick = [];
                input.axes.YAxis(2).Color = [0 0 0];
%                 input.axes.YLim = [1 mx];
                
                hold(input.axes, 'off'); 
                
                yyaxis(input.axes, 'left'); 
            
            end
            
            xlabel(strrep(name, '_', ' '));
            
            input.axes.FontSize = input.font_size;
            
        end
        
    end    
    
end

