classdef QualityChecker < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        head@head.Header; 
        hours@trig.StarHours;
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
        background_aux; % cut out the background aux for calculating the variance outside the search region
        background_timestamps; % timestamps for the background region
        
        extended_flux; % cutout of the flux from the flux_buffer extended around the search region
        extended_aux; % cutout of the aux from the aux_buffer extended around the search region 
        extended_timestamps; % timestamps for the extended batch region
        
        search_start_idx; % starting index for the search region out of the EXTENDED BATCH!
        search_end_idx;  % end index for the search region out of the EXTENDED BATCH!
        
        search_flux; % the flux in the search region
        search_aux; % the aux in the search region
        search_timestamps;  % timestamps for the search region
        
        aux_names; % get the names (and indices below) from DataStore
        aux_indices; % struct with field=number for each of the above aux names
        
        histograms = []; % cut values accumulated over the run, with dim1 the edges defined for each cut type, dim2 the star number and dim3 is the cut type
        hist_edges = []; % edges for each cut (each row is edges for a different cut)
        
        mean_width_values; % track the mean width for all batches in this run
        mean_background_values;  % track the mean background for all batches in this run
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        num_frames;
        num_stars;
        num_cuts;
        
    end
    
    properties(Hidden=true)
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = QualityChecker(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.QualityChecker')
                if obj.debug_bit>1, fprintf('QualityChecker copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('QualityChecker constructor v%4.2f\n', obj.version); end

                obj.hours = trig.StarHours; 
        
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

            obj.pars.subtract_mean_offsets = true; % before doing any calculations on the offsets, remove the mean offsets that also affect the forced photometry centroids

            % do we want to apply all these cuts? 
            obj.pars.use_delta_t = true;
            obj.pars.use_shakes = true;
            obj.pars.use_defocus = true;
            obj.pars.use_slope = true;
            obj.pars.use_near_bad_rows_cols = true;     
            
            obj.pars.use_offset_size = true;
            obj.pars.use_linear_motion = true;
            obj.pars.use_background_intensity = true; 
            
            obj.pars.use_nan_flux = false;
            obj.pars.use_nan_offsets = false;
            obj.pars.use_photo_flag = false; 
            
            obj.pars.use_correlations = true;

            obj.pars.corr_types = {'b', 'x', 'y', 'r', 'w'}; % types of auxiliary we will use for correlations (r is derived from x and y)
            obj.pars.corr_timescales = [10, 25, 50]; % number of frames to run each correlation sum

            obj.pars.thresh_delta_t = 0.3; % events where the difference in timestamps, relative to the mean time-step are disqualified
            obj.pars.thresh_shakes = 5; % events where the mean offset r is larger than this are disqualified
            obj.pars.thresh_defocus = 2; % events with PSF width above this value are disqualified
            obj.pars.thresh_slope = 5; % events where the slope is larger than this value (in abs. value) are disqualified
            obj.pars.thresh_offset_size = 4; % events with offsets above this number are disqualified (after subtracting mean offsets)
            obj.pars.thresh_linear_motion = 2; % events showing linear motion of the centroids are disqualified
            obj.pars.thresh_background_intensity = 10; % events where the background per pixel is above this threshold are disqualified
            
            obj.pars.thresh_correlation = 4; % correlation max/min of flux (with e.g., background) with value above this disqualifies the region

            obj.pars.smoothing_slope = 50; % number of frames to average over when calculating slope
            obj.pars.distance_bad_rows_cols = 5; % how many pixels away from a bad row/column would we still disqualify a star? (on either side, inclusive)
            obj.pars.bad_columns = []; % which columns are considered bad
            obj.pars.bad_rows = []; % which rows are considered bad

            obj.pars.linear_timescale = 25; % timescale for linear_motion cut

            obj.pars.use_dilate = true; % do we want to spread out bad regions a few frames in either direction? 
            obj.pars.dilate_region = 5; % how many frames, on either side, do we flag next to each bad point? 

            obj.pars.num_hist_edges = 200; % including 100 negative and 100 positive values
            
            obj.setupSensor;
            
            obj.reset; 
            
        end
        
        function reset(obj)
            
            obj.hours.reset;
            
            obj.setupCuts;
            
            obj.hist_edges = [];
            obj.histograms = [];
            
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
            obj.background_aux = [];
            obj.background_timestamps = [];

            obj.extended_flux = [];
            obj.extended_aux = [];
            obj.extended_timestamps = [];

            obj.search_start_idx = [];
            obj.search_end_idx = [];

            obj.search_flux = [];
            obj.search_aux = [];
            obj.search_timestamps = [];

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
            
            if util.text.cs(camera, 'Balor')
                obj.pars.bad_columns = [1:20, 1960, 4085:4104]; 
                obj.pars.bad_rows = [1161 1765 3716]; 
            elseif util.text.cs(camera, 'Zyla')
                % I am not sure we have any bad rows/columns in the Zyla...
            else
                error('Unknown camera "%s". Use "Balor" or "Zyla...', camera);
            end
            
        end
        
    end
    
    methods % calculations
        
        function input(obj, store)
            
            if nargin<2 || isempty(store) || ~isa(store, 'trig.DataStore')
                error('Must supply a valid DataStore!'); 
            end
            
            obj.clear;
            
            obj.ingestStore(store); 
            
            %%%%%%%%%%%%% load the data into shorthands %%%%%%%%%%%%%%%%%%%
            
            t = obj.extended_timestamps; 
            f = obj.extended_flux; 
            F = nanmean(f,1); 
            
%             a = obj.extended_aux(:,:,obj.aux_indices.areas); 
            b = obj.extended_aux(:,:,obj.aux_indices.backgrounds); 
            x = obj.extended_aux(:,:,obj.aux_indices.offsets_x); 
            y = obj.extended_aux(:,:,obj.aux_indices.offsets_y); 
            X = obj.extended_aux(:,:,obj.aux_indices.centroids_x); 
            Y = obj.extended_aux(:,:,obj.aux_indices.centroids_y); 
            w = obj.extended_aux(:,:,obj.aux_indices.widths); 
            
            %%%%%%%%%%%%% calculate / prepare the data %%%%%%%%%%%%%%%%%%%%
            
            obj.mean_x = util.vec.weighted_average(x,F.^2,2); % maybe use square of flux instead??
            obj.mean_y = util.vec.weighted_average(y,F.^2,2); % last arg is for dimension 2 (average over stars)
            obj.defocus = util.vec.weighted_average(w,F.^2,2); % get the average PSF width (for focus tests)
            
            obj.mean_width_values = vertcat(obj.mean_width_values, obj.defocus); % keep a log of the focus for the entire run
            obj.mean_background_values = vertcat(obj.mean_background_values, nanmedian(b,2)); % keep a log of the sky background level for the entire run
            
            if obj.pars.subtract_mean_offsets
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
            
%             df = diff(nanmean(f,2)./nanmean(F)); % rate of change of the mean flux
%             obj.slope = filter2(ones(obj.smoothing_slope,1), [df(1); df]); % repeat the first df position because it needs to be the same length
            
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
            
            ff = filter2(ones(obj.pars.linear_timescale,1)./obj.pars.linear_timescale, (f-F)./std(f)); % filter the normalized flux with a smoothing window
            
            obj.linear_motion = sqrt(LX.^2 + LY.^2).*ff; % linear motion is set to be proportional to the filtered flux
            
            obj.background_intensity = b; % just the background level
            
            
            
            aux = zeros(size(f,1),size(f,2),length(obj.pars.corr_types), 'like', f); 
%             aux(:,:,obj.corr_indices.a) = a;
            aux(:,:,obj.corr_indices.b) = b;
            aux(:,:,obj.corr_indices.x) = x;
            aux(:,:,obj.corr_indices.y) = y;
            aux(:,:,obj.corr_indices.r) = r;
            aux(:,:,obj.corr_indices.w) = w;
            
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
                    obj.correlations(:,:,:,ii) = util.series.correlation(f, aux, obj.pars.corr_timescales(ii)).*norm; 
                catch ME
                    disp('here'); 
                    rethrow(ME); 
                end
            end
            
            %%%%%%%%%%%%% find all the bad frames %%%%%%%%%%%%%%%%%%%%%

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
            
            if obj.pars.use_slope
                obj.cut_values_matrix(:,:,obj.cut_indices.slope) = obj.slope.*all_stars;
            end
            
            if obj.pars.use_offset_size
                obj.cut_values_matrix(:,:,obj.cut_indices.offset_size) = obj.offset_size;
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
            
            if obj.pars.use_near_bad_rows_cols
                obj.cut_values_matrix(:,:,obj.cut_indices.near_bad_rows_cols) = obj.near_bad_rows_cols; 
            end
            
            if obj.pars.use_linear_motion
                obj.cut_values_matrix(:,:,obj.cut_indices.linear_motion) = obj.linear_motion; 
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
            
            obj.fillHistograms; 
            
            if obj.pars.use_dilate
                obj.cut_flag_matrix = imdilate(obj.cut_flag_matrix, ones(1+2.*obj.pars.dilate_region,1)); 
            end

            if obj.pars.use_auto_count_hours
                obj.hours.input(obj); % count the star hours 
            end
            
        end
        
        function ingestStore(obj, store)
            
            obj.background_flux = store.background_flux;
            obj.background_aux = store.background_aux;
            obj.background_timestamps = store.background_timestamps;

            obj.extended_flux = store.extended_flux;
            obj.extended_aux = store.extended_aux;
            obj.extended_timestamps = store.extended_timestamps;

            obj.search_start_idx = store.search_start_idx;
            obj.search_end_idx = store.search_end_idx;
            obj.search_flux = store.search_flux;
            obj.search_aux = store.search_aux;
            obj.search_timestamps = store.search_timestamps;

            obj.aux_names = store.aux_names;
            obj.aux_indices = store.aux_indices;
            
        end
        
        function fillHistograms(obj)
            
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
                obj.histograms(:,:,ii) = histcounts2(values, star_edges_rep, obj.hist_edges(ii,:), star_edges-0.5); 
                
            end
            
        end
        
        function makeHistograms(obj)
            
            obj.histograms = zeros(obj.pars.num_hist_edges, obj.num_stars, length(obj.cut_indices), 'single'); 
            
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
        
    end
    
    methods % plotting tools / GUI
        
        function showHistogram(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('type', 1); 
            input.input_var('stars', []); 
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

