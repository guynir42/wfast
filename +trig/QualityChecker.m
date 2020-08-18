classdef QualityChecker < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        hours@trig.StarHours;
        
    end
    
    properties % inputs/outputs
        
        correlations; % 4D matrix: dim1 is frame index, dim2 is star index, dim3 is aux name, dim4 is timescale
        corr_indices; % generate a struct using setupIndices()
        
        delta_t; % difference between average time-step (dt) and average time-step (usually 0.04s). 
        shakes; % the offset size for each star
        defocus; % the size of the PSF
        slope; % the slope at different points in the flux
        offset_size; % the size of the offset in each star
        
        nan_flux; % places where the flux is NaN
        nan_offsets; % places where the offset_x or offset_y is nan
        photo_flag; % places flagged by photometry (mainly offsets and widths that are NaN)
        near_bad_rows_cols; % flag places where the centroids are too close to a bad row/column
        
        mean_x; % the weighted average offset for all stars
        mean_y; % the weighted average offset for all stars
        
        % rename flags to cuts? make it 2D single matrix
        % dim1 of flag is frame index, dim2 is star index 
        % flags show what kind of bad time affects each frame
        % these are logicals and are summed up (ORed) to give "bad_times"
        cut_flag_matrix; % include all flags in this big matrix, dim1 frame number, dim2 star index, dim3 which cut is used
        cut_names = {}; % automatically assemble these at startup
%         {'delta_t', 'shakes', 'defocus', 'slope', 'offset_size', ...
%             'corr_a', 'corr_b', 'corr_x', 'corr_y', 'corr_r', 'corr_w', ...
%             'size_r', 'nan_flux', 'nan_offsets', 'photometry'}; 
        cut_indices;  % generate a struct using setupIndices()
        
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
        
        histograms = [];
        hist_edges = []; 
        
        mean_width_values; % track the mean width for all batches in this run
        
        
    end
    
    properties % switches/controls
        
        subtract_mean_offsets = true; % before doing any calculations on the offsets, remove the mean offsets that also affect the forced photometry centroids
        
        % do we want to apply all these cuts? 
        use_delta_t = true;
        use_shakes = true;
        use_defocus = true;
        use_slope = true;
        use_offset_size = true;
        use_nan_flux = true;
        use_nan_offsets = false;
        use_photo_flag = false;
        use_near_bad_rows_cols = true;         
        use_correlations = true;
        
        corr_types = {'a', 'b', 'x', 'y', 'r', 'w'}; % types of auxiliary we will use for correlations (r is derived from x and y)
        corr_timescales = [25, 50, 100]; % number of frames to run each correlation sum
        
        thresh_delta_t = 0.5; % events where the difference in timestamps, relative to the mean time-step are disqualified
        thresh_shakes = 5; % events where the mean offset r is larger than this are disqualified
        thresh_defocus = 2; % events with PSF width above this value are disqualified
        thresh_slope = 0.25; % events where the slope is larger than this value (in abs. value) are disqualified
        thresh_offset_size = 5; % events with offsets above this number are disqualified (after subtracting mean offsets?)
        thresh_correlation = 5; % correlation max/min of flux (with e.g., background) with value above this disqualifies the region
        
        smoothing_slope = 50; % number of frames to average over when calculating slope
        distance_bad_rows_cols = 5; % how many pixels away from a bad row/column would we still disqualify a star? (on either side, inclusive)
        bad_columns = [];
        bad_rows = []; 
        
        use_dilate = false; % do we want to spread out bad regions a few frames in either direction? 
        dilate_region = 5; % how many frames, on either side, do we flag next to each bad point? 
        
        num_hist_edges = 200; % including 100 negative and 100 positive values
        
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
        
                obj.setupSensor; 
                
                obj.reset;
        
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.hours.reset;
            
            obj.setupIndices;
            
            obj.hist_edges = [];
            obj.histograms = [];
            
            obj.mean_width_values = [];
            
            obj.clear;
            
        end
        
        function setupIndices(obj)
            
            obj.corr_indices = struct;
            for ii = 1:length(obj.corr_types)
                obj.corr_indices.(obj.corr_types{ii}) = ii; 
            end
            
            obj.cut_names = {}; 
            if obj.use_delta_t, obj.cut_names{end+1} = 'delta_t'; end
            if obj.use_shakes, obj.cut_names{end+1} = 'shakes'; end
            if obj.use_defocus, obj.cut_names{end+1} = 'defocus'; end
            if obj.use_slope, obj.cut_names{end+1} = 'slope'; end
            if obj.use_offset_size, obj.cut_names{end+1} = 'offset_size'; end
            if obj.use_nan_flux, obj.cut_names{end+1} = 'nan_flux'; end
            if obj.use_nan_offsets, obj.cut_names{end+1} = 'nan_offsets'; end
            if obj.use_photo_flag, obj.cut_names{end+1} = 'photo_flag'; end
            if obj.use_near_bad_rows_cols, obj.cut_names{end+1} = 'near_bad_rows_cols'; end
            
            for ii = 1:length(obj.corr_types)
                
                for jj = 1:length(obj.corr_timescales)

                    name = sprintf('corr_%s_%d', obj.corr_types{ii}, obj.corr_timescales(jj)); % e.g., corr_a_25 or corr_x_100 
                    obj.cut_names{end+1} = name;
                    
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
        
    end
    
    methods % setters
        
        function setupSensor(obj, camera)
            
            if nargin<2 || isempty(camera)
                camera = 'Balor'; 
            end
            
            obj.bad_columns = [];
            obj.bad_rows = [];
            
            if util.text.cs(camera, 'Balor')
                obj.bad_columns = [1:20, 1960, 4085:4104]; 
                obj.bad_rows = [1161 1765 3716]; 
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
            
            a = obj.extended_aux(:,:,obj.aux_indices.areas); 
            b = obj.extended_aux(:,:,obj.aux_indices.backgrounds); 
            x = obj.extended_aux(:,:,obj.aux_indices.offsets_x); 
            y = obj.extended_aux(:,:,obj.aux_indices.offsets_y); 
            X = obj.extended_aux(:,:,obj.aux_indices.centroids_x); 
            Y = obj.extended_aux(:,:,obj.aux_indices.centroids_y); 
            w = obj.extended_aux(:,:,obj.aux_indices.widths); 
            
            %%%%%%%%%%%%% calculate / prepare the data %%%%%%%%%%%%%%%%%%%%
            
            obj.mean_x = util.vec.weighted_average(x,F,2); % maybe use sqrt of flux instead??
            obj.mean_y = util.vec.weighted_average(y,F,2); % last arg is for dimension 2 (average over stars)
            obj.defocus = util.vec.weighted_average(w,F,2); % get the average PSF width (for focus tests)
            
            obj.mean_width_values = vertcat(obj.mean_width_values, obj.defocus); % keep a log of the focus for the entire run
            
            if obj.subtract_mean_offsets
                x = x - obj.mean_x; 
                y = y - obj.mean_y;
            end

            w = fillmissing(w, 'spline'); 
            
            r = sqrt(x.^2 + y.^2); % size of offsets
            obj.offset_size = r; 
            obj.shakes = util.vec.weighted_average(r,F,2); % the average offset size

            mean_dt = median(diff(t)); % get the cadence (usually this is 0.04 seconds)
            dt_vector = [mean_dt; diff(t)]; % add one mean_dt in the first place to make dt_vector the same size as t
            
            obj.delta_t = (dt_vector - mean_dt)./mean_dt;
            
            df = diff(nanmean(f,2)./nanmean(F)); % rate of change of the mean flux
            obj.slope = filter2(ones(obj.smoothing_slope,1), [df(1); df]); % repeat the first df position because it needs to be the same length
            
            obj.nan_flux = isnan(f); 
            obj.nan_offsets = isnan(x) | isnan(y);  
            obj.photo_flag = obj.extended_aux(:,:,obj.aux_indices.flags); 
            
            obj.near_bad_rows_cols = false(size(X)); % initialize to all zero
            
            for ii = 1:length(obj.bad_rows)
                obj.near_bad_rows_cols = obj.near_bad_rows_cols | abs(Y-obj.bad_rows(ii))>obj.distance_bad_rows_cols;
            end
            
            for ii = 1:length(obj.bad_columns)
                obj.near_bad_rows_cols = obj.near_bad_rows_cols | abs(X-obj.bad_columns(ii))>obj.distance_bad_rows_cols;
            end
            
            aux = zeros(size(f,1),size(f,2),length(obj.corr_types), 'like', f); 
            aux(:,:,obj.corr_indices.a) = a;
            aux(:,:,obj.corr_indices.b) = b;
            aux(:,:,obj.corr_indices.x) = x;
            aux(:,:,obj.corr_indices.y) = y;
            aux(:,:,obj.corr_indices.r) = r;
            aux(:,:,obj.corr_indices.w) = w;

            for ii = 1:length(obj.corr_timescales)
                % correlation of (mean reduced) x and y is x.*y./sqrt(sum(x.^2)*sum(y.^2))
                % the correction term sqrt(N) is to compensate because the 
                % shorter time scales always have higher correlations. 
                obj.correlations(:,:,:,ii) = util.series.correlation(f, aux, obj.corr_timescales(ii))*sqrt(obj.corr_timescales(ii)); 
            end
            
            obj.fillHistograms; 
            
            %%%%%%%%%%%%% find all the bad frames %%%%%%%%%%%%%%%%%%%%%

            obj.cut_flag_matrix = false(obj.num_frames, obj.num_stars, obj.num_cuts); % preallocate the flag matrix with zeros
            
            all_stars = true(obj.num_frames, obj.num_stars); 
            idx = obj.search_start_idx:obj.search_end_idx;
            
            if obj.use_delta_t
                obj.cut_flag_matrix(:,:,obj.cut_indices.('delta_t')) = (abs(obj.delta_t) > obj.thresh_delta_t) .* all_stars; % any large time delay/jump is considered a region with bad timestamps
            end
            
            if obj.use_shakes
                obj.cut_flag_matrix(:,:,obj.cut_indices.('shakes')) = (obj.shakes > obj.thresh_shakes) .* all_stars; 
            end
            
            if obj.use_defocus
                obj.cut_flag_matrix(:,:,obj.cut_indices.('defocus')) = (obj.defocus > obj.thresh_defocus) .* all_stars; 
            end
            
            if obj.use_offset_size
                obj.cut_flag_matrix(:,:,obj.cut_indices.('offset_size')) = obj.offset_size > obj.thresh_offset_size;
            end

            if obj.use_nan_flux
                obj.cut_flag_matrix(:,:,obj.cut_indices.('nan_flux')) = obj.nan_flux;
            end
            
            if obj.use_nan_offsets
                obj.cut_flag_matrix(:,:,obj.cut_indices.('nan_offsets')) = obj.nan_offsets;
            end
            
            if obj.use_photo_flag
                obj.cut_flag_matrix(:,:,obj.cut_indices.('photometry')) = obj.photo_flag; 
            end
            
            if obj.use_correlations
                
                for ii = 1:length(obj.corr_types)

                    for jj = 1:length(obj.corr_timescales)

                        name = sprintf('corr_%s_%d', obj.corr_types{ii}, obj.corr_timescales(jj)); % e.g., corr_a_25 or corr_x_100 

                        obj.cut_flag_matrix(:,:,obj.cut_indices.(name)) = obj.correlations(:,:,ii,jj) > obj.thresh_correlation; 

                    end

                end
                
            end
            
            if obj.use_dilate
                obj.cut_flag_matrix = imdilate(obj.flags_matrix, ones(1+2.*obj.dilate_region,1)); 
            end

            obj.hours.input(obj); % count the star hours 
            
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
                
                name = obj.cut_names{ii};
                
                try
                
                if util.text.cs(name(1:4), 'corr')
                    corr_type = name(6); 
                    timescale = str2double(name(8:end)); 
                    timescale_idx = find(timescale==obj.corr_timescales); 
                    values = obj.correlations(:,:,obj.corr_indices.(corr_type), timescale_idx); 
                else
                    values = obj.(name);
                end
                    
                if size(values,2)==1
                    values = repmat(values, [1, obj.num_stars]); 
                end
                
                if ~isempty(values)
                    values = values(idx,:); 
                    obj.histograms(:,:,ii) = histcounts2(values, star_edges_rep, obj.hist_edges(ii,:), star_edges-0.5); 
                end
                
                catch ME
                    fprintf('Problem found in processing "%s"\n', name); 
                    rethrow(ME); 
                end
                
            end
            
        end
        
        function makeHistograms(obj)
            
            obj.histograms = zeros(obj.num_hist_edges, obj.num_stars, length(obj.cut_indices)); 
            
            for ii = 1:length(obj.cut_names) 
                
                name = obj.cut_names{ii};
                
                if util.text.cs(name(1:4), 'corr')
                    thresh = obj.thresh_correlation*2;
                elseif ~isprop(obj, ['thresh_' name])
                    thresh = 1; 
                else
                    thresh = obj.(['thresh_' name])*2; 
                end
                
                obj.hist_edges(obj.cut_indices.(name),:) = linspace(-thresh, thresh, obj.num_hist_edges+1); 
            
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

            if util.text.cs(name(1:4), 'corr')
                thresh = obj.thresh_correlation;
            else
                if isprop(obj, ['thresh_' name])
                    thresh = obj.(['thresh_' name]); 
                else
                    thresh = 0.5;
                end
            end

            if input.sum
                values = sum(obj.histograms(:,:,input.type),2); 
                bar(input.axes, obj.hist_edges(input.type, 1:end-1), values); 
                mx = max(values); 
                mx = mx.*1.1;
                
                input.axes.YLim = [0 mx]; 
                
                if input.log
                    input.axes.YScale = 'log'; 
                end
                
            else
                
                values = obj.histograms(:,:,input.type);
                imagesc(input.axes, 'XData', obj.hist_edges(input.type, 1:end-1), 'CData', values');
                ylabel(input.axes, 'Star index'); 
                input.axes.XLim = [obj.hist_edges(input.type,1), obj.hist_edges(input.type,end)]; 
                input.axes.YLim = [1 obj.num_stars]; 
                colorbar(input.axes); 
                
                mx = obj.num_stars; 
                
                if input.log
                    input.axes.ColorScale = 'log'; 
                end
                
            end
            
            hold(input.axes, 'on'); 
            plot(input.axes, -thresh.*[1 1], [0 mx], 'r--'); 
            plot(input.axes,  thresh.*[1 1], [0 mx], 'r--'); 
            hold(input.axes, 'off'); 
            

            xlabel(strrep(name, '_', ' '));
            
            input.axes.FontSize = input.font_size;
            
        end
        
    end    
    
end

