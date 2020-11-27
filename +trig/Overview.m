classdef Overview < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        runtime; % number of seconds of processed data
        batch_counter; % number of batches that were processed
        
        ecl_edges; % ecliptic latitude bin edges used in the star_seconds and losses histograms
        snr_edges; % S/N bin edges used in the star_seconds and losses histograms
        size_edges; % stellar size edges used in the star_seconds and losses histograms
        
        star_seconds; % the number of useful seconds accumulated each S/N bin (this is saved after subtracting losses)
        star_seconds_with_losses; % the number of seconds accumulated without excluding anything
        losses_exclusive; % same as histogram, only counting the losses due to each cut (exclusive means times where ONLY this cut was responsible)
        losses_inclusive; % same as histogram, only counting the losses due to each cut (inclusive means times where this cut was ALSO responsible)
        losses_bad_stars; % number of seconds lost when disqualifying stars from the black list 
        
        cut_names = {}; % get these from the QualityChecker
        cut_indices = {}; % get these from the QualityChecker
        cut_thresholds = []; % get these from the QualityChecker
        cut_two_sided = []; % get these from the QualityChecker
        
        cut_histograms; % accumualted values of different cuts for each value, summed over all stars (dim1 is values, dim2 is cut type)
        cut_bin_edges; % a vector with the bin edges for each cut (same dimensions as the histograms)
        
        sim_events; % struct describing simulated events that passed/failed the threshold

        
        num_events_expected; % use simulations to estimate this
        
    end
    
    properties % switches/controls
        
        ecl_bin_width = 2; % ecliptic latitude bins (from -90 to 90)
        snr_bin_width = 1; % S/N bin width
        size_bin_width = 0.1; % stellar size bin width
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Overview(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.Overview')
                if obj.debug_bit>1, fprintf('Overview copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('Overview constructor v%4.2f\n', obj.version); end
                
                obj.reset; 
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.runtime = 0;
            obj.batch_counter = 0;

            obj.ecl_edges = [];
            obj.snr_edges = [];
            obj.size_edges = [];

            obj.star_seconds = [];
            obj.star_seconds_with_losses = [];
            obj.losses_inclusive = [];
            obj.losses_exclusive = [];
            obj.losses_bad_stars = [];

            obj.cut_names = {};
            obj.cut_indices = {}; 
            obj.cut_thresholds = [];
            obj.cut_two_sided = [];

            obj.cut_histograms = []; 
            obj.cut_bin_edges = [];

            obj.sim_events = [];

            obj.num_events_expected = 0; 
            
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations / API
        
        function input(obj, summary) % accepts a summary or an overview object...
            
            if isa(summary, 'trig.RunSummary')
                
                if isempty(summary.head) || isempty(summary.head.ephem) || isempty(summary.head.ephem.ECL_lat)
                    error('Summary object must have a valid header with an Ephemeris object / ecliptic latitude value!'); 
                end
                
                ecl = summary.head.ephem.ECL_lat; 
                snr = summary.snr_bin_edges;
                sizes = summary.size_bin_edges;
                
                inclusive = permute(summary.losses_inclusive, [1,2,4,3]); 
                exclusive = permute(summary.losses_exclusive, [1,2,4,3]); 
                
            elseif isa(summary, 'trig.Overview')
                
                ecl = summary.ecl_edges;
                snr = summary.snr_edges;
                sizes = summary.size_edges;
                
                inclusive = summary.losses_inclusive;
                exclusive = summary.losses_exclusive;
                
            else
                error('Must input a "trig.RunSummary" or "trig.Overview" object to this function! Instead got "%s"...', class(summary)); 
            end
            
            if ~isempty(obj.cut_names) % check both objects have the same cuts
                if ~isequal(obj.cut_names, summary.cut_names) || ...
                    ~isequal(obj.cut_indices, summary.cut_indices) || ...
                    ~isequal(obj.cut_bin_edges, summary.cut_bin_edges) || ...
                    ~isequal(obj.cut_thresholds, summary.cut_thresholds) || ...
                    ~isequal(obj.cut_two_sided, summary.cut_two_sided)
                    error('Cannot combine objects with different cuts!'); 
                end
            end
            
            % these must be the same for all inputs!
            obj.cut_names = summary.cut_names;
            obj.cut_indices = summary.cut_indices;
            obj.cut_bin_edges = summary.cut_bin_edges;
            obj.cut_thresholds = summary.cut_thresholds;
            obj.cut_two_sided = summary.cut_two_sided;
            
            if isempty(obj.cut_histograms)
                obj.cut_histograms = summary.cut_histograms;
            else
                obj.cut_histograms = obj.cut_histograms + summary.cut_histograms;
            end
            
            seconds = summary.star_seconds;
            seconds_with_losses = summary.star_seconds_with_losses;
            bad_stars = summary.losses_bad_stars;
            
            obj.runtime = obj.runtime + summary.runtime;
            obj.batch_counter = obj.batch_counter + summary.batch_counter;
            
            % these commands also allocate space for the histograms
            obj.add_bin_edges_snr(snr);
            obj.add_bin_edges_size(sizes); 
            obj.add_bin_edges_ecl(ecl); 
            
            % turn scalars into bin edges... 
            if isscalar(snr), snr = [snr, snr+1]; end
            if isscalar(sizes), sizes = [sizes, sizes+1]; end
            if isscalar(ecl), ecl = [ecl, ecl+1]; end
            
            for ii = 1:length(snr)-1
                idx1(ii) = find(obj.snr_edges<=snr(ii),1,'last'); 
            end
            
            for jj = 1:length(sizes)-1
                idx2(jj) = find(obj.size_edges<=sizes(jj),1,'last'); 
            end
            
            for kk = 1:length(ecl)-1
                idx3(kk) = find(obj.ecl_edges<=ecl(kk),1,'last'); 
            end
            
            for ii = 1:length(snr)-1
                
                for jj = 1:length(sizes)-1
                    
                    for kk = 1:length(ecl)-1
                        
                        obj.star_seconds(idx1(ii), idx2(jj), idx3(kk)) = obj.star_seconds(idx1(ii), idx2(jj), idx3(kk)) + seconds(ii,jj,kk);
                        obj.star_seconds_with_losses(idx1(ii), idx2(jj), idx3(kk)) = obj.star_seconds_with_losses(idx1(ii), idx2(jj), idx3(kk)) + seconds_with_losses(ii,jj,kk);
                        
                        for c = 1:length(obj.cut_names)
                            obj.losses_inclusive(idx1(ii), idx2(jj), idx3(kk), c) = obj.losses_inclusive(idx1(ii), idx2(jj), idx3(kk), c) + inclusive(ii,jj,kk,c);
                            obj.losses_exclusive(idx1(ii), idx2(jj), idx3(kk), c) = obj.losses_exclusive(idx1(ii), idx2(jj), idx3(kk), c) + exclusive(ii,jj,kk,c);
                        end
                        
                        obj.losses_bad_stars(idx1(ii), idx2(jj), idx3(kk)) = obj.losses_bad_stars(idx1(ii), idx2(jj), idx3(kk)) + bad_stars(ii,jj,kk);
                        
                    end % for kk (ecl)
                    
                end % for jj (size)
                
            end % for ii (snr)
            
            if isempty(obj.sim_events)
                obj.sim_events = summary.sim_events;
            else
                obj.sim_events = [obj.sim_events, summary.sim_events]; 
            end
            
        end
        
    end
    
    methods (Hidden=false) % internal calculations
        
        function add_bin_edges_snr(obj, new_values)
            
            if isempty(new_values)
                return;
            end
            
            mx = ceil(max(new_values(:))./obj.snr_bin_width).*obj.snr_bin_width; 
            if max(obj.snr_edges)>mx
                mx = max(obj.snr_edges); 
            end
            
            new_edges = 0:obj.snr_bin_width:mx;
            
            if ~isequal(obj.snr_edges, new_edges)
            
                % allocate space for the bigger histograms
                new_star_seconds = zeros(length(new_edges)-1, length(obj.size_edges)-1, length(obj.ecl_edges)-1, 'like', obj.star_seconds); 
                new_star_seconds_with = zeros(length(new_edges)-1, length(obj.size_edges)-1, length(obj.ecl_edges)-1, 'like', obj.star_seconds_with_losses); 
                inclusive = zeros(length(new_edges)-1, length(obj.size_edges)-1, length(obj.ecl_edges)-1, length(obj.cut_names), 'like', obj.losses_inclusive); 
                exclusive = zeros(length(new_edges)-1, length(obj.size_edges)-1, length(obj.ecl_edges)-1, length(obj.cut_names), 'like', obj.losses_exclusive); 
                bad_stars = zeros(length(new_edges)-1, length(obj.size_edges)-1, length(obj.ecl_edges)-1, 'like', obj.losses_bad_stars); 

                if ~isempty(obj.snr_edges)
                    
                    start_idx = 1;
                    end_idx = find(new_edges==obj.snr_edges(end))-1; % find the last bin index in the new matrix

                    if isempty(start_idx) || isempty(end_idx)
                        error('the new edges should include the old edges! '); 
                    end

                    new_star_seconds(start_idx:end_idx,:,:) = obj.star_seconds;
                    new_star_seconds_with(start_idx:end_idx,:,:) = obj.star_seconds_with_losses;
                    inclusive(start_idx:end_idx,:,:,:) = obj.losses_inclusive;
                    exclusive(start_idx:end_idx,:,:,:) = obj.losses_exclusive;
                    bad_stars(start_idx:end_idx,:,:) = obj.losses_bad_stars;
                    
                end
                
                obj.star_seconds = new_star_seconds;
                obj.star_seconds_with_losses = new_star_seconds_with;
                obj.losses_inclusive = inclusive;
                obj.losses_exclusive = exclusive;
                obj.losses_bad_stars = bad_stars;
                
            end
            
            obj.snr_edges = new_edges;
            
        end
        
        function add_bin_edges_size(obj, new_values)
            
            if isempty(new_values)
                return;
            end
            
            mx = ceil(max(new_values(:))./obj.size_bin_width).*obj.size_bin_width; 
            if max(obj.size_edges)>mx
                mx = max(obj.size_edges); 
            end
            
            new_edges = 0:obj.size_bin_width:mx;
            
            if ~isequal(obj.size_edges, new_edges)
            
                % allocate space for the bigger histograms
                new_star_seconds = zeros(length(obj.snr_edges)-1, length(new_edges)-1, length(obj.ecl_edges)-1, 'like', obj.star_seconds); 
                new_star_seconds_with = zeros(length(obj.snr_edges)-1, length(new_edges)-1, length(obj.ecl_edges)-1, 'like', obj.star_seconds_with_losses); 
                inclusive = zeros(length(obj.snr_edges)-1, length(new_edges)-1, length(obj.ecl_edges)-1, length(obj.cut_names), 'like', obj.losses_inclusive); 
                exclusive = zeros(length(obj.snr_edges)-1, length(new_edges)-1, length(obj.ecl_edges)-1, length(obj.cut_names), 'like', obj.losses_exclusive); 
                bad_stars = zeros(length(obj.snr_edges)-1, length(new_edges)-1, length(obj.ecl_edges)-1, 'like', obj.losses_bad_stars); 
                
                if ~isempty(obj.size_edges)
                    
                    start_idx = 1;
                    end_idx = find(new_edges==obj.size_edges(end))-1; % find the last bin index in the new matrix

                    if isempty(start_idx) || isempty(end_idx)
                        error('the new edges should include the old edges! '); 
                    end

                    new_star_seconds(:,start_idx:end_idx,:) = obj.star_seconds;
                    new_star_seconds_with(:,start_idx:end_idx,:) = obj.star_seconds_with_losses;
                    inclusive(:,start_idx:end_idx,:,:) = obj.losses_inclusive;
                    exclusive(:,start_idx:end_idx,:,:) = obj.losses_exclusive;
                    bad_stars(:,start_idx:end_idx,:) = obj.losses_bad_stars;

                end
                
                obj.star_seconds = new_star_seconds;
                obj.star_seconds_with_losses = new_star_seconds_with;
                obj.losses_inclusive = inclusive;
                obj.losses_exclusive = exclusive;
                obj.losses_bad_stars = bad_stars;
                
            end
            
            obj.size_edges = new_edges;
            
        end
        
        function add_bin_edges_ecl(obj, new_values)
            
            if isempty(new_values)
                return;
            end
            
            mx = ceil(max(new_values(:))./obj.ecl_bin_width).*obj.ecl_bin_width; 
            if max(obj.ecl_edges)>mx
                mx = max(obj.ecl_edges); 
            end
            
            mn = floor(min(new_values(:))./obj.ecl_bin_width).*obj.ecl_bin_width; 
            
            if min(obj.ecl_edges)<mn
                mn = min(obj.ecl_edges);
            end
            
            pos = 0:obj.ecl_bin_width:mx;
            neg = mn:obj.ecl_bin_width:0;
            
            new_edges = unique([neg, pos]);
            
            if ~isequal(obj.ecl_edges, new_edges)
            
                % allocate space for the bigger histograms
                new_star_seconds = zeros(length(obj.snr_edges)-1, length(obj.size_edges)-1, length(new_edges)-1, 'like', obj.star_seconds); 
                new_star_seconds_with = zeros(length(obj.snr_edges)-1, length(obj.size_edges)-1, length(new_edges)-1, 'like', obj.star_seconds_with_losses); 
                inclusive = zeros(length(obj.snr_edges)-1, length(obj.size_edges)-1, length(new_edges)-1, length(obj.cut_names), 'like', obj.losses_inclusive); 
                exclusive = zeros(length(obj.snr_edges)-1, length(obj.size_edges)-1, length(new_edges)-1, length(obj.cut_names), 'like', obj.losses_exclusive); 
                bad_stars = zeros(length(obj.snr_edges)-1, length(obj.size_edges)-1, length(new_edges)-1, 'like', obj.losses_bad_stars); 
                
                if ~isempty(obj.ecl_edges)
                    
                    start_idx = find(new_edges==obj.ecl_edges(1)); % find the first bin index in the new matrix
                    end_idx = find(new_edges==obj.ecl_edges(end))-1; % find the last bin index in the new matrix

                    if isempty(start_idx) || isempty(end_idx)
                        error('the new edges should include the old edges! '); 
                    end

                    new_star_seconds(:,:,start_idx:end_idx) = obj.star_seconds;
                    new_star_seconds_with(:,:,start_idx:end_idx) = obj.star_seconds_with_losses;
                    inclusive(:,:,start_idx:end_idx,:) = obj.losses_inclusive;
                    exclusive(:,:,start_idx:end_idx,:) = obj.losses_exclusive;
                    bad_stars(:,:,start_idx:end_idx) = obj.losses_bad_stars;

                end
                
                obj.star_seconds = new_star_seconds;
                obj.star_seconds_with_losses = new_star_seconds_with;
                obj.losses_inclusive = inclusive;
                obj.losses_exclusive = exclusive;
                obj.losses_bad_stars = bad_stars;
                
            end
            
            obj.ecl_edges = new_edges;
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function showTotalHours(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('axes', [], 'axis'); 
            input.input_var('font_size', 18); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            bar(input.axes, obj.ecl_edges(1:end-1), squeeze(util.stat.sum2(obj.star_seconds))./3600,1);

            xlabel(input.axes, 'ecliptic latitude [deg]'); 
            ylabel(input.axes, 'total star hours'); 
           
            input.axes.FontSize = input.font_size;
            
        end
        
    end    
    
end

