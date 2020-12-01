classdef Overview < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        runtime; % number of seconds of processed data
        batch_counter; % number of batches that were processed
        
        ecl_edges; % ecliptic latitude bin edges used in the star_seconds and losses histograms
        vel_edges; % transverse velocity bin edges used in the star_seconds and losses histograms
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
        vel_bin_width = 1; % transverse velocity bin width
        vel_max = 31; % transverse velocity bin maximum
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
            obj.vel_edges = [];
            obj.snr_edges = [];
            obj.size_edges = [];

            obj.star_seconds = single([]);
            obj.star_seconds_with_losses = single([]);
            obj.losses_inclusive = single([]);
            obj.losses_exclusive = single([]);
            obj.losses_bad_stars = single([]);

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
        
        function val = getSeconds(obj, varargin) % sum the star_seconds histograms with some constrainsts
            
            import util.text.cs;
            
            if isempty(obj.star_seconds)
                val = 0;
                return;
            end
            
            input = util.text.InputVars;
            input.input_var('source', 'star_seconds'); 
            input.input_var('cut', []); 
            input.input_var('snrs', [], 'signal to noise ratio'); 
            input.input_var('sizes', [], 'stellar sizes'); 
            input.input_var('ecl', [], 'ecliptic latitudes'); 
            input.input_var('vel', [], 'velocity', 'transverse velocity', 'shadow velocity'); 
            input.scan_vars(varargin{:}); 
            
            %%% get the source matrix (star_hours etc). %%%
            
            if cs(input.source, 'star_seconds', 'seconds')
                M = obj.star_seconds;
            elseif cs(input.source, 'star_seconds_with_losses', 'with_losses')
                M = obj.star_seconds_with_losses;
            elseif cs(input.source, 'inclusive')
                M = obj.losses_inclusive;
            elseif cs(input.source, 'exclusive')
                M = obj.losses_exclusive;
            elseif cs(input.source, 'bad_stars')
                M = obj.losses_bad_stars; 
            else
                error('Unknown "source" option "%s". Use "star_seconds" or "inclusive" etc. ', input.source); 
            end
            
            if ndims(M)>4
                if isempty(input.cut)
                    M = nansum(M,5); % take all cuts together! 
                else
                    idx = obj.cut_indices.(input.cut); % choose one of the cuts
                    M = M(:,:,:,:,idx); 
                end
            end
            
            %%% apply cuts on parameters %%%
            
            idx_snr = obj.findIndex(input.snrs, 'snr'); 
            idx_size = obj.findIndex(input.sizes, 'size'); 
            idx_ecl = obj.findIndex(input.ecl, 'ecl'); 
            idx_vel = obj.findIndex(input.vel, 'vel');
            
            M2 = M(idx_snr, idx_size, idx_ecl, idx_vel); 
            
            val = nansum(M2(:)); 
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations / API
        
        function input(obj, summary) % accepts a summary or an overview object...
            
            %%% ingest either Overview or RunSummary object %%%
            
            if isa(summary, 'trig.RunSummary')
                
                if isempty(summary.head) || isempty(summary.head.ephem) || isempty(summary.head.ephem.ECL_lat)
                    error('Summary object must have a valid header with an Ephemeris object / ecliptic latitude value!'); 
                end
                
                ecl = summary.head.ephem.ECL_lat; 
                vel = sqrt(sum(summary.head.ephem.getShadowVelocity.^2));
                snr = summary.snr_bin_edges;
                sizes = summary.size_bin_edges;
                
                inclusive = permute(summary.losses_inclusive, [1,2,5,4,3]); 
                exclusive = permute(summary.losses_exclusive, [1,2,5,4,3]); 
                
            elseif isa(summary, 'trig.Overview')
                
                ecl = summary.ecl_edges;
                vel = summary.vel_edges; 
                snr = summary.snr_edges;
                sizes = summary.size_edges;
                
                inclusive = summary.losses_inclusive;
                exclusive = summary.losses_exclusive;
                
            else
                error('Must input a "trig.RunSummary" or "trig.Overview" object to this function! Instead got "%s"...', class(summary)); 
            end
            
            %%% get the cut values and histogram %%%
            
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
            
            %%% get the star second histograms %%%
            
            seconds = summary.star_seconds;
            seconds_with_losses = summary.star_seconds_with_losses;
            bad_stars = summary.losses_bad_stars;
            
            obj.runtime = obj.runtime + summary.runtime;
            obj.batch_counter = obj.batch_counter + summary.batch_counter;
            
            %%% update the bin edges and allocate space for the histograms %%%
            
            obj.add_bin_edges_snr(snr);
            obj.add_bin_edges_size(sizes); 
            obj.add_bin_edges_ecl(ecl);
            obj.add_bin_edges_vel(vel); 
            
            %%% add the new data into the histograms %%%
            
            % turn scalars into bin edges... 
            if isscalar(snr), snr = [snr, snr+1]; end
            if isscalar(sizes), sizes = [sizes, sizes+1]; end
            if isscalar(ecl), ecl = [ecl, ecl+1]; end
            if isscalar(vel), vel = [vel, vel+1]; end
            
            % precalculate the bin edge indices
            for ii = 1:length(snr)-1
                idx1(ii) = find(obj.snr_edges<=snr(ii),1,'last'); 
            end
            
            for jj = 1:length(sizes)-1
                idx2(jj) = find(obj.size_edges<=sizes(jj),1,'last'); 
            end
            
            for kk = 1:length(ecl)-1
                idx3(kk) = find(obj.ecl_edges<=ecl(kk),1,'last'); 
            end
            
            for mm = 1:length(vel)-1
                idx4(mm) = find(obj.vel_edges<=vel(kk),1,'last'); 
            end
            
            % now do the mega-loop binning
            for ii = 1:length(snr)-1
                
                for jj = 1:length(sizes)-1
                    
                    for kk = 1:length(ecl)-1
                        
                        for mm = 1:length(vel)-1

                            obj.star_seconds(idx1(ii), idx2(jj), idx3(kk), idx4(mm)) = obj.star_seconds(idx1(ii), idx2(jj), idx3(kk), idx4(mm)) + seconds(ii,jj,kk,mm);
                            obj.star_seconds_with_losses(idx1(ii), idx2(jj), idx3(kk), idx4(mm)) = obj.star_seconds_with_losses(idx1(ii), idx2(jj), idx3(kk), idx4(mm)) + seconds_with_losses(ii,jj,kk,mm);

                            for c = 1:length(obj.cut_names)
                                obj.losses_inclusive(idx1(ii), idx2(jj), idx3(kk), idx4(mm), c) = obj.losses_inclusive(idx1(ii), idx2(jj), idx3(kk), idx4(mm), c) + inclusive(ii,jj,kk,mm,c);
                                obj.losses_exclusive(idx1(ii), idx2(jj), idx3(kk), idx4(mm), c) = obj.losses_exclusive(idx1(ii), idx2(jj), idx3(kk), idx4(mm), c) + exclusive(ii,jj,kk,mm,c);
                            end

                            obj.losses_bad_stars(idx1(ii), idx2(jj), idx3(kk), idx4(mm)) = obj.losses_bad_stars(idx1(ii), idx2(jj), idx3(kk), idx4(mm)) + bad_stars(ii,jj,kk,mm);

                        end % for mm (vel)
                        
                    end % for kk (ecl)
                    
                end % for jj (size)
                
            end % for ii (snr)
            
            %%% accumulate the simulated events as well
            
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
                new_star_seconds = zeros(length(new_edges)-1, length(obj.size_edges)-1, length(obj.ecl_edges)-1, length(obj.vel_edges)-1, 'like', obj.star_seconds); 
                new_star_seconds_with = zeros(length(new_edges)-1, length(obj.size_edges)-1, length(obj.ecl_edges)-1, length(obj.vel_edges)-1, 'like', obj.star_seconds_with_losses); 
                inclusive = zeros(length(new_edges)-1, length(obj.size_edges)-1, length(obj.ecl_edges)-1, length(obj.vel_edges)-1, length(obj.cut_names), 'like', obj.losses_inclusive); 
                exclusive = zeros(length(new_edges)-1, length(obj.size_edges)-1, length(obj.ecl_edges)-1, length(obj.vel_edges)-1, length(obj.cut_names), 'like', obj.losses_exclusive); 
                bad_stars = zeros(length(new_edges)-1, length(obj.size_edges)-1, length(obj.ecl_edges)-1, length(obj.vel_edges)-1, 'like', obj.losses_bad_stars); 

                if ~isempty(obj.snr_edges)
                    
                    start_idx = 1;
                    end_idx = find(new_edges==obj.snr_edges(end))-1; % find the last bin index in the new matrix

                    if isempty(start_idx) || isempty(end_idx)
                        error('the new edges should include the old edges! '); 
                    end

                    new_star_seconds(start_idx:end_idx,:,:,:) = obj.star_seconds;
                    new_star_seconds_with(start_idx:end_idx,:,:,:) = obj.star_seconds_with_losses;
                    inclusive(start_idx:end_idx,:,:,:,:) = obj.losses_inclusive;
                    exclusive(start_idx:end_idx,:,:,:,:) = obj.losses_exclusive;
                    bad_stars(start_idx:end_idx,:,:,:) = obj.losses_bad_stars;
                    
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
                new_star_seconds = zeros(length(obj.snr_edges)-1, length(new_edges)-1, length(obj.ecl_edges)-1, length(obj.vel_edges)-1, 'like', obj.star_seconds); 
                new_star_seconds_with = zeros(length(obj.snr_edges)-1, length(new_edges)-1, length(obj.ecl_edges)-1, length(obj.vel_edges)-1, 'like', obj.star_seconds_with_losses); 
                inclusive = zeros(length(obj.snr_edges)-1, length(new_edges)-1, length(obj.ecl_edges)-1, length(obj.vel_edges)-1, length(obj.cut_names), 'like', obj.losses_inclusive); 
                exclusive = zeros(length(obj.snr_edges)-1, length(new_edges)-1, length(obj.ecl_edges)-1, length(obj.vel_edges)-1, length(obj.cut_names), 'like', obj.losses_exclusive); 
                bad_stars = zeros(length(obj.snr_edges)-1, length(new_edges)-1, length(obj.ecl_edges)-1, length(obj.vel_edges)-1, 'like', obj.losses_bad_stars); 
                
                if ~isempty(obj.size_edges)
                    
                    start_idx = 1;
                    end_idx = find(new_edges==obj.size_edges(end))-1; % find the last bin index in the new matrix

                    if isempty(start_idx) || isempty(end_idx)
                        error('the new edges should include the old edges! '); 
                    end

                    new_star_seconds(:,start_idx:end_idx,:,:) = obj.star_seconds;
                    new_star_seconds_with(:,start_idx:end_idx,:,:) = obj.star_seconds_with_losses;
                    inclusive(:,start_idx:end_idx,:,:,:) = obj.losses_inclusive;
                    exclusive(:,start_idx:end_idx,:,:,:) = obj.losses_exclusive;
                    bad_stars(:,start_idx:end_idx,:,:) = obj.losses_bad_stars;

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
                new_star_seconds = zeros(length(obj.snr_edges)-1, length(obj.size_edges)-1, length(new_edges)-1, length(obj.vel_edges)-1, 'like', obj.star_seconds); 
                new_star_seconds_with = zeros(length(obj.snr_edges)-1, length(obj.size_edges)-1, length(new_edges)-1, length(obj.vel_edges)-1, 'like', obj.star_seconds_with_losses); 
                inclusive = zeros(length(obj.snr_edges)-1, length(obj.size_edges)-1, length(new_edges)-1, length(obj.vel_edges)-1, length(obj.cut_names), 'like', obj.losses_inclusive); 
                exclusive = zeros(length(obj.snr_edges)-1, length(obj.size_edges)-1, length(new_edges)-1, length(obj.vel_edges)-1, length(obj.cut_names), 'like', obj.losses_exclusive); 
                bad_stars = zeros(length(obj.snr_edges)-1, length(obj.size_edges)-1, length(new_edges)-1, length(obj.vel_edges)-1, 'like', obj.losses_bad_stars); 
                
                if ~isempty(obj.ecl_edges)
                    
                    start_idx = find(new_edges==obj.ecl_edges(1)); % find the first bin index in the new matrix
                    end_idx = find(new_edges==obj.ecl_edges(end))-1; % find the last bin index in the new matrix

                    if isempty(start_idx) || isempty(end_idx)
                        error('the new edges should include the old edges! '); 
                    end

                    new_star_seconds(:,:,start_idx:end_idx,:) = obj.star_seconds;
                    new_star_seconds_with(:,:,start_idx:end_idx,:) = obj.star_seconds_with_losses;
                    inclusive(:,:,start_idx:end_idx,:,:) = obj.losses_inclusive;
                    exclusive(:,:,start_idx:end_idx,:,:) = obj.losses_exclusive;
                    bad_stars(:,:,start_idx:end_idx,:) = obj.losses_bad_stars;

                end
                
                obj.star_seconds = new_star_seconds;
                obj.star_seconds_with_losses = new_star_seconds_with;
                obj.losses_inclusive = inclusive;
                obj.losses_exclusive = exclusive;
                obj.losses_bad_stars = bad_stars;
                
            end
            
            obj.ecl_edges = new_edges;
            
        end
        
        function add_bin_edges_vel(obj, new_values)
            
            if isempty(new_values)
                return;
            end
            
            mx = ceil(max(new_values(:))./obj.vel_bin_width).*obj.vel_bin_width; 
            if max(obj.vel_edges)>mx
                mx = max(obj.vel_edges); 
            end
            
            new_edges = 0:obj.vel_bin_width:mx;
            
            if ~isequal(obj.vel_edges, new_edges)
            
                % allocate space for the bigger histograms
                new_star_seconds = zeros(length(obj.snr_edges)-1, length(obj.size_edges)-1, length(obj.ecl_edges)-1, length(new_edges)-1, 'like', obj.star_seconds); 
                new_star_seconds_with = zeros(length(obj.snr_edges)-1, length(obj.size_edges)-1, length(obj.ecl_edges)-1, length(new_edges)-1, 'like', obj.star_seconds_with_losses); 
                inclusive = zeros(length(obj.snr_edges)-1, length(obj.size_edges)-1, length(obj.ecl_edges)-1, length(new_edges)-1, length(obj.cut_names), 'like', obj.losses_inclusive); 
                exclusive = zeros(length(obj.snr_edges)-1, length(obj.size_edges)-1, length(obj.ecl_edges)-1, length(new_edges)-1, length(obj.cut_names), 'like', obj.losses_exclusive); 
                bad_stars = zeros(length(obj.snr_edges)-1, length(obj.size_edges)-1, length(obj.ecl_edges)-1, length(new_edges)-1, 'like', obj.losses_bad_stars); 
                
                if ~isempty(obj.vel_edges)
                    
                    start_idx = find(new_edges==obj.vel_edges(1)); % find the first bin index in the new matrix
                    end_idx = find(new_edges==obj.vel_edges(end))-1; % find the last bin index in the new matrix

                    if isempty(start_idx) || isempty(end_idx)
                        error('the new edges should include the old edges! '); 
                    end

                    new_star_seconds(:,:,:,start_idx:end_idx) = obj.star_seconds;
                    new_star_seconds_with(:,:,:,start_idx:end_idx) = obj.star_seconds_with_losses;
                    inclusive(:,:,:,start_idx:end_idx,:) = obj.losses_inclusive;
                    exclusive(:,:,:,start_idx:end_idx,:) = obj.losses_exclusive;
                    bad_stars(:,:,:,start_idx:end_idx) = obj.losses_bad_stars;

                end
                
                obj.star_seconds = new_star_seconds;
                obj.star_seconds_with_losses = new_star_seconds_with;
                obj.losses_inclusive = inclusive;
                obj.losses_exclusive = exclusive;
                obj.losses_bad_stars = bad_stars;
                
            end
            
            obj.vel_edges = new_edges;
            
        end
        
        function idx = findIndex(obj, values, edge_name)
            
            edges = obj.([edge_name, '_edges']); 
            
            if isempty(values)
                idx = 1:length(edges)-1; 
            elseif isscalar(values)
                idx = find(edges<=values, 1, 'last'); 
            elseif length(values)==2
                idx1 = find(edges<=values(1), 1, 'last'); 
                idx2 = find(edges>values(2), 1, 'first')-1; 
                idx = idx1:idx2; 
            else
                error('The values given to "%s" must be empty, or scalar, or a [min,max] vector', edge_name); 
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function showTotalHours(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('velocity', false); 
            input.input_var('axes', [], 'axis'); 
            input.input_var('font_size', 18); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            if input.velocity
                imagesc(input.axes, obj.ecl_edges(1:end-1)+obj.ecl_bin_width/2, obj.vel_edges(1:end-1), squeeze(nansum(nansum(obj.star_seconds,1),2))'./3600);
                c = colormap(input.axes);
                c(1,:) = [1 1 1]; 
                colormap(input.axes, c); 
                
                xlabel(input.axes, 'ecliptic latitude [deg]'); 
                ylabel(input.axes, 'transverse velocity [FSU/s]'); 
                
                hc = colorbar(input.axes); 
                ylabel(hc, 'star hours'); 
                hc.FontSize = input.font_size;
                
                input.axes.ColorScale = 'log';
                
            else
                bar(input.axes, obj.ecl_edges(1:end-1)+obj.ecl_bin_width/2, squeeze(nansum(nansum(nansum(obj.star_seconds,1),2),4))./3600,1);
                xlabel(input.axes, 'ecliptic latitude [deg]'); 
                ylabel(input.axes, 'total star hours'); 
            end
            
            input.axes.FontSize = input.font_size;
            
        end
        
        function showHours(obj, varargin)
            
            import util.text.cs;
            
            input = util.text.InputVars;
            input.input_var('snrs', [], 'signal to noise ratio'); 
            input.input_var('sizes', [], 'stellar sizes'); 
            input.input_var('ecl', [], 'ecliptic latitudes'); 
            input.input_var('vel', [], 'velocity', 'transverse velocity', 'shadow velocity'); 
            input.input_var('cuts', {}); % names of cuts to apply
            input.input_var('exclusive', false); % by default show inclusive losses
            input.input_var('log', false, 'logarithm'); 
            input.input_var('axes', [], 'axis'); 
            input.input_var('font_size', 18); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.snrs) && isempty(input.sizes) && isempty(input.ecl) && isempty(input.vel)
                input.ecl = obj.ecl_edges; % default behavior is to show star-hours along ecliptic latitude
            end
            
            axes_names = {};
            if length(input.snrs)>2, axes_names{end+1} = 'snr'; end
            if length(input.sizes)>2, axes_names{end+1} = 'size'; end
            if length(input.ecl)>2, axes_names{end+1} = 'ecl'; end
            if length(input.vel)>2, axes_names{end+1} = 'vel'; end
            
            if length(axes_names)>2
                error('Cannot plot more than two dimensions!'); 
            end
            
            if ischar(input.cuts)
                input.cuts = {input.cuts};
            end
            
            % temporary variables to pass to getSeconds()
            snrs = input.snrs;
            sizes = input.sizes;
            ecl = input.ecl;
            vel = input.vel;
            
            %%% get the relevant data into a matrix %%%
            
            if length(axes_names)==1 % one dimensional bar chart
                
                edges1 = obj.([axes_names{1} '_edges']); 
                width1 = obj.([axes_names{1}, '_bin_width']); 
                
                hours = zeros(length(edges1)-1, 1, 'like', obj.star_seconds); 
                hours_all = zeros(size(hours), 'like', obj.star_seconds_with_losses); 
                hours_lost = zeros(length(edges1)-1, length(input.cuts)); 
                
                exc_or_inc = 'inclusive';
                if input.exclusive, exc_or_inc = 'exclusive'; end
                
                for ii = 1:length(edges1)-1
                    hours(ii,1) = obj.getSeconds('snrs', snrs, 'sizes', sizes, 'ecl', ecl, 'vel', vel,...
                        axes_names{1}, edges1(ii))/3600; % the last input overrides one of the previous four inputs
                    hours_all(ii,1) = obj.getSeconds('source', 'with_losses', 'snrs', snrs, 'sizes', sizes, 'ecl', ecl, 'vel', vel,...
                        axes_names{1}, edges1(ii))/3600; % the last input overrides one of the previous four inputs
                    
                    for jj = 1:length(input.cuts)
                        
                        if cs(input.cuts{jj}, 'bad_stars')
                            hours_lost(ii,jj) = obj.getSeconds('source', 'bad_stars', ...
                                'snrs', snrs, 'sizes', sizes, 'ecl', ecl, 'vel', vel,...
                                axes_names{1}, edges1(ii))/3600; % the last input overrides one of the previous four inputs
                        else
                            hours_lost(ii,jj) = obj.getSeconds('source', exc_or_inc, 'cut', input.cuts{jj}, ...
                                'snrs', snrs, 'sizes', sizes, 'ecl', ecl, 'vel', vel,...
                                axes_names{1}, edges1(ii))/3600; % the last input overrides one of the previous four inputs
                        end
                        
                    end
                    
                end
                
            elseif length(axes_names)==2 % two dimensional image
                
                edges1 = obj.([axes_names{1} '_edges']); 
                edges2 = obj.([axes_names{2} '_edges']); 
                
                width1 = obj.([axes_names{1}, '_bin_width']); 
                width2 = obj.([axes_names{2}, '_bin_width']); 
                
                
            end
            
            
            %%% now plot the data! %%%
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            if length(axes_names)==1 % one dimensional bar chart
                
                cla(input.axes); 
                
                input.axes.NextPlot = 'replace'; 
                
                x = edges1(1:end-1)+width1/2;
                
                h = bar(input.axes, x, hours, 1); 
                
                input.axes.NextPlot = 'add'; 
                
                hours_all(hours==0) = NaN;
                h2 = plot(input.axes, x, hours_all, '^', 'Color', h.FaceColor, 'MarkerSize', 8, 'LineWidth', 3); 
                
                if ~isempty(hours_lost)
                    h3 = bar(input.axes, x, hours_lost, 0.8); 
                    legend(input.axes, strrep({'star hours', 'without losses', input.cuts{:}}, '_', ' '), 'Location', 'NorthWest'); 
                end
                
                input.axes.NextPlot = 'replace'; 
                
                ylabel(input.axes, 'Star hours'); 
                
                if input.log
                    input.axes.YScale = 'log';
                else
                    input.axes.YScale = 'linear';
                end
                
            elseif length(axes_names)==2 % two dimensional image
                
            end
            
            if cs(axes_names{1}, 'snr')
                xlabel(input.axes, 'Signal to Noise ratio'); 
            elseif cs(axes_names{1}, 'size')
                xlabel(input.axes, 'Stellar radius [FSU]');
            elseif cs(axes_names{1}, 'ecl')
                xlabel(input.axes, 'Ecliptic latitude [deg]');
            elseif cs(axes_names{1}, 'vel')
                xlabel(input.axes, 'Transverse velocity [FSU/s]');
            end
            
            input.axes.FontSize = input.font_size;
            
        end
        
        function printReport(obj, varargin)
            
            import util.text.cs;
            
            input = util.text.InputVars;
            input.input_var('units', 'hours'); % can also choose "hours" or "percent"
            input.input_var('format', 'text'); % can aso choose "latex"
            input.scan_vars(varargin{:}); 
            
            if cs(input.units, 'seconds')
                unit_str = 's'; 
            elseif cs(input.units, 'hours')
                unit_str = 'h'; 
            elseif cs(input.units, 'percent', '%')
                unit_str = '%'; 
            else
                error('Unknown option "%s" to "units" argument. Try using "seconds" or "hours" or "percent"...', input.units); 
            end
            
            if cs(input.format, 'text')
                sep = '|';
                line = '--------------------------+--------------+------------- \n';
                ending = newline;
                code = @(str) str;
            elseif cs(input.format, 'latex')
                sep = '&';
                line = ['\\hline' newline]; 
                ending = ['\\' newline]; 
                code = @(str) sprintf('\\code{%s}', str); 
            else
                error('Unknown "format" option "%s". Use "text" or "latex" instead...', input.format); 
            end
            
            fprintf('%-25s %s inclusive[%s] %s exclusive[%s] %s', 'Cut name', sep, unit_str, sep, unit_str, ending); 
            fprintf(line); 
            
            M = obj.star_seconds; 
            useful = nansum(M(:)); 
            
            M = obj.star_seconds_with_losses;
            total = nansum(M(:)); 
            
            % get the total star hours, not including losses
            
            for ii = 1:length(obj.cut_names)
                
                M = obj.losses_inclusive(:,:,:,:,ii);
                inc = nansum(M(:)); 
                
                if unit_str=='h'
                    inc = inc/3600;
                elseif unit_str=='%'
                    inc = inc/total.*100; 
                end
                
                M = obj.losses_exclusive(:,:,:,:,ii);
                exc = nansum(M(:)); 
                
                if unit_str=='h'
                    exc = exc/3600;
                elseif unit_str=='%'
                    exc = exc/total.*100; 
                end
                
                fprintf('%-25s %s %12.2f %s %12.2f %s', code(obj.cut_names{ii}), sep, inc, sep, exc, ending);
                
                
            end
            
            M = obj.losses_bad_stars;
            stars = nansum(M(:)); 
            
            if unit_str=='h'
                stars= stars/3600;
            elseif unit_str=='%'
                stars = stars/total.*100; 
            end
            
            fprintf('%-25s %s %12.2f %s         ---  %s', code('Bad stars'), sep, stars, sep, ending); 
            
            fprintf(line); 
            
            fprintf('%-25s %s %11.1fh %s %11.1fh %s', 'Good / total time', sep, useful/3600, sep, total/3600, ending); % always show it as hours!  

            
        end
        
    end    
    
end

