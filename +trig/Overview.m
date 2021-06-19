classdef Overview < handle
% Calculate some statistics on data collected from multiple runs. 
% To generate this data, use the input() method, which accepts either 
% Overview or RunSummary objects. It aggregates the star hours and simulated
% events saved in such objects and builds multi-dimensional histograms. 
%
% To use these data, the following functions are most useful: 
%   -showHours: show the total star hours as function of some of the parameters
%               they were observed at, such as stellar size, ecliptic latitude, 
%               transverse velocity, etc. 
%   -numDetections: get an estimate for the number of occultations that would
%                   be seen, based on the model, the efficiency, and the 
%                   number of star hours that were collected. 
%   -showCoverage: plot the occulter model and the coverage as a function 
%                  of the occulter radius. 
% The last two functions use calcCoverage() which uses calcEfficiency(). 
%
% calcEfficiency() compares the number of simulated (injected) events that
% passed the threshold with the total number of simulated events, as a 
% function of the occultation parameters. 
% The results are summed over some parameters, such as stellar radius, as we 
% can assume their distributions are fairly similar across runs. 
% The output is an efficiency fraction as function of velocity and occulter
% radius, which are not integrated as they are connected to the star hours 
% and the occulter size model. 
% 
% calcCoverage() takes the efficiency and multiplies it by the star hours 
% taken at different runs, summing over velocities (with different efficiency
% for each velocity bin) and returns a matrix of coverage. 
% Coverage is the number of square degress that were "scanned" by the survey, 
% and is the reciprocal of the number of objects per square degree. 
% This is given as a function of occulter radius and ecliptic latitude, 
% because different star hours are collected at different latitudes, 
% and because the coverage is higher for larger occulters. 
% The coverage can be compared to the number of objects per square degree
% as predicted by the occulter model. 
%
% EXAMPLES: 
% % ecliptic latitude range in degrees, distance to occulters, in AU:
% >> numDetections('ecliptic latitude', [-5 5], 'distance', 40); 
% >> showCoverage('ecliptic latitude', [-5 5], 'distance', 40, 'log', 1, 'ax', gca); 
% 
% 
% This object also contains a couple of CometModel objects to describe the 
% size distribution of KBOs or Oort cloud objects. 
% The "folders" object contains a vector of the RunFolder objects used to 
% build up the Overview (if created using Scanner's calcOverview()). 
% Two Overview objects can be combined using the input() method on one of them. 

    
    properties(Transient=true)
        
    end
    
    properties % objects
        
        kbos@trig.CometModel; 
        hills@trig.CometModel; 
        oort@trig.CometModel; 
        folders@trig.RunFolder; 
        
    end
    
    properties % inputs/outputs
        
        runtime; % number of seconds of processed data
        batch_counter; % number of batches that were processed
        
        ecl_edges; % ecliptic latitude bin edges used in the star_seconds and losses histograms (degrees)
        vel_edges; % transverse velocity bin edges used in the star_seconds and losses histograms (km/s)
        snr_edges; % S/N bin edges used in the star_seconds and losses histograms
        size_edges; % stellar size edges used in the star_seconds and losses histograms (Fresnel units, estimated at 40 AU)
%         r_edges; % occulter radius edges used in efficiency estimates (km)
        
        b_max; % the maximumum impact parameter used in simulations. This tells us we can integrate -b_max to b_max to get the coverage
        
        star_seconds; % the number of useful seconds. Dims 1- S/N, 2- R, 3- ECL 4- vel
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
        
        % defining the width and maximal values for the histograms
        ecl_bin_width = 2; % ecliptic latitude bins (from -90 to 90)
        vel_bin_width = 1; % transverse velocity bin width
        vel_max = 31; % transverse velocity bin maximum
        snr_bin_width = 1; % S/N bin width
        size_bin_width = 0.1; % stellar size bin width
           
        % these are used to calculate the coverage/number of detections
        lambda_nm = 500; % to calculate the Fresnel scale
        distance_au = 40; % to calculate the Fresnel scale
        
        num_occulters_above_r = 4.4e6; % objects per deg^2 
        r_normalization_km = 0.25; % the radius at which the above is given
        power_law_index = -3.8; % slow of the differential power law
        
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
        
        default_r_edges_fsu = 0.2:0.1:2;
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = Overview(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.Overview')
                if obj.debug_bit>1, fprintf('Overview copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('Overview constructor v%4.2f\n', obj.version); end
                
                obj.kbos = trig.CometModel; 
                
                obj.hills = trig.CometModel('hills'); 
                
                obj.oort = trig.CometModel('oort'); 
                
                obj.reset; 
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.folders = trig.RunFolder.empty;
            
            obj.runtime = 0;
            obj.batch_counter = 0;

            obj.ecl_edges = [];
            obj.vel_edges = [];
            obj.snr_edges = [];
            obj.size_edges = [];

            obj.b_max = []; 
            
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
        
        function val = km2fsu(obj, dist_au)
            
            if nargin<2 || isempty(dist_au)
                dist_au = obj.distance_au; 
            end
           
            val = 1./sqrt(obj.lambda_nm.*1e-12.*dist_au.*1.496e+8/2); 
            
        end
        
        function val = fsu2deg2(obj, dist_au)

            if nargin<2 || isempty(dist_au)
                dist_au = obj.distance_au; 
            end
            
            val = obj.lambda_nm.*1e-12./(dist_au.*1.496e+8)/2*(180/pi).^2; 

        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations / API
        
        function input(obj, summary) % accepts a summary or an overview object...
            
            %%% ingest either Overview or RunSummary object %%%

            if isa(summary, 'trig.RunSummary') && ~isempty(summary)
                
                if isempty(summary.head) || isempty(summary.head.ephem) || isempty(summary.head.ephem.ECL_lat)
                    error('Summary object must have a valid header with an Ephemeris object / ecliptic latitude value!'); 
                end
                
                ecl = summary.head.ephem.ECL_lat; 
                vel = sqrt(sum(summary.head.ephem.getShadowVelocity.^2)); % in km/s
                snr = summary.snr_bin_edges;
                sizes = summary.size_bin_edges;
                
                inclusive = permute(summary.losses_inclusive, [1,2,5,4,3]); 
                exclusive = permute(summary.losses_exclusive, [1,2,5,4,3]); 
                
            elseif isa(summary, 'trig.Overview') && ~isempty(summary)
                
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

%                             for c = 1:length(obj.cut_names)
                                obj.losses_inclusive(idx1(ii), idx2(jj), idx3(kk), idx4(mm), :) = obj.losses_inclusive(idx1(ii), idx2(jj), idx3(kk), idx4(mm), :) + inclusive(ii,jj,kk,mm,:);
                                obj.losses_exclusive(idx1(ii), idx2(jj), idx3(kk), idx4(mm), :) = obj.losses_exclusive(idx1(ii), idx2(jj), idx3(kk), idx4(mm), :) + exclusive(ii,jj,kk,mm,:);
%                             end

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
        
        function [N_total, N_passed] = calcEfficiency(obj, r_edges, distance_au) % a 2D matrix telling what fraction of events we got from injected events
        % Get the total number of simulated events and those that passed.
        % INPUTS: 
        %   -r_edges: bin edges for occulter radius in km. 
        %   -distance_au: the distance to the occulters, in AU. 
        %
        % OUTPUTS:
        %   -N_total: number of simulated events. 
        %   -N_passed: number of events the were detected. 
        %
        % NOTE: output dim1 is velocity and dim2 is occulter radius. 
        %
        % Since we have injected events with the "true" stellar radius 
        % S/N and impact parameter (b), we can assume the detection 
        % efficiency already takes into account these distributions.
        % Thus we can sum all the star seconds over R, b and S/N. 
        %
        % The efficiency matrix should have the following dimensions:
        % 1) transverse velocity (vel): we need to weigh the different star
        %    hours based on events with different velocities.
        % 2) occulter radius (r): we want to say what fraction of events
        %    with r in this range have passed. 
        % 
        % The star hours should be summed along R, and S/N, and then
        % permuted to have the same dimensions, with a singleton 2nd
        % dimension (star seconds are agnostic to occulter radii) and vel
        % on the 1st dimension and ecl on the 3rd dimension.         
        % Then, the product of star-seconds and efficiency can be
        % integrated to give the effective star-time. 
        % That will still need to be multiplied by the scanning speed and 
        % the KBO size distribution, to get the total coverage. 
        
            import util.text.cs;
        
            if nargin<3 || isempty(distance_au)
                distance_au = 40; 
            end
            
            if ischar(distance_au)
                
                if cs(distance_au, 'kbos', 'kuiper belt objects')
                    distance_au = 40; 
                elseif cs(distance_au, 'hills cloud', 'inner oort')
                    distance_au = 3000; 
                elseif cs(distance_au, 'oort cloud')
                    distance_au = 10000; 
                else
                    error('Unknown "distance" option "%s". Use a numeric value in AU, or "KBOs", "Hills" or "Oort"', distance_au); 
                end
                
            end
            
            if nargin<2 || isempty(r_edges)                
                r_edges = obj.default_r_edges_fsu./obj.km2fsu(distance_au); 
            end
            
            r_edges = obj.km2fsu(distance_au).*r_edges; % convert to FSU
            v_edges = obj.km2fsu(distance_au).*obj.vel_edges; % convert to FSU
            
            obj.b_max = nanmax([obj.sim_events.b]); 
            
            N_total = zeros(length(v_edges)-1, length(r_edges)-1, 'like', obj.star_seconds); 
            
            N_passed = zeros(size(N_total), 'like', N_total); 
            
            for ii = 1:length(obj.sim_events)
                
                ev = obj.sim_events(ii); 
                
                if ev.D~=distance_au
                    continue;
                end

                idx_v = find(ev.v>v_edges, 1, 'last'); 
                if idx_v>=length(v_edges)
%                         idx_v = length(v_edges) - 1;
                    continue;
                end

                idx_r = find(ev.r>r_edges, 1, 'last'); 
                if idx_r>=length(r_edges)
%                         idx_r = length(r_edges) - 1;
                    continue;
                end

                % how many events, in total, were simulated
                N_total(idx_v, idx_r) = N_total(idx_v, idx_r) + 1;

                % how many events passed the cut
                % (we may want to replace these with the full candidates
                % that passed classification, too)
                if ev.passed
                    N_passed(idx_v, idx_r) = N_passed(idx_v, idx_r) + 1;
                end
                
            end
            
        end
        
        function [coverage, cov_lower, cov_upper, E, E_l, E_u, N_total, N_passed, v_weights] = calcCoverage(obj, r_edges, distance_au)
        % Get the coverage: the number of square degrees scanned out of the 
        % entire sky, for a given occulter radius and ecliptic latitude. 
        % 
        % INPUTS: 
        %   -r_edges: bin edges for occulter radius in km. 
        %   -distance_au: the distance to the occulters, in AU. 
        %
        % OUTPUTS: 
        %   -coverage: the number of scanned square degrees, reciprocal  
        %              of the number of objects per square degree. 
        %   -cov_lower and cov_upper: limits based on the Poisson errors in
        %                             the efficiency estimate. 
        %   -E, E_l and E_u: efficiency and lower/upper limits. 
        %                    dim1 velocity and dim2 occulter radius.
        %   -N_total and N_passed: outputs of calcEfficiency. 
        %   -v_weights: the number of seconds in each velocity bin. 
        %               dim 1 is velocity, dim2 is ecliptic latitude. 
        % 
        % NOTE: the coverage has 2 dimensions: dim1 is occulter radius and 
        %       dim2 is ecliptic latitude. 
        
        
            if isempty(obj.star_seconds)
                disp('No star seconds have been accumulated yet!'); 
                coverage = []; 
                return;
            end
            
            if nargin<3 || isempty(distance_au)
                distance_au = 40; 
            end
            
            if ischar(distance_au)
                
                if cs(distance_au, 'kbos', 'kuiper belt objects')
                    distance_au = 40; 
                elseif cs(distance_au, 'hills cloud', 'inner oort')
                    distance_au = 3000; 
                elseif cs(distance_au, 'oort cloud')
                    distance_au = 10000; 
                else
                    error('Unknown "distance" option "%s". Use a numeric value in AU, or "KBOs", "Hills" or "Oort"', distance_au); 
                end
                
            end
            
            if nargin<2 || isempty(r_edges)                
                r_edges = obj.default_r_edges_fsu./obj.km2fsu(distance_au); 
            end
            
            % the result should have 2D (dim1 velocity and dim2 occulter radius)
            [N_total, N_passed] = obj.calcEfficiency(r_edges, distance_au);
            
            E = N_passed./N_total; 
            E(N_total==0) = NaN; % where there are no events at all, efficiency is NaN
%             E = fillmissing(E, 'next'); 
            
            
            [N_passed_lower, N_passed_upper] = util.stat.poisson_errors(N_passed, .32); 
            
            E_l = N_passed_lower./N_total; 
            E_l(N_total==0) = NaN; % where there are no events at all, efficiency is NaN
%             E_l = fillmissing(E_l, 'next'); 
            
            
            E_u = N_passed_upper./N_total; 
            E_u(N_total==0) = NaN; % where there are no events at all, efficiency is NaN
%             E_u = fillmissing(E_u, 'next'); 
            
            
%             E = ones(size(E)); % debugging only

            T = nansum(nansum(obj.star_seconds,1),2); % star-seconds integrated over all S/N and stellar sizes
            T = permute(T, [4,1,3,2]); % arrange velocity into the 1st dim, leave dim 2 as scalar (for r) and dim 3 for ecl 
            
            v_shift = 0; % add a shift for e.g., KBO intrinsic orbital velocity (in km/s)
            v = obj.km2fsu(distance_au).*(obj.vel_edges(1:end-1)+obj.vel_bin_width/2 - v_shift); % velocity bin centers, translated from km/s to FSU/s
            v = abs(util.vec.tocolumn(v)); % put it on the same dimension as the matrices and take the abs() in case we got negative velocities
            
            b = 2.*obj.b_max; % the range of impact parameters, integrated from -b_max to +b_max
            
            coverage = nansum(fillmissing(E, 'next').*T.*b.*v, 1); % integral of efficiency and time and impact parameter, for each velocity
            coverage = permute(coverage, [2,3,1]); % remove the velocity dimension we've integrated on, and leave radius and ecliptic latitute
            coverage = coverage.*obj.fsu2deg2(distance_au); 
            
            cov_lower = nansum(fillmissing(E_l, 'next').*T.*b.*v, 1); % integral of efficiency and time and impact parameter, for each velocity
            cov_lower = permute(cov_lower, [2,3,1]); % remove the velocity dimension we've integrated on, and leave radius and ecliptic latitute
            cov_lower = cov_lower.*obj.fsu2deg2(distance_au); 
            
            cov_upper = nansum(fillmissing(E_u, 'next').*T.*b.*v, 1); % integral of efficiency and time and impact parameter, for each velocity
            cov_upper = permute(cov_upper, [2,3,1]); % remove the velocity dimension we've integrated on, and leave radius and ecliptic latitute
            cov_upper = cov_upper.*obj.fsu2deg2(distance_au); 
            
            v_weights = permute(T, [1,3,2]); % dim 1 is velocity, dim 2 is ecliptic latitude
            
        end
        
        function [N, N_err] = numDetections(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ecl', [-5,5], 'ecliptic latitude', 'ecliptic limits'); 
            input.input_var('r_edges', []);
            input.input_var('distance', 40, 'distance_au', 'dist_au'); 
            input.scan_vars(varargin{:}); 
            
            if length(input.ecl)~=2
                error('Give the ecliptic latitude parameter as [min,max] values'); 
            end
            
            ecl_idx1 = find(obj.ecl_edges>=input.ecl(1), 1, 'first');
            ecl_idx2 = find(obj.ecl_edges<input.ecl(2), 1, 'last') - 1;
            
            if isempty(ecl_idx1) || isempty(ecl_idx2) || ecl_idx1>ecl_idx2
                error('No data for the ecliptic latitude range %s', util.text.print_vec(input.ecl, ' to ')); 
            end
            
            r = util.vec.tocolumn(input.r_edges);
            r = r(1:end-1) + diff(r)/2; 
            
            if ischar(input.distance)
                
                if cs(input.distance, 'kbos', 'kuiper belt objects')
                    input.distance = 40; 
                elseif cs(input.distance, 'hills cloud', 'inner oort')
                    input.distance = 3000; 
                elseif cs(input.distance, 'oort cloud')
                    input.distance = 10000; 
                else
                    error('Unknown "distance" option "%s". Use a numeric value in AU, or "KBOs", "Hills" or "Oort"', input.distance); 
                end
                
            end
            
            if input.distance<=200
                comets = obj.kbos; 
            elseif input.distance<=5000
                comets = obj.hills;
            else
                comets = obj.oort; 
            end
            
            if isempty(input.r_edges)                
                input.r_edges = obj.default_r_edges_fsu./obj.km2fsu(input.distance); 
            end
            
            [coverage, cov_lower, cov_upper] = obj.calcCoverage(input.r_edges, input.distance);
            
            % integrate over the relevant ecliptic latitudes
            C = nansum(coverage(:,ecl_idx1:ecl_idx2),2);
            C_l = nansum(cov_lower(:,ecl_idx1:ecl_idx2),2);
            C_u = nansum(cov_upper(:,ecl_idx1:ecl_idx2),2); 
            
            [C_l,C,C_u]
            
%             C_p = C_u - C; % coverage plus 
%             C_m = C - C_l; % coverage minus
            
            % get the KBO abundance in each r interval
            [n, n_l, n_u] = comets.numDensityIntervals(input.r_edges);
            
%             n_p = n_u - n; % number density plus
%             n_m = n - n_l; % number density minus
            
            N = nansum(C.*n, 1); 
            
            N_err = [nansum(C_l.*n_l,1), nansum(C.*n_l,1), nansum(C_u.*n_l,1);
                     nansum(C_l.*n,1), nansum(C.*n,1), nansum(C_u.*n,1);
                     nansum(C_l.*n_u,1), nansum(C.*n_u,1), nansum(C_u.*n_u,1)];
            
            
%             N_p = sqrt(nansum(C_p.*n).^2 + nansum(C.*n_p).^2);
%             N_m = sqrt(nansum(C_m.*n).^2 + nansum(C.*n_m).^2);
%             
%             N_l = N - N_m; 
%             N_u = N + N_p;
            
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
                    legend(input.axes, strrep({'star hours', 'without losses', input.cuts{:}}, '_', ' '), 'Location', 'NorthEast'); 
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
        
        function showEfficiency(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ecl', [-5,5], 'ecliptic latitude', 'ecliptic limits'); 
            input.input_var('r_edges', []);
            input.input_var('distance', 40, 'distance_au', 'dist_au'); 
            input.input_var('parent', [], 'figure'); 
            input.input_var('font_size', 18); 
            input.scan_vars(varargin{:}); 
            
            if length(input.ecl)~=2
                error('Give the ecliptic latitude parameter as [min,max] values'); 
            end
            
            ecl_idx1 = find(obj.ecl_edges>=input.ecl(1), 1, 'first');
            ecl_idx2 = find(obj.ecl_edges<input.ecl(2), 1, 'last') - 1;
            
            if isempty(ecl_idx1) || isempty(ecl_idx2) || ecl_idx1>ecl_idx2
                error('No data for the ecliptic latitude range %s', util.text.print_vec(input.ecl, ' to ')); 
            end
            
            if ischar(input.distance)
                
                if cs(input.distance, 'kbos', 'kuiper belt objects')
                    input.distance = 40; 
                elseif cs(input.distance, 'hills cloud', 'inner oort')
                    input.distance = 3000; 
                elseif cs(input.distance, 'oort cloud')
                    input.distance = 10000; 
                else
                    error('Unknown "distance" option "%s". Use a numeric value in AU, or "KBOs", "Hills" or "Oort"', input.distance); 
                end
                
            end
            
            if isempty(input.r_edges)
                input.r_edges = obj.default_r_edges_fsu./obj.km2fsu(input.distance); 
            end
            
            r = util.vec.tocolumn(input.r_edges);
            r = r(1:end-1) + diff(r)/2; 
            
            [~,~,~,E, E_l, E_u, N_total, N_passed, v_weights] = obj.calcCoverage(input.r_edges, input.distance);
            
            v_weights = nansum(v_weights(:,ecl_idx1:ecl_idx2),2); % sum over all runs in the ecliptic range
            
            E1 = nansum(N_passed)./nansum(N_total); % simple sum over all velocities
            E1 = fillmissing(E1, 'constant', 0);
            
            E2 = util.vec.weighted_average(fillmissing(E, 'next'),v_weights,1); % sum according to how many hours were spent in each velocity
            
            [N_lower, N_upper] = util.stat.poisson_errors(nansum(N_passed), 0.32); 
            
            E1_l = N_lower./nansum(N_total); % simple sum over all velocities
            E1_l = fillmissing(E1_l, 'constant', 0);
            
%             E2_l = util.vec.weighted_average(fillmissing(E_l, 'next'),v_weights,1); % sum according to how many hours were spent in each velocity
            
            E1_u = N_upper./nansum(N_total); % simple sum over all velocities
            E1_u = fillmissing(E1_u, 'constant', 0);
            E1_u(isinf(E1_u)) = 0; 
            
%             E2_u = util.vec.weighted_average(fillmissing(E_u, 'next'),v_weights,1); % sum according to how many hours were spent in each velocity
            
            if isempty(input.parent)
                input.parent = gcf;
            end
            
            delete(input.parent.Children); 
            
            ax1 = axes('Parent', input.parent); 
            
%             errorbar(ax1, r, E1*100, (E1-E1_l)*100, (E1_u-E1)*100, '+-', 'LineWidth', 2, 'DisplayName', 'average efficiency'); 
            [h_line, h_fill] = util.plot.shaded(r, E1'*100, [(E1-E1_l)'*100, (E1_u-E1)'*100], 'ax', ax1, 'positive', 1); 

            hold(ax1, 'on'); 
            
            h_line.DisplayName = 'average efficiency';
            h_fill.DisplayName = '68% confidence'; 
            
            plot(ax1, r, E2*100, 'go', 'LineWidth', 2, 'DisplayName', 'weighted by velocities'); 
            
            for ii = 1:length(E1)
                text(ax1, mean([r(ii),input.r_edges(ii)]), double(E1(ii)*100), sprintf('%d/%d', sum(N_passed(:,ii),1), sum(N_total(:,ii),1)), ...
                    'FontSize', 14, 'Rotation', 90, 'HorizontalAlignment', 'left', 'Color', 'r'); 
            end
            
            hold(ax1, 'off'); 
            
            box(ax1, 'on'); 
            
            ytickformat(ax1, '%d%%'); 
            xtickformat(ax1, '%3.1f km'); 
            
            xlabel(ax1, 'Occulter radius'); 
            ylabel(ax1, 'Detection efficiency'); 
            
            if ax1.XTick(1)==0
                ax1.XTick = ax1.XTick(2:end);
            end
            
            ax1.YLim = [0 100]; 
            grid(ax1, 'on'); 
            
            ax1.FontSize = input.font_size; 
            
            legend(ax1, 'Location', 'SouthEast'); 
            
            
            %%%%% velocity distribution %%%%%
            
            ax2 = axes('Parent', input.parent, 'Position', [0.175 0.55 0.3 0.3]); 
            
            bar(ax2, obj.vel_edges(1:end-1)-obj.vel_bin_width/2, v_weights/3600); 
            
            ylabel('Star hours'); 
            xlabel('Velocity [km/s]'); 
            
            ax2.FontSize = input.font_size - 2;
            
        end
        
        function showCoverage(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ecl', [-5,5], 'ecliptic latitude', 'ecliptic limits'); 
            input.input_var('r_edges', []);
            input.input_var('distance', 40, 'distance_au', 'dist_au'); 
            input.input_var('log', true, 'logarithm'); 
            input.input_var('axes', [], 'axis'); 
            input.input_var('font_size', 18); 
            input.scan_vars(varargin{:}); 
            
            if length(input.ecl)~=2
                error('Give the ecliptic latitude parameter as [min,max] values'); 
            end
            
            ecl_idx1 = find(obj.ecl_edges>=input.ecl(1), 1, 'first');
            ecl_idx2 = find(obj.ecl_edges<input.ecl(2), 1, 'last') - 1;
            
            if isempty(ecl_idx1) || isempty(ecl_idx2) || ecl_idx1>ecl_idx2
                error('No data for the ecliptic latitude range %s', util.text.print_vec(input.ecl, ' to ')); 
            end
            
            if ischar(input.distance)
                
                if cs(input.distance, 'kbos', 'kuiper belt objects')
                    input.distance = 40; 
                elseif cs(input.distance, 'hills cloud', 'inner oort')
                    input.distance = 3000; 
                elseif cs(input.distance, 'oort cloud')
                    input.distance = 10000; 
                else
                    error('Unknown "distance" option "%s". Use a numeric value in AU, or "KBOs", "Hills" or "Oort"', input.distance); 
                end
                
            end
            
            if input.distance<=200
                comets = obj.kbos; 
            elseif input.distance<=5000
                comets = obj.hills;
            else
                comets = obj.oort; 
            end
            
            if isempty(input.r_edges)
                input.r_edges = obj.default_r_edges_fsu./obj.km2fsu(input.distance); 
            end
            
            r = util.vec.tocolumn(input.r_edges);
            r = r(1:end-1) + diff(r)/2; 
            
            [coverage, cov_lower, cov_upper] = obj.calcCoverage(input.r_edges, input.distance);
            
            star_hours = obj.star_seconds(:,:,ecl_idx1:ecl_idx2,:); 
            star_hours = nansum(star_hours(:))./3600; 
            
            C = nansum(coverage(:,ecl_idx1:ecl_idx2),2); 
            C_l = nansum(cov_lower(:,ecl_idx1:ecl_idx2),2); 
            C_u = nansum(cov_upper(:,ecl_idx1:ecl_idx2),2); 
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
%             errorbar(input.axes, r, 1./C, 1./C-1./C_l, 1./C_u-1./C); 
            
            hold_state = input.axes.NextPlot;

            h_cov = plot(input.axes, r, 1./C, '-', 'LineWidth', 3); 
            h_cov.DisplayName = sprintf('Inverse coverage (%d star hours)', round(star_hours));
            
            input.axes.NextPlot = 'add';
            
            h_c_u = plot(input.axes, r, 1./C_u, '--', 'LineWidth', 1.5); 
            h_c_u.DisplayName = '1\sigma confidence interval'; 
            
            h_c_l = plot(input.axes, r, 1./C_l, '--', 'LineWidth', 1.5, 'Color', h_c_u.Color); 
            h_c_l.HandleVisibility = 'off'; 
            
            comets.show('r_edges', input.r_edges); 
            
            input.axes.NextPlot = hold_state;
            
            if input.log
                input.axes.YScale = 'log';
            else
                input.axes.YScale = 'linear';
            end
            
            xlabel(input.axes, 'Occulter radius [km]'); 
            ylabel(input.axes, 'Number density [deg^{-2}]'); 
            
            input.axes.FontSize = input.font_size; 
            
            legend(input.axes, 'Location', 'NorthEast'); 
            
        end
        
        function printReport(obj, varargin)
            
            import util.text.cs;
            
            input = util.text.InputVars;
            input.input_var('format', 'text'); % can aso choose "latex"
            input.scan_vars(varargin{:}); 
            
            fprintf('\n'); 
            
            if cs(input.format, 'text')
                sep = '|';
                line = '-----------------------------+-------------------+------------------- \n';
                ending = newline;
                code = @(str) str;
                percent = '%'; 
            elseif cs(input.format, 'latex')
                sep = '&';
                line = ['\\hline' newline]; 
                ending = ['\\' newline]; 
                code = @(str) sprintf('\\code{%s}', str); 
                percent = '\%'; 
            else
                error('Unknown "format" option "%s". Use "text" or "latex" instead...', input.format); 
            end
            
            fprintf('%-28s %s   inclusive[h]    %s   exclusive[h]    %s', 'Cut name', sep, sep, ending); 
            fprintf(line); 
            
            M = obj.star_seconds; 
            useful = nansum(M(:)); 
            
            M = obj.star_seconds_with_losses;
            total = nansum(M(:)); 
            
            % get the total star hours, not including losses
            
            for ii = 1:length(obj.cut_names)
                
                M = obj.losses_inclusive(:,:,:,:,ii);
                inc = nansum(M(:)); 
                inc_h = inc/3600;
                inc_p = inc/total*100; 
                
                M = obj.losses_exclusive(:,:,:,:,ii);
                exc = nansum(M(:)); 
                exc_h = exc/3600;
                exc_p = exc/total*100; 
                
                fprintf('%-28s %s %8.2f (%5.2f%s) %s %8.2f (%5.2f%s) %s', code(obj.cut_names{ii}), sep, inc_h, inc_p, percent, sep, exc_h, exc_p, percent, ending);
                
                
            end
            
            M = obj.losses_bad_stars;
            stars = nansum(M(:)); 
            stars_h = stars/3600;
            stars_p = stars/total*100;
            
            fprintf('%-28s %s %8.2f (%5.2f%s) %s         ---       %s', code('Bad stars'), sep, stars_h, stars_p, percent, sep, ending); 
            
            fprintf(line); 
            
            fprintf('%-28s %s %8.1fh (%5.2f%s) out of %8.1fh   %s', 'good/total star hours', sep, useful/3600, useful/total.*100, percent, total/3600, ending); % always show it as hours!  

            
        end
        
    end    
    
end

