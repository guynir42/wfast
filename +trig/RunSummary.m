classdef RunSummary < handle
% Lightweight object that summarizes the parameters, the data quality and 
% the key results from an occultation scanning run. 
% It is produced by the EventFinder using the produceSummary() method. 
% 
% The parameters for each relevant object are saved as structs, "finder_pars", 
% "store_pars", and "checker_pars". The header info is save in "head" as usual. 
% 
% An important result is the histogram of star-time per S/N and R_star. 
% This is stored without separating into bins for each star, instead we 
% histogram over each S/N value and each star's Fresnel size. 
% This is done so we can stack the total hours or bin them by e.g., 
% ecliptic latitude later on. 
% 
% Other interesting results are the "snr_values" saved for each batch in the 
% run, the cut values (again, histogrammed by values only, not by star), and 
% the parameters of failed/passed simulated events. 
% That last one can be used to find what occultation parameters manage to 
% pass our threshold (and it can be compared to the theoretical values). 
% accumulating these events can help define the real thresholds for different
% types of occultations. 
% 
    
    properties % objects
        
        head;
        
        finder_pars; 
        store_pars;
        checker_pars;
        
    end
    
    properties % inputs/outputs
        
        snr_values; % best S/N value for each batch
        runtime; % number of seconds of processed data
        batch_counter; % number of batches that were processed
        total_batches; % total batches in the run (optional)
        
        snr_bin_edges; % the edges used in the star_seconds and losses histograms
        size_bin_edges; % the edges used in the star_seconds and losses histograms
        fwhm_hist; % histogram of the number of seconds spent at each FWHM value (in steps of 0.1") 
        fwhm_edges; % histogram edges
        seeing_log; % FWHM measurements from the model PSF, in arcsec, for each batch
        airmass_log; % calculate the airmass based on header data and julian dates
        background_log; % get the average background value for each batch
        juldates_log; % julian date of the middle of each batch
        
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
        
        black_list_stars; % a vector of star numbers (from the full list of stars?)
        good_stars; % a list of the stars that passed the threshold after the burn-in period
        
        % structs describing simulated events that passed/failed the threshold
        sim_events;
        
        star_snr; 
        star_sizes;
        size_snr_coeffs;
        
        num_events_expected; % use simulations to estimate this
        
    end
    
    properties % switches/controls
        
        % these are used to histogram stars into Fresnel scale bins
        % We use the Fresnel scale calculated at 40 AU (KBOs)
        R_bin_width = 0.1;
        R_bin_max = 5; % stars bigger than this are not saved in the summary... 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
        
        version = 1.00;
        
    end
    
    properties(Hidden=true, Transient=true)
        
        results; 
        axes_names;
        axes_indices; 
        type_result=''; 
        R_hist;
        R_edges;
        r_hist;
        r_edges;
        b_hist;
        b_edges;
        v_hist;
        v_edges;
        
        
    end
    
    methods % constructor
        
        function obj = RunSummary(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.RunSummary')
                if obj.debug_bit>1, fprintf('RunSummary copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('RunSummary constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
       
        function inputHours(obj, hours, stellar_sizes)
            
            if nargin<2 || isempty(hours) || ~isa(hours, 'trig.StarHours')
                error('Must supply a StarHours object to this method'); 
            end
            
            if nargin<3 || isempty(stellar_sizes) || length(stellar_sizes)~=size(hours.histogram,2)
                error('stellar_sizes input must match the number of stars in the star-hour histograms!'); 
            end
            
            obj.snr_bin_edges = hours.snr_bin_edges;
            obj.size_bin_edges = 0:obj.R_bin_width:obj.R_bin_max; 
            
            for ii = 1:length(obj.size_bin_edges)-1
                
                % indices of stars in this R bin
                if ii==length(obj.size_bin_edges)-1
                    idx = ~(stellar_sizes<obj.size_bin_edges(ii)); % includes everything in the last bin, everything bigger than that, and all NaN sizes
                else
                    idx = stellar_sizes>=obj.size_bin_edges(ii) & stellar_sizes<obj.size_bin_edges(ii+1);
                end
                
                obj.star_seconds(:,ii) = nansum(hours.histogram(:,idx),2);
                obj.star_seconds_with_losses(:,ii) = nansum(hours.histogram_with_losses(:,idx),2);
                obj.losses_exclusive(:,ii,:) = nansum(hours.losses_exclusive(:,idx,:), 2);
                obj.losses_inclusive(:,ii,:) = nansum(hours.losses_inclusive(:,idx,:), 2);
                obj.losses_bad_stars(:,ii) = nansum(hours.losses_bad_stars(:,idx), 2);
                obj.runtime = hours.runtime; 
            
            end
            
        end
        
        function val = calcSimStats(obj, varargin) % show the detection statistics for simulated events
            
            import util.text.cs; 
            
            input = util.text.InputVars;
            input.input_var('result', 'fraction'); % can choose "fraction" or "snr"
            input.input_var('distance', 'kbos'); % which population we want to probe? KBOs? Hills (inner Oort) or Oort?
            input.input_var('star_size_bin', 0.5); % size of bin for R (star radius)
            input.input_var('radius_bin', 0.1); % size of bin for r (occulter radius)
            input.input_var('impact_bin', 0.2); % size of bin for b (impact parameter)
            input.input_var('velocity_bin', 5); % size of bin for v (velocity)
            input.input_var('x_axis', 'r'); % which parameter to put in the x-axis of each subframe (dim 2)
            input.input_var('y_axis', 'b'); % which parameter to put in the y-axis of each subframe (dim 1)
            input.input_var('horizontal', 'v'); % which parameter to put in the horizontal spread of frames (dim 3)
            input.input_var('vertical', 'R'); % which parameter to put in the vertical spread of frames (dim 4)
            input.scan_vars(varargin{:}); 
            
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
            
            val = ''; 
            
            %%%%%%%%%%% find the bin edges and distributions of each variable %%%%%%%%%%%%%%%%%
            events = obj.sim_events;
            events = events([events.D]==input.distance); 
            
            if cs(input.result, 'fraction')
                events = events; % mix the two together (the failed events have detect_snr = []); 
                mode = 1; 
            elseif cs(input.result, 'snr')
                events = events(logical([events.passed])); % take only the passing events for S/N calculation 
                mode = 2; 
            else
                error('Unkown option "%s" to input "result". Try using "fraction" or "snr". ', input.result); 
            end
            
            [obj.R_hist, obj.R_edges] = histcounts([events.R]+0.001, 'BinWidth', input.star_size_bin);
            [obj.r_hist, obj.r_edges] = histcounts([events.r]+0.001, 'BinWidth', input.radius_bin);
            [obj.b_hist, obj.b_edges] = histcounts([events.b]+0.001, 'BinWidth', input.impact_bin);
            [obj.v_hist, obj.v_edges] = histcounts([events.v]+0.001, 'BinWidth', input.velocity_bin);

            %%%%%%%%%% setup the size of the matrix for the result %%%%%%%%%%%%%%%%%
            
            % name of the dimension (single letter)
            dim1 = input.y_axis;
            dim2 = input.x_axis;
            dim3 = input.horizontal;
            dim4 = input.vertical;
            
            obj.axes_names = {dim1,dim2,dim3,dim4}; 
            obj.axes_indices.dim1 = dim1;
            obj.axes_indices.dim2 = dim2;
            obj.axes_indices.dim3 = dim3;
            obj.axes_indices.dim4 = dim4;
            
            % edge vectors
            e1 = obj.([dim1 '_edges']); 
            e2 = obj.([dim2 '_edges']); 
            e3 = obj.([dim3 '_edges']); 
            e4 = obj.([dim4 '_edges']); 
            
            % adjust the bins just so they don't fall on the uniformly sampled grid:
            e1 = e1 - 0.001;
            e2 = e2 - 0.001;
            e3 = e3 - 0.001;
            e4 = e4 - 0.001;
            
            % number of bins (equal to number of edges minus one)
            Ndim1 = length(e1)-1; 
            Ndim2 = length(e2)-1; 
            Ndim3 = length(e3)-1; 
            Ndim4 = length(e4)-1; 
            
            val = zeros(Ndim1, Ndim2, Ndim3, Ndim4, 'single'); 
            
            %%%%%%%%%% start binning the data %%%%%%%%%%%%%%%%%%

            for ii = 1:Ndim1
                
                subset1 = events([events.(dim1)]>=e1(ii) & [events.(dim1)]<e1(ii+1)); % choose only events inside this range
                
                for jj = 1:Ndim2
                    
                    subset2 = subset1([subset1.(dim2)]>=e2(jj) & [subset1.(dim2)]<e2(jj+1)); % choose only events inside this range
                    
                    for kk = 1:Ndim3
                        
                        subset3 = subset2([subset2.(dim3)]>=e3(kk) & [subset2.(dim3)]<e3(kk+1)); % choose only events inside this range
                    
                        for mm = 1:Ndim4
                            
                            subset4 = subset3([subset3.(dim4)]>=e4(mm) & [subset3.(dim4)]<e4(mm+1)); % choose only events inside this range

                            if mode==1 % fraction mode
                                
                                if isempty(subset4)
                                    val(ii,jj,kk,mm) = NaN; % don't have a better idea what to put here
                                else
                                    passed = nnz([subset4.passed]); 
                                    total = numel(subset4); 

                                    val(ii,jj,kk,mm) = passed./total; 
                                end
                                
                            elseif mode==2 % snr mode
                                
                                if isempty(subset4)
                                    val(ii,jj,kk,mm) = 0; % don't have a better idea what to put here
                                else
                                    val(ii,jj,kk,mm) = mean([subset4.detect_snr]./[subset4.star_snr].*10); % adjust all S/N values to a standard star with S/N=10
                                end
                                
                            end % mode

                        end % for mm
                    
                    end % for kk
                    
                end % for jj
                
            end % for ii
            
            obj.results = val; 
            obj.type_result = input.result; % keep track of what we calculated
            
            if nargin==0
                clear val;
            end
                                    
        end
        
        function func = getRateFunction(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('R_bin', 0.25); 
            input.input_var('SNR_bin', 1); 
            input.input_var('plot', 0); 
            input.input_var('axes', [], 'axis'); % axes to draw to if plotting is enabled
            input.input_var('font_size', 18); % fonts on the axes
            input.scan_vars(varargin{:}); 
            
            et = obj.sim_events; % events total
            ep = obj.sim_events([obj.sim_events.passed]); 
            
            [N_total, R_edges, SNR_edges] = histcounts2([et.R], [et.star_snr], 'BinWidth', [input.R_bin, input.SNR_bin]); 
            
            N_passed = histcounts2([ep.R], [ep.star_snr], R_edges, SNR_edges);
            
            x = R_edges(1:end-1) + 0.5*input.R_bin; 
            y = SNR_edges(1:end-1) + 0.5*input.SNR_bin; 
            
            [X,Y] = meshgrid(x, y);
            
            A = [ones(numel(X),1) X(:) X(:).^2 X(:).*Y(:) Y(:) Y(:).^2];
            
            R = (N_passed./N_total)'; % the detection rate
            B = (R); % the rate is the result we want to fit to (flip the axes to put R on x and SNR on y)
            B = B(:); 
            
            w = sqrt(N_total)'; % weights are proportional to the number of detections
            w = w(:); 
            
            idx = w>0; 
            
            B = B(idx); 
            A = A(idx,:); 
            w = w(idx); 
            
            coeffs = lscov(A, B, w); 
            
%             func = str2func(sprintf('@(R,S) (%10.8f %+10.8f.*R   %+10.8f.*S).^2', coeffs(1), coeffs(2), coeffs(3))); 
            func = str2func(sprintf('@(R,S) %10.8f %+10.8f.*R %+10.8f.*R.^2 %+10.8f.*R.*S %+10.8f.*S %+10.8f.*S.^2', ...
                coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6))); 
            
            if input.plot
                
                
                if isempty(input.axes)
                    input.axes = gca;
                end
                
                util.plot.show(func(X,Y), 'fancy', 'on', 'xvalues', x, 'yvalues', y); 
                
                axis(input.axes, 'normal', 'xy'); 
                
                hold(input.axes, 'on'); 
                
                scatter(input.axes, X(idx), Y(idx), w(:), R(idx), 'filled'); % show the measurements as filled circles (color indicates rate, size indicates weight)
                
                hold(input.axes, 'off'); 
                
            end
            
        end
        
        function coverage = getCoverageDegSquare(obj, varargin) % how many deg^2 did we cover (multiply this by angular surface density to get detection number)
        % Some notes: the first value of radius_axis defines the size where
        % the coverage is normalized, with the assumption the radius goes
        % down as the power_law input. 
            
            input = util.text.InputVars;
            input.input_var('radius_axis', 0.5:0.1:2); % occulter radius values
            input.input_var('power_law', 3.5); % index of the occulter radius values
            input.input_var('impact_axis', 0:0.1:2); % impact parameter values
            input.input_var('velocity', []); % shadow veclocity in km/sec! If empty, will use the header data
            input.input_var('range_v', 3); % include values above/below the nominal velocity in this range
            input.input_var('distance_au', 40);
            input.input_var('wavelength_nm', 500); 
            input.input_var('filename', 'sim_results_generic'); 
            input.scan_vars(varargin{:}); 
            
            if ~isempty(input.velocity)
                v = input.velocity; 
            else
                v = obj.head.ephem.getShadowVelocity; 
                v = sqrt(sum(v.^2)); 
            end
            
            % get the simulation results as a table
            res = load(fullfile(getenv('DATA'), 'WFAST/occultations', input.filename)); 
            res = res.results.sim_events; % table! 
            
            % Fresnel Scale Units: 
            FSU = sqrt(input.wavelength_nm.*1e-12*input.distance_au.*150e6./2); % physical Fresnel scale, in km!
            fsu = sqrt(input.wavelength_nm.*1e-12./(2.*input.distance_au.*150e6))*180/pi; % angular Fresnel scale, in degrees!
            fsu_default = sqrt(500e-9./(2.*40.*150e9))*180/pi;
            
            % this is the only result from the actual run we are summarizing
            H = obj.star_seconds; % histogram of star-seconds per S/N and R
            
            % the S/N axis is defined by the existing histograms
            S_axis = obj.snr_bin_edges; 
            dS = median(diff(S_axis)); 
%             S_axis = S_axis(1:end-1)+dS/2; 
            
            % the stellar size is defined by the existing histograms
            R_correction = fsu./fsu_default; % in the default case this is 1 (assume R was calculated for each star using the defaults)
            R_axis = R_correction.*obj.size_bin_edges; 
            dR = median(diff(R_axis)); 
%             R_axis = R_axis(1:end-1) + dR/2; 
            
            % the occulter radius distribution is given by the user
            r_axis = input.radius_axis; 
            
            % the impact parameter distribution is given by the user
            b_axis = input.impact_axis; 
            
            v_axis = (v + [-input.range_v, input.range_v])./FSU;  % velocity range in Fresnel scale units
            v_deg = v./(input.distance_au.*150e6).*180./pi; % scanning speed, in deg/sec
            
            prog = util.sys.ProgressBar;
            
            if obj.debug_bit, disp('Calculating coverage!'); end
            
            coverage = 0; 
            
            for ii = 1:length(v_axis)-1
                
                values_v = res(res.v>=v_axis(ii) & res.v<v_axis(ii+1),:); 
                if isempty(values_v), continue; end
                
                prog.start(length(b_axis)-1);
                
                for jj = 1:length(b_axis)-1
                    
                    values_b = values_v(values_v.b>=b_axis(jj) & values_v.b<b_axis(jj+1),:); 
                    if isempty(values_b), continue; end
                    
                    for kk = 1:length(r_axis)-1
                        
                        values_r = values_b(values_b.r>=r_axis(kk) & values_b.r<r_axis(kk+1),:); 
                        if isempty(values_r), continue; end
                        
                        for mm = 1:length(R_axis)-1
                            
                            values_R = values_r(values_r.R>=R_axis(mm) & values_r.R<R_axis(mm+1),:); 
                            if isempty(values_R), continue; end
                            
                            for nn = 1:length(S_axis)-1
                                
                                values = values_R(values_R.star_snr>=S_axis(nn) & values_R.star_snr<S_axis(nn+1),:); 
                                if isempty(values), continue; end
                                
                                N_total = height(values); 
                                N_passed = height(values(values.passed,:));
                                
                                dt = H(nn,mm); % star-seconds that were observed on this S/N and R
                                
                                r_rate = (r_axis(1)./r_axis(kk)).^input.power_law; % occurance rate of occulters of this radius, relative to the minimal size! 
                                
                                db = (b_axis(jj+1)-b_axis(jj))*fsu; % size of impact parameter bin, in degrees
                                
                                coverage = coverage + N_passed/N_total*r_rate*v_deg*db*dt; % units of deg^2
                                
                            end % for nn (S)
                            
                        end % for mm (R)
                        
                    end % for kk (r)
                    
                    if obj.debug_bit, prog.showif(jj); end
                    
                end % for jj (b)
                
            end % for ii (v)
            
        end
            
        function val = getNumDetections(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('coverage', []); % if this is left empty, will call getCoverageDegSquare with the other parameters
            input.input_var('density', 1.1e7); % density of objects per deg^2. from ref: https://ui.adsabs.harvard.edu/abs/2012ApJ...761..150S/abstract
            input.input_var('norm_radius', 0.25); % normalization radius in km
            input.input_var('radius_axis', 0.5:0.1:2); % occulter radius values
            input.input_var('power_law', 3.5); % index of the occulter radius values
            input.input_var('impact_axis', 0:0.1:2); % impact parameter values
            input.input_var('velocity', []); % shadow veclocity in km/sec! If empty, will use the header data
            input.input_var('range_v', 3); % include values above/below the nominal velocity in this range
            input.input_var('distance_au', 40);
            input.input_var('wavelength_nm', 500); 
            input.input_var('filename', 'sim_results_generic'); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.coverage)
                input.coverage = obj.getCoverageDegSquare(varargin{:}); 
            end
            
            FSU = sqrt(input.wavelength_nm.*1e-12*input.distance_au.*150e6./2); % physical Fresnel scale, in km!

            r_norm = input.norm_radius./FSU; % translate this size to FSU
            
            density = input.density.*(r_norm./input.radius_axis(1)).^(input.power_law-1); % re-normalize the density using the minimal radius in the simulation
            
            val = density.*input.coverage; 
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function showSimulations(obj, varargin)
           
            input = util.text.InputVars;
            input.input_var('distance', 'kbos'); % which population we want to probe? KBOs? Hills (inner Oort) or Oort?
            input.input_var('parent', []); % parent can be a figure or a panel (default is gcf())
            input.input_var('font_size', 18); % fonts on the axes
            input.scan_vars(varargin{:}); 
            
            if isempty(input.parent)
                input.parent = gcf;
            end
            
            delete(input.parent.Children);
            
            obj.calcSimStats(varargin{:}); % pass through all the relevant inputs
            
            % get shorthands for the edges and names of the variables
            name1 = obj.axes_names{1}; % b
            edges1 = obj.([name1 '_edges']); % b edges
            
            name2 = obj.axes_names{2}; % r
            edges2 = obj.([name2 '_edges']); % r edges
            
            name3 = obj.axes_names{3}; % v
            edges3 = obj.([name3 '_edges']); % v edges
            
            name4 = obj.axes_names{4}; % R
            edges4 = obj.([name4 '_edges']); % R edges
            
            % parameters for the axes grid
            gap = 0.01; 
            margin_x = 0.08;
            margin_y = 0.1;
            Nx = length(edges3)-1;
            Ny = length(edges4)-1;
            width =  (1 - margin_x)/Nx - gap;
            height = (1 - margin_y)/Ny - gap;
            
            ax = {}; 
            
            for ii = 1:Ny % run over R values
                
                for jj = 1:Nx % run over v values
                    
                    ax{ii,jj} = axes('Parent', input.parent, ...
                        'Position', [margin_x + (jj-1)*(width+gap) margin_y + (ii-1)*(height+gap) width height]); 
                    
                    if strcmp(obj.type_result, 'fraction')
                        util.plot.show(obj.results(:,:,jj,ii), 'monochrome', 1,'fancy', 'on', ...
                            'x_values', edges2(1:end-1), 'y_values', edges1(1:end-1)); 
                        
                        colorbar(ax{ii,jj}, 'off'); 
                        axis(ax{ii,jj}, 'normal'); 
                        ax{ii,jj}.YDir = 'normal';
                        title(ax{ii,jj}, ''); 
                    elseif strcmp(obj.type_result, 'snr')
                        
                        contour(ax{ii,jj}, edges2(1:end-1),edges1(1:end-1), obj.results(:,:,jj,ii), [7.5 10 12.5 15 17.5 20], 'k', 'ShowText', 'on'); 
                        
                    end
                    
                    if ii==1 % bottom row of axes
                        
                        xlabel(ax{ii,jj}, name2); 
                        ax{ii,jj}.XTick = [edges2(2) edges2(floor(length(edges2)/2)) edges2(end-3)]; 
                        
                    else
                        ax{ii,jj}.XTick = []; 
                    end
                    
                    if ii==Ny % top row of axes
%                         [n,u] = obj.getFullName(name2);                        
%                         xlabel(ax{ii,jj}, sprintf('%s \n[%s]', n,u));
                    else % non-top row
                        
                        
                        if ~isempty(ax{ii,jj}.YTick)
                            ax{ii,jj}.YTick(end) = [];
                        end
                        
                    end
                    
                    if jj==1 % left-most column of axes
%                         [n,u] = obj.getFullName(name1);
%                         ylabel(ax{ii,jj}, sprintf('%s [%s]', n,u));
                        ylabel(ax{ii,jj}, name1); 
                    else
                        ax{ii,jj}.YTick = [];
                    end
                    
                    if jj==Nx % right-most column of axes
                        
                    else
                        if ~isempty(ax{ii,jj}.XTick)
%                             ax{ii,jj}.XTick(end) = [];
                        end
                    end
                    
                    
                    ax{ii,jj}.FontSize = input.font_size; 
                    
                    util.plot.inner_title(sprintf('%s= %4.2f, %s= %4.2f', name4, edges4(ii), name3, edges3(jj)), ...
                        'Position', 'NorthWest', 'FontSize', input.font_size-2, 'margin', 0.1); 
                    
                end % for jj (v)
                
            end % for ii (R)
                        
        end
        
        function [name, units] = getFullName(obj, letter)
            
            if letter=='R'
                name = 'stellar size'; 
            elseif letter=='r'
                name = 'radius';
            elseif letter=='b'
                name = 'impact parameter'; 
            elseif letter=='v'
                name = 'velocity';
            end
            
            if nargout>1
                if letter=='R'
                    units = 'FSU'; 
                elseif letter=='r'
                    units = 'FSU';
                elseif letter=='b'
                    units = 'FSU'; 
                elseif letter=='v'
                    units = 'FSU/s';
                end 
            end
            
        end
        
        function showDetectionRate(obj, varargin)
            
            import util.text.cs;
            
            input = util.text.InputVars;
            input.input_var('distance', 'kbos'); % which population we want to probe? KBOs? Hills (inner Oort) or Oort?
            input.input_var('parent'); % parent can be figure or panel (default is gcf())
            input.input_var('font_size', 18); % fonts on the axes
            input.scan_vars(varargin{:}); 
            
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
            
            if isempty(input.parent)
                input.parent = gcf;
            end
            
            delete(input.parent.Children);
            
            height = 0.35;
            width = 0.33;
            margin = 0.13;
                        
            ax_R = axes('Parent', input.parent, 'Position', [margin+(margin+width)*0 margin+(margin+height)*0 width height]); 
            ax_r = axes('Parent', input.parent, 'Position', [margin+(margin+width)*1 margin+(margin+height)*0 width height]); 
            ax_b = axes('Parent', input.parent, 'Position', [margin+(margin+width)*0 margin+(margin+height)*1 width height]); 
            ax_v = axes('Parent', input.parent, 'Position', [margin+(margin+width)*1 margin+(margin+height)*1 width height]); 
            
            et = obj.sim_events([obj.sim_events.D]==input.distance); % events total
            ep = obj.sim_events([obj.sim_events.D]==input.distance & [obj.sim_events.passed]); 
            
            %%%%%%%% stellar radius R %%%%%%%%%%%%%%%%%%%%%
            bin_size = max([et.R])/25;
            [N_R,E_R] = histcounts([et.R], 'BinWidth', bin_size); 
            bar(ax_R, E_R(1:end-1)+bin_size/2, N_R); 
            hold(ax_R, 'on'); 
            N_R_passed = histcounts([ep.R], 'BinEdges', E_R); 
            bar(ax_R, E_R(1:end-1)+bin_size/2, N_R_passed); 
            hold(ax_R, 'off'); 
            xlabel(ax_R, 'stellar radius R [FSU]'); 
            ylabel(ax_R, 'number of events'); 
            ax_R.YScale = 'log';
            
            yyaxis(ax_R, 'right'); 
            plot(ax_R, E_R(1:end-1)+bin_size/2, N_R_passed./N_R*100, '-*', 'LineWidth', 2); 
            ax_R.YLim = [0.1 100]; 
            grid(ax_R, 'on'); 
            ytickformat(ax_R, '%d%%'); 
            
            yyaxis(ax_R, 'left'); 
%             legend(ax_R, {'all events', 'passed events', 'percent'}, 'Location', 'NorthEast'); 
            ax_R.FontSize = input.font_size;
            
            %%%%%%%% occulter radius r %%%%%%%%%%%%%%%%%%%%%
            bin_size = 0.25;
            [N_r,E_r] = histcounts([et.r], 'BinWidth', bin_size); 
            bar(ax_r, E_r(1:end-1)+bin_size/2, N_r); 
            hold(ax_r, 'on'); 
            N_r_passed = histcounts([ep.r], 'BinEdges', E_r); 
            bar(ax_r, E_r(1:end-1)+bin_size/2, N_r_passed); 
            hold(ax_r, 'off'); 
            xlabel(ax_r, 'occulter radius r [FSU]'); 
            ylabel(ax_r, 'number of events'); 
            ax_r.YScale = 'log';
            
            yyaxis(ax_r, 'right'); 
            plot(ax_r, E_r(1:end-1)+bin_size/2, N_r_passed./N_r*100, '-*', 'LineWidth', 2); 
            ax_r.YLim = [0.1 100]; 
            grid(ax_r, 'on'); 
            ytickformat(ax_r, '%d%%'); 
            
            yyaxis(ax_r, 'left'); 
%             legend(ax_r, {'all events', 'passed events', 'percent'}, 'Location', 'NorthEast'); 
            ax_r.FontSize = input.font_size;
            
            
            %%%%%%%% impact parameter b %%%%%%%%%%%%%%%%%%%%%
            bin_size = 0.25;
            [N_b,E_b] = histcounts([et.b], 'BinWidth', bin_size); 
            bar(ax_b, E_b(1:end-1)+bin_size/2, N_b); 
            hold(ax_b, 'on'); 
            N_b_passed = histcounts([ep.b], 'BinEdges', E_b); 
            bar(ax_b, E_b(1:end-1)+bin_size/2, N_b_passed); 
            hold(ax_b, 'off'); 
            xlabel(ax_b, 'impact parameter b [FSU]'); 
            ylabel(ax_b, 'number of events'); 
            ax_b.YScale = 'log';
            
            yyaxis(ax_b, 'right'); 
            plot(ax_b, E_b(1:end-1)+bin_size/2, N_b_passed./N_b*100, '-*', 'LineWidth', 2); 
            ax_b.YLim = [0 100]; 
            grid(ax_b, 'on'); 
            ytickformat(ax_b, '%d%%'); 
            
            yyaxis(ax_b, 'left'); 
%             legend(ax_b, {'all events', 'passed events', 'percent'}, 'Location', 'NorthEast'); 
            ax_b.FontSize = input.font_size;
            
            
            %%%%%%%% velocity v %%%%%%%%%%%%%%%%%%%%%
%             bin_size = 2.5;
            bin_size = max([et.v])/15;
            [N_v,E_v] = histcounts([et.v], 'BinWidth', bin_size); 
            bar(ax_v, E_v(1:end-1)+bin_size/2, N_v); 
            hold(ax_v, 'on'); 
            N_v_passed = histcounts([ep.v], 'BinEdges', E_v); 
            bar(ax_v, E_v(1:end-1)+bin_size/2, N_v_passed); 
            hold(ax_v, 'off'); 
            xlabel(ax_v, 'velocity v [FSU/s]'); 
            ylabel(ax_v, 'number of events'); 
            ax_v.YScale = 'log';
            
            yyaxis(ax_v, 'right'); 
            plot(ax_v, E_v(1:end-1)+bin_size/2, N_v_passed./N_v*100, '-*', 'LineWidth', 2); 
            ax_v.YLim = [0 100]; 
            grid(ax_v, 'on'); 
            ytickformat(ax_v, '%d%%'); 
            
            yyaxis(ax_v, 'left'); 
%             legend(ax_v, {'all events', 'passed events', 'percent'}, 'Location', 'NorthEast'); 
            ax_v.FontSize = input.font_size;
            
        end
        
    end    
    
end

