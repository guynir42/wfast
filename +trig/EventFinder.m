classdef EventFinder < handle

    properties(Transient=true)

        gui; % class added later
        bank@occult.ShuffleBank; % randomly picked filter kernels used for matched-filtering (loaded from file when needed)
        bank_small@occult.ShuffleBank; % a smaller filter bank used to weed out only the good stars for full-filtering
        sim_bank@occult.FilterBank; % simulated occultation templates in predefined intervals for checking detection S/N of injected events
        
    end
    
    properties % objects
        
        head@head.Header; % link back to parent object, keep a lot of useful metadata
        cat@head.Catalog; % link back to catalog object, hopefully containing data on the star properties
        
        store@trig.DataStore; 
        
        psd@trig.CorrectPSD; % calculates the Power Spectra Density using Welch's method, then dereddens the flux
        
        var_buf@util.vec.CircularBuffer; % circular buffers for each star/kernel, holding a few variance measurements from past batches
        
        cand@trig.Candidate;
        latest_candidates@trig.Candidate; 
        
    end
    
    properties % inputs/outputs
        
        snr_values; % best S/N of the filtered fluxes, saved for each batch (use this to determine the ideal threshold)
        
        var_values; % keep track of the variance of each star/kernel over the length of the run (only when use_keep_variances=1)
        
        corrected_fluxes; % extended flux after correction (by PSD or by removing linear fit)
        corrected_stds; % standard deviation of the corrected fluxes
        
        filtered_fluxes; % fluxes after matched-filtering with the bank of templates
        bg_filtered_std; % measured background standard deviation after filtering
        
%         star_indices; % a list of the more promising stars that should be turned over to full filtering (e.g., from the pre-filter)
%         star_indices_sim; % same thing, only for stars that have been injected with a simulated event
        
        black_list_stars; % a slowly growing list of stars that display too many events, so that all their events are marked as false
        black_list_batches; % a slowly growing list of batches that display too many events, so that all their events are marked as false
        
    end
    
    properties % switches/controls
        
        threshold = 7.5; % threshold (in units of S/N) for peak of event 
        time_range_thresh = -2.5; % threshold for including area around peak (in continuous time)
        kern_range_thresh = -1; % area threshold (in kernels, discontinuous) NOTE: if negative this will be relative to "threshold"
        star_range_thresh = -1; % area threshold (in stars, discontinuous) NOTE: if this is higher than "threshold" there will be no area around peak
        min_time_spread = 4; % how many frames around the peak to keep, in both directions, even if the flux drops below threshold 
        % (min event duration is 1+2*min_time_spread, unless at the edge) 
        
        max_events = 5; % how many events can we have triggered on the same 2-batch window?
        max_stars = 5; % how many stars can we afford to have triggered at the same time? 
        
        num_hits_black_list = 4; % how many repeated events can we allow before including star or batch in black list
        
        use_psd_correction = 1; % use welch on a flux buffer to correct red noise
        use_std_filtered = 1; % normalize variance of each filtered flux to the average of previous batches (averaging size is set by length of "var_buf")
        
        use_prefilter = 1; % filter all stars on a smaller bank, using a lower threshold, and then vet the survivors only with the full filter bank
        pre_threshold = 5; % threshold for the pre-filter
        
        use_keep_variances = 0; % keep a copy of the variance of each star/kernel for the entire run (this may be a bit heavy on the memory)
        use_conserve_memory = 1; % throw away large image data for disqualified events (reduces size of MAT files, and these data can be restored from file)
        
        use_sim_sporadic = 0; % add simulated events in a few batches randomly spread out in the whole run
        num_sim_events_per_batch = 0.01; % fractional probability to add a simulated event into each new batch. 
        
        use_sim_full = 0; % if true, will run all filter banks on each batch, to measure detectability of different kernels. 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
        
        display_candidate_idx = []; % GUI parameter, which event is to be shown now
        use_display_kept_candidates = 0; % GUI parameter, to show only kept events or all events (default)
        use_show_truth = 0; % if true, will show the field name/coords and for each candidate will reveal if it is simulated
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = EventFinder(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.EventFinder')
                if obj.debug_bit>1, fprintf('EventFinder copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('EventFinder constructor v%4.2f\n', obj.version); end
            
                obj.psd = trig.CorrectPSD;
            
                obj.store = trig.DataStore;
                
                obj.var_buf = util.vec.CircularBuffer; 
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.cand = trig.Candidate.empty;
            obj.latest_candidates = trig.Candidate.empty;
            
            obj.black_list_stars = [];
            obj.black_list_batches = [];
            
            obj.display_candidate_idx = [];
            
            obj.store.reset;
            obj.var_buf.reset;
            obj.psd.reset;
            
            obj.snr_values = [];
            obj.var_values = [];
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            obj.psd.clear;
            obj.store.clear;
            
            obj.corrected_fluxes = [];
            obj.corrected_stds = [];
            
            obj.filtered_fluxes = [];
            obj.bg_filtered_std = [];
            
            obj.latest_candidates = trig.Candidate.empty;
            
        end
        
    end
    
    methods % getters
        
        function val = getTimeThresh(obj)
            
            if obj.time_range_thresh>0
                val = obj.time_range_thresh;
            else
                val = obj.threshold + obj.time_range_thresh;
            end
            
        end
        
        function val = getKernThresh(obj)
            
            if obj.time_range_thresh>0
                val = obj.kern_range_thresh;
            else
                val = obj.threshold + obj.kern_range_thresh;
            end
            
        end
        
        function val = getStarThresh(obj)
            
            if obj.star_range_thresh>0
                val = obj.star_range_thresh;
            else
                val = obj.threshold + obj.star_range_thresh;
            end
            
        end
        
        function val = num_batches_scanned(obj)
            
            val = length(obj.snr_values);
            
        end
        
        function val = num_candidates(obj)
            
            val = length(obj.cand);
            
        end
        
        function val = kept(obj)
            
            val = obj.cand([obj.cand.keep]==1);
            
        end
        
        function val = num_kept(obj)
            
            val = length(obj.kept); 
            
        end
        
        function val = obs_pars_str(obj)
            
            if isempty(obj.head)
                val = '';
            else
                if isempty(obj.head.STARTTIME)
                    date_str = '';
                else
                    date_str = obj.head.STARTTIME(1:10); % get only date, no hours/minutes/seconds
                end
                
                if isempty(obj.head.RA)
                   
                else
                    
                end
                
                val = sprintf('date: %s | RA= %s | DE= %s | airmass= %4.2f | wind= %4.2f km/h | moon: %d%%, %d deg',...
                    date_str, obj.head.RA, obj.head.DEC, obj.head.AIRMASS, obj.head.WIND_SPEED, round(100*obj.head.MOONILL), round(obj.head.MOONDIST));
            end
            
        end
        
        function val = phot_pars_str(obj)
            
            if isempty(obj.head.PHOT_PARS)
                val = '';
            else
                
                val = sprintf('aperture: %s pix | annulus: %s pix | iterations: %d | var_buf: %d | PSD: %d ', ...
                    util.text.print_vec(obj.head.PHOT_PARS.aperture_radius),...
                    util.text.print_vec(obj.head.PHOT_PARS.annulus_radii),...
                    obj.head.PHOT_PARS.iterations, obj.use_var_buf, obj.use_psd_correction);
                
            end
            
        end
                
    end
    
    methods % setters
        
        function set.display_candidate_idx(obj, val)
            
            obj.display_candidate_idx = val;
            if ~isempty(obj.display_candidate_idx) && ~isempty(obj.gui) && obj.gui.check
                obj.cand(obj.display_candidate_idx).show('parent', obj.gui.panel_image);
            end

        end
        
    end
    
    methods % calculations
        
        function input(obj, varargin)
            
            obj.clear;
            
            obj.store.input(varargin{:});
            
            if isempty(obj.bank)
                obj.loadFilterBank;
            end
            
            if obj.use_prefilter && isempty(obj.bank_small)
                obj.loadFilterBankSmall;
            end
            
            if obj.use_sim_full
                obj.loadSimulationBank;
            end
            
            if obj.store.is_done_burn % do not do any more calculations until done with "burn-in" period
                
                if obj.use_psd_correction
                    obj.psd.calcPSD(obj.store.flux_buffer, obj.store.timestamps_buffer, obj.store.length_extended, obj.store.length_background*2); % first off, make the PSD correction for the current batch
                end
                                
                if obj.var_buf.is_empty % check if we haven't used this since the last reset() call
                    N = ceil(obj.store.length_background./obj.store.length_search); % how many batches we want to keep as "background" sample 
                    obj.var_buf.reset(N); 
                end

                if obj.use_sim_full
                    
                else
                     
                    [obj.corrected_fluxes, obj.corrected_stds] = obj.correctFluxes(obj.store.extended_flux, 1);
                    
                    if obj.use_prefilter % use a small filter bank and only send the best stars to the big bank
                        
                        obj.bank_small.input(obj.corrected_fluxes, obj.corrected_stds); 
                        
                        if obj.use_std_filtered
                            
                            V = nanvar(obj.bank_small.fluxes_filtered);

                            obj.var_buf.input(V); 

                            if obj.use_keep_variances
                                obj.var_values = vertcat(obj.var_values, V); 
                            end
                            
                            bg_filtered_std = sqrt(obj.var_buf.mean); 

                        else
                            bg_filtered_std = 1; % do not correct for filtered flux background noise
                        end
                        
                        pre_snr = squeeze(util.stat.max2(obj.bank_small.fluxes_filtered./bg_filtered_std)); % signal to noise peak for each star at any time and any filter kernel
                        
                        star_indices = find(pre_snr>=obj.pre_threshold); 
                        
                    else % do not prefilter, just put all stars into the big filter bank
                        
                        star_indices = 1:size(obj.corrected_fluxes,2); % just list all the stars
                        
                    end
                    
                    if ~isempty(star_indices) % if no stars passed the pre-filter, we will just skip to the next batch
                        
                        obj.filtered_fluxes = obj.bank.input(obj.corrected_fluxes(:,star_indices), obj.corrected_stds(1,star_indices)); % filtered flux for each star/kernel
                        
                        if obj.use_std_filtered % get the noise in the filtered flux background region

                            if obj.use_prefilter % if we used a pre-filter we must also get the background noise for each chosen star, and we don't have var_buf for this!

                                [background_flux, background_std] = obj.correctFluxes(obj.store.background_flux(:,star_indices), star_indices); % run the PSD or linear fit removal, but this time on the background area

                                background_ff = obj.bank.input(background_flux, background_std); % filter these background fluxes also

                                obj.bg_filtered_std = nanstd(background_ff); % this reflects the noise in a few previous batches, on the filtered fluxes, including only the chosen stars

                            else % we can just use the var_buf for this
                                obj.var_buf.input(nanvar(obj.filtered_fluxes)); 
                                obj.bg_filtered_std = sqrt(obj.var_buf.mean); 
                            end

                        else
                            obj.bg_filtered_std = 1; % just assume the filters have unit standard deviation by construction
                        end
                        
                        obj.filtered_fluxes = obj.filtered_fluxes./obj.bg_filtered_std; % normalize the filtered flux by the background std
                        
                        obj.latest_candidates = obj.searchForCandidates(obj.store.extended_flux, obj.corrected_fluxes, obj.filtered_fluxes, star_indices); 
                        
                    end
                    
                    if obj.use_sim_sporadic

                        if rand<obj.num_sim_events_per_batch

                            candidate_stars = false(1,size(obj.corrected_fluxes,2)); % logical vector with 1s where any star has a chance of already having an event

                            if length(star_indices)==size(obj.corrected_fluxes,2) % when not using pre-filter, or if by some coincidence the pre-filter triggered all stars! 
                                candidate_stars([obj.latest_candidates.star_index]) = true; % choose only stars that have a recent candidate
                            else
                                candidate_stars(star_indices) = true; % just mark all the stars chosen by the pre-filter
                            end

                            star_indices = find(~candidate_stars); % all the other stars are good choices for simulations
                            
                            if isscalar(star_indices)
                                star_index_sim = star_indices;
                            else
                                star_index_sim = star_indices(randperm(length(star_indices),1)); % this should now contain exactly one star
                            end

                            f = obj.store.extended_flux; 
                            
                            [f_sim, sim_pars] = obj.addSimulatedEvent(f(:,star_index_sim)); % add the simulated occultation to the raw fluxes
                            
                            f(:,star_index_sim) = f_sim;
                            
                            [f_corr, std_corr] = obj.correctFluxes(f); 
                            
                            f_filt = obj.bank.input(f_corr(:,star_index_sim), std_corr(star_index_sim)); % filter only the star with the simulated event on it
                            
                            % get an estimate for the background of this flux
                            
                            if obj.use_std_filtered % need to correct the filtered fluxes by their measured noise
                                
                                if obj.use_prefilter % need to draw from the background region an calculate the filtered flux noise
                                    [bg_flux, bg_std] = obj.correctFluxes(obj.store.background_flux(:,star_index_sim), star_index_sim); % run correction on the background region
                                    bg_ff = obj.bank.input(bg_flux, bg_std); % filter these background fluxes also
                                    bg_ff_std = nanstd(bg_ff); 
                                else % we have a var_buf for every star, we can just get the noise from that
                                    bg_ff_std = sqrt(obj.var_buf.mean);
                                    bg_ff_std = bg_ff_std(star_index_sim); 
                                end
                                
                                f_filt = f_filt./bg_ff_std; 
                                
                            end
                            
                            new_event_sim = obj.searchForCandidates(f, f_corr, f_filt, star_index_sim, 1); 
                            
                            if ~isempty(new_event_sim)
                                new_event_sim.is_simulated = 1; 
                                new_event_sim.sim_pars = sim_pars; 
                                obj.latest_candidates = vertcat(obj.latest_candidates, new_event_sim); 
                            end
                            
                        end
                        
                    end
                    
                    obj.cand = vertcat(obj.cand, obj.latest_candidates); % store all the candidates that were found in the last batch
                    
                end
                
            end
            
        end
        
        function new_candidates = searchForCandidates(obj, fluxes_raw, fluxes_corrected, fluxes_filtered, star_indices, max_num_events)

            if nargin<6 || isempty(max_num_events)
                max_num_events = obj.max_events;
            end
            
            new_candidates = trig.Candidate.empty;
            
            time_indices = obj.store.search_start_idx:obj.store.search_end_idx;

            for ii = 1:length(max_num_events)

                [mx, idx] = util.stat.maxnd(abs(fluxes_filtered(time_indices,:,:))); % find the maximum point in the filtered fluxes (after correcting for background variance)

                idx(1) = time_indices(idx(1)); % adjust the idx(1) time index to be inside the extended region, instead of in the search region
                
                if length(idx)<3
                    idx(3) = 1;
                end
                
                if ii==1
                    obj.snr_values(end+1) = fluxes_filtered(idx(1),idx(2),idx(3)); % keep track of the highest S/N in this batch (including if it is positive OR negative)
                end

                if mx>=obj.threshold % at least one part of the filtered lightcurves passed the threshold

                    c = obj.makeNewCandidate(fluxes_raw, fluxes_corrected, fluxes_filtered, idx, star_indices); 

                    % add the new event and remove its footprint from the fluxes to search
                    new_candidates = vertcat(new_candidates, c); 

                    fluxes_filtered(c.time_range, :, :) = NaN; % don't look at the same region twice

                end

            end
            
        end
        
        function c = makeNewCandidate(obj, fluxes_raw, fluxes_corrected, fluxes_filtered, idx, star_indices)

            c = trig.Candidate;
            
            c.serial = length(obj.cand) + length(obj.latest_candidates) + 1; 
            
            c.snr = abs(fluxes_filtered(idx(1),idx(2),idx(3))); % note this is positive even for negative filter responses! 
            c.is_positive = fluxes_filtered(idx(1),idx(2),idx(3))>0; 
            c.threshold = obj.threshold; 
            
            c.time_index = idx(1); % time index inside of the extended region
            c.kern_index = idx(2); % which kernel triggered
            c.star_index = star_indices(idx(3)); % which star (out of all stars that passed the burn-in)
            c.star_index_global = obj.store.star_indices(star_indices(idx(3))); % get the star index in the original list of stars before pre-filter and before burn-in

            c.kernel = obj.bank.kernels(:,c.kern_index); 
            c.kern_props = obj.bank.pars(c.kern_index); 
            
            c.time_range = obj.findTimeRange(fluxes_filtered, idx(1), idx(2), idx(3)); % find continuous area (in the extended region) that is above time_range_thresh 
            c.thresh_time = obj.getTimeThresh;
            
            c.kern_extra = find(max(abs(fluxes_filtered(c.time_range, :, idx(3))))>obj.getKernThresh); % find (non-continuous) list of kernels that also show some response at the same time/star
            c.thresh_kern = obj.getKernThresh;
            
            c.star_extra = star_indices(find(max(abs(fluxes_filtered(c.time_range, idx(2), :)))>obj.getStarThresh)); % find (non-continuous) list of stars that also show some response at the same time/kernel
            c.thresh_star = obj.getStarThresh; 
            
            % save the timing data
            c.timestamps = obj.store.extended_timestamps;
            c.juldates = obj.store.extended_juldates;
            c.search_start_idx = obj.store.search_start_idx;
            c.search_end_idx = obj.store.search_end_idx;
            
            % save the filenames and frame numbers
            c.filenames = obj.store.extended_filenames;
            c.frame_numbers = obj.store.extended_frame_num;
            c.frame_index = obj.store.extended_frame_num(c.time_index); 

            % keep a copy of the cutouts, but only for the relevant star!
            c.cutouts = obj.store.cutouts(:,:,:,obj.star_index); 
            
            % store the different types of flux
            c.flux_raw = fluxes_raw(:,c.star_index); 
            c.flux_corrected = fluxes_corrected(:,c.star_index); 
            c.flux_filtered = fluxes_filtered(:,idx(2),idx(3));

            c.flux_raw_all = fluxes_raw; 
            c.flux_corrected_all = fluxes_corrected; 
            
            % store some statistics on this star's raw flux
            c.flux_mean = nanmean(obj.store.background_flux(:,c.star_index)); 
            c.flux_std = nanstd(obj.store.background_flux(:,c.star_index)); 
            
            % save the auxiliary data (like background and offsets)
            c.auxiliary = obj.store.extended_aux; 
            c.aux_names = obj.store.aux_names;
            c.aux_indices = obj.store.aux_indices;

            F = nanmean(obj.store.extended_flux,1);
            DX = c.auxiliary(:,:, c.aux_indices.offsets_x); % the offsets_x for all stars
            dx = DX(:,c.star_index) - util.vec.weighted_average(DX,F,2); % reduced the mean offsets_x
            DY = c.auxiliary(:,:, c.aux_indices.offsets_y); % the offsets_y for all stars
            dy = DY(:,c.star_index) - util.vec.weighted_average(DY,F,2); % reduced the mean offsets_y
            
            c.relative_dx = dx; 
            c.relative_dy = dy; 
            
            % keep the correlations for this star ONLY!
            c.correlations = obj.store.checker.correlations(:,c.star_index,:); % linearize the last dimension so instead of dim3=types and dim4=timescales we have dim3=names on one list
            c.corr_names = obj.store.checker.getCorrNames;
            
            c.corr_indices = struct; 
            for ii = 1:length(c.corr_names)
                c.corr_indices.(c.corr_names{ii}) = ii; 
            end
            
            % keep track of what parameters were used in this analysis
            c.used_psd_corr = obj.use_psd_correction;
            c.used_filt_std = obj.use_std_filtered; 

            % get the observational parameters and the star parameters from astrometry/GAIA
            c.head = obj.head;
            if ~isempty(obj.cat) && obj.cat.success
                c.star_props = obj.cat.data(c.star_index,:); % copy a table row from the catalog
            end

            % disqualify event based on any cut flag
            M = obj.store.checker.cut_flag_matrix; % dim1 is time, dim2 is stars, dim3 is type of cut
            cut = obj.store.checker.cut_names; % cell array of names of each cut in the matrix

            for ii = 1:size(M,3)
                
                if any(M(c.time_range, c.star_index, ii))
                    
                    if strcmpi(cut{ii}(1:4), 'corr')
                        
                        % the name of the cut is corr_<letter>_<timescale number> 
                        % so we need to extract the index of the correlation type and timescale
                        corr_index = obj.store.checker.corr_indices.(cut{ii}(6)); 
                        timescale_index = find(obj.store.checker.corr_timescales==str2double(cut{ii}(8:end))); 
                        cut_values = obj.store.checker.correlations(:,:,corr_index,timescale_index); 
                        thresh = obj.store.checker.thresh_correlation; 
                        
                    elseif isprop(obj.store.checker, cut{ii})
                        cut_values = obj.store.checker.(cut{ii});
                        thresh = obj.store.checker.(['thresh_' cut{ii}]); 
                    end
                    
                    if size(cut_values,2)>1
                        cut_vector = cut_values(:, c.star_index); % extract only the vector pertaining to this star
                    else
                        cut_vector = cut_values;
                    end
                    
                    mx = max(cut_vector(c.time_range)); % the maximum cut value inside the event time
                    
                    str = sprintf('Cut "%s" failed with value %4.2f', cut{ii}, mx);
                    c.addCutString(str); 
                    c.addCutValue(mx); 
                    c.addCutName(cut{ii}); 
                    c.addCutRegion(M(:,c.star_index,ii)); 
                    c.addCutVector(cut_vector); % also keep a vector of the values that didn't pass the cut
                    c.addCutThreshold(thresh); 
                    
                    c.kept = 0; 
                    
                    c.setRunNameDate;
                    
                end
                
            end
            
        end
        
        function time_range = findTimeRange(obj, ff, time_index, kern_index, star_index)
            
            N = size(ff,1); % time length
            
            thresh = obj.getTimeThresh;
            time_range = [];
            
            kernel_min_spread = ceil(nnz(abs(obj.bank.kernels(:,kern_index)>=0.02))/2); % minimal time range cannot be smaller than "active" part of kernel (i.e., above 2% deviation)

            for jj = 0:N % go backward in time

                idx = time_index - jj;

                if idx<1, break; end

%                 if any(abs(ff(idx, :, star_index))>=thresh)
                if abs(ff(idx, kern_index, star_index))>=thresh || jj<=obj.min_time_spread || jj<=kernel_min_spread
                    time_range = [time_range, idx];
                else
                    break;
                end

            end

            time_range = flip(time_range);

            for jj = 1:N % go forward in time

                idx = time_index + jj;

                if idx>N, break; end

%                 if any(abs(ff(idx, :, star_index))>=thresh)
                if abs(ff(idx, kern_index, star_index))>=thresh || jj<=obj.min_time_spread || jj<=kernel_min_spread
                    time_range = [time_range, idx];
                else
                    break;
                end

            end

        end
        
        function [flux_out, std_out] = correctFluxes(obj, fluxes, star_indices) % only supply the star_indices if you are running large matrices like the b/g calculation after a pre-filter
            
            if nargin<3 || isempty(star_indices)
                star_indices = 1:size(fluxes,2);
            end
            
            if obj.use_psd_correction
                flux_out = obj.psd.input(fluxes, 1, star_indices); 
            else
                flux_out = obj.removeLinearFit(fluxes, obj.store.extended_timestamps); 
                flux_out = fillmissing(flux_out, 'linear'); 
            end

            std_out = nanstd(flux_out); 

        end
        
        function flux_rem = removeLinearFit(obj, flux, timestamps)
            
            s = @(x) nansum(x);
            f = flux;
            t = timestamps;
            n = size(f,1);
            
            denom = n.*s(t.^2)-s(t).^2;
            
            a = (n.*s(t.*f)-s(t).*s(f))./denom;
            b = (s(f).*s(t.^2)-s(t).*s(t.*f))./denom;
            
            flux_rem = flux - (a.*t + b); 
            
        end
        
        function [flux, sim_pars] = addSimulatedEvent(obj, flux)
        
            obj.loadSimulationBank; % lazy load the sim_bank
            
            [lc, sim_pars] = obj.sim_bank.randomLC(0.5); % let's set the stellar radius R=0.5 Fresnel units
            
            margin = length(flux) - length(lc);
            
            t = (1:length(flux))';
            
            fr = util.fit.polyfit(t, flux, 'order', 2); 
            
            flux_mean = fr.ym; 
            flux_noise = flux - flux_mean; % separate the noise from the mean flux
            
            if margin>0 % if the lc is similar in size to the search region, this push ensures the peak is in that region
                push = randi(margin); 
                lc = vertcat(lc, ones(margin,1)); 
                lc = circshift(lc, push); 
            elseif margin<0
                error('need to handle this case at some point...'); 
            else
                 % do nothing?
            end
            
            flux_mean = flux_mean.*lc;
            
            flux = flux_mean + flux_noise; 
            
        end
            
        function loadFilterBank(obj)

            f = fullfile(getenv('DATA'), '/WFAST/saved/FilterBankShuffle.mat');
            if exist(f, 'file')
                load(f, 'bank');
                obj.bank = bank;
            else
                error('Cannot load kernels from ShuffleBank object'); 
            end

        end
        
        function loadFilterBankSmall(obj)

            f = fullfile(getenv('DATA'), '/WFAST/saved/FilterBankShuffleSmall.mat');
            if exist(f, 'file')
                load(f, 'bank');
                obj.bank_small = bank;
            else
                error('Cannot load kernels from ShuffleBank object'); 
            end

        end
        
        function loadSimulationBank(obj)

            if isempty(obj.sim_bank)
                obj.sim_bank = occult.FilterBank; % simulations require we produce a FilterBank to make injection events
            end

            if isempty(obj.sim_bank.bank) % simulations require we produce a FilterBank to make injection events
                disp('making a sim_bank!'); 
                obj.sim_bank.makeBank;
            end

%             if isempty(obj.sim_bank.filtered_bank) || numel(obj.sim_bank.filtered_bank)~=numel(obj.bank.kernels)*obj.sim_bank.num_pars
% 
%                 if obj.debug_bit, fprintf('Cross filtering all kernels in ShuffleBank (%d) with all templates in FilterBank (%d)...\n', size(obj.bank.kernels,2), obj.sim_bank.num_pars); end
% 
%                 obj.sim_bank.filtered_bank = util.vec.convolution(obj.bank.kernels, permute(obj.sim_bank.bank-1, [1,6,2,3,4,5])); % making pre-filtered template kernels, these are simply added to filtered data
% 
%             end

        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

