classdef EventFinder < handle
% This class contains the full pipeline for detecting occultation events. 
% 
% To use, make an object and make sure to give it the header and catalog:
% >> finder = trig.EventFinder;
% >> finder.head = header;
% >> finder.cat = catalog;
% Then each time you load some data into an img.Photometry object you can 
% give that object over to the finder for processing the data. 
% >> phot.input(calibrated_cutouts, ...); 
% >> finder.input(phot); 
% After looping over all batches of data, make sure to do:
% >> finder.finishup(); 
% 
% The user-defined parameters of this class are saved as a struct "pars". 
% These parameters are defined in the "reset/clear" methods block, inside
% the resetPars() function. Refer to the inline docs for info on each parameter. 
%
% The finder does several things with the data:
% 1) inputs them into the DataStore object "store". Which in turns feeds the
%    data into a QualityChecker object "checker" inside of itself, which 
%    feeds the data into a StarHours object "hours" inside the "checker". 
%    This branch of the data analysis stores the data in appropriate buffers, 
%    runs different checks and disqualifies some parts of the data, and then 
%    counts the number of useful star hours (each is done in turn by the 
%    abovementioned objects). 
%    This chain of processing can be used independently of the finder object. 
%    If the only requirement is to assess the data quality and quantity, 
%    then you should generate a DataStore object and feed the photometry to 
%    it directly. 
% 2) Data from the store is used to calculate the Power Spectral Density 
%    (PSD) of the data using the CorrectPSD object "psd".
%    This is applied to the fluxes in before searching for occultations. 
%    To skip using the PSD calculation (which is slow) set 
%    pars.use_psd_correction=0. By default we use the PSD correction. 
% 3) Feed the corrected fluxes and auxiliary data into a search algorithm. 
%    This uses a matched-filter search using two filter banks saved as 
%    occult.ShuffleBank objects called "bank" and "bank_small". 
%    Another filter bank called "bank_oort" is used to store wider
%    templates for finding Oort cloud objects (set use_oort=1). 
%
%    Set pars.use_sim = 0 to skip using simulated events (default is 1). 
%    To determine how often simulated events are injected, set
%    "pars.num_sim_events_per_batch" which can be a small fraction (the 
%    default is 0.01, so only a few simulated events per run). 
% 
% The output of the search is two-fold:
% a) The data quality cuts values and the amount of accumulated star hours 
%    are histogrammed and saved as a measure of the amount and quality of 
%    the data over the run. You can get a summary of this information, along
%    with header data and the parameters used in this run by each object by:
%    >> finder.produceSummary; 
%    which produces a RunSummary object. These objects can be collected for 
%    multiple runs and used to, e.g., draw statistics on how much data has 
%    been scanned for occultations. 
% b) The occultation candidates that surpassed the threshold during the event
%    finding branch of the code. 
%    These are returned in a vector of Candidate objects. If it is empty, no
%    events have been found. 
%    Each event has a "kept" property spelling if it passed all data quality
%    checks. Unkept events are saved to give the user a chance to check that
%    the quality cuts are not to harsh, and identify artefacts that somehow
%    pass the cuts (by looking for similar events that are more clearly 
%    marked as artefacts). 
%    The kept events should be classified and saved to an external database. 
% 
    
    
    properties(Transient=true)

        gui@trig.gui.EvFinderGUI; % class added later
        
        bank@occult.ShuffleBank; % randomly picked filter kernels used for matched-filtering (loaded from file when needed)
        bank_small@occult.ShuffleBank; % a smaller filter bank used to weed out only the good stars for full-filtering
        bank_oort@occult.ShuffleBank; % a filter bank for Oort cloud occultations
        
    end
    
    properties % objects
        
        head@head.Header; % link back to parent object, keep a lot of useful metadata
        cat@head.Catalog; % link back to catalog object, hopefully containing data on the star properties
        
        store@trig.DataStore; 
        
        psd@trig.CorrectPSD; % calculates the Power Spectra Density using Welch's method, then dereddens the flux
        
        var_buf@util.vec.CircularBuffer; % circular buffers for each star/kernel, holding a few variance measurements from past batches
        var_buf_oort@util.vec.CircularBuffer; % another buffer for the Oort cloud templates
        
        cand@trig.Candidate;
        latest_candidates@trig.Candidate; 
        
        pars; % a struct with all the user-defined parameters
        
        summary@trig.RunSummary; 
        
    end
    
    properties % inputs/outputs
        
        snr_values; % best S/N of the filtered fluxes, saved for each batch (use this to determine the ideal threshold)
        
        var_values; % keep track of the variance of each star/kernel over the length of the run (only when use_keep_variances=1)
        
        corrected_fluxes; % extended flux after correction (by PSD or by removing linear fit)
        corrected_stds; % standard deviation of the corrected fluxes
        
        filtered_fluxes; % fluxes after matched-filtering with the bank of templates
        bg_filtered_std; % measured background standard deviation after filtering
        
        black_list_stars; % a slowly growing list of stars that display too many events, so that all their events are marked as false
%         black_list_batches; % a slowly growing list of batches that display too many events, so that all their events are marked as false
        
        % simulated events are saved (a short struct describing the event). 
        % The passed events have been detected and the failed events have
        % been rejected by the system. Note that simulated events are placed
        % outside of the data cuts, so most of them are not rejected by them, 
        % but only by the S/N threshold. 
        sim_events; 
    
        total_batches = []; % optionally give the finder the number of files/batches in this run
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
        
        background_ff; % filtered flux background data (either flux values or the variance from the var_buf)
        
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
                obj.var_buf_oort = util.vec.CircularBuffer; 
                
                obj.resetPars; 
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function resetPars(obj)
            
            obj.pars = struct;
            
            obj.pars.threshold = 7.5; % threshold (in units of S/N) for peak of event 
            obj.pars.time_range_thresh = -2.5; % threshold for including area around peak (in continuous time)
            obj.pars.kern_range_thresh = -1; % area threshold (in kernels, discontinuous) NOTE: if negative this will be relative to "threshold"
            obj.pars.star_range_thresh = -1; % area threshold (in stars, discontinuous) NOTE: if this is higher than "threshold" there will be no area around peak
            obj.pars.min_time_spread = 4; % how many frames around the peak to keep, in both directions, even if the flux drops below threshold 
            % (min event duration is 1+2*min_time_spread, unless at the edge) 

            obj.pars.max_events = 5; % we will search iteratively up to this many times for new candidates in each batch
            
            obj.pars.use_psd_correction = 1; % use welch on a flux buffer to correct red noise
            obj.pars.use_std_filtered = 1; % normalize variance of filtered flux of each kernel to the average of previous batches (averaging size is set by length_background in the store)

            obj.pars.use_prefilter = 1; % filter all stars on a smaller bank, using a lower threshold, and then only vet the survivors with the full filter bank
            obj.pars.pre_threshold = 5; % threshold for the pre-filter
            
            obj.pars.use_oort = false; % use the Oort cloud template bank as well
            
            obj.pars.filter_bank_full_filename = '/WFAST/occultations/TemplateBankKBOs.mat'; % filename where the filter bank was taken from (relative to the DATA folder)
            obj.pars.filter_bank_small_filename = '/WFAST/occultations/TemplateBankKBOs_small.mat'; % filename where the smaller filter bank was taken from (relative to the DATA folder)
            obj.pars.filter_bank_oort_filename = '/WFAST/occultations/TemplateBankOort.mat'; % filename where the Oort cloud templates are taken from (relative to the DATA folder)
            
            obj.pars.limit_events_per_batch = 5; % too many events in one batch will mark all events as black listed! 
            obj.pars.limit_events_per_star = 5; % too many events on the same star will mark all events on that star as black listed! 

            obj.pars.use_keep_variances = 1; % keep a copy of the variance of each star/kernel for the entire run (this may be a bit heavy on the memory)

            obj.pars.use_sim = 1; % add simulated events in a few batches randomly spread out in the whole run
            obj.pars.num_sim_events_per_batch = 0.01; % fractional probability to add a simulated event into each new batch. 
            obj.pars.use_keep_simulated = true; % if false, the simulated events would not be kept with the list of detected candidates (for running massive amount of simualtions)
            obj.pars.sim_max_R = 3; % maximum value of stellar size for simulated events
            
            obj.reset;
            
        end
        
        function reset(obj)
            
            obj.cand = trig.Candidate.empty;
            obj.latest_candidates = trig.Candidate.empty;
            
            obj.black_list_stars = [];
            
            obj.display_candidate_idx = [];
            
            obj.store.reset;
            obj.var_buf.reset;
            obj.psd.reset;
            
            obj.snr_values = [];
            obj.var_values = [];
            
            obj.sim_events = [];
            
            obj.total_batches = [];
            
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
        
        function val = getTimeThresh(obj) % get the threshold for timestamps adjacent (continuous) to the peak (translate the relative pars.time_range_thresh to absolute threshold)
            
            if obj.pars.time_range_thresh>0
                val = obj.pars.time_range_thresh;
            else
                val = obj.pars.threshold + obj.pars.time_range_thresh;
            end
            
        end
        
        function val = getKernThresh(obj) % get the threshold for other kernels during the occultation time (translate the relative pars.kern_range_thresh to absolute threshold)
            
            if obj.pars.time_range_thresh>0
                val = obj.pars.kern_range_thresh;
            else
                val = obj.pars.threshold + obj.pars.kern_range_thresh;
            end
            
        end
        
        function val = getStarThresh(obj) % get the threshold for other stars during the occultation time (translate the relative pars.star_range_thresh to absolute threshold)
            
            if obj.pars.star_range_thresh>0
                val = obj.pars.star_range_thresh;
            else
                val = obj.pars.threshold + obj.pars.star_range_thresh;
            end
            
        end
        
        function val = batch_counter(obj) % how many batches were processed and searched for occultations (not in the burn-in period)
            
            val = length(obj.snr_values);
            
        end
        
        function val = num_candidates(obj) % how many candidates did we find in this run
            
            val = length(obj.cand);
            
        end
        
        function val = kept(obj) % filter only the kept candidates
            
            val = obj.cand([obj.cand.kept]==1);
            
        end
        
        function val = num_kept(obj) % how many candidates were kept
            
            val = length(obj.kept); 
            
        end
        
        function val = obs_pars_str(obj) % summarize the observational parameters
            
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
        
        function val = phot_pars_str(obj) % summarize the photometric parameters used 
            
            if isempty(obj.head.PHOT_PARS)
                val = '';
            else
                
                val = sprintf('aperture: %s pix | annulus: %s pix | iterations: %d | var_buf: %d | PSD: %d ', ...
                    util.text.print_vec(obj.head.PHOT_PARS.aperture_radius),...
                    util.text.print_vec(obj.head.PHOT_PARS.annulus_radii),...
                    obj.head.PHOT_PARS.iterations, obj.pars.use_var_buf, obj.pars.use_psd_correction);
                
            end
            
        end
                
        function val = runtime_hours(obj)
            
            val = obj.store.checker.hours.runtime/3600; 
            
        end
        
        function val = runtime_batches_str(obj)
            
            val = sprintf('%5.1fh, %d batches', round(obj.runtime_hours,1), obj.batch_counter); 
            
        end
        
        function val = star_hours(obj)
            
            val = util.stat.sum2(obj.store.checker.hours.histogram)/3600;
            
        end
        
        function val = star_hours_total(obj)
            
            val = util.stat.sum2(obj.store.checker.hours.histogram_with_losses)/3600;
            
        end
        
        function val = star_hours_str(obj)
            
            val = sprintf('star hours= %d/%d', round(obj.star_hours), round(obj.star_hours_total)); 
            
        end
        
        function val = num_stars_str(obj)
            
            val = sprintf('stars= %d/%d', length(obj.store.star_indices), length(obj.store.star_snr)); 
            
        end
        
    end
    
    methods % setters
        
        function set.head(obj, val) % setting the header of this object also cascades down to the store and its sub-objects
            
            obj.head = val;
            
            obj.store.head = val;
            
        end
        
    end
    
    methods % calculations
        
        function input(obj, varargin) % provide the photometric data as varargin pairs, or as a util.text.InputVar object, or as an img.Photometry object
            
            obj.clear;
            
            obj.store.input(varargin{:}); % the store does the actual parsing and organizing of data into buffers
            
            obj.pars.analysis_time = util.text.time2str('now'); % keep a record of when the analysis was done
            
            if isempty(obj.bank) % lazy load the filter bank from file
                obj.loadFilterBank; 
            end
            
            if obj.pars.use_prefilter && isempty(obj.bank_small) % lazy load the small filter bank from file
                obj.loadFilterBankSmall; 
            end
            
            if obj.pars.use_oort && isempty(obj.bank_oort) % lazy load the small filter bank from file
                obj.loadFilterBankOort; 
            end
            
            if isempty(obj.head) % must have a header to uniquely identify candidates! 
                error('Cannot find events without a header object...'); 
            end
            
            if isempty(obj.head.run_identifier) % must have a run identifier to tag each event
                d = fileparts(obj.store.filename_buffer{1}); 
                [d, run_name] = fileparts(d); 
                [d, run_date] = fileparts(d); 
                obj.head.run_identifier = fullfile(run_date, run_name); 
            end
            
            if obj.store.is_done_burn % do not do any more calculations until done with "burn-in" period
                
                if obj.pars.use_psd_correction
                    obj.psd.calcPSD(obj.store.flux_buffer, obj.store.timestamps_buffer, obj.store.pars.length_extended, obj.store.pars.length_background*2); % first off, make the PSD correction for the current batch
                end
                                
                if obj.var_buf.is_empty % check if we haven't used this since the last reset() call
                    N = ceil(obj.store.pars.length_background./obj.store.pars.length_search); % how many batches we want to keep as "background" sample 
                    obj.var_buf.reset(N); 
                end
 
                [obj.corrected_fluxes, obj.corrected_stds] = obj.correctFluxes(obj.store.extended_flux, 1); % runs either PSD correction or simple removal of linear fit to each lightcurve

                if obj.pars.use_prefilter % use a small filter bank and only send the best stars to the big bank

                    obj.bank_small.input(obj.corrected_fluxes, obj.corrected_stds); % matched-filtering produces a 3D matrix inside bank_small.fluxes_filtered

                    if obj.pars.use_std_filtered % if true, will apply a correction based on the STD of each filtered flux individually

                        V = nanvar(obj.bank_small.fluxes_filtered);

                        obj.var_buf.input(V); % variance buffer saves the variance of several batches in the past

                        if obj.pars.use_keep_variances % keep a log of the variance of each star and each filter (this accumulates to a big chunk of memory)
                            obj.var_values = vertcat(obj.var_values, V); 
                        end

                        bg_filtered_std = sqrt(obj.var_buf.mean); % get the STD of the filtered flux (from previous batches - i.e., background region)

                    else
                        bg_filtered_std = 1; % do not correct for filtered flux background noise
                    end

                    pre_snr = squeeze(util.stat.max2(obj.bank_small.fluxes_filtered./bg_filtered_std)); % signal to noise peak for each star at any time and any filter kernel

                    star_indices = find(pre_snr>=obj.pars.pre_threshold); % only stars that have shown promising response to the small filter bank are on this list 

                else % do not prefilter, just put all stars into the big filter bank

                    star_indices = 1:size(obj.corrected_fluxes,2); % just list all the stars

                end

                best_snr = []; % need to keep track what is the best S/N for this batch, regardless of which filter and regardless of if there were any candidates detected. 
                
                if ~isempty(star_indices) % if no stars passed the pre-filter, we will just skip to the next batch (if no pre-filter is used, all stars will pass)

                    obj.filtered_fluxes = obj.bank.input(obj.corrected_fluxes(:,star_indices), obj.corrected_stds(1,star_indices)); % filtered flux for each star/kernel

                    if obj.pars.use_std_filtered % get the noise in the filtered flux background region

                        if obj.pars.use_prefilter % if we used a pre-filter we must also get the background noise for each chosen star, and we don't have var_buf for this!

                            [background_flux, background_std] = obj.correctFluxes(obj.store.background_flux(:,star_indices), star_indices); % run the PSD or linear fit removal, but this time on the background area

                            obj.background_ff = obj.bank.input(background_flux, background_std); % filter these background fluxes also

                            obj.bg_filtered_std = nanstd(obj.background_ff); % this reflects the noise in a few previous batches, on the filtered fluxes, including only the chosen stars

                        else % not using pre-filter, means we have a var_buf with background variance for each star/filter
                            obj.var_buf.input(nanvar(obj.filtered_fluxes)); 
                            obj.background_ff = obj.var_buf.data_ordered; 
                            obj.bg_filtered_std = sqrt(nanmean(obj.background_ff)); 
                        end

                    else
                        obj.bg_filtered_std = 1; % just assume the filters have unit standard deviation by construction
                    end

                    obj.filtered_fluxes = obj.filtered_fluxes./obj.bg_filtered_std; % normalize the filtered flux by the background std

                    [obj.latest_candidates, best_snr] = obj.searchForCandidates(obj.store.extended_flux, obj.corrected_fluxes, obj.filtered_fluxes, star_indices); % loop over the normalized filtered flux and find multiple events

                    if length(obj.latest_candidates)>=obj.pars.limit_events_per_batch % this batch has too many events, need to mark them as black-listed! 
                        
                        for ii = 1:length(obj.latest_candidates)
                            obj.latest_candidates(ii).notes{end+1} = sprintf('Batch is black listed with %d events', length(obj.latest_candidates)); 
                            obj.latest_candidates(ii).kept = 0;
                        end
                        
                        obj.store.saveHours(1); % the parameter 1 is used to mark this as a bad batch...
                        
                    else
                        obj.store.saveHours; % with no arguments it just counts the hours in this batch as good times
                    end
                    
                    if obj.pars.use_sim % simulated events are injected into the data and treated like real events
                        
                        star_indices = obj.findStarsForSim(star_indices); % get stars that do not include any possible real candidates
                        
                        for ii = 1:obj.getNumSimulations % this number can be more than one, or less (then we randomly decide if to include a simulated event in this batch)
                            obj.simulate(star_indices); % each call to this function tries to add a single simulated event to a random star from the list                            
                        end
                        
                    end
                    
                    obj.cand = vertcat(obj.cand, obj.latest_candidates); % store all the candidates that were found in the last batch
                    
                end % if no stars passed the pre-filter, we will just skip to the next batch

                if obj.pars.use_oort % run the fluxes through the Oort cloud template bank as well 
                    
                    star_indices = 1:size(obj.corrected_fluxes,2); % just list all the stars
                    
                    oort_filtered_fluxes = obj.bank.input(obj.corrected_fluxes(:,star_indices), obj.corrected_stds(1,star_indices)); % filtered flux for each star/kernel

                    if obj.pars.use_std_filtered % get the noise in the filtered flux background region

                        obj.var_buf_oort.input(nanvar(oort_filtered_fluxes)); 
                        oort_background_ff = obj.var_buf_oort.data_ordered; 
                        oort_bg_filtered_std = sqrt(nanmean(oort_background_ff)); 

                    else
                        oort_bg_filtered_std = 1; % just assume the filters have unit standard deviation by construction
                    end

                    oort_filtered_fluxes = oort_filtered_fluxes./oort_bg_filtered_std; % normalize the filtered flux by the background std

                    oort_latest_candidates = obj.searchForCandidates(obj.store.extended_flux, obj.corrected_fluxes, oort_filtered_fluxes, star_indices); % loop over the normalized filtered flux and find multiple events

                    for ii = 1:length(oort_latest_candidates)
                        oort_latest_candidates(ii).oort_filter = true; 
                    end
                    
                    obj.cand = vertcat(obj.cand, oort_latest_candidates);
                    
                end
                
                if ~isempty(obj.latest_candidates) % if candidates are found, check them for highest S/N
                    obj.snr_values(end+1) = max(abs([obj.latest_candidates.snr])); 
                elseif ~isempty(best_snr) % we managed to find a good S/N value from the large filter bank
                    obj.snr_values(end+1) = best_snr; % no candidates, instead just return the best S/N found when searching for candidates
                elseif obj.pars.use_prefilter && ~isempty(pre_snr) % we used a pre-filter and didn't get any stars that passed
                    obj.snr_values(end+1) = max(pre_snr); % take the best S/N from the smaller filter bank
                else
                    error('This shouldn''t happen!'); 
                end
                
            end % do not do any more calculations until done with "burn-in" period
            
        end
        
        function finishup(obj) % call this at the end of the run 
            
            if obj.pars.limit_events_per_star % check if any stars need to be black listed
                
                real_events = ~[obj.cand.is_simulated]; % don't include simulated events in the black list! 

                N = histcounts([obj.cand(real_events).star_index], 'BinEdges', 1:size(obj.store.extended_flux,2)+1);

                idx = N>=obj.pars.limit_events_per_star;

                obj.black_list_stars = find(idx); 

                obj.store.checker.hours.removeStars(idx); 

                for ii = 1:length(obj.cand)

                    star = obj.cand(ii).star_index; 

                    if ismember(star, obj.black_list_stars)
                        str = sprintf('Star %d blacklisted with %d events', star, N(star)); 
                        str_idx = strcmp(str, obj.cand(ii).notes);
                        if isempty(str_idx)
                            obj.cand(ii).notes{end+1} = str; 
                        end
                        obj.cand(ii).kept = 0; 
                    end

                end

            end % check if any stars need to be black listed
            
            obj.summary = obj.produceSummary; 
            
        end
        
        function s = produceSummary(obj) % makes a RunSummary with all relevant data on the quality and quantity of star hours, etc
            
            % add varargin? 
            
            s = trig.RunSummary; % generate a new object
            s.head = util.oop.full_copy(obj.head); % I want to make sure this summary is distinct from the original finder
            
            % load the content of the finder
            s.finder_pars = obj.pars;
            s.snr_values = obj.snr_values;
            s.total_batches = obj.total_batches; 
            s.sim_events = obj.sim_events;
            s.black_list_stars = obj.black_list_stars;
            
            % load the content of the store
            s.store_pars = obj.store.pars;
            s.good_stars = obj.store.star_indices;
            
            % load the content of the checker
            s.checker_pars = obj.store.checker.pars;
            s.cut_names = obj.store.checker.cut_names;
            s.cut_indices = obj.store.checker.cut_indices; 
            s.cut_thresholds = obj.store.checker.cut_thresholds;
            s.cut_two_sided = obj.store.checker.cut_two_sided;
            s.cut_histograms = permute(nansum(obj.store.checker.histograms,2), [1,3,2]);
            s.cut_bin_edges = obj.store.checker.hist_edges;
            
            % load the content of the star hours
            if ~ismember('FresnelSize', obj.cat.data.Properties.VariableNames)
                obj.cat.addStellarSizes;
            end
            
            s.inputHours(obj.store.checker.hours, obj.cat.data.FresnelSize(obj.store.star_indices));
            
%             s.snr_bin_edges = obj.store.checker.hours.snr_bin_edges;
%             s.star_seconds = permute(nansum(obj.store.checker.hours.histogram,2),[1,3,2]);
%             s.star_seconds_with_losses = permute(nansum(obj.store.checker.hours.histogram_with_losses,2),[1,3,2]);
%             s.losses_exclusive = permute(nansum(obj.store.checker.hours.losses_exclusive, 2), [1,3,2]);
%             s.losses_inclusive = permute(nansum(obj.store.checker.hours.losses_inclusive, 2), [1,3,2]);
%             s.losses_bad_stars = permute(nansum(obj.store.checker.hours.losses_bad_stars, 2), [1,3,2]);
%             s.runtime = obj.store.checker.hours.runtime; 
            
            obj.summary = s; 
            
        end
        
    end
    
    methods(Hidden=true) % internal methods
        
        function [new_candidates, best_snr] = searchForCandidates(obj, fluxes_raw, fluxes_corrected, fluxes_filtered, star_indices, max_num_events) % loop over filtered fluxes, find above-threshold peaks, remove that area, then repeat

            if nargin<6 || isempty(max_num_events)
                max_num_events = obj.pars.max_events;
            end
            
            new_candidates = trig.Candidate.empty;
            
            time_indices = obj.store.search_start_idx:obj.store.search_end_idx;

            for ii = 1:length(max_num_events)

                [mx, idx] = util.stat.maxnd(abs(fluxes_filtered(time_indices,:,:))); % find the maximum point in the filtered fluxes (after correcting for background variance)

                idx(1) = time_indices(idx(1)); % adjust the idx(1) time index to be inside the extended region, instead of in the search region
                
                if length(idx)<3
                    idx(3) = 1;
                end
                
                if ii==1, best_snr = mx; end % return the best S/N from this batch
                
                if mx>=obj.pars.threshold % at least one part of the filtered lightcurves passed the threshold

                    c = obj.makeNewCandidate(fluxes_raw, fluxes_corrected, fluxes_filtered, idx, star_indices); 

                    c.serial = length(obj.cand) + length(obj.latest_candidates) + 1; % new event gets a serial number
                    new_candidates = vertcat(new_candidates, c); % add this event to the list of new candidates
                    
                    % remove the new event's footprint from the fluxes to search
                    fluxes_filtered(c.time_range, :, :) = NaN; % don't look at the same region twice

                    % do we have to ALSO reduce this region from the star hours?
                    
                end

            end
            
        end
        
        function c = makeNewCandidate(obj, fluxes_raw, fluxes_corrected, fluxes_filtered, idx, star_indices) % generate a Candidate object and fill it with data

            c = trig.Candidate;
            
            c.snr = abs(fluxes_filtered(idx(1),idx(2),idx(3))); % note this is positive even for negative filter responses! 
            c.is_positive = fluxes_filtered(idx(1),idx(2),idx(3))>0; % negative events are this where the filter response is large negative (e.g., flares)
                        
            c.time_index = idx(1); % time index inside of the extended region
            c.kern_index = idx(2); % which kernel triggered
            c.star_index = star_indices(idx(3)); % which star (out of all stars that passed the burn-in)
            c.star_index_global = obj.store.star_indices(star_indices(idx(3))); % get the star index in the original list of stars before pre-filter and before burn-in

            c.kernel = obj.bank.kernels(:,c.kern_index); % the lightcurve of the kernel used (zero centered and normalized)
            c.kern_props = obj.bank.pars(c.kern_index); % properties of the occultation that would generate this kernel
            c.kern_props.inverted = ~c.is_positive; % if the filtered flux was negative, it means we had to invert the kernel
            
            % find continuous area (in the extended region) that is above time_range_thresh 
            c.time_range = obj.findTimeRange(fluxes_filtered, idx(1), idx(2), idx(3)); 
            
            % find (non-continuous) list of kernels that also show some response at the same time/star
            c.kern_extra = find(max(abs(fluxes_filtered(c.time_range, :, idx(3))))>obj.getKernThresh); 
            
            % find (non-continuous) list of stars that also show some response at the same time/kernel
            c.star_extra = star_indices(find(max(abs(fluxes_filtered(c.time_range, idx(2), :)))>obj.getStarThresh));
            
            % save the timing data
            c.timestamps = obj.store.extended_timestamps; % internal camera clock timestamps
            c.juldates = obj.store.extended_juldates; % translated to Julian date
            c.search_start_idx = obj.store.search_start_idx; % index in the extended region where the search region starts
            c.search_end_idx = obj.store.search_end_idx; % index in the extended region where the search region ends
            
            % save the filenames and frame numbers
            c.filenames = obj.store.extended_filenames; % cell array with full path filenames of each frame (mostly there are only 2 unique filenames)
            c.frame_numbers = obj.store.extended_frame_num; % frame number inside the files named above
            c.frame_index = obj.store.extended_frame_num(c.time_index); % the frame index of the event peak, inside the specific file where the event occured

            % keep a copy of the cutouts, but only for the relevant star!
            c.cutouts = obj.store.cutouts(:,:,:,c.star_index); 
            
            % store the different types of flux
            c.flux_raw = fluxes_raw(:,c.star_index); 
            c.flux_corrected = fluxes_corrected(:,c.star_index); 
            c.flux_filtered = fluxes_filtered(:,idx(2),idx(3));

            % store data for all stars as reference
            c.flux_raw_all = fluxes_raw; 
            c.flux_corrected_all = fluxes_corrected; 
            c.auxiliary_all = obj.store.extended_aux; 
            
            idx = find(star_indices==c.star_index); % the index of the star inside the list of stars that passed the pre-filter
            c.filtered_flux_past_values = obj.background_ff(:,idx); % because we only have this background_ff for the subset of stars
            c.flux_buffer = obj.store.flux_buffer(:,c.star_index); % flux history for this star
            c.timestamps_buffer = obj.store.timestamps_buffer; % timestamps for that duration
            c.psd = obj.psd.power_spectrum(:,c.star_index); % the Power Spectral Density (PSD) for this star
            c.freq_psd = obj.psd.freq; % frequency axis for the PSD
            
            % store some statistics on this star's raw flux
            A = obj.store.background_aux(:,c.star_index, obj.store.aux_indices.areas); 
            B = obj.store.background_aux(:,c.star_index, obj.store.aux_indices.backgrounds); 
            c.flux_mean = nanmean(obj.store.background_flux(:,c.star_index)-A.*B); 
            c.flux_std = nanstd(obj.store.background_flux(:,c.star_index)); 
            
            % save the auxiliary data (like background and offsets)
            c.auxiliary = permute(obj.store.extended_aux(:,c.star_index,:), [1,3,2]); 
            c.aux_names = obj.store.aux_names;
            c.aux_indices = obj.store.aux_indices;

            DX = c.auxiliary(:,c.aux_indices.offsets_x); % the offsets_x for this star
            dx = DX - obj.store.checker.mean_x; % reduce the mean offsets_x
            DY = c.auxiliary(:,c.aux_indices.offsets_y); % the offsets_y for this star
            dy = DY - obj.store.checker.mean_y; % reduce the mean offsets_y
            
            % store the relative offsets for this star
            c.relative_dx = dx; 
            c.relative_dy = dy; 
            
            % get the observational parameters and the star parameters from astrometry/GAIA
            c.head = obj.head;
            if ~isempty(obj.cat) && obj.cat.success
                c.star_props = obj.cat.data(c.star_index,:); % copy a table row from the catalog
            end

            % save the parameters used by the finder, store and quality-checker
            c.finder_pars = obj.pars; 
            c.store_pars = obj.store.pars;
            c.checker_pars = obj.store.checker.pars;
            % add StarHours parameters?? 
            
            % find the time-frames, the kernels and the stars that also triggered at the same time
            c.finder_pars.time_range_thresh = obj.getTimeThresh;
            c.finder_pars.kern_range_thresh = obj.getKernThresh;
            c.finder_pars.star_range_thresh = obj.getStarThresh; 
            c.finder_pars.total_batches = obj.total_batches; 
            c.batch_number = obj.batch_counter; 

            % store the name and date of the run (from the filenames/header)
            c.run_identifier = obj.head.run_identifier; 

            c.checkIfPeakIsIncluded; % set kept=0 for candidates that have a higher peak outside the search region than inside it (duplicate events)
            
            % disqualify event based on any cut flag
            c.cut_matrix = permute(obj.store.checker.cut_values_matrix(:,c.star_index,:), [1,3,2]); % dim1 is time, dim2 is type of cut
            c.cut_names = obj.store.checker.cut_names; % cell array of names of each cut in the matrix
            c.cut_indices = obj.store.checker.cut_indices; % struct with the indices of each cut in each field
            c.cut_thresh = obj.store.checker.cut_thresholds; 
            c.cut_two_sided = obj.store.checker.cut_two_sided; 
            
            for ii = 1:length(c.cut_names) % go oveer each cut and check if it applies to this event
            
                if obj.store.checker.cut_flag_matrix(c.time_index, c.star_index, ii) % the event occurs on one of the flagged frames/stars
                    
                    if obj.store.checker.pars.use_dilate % must look in the entire dilated area to look for the biggest value of this cut
                        indices = (-obj.store.checker.pars.dilate_region:obj.store.checker.pars.dilate_region) + c.time_index; % region around peak that can be excluded by this cut
                        if c.cut_two_sided(ii)
                            [~, idx] = max(abs(c.cut_matrix(indices,ii))); % find the absolute max, then...
                            mx = max(c.cut_matrix(indices(idx),ii)); % get the value with sign
                        else
                            mx = max(c.cut_matrix(indices,ii)); % maximum value of this cut in this region
                        end
                        
                    else
                        mx = c.cut_matrix(c.time_index,ii); % no dilation, just keep the cut value at peak
                    end
                    
                    if isfield(obj.store.checker.pars, ['thresh_' c.cut_names{ii}]) || strcmp(c.cut_names{ii}(1:4), 'corr') % this is a "threshold" type test
                        c.cut_string{end+1,1} = sprintf('Cut "%s" failed with value %4.2f', c.cut_names{ii}, mx);
                    else % this is a pass/fail type cut
                        c.cut_string{end+1,1} = sprintf('Cut "%s" failed', c.cut_names{ii});
                    end
                    
                    c.cut_value(end+1) = mx; % store the value of the cut that this event overlapped
                    c.cut_hits(end+1) = ii; % mark which cut number has triggered (could be more than one!)
                    c.kept = 0; % mark the event as not real
                    
                end % the event occurs on one of the flagged frames/stars
                
            end % go oveer each cut and check if it applies to this event
            
        end
        
        function time_range = findTimeRange(obj, ff, time_index, kern_index, star_index) % find all continuous time frames around peak that have surpassed a (lower) threshold, or just mark a minimal event region
            
            N = size(ff,1); % time length
            
            thresh = obj.getTimeThresh;
            time_range = [];
            
            kernel_min_spread = ceil(nnz(abs(obj.bank.kernels(:,kern_index)>=0.02))/2); % minimal time range cannot be smaller than "active" part of kernel (i.e., above 2% deviation)

            for jj = 0:N % go backward in time

                idx = time_index - jj;

                if idx<1, break; end

%                 if any(abs(ff(idx, :, star_index))>=thresh)
                if abs(ff(idx, kern_index, star_index))>=thresh || jj<=obj.pars.min_time_spread || jj<=kernel_min_spread
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
                if abs(ff(idx, kern_index, star_index))>=thresh || jj<=obj.pars.min_time_spread || jj<=kernel_min_spread
                    time_range = [time_range, idx];
                else
                    break;
                end

            end

        end
        
        function [flux_out, std_out] = correctFluxes(obj, fluxes, star_indices) % only supply the star_indices if you are running large matrices like the b/g calculation after a pre-filter
            
            if nargin<3 || isempty(star_indices) % by default run the correction for all stars
                star_indices = 1:size(fluxes,2);
            end
            
            if obj.pars.use_psd_correction
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
        
        function loadFilterBank(obj)

            f = fullfile(getenv('DATA'), obj.pars.filter_bank_full_filename);
            if exist(f, 'file')
                load(f, 'bank');
                obj.bank = bank;
            else
                error('Cannot load kernels from ShuffleBank object'); 
            end

        end
        
        function loadFilterBankSmall(obj)

            f = fullfile(getenv('DATA'), obj.pars.filter_bank_small_filename);
            if exist(f, 'file')
                load(f, 'bank');
                obj.bank_small = bank;
            else
                error('Cannot load kernels from ShuffleBank object'); 
            end

        end
        
        function loadFilterBankOort(obj)

            f = fullfile(getenv('DATA'), obj.pars.filter_bank_oort_filename);
            if exist(f, 'file')
                load(f, 'bank');
                obj.bank_oort = bank;
            else
                error('Cannot load kernels from ShuffleBank object'); 
            end

        end
        
    end
    
    methods (Hidden=true) % simulations
        
        function star_indices = findStarsForSim(obj, star_indices) % take a list of stars from the prefilter (or all stars) and get a list of good stars for a simulated event
           
            stars_not_for_sim = false(1,size(obj.corrected_fluxes,2)); % logical vector with 1s where any star has a chance of already having an event

            if length(star_indices)==size(obj.corrected_fluxes,2) % when not using pre-filter, or if by some coincidence the pre-filter triggered all stars! 
                stars_not_for_sim([obj.latest_candidates.star_index]) = true; % choose only stars that have a recent candidate
            else
                stars_not_for_sim(star_indices) = true; % just mark all the stars chosen by the pre-filter
            end

            % make sure to take only stars with good frames! 
            B = obj.store.checker.bad_times;                            
            stars_not_for_sim = stars_not_for_sim | all(B(obj.store.search_start_idx:obj.store.search_end_idx,:)); % also disqualify stars that are completely ruled out by the checker

            star_indices = find(~stars_not_for_sim); % all the other stars are good choices for simulations
            
        end
        
        function val = getNumSimulations(obj) % if num_sim_events_per_batch is < 1, randomly choose 0 or 1 based on uniform distribution and the given fraction. If >1, use that number in Poisson distribution to choose a number of events
            
            if obj.pars.num_sim_events_per_batch<1
                val = double(rand<obj.pars.num_sim_events_per_batch);
            else
                val = poissrnd(obj.pars.num_sim_events_per_batch); 
            end
            
        end
        
        function simulate(obj, star_indices) % choose a star from the list of star_indices and add a randomly chosen occultation to it, then search for that event
            
            if isempty(star_indices) % no stars, no simulation! 
                return;
            elseif isscalar(star_indices) % single star is given, just use it
                star_index_sim = star_indices;
            else
                star_index_sim = star_indices(randperm(length(star_indices),1)); % randomly choose a single star from the list
            end

%             f = obj.store.extended_flux; 

            [f_sim, sim_pars] = obj.addSimulatedEvent(star_index_sim); % add the simulated occultation to the raw fluxes

%             f(:,star_index_sim) = f_sim;

            [f_corr, std_corr] = obj.correctFluxes(f_sim); % PSD or linear fit correction (on top of the simulated event!)

            f_filt = obj.bank.input(f_corr(:,star_index_sim), std_corr(star_index_sim)); % filter only the star with the simulated event on it

            % get an estimate for the background of this flux:
            if obj.pars.use_std_filtered % need to correct the filtered fluxes by their measured noise

                if obj.pars.use_prefilter % need to draw from the background region and calculate the filtered flux noise
                    [bg_flux, bg_std] = obj.correctFluxes(obj.store.background_flux(:,star_index_sim), star_index_sim); % run correction on the background region
                    obj.background_ff = obj.bank.input(bg_flux, bg_std); % filter these background fluxes also
                    bg_ff_std = nanstd(obj.background_ff); 
                else % we have a var_buf for every star, we can just get the noise from that
                    obj.background_ff = obj.var_buf.data_ordered; 
                    bg_ff_std = sqrt(nanmean(obj.background_ff));
                    bg_ff_std = bg_ff_std(star_index_sim); 
                end

                f_filt = f_filt./bg_ff_std; 

            end % need to correct the filtered fluxes by their measured noise

            new_event_sim = obj.searchForCandidates(f, f_corr, f_filt, star_index_sim, 1); % try to find a single event on this single star

            if ~isempty(new_event_sim) % we recovered the injected event! 
                
                sim_pars.detect_snr = new_event_sim.snr; % keep track of the detection S/N for this event
                sim_pars.passed = true; 
                
                new_event_sim.is_simulated = 1; 
                new_event_sim.sim_pars = sim_pars; 
                
                if obj.pars.use_keep_simulated
                    obj.latest_candidates = vertcat(obj.latest_candidates, new_event_sim); % add the simulated events to the list of regular events
                end
                
                % add this parameter struct to the list of passed events
                if isempty(obj.sim_events)
                    obj.sim_events = sim_pars;
                else
                    obj.sim_events(end+1) = sim_pars;
                end
                
            else % no events were recovered
                
                sim_pars.detect_snr = 0; % event was not detected, so S/N is zero
                sim_pars.passed = false; 
                
                % add this parameter struct to the list of failed events
                if isempty(obj.sim_events)
                    obj.sim_events = sim_pars;
                else
                    obj.sim_events(end+1) = sim_pars;
                end
                
            end

        end % we recovered the injected event!
        
        function [flux, sim_pars] = addSimulatedEvent(obj, star_idx) % take a flux matrix and an index for a star and add a random occultation event on top of it
        
            if ~isempty(obj.cat) && obj.cat.success
                
                if ~ismember('FresnelSize', obj.cat.data.Properties.VariableNames)
                    obj.cat.addStellarSizes;
                end
                
                R = obj.cat.data.FresnelSize;
                if isempty(R)
                    obj.cat.addStellarSizes;
                end
                
                if isnan(R(star_idx)) % no stellar radius known from GAIA, just randomly pick one (from the distribution of R of other stars)
                    R_star = util.stat.inverseSampling(R, 'max', 3);  
                elseif R(star_idx)>3
                    R_star = 3; 
                else
                    R_star = R(star_idx);
                end
                
                [lc, sim_pars] = obj.bank.gen.randomLC('stellar_size', R_star); % let's set the stellar radius R=0.5 (Fresnel units)
                
            end
            
            flux_all = obj.store.extended_flux; 
            flux = flux_all(:,star_idx); % pick out the one star
            bg = obj.store.extended_aux(:,star_idx,obj.store.aux_indices.backgrounds).*...
                obj.store.extended_aux(:,star_idx,obj.store.aux_indices.areas);
            
            flux = flux - bg; % don't forget to add the background at the end
            
            margin = length(flux) - length(lc); % margins on both sides of the template, combined
            
            t = (1:length(flux))'; % fake timestamps for the fit
            
            fr = util.fit.polyfit(t, flux, 'order', 2); % fit the original lightcurve to a second order polynomial
            
            flux_mean = fr.ym; 
            flux_noise = flux - flux_mean; % separate the noise from the mean flux
            
            B = obj.store.checker.bad_times; % bad times matrix
            
            good_times = obj.store.search_start_idx + find(~B(obj.store.search_start_idx:obj.store.search_end_idx,star_idx)); % mark regions of the star's flux that are not disqualified by the checker
            
            if isempty(good_times) 
                error('this shouldn''t happen!'); % we chose a star that has good times in its lightcurve! 
            end
            
            if margin>0
                
                if length(good_times)>1 % choose a random point inside the good_times regions
                    peak_idx = randperm(length(good_times),1); 
                else
                    peak_idx = 1; % only one option, no need for random choice
                end
                
                peak = good_times(peak_idx); % position of the simulated event's center
                
                lc2 = vertcat(lc, ones(margin,1)); % make the template long enough to fit the real lightcurve
                lc2 = circshift(lc2, peak - floor(length(lc)/2) - 1); % shift it until the middle of the template is at "peak"
                
            elseif margin<0
                error('need to handle this case at some point...'); % this would be a problem if we inject events with W>=8 seconds
            else
                error('need to handle this case at some point...'); 
            end
            
            flux_mean = flux_mean.*lc2; % the raw flux is multiplied by the template (that should be 1 outside the occulation region)
            
            flux = flux_mean + flux_noise + bg; % now we can re-attach the noise and the background we extracted before
            
            % add some parameters known from the simulation
            sim_pars.time_index = peak; % what was the real peak time of this event
            sim_pars.star_index = star_idx; % what was the star index
            sim_pars.star_snr = obj.store.star_snr(obj.store.star_indices(star_idx)); % translate to the global star index to get the S/N from the store
            sim_pars.calc_snr = sqrt(nansum((lc-1).^2)).*sim_pars.star_snr; % estimate the event S/N with analytical formula (and assuming the star has white, gaussian noise)
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = trig.gui.EvFinderGUI(obj); 
            end
            
            obj.gui.make; 
            
        end
        
        function showCandidates(obj, varargin)
            
            if isempty(obj.cand)
                return;
            end
            
            if ~isempty(obj.gui) && obj.gui.check
                args = horzcat({'parent', obj.gui.panel_image}, varargin); 
            else
                args = horzcat({'parent', gcf}, varargin); 
            end
            
            obj.cand.show(args{:}); 
            
        end
        
        function showLatest(obj) % call this from e.g., the Analysis to update the plot if there is a new candidate
            
            if isempty(obj.cand)
                return;
            end
            
            if isempty(obj.gui.panel_image.Children) || isempty(obj.gui.panel_image.UserData) || ... 
                    ~isstruct(obj.gui.panel_image.UserData) || ~isfield(obj.gui.panel_image.UserData, 'index') % no proper plotting has been done so far
                
                obj.showCandidates('index', obj.num_candidates); 
                
            elseif obj.gui.panel_image.UserData.show_kept==0 && obj.gui.panel_image.UserData.index<obj.num_candidates % not showing the last candidate in the list
                obj.showCandidates('index', obj.num_candidates); 
            elseif obj.gui.panel_image.UserData.show_kept % if we are only showing the kept events, we must check if there is a more recent kept event to show:

                index = 0; 
                
                for ii = obj.gui.panel_image.UserData.index+1:obj.num_candidates

                    if obj.cand(ii).kept
                        index = ii; 
                    end

                end
            
                if index
                    obj.showCandidates('index', index); 
                end
                
            end
            
        end
        
        function popupStars(obj)
            
            f = util.plot.FigHandler('Star S/N'); 
            f.width = 26;
            f.height = 16;
            f.clear;
            
            ax = axes('Parent', f.fig); 
            
            obj.showStars('ax', ax, 'mag', 1); 
            
        end
        
        function showStars(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('threshold', true); % show the threshold used by the DataStore
            input.input_var('magnitudes', false); % show each star's magnitude
            input.input_var('axes', [], 'axis'); % which axes to plot to? default is gca()
            input.input_var('font_size', 20); % fonts on the axes
            input.scan_vars(varargin{:}); 
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            plot(input.axes, obj.store.star_snr, 'DisplayName', 'S/N per sample'); 
            
            if input.threshold
                
                N = length(obj.store.star_snr); 
                
                hold_state = input.axes.NextPlot; 
                
                hold(input.axes, 'on'); 
                
                plot(input.axes, 1:N, ones(N,1).*obj.store.pars.threshold, '--g', ...
                    'DisplayName', sprintf('threshold= %4.2f', obj.store.pars.threshold)); 
                
                input.axes.NextPlot = hold_state; 
                
            end
            
            if input.magnitudes && ~isempty(obj.cat) && obj.cat.success
                
                yyaxis(input.axes, 'right'); 
                
                plot(obj.cat.magnitudes, 'DisplayName', 'GAIA B_p mag'); 
                
                input.axes.YDir = 'reverse'; 
                
                ylabel(input.axes, 'magnitudes'); 
                
                legend(input.axes); 
                
                yyaxis(input.axes, 'left'); 
                
            end
            
            ylabel(input.axes, 'Star S/N per sample');
            xlabel(input.axes, 'Star index'); 
            
            input.axes.YLim = [0 max(obj.store.star_snr).*1.2]; 
            
            input.axes.FontSize = input.font_size; 
            
        end
        
        function popupSizes(obj)
            
            f = util.plot.FigHandler('Star sizes'); 
            f.width = 26;
            f.height = 16;
            f.clear;
            
            ax = axes('Parent', f.fig); 
            
            obj.showSizes('ax', ax, 'lines', [0.5 3]); 
            
        end
        
        function showSizes(obj, varargin)
            
            if isempty(obj.cat) || obj.cat.success==0
                error('No catalog was matched to this run!'); 
            end
            
            input = util.text.InputVars;
            input.input_var('distance_au', 40); % the distance at which to calculate the Fresnel scale... 
            input.input_var('log', true); % show the log10 of the Fresnel sizes. 
            input.input_var('lines', []); % add vertical lines showing some limits (give a scalar/vector of values)
            input.input_var('axes', [], 'axis'); % which axes to plot to? default is gca()
            input.input_var('font_size', 20); % fonts on the axes
            input.scan_vars(varargin{:}); 
            
            obj.cat.addStellarSizes;
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            S = obj.cat.data.FresnelSize(obj.store.star_indices);
            
            if input.log
                S = log10(S); 
            end
            
            histogram(input.axes, S, 'BinWidth', 0.05); 
            
%             input.axes.YScale = 'log'; 
            
            if input.log
                xlabel(input.axes, 'log10(Fresnel size)'); 
                input.axes.XLim = [-1 1.5]; 
            else
                xlabel(input.axes, 'Fresnel size'); 
                input.axes.XLim = [0 20]; 
            end
            
            [N,E] = histcounts(S, 'BinWidth', 0.05); 
                
            input.axes.YLim = [1, max(N)*1.1]; 
            
            if ~isempty(input.lines)
                
                if input.log
                    lines = log10(input.lines);
                else
                    lines = input.lines;
                end
                
                hold_state = input.axes.NextPlot;
                hold(input.axes, 'on'); 
                
                for ii = 1:length(input.lines)
                    
                    h = plot(input.axes, lines(ii).*[1 1], input.axes.YLim, '--r'); 
                    
%                     [~, idx] = min(abs(lines(ii)-E)); % find the bin index of this line
                    
                    text(input.axes, lines(ii)+0.05, input.axes.YLim(2), ...
                        sprintf('R= %4.2f ', input.lines(ii)), 'HorizontalAlignment', 'right', ...
                        'Color', h.Color, 'FontSize', input.font_size-2, 'Rotation', 90); 
                end
                
                hold(input.axes, 'off'); 
                
            end
            
            ylabel(input.axes, 'Number of stars in analysis'); 
            
            input.axes.NextPlot = hold_state; 
            input.axes.FontSize = input.font_size; 
            
        end
        
        function popupSNR(obj)
            
            f = util.plot.FigHandler('batch maximum S/N'); 
            f.width = 26;
            f.height = 16;
            f.clear;
            
            ax = axes('Parent', f.fig); 
            
            obj.showSNR('ax', ax); 
            
        end
        
        function showSNR(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('axes', [], 'axis'); % which axes to plot to? default is gca()
            input.input_var('font_size', 20); % fonts on the axes
            input.scan_vars(varargin{:}); 
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            h = histogram(input.axes, obj.snr_values, 'BinWidth', 0.5); 
            
            hold_state = input.axes.NextPlot; 
            
            input.axes.YScale = 'log'; 
            input.axes.YLim(1) = 0.5; 
            
            hold(input.axes, 'on'); 
            
            h2 = plot(input.axes, obj.pars.threshold.*[1 1], input.axes.YLim, '--r'); 
            text(input.axes, obj.pars.threshold+1.5, input.axes.YLim(2), ...
                sprintf(' threshold= %4.2f', obj.pars.threshold), ...
                'Color', h2.Color, 'FontSize', input.font_size, 'Rotation', -90); 
            
            N_total = numel(obj.snr_values);
            N_thresh = nnz(obj.snr_values>=obj.pars.threshold);
            
            text(input.axes, input.axes.XLim(2)-2, 10, sprintf('Number of batches= %d', N_total), ...
                'FontSize', input.font_size, 'Color', 'black', 'HorizontalAlignment', 'right'); 
            
            text(input.axes, input.axes.XLim(2)-2, 3, sprintf('batches above threshold= %d', N_thresh), ...
                'FontSize', input.font_size, 'Color', 'black', 'HorizontalAlignment', 'right'); 
            
            xlabel(input.axes, 'Best S/N for each batch'); 
            ylabel(input.axes, 'Number of batches'); 
            
            input.axes.NextPlot = hold_state; 
            input.axes.FontSize = input.font_size; 
            
        end
        
        function popupRunQuality(obj)
            
            f = util.plot.FigHandler('Run Quality'); 
            f.width = 26;
            f.height = 16;
            f.clear;
            
            ax = axes('Parent', f.fig); 
            
            obj.showRunQuality('ax', ax); 
            
        end
        
        function showRunQuality(obj, varargin)
            
            import util.series.binning; 
            
            input = util.text.InputVars;
            input.input_var('smooth', 100); % what factor to smooth the data?
            input.input_var('line', 3); % line width
            input.input_var('axes', [], 'axis'); % which axes to plot to? default is gca()
            input.input_var('font_size', 20); % fonts on the axes
            input.scan_vars(varargin{:}); 
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            j = obj.store.extended_juldates(end) - double(obj.store.timestamps - obj.store.extended_timestamps(end))/24/3600; 
            a = celestial.coo.airmass(j, obj.head.RA, obj.head.DEC, [obj.head.longitude, obj.head.latitude]./180.*pi);
            
            w = obj.store.checker.mean_width_values.*2.355.*obj.head.SCALE;
            b = obj.store.checker.mean_background_values;
            
            t = datetime(j, 'convertFrom', 'juliandate'); 
            t2 = t(1+obj.store.pars.length_burn_in:end); 
            
            plot(input.axes, t, a, 'DisplayName', 'airmass', 'LineWidth', input.line); 
            
            hold_state = input.axes.NextPlot; 
            
            hold(input.axes, 'on'); 
            
            if ~isempty(b)
                plot(input.axes, t2, b, '-', 'DisplayName', 'background [count/pix]', 'LineWidth', input.line); 
                plot(input.axes, binning(t2, input.smooth), binning(b, input.smooth), '-k',...
                    'DisplayName', 'background smoothed', 'LineWidth', input.line, 'HandleVisibility', 'off'); 
            end
            
            if ~isempty(w)
                
                yyaxis(input.axes, 'right');
                
                hold(input.axes, 'on'); 
            
                plot(input.axes, t2, w, '-', 'DisplayName', 'seeing ["]', 'LineWidth', input.line); 
                plot(input.axes, binning(t2, input.smooth), binning(w, input.smooth),'-k',...
                    'DisplayName', 'seeing smoothed', 'LineWidth', input.line, 'HandleVisibility', 'off'); 

                ylabel(input.axes, 'seeing FWHM ["]'); 
                
                input.axes.NextPlot = hold_state; 
                
                yyaxis(input.axes, 'left'); 
                
            end
            
            hl = legend(input.axes, 'Location', 'NorthWest'); 
            hl.FontSize = input.font_size;
            
            ylabel(input.axes, 'airmass / background'); 
            
            input.axes.NextPlot = hold_state; 
            input.axes.FontSize = input.font_size; 
            
        end
        
        function popupCuts(obj)
            
            f = util.plot.FigHandler('Quality cuts'); 
            f.width = 26;
            f.height = 16;
            f.clear;
            
            ax = axes('Parent', f.fig); 
            
            obj.showCuts('ax', ax); 
            
        end
        
        function showCuts(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('axes', [], 'axis'); % which axes to plot to? default is gca()
            input.input_var('font_size', 20); % fonts on the axes
            input.scan_vars(varargin{:}); 
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            input.axes.FontSize = input.font_size; 
            
        end
        
        function popupHours(obj)
            
            f = util.plot.FigHandler('Star hours'); 
            f.width = 26;
            f.height = 16;
            f.clear;
            
            ax = axes('Parent', f.fig); 
            
            obj.showHours('ax', ax); 
            
        end
        
        function showHours(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('axes', [], 'axis'); % which axes to plot to? default is gca()
            input.input_var('font_size', 20); % fonts on the axes
            input.scan_vars(varargin{:}); 
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            input.axes.FontSize = input.font_size; 
            
        end
        
        function popupPSD(obj)
            
            f = util.plot.FigHandler('Power Spectral Density'); 
            f.width = 26;
            f.height = 16;
            f.clear;
            
            ax = axes('Parent', f.fig); 
            
            obj.showPSD('ax', ax); 
            
            
        end
        
        function showPSD(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('stars', [1 10 30 100 300 1000]); % which stars are we going to use? 
            input.input_var('global_index', false); % use a global index (star indices are from the total stars found, not those that passed a threshold)
            input.input_var('magnitudes', true); % show each star's magnitude
            input.input_var('sqrt', false); % show the sqrt of the PSD
            input.input_var('half', true); % show only half the frequencies
            input.input_var('log_x', true); % show the x-axis on a log-scale (y-axis is always in log!)
            input.input_var('axes', [], 'axis'); % which axes to plot to? default is gca()
            input.input_var('font_size', 20); % fonts on the axes
            input.scan_vars(varargin{:}); 
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            %%%%%%%%%%%%%% parse star list %%%%%%%%%%%%%%%%%%%%%%%
            
            if input.global_index
                N_stars = length(obj.store.star_snr); % use all stars, including those that failed to pass
            else
                N_stars = length(obj.store.star_indices); % the number of stars are those that passed the threshold
            end
            
            if isempty(input.stars)
                error('Must input a non-empty "stars" input!'); 
            elseif ischar(input.stars)
                if util.text.cs(input.stars, 'first')
                    stars = 1;
                elseif util.text.cs(input.stars, 'all', ':')
                    stars = 1:N_stars; 
                else
                    error('Unknown option "%s" to input "stars". Use "first" or "all". ', input.stars); 
                end
            elseif isnumeric(input.stars) && isvector(input.stars)
                stars = input.stars;
            else
                error('Must input "stars" argument as numeric vector or string command ("first" or "all"). Instead got "%s" class input...', class(input.stars)); 
            end
            
            stars = stars(stars<N_stars); % remove stars out of range of the total number of stars we have
            
            if input.global_index % translate the star global index to the internal star index (those that passed the threshold)
                
                internal_stars = []; 
                
                for ii = 1:length(stars)
                    
                    idx = find(stars(ii)==obj.store.star_indices,1);
                    
                    if ~isempty(idx)
                        internal_stars(end+1) = idx; 
                    end
                    
                end
                
                stars = internal_stars; % these stars are the indices inside the PSD that match each one of the stars the user gave, but only if the passed the cut! 
            
            end

            %%%%%%%%%%%%%%%% prepare the data %%%%%%%%%%%%%%%%%%%%%%
            
            N = length(obj.psd.freq); 
            
            if input.half
                N = ceil(N/2); 
            end
            
            x = obj.psd.freq(1:N); 
            y = obj.psd.power_spectrum(1:N, stars); 
            
            if input.sqrt
                y = sqrt(y); 
            end
            
            %%%%%%%%%%%%% plot and adjust the settings %%%%%%%%%%%%%%%
            
            h = semilogy(input.axes, x, y);
            
            xlabel(input.axes, 'Frequency [Hz]'); 
            ylabel(input.axes, 'Power Spectral Density'); 
            
            if input.sqrt
                ylabel(input.axes, 'sqrt(PSD)'); 
            end
            
            if input.log_x
                input.axes.XScale = 'log'; 
            end
            
            input.axes.XLim = [x(1), x(end)]; 
            input.axes.YLim = [1 y(1).*10];
            
            input.axes.FontSize = input.font_size; 
            
            %%%%%%%%%%%%%% add magnitudes? %%%%%%%%%%%%%%%%%%%%%
            
            if input.magnitudes && ~isempty(obj.cat) && obj.cat.success
                
                if input.global_index
                    global_stars = input.stars; % the original input to this function IS the global star index! 
                else
                    global_stars = obj.store.star_indices(stars); % translate the internal star index to the global star index
                end
                
                for ii = 1:length(stars)
                    
                    y_factor = 2; 
                    if ii==length(stars)
                        y_factor = 1/2;
                    end
                    
                    text(input.axes, double(x(2)*2), double(y(1,ii))*y_factor, ...
                        sprintf('B_p= %4.1f', obj.cat.data.Mag_BP(global_stars(ii))), ...
                        'Color', h(ii).Color, 'FontSize', input.font_size-2); 
                end
                
            end
            
        end
        
        function popupSim(obj)
            
            if isempty(obj.summary)
                obj.produceSummary;
            end
            
            f = util.plot.FigHandler('Simulated events'); 
            f.width = 30;
            f.height = 20;
            f.clear;
            
            obj.summary.showSimulations('parent', f.fig); 
            
            f2 = util.plot.FigHandler('Simulation detection rate');
            f2.width = 30;
            f2.height = 20;
            f2.clear;
            
            obj.summary.showDetectionRate('parent', f2.fig); 
            
        end
        
    end 
    
end

