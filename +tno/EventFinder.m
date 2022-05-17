classdef EventFinder < handle
% This class contains the full pipeline for detecting occultation events. 
% 
% To use, make an object and make sure to give it the header and catalog:
% >> finder = tno.EventFinder;
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
% These parameters are defined in the "reset/clear" methods block, in the 
% resetPars() function. Refer to the inline docs for info on each parameter. 
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
%    This is applied to the fluxes before searching for occultations. 
%    To skip using the PSD calculation (which is slow) set 
%    pars.use_psd_correction=0. By default we use the PSD correction. 
% 3) Feed the corrected fluxes and auxiliary data into a search algorithm. 
%    This uses a matched-filter search using two filter banks saved as 
%    occult.ShuffleBank objects called "bank" and "bank_small". 
%    Additional filter banks called "bank_oort" and "bank_oort_small" are
%    used to store wider templates for finding Oort cloud objects 
%    (to use this, set use_oort=1). 
% 4) Inject simulated events into real data to check the recovery
%    efficiency of the pipeline and human vetters. 
%    Set pars.use_sim = 0 to skip using simulated events (default is 1). 
%    To determine how often simulated events are injected, set
%    "pars.num_sim_events_per_batch" which can be a small fraction (the 
%    default is 0.25, which means many simulated events are injected, but 
%    only a small fraction are saved). 
% 
% The output of the search is two-fold:
% a) The data quality cuts values and the amount of accumulated star hours 
%    are histogrammed and saved as a measure of the amount and quality of 
%    the data over the run. You can get a summary of this information, along
%    with header data and the parameters used in this run by each object by:
%    >> finder.produceSummary; 
%    which produces a Summary object. These objects can be collected for 
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

        gui@tno.gui.EvFinderGUI; % class added later
        
        bank@occult.ShuffleBank; % randomly picked filter kernels used for matched-filtering (loaded from file when needed)
        bank_small@occult.ShuffleBank; % a smaller filter bank used to weed out only the good stars for full-filtering
        bank_oort@occult.ShuffleBank; % a vector of filter banks for Oort cloud occultations at different distances
        bank_oort_small@occult.ShuffleBank; % a vector of filter banks for Oort cloud occultations with lower threshold for initial filtering
        
    end
    
    properties % objects
        
        head@head.Header; % link back to parent object, keep a lot of useful metadata
        cat@head.Catalog; % link back to catalog object, hopefully containing data on the star properties
        
        store@tno.DataStore; 
        
        psd@tno.CorrectPSD; % calculates the Power Spectra Density using Welch's method, then dereddens the flux
        
        var_buf@util.vec.CircularBuffer; % circular buffers for each star/kernel, holding a few variance measurements from past batches
        var_buf_oort@util.vec.CircularBuffer; % a vector of buffers for the Oort cloud templates at different distances
        
        cand@tno.Candidate;
        latest_candidates@tno.Candidate; 
        
        pars; % a struct with all the user-defined parameters
        
        summary@tno.Summary; 
        
        monitor@tno.StarMonitor; 
        
    end
    
    properties % inputs/outputs
        
        snr_values; % best S/N of the filtered fluxes, saved for each batch (use this to determine the ideal threshold)
        kernel_scores; % table of batch number, star S/N and kernel best score for each batch for stars passing the prefilter. 
        var_values; % keep track of the variance of each star/kernel over the length of the run (only when use_keep_variances=1)
        
        corrected_fluxes; % extended flux after correction (by PSD or by removing linear fit)
        corrected_stds; % standard deviation of the corrected fluxes
        
%         filtered_fluxes; % fluxes after matched-filtering with the bank of templates
%         bg_filtered_std; % measured background standard deviation after filtering
        
        black_list_stars; % a slowly growing list of stars that display too many events, so that all their events are marked as false
%         black_list_batches; % a slowly growing list of batches that display too many events, so that all their events are marked as false
        
        % simulated events are saved (a short struct describing the event). 
        % The passed events have been detected and the failed events have
        % been rejected by the system. Note that simulated events are placed
        % outside of the data cuts, so most of them are not rejected by them, 
        % but only by the S/N threshold. 
        sim_events; 
    
        flux_histograms = [];
        flux_histograms_log = [];
        flux_edges = [];
        flux_edges_log = [];
        
        total_batches = []; % optionally give the finder the number of files/batches in this run
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
        
%         background_ff; % filtered flux background data (either flux values or the variance from the var_buf)
        
        display_candidate_idx = []; % GUI parameter, which event is to be shown now
        use_display_kept_candidates = 0; % GUI parameter, to show only kept events or all events (default)
        use_show_truth = 0; % if true, will show the field name/coords and for each candidate will reveal if it is simulated
        
        version = 1.03;
        
    end
    
    methods % constructor
        
        function obj = EventFinder(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'tno.EventFinder')
                if obj.debug_bit>1, fprintf('EventFinder copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                
                if obj.debug_bit>1, fprintf('EventFinder constructor v%4.2f\n', obj.version); end
            
                obj.psd = tno.CorrectPSD;
            
                obj.store = tno.DataStore;
                
                obj.var_buf = util.vec.CircularBuffer; 
                
                obj.monitor = tno.StarMonitor; 
                
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
            obj.pars.max_time_spread = 40; % maximal filter cannot be larger than 30 frames in either direction (to prevent edge effects leaking into search region)
            
            obj.pars.max_events = 5; % we will search iteratively up to this many times for new candidates in each batch
            
            obj.pars.use_psd_correction = 1; % use welch on a flux buffer to correct red noise
            obj.pars.use_std_filtered = 1; % normalize variance of filtered flux of each kernel to the average of previous batches (averaging size is set by length_background in the store)
            obj.pars.use_detrend_after_psd = 1; % apply another linear detrend after PSD correction
            
            obj.pars.use_prefilter = 1; % filter all stars on a smaller bank, using a lower threshold, and then only vet the survivors with the full filter bank
            obj.pars.pre_threshold = 5.0; % threshold for the pre-filter
            
            obj.pars.use_oort = false; % use the Oort cloud template bank as well
            obj.pars.oort_distances = [3000]; % maybe add 1000 or 10,000 AU as well? 
            
            obj.pars.limit_events_per_batch = 5; % too many events in one batch will mark all events as black listed! 
            obj.pars.limit_events_per_star = 5; % too many events on the same star will mark all events on that star as black listed! 

            obj.pars.use_keep_variances = 1; % keep a copy of the variance of each star/kernel for the entire run (this may be a bit heavy on the memory)

            obj.pars.use_save_histograms = true;
            obj.pars.min_flux_histogram = 1e4;
            obj.pars.flux_binning_factors = [1 3 5];
        
            obj.pars.use_sim = 1; % add simulated events in a few batches randomly spread out in the whole run
            obj.pars.num_sim_events_per_batch = 0.25; % fractional probability to add a simulated event into each new batch. 
            obj.pars.use_keep_simulated = true; % if false, the simulated events would not be kept with the list of detected candidates (for running massive amount of simualtions)
            obj.pars.sim_max_R = 10; % maximum value of stellar size for simulated events
            obj.pars.sim_power_law = 3.5; % KBO radius power law (use positive values but produce negative slope)
            obj.pars.sim_rescale_noise = true;
            
            obj.reset;
            
        end
        
        function deleteBanks(obj)
            
            obj.bank = occult.ShuffleBank.empty;
            obj.bank_small = occult.ShuffleBank.empty;
            obj.bank_oort = occult.ShuffleBank.empty;
            obj.bank_oort_small = occult.ShuffleBank.empty;
            
        end
        
        function reset(obj)
            
            obj.cand = tno.Candidate.empty;
            obj.latest_candidates = tno.Candidate.empty;
            
            obj.black_list_stars = [];
            
            obj.display_candidate_idx = [];
            
            obj.store.reset;
            obj.var_buf.reset;
            obj.var_buf_oort = util.vec.CircularBuffer.empty; % this is initialized after loading the Oort template banks
            obj.psd.reset;
            
            obj.snr_values = [];
            obj.var_values = [];
            
            obj.sim_events = [];
            
            obj.total_batches = [];
            
            obj.flux_histograms = [];
            obj.flux_histograms_log = [];
            obj.flux_edges = [];
            obj.flux_edges_log = [];
            
            obj.kernel_scores = []; 
            
            obj.monitor.reset;
            
            obj.deleteBanks;
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            obj.psd.clear;
            obj.store.clear;
            
            obj.corrected_fluxes = [];
            obj.corrected_stds = [];
            
            obj.latest_candidates = tno.Candidate.empty;
            
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
                
                val = sprintf('aperture: %s pix | annulus: %s pix | iterations: %d | PSD: %d ', ...
                    util.text.print_vec(obj.head.PHOT_PARS.aperture_radius),...
                    util.text.print_vec(obj.head.PHOT_PARS.annulus_radii),...
                    obj.head.PHOT_PARS.iterations, obj.pars.use_psd_correction);
                
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
        
        function printNumDetections(obj)
            
            N = numel(obj.sim_events); 
            P = sum([obj.sim_events.passed]);
            [low, high] = util.stat.poisson_errors(P); % two sigma
            
            fprintf('Detected %d [%4.1f - %4.1f] out of %d events.\n', ...
                P, low, high, N); 
        end
        
    end
    
    methods % setters
        
        function set.head(obj, val) % setting the header of this object also cascades down to the store and its sub-objects
            
            obj.head = val;
            
            obj.store.head = val;
            
        end
        
    end
    
    methods % interface methods
        
        function input(obj, varargin) % provide the photometric data as varargin pairs, or as a util.text.InputVar object, or as an img.Photometry object
            
            obj.clear;
            
            if isempty(obj.store.star_sizes) && ~isempty(obj.cat) && ~isempty(obj.cat.success) && obj.cat.success
                obj.cat.addStellarSizes; 
                obj.store.star_sizes = obj.cat.data.FresnelSize;
            end
            
            obj.store.input(varargin{:}); % the store does the actual parsing and organizing of data into buffers
            
            obj.pars.analysis_time = util.text.time2str('now'); % keep a record of when the analysis was done
            
            if isempty(obj.bank) % lazy load the KBOs filter bank from file
                obj.loadFilterBankKBOs; 
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
            
            if obj.pars.use_save_histograms
                obj.calcHistograms;
            end
            
            if obj.pars.use_psd_correction
                    
                flux = obj.store.extended_detrend;
                aperture_index = obj.store.aperture_index;
                if isempty(aperture_index) || aperture_index > size(flux,3)
                    % choose last flux aperture by default, 
                    % or when flux is already only 2D, set aperture_index=1
                    aperture_index = size(flux,3); 
                end
                
                if size(flux,1)>=obj.store.pars.length_extended % only start adding data to the PSD after 2 batches are loaded

                    flux = flux(1:obj.store.pars.length_search,:,aperture_index); % grab the first frames from the 2-batch extended region

                    obj.psd.frame_rate = 1./nanmedian(diff(obj.store.extended_timestamps)); 
%                     obj.psd.window_size = ceil(obj.store.pars.length_extended); 
                    obj.psd.num_points = obj.store.pars.length_background*2;

                    % remove stars that didn't pass the threshold
                    if size(flux,2) < size(obj.psd.flux_buffer,2)
%                         obj.psd.flux_buffer = obj.psd.flux_buffer(:,obj.store.star_indices); 
                        start_idx = size(obj.psd.flux_buffer,1) - obj.store.pars.length_psd;
                        if start_idx < 1
                            start_idx = 1;
                        end
                        obj.psd.flux_buffer = obj.store.flux_buffer(start_idx:end,:); 
                        obj.psd.power_spectrum = obj.psd.power_spectrum(:,obj.store.star_indices); % only keep PSDs for current stars
                        obj.psd.flux_buffer = obj.psd.removeOutliers(obj.psd.flux_buffer); % inpaint outliers with current PSD
                        obj.psd.calcPSD; % update PSD using new flux aperture
                    else
                        obj.psd.addToBuffer(flux, obj.store.pars.length_psd); % also calculates the PSD
                    end

                end
                
            end
            
            if obj.store.is_done_burn % do not do any more calculations until done with "burn-in" period
                
                %%%%%%%%%%%%%% PSD CORRECTION %%%%%%%%%%%%%%%
                
%                 if obj.pars.use_psd_correction
%                     
%                     try
%                         obj.psd.calcPSD; % first off, make the PSD correction for the current batch
%                     catch ME
%                         if strcmp(ME.identifier, 'MATLAB:nomem')
%                             util.text.date_printf('Out of memory while caclulating PSD! Trying again...');
%                             pause(30);
%                             obj.psd.calcPSD; % try again to make the PSD correction for the current batch
%                         else
%                             rethrow(ME);
%                         end
%                     end
% 
%                 end
                
                [obj.corrected_fluxes, obj.corrected_stds] = obj.correctFluxes(obj.store.extended_detrend); % runs either PSD correction or simple removal of linear fit to each lightcurve

                %%%%%%%%%%% FILTER FOR KBOS %%%%%%%%%%%%%%
                
                [obj.latest_candidates, star_indices, best_snr, new_scores] = obj.applyTemplateBank('kbos'); % by default use corrected_fluxes and correctd_stds given above
                
                obj.kernel_scores = vertcat(obj.kernel_scores, new_scores); 
                
                obj.cand = vertcat(obj.cand, obj.latest_candidates); % store all the candidates that were found in the last batch
                
                %%%%%%%%%%%%%% STAR HOURS %%%%%%%%%%%%%%
                
                % save the star hours, but also check if this batch is black-listed
                if length(obj.latest_candidates(~[obj.latest_candidates.is_simulated]))>=obj.pars.limit_events_per_batch % this batch has too many events, need to mark them as black-listed! 

                    for ii = 1:length(obj.latest_candidates)
                        obj.latest_candidates(ii).notes{end+1} = sprintf('Batch is black listed with %d events', length(obj.latest_candidates)); 
                        obj.latest_candidates(ii).kept = 0;
                    end

                    obj.store.saveHours(1); % the parameter 1 is used to mark this as a bad batch...

                else
                    obj.store.saveHours; % with no arguments it just counts the hours in this batch as good times
                end

                %%%%%%%% MONITOR SPECIFIC STARS %%%%%%%%%%%%
                
                obj.monitor.input(obj.store.extended_timestamps, obj.store.extended_juldates, ... % track additional data on specific star indices
                    obj.store.extended_flux, obj.store.extended_aux, obj.store.cutouts, obj.latest_candidates);
                
                %%%%%%%%% KEEP TRACK OF VARIANCE VALUES %%%%%%%%
                
                if obj.pars.use_keep_variances && ~obj.var_buf.is_empty % just like the var buffer, only for entire run (this accumulates to a big chunk of memory)
                    obj.var_values = vertcat(obj.var_values, obj.var_buf.data(obj.var_buf.idx,:,:)); 
                end
                
                %%%%%%%% OORT CLOUD %%%%%%%%%%%%%
                
                if isempty(obj.bank_oort) || isempty(obj.bank_oort_small) || isempty(obj.var_buf_oort) 
                    obj.loadFilterBankOort; % lazy load the Oort filter bank from file
                end

                if obj.pars.use_oort && ~isempty(obj.bank_oort)
                
                    [oort_candidates, oort_star_indices] = obj.applyTemplateBank('oort'); % by default use corrected_fluxes and correctd_stds given above
                                        
                    star_indices = unique([star_indices; oort_star_indices]); % add the stars selected by the Oort prefilter / filter

                    obj.cand = vertcat(obj.cand, oort_candidates);
                    
                end
                
                %%%%%%%%%%%% SIMULATIONS %%%%%%%%%%%%%%%
                
                if obj.pars.use_sim % simulated events are injected into the data and treated like real events

                    star_indices_sim = obj.findStarsForSim(star_indices); % get stars that do not include any possible real candidates

                    for ii = 1:obj.getNumSimulations % this number can be more than one, or less (then we randomly decide if to include a simulated event in this batch)
                        
                        [sim_cand, sim_ev] = obj.simulateSingleEvent(star_indices_sim); % each call to this function tries to add a single simulated event to a random star from the list                            
                        
                        if obj.pars.use_keep_simulated
                            obj.cand = vertcat(obj.cand, sim_cand); % add the simulated events to the list of regular events
                        end
                        
                        if ~isempty(sim_ev)
                            
                            % add this sim event struct to the list of events
                            if isempty(obj.sim_events)
                                obj.sim_events = sim_ev;
                            else
                                obj.sim_events(end+1) = sim_ev;
                            end

                        end
                        
                    end

                end % use sim
                
                obj.snr_values(end+1) = best_snr; % this used to be before the Oort cloud and sim, but then the batch counter was wrong! 
                
            end % do not do any more calculations until done with "burn-in" period
            
        end
        
        function finishup(obj) % call this at the end of the run 
            
            if obj.pars.limit_events_per_star % check if any stars need to be black listed
                
                real_events = ~[obj.cand.is_simulated] & ~[obj.cand.oort_template]; % don't include simulated events or Oort cloud triggers in the black list! 

                N = histcounts([obj.cand(real_events).star_index], 'BinEdges', 1:size(obj.store.flux_buffer,2)+1); % star index among the subset of stars in the analysis

                idx = N>=obj.pars.limit_events_per_star;

                obj.black_list_stars = find(idx); 

                obj.store.checker.hours.removeStars(idx); 

                for ii = 1:length(obj.cand)

                    star = obj.cand(ii).star_index; 

                    if ismember(star, obj.black_list_stars)
                        str = sprintf('Star %d blacklisted with %d events', star, N(star)); 
                        str_idx = strcmp(str, obj.cand(ii).notes);
                        if nnz(str_idx)==0 % none of the notes includes this black list note
                            obj.cand(ii).notes{end+1} = str; 
                        end
                        obj.cand(ii).kept = 0; 
                    end

                    obj.cand(ii).best_snrs = obj.snr_values; % keep track of the best S/N during the whole run
                    
                end

            end % check if any stars need to be black listed
            
            obj.summary = obj.produceSummary; 
            
            if obj.store.is_done_burn
                obj.summary.num_events_expected = obj.summary.getNumDetections; % use simulations to estimate the detection rate
            end
            
        end
        
        function s = produceSummary(obj) % makes a RunSummary with all relevant data on the quality and quantity of star hours, etc
            
            % add varargin? 
            
            if ~obj.store.is_done_burn
                s = tno.Summary.empty;
                return;
            end
            
            s = tno.Summary; % generate a new object
            s.head = util.oop.full_copy(obj.head); % I want to make sure this summary is distinct from the original finder
            
            % load the content of the finder
            s.finder_pars = obj.pars;
            s.snr_values = obj.snr_values;
            s.batch_counter = obj.batch_counter; 
            s.total_batches = obj.total_batches; 
            s.sim_events = obj.sim_events;
            s.black_list_stars = obj.black_list_stars;
            
            % load the content of the store
            s.store_pars = obj.store.pars;
            s.good_stars = obj.store.star_indices;
            s.star_snr = obj.store.star_snr;
            s.star_sizes = obj.store.star_sizes;
            s.size_snr_coeffs = obj.store.size_snr_coeffs;
            
            
%             FWHM = obj.store.checker.defocus_log*2.355*obj.head.SCALE;
% 
%             if ~isempty(FWHM)
%                 s.fwhm_edges = 0:0.1:round(nanmax(FWHM)*10)/10;
%                 s.fwhm_hist = histcounts(FWHM, 'BinEdges', s.fwhm_edges);             
%                 s.fwhm_hist = s.fwhm_hist.*nanmedian(diff(obj.store.checker.juldate_log))*24*3600;
%             end
            
            s.seeing_log = obj.store.fwhm_log; 
            s.juldates_log = obj.store.juldates_log; 
            s.airmass_log = celestial.coo.airmass(s.juldates_log, obj.head.RA, obj.head.DEC, [obj.head.longitude, obj.head.latitude]./180.*pi);
            b = obj.store.checker.mean_background_values;
            b = [NaN(obj.store.pars.length_burn_in,1); b]; 
            s.background_log = util.series.binning(b, obj.store.pars.length_search); % should be the same size as the other logs, with some NaNs in the beginning for the burn in
            
            s.flux_histograms = obj.flux_histograms;
            s.flux_edges = obj.flux_edges;
            s.flux_histograms_log = obj.flux_histograms_log;
            s.flux_edges_log = obj.flux_edges_log;
            s.flux_binning_factors = obj.pars.flux_binning_factors;
            
            % load the content of the checker
            s.checker_pars = obj.store.checker.pars;
            s.cut_names = obj.store.checker.cut_names;
            s.cut_indices = obj.store.checker.cut_indices; 
            s.cut_thresholds = obj.store.checker.cut_thresholds;
            s.cut_two_sided = obj.store.checker.cut_two_sided;
            s.cut_histograms = permute(nansum(obj.store.checker.histograms,2), [1,3,2]);
            s.cut_bin_edges = obj.store.checker.hist_edges;
            
            % make sure the catalog has Fresnel sizes
            if ~isempty(obj.cat) && obj.cat.success
                if ~ismember('FresnelSize', obj.cat.data.Properties.VariableNames)
                    obj.cat.addStellarSizes;
                end
            end
            
            % load the content of the star hours
            s.inputHours(obj.store.checker.hours, obj.cat.data.FresnelSize(obj.store.star_indices));
            
            obj.summary = s; 
            
        end
        
        function val = checkBatchGood(obj)
            
            val = obj.store.checker.checkBatchGood; 
            
        end
        
    end
    
    methods(Hidden=true) % internal calculation methods
        
        function calcHistograms(obj)
            
            f = obj.store.this_input.fluxes(:,:,end);
            F = nanmean(f); 

            idx = F > obj.pars.min_flux_histogram;

            if isempty(obj.flux_edges)
                obj.flux_edges = -5:0.01:5; 
            end
            
            if isempty(obj.flux_edges_log)
                obj.flux_edges_log = -3:0.01:3; 
            end

            f_norm = (f(:,idx) - F(idx)) / F(idx); 
            f_ratio = f(:,idx) / F(idx); 
            
            b = obj.pars.flux_binning_factors;

            if isempty(b)
                b = 1;
            end

            if isempty(obj.flux_histograms)
                obj.flux_histograms = zeros(length(obj.flux_edges)-1, length(b), 'uint64'); 
            end
            
            if isempty(obj.flux_histograms_log)
                obj.flux_histograms_log = zeros(length(obj.flux_edges_log)-1, length(b), 'uint64'); 
            end

            for ii = 1:length(b)
                
                if b(ii) > 1
                    ff = filter2(ones(b(ii),1)/b(ii), f_norm);
                else
                    ff = f_norm;
                end
                
                N = histcounts(ff, obj.flux_edges); 
                
                obj.flux_histograms(:,ii) = obj.flux_histograms(:,ii) + uint64(N'); 
                
                if b(ii) > 1
                    ff = filter2(ones(b(ii),1)/b(ii), f_ratio);
                else
                    ff = f_ratio;
                end
                
                N = histcounts(log2(abs(ff)+1e-10), obj.flux_edges_log); 
                
                obj.flux_histograms_log(:,ii) = obj.flux_histograms_log(:,ii) + uint64(N'); 
                
            end
            
        end
        
        function [flux_out, std_out] = correctFluxes(obj, fluxes, star_indices) % star indices tells the PSD which stars were given in fluxes (default is all stars)
            
            if nargin<3 || isempty(star_indices) % by default run the correction for all stars
                star_indices = 1:size(fluxes,2);
            end
            
            if obj.pars.use_psd_correction
                
                flux_out = obj.psd.input(fluxes, 1, star_indices); 
                
                if obj.pars.use_detrend_after_psd
                    flux_out = obj.removeTrendByBatches(flux_out); 
                end
                
            else
                flux_out = fluxes; 
%                 flux_out = obj.removeLinearFit(fluxes, obj.store.extended_timestamps); 
%                 flux_out = fillmissing(flux_out, 'linear'); 
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
        
        function flux_out = removeTrendByBatches(obj, flux)
            
            flux_out = flux; 
            
            idx1 = find(obj.store.extended_frame_num==1); % indices where a new batch begins
            idx2 = find(obj.store.extended_frame_num==obj.head.NAXIS3); % indices where the batch ends

            for ii = 1:length(idx1)

                if length(idx2)<ii
                    f = flux(idx1(ii):end,:); 
                else
                    f = flux(idx1(ii):idx2(ii),:); 
                end

                f = util.series.detrend(f, 'iterations', 2); % remove linear fit to the data

                if length(idx2)<ii
                    flux_out(idx1(ii):end,:) = f; 
                else
                    flux_out(idx1(ii):idx2(ii),:) = f; 
                end

            end
            
        end
        
        function [candidates_all, star_indices_all, best_snr_all, scores_all] = applyTemplateBank(obj, bank_name, fluxes_det, fluxes_corr, stds)
        % do a prefilter (optional) and then a filter and then find candidates. 
        % INPUTS: 
        %   -bank_name: a string: "kbos" or "oort". 
        %   -fluxes: the corrected flux (after removing linear fit or PSD
        %            correction. Default is obj.corrected_fluxes. 
        %   -std: the standard deviation of the fluxes. Default is 
        %         obj.corrected_stds. 
        %
        % OUTPUTS: 
        %   -candidates_all: a vector of tno.Candidate objects, containing
        %                all successfull triggers. If none are found, returns
        %                an empty vector. 
        % 
        %   -star_indices_all: a vector containing all the stars that passed
        %                  either the pre-filter or (if all stars passed)
        %                  only the stars that triggered the full filter. 
        %   -best_snr_all: maximal S/N value measured from the filtered fluxes.
        %   -scores_all: a table containing each kernel's best response for
        %                each star that passed the prefilter.  
        %   
        % NOTE: if there are multiple template banks, it will run them in
        %       sequence and add up the results. The best S/N is max of all
        %       results. 
        
            candidates_all = tno.Candidate.empty;
            star_indices_all = [];
            best_snr_all = 0;
            scores_all = []; 
            
            if nargin<2 || isempty(bank_name)
                bank_name = 'kbos';
            end
            
            if nargin<3 || isempty(fluxes_det)
                fluxes_det = obj.store.extended_detrend;
            end
            
            if nargin<4 || isempty(fluxes_corr)
                fluxes_corr = obj.corrected_fluxes;
            end
            
            if nargin<5 || isempty(stds)
                stds = obj.corrected_stds;
            end
            
            if util.text.cs(bank_name, 'kbos')
                bank = obj.bank;
                bank_small = obj.bank_small;
                var_buf = obj.var_buf; 
            elseif util.text.cs(bank_name, 'oort')
                bank = obj.bank_oort;
                bank_small = obj.bank_oort_small;
                var_buf = obj.var_buf_oort; 
            else
                error('Unknown "bank_name" option "%s". Use "KBOs" or "Oort". ', bank_name); 
            end
            
            for ii = 1:length(bank)
            
                candidates = tno.Candidate.empty;
                scores = [];
                
                if obj.pars.use_prefilter % use a small filter bank and only send the best stars to the big bank

                    fluxes_filtered_small = obj.filterFluxes(fluxes_corr, stds, bank_small(ii), [], var_buf(ii)); 

                    pre_snr = squeeze(util.stat.max2(fluxes_filtered_small)); % signal to noise peak for each star at any time and any filter kernel

                    star_indices = find(pre_snr>=obj.pars.pre_threshold); % only stars that have shown promising response to the small filter bank are on this list 

                else % do not prefilter, just put all stars into the big filter bank

                    star_indices = util.vec.tocolumn(1:size(fluxes_corr,2)); % just list all the stars
                    pre_snr = 0; 
                    
                end

                best_snr = nanmax(pre_snr); % need to keep track what is the best S/N for this batch, regardless of which filter and regardless of if there were any candidates detected. 

                if ~isempty(star_indices) % if no stars passed the pre-filter, we will just skip to the next batch (if no pre-filter is used, all stars will pass)

                    if obj.pars.use_prefilter
                        var_buf_large_filter = []; % each batch we will have different stars at this stage, cannot use var_buf
                    else
                        var_buf_large_filter = var_buf(ii); % always use all stars, no problem applying the var buffer
                    end

                    filtered_fluxes = obj.filterFluxes(fluxes_corr, stds, bank(ii), star_indices, var_buf_large_filter); 

                    [candidates, best_snr] = obj.searchForCandidates(obj.store.extended_flux, fluxes_det, fluxes_corr, filtered_fluxes, star_indices, bank(ii)); % loop over the normalized filtered flux and find multiple events

                    if length(star_indices)==size(fluxes_corr,2) % if all stars passed prefilter (or when skipping pre-filter)
                        star_indices = [candidates.star_index]; % only keep stars that triggered
                    end
                    
                    if nargout >= 4
                        % the kernel response
                        time_indices = obj.store.search_start_idx:obj.store.search_end_idx;
                        kernel_best = nanmax(filtered_fluxes(time_indices, :, :), [], 1);
                        kernel_best = permute(kernel_best, [3,2,1]); % get rid of reduced first index
                        
                        % the star S/N
                        A = obj.store.background_aux(:,star_indices, obj.store.aux_indices.areas); 
                        B = obj.store.background_aux(:,star_indices, obj.store.aux_indices.backgrounds); 
                        flux = obj.store.background_flux(:,star_indices)-A.*B; % background corrected flux
                        flux_mean = nanmedian(flux, 1); 
                        flux_std = mad(flux, 1)./0.6745; % match the MAD to STD using sqrt(2)*erfinv(0.5)
                        star_snr = (flux_mean./flux_std)';
                        
                        % the batch number
                        batch_number = ones(length(star_snr), 1) * obj.batch_counter; 
                        
                        scores = table(batch_number, star_snr, kernel_best, ...
                            'VariableNames', {'batch_number', 'star_snr', 'kernel_scores'});
                        
                    end

                end % if no stars passed the pre-filter, we will just skip to the next batch
            
                candidates_all = vertcat(candidates_all, candidates); 
                star_indices_all = unique([star_indices_all; star_indices]); 
                best_snr_all = nanmax(best_snr_all, best_snr); 
                scores_all = vertcat(scores_all, scores);
                
            end
            
        end
        
        function ff = filterFluxes(obj, fluxes, stds, bank, star_indices, var_buf, bg_fluxes) % filter the fluxes with the template bank
        % Run the input fluxes through the template bank. 
        % If given a var buffer, it will also update it with the latest
        % variance values measured from the filtered fluxes. 
        %
        % INPUTS:
        %   -fluxes: should be corrected (for PSD etc). Can include all
        %            stars or a subset given by star_indices. 
        %            2D matrix, dim 1 is time, dim 2 is star index. 
        %   -stds: standard deviations for the fluxes as given. 
        %          Row vector, one value per star. 
        %   -bank: an occult.ShuffleBank object with templates. 
        %   -star_indices: optional argument, specifying which stars were
        %                  chosen (e.g., by prefilter). Default is all stars. 
        %   -var_buf: a util.vec.CircularBuffer object with variance per
        %             star from the previous batches. If not given, or if
        %             it is empty, the STD correction will be calculated by
        %             filtering the background flux. 
        %   -bg_fluxes: optionally give the flux of the background region.
        %               If not given, will just take the background_flux
        %               from the data store. 
        % 
            
            N = ceil(obj.store.pars.length_background./obj.store.pars.length_search); % how many batches we want to keep as "background" sample in the var_buf

            if nargin<4 || isempty(bank) || ~isa(bank, 'occult.ShuffleBank')
                error('Must supply a valid ShuffleBank object!');  
            end
            
            if nargin<5 || isempty(star_indices)
                star_indices = 1:size(fluxes,2); % just use all stars
            end
            
            if nargin<6 || isempty(var_buf)
                var_buf = []; % no buffer is given, must calculate variance on background fluxes
            else
                if var_buf.is_empty % we haven't used this since the last reset() call
                    var_buf.reset(N); 
                end
            end
            
            % first calculate the STD correction for the filtered flux
            if obj.pars.use_std_filtered
                
                if isempty(var_buf) || var_buf.counter<N  % no buffer (or the buffer is not yet full), must calculate its own variance on the background
                    
                    if nargin<7 || isempty(bg_fluxes) % if we were not given a background flux, use the store
                        [bg_fluxes, bg_stds] = obj.correctFluxes(obj.removeOutliers(obj.store.background_detrend(:,star_indices)), star_indices); 
                    end
                    
                    try
                        bg_filtered_fluxes = bank.input(bg_fluxes, bg_stds); 
                    catch ME
                        if strcmp(ME.identifier, 'MATLAB:nomem')
                            util.text.date_printf('Out of memory while filtering fluxes! Trying again...'); 
                            pause(30); 
                            bg_filtered_fluxes = bank.input(bg_fluxes, bg_stds); 
                        else
                            rethrow(ME); 
                        end

                    end
                    
                    ff_std = nanstd(obj.removeOutliers(bg_filtered_fluxes)); % STD correction from filtered background fluxes
                    
                else
                    ff_std = sqrt(var_buf.mean); % STD correction from var buffer
                end
                
            else
                ff_std = 1; % no STD correction
            end
            
            try
                ff = bank.input(fluxes(:,star_indices), stds(1,star_indices));
            catch ME
                if strcmp(ME.identifier, 'MATLAB:nomem')
                    util.text.date_printf('Out of memory while filtering fluxes! Trying again...'); 
                    pause(30); 
                    ff = bank.input(fluxes(:,star_indices), stds(1,star_indices));
                else
                    rethrow(ME); 
                end

            end
            
            if ~isempty(var_buf) % if we are using a var buffer, load the new results into it
                
                var_buf.input(nanvar(obj.removeOutliers(ff))); 
                
            end
            
            ff = ff./ff_std; % apply the STD correction
            
        end
        
        function flux_out = removeOutliers(obj, flux, sigma)
            
            if nargin<3 || isempty(sigma)
                sigma = 5;
            end
            
            idx = abs(flux - nanmean(flux)) > sigma*nanstd(flux);
            flux_out = flux;
            flux_out(idx) = NaN; 
            
        end
            
        function [new_candidates, best_snr] = searchForCandidates(obj, fluxes_raw, fluxes_detrended, fluxes_corrected, fluxes_filtered, star_indices, bank, max_num_events) % loop over filtered fluxes, find above-threshold peaks, remove that area, then repeat

            if nargin<8 || isempty(max_num_events)
                max_num_events = obj.pars.max_events;
            end
            
            new_candidates = tno.Candidate.empty;
            
            time_indices = obj.store.search_start_idx:obj.store.search_end_idx;

            for ii = 1:length(max_num_events)

                [mx, idx] = util.stat.maxnd(abs(fluxes_filtered(time_indices,:,:))); % find the maximum point in the filtered fluxes (after correcting for background variance)

                idx(1) = time_indices(idx(1)); % adjust the idx(1) time index to be inside the extended region, instead of in the search region
                
                if length(idx)<3
                    idx(3) = 1;
                end
                
                if ii==1, best_snr = mx; end % return the best S/N from this batch
                
                if mx>=obj.pars.threshold % at least one part of the filtered lightcurves passed the threshold

                    c = obj.makeNewCandidate(fluxes_raw, fluxes_detrended, fluxes_corrected, fluxes_filtered, idx, star_indices, bank); 

                    c.serial = length(obj.cand) + length(obj.latest_candidates) + 1; % new event gets a serial number
                    new_candidates = vertcat(new_candidates, c); % add this event to the list of new candidates
                    
                    % remove the new event's footprint from the fluxes to search
                    fluxes_filtered(c.time_range, :, :) = NaN; % don't look at the same region twice

                    % do we have to ALSO reduce this region from the star hours?
                    
                end

            end
            
        end
        
        function c = makeNewCandidate(obj, fluxes_raw, fluxes_detrended, fluxes_corrected, fluxes_filtered, idx, star_indices, bank) % generate a Candidate object and fill it with data

            c = tno.Candidate;
            
            c.snr = abs(fluxes_filtered(idx(1),idx(2),idx(3))); % note this is positive even for negative filter responses! 
            c.is_positive = fluxes_filtered(idx(1),idx(2),idx(3))>0; % negative events are this where the filter response is large negative (e.g., flares)
                        
            c.time_index = idx(1); % time index inside of the extended region
            c.kern_index = idx(2); % which kernel triggered
            c.star_index = star_indices(idx(3)); % which star (out of all stars that passed the burn-in)
            c.star_index_global = obj.store.star_indices(star_indices(idx(3))); % get the star index in the original list of stars before pre-filter and before burn-in
            c.star_snr = obj.store.star_snr(c.star_index_global); % get the S/N that was recorded at end of burn in
            
            c.template_bank = bank.getBankName; % keep a record of which template bank was used
            c.oort_template = bank.D_au>=300; % arbitrary cutoff to include template banks for the oort cloud
            c.kernel = bank.kernels(:,c.kern_index); % the lightcurve of the kernel used (zero centered and normalized)
            c.kern_props = bank.pars(c.kern_index); % properties of the occultation that would generate this kernel
            c.kern_props.inverted = ~c.is_positive; % if the filtered flux was negative, it means we had to invert the kernel
            
            % find continuous area (in the extended region) that is above time_range_thresh 
            c.time_range = obj.findTimeRange(fluxes_filtered, idx(1), idx(2), idx(3), bank); 
            
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
            c.flux_detrended = fluxes_detrended(:,c.star_index);
            c.flux_corrected = fluxes_corrected(:,c.star_index); 
            c.flux_filtered = fluxes_filtered(:,idx(2),idx(3));
            c.flux_extra = obj.store.extended_fluxes_extra(:,c.star_index,:); % keep track of alternative apertures (e.g., unforced or gaussian)
            c.extra_aperture_indices = obj.store.extra_fluxes_indices; % which photometry indices were saved 
            
            % store data for all stars as reference
            c.flux_raw_all = fluxes_raw; 
            c.flux_detrended_all = fluxes_detrended; 
            c.flux_corrected_all = fluxes_corrected; 
            c.auxiliary_all = obj.store.extended_aux; 
            
            % keep track of each star's average position
            x = nanmean(obj.store.extended_aux(:,:,obj.store.aux_indices.centroids_x),1)'; 
            y = nanmean(obj.store.extended_aux(:,:,obj.store.aux_indices.centroids_y),1)'; 
            
            c.positions = [x,y];
            
            idx = find(star_indices==c.star_index); % the index of the star inside the list of stars that passed the pre-filter
%             c.filtered_flux_past_values = obj.background_ff(:,idx); % because we only have this background_ff for the subset of stars
            c.flux_buffer = obj.store.flux_buffer(:,c.star_index); % flux history for this star
            c.detrend_buffer = obj.store.detrend_buffer(:,c.star_index); % detrended flux history for this star
            c.timestamps_buffer = obj.store.timestamps_buffer; % timestamps for that duration
            
            if ~isempty(obj.psd.power_spectrum)
                c.psd = obj.psd.power_spectrum(:,c.star_index); % the Power Spectral Density (PSD) for this star
                c.freq_psd = obj.psd.freq; % frequency axis for the PSD
            end
            
            % store some statistics on this star's raw flux
            A = obj.store.background_aux(:,c.star_index, obj.store.aux_indices.areas); 
            B = obj.store.background_aux(:,c.star_index, obj.store.aux_indices.backgrounds); 
            c.flux_mean = nanmedian(obj.store.background_flux(:,c.star_index)-A.*B); 
            c.flux_std = mad(obj.store.background_flux(:,c.star_index), 1)./0.6745; % match the MAD to STD using sqrt(2)*erfinv(0.5)
            
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
            c.average_offsets = obj.store.extended_average_offsets; % keep track of forced photometry offsets
%             c.fwhm = obj.store.checker.fwhm;
            
            % get the observational parameters and the star parameters from astrometry/GAIA
            c.head = util.oop.full_copy(obj.head);
            if ~isempty(obj.cat) && obj.cat.success
                c.star_props = obj.cat.data(c.star_index_global,:); % copy a table row from the catalog
            end

            % save the parameters used by the finder, store and quality-checker
            c.finder_pars = obj.pars; 
            c.store_pars = obj.store.pars;
            c.store_pars.aperture_index = obj.store.aperture_index; 
            c.store_pars.star_indices = obj.store.star_indices; 
            c.store_pars.star_sizes = obj.store.star_sizes;
            c.store_pars.star_snrs = obj.store.star_snr;
            c.global_star_indices = obj.store.star_indices; 
            c.aperture_index = obj.store.aperture_index; 
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
            
            for ii = 1:length(c.cut_names) % go over each cut and check if it applies to this event
            
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
                    
                    if isfield(obj.store.checker.pars, ['thresh_' c.cut_names{ii}]) ...
                            || strcmp(c.cut_names{ii}(1:4), 'corr') || strcmp(c.cut_names{ii}(1:9), 'flux_corr') % this is a "threshold" type test
                        c.cut_string{end+1,1} = sprintf('Cut "%s" failed with value %4.2f', c.cut_names{ii}, mx);
                    
                    else % this is a pass/fail type cut
                        c.cut_string{end+1,1} = sprintf('Cut "%s" failed', c.cut_names{ii});
                    end
                    
                    c.cut_value(end+1) = mx; % store the value of the cut that this event overlapped
                    c.cut_hits(end+1) = ii; % mark which cut number has triggered (could be more than one!)
                    c.kept = 0; % mark the event as not real
                    
                end % the event occurs on one of the flagged frames/stars
                
            end % go oveer each cut and check if it applies to this event
            
            c.calcTrackingErrorValue; % this gives the sum of the first three, I will instead focus on the 3rd star
            
%             if c.corr_flux(3)>obj.store.checker.pars.thresh_tracking_error % too much correlations with other stars
                str = sprintf('Correlation 3rd highest star: %4.2f', c.corr_flux(3)); 
                str_idx = strcmp(str, c.notes);
                if nnz(str_idx)==0
                    c.notes{end+1} = str; 
                end
%                 c.kept = 0; % I don't want to disqualify on this, but I
%                 do want the scanner to see a warning, and judge for
%                 themselves. 
%             end
            
            c.analysis_time = obj.pars.analysis_time;
            
        end
        
        function time_range = findTimeRange(obj, ff, time_index, kern_index, star_index, bank) % find all continuous time frames around peak that have surpassed a (lower) threshold, or just mark a minimal event region
            
            N = size(ff,1); % time length
            
            thresh = obj.getTimeThresh;
            time_range = [];
            
            kernel_min_spread = ceil(nnz(abs(bank.kernels(:,kern_index))>=0.02)/2); % minimal time range cannot be smaller than "active" part of kernel (i.e., above 2% deviation)
            
            for jj = 0:obj.pars.max_time_spread % go backward in time

                idx = time_index - jj;

                if idx<1, break; end

                if abs(ff(idx, kern_index, star_index))>=thresh || jj<=obj.pars.min_time_spread || jj<=kernel_min_spread
                    time_range = [time_range, idx];
                else
                    break;
                end

            end

            time_range = flip(time_range);

            for jj = 1:obj.pars.max_time_spread % go forward in time

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
        
        function loadFilterBankKBOs(obj, throw_error)
            
            if nargin<2 || isempty(throw_error)
                throw_error = 1;
            end
            
            frame_rate = floor(obj.head.FRAME_RATE); 
            
            obj.bank = occult.ShuffleBank.empty; 
            obj.bank_small = occult.ShuffleBank.empty; 
            
            f = fullfile(getenv('DATA'), sprintf('WFAST/occultations/templates_40AU_%dHz.mat', frame_rate));
            if exist(f, 'file')
                load(f, 'bank');
                obj.bank = bank;
            elseif throw_error
                error('Cannot load template bank from file "%s"', f); 
            end

            f = fullfile(getenv('DATA'), sprintf('WFAST/occultations/templates_40AU_%dHz_small.mat', frame_rate));
            if exist(f, 'file')
                load(f, 'bank');
                obj.bank_small = bank;
            elseif throw_error
                error('Cannot load template bank from file "%s"', f); 
            end

        end
        
        function loadFilterBankOort(obj, throw_error)

            if nargin<2 || isempty(throw_error)
                throw_error = 0;
            end
            
            frame_rate = floor(obj.head.FRAME_RATE); 
            
            obj.bank_oort = occult.ShuffleBank.empty; 
            obj.bank_oort_small = occult.ShuffleBank.empty; 
            
            for ii = 1:length(obj.pars.oort_distances)

                f = fullfile(getenv('DATA'), sprintf('WFAST/occultations/templates_%dAU_%dHz.mat', ...
                    floor(obj.pars.oort_distances(ii)), frame_rate));
                if exist(f, 'file')
                    load(f, 'bank');
                    obj.bank_oort(ii) = bank;
                elseif throw_error
                    error('Cannot load template bank from file "%s"', f); 
                end

                f = fullfile(getenv('DATA'), sprintf('WFAST/occultations/templates_%dAU_%dHz_small.mat', ...
                    floor(obj.pars.oort_distances(ii)), frame_rate));
                if exist(f, 'file')
                    load(f, 'bank');
                    obj.bank_oort_small(ii) = bank;
                elseif throw_error
                    error('Cannot load template bank from file "%s"', f); 
                end

            end
            
            for ii = 1:length(obj.bank_oort)
                obj.var_buf_oort(ii) = util.vec.CircularBuffer;
            end

        end
        
    end
    
    methods (Hidden=true) % simulations
        
        function star_indices = findStarsForSim(obj, star_indices) % take a list of stars from the prefilter (or all stars) and get a list of good stars for a simulated event
           
            stars_not_for_sim = false(1,size(obj.corrected_fluxes,2)); % logical vector with 1s where any star has a chance of already having an event
            stars_not_for_sim(star_indices) = true; % choose only stars that have a recent candidate
            
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
        
        function [all_cand, sim_pars] = simulateSingleEvent(obj, star_indices) % choose a star from the list of star_indices and add a randomly chosen occultation to it, then search for that event
            
            all_cand = tno.Candidate.empty;
            sim_pars = []; 
            
            if isempty(star_indices) % no stars, no simulation! 
                return;
            elseif isscalar(star_indices) % single star is given, just use it
                star_index_sim = star_indices;
            else
                star_index_sim = star_indices(randperm(length(star_indices),1)); % randomly choose a single star from the list
            end

            if obj.pars.use_oort
                all_banks = [obj.bank, obj.bank_oort]; % append all template banks, KBOs and Oort
            else
                all_banks = obj.bank; % only KBO banks
            end
            
            bank_idx = randi(length(all_banks)); % index of bank that was used to produce this event
            bank = all_banks(bank_idx); % choose a random bank to simulate from
            
            [flux_sim, detrend_sim, sim_pars] = obj.getFluxWithSimulatedLightcurve(star_index_sim, bank); % add the simulated occultation to the raw fluxes

            for ii = 1:length(all_banks) % check if event triggers on all filter banks
            
                bank = all_banks(ii); % which bank is used to filter
            
                [f_corr, std_corr] = obj.correctFluxes(detrend_sim, star_index_sim); % PSD or linear fit correction (on top of the simulated event!)

                f_filt = bank.input(f_corr(:,star_index_sim), std_corr(1,star_index_sim)); % filter only the star with the simulated event on it

                % get an estimate for the background of this flux:
                if obj.pars.use_std_filtered % need to correct the filtered fluxes by their measured noise

                    [bg_flux, bg_std] = obj.correctFluxes(obj.removeOutliers(obj.store.background_detrend(:,star_index_sim)), star_index_sim); % run correction on the background region
                    background_ff = bank.input(bg_flux, bg_std); % filter the background flux also
                    bg_ff_std = nanstd(background_ff); 

                    f_filt = f_filt./bg_ff_std; 

                end % need to correct the filtered fluxes by their measured noise

                [candidate, best_snr(ii)] = obj.searchForCandidates(flux_sim, detrend_sim, f_corr, f_filt, star_index_sim, all_banks(ii), 1); % try to find a single event on this single star

                all_cand = vertcat(all_cand, candidate);
                
            end
            
            N_triggers = length(all_cand); 
            
            trig_names = {all_cand.template_bank}; % append all the bank names that triggered
            
            if N_triggers>0 % we recovered the injected event! 
            
                [~,idx] = max([all_cand.snr]); % find the best candidate (highest S/N)    
                candidate = all_cand(idx); % get the detection S/N from the best candidate/template bank
                
                sim_pars.detect_snr = candidate.snr; % keep track of the detection S/N for this event
                sim_pars.passed = true; 
                
                sim_pars.num_triggers = N_triggers; % how many template banks triggered
                sim_pars.bank_names = trig_names; % what the names of those banks were
                
                for ii = 1:length(all_cand)
                    all_cand(ii).is_simulated = 1; 
                    all_cand(ii).sim_pars = sim_pars; 
                    all_cand(ii).applySimToExtraFlux(obj.store); 
                end
                
            else % no events were recovered
                
                r = sim_pars.r;
                R = sim_pars.R;
                b = sim_pars.b;
                
%                 if r>1.5 && r>R && r>b
%                     fprintf('failed to detect event with r= %4.2f | R= %4.2f | b= %4.2f\n', r, R, b); 
%                 end
                
                sim_pars.passed = false; 
%                 sim_pars.detect_snr = 0; % event was not detected, so S/N is zero
                sim_pars.detect_snr = best_snr(bank_idx); % this does not guarantee that the S/N is not from just random noise triggering at an unrelated filter/peak time
                
                sim_pars.num_triggers = 0; % how many template banks triggered
                sim_pars.bank_names = {}; % what the names of those banks were
                
            end

        end % we recovered the injected event!
        
        function [flux_all, detrend_all, sim_pars] = getFluxWithSimulatedLightcurve(obj, star_idx, bank) % take a flux matrix and an index for a star and add a random occultation event on top of it
        
            if isempty(obj.cat) || isempty(obj.cat.success) || obj.cat.success==0
                error('Cannot simulate events without stellar sizes in the catalog!'); 
            end

            %%%% get the occultation parameters %%%%
            
            R = obj.cat.data.FresnelSize;

            star_idx_global = obj.store.star_indices(star_idx);

            if isnan(R(star_idx_global)) % no stellar radius known from GAIA, try to estimate it 
                R_star = obj.estimateR(star_idx_global); 
            elseif R(star_idx_global)>obj.pars.sim_max_R % do not simulate stars bigger than this 
                R_star = obj.pars.sim_max_R; 
            else
                R_star = R(star_idx_global);
            end

            R_star = R_star*sqrt(bank.D_au./40); % adjust the stellar size in case the bank is for Hills/Oort cloud

            r_occulter = util.stat.power_law_dist(obj.pars.sim_power_law, 'min', bank.r_range(1), 'max', bank.r_range(2)); % occulter radius drawn from power law distribution
            
            b_par = bank.b_range(1) + rand*(diff(bank.b_range)); % impact parameter drawn from uniform distribution
            vel = bank.v_range(1) + rand*(diff(bank.v_range)); % velocity drawn from uniform distribution

            %%%% make the template %%%%
            
            [template, sim_pars] = bank.gen.randomLC('stellar_size', R_star, ...
                'occulter', r_occulter, 'impact', b_par, 'velocity', vel); % get a random occultation
            sim_pars.D = bank.D_au; 

            %%%% clean up the flux %%%% 
            
            flux_all = obj.store.extended_flux; % flux as it was measured
            detrend_all = obj.store.extended_detrend; % flux after removing linear trend
            
            flux_raw = flux_all(:,star_idx); % pick out the one star
            bg = obj.store.extended_aux(:,star_idx,obj.store.aux_indices.backgrounds).*...
                obj.store.extended_aux(:,star_idx,obj.store.aux_indices.areas);
            
            bg = nanmedian(bg); % prefer the median value to individual measurements, that could be outliers
            
            F = nanmean(flux_all(:,star_idx),1) - bg; % the mean flux, excluding the background
            
            t = (1:size(flux_all,1))'; % fake timestamps for the fit
            fr = util.fit.polyfit(t, detrend_all(:,star_idx), 'order', 2, 'double', 1); % fit the detrended flux to a second order polynomial
            
            flux_smoothed = fr.ym + F; % any residual 2nd order polynomial from the detrended flux, added to the mean flux level (this does not include the b/g)
            flux_noise = detrend_all(:,star_idx) - fr.ym; % separate the noise from the polynomial fit
            
            %%%% choose a good point to put the peak of the template %%%%
            
            B = obj.store.checker.bad_times; % bad times matrix
            
            good_times = obj.store.search_start_idx + find(~B(obj.store.search_start_idx:obj.store.search_end_idx,star_idx)); % mark regions of the star's flux that are not disqualified by the checker
            
            if isempty(good_times) 
                error('this shouldn''t happen!'); % we chose a star that has good times in its lightcurve! 
            end
            
            template = util.img.crop2size(template, size(flux_all,1)); 
            
            margin = size(flux_all,1) - length(template); % margins on both sides of the template, combined
            
            if margin>=0
                
                if length(good_times)>1 % choose a random point inside the good_times regions
                    peak_idx = randperm(length(good_times),1); 
                else
                    peak_idx = 1; % only one option, no need for random choice
                end
                
                peak = good_times(peak_idx); % position of the simulated event's center
                
                template_shift = vertcat(template, ones(margin,1)); % make the template long enough to fit the real lightcurve
                template_shift = circshift(template_shift, peak - floor(length(template)/2) - 1); % shift it until the middle of the template is at "peak"
                
            elseif margin<0
                error('need to handle this case at some point...'); % this would be a problem if we inject events with W>=8 seconds
            else
                error('need to handle this case at some point...'); 
            end
            
            %%%% add the template on the flux %%%%
            
            flux_smoothed_sim = flux_smoothed.*template_shift; % the mean flux is multiplied by the template (which is 1 outside the occulation region)
            
            if obj.pars.sim_rescale_noise
                
                background_var = nanmedian(obj.store.extended_aux(:,star_idx,obj.store.aux_indices.variances).*...
                    obj.store.extended_aux(:,star_idx,obj.store.aux_indices.areas)); % the read noise etc, which is independent of the source flux
                
                v1 = flux_smoothed + background_var; % variance taking into account only read-noise+bg and source noise before template is added
                v2 = flux_smoothed_sim + background_var; % variance taking into account only read-noise+bg and source noise after template is added
                
                flux_noise_corrected = (flux_noise - nanmean(flux_noise)).*sqrt(v2./v1) + nanmean(flux_noise); 
                
            end
            
            flux_final = flux_smoothed_sim + flux_noise_corrected + bg; % now we can re-attach the noise and the background we extracted before
            
            flux_all(:,star_idx) = flux_final; 
            detrend_all(:,star_idx) = flux_final - F - bg; % the detrended flux does not include the mean flux level!  
            detrend_all(:,star_idx) = obj.removeTrendByBatches(detrend_all(:,star_idx)); % remove trends after adding the simulated event (which is more realistic!)
            
            % add some parameters known from the simulation
            sim_pars.time_index = peak; % what was the real peak time of this event
            sim_pars.star_index = star_idx; % what was the star index
            sim_pars.star_snr = obj.store.star_snr(obj.store.star_indices(star_idx)); % translate to the global star index to get the S/N from the store
            sim_pars.calc_snr = sqrt(nansum((template-1).^2)).*sim_pars.star_snr; % estimate the event S/N with analytical formula (and assuming the star has white, gaussian noise)
            
            sim_pars.fluxes.raw_flux = flux_raw;
            sim_pars.fluxes.smoothed_flux = flux_smoothed;
            sim_pars.fluxes.sim_smoothed_flux = flux_smoothed_sim;
            sim_pars.fluxes.background_mean = bg;
            sim_pars.fluxes.background_var = background_var;
            sim_pars.fluxes.template = template_shift;
            sim_pars.fluxes.noise_flux = flux_noise;
            sim_pars.fluxes.noise_flux_corr = flux_noise_corrected; 
            sim_pars.fluxes.final_flux = flux_final; 
            
        end
        
        function val = estimateR(obj, idx) % get the star index and return an estimate for the star size (in FSU)
            
            S = obj.store.star_snr(idx); % the index must be given from the global list of stars! 
            
            if isempty(obj.store.size_snr_coeffs) || S<obj.store.pars.threshold % the second conditions should not occur... 
                val = util.stat.inverseSampling(obj.store.star_sizes, 'max', 3, 'width', 0.02);  
            else
                
                val = 0;
                
                for ii = 1:length(obj.store.size_snr_coeffs)
                    val = val + obj.store.size_snr_coeffs(ii).*S.^(ii-1); 
                end
                
            end

        end
        
    end
    
    methods % plotting tools / GUI
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = tno.gui.EvFinderGUI(obj); 
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
            
            obj.cat.addStellarSizes(input.distance_au);
            
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
            
            j = obj.store.extended_juldates(end) + double(obj.store.timestamps - obj.store.extended_timestamps(end))/24/3600; 
            a = celestial.coo.airmass(j, obj.head.RA, obj.head.DEC, [obj.head.longitude, obj.head.latitude]./180.*pi);
            
            t = datetime(j, 'convertFrom', 'juliandate'); 
            
            t2 = t(1+obj.store.pars.length_burn_in:end); 
            
%             w = obj.store.checker.mean_width_values.*2.355.*obj.head.SCALE;
%             w = obj.store.checker.defocus_log.*2.355.*obj.head.SCALE;
            
            b = obj.store.checker.mean_background_values;
            
            w = obj.store.fwhm_log;
            t3 = datetime(obj.store.juldates_log, 'ConvertFrom', 'JulianDate');
            
            
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
            
                plot(input.axes, t3, w, '-', 'DisplayName', 'seeing ["]', 'LineWidth', input.line); 
                
%                 plot(input.axes, binning(t2,obj.store.pars.length_search), w, '-', 'DisplayName', 'seeing ["]', 'LineWidth', input.line); 
                
%                 plot(input.axes, binning(t2, input.smooth*100), binning(w, input.smooth),'-k',...
%                     'DisplayName', 'seeing smoothed', 'LineWidth', input.line, 'HandleVisibility', 'off'); 

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
        
        function popupHours(obj)
            
            f = util.plot.FigHandler('Star hours'); 
            f.width = 30;
            f.height = 18;
            f.clear;
            
            obj.store.checker.hours.viewer('Parent', f.fig); 
            
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
            input.axes.YLim = [util.stat.min2(y)/2 util.stat.max2(y).*2];
            
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
        
        function popupSimEvents(obj)
            
            if isempty(obj.summary)
                obj.produceSummary;
            end
            
            f = util.plot.FigHandler('Simulated events'); 
            f.width = 30;
            f.height = 20;
            f.clear;
            
            obj.summary.showSimulations('parent', f.fig); 
            
        end
        
        function popupSimRates(obj)
            
            f = util.plot.FigHandler('Simulation detection rate');
            f.width = 30;
            f.height = 20;
            f.clear;
            
            obj.summary.showDetectionRate('parent', f.fig); 
            
        end
        
    end 
    
end

