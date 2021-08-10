classdef Candidate < handle
% A container for a single occultation event candidate. 
% This is one of the main outputs of the occultation pipeline. 
% Each candidate contains all data needed to make a decision if it is real
% or just an artefact / false positive. 
% The candidate object contains a history of the parameters of all objects
% used in the event finding pipeline, including the header data that also
% holds the photometric parameters (e.g., what aperture was used). 
% The objects used to keep this metadata are:
% "head", "kern_props", "star_props", "sim_pars", "finder_pars", 
% "store_pars", and "checker_pars". 
% 
% Other data saved with the candidate includes the time_index, star_index
% and kern_index (specifying the time, place, and best filter template) 
% that triggered the event. 
% There are fluxes, auxiliary data and cutouts saved for the time surrounding
% the event peak (the extended region as defined in the DataStore). 
% Additional parameters are described in the properties blocks, in the 
% inline comments.  
% 
% Some events are marked as kept=1 and some kept=0. The unkept events have
% been disqualified based on various conditions, mostly quality cuts and 
% black lists. See the EventFinder and QualityChecker objects for more info. 
% Each unkept candidate keeps a list of the cuts that triggered it's 
% removal, with the name and value of those cuts in "cut_hits", "cut_string"
% and "cut_value". In addition, the user can add notes to the event, and
% classify it into one of several categories (see getListOfClasses()). 
%
% Kept candidates should be vetted and classified by a human scanner. 
% Once they are classified they can be saved into some database for future
% reference and for applying statistical tools on multiple detections. 
% 
% Some candidates are simulated events that are injected into real data, and 
% are meant to be classified along with the real candidates. The idea is that
% during the classification / scanning phase, simulated events would be 
% indistinguishable from regular events, and a census of the number of such 
% events that were accepted/rejected by the scanners also helps quantify 
% the efficiency of the scanners. 
% Once classified, the simulated candidates should be skipped when running 
% statistics on real occultations. 
% 
% Plotting tools: 
% To make candidate vetting easier, there is a show() method that can be
% called on a vector of candidates, e.g., the output from the finder:
% >> finder.cand.show;
% This generates a simple GUI that shows all relevant data on the candidate, 
% and makes it easy to view, vet, and classify candidates. 
% The optional arguments to this function are:
% -index: which candidate from the vector given, should be plotted. This 
%         defaults to the previous index plotted in this figure/panel, but
%         will be set to 1 if the previous number doesn't exist or if it 
%         exceeds the number of candidates in the vector. 
% -kept: set the GUI to filter only kept events when pushing the prev/next
%        buttons to jump the index back/forward. 
% -cuts: add a list of cut values to display, on top of any triggered cuts. 
%        Triggered cuts are those that caused the candidate to be disqualified. 
%        Adding extra cuts helps figure out if some cut is borderline or if
%        the cut thresholds should be adjusted to handle more false positives. 
% -parent: the figure or panel into which we plot the different panels of 
%          this GUI. Default is gcf(). 
% -font_size: for axes labels etc. Default is 18. 
% 
% NOTE: the preferences for display are not saved in the candidate objects
%       themselves, but as a struct in the UserData field of the parent. 
%       Make sure not to overwrite this struct to maintain consistent plots. 
%       When closing the figure, it also resets these preferences. 
%
% The "reveal" button is only available if "use_show_secrets" is true. It 
% makes a pop-up window that shows the most important (and hidden) data
% about the event. This includes the ecliptic latitude, the shadow velocity
% and if the event is simulated. 
% To prevent biasing the scanners, we hide all this information. 
% E.g., a scanner might disqualify more events on high ecliptic latitude, 
% since KBOs are not expected to be found there, or disqualify events with 
% a velocity that doesn't match the observing direction. 
% When debugging is done, the unblinded sample should be scanned with the
% use_show_secrets=0. 

    properties % objects
        
        head; % header information
        kern_props; % struct with the data about the best kernel
        star_props; % table with one row from astrometry/GAIA
        sim_pars = struct; % struct with the simulation parameters (for simulated events only)

        finder_pars; % parameters defining the EventFinder object
        store_pars; % parameters defining the DataStore object
        checker_pars; % parameters defining the QualityChecker object
        
        mcmc_wide_v@occult.MCMC; % results of an MCMC run with a wide velocity prior
        mcmc_narrow_v@occult.MCMC; % results of an MCMC run with a narrow velocity prior around the expected star's velocity
        
    end
    
    properties % inputs/outputs
        
        serial; % each candidates's id in the run
        
        analysis_time = ''; % timestring when this event was identified by analysis (UTC)
        
        notes = {}; % anything added to the event, like failed cuts
        
        cutouts; % cutouts of the selected star, for the duration of the extended region
        kernel; % the lightcurve of the matched-filter kernel used 
        
        template_bank = ''; % note which template bank was used to trigger this event
        oort_template = 0; % was this candidate produced using an oort cloud template? 
        
        time_index; % out of the extended region
        kern_index; % from the full filter bank
        star_index; % from the stars that passed the initial burn-in (not from the subset that survived the pre-filter
        star_index_global; % from all the stars in this run (compare this to the catalog, cutouts, etc)
        star_snr; % from the store's burn in period
        frame_index; % the frame of the peak inside the original batch/file
        
        time_range; % time indices around time_index that are considered "part of the event"
        flux_mean; % mean of the raw flux (calculated over the background region)
        flux_std; % standard deviation of the raw flux (calculated over the background region)
        
        filenames; % cell array with the filenames for each frame in the extended region
        frame_numbers; % frame numbers for each frame inside its respective file
        
        % all of these are for the duration of the extended region
        timestamps; 
        juldates;
        flux_raw; % all fluxes for all stars, without any processing
        flux_detrended; % fluxes for all stars, after removing linear fit from each batch
        flux_corrected; % fluxes for all stars, corrected by PSD or by removing linear fit
        flux_filtered; % flux after matched-filtering, for the peak star and kernel
        aperture_index; % which aperture was used in the original photometry analysis (the index is for all types of photometry, including gaussian, as listed in head.PHOT_PARS)
        aperture_radius; % the size (in pixels) of the forced aperture chosen in the original analysis
        
        flux_unforced; % flux with aperture photometry used instead of forced, for the trigger star and nearest aperture to what was used in forced (loaded on demand)
        flux_unforced_all; % flux with aperture photometry used instead of forced, for all stars and multiple apertures (loaded on demand)
        unforced_aperture_sizes; % the radius of the apertures used, for the recalculated non-forced apertures, in units of pixels
        unforced_aperture_index; % which of the above apertures was chosen (closest match to the original, forced aperture size)
        unforced_aperture_radius; % the radius (in pixels) of the chosen aperture
        
        psd; % power spectral density for this star only
        freq_psd; % frequency values for the PSD
        
        auxiliary; % values of auxiliary measurements (e.g., background, offsets) for all stars and all frames in the extended region
        aux_names; % cell array with the names of each auxiliary
        aux_indices; % struct with a field for each auxiliary name, containing the index in the big matrix
        aux_unforced; % auxiliary data using aperture photometry for the specific star and aperture size (loaded on demand)
        aux_unforced_all; % auxiliary data using aperture photometry for all stars and apertures (loaded on demand)
        
        relative_dx; % the offsets_x for this star, minus the flux weighted mean offsets_x of all stars
        relative_dy; % the offsets_x for this star, minus the flux weighted mean offsets_x of all stars
        
        cut_matrix; % the different cuts from the QualityChecker (for this star, on the extended region)
        cut_names; % a cell array with the names of the cuts, e.g., 'shakes', 'offset_size', 'corr_x_25',...
        cut_indices; % a struct with each cut name as a field containing the index of that cut in the matrix's 2nd dim
        cut_thresh; % which threshold is used for each cut (a row vector)
        cut_two_sided; % true when the cut value should be compared in abs() to the threshold (a row vector)
        
        cut_string = {}; % cell string describing the cuts that this event hit
        cut_hits = []; % which cuts coincide with the event peak
        cut_value = []; % what is the value that triggered the cut (maximum within range of the peak)
        
        corr_flux = []; % correlation values, corrected for length of event region (normalized to 100 frames), for a few stars with the highest correlation
        corr_flux_stars = []; % indices of stars with highest flux correlation
        
        batch_number;
        
        snr;
        
    end
    
    properties % switches/controls
        
        kept = 1; % by default any new event is considered real (kept)
        
        classification = ''; 
        
        display_cuts = []; % user can add more cuts to display even if the event did not happen to cross that cut
        
        debug_bit = 1;
        
    end
    
    properties(Hidden=true)
        
        use_show_secrets = false; % this should be turned on only when debuggin (on the un-blinded set)
        
        run_identifier = ''; % use this to find the files in the database. Format: <date folder>\<run_name_folder> e.g., 2020-09-07\ecliptic_run1

        folder = ''; % use this when loading the candidates from file, so you know where to save them to after classification is done
        
        is_positive; % if the filtered flux peak is positive. If not, it means the flux is opposite the template (e.g., a flare with an occultation template)
        
        is_simulated = 0; % by default candidates are not simulated. If true, this candidate is an injected simulated event.  

        flux_buffer; % flux buffer going back as far as possible, for this star only
        detrend_buffer; % flux after removing a linear fit from each batch
        timestamps_buffer; % timestamps for the above
%         filtered_flux_past_values; % normalization data for this star and this kernel only. Can be either raw filtered_flux values for the last "backround" period, or the variance values from the var_buf
        
        flux_raw_all; % the raw flux (over the extended region) for all stars
        flux_detrended_all; % the flux after removing a linear fit from each batch, for all stars
        flux_corrected_all; % the corrected flux (over the exteneded region) for all stars
        auxiliary_all; % the auxiliary data (for the extended region) for all stars
        cutouts_all; % the cutouts for all stars (in the extended region)
        global_star_indices; % the indices of the stars chosen by the DataStore to be above some threshold
        
        search_start_idx; % the frame index inside the extended region, where the search region starts (often this is 51)
        search_end_idx; % the frame index inside the extended region, where the search region ends (often this is 150)
        
        kern_extra; % any other kernels that passed the lower threshold for kernels
        star_extra; % any other stars that passed the lower threshold for stars
        
        version = 1.02;
        
    end
    
    properties(Hidden=true, Transient=true) % lazy loaded heavy stuff like the stack image
        
        stack; 
        calibration;
        
    end
    
    methods % constructor
        
        function obj = Candidate(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.Candidate')
                if obj.debug_bit>1, fprintf('Candidate copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('Candidate constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % getters
        
        function classes = getListOfClasses(obj)
            
            classes = {'occultation certain', 'occultation possible', 'satellite', 'flare', 'cosmic ray', 'bad pixel', 'edge effect', 'tracking error', 'artefact', 'mystery'};
            
        end
        
        function val = kernel_timestamps(obj)
            
            dt = median(diff(obj.timestamps));
            
            val = (-floor(length(obj.kernel)/2):floor(length(obj.kernel)/2)).*dt;
            
        end
        
        function val = kernel_lightcurve(obj, flip)
            
            if nargin<2 || isempty(flip)
                flip = 0;
            end
            
            if isempty(obj.kernel) || isempty(obj.kern_props)
                val = [];
            else
                val = (-1).^flip.*obj.kernel.*obj.kern_props.norm + 1; 
            end
            
        end
        
        function val = star_current_snr(obj)
            
            if isempty(obj.flux_mean) || isempty(obj.flux_std)
                val = [];
            else
                val = obj.flux_mean./obj.flux_std;
            end
            
        end
        
        function val = time_start(obj)
            
            val = obj.timestamps(obj.search_start_idx);
            
        end
        
        function val = time_end(obj)
            
            val = obj.timestamps(obj.search_end_idx);
            
        end
        
        function val = getFilename(obj)
            
            if isempty(obj.filenames) || isempty(obj.time_index)
                val = '';
            else
                val = obj.filenames{obj.time_index}; 
            end
            
        end
        
        function str = printSummary(obj)
            
            J = mean(obj.juldates); 
            
            t = datetime(J, 'convertFrom', 'juliandate', 'TimeZone', 'UTC'); 
            
            str = sprintf('id: %d | star: %d | frame= %d | batch= %d | time: %02d:%02d:%02d | event S/N= %4.2f | star S/N= %4.2f (%4.2f) | FWHM= %4.1f" | x/y= %d, %d ', ...
                    obj.serial, obj.star_index, obj.time_index, obj.batch_number, hour(t), minute(t), round(second(t)), obj.snr, obj.star_current_snr, obj.star_snr, ...
                    obj.fwhm, round(nanmean(obj.auxiliary(:,obj.aux_indices.centroids_x))), round(nanmean(obj.auxiliary(:,obj.aux_indices.centroids_y))));
            
        end
        
        function val = fwhm(obj)
            
            if isempty(obj.cut_matrix) || isempty(obj.cut_indices) || ~isfield(obj.cut_indices, 'fwhm')
                val = [];
            else
                val = obj.cut_matrix(1,obj.cut_indices.fwhm);
            end
            
        end
        
        function str = printNotes(obj)
            
            str = {};
            
            if ~isempty(obj.cut_string)
                str = vertcat(str, util.vec.tocolumn(obj.cut_string));
            end
            
            if ~isempty(obj.classification)
                str = vertcat(str, {sprintf('classified as "%s"', obj.classification)}); 
            end
            
            if ~isempty(obj.notes)
                str = vertcat(str, util.vec.tocolumn(obj.notes));
            end
            
            str = strjoin(str, newline); 
            
        end
        
        function str = printKernelProps(obj)
            
            if isempty(obj.kern_props)
                str = '';
            else
                str = sprintf('R= %4.2f | r= %4.2f | b= %4.2f | v= %4.1f | D= %d | invert= %d', ...
                    obj.kern_props.R, obj.kern_props.r, obj.kern_props.b, obj.kern_props.v, obj.kern_props.D, ~obj.is_positive);
            end 
            
        end
        
        function str = printStellarProps(obj)
            
            if isempty(obj.star_props)
                str = '';
            else
                str = sprintf('%4.2fmag, %dK, %3.1f FSU', obj.star_props.Mag_BP, round(obj.star_props.Teff), obj.getStellarSize); 
            end
            
        end
        
        function val = getStellarSize(obj)
            
            if ~isempty(obj.star_props) && ismember('Mag_BP', obj.star_props.Properties.VariableNames) && ~isnan(obj.star_props.Mag_BP) && ...
                    ismember('Teff', obj.star_props.Properties.VariableNames) && ~isnan(obj.star_props.Teff) % catalog entry exists and has non NaN mag/temperature
                
                s = util.ast.star_sizes_gaia(obj.star_props); 
                
                val = util.ast.fresnel_size(s, 40, 500); 
                
%                 if ~ismember('BolMag', obj.star_props.Properties.VariableNames) || isnan(obj.star_props.BolMag) || ...
%                     ~ismember('BolTemp', obj.star_props.Properties.VariableNames) || ~isnan(obj.star_props.BolTemp) % no bolometric data, need to calculate it now
%                 
%                     [C,T] = util.ast.bol_corr(obj.star_props.Teff, obj.star_props.Mag_BP, obj.star_props.Mag_RP, 'BP', 'RP', 'GAIA', 0);
%                 
%                     bol_mag = obj.star_props.Mag_BP + C; 
%                 
%                     obj.star_props.BolMag = bol_mag;
%                     obj.star_props.BolTemp = T; 
%                 end
%                 
%                 val = util.ast.stellar_size(obj.star_props.BolMag, obj.star_props.BolTemp, 'units', 'FSU', 'dist', 40/206000); % distance is given in parsec (40AU)
                
            else
                val = NaN; 
            end
            
        end
        
        function str = printPhotometricPars(obj)
            
            if isempty(obj.head) || isempty(obj.head.PHOT_PARS)
                str = '';
            else
                str = sprintf(); % fill this at some point? 
            end
            
        end
        
        function str = printRunData(obj)
           
            RA = ''; 
            Dec = '';
            ECL = [];
            
            if ~isempty(obj.head)
                RA = obj.head.OBSRA;
                Dec = obj.head.OBSDEC; 
                ECL = obj.head.ephem.ECL_lat;
            end
            
            v = obj.head.ephem.getShadowVelocity;
            v = sqrt(sum(v.^2)); 
            
            str = sprintf('Run: "%s" | ecl lat= %4.1f deg \ncoords: %s %s | V= %d', strrep(obj.run_identifier, '\', '/'), ECL, RA, Dec, round(v)); 
            
        end
        
        function str = printSimulationData(obj)
            
            if ~obj.is_simulated
                str = '';
            else
                str = sprintf('SIM: R= %4.2f | r= %4.2f | b= %4.2f | v= %4.1f', ...
                    obj.sim_pars.R, obj.sim_pars.r, obj.sim_pars.b, obj.sim_pars.v); 
            end
            
        end
        
        function dir = findAnalysisFolder(obj)
            
            d = util.sys.WorkingDirectory(fullfile(getenv('DATA'), 'WFAST', obj.run_identifier(1:4), obj.run_identifier)); % the run folder
            
            if ~exist(d.pwd, 'dir')
                error('Could not find folder "%s"', d.pwd); 
            end
            
            list = sort(d.match_folders('analysis*'));
            
            if isempty(list)
                error('Could not find any analysis folders in "%s"', d.pwd); 
            end
            
            dir = list{end}; 
            
        end
        
        function summary = getSummary(obj)
            
            d = obj.findAnalysisFolder; 
            
            load(fullfile(d, 'summary.mat')); 
            
            if ~exist('summary', 'var')
                summary = [];
            end
            
        end
        
        function flux = getNormalizedFlux(obj, margin) % get a centered lightcurve, normalized to 1, with NaNs outside the event region (+margin)
            
            if nargin<2 || isempty(margin)
                margin = 5;
            end
            
            if isempty(obj.flux_raw)
                flux = [];
            else
                flux = obj.flux_raw; 
                flux = flux./obj.flux_mean; 
                idx = false(size(flux)); 
                idx(obj.time_range) = true; 
                idx = imdilate(idx, ones(margin)); 
                flux(~idx) = NaN; 
                flux = circshift(flux, floor(size(flux,1)/2) + 1 - obj.time_index); 
            end
            
        end
        
        function val = getEclipticLatitutde(obj)
            
            if isempty(obj.star_props) || ~isa(obj.star_props, 'table') ...
                    || ~ismember('RA', obj.star_props.Properties.VariableNames) ...
                    || ~ismember('Dec', obj.star_props.Properties.VariableNames) 
                val = [];
            else
                out_coord = celestial.coo.coco([obj.star_props.RA, obj.star_props.Dec], 'J2000', 'e', 'd', 'd');
                val = out_coord(2); 
            end
            
        end
        
        function val = getShadowVelocity(obj)
            
            if isempty(obj.star_props) || ~isa(obj.star_props, 'table') ...
                    || ~ismember('RA', obj.star_props.Properties.VariableNames) ...
                    || ~ismember('Dec', obj.star_props.Properties.VariableNames) 
                val = [];
            else
                e = head.Ephemeris; 
                e.RA_deg = obj.star_props.RA;
                e.Dec_deg = obj.star_props.Dec;
                val = e.getShadowVelocity; 
                val = sqrt(sum(val.^2));                 
            end
            
        end
        
    end
    
    methods % setters
        
        function addNote(obj, str) % free form string comment by users
            
            if isempty(obj.notes)
                obj.notes{1} = str;
            else
                obj.notes{end+1} = str;
            end
            
        end
        
    end
    
    methods % calculations
        
        function checkIfPeakIsIncluded(obj) % check if the time_index is really the highest point in the filtered flux (return false if the real peak is outside the search region)
            
            f = obj.flux_filtered;
            
            idx = obj.time_index; 
            
            mx = 0; % best maximum outside the search region
            real_peak_idx = []; 
            
            % go forward
            for ii = 1:length(f)
                
                idx = idx + 1;
                
                if idx>length(f)
                    break; % out of bounds!
                end
                
                if abs(f(idx))<obj.finder_pars.time_range_thresh
                    break; % reached the end of the event on this side
                end
                
                if idx>obj.search_end_idx % we are now outside the search region
                
                    if abs(f(idx))>abs(mx)
                        mx = f(idx);
                        real_peak_idx = idx; 
                    end
                    
                end
                
            end
            
            idx = obj.time_index;
             
            % go backwards
            for ii = 1:length(f)
                
                idx = idx - 1;
                
                if idx<1
                    break; % out of bounds!
                end
                
                if abs(f(idx))<obj.finder_pars.time_range_thresh
                    break; % reached the end of the event on this side
                end
                
                if idx<obj.search_start_idx % we are now outside the search region
                
                    if abs(f(idx))>abs(mx)
                        mx = f(idx);
                        real_peak_idx = idx; 
                    end
                    
                end
                
            end
            
            if abs(mx)>abs(f(obj.time_index))
                obj.notes{end+1} = sprintf('Peak outside, frame= %d, S/N= %4.2f', real_peak_idx, mx); 
                obj.kept = 0;
            end
            
        end
        
        function [star_idx, corr] = findHighestCorrelations(obj, number, frames)
        % check the correlations to other stars' fluxes, with +-"frames"
        % around the time_index, and give the indices of the most
        % correlated stars (find "number" such matches). 
        
            if isempty(obj.flux_raw_all)
                error('Must have "flux_raw_all" to calculate correlation with other star fluxes');
            end
        
            if nargin<2 || isempty(number)
                number = 3; 
            end
            
            if nargin<3 || isempty(frames)
                frames = 2*length(obj.time_range); 
            end
            
            frame_idx = (-frames:frames) + obj.time_index; % plus/minus number of frames around peak
            frame_idx(frame_idx<1 | frame_idx>size(obj.flux_raw,1)) = []; 
            f = obj.flux_raw(frame_idx);
            F = obj.flux_raw_all(frame_idx,:); 
            
            % normalize both fluxes
            f = (f-nanmean(f)); 
            F = (F-nanmean(F)); 
%             
%             C = nansum(f.*F)./length(idx); % normalized fluxes just give the correlation (up to number of samples)
            
            C = nansum(f.*F)./sqrt(nansum(f.^2).*nansum(F.^2)).*sqrt(length(frame_idx)); 
            C(isnan(C)) = 0; 
            [~, sort_idx] = sort(C); % sort, and get a list of the sorted indices
            
            star_idx = flip(sort_idx(end-number:end-1));
            
            corr = C(:,star_idx); 
            
        end
        
        function val = calcTrackingErrorValue(obj, N_stars)
            
            if nargin<2 || isempty(N_stars)
                N_stars = 3; 
            end
            
%             N_frames = 2*length(obj.time_range);
            
            if isempty(obj.corr_flux)
                [obj.corr_flux_stars, obj.corr_flux] = obj.findHighestCorrelations(N_stars);
%                 obj.corr_flux = corr.*sqrt(N_frames)./10; % normalize by number of frames relative to 100 frames
            end
            
            val = sum(obj.corr_flux); 
            
        end
        
        function clearExtraData(obj_vec) % remove large arrays (e.g., before saving to disk)
            
            for ii = 1:length(obj_vec)

                obj_vec(ii).flux_raw_all = [];
                obj_vec(ii).flux_detrended_all = []; 
                obj_vec(ii).flux_corrected_all = [];
                obj_vec(ii).auxiliary_all = [];
                obj_vec(ii).cutouts_all = [];

            end
            
        end
        
        function saveClassified(obj_vec, filename)
            
            if nargin<2 || isempty(filename)
                filename = 'classified.mat';
            end
            
            candidates = obj_vec(~cellfun(@isempty, {obj_vec.classification}));
            
            candidates = util.oop.full_copy(candidates); 
            
            candidates.clearExtraData;
            
            save(filename, 'candidates', '-v7.3'); 
            
        end
        
        function val = isSameEvent(obj, other)
            
            val = 0; % assume the candidates are not the same underlying event
            
            if obj.star_index~=other.star_index, return; end % must be on the same star
            
            if obj.batch_number~=other.batch_number, return; end % must be in the same batch
            
            if abs(obj.frame_index-other.frame_index)>10, return; end % must have peaks at least 10 frames from each other (allow for kernel shape to move the peak)
            
            val =1; % if we didn't short circuit anywhere, the events are the same
            
        end
        
        function runMCMC(obj, varargin)
            
            input = util.text.InputVars; 
            input.input_var('velocity_prior', 'wide'); % can also choose "narrow"
            input.input_var('chains', 10, 'num_chains', 'number_chains'); 
            input.input_var('plot', false); % load up an MCMC GUI for the run
            input.scan_vars(varargin{:}); 
            
            if util.text.cs(input.velocity_prior, 'wide')
                use_wide = 1;
                wide_str = 'wide'; 
            elseif util.text.cs(input.velocity_prior, 'narrow')
                use_wide = 0;
                wide_str = 'narrow'; 
            else
                error('Unknown velocity prior option "%s". Use "wide" or "narrow"...', input.velocity_prior);
            end
            
            mcmc = occult.MCMC; 
            mcmc.input_flux = obj.flux_raw./obj.flux_mean;
            mcmc.input_flux = circshift(mcmc.input_flux, 100-obj.time_index); % center the peak
            mcmc.input_errors = obj.flux_std./obj.flux_mean; 
            mcmc.input_R = obj.star_props.FresnelSize; 
            mcmc.input_v = obj.getShadowVelocity; 
            
            mcmc.setupDeepScan; 
            mcmc.num_chains = input.chains; 
            
            if use_wide
                mcmc.setupWideVelocityPrior;
            else
                mcmc.setupNarrowVelocityPrior; 
            end
            
            if input.plot
                mcmc.makeGUI;
                drawnow;
            end
            
            if obj.debug_bit, fprintf('Running MCMC with %s velocity prior.\n', wide_str); end
            
            mcmc.run;
            
            if use_wide
                obj.mcmc_wide_v = mcmc;
            else
                obj.mcmc_narrow_v = mcmc; 
            end
            
        end
        
        function loadUnforcedPhotometry(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('apertures', obj.head.PHOT_PARS.aperture_radius, 'radius', 'radii'); 
            input.input_var('annulus', obj.head.PHOT_PARS.annulus_radii)
            input.input_var('folder', [], 'root'); 
            input.scan_vars(varargin{:}); 
            
            %%% verify the filenames with the right path %%%
            
            f = unique(obj.filenames); % get the unique set of filenames where the data is saved
            
            if isempty(input.folder) 
                
                if ~exist(fileparts(f{1}), 'dir')
                    error('Cannot find folder "%s". Use the "folder" optional input to specify the root folder', fileparts(f{1})); 
                end
                
            else % the base folder is not the one in the saved filenames (e.g., the root folder was moved)
               
                if ~ischar(input.folder) || ~exist(input.folder, 'dir')
                    error('Illegal input to "folder" option. Must be a string with a valid path'); 
                end
                    
                for ii = 1:length(f)
                    
                    [path, file_name, ext] = fileparts(f{ii}); 
                    [path, run_name] = fileparts(path); 
                    [path, run_date] = fileparts(path); 
                    
                    f{ii} = fullfile(input.folder, run_date, run_name, [file_name ext]); 
                    
                end
                
            end
            
            if obj.version<1.02 % addresses bug where we saved the wrong filenames (for file numbers N-1 and N we saved numbers N and N+1)
                for ii = 1:length(f)-1
                    f{ii+1} = f{ii}; 
                end
                
                d = util.sys.WorkingDirectory(fileparts(f{1})); 
                file_list = d.files([], 1); % get the files with the full path
                idx = find(strcmp(file_list, f{1})); 
                f{1} = file_list{idx-1}; 
            end
            
            
            %%% figure out the original aperture %%% 
            
            if isempty(obj.global_star_indices) || isempty(obj.aperture_index) || isempty(obj.aperture_radius) || isempty(obj.store_pars) || ~isstruct(obj.store_pars)
                % reload the store pars from the summary file
                t = datetime(obj.analysis_time(1:10), 'Format', 'yyyy-MM-dd'); 
                base_dir = fileparts(f{end});
                
                for ii = 1:100
                    
                    test_dir = fullfile(base_dir, sprintf('analysis_%s', t));
                    
                    if exist(test_dir, 'dir')
                        break;
                    end
                    
                    t = t - 1; % remove one day
                    
                end
                
                if ~exist(test_dir, 'dir')
                    error('Could not load analysis folder for recovering the "store_pars"'); 
                end
                    
                finder_filename = fullfile(test_dir, 'finder.mat');
                load(finder_filename); 
                obj.store_pars = finder.store.pars; 
                obj.store_pars.aperture_index = finder.store.aperture_index;
                obj.store_pars.star_sizes = finder.store.star_sizes;
                obj.store_pars.star_snrs = finder.store.star_snr;
                
                obj.global_star_indices = finder.store.star_indices; 
                obj.aperture_index = finder.store.aperture_index;
                obj.aperture_radius = util.text.extract_numbers(obj.head.PHOT_PARS.types{obj.aperture_index}); 
                obj.aperture_radius = obj.aperture_radius{1}; % extract_numbers() gives a cell array
                
            end
            
            
            %%% now load the data and run the photometry %%%
            
            if isempty(obj.calibration)
                obj.calibration = img.Calibration; 
                obj.calibration.loadByDate(obj.head.RUNSTART(1:10), obj.head.INST, obj.head.PROJECT); 
            end
            
            F = []; % raw flux
            A = []; % auxiliary data
            
            for ii = 1:length(f) % go over files
                
                C = h5read(f{ii}, '/cutouts'); 
                P = h5read(f{ii}, '/positions'); 
                C = obj.calibration.input(C, 'pos', P); 
                
                if isempty(input.apertures) || isempty(input.annulus)
                    error('Could not find the original aperture or annulus radii'); 
                end
                
                s = util.img.photometry2(C, 'use_aperture', 1, 'resolution', 2, 'annulus', input.annulus, 'aperture', input.apertures); 
                S = s.apertures_photometry; 
%                 S = s.forced_photometry; % debug only! 
                
                F = vertcat(F, S.flux); 
                
                A_new = [];
                for jj = 1:length(obj.aux_names)
                    if util.text.cs('offsets_', obj.aux_names{jj})
                        A_new = cat(3, A_new, permute(S.(['offset_' obj.aux_names{jj}(end)]), [1,2,4,3])); 
                    elseif strcmp(obj.aux_names{jj}, 'centroids_x')
                        A_new = cat(3, A_new, P(:,1)' + permute(S.offset_x, [1,2,4,3])); 
                    elseif strcmp(obj.aux_names{jj}, 'centroids_y')
                        A_new = cat(3, A_new, P(:,2)' + permute(S.offset_y, [1,2,4,3])); 
                    elseif isfield(S, obj.aux_names{jj})
                        A_new = cat(3, A_new, permute(S.(obj.aux_names{jj}), [1, 2, 4, 3]));
                    elseif isfield(S, obj.aux_names{jj}(1:end-1))
                        A_new = cat(3, A_new, permute(S.(obj.aux_names{jj}(1:end-1)), [1, 2, 4, 3]));                        
                    elseif isfield(S, obj.aux_names{jj}(1:end-2))
                        A_new = cat(3, A_new, permute(S.(obj.aux_names{jj}(1:end-2)), [1, 2, 4, 3]));
                    else
                        error('Could not match name "%s" with any fields in the photometry output', obj.aux_names{jj}); 
                    end
                end
                
                A = vertcat(A, A_new); 
                
            end
            
            new_radii = s.parameters.aperture_radius; 
            obj.unforced_aperture_sizes = new_radii; 
            
            [~, obj.unforced_aperture_index] = min(abs(obj.aperture_radius - new_radii)); % find the closest aperture to what was used in the forced photometry
            
            obj.flux_unforced_all = F; 
            obj.flux_unforced = F(:,obj.star_index_global, obj.unforced_aperture_index); 
            
            obj.aux_unforced_all = A; 
            obj.aux_unforced = A(:,obj.star_index_global, :, obj.unforced_aperture_index); 
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj_vec, varargin) % generate a GUI interface for a vector of Candidate objects
            
            input = util.text.InputVars;
            input.input_var('index', []); % which event in the vector to plot (default is 1 or what was plotted previously)
            input.input_var('kept', [], 'show_kept'); % if true, filter only the kept events when pressing the prev/next buttons
            input.input_var('cuts', []); % additional cut types to show on the flux plot
            input.input_var('duplicates', []); % skip duplicate events when classifying
            input.input_var('parent', []); % parent graphic object to plot to (figure or panel, default is gcf())
            input.input_var('font_size', 18);
            input.scan_vars(varargin{:});
            
            if isempty(input.parent)
                input.parent = gcf;
            end
            
            if isempty(input.index) % default value
                if ~isempty(input.parent.UserData) && isfield(input.parent.UserData, 'index')
                    input.index = input.parent.UserData.index; 
                else
                    input.index = 1;
                end
            end
            
            if isempty(input.kept) % default value
                if ~isempty(input.parent.UserData) && isfield(input.parent.UserData, 'show_kept')
                    input.kept = input.parent.UserData.show_kept; 
                else
                    input.kept = false;
                end
            end
            
            if isempty(input.cuts) % default value
                if ~isempty(input.parent.UserData) && isfield(input.parent.UserData, 'show_cuts')
                    input.cuts = input.parent.UserData.show_cuts; 
                else
                    input.cuts = false;
                end
            end
            
            if isempty(input.duplicates) % default value
                if ~isempty(input.parent.UserData) && isfield(input.parent.UserData, 'skip_duplicates')
                    input.duplicates = input.parent.UserData.skip_duplicates; 
                else
                    input.duplicates = true;
                end
            end
            
            if isempty(input.parent.UserData)
                input.parent.UserData = struct('number', length(obj_vec), 'index', input.index, ...
                    'show_kept', input.kept, 'show_cuts', input.cuts, 'skip_duplicates', input.duplicates); % add other state variables here
            else
                input.parent.UserData.number = length(obj_vec); 
                input.parent.UserData.index = input.index;
                input.parent.UserData.show_kept = input.kept;
                input.parent.UserData.show_cuts = input.cuts; 
                % and add them here too
            end
            
            idx = input.parent.UserData.index;
            if idx>length(obj_vec)
                input.parent.UserData.index = 1;
                idx = 1;
            end
            
            input.parent.UserData.num_kept = nnz([obj_vec.kept]);
            
            obj = obj_vec(idx); 
            
            delete(input.parent.Children);
            
            margin_left = 0.05;
            
            ax1 = axes('Parent', input.parent, 'Position', [margin_left 0.55 0.55 0.35]);
            obj.showRawFlux('ax', ax1, 'on_top', 1, 'equal', 0, 'title', 0);
            
            ax2 = axes('Parent', input.parent, 'Position', [margin_left 0.2 0.55 0.35]);
            obj.showFilteredFlux('ax', ax2, 'title', 0, 'cuts', input.cuts); 
            
            ax3 = axes('Parent', input.parent, 'Position', [0.68 0.2 0.3 0.5]);
            obj.showCutouts('ax', ax3);
            
            %%%%%%%%%%%%%%%%%%%%%% popup panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            popup_panel = uipanel(input.parent, 'Units', 'Normalized', 'Position', [0.68 0.11 0.3 0.08], 'title', 'popups'); 
            
            button = uicontrol(popup_panel, 'Style', 'pushbutton', 'string', 'other stars', ...
                'Units', 'Normalized', 'Position', [0.0 0.0 0.25 1], 'FontSize', 14, ...
                'Callback', @obj.popupOtherStars, 'UserData', input.parent, ...
                'Tooltip', 'Show the flux of other stars with high correlation to this star'); 
            
            % turn button red if there are high correlations! 
            if length(obj.corr_flux)>=3 && obj.corr_flux(3)>obj.checker_pars.thresh_tracking_error
                button.BackgroundColor = 'r'; 
            end
            
            button = uicontrol(popup_panel, 'Style', 'pushbutton', 'string', 'stack image', ...
                'Units', 'Normalized', 'Position', [0.25 0.0 0.25 1], 'FontSize', 14, ...
                'Callback', @obj.popupStack, 'UserData', input.parent, ...
                'Tooltip', 'Show the stack image for the event'); 
            
            button = uicontrol(popup_panel, 'Style', 'pushbutton', 'string', 'cutouts', ...
                'Units', 'Normalized', 'Position', [0.50 0.0 0.25 1], 'FontSize', 14, ...
                'Callback', @obj.popupCutouts, 'UserData', input.parent, ...
                'Tooltip', 'show more cutouts around the event center'); 
            
            button = uicontrol(popup_panel, 'Style', 'pushbutton', 'string', '', ...
                'Units', 'Normalized', 'Position', [0.75 0.0 0.25 1], 'FontSize', 14, ...
                'Callback', '', 'UserData', input.parent, ...
                'Tooltip', 'placeholder'); 
            
            %%%%%%%%%%%%%%%%%%%%%% panel info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            info_panel = uipanel(input.parent, 'Units', 'Normalized', 'Position', [margin_left 0.9 1-margin_left*2 0.1], 'title', 'info'); 
            
            button = uicontrol(info_panel, 'Style', 'text', 'string', strjoin({obj.printSummary, ['Star: ' obj.printStellarProps ' | Kernel: ' obj.printKernelProps]}, newline), ...
                'Units', 'Normalized', 'Position', [0.02 0 0.88 1], 'FontSize', 14, 'HorizontalAlignment', 'Left'); 
            
            if obj.use_show_secrets
                button = uicontrol(info_panel, 'Style', 'pushbutton', 'string', 'reveal', ...
                    'Units', 'Normalized', 'Position', [0.9 0.1 0.1 0.9], 'FontSize', 14,...
                    'Callback', @obj.callback_show_secrets, 'UserData', input.parent, ...
                    'Tooltip', 'Show the run name, coordinates, and simulation status for this event');
            end
            
            %%%%%%%%%%%%%%%%%%%%%% panel notes %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            notes_panel = uipanel(input.parent, 'Units', 'Normalized', 'Position', [0.68 0.7 0.3 0.2], 'Title', 'notes'); 
            
            button = uicontrol(notes_panel, 'Style', 'text', 'string', obj.printNotes, ...
                'Units', 'Normalized', 'Position', [0.05 0 0.9 1], 'FontSize', 12, 'HorizontalAlignment', 'Left'); 
            
            
            %%%%%%%%%%%%%%%%%%%%%% panel cuts %%%%%%%%%%%%%%%%%%%%%%%%%
            
            top_margins = 0.05; 
            pos = margin_left;
            
            cuts_panel = uipanel(input.parent, 'Units', 'Normalized', 'Position', [pos 0.0 0.2 0.1], 'title', 'cuts'); 
            
            if input.cuts
                str = 'with cuts';
            else
                str = 'no cuts';
            end
            
            button = uicontrol(cuts_panel, 'Style', 'pushbutton', 'String', str, ...
                'Units', 'Normalized', 'Position', [0.05 top_margins 0.25 1-2*top_margins], 'FontSize', 12, ...
                'Callback', @obj_vec.callback_cuts, 'UserData', input.parent);
            
            if input.cuts
                button.TooltipString = 'show cuts on top of the filtered flux';
            else
                button.TooltipString = 'show only the filtered flux without cuts'; 
            end
            
            button = uicontrol(cuts_panel, 'Style', 'edit', 'String', ['Cuts= ' util.text.print_vec([obj.cut_hits, obj.display_cuts], ', ')], 'FontSize', 18, ...
                'Units', 'Normalized', 'Position', [0.35 top_margins 0.6 1-2*top_margins], ...
                'Callback', @obj_vec.callback_parse_cuts, 'UserData', input.parent);
            
            button.TooltipString = sprintf('Showing the cuts for: %s', strjoin(obj.cut_names([obj.cut_hits, obj.display_cuts]), ', '));  
            
            %%%%%%%%%%%%%%%%%%%%%% panel indexing %%%%%%%%%%%%%%%%%%%%%%%%%
            
            pos = pos + 0.05 + cuts_panel.Position(3); 
            
            index_panel = uipanel(input.parent, 'Units', 'Normalized', 'Position', [pos 0.0 0.4 0.1], 'title', 'controls'); 
            
            if input.parent.UserData.show_kept
%                 kept_string = sprintf('kept (%d)', input.parent.UserData.num_kept);
                kept_string = 'kept'; 
            else
                kept_string = 'all'; 
            end
            
            button = uicontrol(index_panel, 'Style', 'pushbutton', 'String', kept_string, 'FontSize', 20, ...
                'Units', 'Normalized', 'Position', [0.01 top_margins 0.18 1-2*top_margins], ...
                'Callback', @obj_vec.callback_kept, 'UserData', input.parent); 
            
            if input.parent.UserData.show_kept
                button.TooltipString = 'Showing only the kept candidates';
            else
                button.TooltipString = 'Showing all candidates';
            end
            
            button = uicontrol(index_panel, 'Style', 'pushbutton', 'String', '-', 'FontSize', 20, ...
                'Units', 'Normalized', 'Position', [0.2 top_margins 0.18 1-2*top_margins], ...
                'Callback', @obj_vec.callback_prev, 'UserData', input.parent, ...
                'Tooltip', 'go back to previous candidate'); 
            
            button = uicontrol(index_panel, 'Style', 'edit', 'String', sprintf('index= %d / %d', input.parent.UserData.index, length(obj_vec)), 'FontSize', 16, ...
                'Units', 'Normalized', 'Position', [0.39 top_margins 0.4 1-2*top_margins], ...
                'Callback', @obj_vec.callback_index, 'UserData', input.parent, ...
                'Tooltip', 'jump to any candidate index'); 
            
            button = uicontrol(index_panel, 'Style', 'pushbutton', 'String', '+', 'FontSize', 20, ...
                'Units', 'Normalized', 'Position', [0.8 top_margins 0.18 1-2*top_margins], ...
                'Callback', @obj_vec.callback_next, 'UserData', input.parent, ...
                'Tooltip', 'go forward to next candidate'); 
            
            %%%%%%%%%%%%%%%%%%%%%% panel classify %%%%%%%%%%%%%%%%%%%%%%%%%
            
            pos = pos + 0.05 + index_panel.Position(3); 
            
            classify_panel = uipanel(input.parent, 'Units', 'Normalized', 'Position', [pos 0.0 0.2 0.1], 'title', 'classify'); 
            
            button = uicontrol(classify_panel, 'Style', 'pushbutton', 'String', 'Classify', 'FontSize', 18, ...
                'Units', 'Normalized', 'Position', [0.025 top_margins 0.45 1-2*top_margins], ...
                'Callback', @obj_vec.callback_classify, 'UserData', input.parent); 
            
            button.TooltipString = 'choose a classification for this candidate'; 
            
            button = uicontrol(classify_panel, 'Style', 'pushbutton', 'String', 'SAVE', 'FontSize', 20, ...
                'Units', 'Normalized', 'Position', [0.525 top_margins 0.45 1-2*top_margins], ...
                'Callback', @obj_vec.callback_save, 'UserData', input.parent); 
            
            kept = obj_vec(logical([obj_vec.kept])); 
            
            N_class = nnz(~cellfun(@isempty, {kept.classification}));
            N_total = numel(kept); 
            
            if isprop(button, 'Tooltip')
                tool_tip_name = 'Tooltip';
            elseif isprop(button, 'TooltipString')
                tool_tip_name = 'TooltipString';
            end
            
            if N_class==N_total
                button.String = 'SAVE'; 
                button.(tool_tip_name)= sprintf('Save %d candidates', N_total); 
            else                
                button.String = sprintf('%d/%d', N_class, N_total); 
                button.(tool_tip_name) = sprintf('Cannot save. Only %d candidates have been classified out of %d', N_class, N_total); 
            end
            
            % TODO: notes section + add note button
            % TODO: popup stack image (+load from file) 
            % TODO: popup cutouts viewer
            
        end
        
        function showRawFlux(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('timestamps', false); % show frame indices or timestamps
            input.input_var('title', true); 
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('font_size', 16);
            input.input_var('equal_limits', false); 
            input.input_var('on_top', false); 
            input.scan_vars(varargin{:});
            
            if isempty(input.ax), input.ax = gca; end
            
            if input.timestamps
                x = obj.timestamps; 
                xk = obj.kernel_timestamps+obj.peak_timestamp;
            else
                x = 1:length(obj.timestamps); 
                xk = 1:length(obj.kernel);
                xk = xk + obj.time_index - length(obj.kernel)/2;    
            end
            
            yyaxis(input.ax, 'left'); 
            
            input.ax.NextPlot = 'replace';
            
            f = obj.flux_raw;
            a = obj.auxiliary(:,obj.aux_indices.areas); 
            b = obj.auxiliary(:,obj.aux_indices.backgrounds); 
            f2 = f - nanmedian(a.*b);
            
            h1 = plot(input.ax, x, f2, '-', 'LineWidth', 2, 'Color', [0.929 0.694 0.125]);
            h1.DisplayName = 'raw flux';
            
            input.ax.NextPlot = 'add';
            h2 = plot(input.ax, xk, obj.flux_mean*obj.kernel_lightcurve(~obj.is_positive), ':', 'LineWidth', 2, 'Color', input.ax.Colormap(4,:));
            h2.DisplayName = 'best kernel';
            
            input.ax.NextPlot = 'replace';
            
            input.ax.FontSize = input.font_size;
            input.ax.YAxis(1).Color = [0 0 0];
            
            if input.ax.YLim(2)>0
                if nanmin(f2)<0
                    input.ax.YLim(1) = nanmin(f2).*1.1;
                else
                    input.ax.YLim(1) = 0;
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            yyaxis(input.ax, 'right'); 
                        
            input.ax.NextPlot = 'replace';
            input.ax.ColorOrderIndex = 5;
            
            h3 = plot(input.ax, x, obj.auxiliary(:,obj.aux_indices.backgrounds), '--', 'Color', [0.2 0.6 0.2]); 
            h3.DisplayName = 'background';
            
            aux = obj.auxiliary(:,obj.aux_indices.backgrounds);
            
            input.ax.NextPlot = 'add';
            
            w = obj.auxiliary(:,obj.aux_indices.widths);
            
            w = w.*2.355; 
            
%             if ~isempty(obj.head) && ~isempty(obj.head.SCALE)
%                 w = w.*obj.head.SCALE;
%             end
            
            h4 = plot(input.ax, x, w, 'p', 'MarkerSize', 4, 'Color', [0.2 0.5 1]);
            h4.DisplayName = 'PSF FWHM';
            
            aux = horzcat(aux, w); 
            
            h5 = plot(input.ax, x, obj.relative_dx, 'x', 'MarkerSize', 7, 'Color', [1 0.3 0.2]);
            h5.DisplayName = 'relative dx';
            
            aux = horzcat(aux, obj.relative_dx); 
            
            h6 = plot(input.ax, x, obj.relative_dy, '+', 'MarkerSize', 6, 'Color', [0.8 0.3 0.4]);
            h6.DisplayName = 'relative dy';
            
            aux = horzcat(aux, obj.relative_dy); 
            
            if obj.is_positive
                leg_pos_y = 'South'; 
            else
                leg_pos_y = 'North'; 
            end
            
            if obj.time_index>=floor(length(obj.timestamps)/2)
                leg_pos_x = 'West'; 
            else
                leg_pos_x = 'East'; 
            end
            
            lh = legend(input.ax, 'Location', [leg_pos_y leg_pos_x], 'Orientation', 'Vertical');
            lh.FontSize = input.font_size-4;
%             lh.NumColumns = 3;

            if input.timestamps
                
                xlabel(input.ax, 'timestamp (seconds)');
                
                if obj.timestamps(1)<obj.timestamps(end)
                    input.ax.XLim = [obj.timestamps(1) obj.timestamps(end)];
                end
                
            else
                xlabel(input.ax, 'frame index');
                input.ax.XLim = [x(1) x(end)]; 
            end
            
%             title(input.ax, strjoin(obj.notes, ', '), 'FontSize', input.font_size);
            
            if input.title
                util.plot.inner_title(obj.printSummary, 'ax', input.ax, 'Position', 'NorthWest', 'FontSize', input.font_size);
            end
            
            mx = util.stat.max2(aux(obj.time_range,:));
            mn = util.stat.min2(aux(obj.time_range,:)); 
            
            if ~isnan(mn) && ~isnan(mx)
                input.ax.YLim = [mn-0.25.*abs(mn) mx+0.25.*abs(mx)]; 
            end
            
%             ylabel(input.ax, 'pixels or count/pixel');
            ylabel(input.ax, 'auxiliary data');
    
            input.ax.YAxis(2).Color = [0 0 0];
            
            input.ax.NextPlot = 'replace';
            
            input.ax.FontSize = input.font_size;
            
            yyaxis(input.ax, 'left'); % go back to the left axis by default
            
            ylabel(input.ax, 'raw flux [counts]');
            
            if input.equal_limits
                mx = max(abs(obj.flux_raw - obj.flux_mean)); 
                input.ax.YLim = obj.flux_mean + [-1 1].*1.2.*mx;
            end
            
            if input.on_top
                input.ax.XTick = [];
%                 input.ax.YTick = input.ax.YTick(2:end); 
            end
            
        end
        
        function showFilteredFlux(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('timestamps', false); % show frame indices or timestamps
            input.input_var('title', true); 
            input.input_var('cuts', false);
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('font_size', 16);
            input.input_var('equal_limits', false); 
            input.input_var('on_top', false); 
            input.scan_vars(varargin{:});
            
            if isempty(input.ax), input.ax = gca; end
            
            if input.timestamps
                x = obj.timestamps; 
                xk = obj.kernel_timestamps+obj.peak_timestamp;
            else
                x = 1:length(obj.timestamps); 
                xk = 1:length(obj.kernel);
                xk = xk + obj.time_index - length(obj.kernel)/2;                
            end
            
            input.ax.NextPlot = 'replace';
            h1 = plot(input.ax, x, obj.flux_filtered, 'LineWidth', 2);
            h1.DisplayName = 'filtered flux';
            
            input.ax.NextPlot = 'add';
            range = obj.time_range;
            if ~isempty(range)
                h2 = plot(input.ax, x(range), obj.flux_filtered(range), 'LineWidth', 2);
                h2.DisplayName = 'time range';
            end
            
            f = obj.flux_corrected;
            s = nanstd(obj.flux_corrected); 
            
            input.ax.ColorOrderIndex = 3;
            h3 = plot(input.ax, x, f./s, '-');
            h3.DisplayName = 'corrected flux';
            
            sign = 1;
            if obj.is_positive==0
                sign = -1;
            end

            h4 = plot(input.ax, xk, obj.kernel*5*sign, '--');
            h4.DisplayName = 'best kernel';
            
            if obj.is_positive
                leg_pos_y = 'North'; 
            else
                leg_pos_y = 'South'; 
            end
            
            if obj.time_index>=floor(length(obj.timestamps)/2)
                leg_pos_x = 'West'; 
            else
                leg_pos_x = 'East'; 
            end
            
            lh = legend(input.ax, 'Location', [leg_pos_y leg_pos_x], 'Orientation', 'Vertical');
            lh.FontSize = input.font_size-4;
%             lh.NumColumns = 3;

            if input.equal_limits
                mx = max(max(abs(obj.flux_filtered)), max(abs(f./s)));             
                input.ax.YLim = [-1 1].*1.2.*mx;
            else
                ff_ext = [ -abs(nanmin(obj.flux_filtered)).*1.3 abs(nanmax(obj.flux_filtered).*1.3)]; % fluxes_filtered extrema
                fc_ext = [ -abs(nanmin(f./s)).*1.3 abs(nanmax(f./s).*1.3)]; % fluxes_corrected extrema
                input.ax.YLim = [min(fc_ext(1), ff_ext(1)) max(fc_ext(2), ff_ext(2))];
            end
            
            if input.timestamps
                
                xlabel(input.ax, 'timestamp (seconds)');
                
                if obj.timestamps(1)<obj.timestamps(end)
                    input.ax.XLim = [obj.timestamps(1) obj.timestamps(end)];
                end
                
            else
                xlabel(input.ax, 'frame index');
                input.ax.XLim = [x(1) x(end)]; 
            end
            
            ylabel(input.ax, 'flux [S/N]');
            
            if input.on_top
                input.ax.XTick = [];
%                 input.ax.YTick = input.ax.YTick(2:end); 
            end
            
%             title(input.ax, strjoin(obj.notes, ', '), 'FontSize', input.font_size);
            
            if input.title
                util.plot.inner_title(obj.printSummary, 'ax', input.ax, 'Position', 'NorthWest', 'FontSize', input.font_size);
            end
            
            input.ax.NextPlot = 'replace';
            input.ax.FontSize = input.font_size;
            
            hits = obj.cut_hits;
            hits = unique([hits, obj.display_cuts]); 
                
            if input.cuts && ~isempty(hits)
                                
                yyaxis(input.ax, 'right');
                input.ax.NextPlot = 'replace'; 
                
                vectors = [];
                
                colors = {'g', 'b', 'y'}; 
                
                for ii = 1:length(hits)
                    
                    vec = obj.cut_matrix(:,hits(ii));
                    
                    
                    if obj.cut_two_sided(hits(ii))
                        region = single(abs(vec)>=obj.cut_thresh(hits(ii)));
                    else
                        region = single(vec>=obj.cut_thresh(hits(ii)));
                        vec(vec<0) = 0; 
                    end
                    
                    col_idx = mod(ii-1,length(colors))+1; 
                    
                    h = plot(input.ax, x, vec, '-', 'LineWidth', 1.5, 'Color', colors{col_idx}); 
                    
                    cut_idx = find(obj.cut_hits==hits(ii));
                    
                    if ~isempty(cut_idx) % if this was a triggered cut 
                        val = obj.cut_value(cut_idx);
                    else % or a user added cut
                        val = obj.cut_matrix(obj.time_index, hits(ii));
                    end
                    
                    h.DisplayName = strrep(sprintf('%s=%4.2f', obj.cut_names{hits(ii)}, val), '_', ' '); 
                    
                    vectors = horzcat(vectors, vec); % keep track of all cuts to set the bounds
                    
                    input.ax.NextPlot = 'add'; 
                    
                    region(~region) = NaN; 
                    
                    plot(input.ax, x, region.*vec, '*', 'Color', 'm', 'HandleVisibility', 'off'); 
                    plot(input.ax, x, ones(length(x),1).*obj.cut_thresh(hits(ii)), ':', 'Color', 'm', 'HandleVisibility', 'off'); 
                    
                end
                
                input.ax.NextPlot = 'replace'; 
                
                mn = util.stat.min2(vectors); 
                mx = util.stat.max2(vectors);
                
                if mn~=0 || mx~=0
                    input.ax.YLim = [mn - 0.2*abs(mn), mx + 0.2*abs(mx)]; 
                end
                
                input.ax.YAxis(2).Color = [0 0 0];
                ylabel(input.ax, 'cut values'); 
            
                yyaxis(input.ax, 'left');
                
            end
            
        end
        
        function showCutouts(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('crosshair', true); 
            input.input_var('parent', []); 
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('position', []);
            input.input_var('number', 9);
            input.input_var('bias', []);
            input.input_var('dynamic_range', []);
            input.input_var('font_size', 10);
            input.scan_vars(varargin{:});
            
            if isempty(obj.cutouts)
                return;
            end
            
            % deal with the rest of this later! 
            
            if isempty(input.parent) && isempty(input.ax)
                    
                input.parent = gcf;
%                     delete(input.parent.Children);
                if ~isempty(input.position)
                    panel = util.plot.stretchy_panel('Parent', input.parent, 'Position', input.position);
                else
                    panel = util.plot.stretchy_panel('Parent', input.parent);
                end

            elseif isempty(input.parent) && ~isempty(input.ax)

                pos = input.ax.Position;

                parent = input.ax.Parent;

                panel = util.plot.stretchy_panel('Position', pos, 'Parent', parent);

                delete(input.ax);

            elseif ~isempty(input.parent)
                panel = input.parent;
            end

            idx = obj.time_index; 

            idx_start = idx - floor(input.number/2); 
            idx_end = idx + ceil(input.number/2) - 1;

            if idx_start<1
                idx_end = idx_end + 1 - idx_start;
                idx_start = 1;
            end

            if idx_end>size(obj.cutouts,3)
                idx_start = idx_start + size(obj.cutouts,3) - idx_end;
                idx_end = size(obj.cutouts,3);
            end

            rad = [];
%             if ~isempty(obj.relative_dx) && ~isempty(obj.relative_dx)
% 
%                 cen = floor([size(obj.cutouts,2), size(obj.cutouts,1)]/2)+1;
%                 cen = cen + [obj.relative_dx(idx_start:idx_end) obj.relative_dx(idx_start:idx_end)];
% 
%                 if ~isempty(obj.head) && ~isempty(obj.head.PHOT_PARS) && ~isempty(obj.head.PHOT_PARS.aperture_radius)
%                     rad = obj.head.PHOT_PARS.aperture_radius;
%                     str = sprintf('ap= %4.2f', rad(end)); 
%                     col = 'green';
%                 end
% 
%             end

            Nrows = ceil(sqrt(input.number));
            Ncols = Nrows;

            if isempty(input.bias) || isempty(input.dynamic_range)
                dyn = util.img.autodyn(obj.cutouts(:,:,obj.time_range));
%                 dyn = [0, nanmax(squeeze(util.stat.max2(obj.cutouts(:,:,obj.time_range))))];
                if ~isempty(input.bias)
                    dyn(1) = input.bias;
                elseif ~isempty(input.dynamic_range)
                    dyn(2) = input.dynamic_range;
                end

            else
                dyn = [input.bias, input.dynamic_range];
            end

            for ii = 1:input.number

                x = mod(ii-1, Nrows);
                y = floor((ii-1)/Nrows);

                ax{ii} = axes('Position', [x/Ncols y/Nrows 1/Ncols 1/Nrows], 'Parent', panel);

                imagesc(ax{ii}, obj.cutouts(:,:,idx_start+ii-1));

                if dyn(1)>=dyn(2)
                    dyn(2) = dyn(1)+1;
                end
                
                ax{ii}.CLim = dyn;
                axis(ax{ii}, 'image');
                ax{ii}.XTick = [];
                ax{ii}.YTick = [];

                if input.crosshair

                    ax{ii}.NextPlot = 'add'; 

                    x0 = ceil(size(obj.cutouts,1)/2); 

                    plot(ax{ii}, [x0 x0], [0.5 x0-1], '-m'); 
                    plot(ax{ii}, [0.5 x0-1], [x0 x0], '-m'); 

                    ax{ii}.NextPlot = 'replace'; 

                end

                if idx_start+ii-1==idx
                    util.plot.inner_title([num2str(idx_start+ii-1) '*'], 'Position', 'NorthWest', 'Color', 'red', 'FontSize', input.font_size, 'ax', ax{ii});
                else
                    util.plot.inner_title(num2str(idx_start+ii-1), 'Position', 'NorthWest', 'FontSize', input.font_size, 'ax', ax{ii});
                end

                if ~isempty(rad)
                    viscircles(ax{ii}, cen(ii,:), rad(end), 'EdgeColor', col);
                    if idx_start+ii-1==idx
                        util.plot.inner_title(str, 'Position', 'bottom', 'Color', col, 'FontSize', input.font_size, 'ax', ax{ii});
                    end
                end

            end

            panel.Children = panel.Children([5,1,2,3,4,6,7,8,9]); % move the center cutout to top (for the inner titles)
%                 
%                 clim = ax{idx-idx_start+1}.CLim;
%                 for ii = 1:length(ax)
%                     ax{ii}.CLim = clim;
%                 end
%                     
            
        end
        
        function popupSecrets(obj)
            
            str1 = obj.printRunData;
            str2 = obj.printSimulationData;
            
            options.WindowStyle = 'non-modal'; 
            options.Interpreter = 'tex'; 
            
            if isempty(str2)
                str = [str1 newline 'This event is not simulated!']; 
            else
                str = (strjoin({str1, str2}, newline));
            end
            
            str = strrep(str, '_', ' '); 
            str = ['\fontsize{16}' str];
            
            msgbox(str, 'secret info', options); 
            
        end
        
        function popupClassifier(obj)
            
%             str = 'Classify this candidate:';
%             
%             rep = listdlg('ListString', classes, 'PromptString', str, 'ListSize',[150,250],...
%                 'InitialValue',length(classes), 'SelectionMode', 'single', 'Name','classification');
%             
%             if ~isempty(rep)
%                 if strcmp(classes{rep}, 'unclassify')
%                     obj.classification = '';
%                 else
%                     obj.classification = classes{rep}; 
%                 end
%             end
            
            classes = obj.getListOfClasses;
            N = ceil(length(classes)/2)+1; % number of rows
            
            f = util.plot.FigHandler('Classifications'); 
            f.clear;
            f.fig.WindowStyle = 'modal';
            
            for ii = 1:N-1
                
                uicontrol(f.fig, 'Style', 'pushbutton', 'String', classes{ii*2-1}, ...
                    'Units', 'Normalized', 'Position', [0 (N-ii)/N 0.5 1/N], ...
                    'Callback', @obj.callback_classification_button, 'FontSize', 18); 
                
                uicontrol(f.fig, 'Style', 'pushbutton', 'String', classes{ii*2}, ...
                    'Units', 'Normalized', 'Position', [0.5 (N-ii)/N 0.5 1/N], ...
                    'Callback', @obj.callback_classification_button, 'FontSize', 18); 
                
            end
            
            uicontrol(f.fig, 'Style', 'pushbutton', 'String', 'Unclassify', ...
                    'Units', 'Normalized', 'Position', [0 0 1 1/N], ...
                    'Callback', @obj.callback_unclassify_button, 'FontSize', 18); 
            
            uiwait(f.fig); 
                
        end
        
        function popupPSD(obj)
            
            f = figure('Name', 'Power Spectral Density'); 
            
            ax = axes('Parent', f); 
            
            obj.showPSD('axes', ax); 
                        
        end
        
        function showPSD(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('sqrt', false); % show the sqrt of the PSD
            input.input_var('half', true); % show only half the frequencies
            input.input_var('log_x', true); % show the x-axis on a log-scale (y-axis is always in log!)
            input.input_var('axes', [], 'axis'); % which axes to plot to? default is gca()
            input.input_var('font_size', 20); % fonts on the axes
            input.scan_vars(varargin{:}); 
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            N = length(obj.freq_psd); 
            
            if input.half
                N = ceil(N/2); 
            end
            
            x = obj.freq_psd(1:N); 
            y = obj.psd(1:N); 
            
            if input.sqrt
                y = sqrt(y); 
            end
            
            semilogy(input.axes, x, y);
            
            xlabel(input.axes, 'Frequency [Hz]'); 
            ylabel(input.axes, 'Power Spectral Density'); 
            
            if input.sqrt
                ylabel(input.axes, 'sqrt(PSD)'); 
            end
            
            if input.log_x
                input.axes.XScale = 'log'; 
            end
            
            input.axes.XLim = [x(1), x(end)]; 
            
            input.axes.FontSize = input.font_size; 
            
        end
        
        function popupOtherStars(obj, ~, ~)
            
            [idx, corr] = obj.findHighestCorrelations; % indices of stars that have highest correlations to this star
            
            flux = obj.flux_raw_all(:,idx); 
            
            f = util.plot.FigHandler('Other stars, flux and position'); 
            f.clear;
            ax1 = axes('Parent', f.fig, 'Position', [0.1 0.2 0.88 0.75]); 
            
            mx = nanmax([nanmax(obj.flux_raw), nanmax(flux)]).*2;
            
            area(ax1, obj.time_range, ones(length(obj.time_range),1).*mx, 'FaceColor', 'g', ...
                'EdgeColor', 'none', 'FaceAlpha', 0.25, 'DisplayName', 'Event region'); 
            
            hold(ax1, 'on');
            
            h1 = plot(ax1, obj.flux_raw, 'LineWidth', 3);
            h1.DisplayName = sprintf('Candidate star, idx= %d', obj.star_index); 
            
            h2 = plot(ax1, flux, 'LineWidth', 1.5); 
            
            for ii = 1:length(h2)
                h2(ii).DisplayName = sprintf('idx= %d | corr= %4.2f | trig= %d', idx(ii), corr(ii), ismember(idx(ii), obj.star_extra)); 
            end
            
            xlabel(ax1, 'Frame index'); 
            ylabel(ax1, 'Raw flux'); 
            
            ax1.FontSize = 14;
            ax1.YScale = 'log'; 
            
            hold(ax1, 'off'); 
            
            legend(ax1, 'Location', 'NorthEastOutside'); 
            
            ax2 = axes('Parent', f.fig, 'Position', [0.73 0.1 0.25 0.5]); 
            
            x = nanmean(obj.auxiliary_all(:,:,obj.aux_indices.centroids_x)); 
            y = nanmean(obj.auxiliary_all(:,:,obj.aux_indices.centroids_y)); 
            
            plot(ax2, x(obj.star_index), y(obj.star_index), 'x', 'Color', h1.Color, 'MarkerSize', 15); 
            
            hold(ax2, 'on'); 
            
            for ii = 1:length(idx)
                
                plot(ax2, x(idx(ii)), y(idx(ii)), 'o', 'Color', h2(ii).Color, 'MarkerSize', 10); 
                
            end
            
            hold(ax2, 'off'); 
            
            ax2.XLim = [1 obj.head.NAXIS1]; 
            ax2.YLim = [1 obj.head.NAXIS2]; 
            
            ax2.FontSize = 12;
            
            uicontrol(f.fig, 'Style', 'pushbutton', 'String', 'CLOSE WINDOW', ...
                'Units', 'Normalized', 'Position', [0.5 0.02 0.2 0.08], ...
                'Callback', @obj.callback_close_window, 'FontSize', 12); 
            
        end
        
        function popupStack(obj, ~, ~)
            
            if isempty(obj.stack)
                filename = obj.filenames{obj.time_index}; 
                obj.stack = h5read(filename, '/stack'); 
            end
            
            f = util.plot.FigHandler('Stack image'); 
            f.width = 18;
            f.height = 16;
            f.clear;
            ax = axes('Parent', f.fig); 
            
            util.plot.show(obj.stack, 'auto', 1); 
            title(ax, ''); 
            
            hold(ax, 'on'); 
            
            x = nanmedian(obj.auxiliary(:,obj.aux_indices.centroids_x)); 
            y = nanmedian(obj.auxiliary(:,obj.aux_indices.centroids_y)); 
            
            plot(ax, x, y, 'go', 'MarkerSize', 15); 
            
            hold(ax, 'off'); 
            
        end
        
        function popupCutouts(obj, ~, ~)
            
            f = util.plot.FigHandler('Cutout viewer'); 
            f.width = 25;
            f.height = 24;
            f.clear;
            
            N = 2*length(obj.time_range); 
            
            ax = util.plot.show_cutouts(obj.cutouts, 'frame', obj.time_index, 'number', N); 
            
            mn = 0;
            mx = 1;
            
            for ii = 1:length(ax)
                if ax{ii}.CLim(1)<mn 
                    mn = ax{ii}.CLim(1);
                end
                if ax{ii}.CLim(2)>mx 
                    mx = ax{ii}.CLim(2); 
                end
            end
            
            for ii = 1:length(ax)
                ax{ii}.CLim = [mn,mx]; 
            end
            
        end
        
        function plotApertureDifferences(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis'); 
            input.input_var('font_size', 16); 
            input.scan_vars(varargin{:}); 
            
            if isempty(obj.flux_unforced)
                obj.loadUnforcedPhotometry; % this may fail if not given the root folder explicitely (call this function with arguments, then plot!)
            end
                
            if isempty(input.ax)
                input.ax = gca;
            end
            
            t = obj.timestamps;
            f1 = obj.flux_raw; 
            f2 = obj.flux_unforced;
            
            plot(input.ax, t, f1, t, f2, 'LineWidth', 1.5); 
            
            xlabel(input.ax, 'Time [s]'); 
            ylabel(input.ax, 'Raw flux [counts]'); 
            
            input.ax.FontSize = input.font_size; 
            
            legend(input.ax, {'forced', 'aperture'}, 'Location', 'SouthEast'); 
            
        end
        
    end
    
    methods % callbacks to the show() method
        
        function callback_prev(obj, hndl, ~)
            
            num = hndl.UserData.UserData.number;
            idx = hndl.UserData.UserData.index;
            use_kept = hndl.UserData.UserData.show_kept;
            
            if use_kept
                N = num+1;
            else
                N = 1;
            end
            
            for ii = 1:N
                
                idx = idx - 1; 
                if idx<1
                    idx = num;
                end
                
                if use_kept && obj(idx).kept % when we stumble on a kept event
                    break;
                end
                
            end
            
            obj.show('parent', hndl.UserData, 'index', idx); 
            
        end
        
        function callback_index(obj, hndl, ~)
            
            num = hndl.UserData.UserData.number;
            idx = hndl.UserData.UserData.index;
            
            val = util.text.extract_numbers(hndl.String);
            
            if ~isempty(val)
                val = val{1}; 
            end
            
            if ~isempty(val)
                idx = val(1);
            end
            
            if idx<1 || idx>num
                error('Index is outside range of candidates (1 to %d)', num); 
            end
            
            obj.show('parent', hndl.UserData, 'index', idx);
            
        end
        
        function idx = getNextIndex(obj, hndl, num_steps)
            
            if nargin<3 || isempty(num_steps)
                num_steps = 1;
            end
            
            num = hndl.UserData.UserData.number;
            idx = hndl.UserData.UserData.index;
            use_kept = hndl.UserData.UserData.show_kept;
            
            if use_kept
                N = num+1;
            else
                N = 1;
            end
            
            for s = 1:num_steps

                for ii = 1:N

                    idx = idx + 1; 
                    if idx>num
                        idx = 1;
                    end

                    if use_kept && obj(idx).kept % when we stumble on a kept event
                        break;
                    end

                end

            end
            
        end
        
        function callback_next(obj, hndl, ~)
            
            idx = obj.getNextIndex(hndl); 
            obj.show('parent', hndl.UserData, 'index', idx); 
            
        end
        
        function callback_kept(obj, hndl, ~)
            
            use_kept = ~hndl.UserData.UserData.show_kept;
            
            hndl.UserData.UserData.show_kept = use_kept;
            
            hndl.Value = use_kept;
            
            if use_kept
%                 hndl.String = sprintf('kept (%d)', hndl.UserData.UserData.num_kept);
                hndl.String = 'kept'; 
            else
                hndl.String = 'all';
            end
            
        end
        
        function callback_cuts(obj, hndl, ~)
            
            use_cuts = ~hndl.UserData.UserData.show_cuts;
            
            hndl.UserData.UserData.show_cuts = use_cuts;
            
            hndl.Value = use_cuts;
            
            if use_cuts
                hndl.String = 'show cuts';
            else
                hndl.String = 'no cuts';
            end
            
            obj.show('parent', hndl.UserData);
            
        end
        
        function callback_show_secrets(obj, hndl, ~)
            
            obj.popupSecrets; 
            
        end
        
        function callback_parse_cuts(obj, hndl, ~)
            
            idx = hndl.UserData.UserData.index;
            str = hndl.String;
            
            c = util.text.parse_inputs(str);
            
            obj(idx).display_cuts = [];
            
            for ii = 1:length(c)
                
                arg = util.text.parse_value(c{ii}); 
                
                if iscell(arg)
                    arg = strjoin(arg, '_');
                end
                
                if ischar(arg) 
                    if isfield(obj(idx).cut_indices, arg)
                        value = obj(idx).cut_indices.(arg); 
                    else
                        value = [];
                    end
                elseif isnumeric(arg) && isscalar(arg)
                    value = arg;
                else
                    value = [];
                end
                
                if ~isempty(value) && ~ismember(value, obj(idx).cut_hits) && ~ismember(value, obj(idx).display_cuts)
                    obj(idx).display_cuts = horzcat(obj(idx).display_cuts, value); 
                end
                
            end
            
            obj.show('parent', hndl.UserData); 
            
        end
        
        function callback_classify(obj, hndl, ~)
            
            idx = hndl.UserData.UserData.index;
            
            obj(idx).popupClassifier;
            
            if isempty(obj(idx).classification) % no classification, leave the same index
                obj.show('parent', hndl.UserData);
            else
                
                idx_new = obj.getNextIndex(hndl); 
                
                if hndl.UserData.UserData.skip_duplicates && obj(idx).isSameEvent(obj(idx_new))
                    obj(idx_new).classification = obj(idx).classification; % copy the same classification as this candidate
                    idx_new = obj.getNextIndex(hndl, 2); % skip to the next one after that
                end
                
                obj.show('parent', hndl.UserData, 'index', idx_new);
                
            end
            
        end
        
        function callback_classification_button(obj, hndl, ~)
            
            obj.classification = hndl.String; 
            
            delete(hndl.Parent); 
            
        end
        
        function callback_unclassify_button(obj, hndl, ~)
            
            obj.classification = ''; 
            
            delete(hndl.Parent); 
            
        end
        
        function callback_save(obj_vec, ~, ~)
            
            dir = obj_vec(1).folder; 

            if isempty(dir)
                dir = obj_vec(1).findAnalysisFolder;
            end

            if ~exist(dir, 'dir')
                dir = obj_vec(1).findAnalysisFolder; % try to re-calculate the folder position
            end

            filename = fullfile(dir, 'classified.mat'); 

            f = util.plot.FigHandler('saving candidates'); 
            f.clear;
            f.width = 20;
            f.height = 13;
            
            f.fig.WindowStyle = 'modal';
            
            kept = obj_vec(logical([obj_vec.kept])); 
            
            classifications = {kept.classification};
            
            % list of classifications, with occultations first
            classes = unique(classifications);
            classes(find(strcmp(classes, 'occultation certain'))) = [];
            classes(find(strcmp(classes, 'occultation possible'))) = [];
            classes = ['occultation certain', 'occultation possible', classes]; 
            
            N_class = nnz(~cellfun(@isempty, classifications));
            N_kept = numel(kept); 
            
            if N_class==N_kept
                
                buttons = struct; 
                
                candidates = obj_vec(~cellfun(@isempty, {obj_vec.classification}));
                
                %%%%%%%%%%% total number of candidates %%%%%%%%%%%%%
                
                str = sprintf('A total of %d candidates were classified', length(candidates));
                
                if length(candidates)~=N_kept
                    str = sprintf('%s\nout of %d kept candidates', str, N_kept);
                end
                
                buttons.total = uicontrol(f.fig, 'Style', 'text', 'String', str, ...
                    'Units', 'Normalized', 'Position', [0 0.8 1 0.2], ...
                    'HorizontalAlignment', 'Center', 'FontSize', 18, 'UserData', N_kept);
                
                %%%%%%%%%%% classification break-down %%%%%%%%%%%%%%%
                
                str = ['  *total classifications*  ' newline];
                
                for ii = 1:length(classes)
                    str = sprintf('%s "%s": %d\n', str, classes{ii}, nnz(strcmp(classes{ii}, classifications))); 
                end
                
                buttons.classifications = uicontrol(f.fig, 'Style', 'text', 'String', str, ...
                    'Units', 'Normalized', 'Position', [0 0.3 0.5 0.5], ...
                    'HorizontalAlignment', 'Left', 'FontSize', 18);
                
                %%%%%%%%%%% Simulation results %%%%%%%%%%%%%%%%%%%
                
                str = ''; 
                alt_str = '';
                
                simulated = candidates(logical([candidates.is_simulated])); 
                N_sim = length(simulated); 
                
                if N_sim
                    
                    alt_str = ['   *simulated events*  ' newline]; 
                    
%                     str = sprintf('%s Total simulated: %d\n', str, N_sim); 
                    
                    N_certain = nnz([candidates.is_simulated] & strcmp('occultation certain', {candidates.classification}));
                    alt_str = sprintf('%s "occultation certain": %d\n', alt_str, N_certain); 
                    
                    N_possible = nnz([candidates.is_simulated] & strcmp('occultation possible', {candidates.classification}));
                    alt_str = sprintf('%s "occultation possible": %d\n', alt_str, N_possible);
                    
                    alt_str = sprintf('%s all other classes: %d\n', alt_str, N_sim - N_certain - N_possible); 
                    
                    points_certain = 100/N_sim; % only get full points if correctly classify everything! 
                    points_possible = points_certain/2; % only get half points for possible occultation
                    points_wrong = -points_possible; % any simulated event missed by classifier gives negative points! 
                    
%                     str = sprintf('%s points: %4.2f and %4.2f\n', str, points_certain, points_possible); 
                    
                    total_points = points_certain.*N_certain + points_possible.*N_possible - points_wrong*(N_sim - N_certain - N_possible);

                    alt_str = sprintf('%s\n Classifier grade: %d%%\n', alt_str, round(total_points)); 
                    
                end
                
                buttons.simulations = uicontrol(f.fig, 'Style', 'text', 'String', str, ...
                    'Units', 'Normalized', 'Position', [0.5 0.3 0.5 0.5], ...
                    'HorizontalAlignment', 'Left', 'FontSize', 18, 'UserData', alt_str);
                
                %%%%%%%%%%%% filename confirmation %%%%%%%%%%%%%%%%
                
%                 str = sprintf('Ready to save %d candidates to\n"%s"', N_total, filename); 
                str = ''; % I don't want to reveal the run name until candidates are saved! 
                
                buttons.filename = uicontrol(f.fig, 'Style', 'text', 'String', str, ...
                    'Units', 'Normalized', 'Position', [0.0 0.1 1 0.2], ...
                    'HorizontalAlignment', 'Left', 'FontSize', 12, 'UserData', filename);
                
                %%%%%%%%%%%%%%%%% Save / Cancel %%%%%%%%%%%%%%%%%%%%%%
                
                buttons.save = uicontrol(f.fig, 'Style', 'pushbutton', 'String', 'Save candidates', ...
                    'Units', 'Normalized', 'Position', [0.0 0.0 0.5 0.1], ...
                    'HorizontalAlignment', 'Center',  'FontSize', 18, ...
                    'Callback', @candidates.callback_confirm_save);
                
                buttons.cancel = uicontrol(f.fig, 'Style', 'pushbutton', 'String', 'Cancel', ...
                    'Units', 'Normalized', 'Position', [0.5 0.0 0.5 0.1], ...
                    'HorizontalAlignment', 'Center',  'FontSize', 18, ...
                    'Callback', @obj_vec.callback_close_window);
                
                uicontrol(buttons.cancel); 
                
                buttons.save.UserData = buttons; 
                
            else                
                
                button = uicontrol(f.fig, 'Style', 'text', 'String', sprintf('\n\nOnly %d out of %d candidates were classified. \n\n Classify all candidates before saving!', N_class, N_kept), ...
                    'Units', 'Normalized', 'Position', [0 0.1 1 0.9], ...
                    'HorizontalAlignment', 'Center', 'FontSize', 18);
                    
                button = uicontrol(f.fig, 'Style', 'pushbutton', 'String', 'Cancel', ...
                    'Units', 'Normalized', 'Position', [0 0.0 1 0.1], ...
                    'HorizontalAlignment', 'Center',  'FontSize', 18, ...
                    'Callback', @obj_vec.callback_close_window);
                
                uicontrol(button); 
                
            end
            
            uiwait(f.fig); 
            
        end
        
        function callback_confirm_save(obj_vec, hndl, ~)
            
            buttons = hndl.UserData; % recover handles to all buttons
            
            filename = buttons.filename.UserData; 
            
            fprintf('Saving %d candidates to "%s"\n', length(obj_vec), filename);
            
            buttons.filename.String = sprintf(' Successfully saved %d candidates to \n "%s"', buttons.total.UserData, filename); 
            
            obj_vec.saveClassified(filename); 
            
            buttons.simulations.String = buttons.simulations.UserData; % load the alt-text
            
            buttons.save.Callback = @obj_vec.callback_close_window; 
            buttons.save.String = 'Close figure'; 
            buttons.save.Position(3) = 1; % stretch it over the "cancel" button
            
            delete(buttons.cancel); 
            
        end
        
        function callback_close_window(~, hndl, ~)
            
            delete(hndl.Parent); 
            
        end
        
    end
    
end

