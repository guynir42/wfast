classdef DataStore < handle
% Storage for the flux and auxiliary data from the photometry object, 
% keeping track of buffers of past data, organizing the auxiliary in a more
% reasonable way and providing the other objects in the +tno package easy
% access to the data they need on the timescales they need. 
%
% To input data you can use varargin pairs, or supply your own InputVar object
% or simply give it an img.Photometry object like this:
% >> phot = img.Photometry;
% >> phot.input(calibrated_cutouts, ...); 
% >> store.input(phot); 
% The last option is the prefered method to input data to the store. 
%
% The user-defined parameters of this class are saved as a struct "pars". 
% These parameters are defined in the "reset/clear" methods block, inside
% the resetPars() function. Refer to the inline docs for info on each parameter. 
%
% Some notable parameters are:
% -length_XXXX: see the definitions below. 
% -threshold and use_threshold: this defines a cut on the stars to be kept 
%  in the store after the burn-in period (see below). Default: 3 and true. 
% -use_remove_cosmic_rays and cosmic_rays_threshold are used to get rid of
%  cosmic ray spikes in the flux before even giving the data to any of the 
%  other objects that use it (e.g., the PSD or the event calculcator). 
%  In short it looks for strong, isolated peaks in the flux. See the 
%  removeCosmicRays() function for more info. 
%
% Some definitions of time-scales:
% -buffers: refer to the longest backlog of data. It is mostly used to get
%  the Power Spectra Density (PSD) of each lightcurve. By default this buffer
%  has a length of 5000 to 10,000 frames. 
% -background region is not as long as the full buffer, but contains enough
%  samples to calculate properties such as the mean and STD of each time series. 
%  By default this includes the 2000 frames before (not including) the current
%  batches that are being analyzed (see below:)
% -extended region: usually two adjacent batches (200 frames) that are now
%  being analyzed as one continuous time-span. Filtering and corrections 
%  and correlations are all calculated on this block of data. 
% -search region: a subset of the extended region, including (usually) the 
%  central 100 frames. This region is used (after filtering on the extended
%  region) to find actual event candidates. 
%  These two time-scales allow filtering to be done on a wider region, with 
%  enough margin from the central search region to prevent edge effects. 
%  It means the filtering is always done twice on each time-frame, but each 
%  frame is only ever scanned once for events. 
% -burn-in period: a few thousand frames (default 5000) are used to characterize
%  the different stars and make sure there is a minimal sample of flux in 
%  the buffer to calculate a reliable PSD. During this period no events are 
%  searched for, and all stars are used. 
%  Once the period ends, we disqualify all stars with low S/N (default is 3) 
%  based on their flux during this period. 
%  If multiple apertures have been used in the photometry, we choose the best
%  one based on which aperture left us with the most stars. See note on 
%  multiple apertures below. 
%
% Some definitions of data products:
% -flux: the raw flux from the photometry (often using aperture/forced). 
%        The flux is NOT background subtracted! 
% -aux or auxiliary: other measurements done during the photometry, like
%     background, offsets, width, number of bad pixels, etc. 
%     The nature of the different photometric auxiliary measurements are 
%     described in the util.img.photometry2 documentation. 
% -timestamps: just the internal count from the camera's clock, in seconds. 
% -juldates: translate the local timestamp to global Julian date. 
% -filename: a cell array with the full-path filenames for each frame. This
%            has a lot of redundant info, because frames often belong to one
%            of two files. But we do this to allow different number of batches
%            to be loaded (and thus always have a history of where the data
%            came from no matter how we cut it). 
% -frame_num: what is the corresponding frame index inside each of thos files. 
%             Usually this runs from 1 to 100 and then back to 1 when a new 
%             file is loaded (and the "filename" changes as well). 
% -cutouts: the small images around each star (after calibration) are kept
%           only for the extended region. They are used mainly as visual 
%           aid in inspecting event candidates. 
%
% The flux data is a 2D matrix with time on dim1 and star index on dim2. 
% The aux data is a 3D matrix with dim1 and 2 like fluxes, but with an added
% dim3 that tracks which type of auxiliary it is. The types of aux data are 
% listed in the "aux_names" cell array, and can be indexed using "cut_indices":
% >> widths=store.auxiliary(:,:,store.aux_indices.widths); 
% This organization of the auxiliary data is kept throughout this package. 
% 
% The threshold and the burn-in period: 
% At the end of the burn-in period we check the S/N of each star based on 
% the flux and aux in the buffer. The threshold=3 is used to dump bad stars. 
% The indices of stars that are kept are store in "star_indices" in order
% to recover the original star index from the list of stars that have been 
% loaded from file (e.g., to match a star to its catalog entry). 
% The fluxes, aux, and other data in the store are filtered and only the 
% data for the chosen stars is used from here on out. 
% Stars that later improve or deteriorate are not included/disqualified 
% after the initial calculation is done at the end of the burn-in. 
% 
% A note on multiple apertures:
% If the photometry is done on multiple aperture radii, or using different
% methods like individually aligned apertures vs. forced position apertrues, 
% the store will ingest these additional fluxes and auxiliary with an 
% extra dimension added in the end (dim3 for flux, dim4 for aux). 
% This extra dimension is only used for the burn-in. After the burn-in is 
% finished, when choosing the best stars, the store also chooses the best
% aperture that preserves the most stars above threshold. The index of the 
% chosen aperture is saved in "aperture_index" and defaults to 1 if only a 
% single aperture is ever provided. 
% Like the list of good stars, the choice of aperture does not change during
% the run, after it is decided at the end of the burn-in period. 
% An additional fluxes_extra is used to save additional flux info needed for 
% vetting candidates. Usually this keeps non-forced aperture photometry 
% results (with the same size aperture as the forced) to compare the two
% and rule out cases where false events happen due to aperture misalignment. 
%
% Plotting: use the plot() function to show the data in the buffer.
% The parameters for this function are:
%   -stars: which stars to show (default is 1:5). Stars are counted out of 
%           the stars that passed the cut (after the burn in period). 
%   -aux: specify which kind of aux data you want to see. The default is []
%         which means show the flux and not the aux. 
%         Use either the aux index or name. 
%   -axes: which axes to plot into. 
    
    properties(Transient=true)
        
    end
    
    properties % objects
        
        head@head.Header; 
        checker@tno.QualityChecker; % check the quality of the data
        this_input@util.text.InputVars; % keep a record of the last inputs parsed by this object
        pars; % a struct with all the user-defined parameters
        
    end
    
    properties % inputs/outputs
        
        timestamps; % keep a journal of the timestamps for the entire run
        
        cutouts; % keep cutouts ONLY for the extended region
        
        flux_buffer; % biggest flux buffer used for PSD (we cut from it smaller buffers)
        detrend_buffer; % fluxes after removing a linear fit to each batch
        aux_buffer; % derive all other aux buffers from this (by default this buffer is just equal to the background_aux buffer, but we may change that in the future)
        timestamps_buffer; % timestamps for the flux buffer
        juldates_buffer; % julian dates for each flux measurement
        filename_buffer = {}; % full path + filename for each measurement in the buffer
        frame_num_buffer = []; % frame inside each file
        
        background_flux; % cut out the background flux for calculating the variance outside the search region
        background_detrend; % flux without the linear fit
        background_aux; % cut out the background aux for calculating the variance outside the search region
        background_timestamps; % timestamps for the background region
        background_juldates; % for each timestamp calculate the julian date
        background_filenames = {}; % a cell array the same length as the timestamps, with the filename of the source of each measurement
        background_frame_num = []; % frame inside each file
        
        extended_flux; % cutout of the flux from the flux_buffer extended around the search region
        extended_detrend; % flux without the linear fit
        extended_aux; % cutout of the aux from the aux_buffer extended around the search region 
        extended_timestamps; % timestamps for the extended batch region
        extended_juldates; % for each timestamp calculate the julian date
        extended_filenames = {}; % a cell array the same length as the timestamps, with the filename of the source of each measurement
        extended_frame_num = []; % frame inside each file
        extended_fluxes_extra; % any other apertures we want to save (e.g., unforced aperture, gaussian) or other sizes of apertures, for quality assurance of events
        extended_average_offsets; % the (flux weighted) x/y offsets used for forced photometery
        
        search_start_idx; % starting index for the search region out of the EXTENDED BATCH!
        search_end_idx;  % end index for the search region out of the EXTENDED BATCH!
        
        search_flux; % the flux in the search region
        search_detrend; % flux without the linear fit
        search_aux; % the aux in the search region
        search_timestamps;  % timestamps for the search region
        search_juldates; % juldates for the search region
        search_filenames = {}; % filename cell array for each measurement in the search region
        search_frame_num = [];  % frame inside each file
        
        aux_names = {'errors', 'areas', 'backgrounds', 'variances', 'offsets_x', 'offsets_y', 'centroids_x', 'centroids_y', 'widths', 'bad_pixels', 'flags'}; % add more aux measurements if you want! 
        aux_indices; % struct with field=number for each of the above aux names
        
        frame_counter = 0; % how many frames were entered into the data store since the last reset()
        
        star_indices = []; 
        star_snr = [];
        star_sizes = []; 
        star_rates =[]; 
        
        aperture_index = [];
        
        fwhm_log; % in arcsec!
        juldates_log; % the julian dates for the middle of each batch
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
        
        size_snr_coeffs; % coefficients for a fit of stellar size to the stellar S/N (R = C(0) + C(1).*S + C(2).*S.^2 ...
        
        extra_fluxes_indices; % vector of indices that we want to keep as extras
        
        bad_ratios; % if there is a big difference between the flux in different apertures we disqualify those stars too (e.g., binaries)
        
        version = 1.03;
        
    end
    
    methods % constructor
        
        function obj = DataStore(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'tno.DataStore')
                if obj.debug_bit>1, fprintf('DataStore copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('DataStore constructor v%4.2f\n', obj.version); end
            
                obj.checker = tno.QualityChecker;
                
                obj.resetPars;

            end
            
        end
        
    end
    
    methods % reset/clear
        
        function resetPars(obj)
           
            obj.pars = struct;
            
            % number of frames to keep for each type of calculation
            obj.pars.length_burn_in = 2000; % this mininal number of frames is needed to start event finding! 
            obj.pars.length_psd = 2000; % this many flux measurements used for calculating the PSD
            obj.pars.length_background = 2000; % this many flux/aux measurements used for calculating the variance of the background (outside the search region)
            obj.pars.length_extended = 200; % extended region around the search region, used for overlap, filter edges, etc. 
            obj.pars.length_search = 100; % the area searched in each batch (ideally equal to each batch so there's no need to search the same region twice)

            obj.pars.use_threshold = true; % use the minimal S/N to keep out bad stars and not even store them
            obj.pars.threshold = 3; % stars with S/N lower than this are disqualified after the burn-in period!

            obj.pars.extra_flux_types = {'aperture'}; 
            
            obj.pars.use_remove_cosmic_rays = true; % get rid of intense cosmic rays before even inputting the flux into the buffers
            obj.pars.cosmic_ray_threshold = 8; % in units of S/N

            obj.pars.use_reject_gaussian_photometry = true; 
            obj.pars.use_reject_aperture_photometry = true; 
            
            obj.pars.use_fwhm_per_star = true;
            
            obj.reset;
            
        end
        
        function reset(obj)
            
            obj.checker.reset;
            
            obj.timestamps = [];
            obj.flux_buffer = [];
            obj.detrend_buffer = [];
            obj.aux_buffer = [];
            obj.timestamps_buffer = [];
            obj.juldates_buffer = [];
            obj.filename_buffer = {};
            obj.frame_num_buffer = []; 
            obj.frame_counter = 0; 
            obj.extended_fluxes_extra = [];
            obj.calcAuxIndices;
            
            obj.star_indices = []; 
            obj.star_snr = [];
            obj.star_sizes = [];
            obj.star_rates  =[]; 
            
            obj.aperture_index = []; 
            
            obj.cutouts = [];
            
            obj.juldates_log = [];
            obj.fwhm_log = [];
            
            obj.extra_fluxes_indices = []; 
            
            obj.clear;
            
        end
        
        function clear(obj)
        
            obj.background_flux = [];
            obj.background_detrend = [];
            obj.background_aux = [];
            obj.background_timestamps = [];
            obj.background_juldates = []; 
            obj.background_filenames = {};
            obj.background_frame_num = [];
            
            obj.extended_flux = [];
            obj.extended_detrend = [];
            obj.extended_aux = [];
            obj.extended_timestamps = [];
            obj.extended_juldates = [];
            obj.extended_filenames = {}; 
            obj.extended_frame_num = [];
            
            obj.search_start_idx = [];
            obj.search_end_idx = [];
            
            obj.search_flux = [];
            obj.search_detrend = [];
            obj.search_aux = []; 
            obj.search_timestamps = [];
            obj.search_juldates = [];
            obj.search_filenames = {}; 
            obj.search_frame_num = []; 
            
        end
            
        function calcAuxIndices(obj)
            
            obj.aux_indices = struct;
            
            for ii = 1:length(obj.aux_names)
                obj.aux_indices.(obj.aux_names{ii}) = ii; 
            end
            
        end
        
    end
    
    methods % getters
        
        function val = is_done_burn(obj)
            
            val = obj.frame_counter>obj.pars.length_burn_in;
            
        end
        
        function val = is_reduced_stars(obj)
            
            if isempty(obj.flux_buffer)
                val = 0;
            elseif isempty(obj.star_indices)
                val = 0;
            elseif size(obj.flux_buffer,2)==length(obj.star_indices)
                val = 1;
            else 
                val = 0; 
            end
            
        end
        
        function val = aperture_radius(obj) % to be depricated
            
            if isempty(obj.head.PHOT_PARS)
                val = [];
            else
                val = obj.head.PHOT_PARS.aperture_radius(obj.aperture_index); % this is no longer used, the aperture index is from all types of aperture, including gaussian (see aperture_type())
            end
            
        end
        
        function val = aperture_type(obj)
            
            if isempty(obj.head.PHOT_PARS) || ~isfield(obj.head.PHOT_PARS, 'types') ||isempty(obj.head.PHOT_PARS.types) || isempty(obj.aperture_index)
                val = '';
            else
                val = obj.head.PHOT_PARS.types{obj.aperture_index}; 
            end
            
        end
        
    end
    
    methods % setters
        
        function set.head(obj, val) % setting the header cascades down to sub-objects
            
            obj.head = val;
            
            obj.checker.head = val;
            
        end
        
    end
    
    methods % calculations
        
        function val = makeInputVars(obj) % make an input object that can parse all photometric products and meta-data
            
            val = util.text.InputVars;
            val.use_ordered_numeric = 1; % can just give the numeric values in the correct order, without naming each one
            
            % these inputs are used for explicitely giving the photometric products
            val.input_var('timestamps', []); 
            val.input_var('fluxes', []);
            val.input_var('detrend', []); 
            val.input_var('cutouts', []); 
            
            for ii = 1:length(obj.aux_names)
                if ~strcmp(obj.aux_names{ii}, 'fluxes')
                    val.input_var(obj.aux_names{ii}); 
                end
            end
            
            val.input_var('offset_averages', []); 
            val.input_var('filename', ''); 
            val.input_var('juldates', [], 'julian_dates'); 
            val.input_var('pars_struct', []); 
            
        end
        
        function loadFile(obj, filename) % load data from a Lightcurves file
            
            input = obj.makeInputVars; 
            
            list = input.list_added_properties;
            
            for ii = 1:length(list)
            
                input.(list{ii}) = h5read(filename, ['/' list{ii}]); 
                
            end
            
            obj.input(input); 
            
        end
        
        function input(obj, varargin)
            
            if isempty(varargin)
                error('Must provide photometric products like timestamps, fluxes and so on'); 
            elseif isa(varargin{1}, 'util.text.InputVars')
                obj.this_input = util.oop.full_copy(varargin{1}); % just take the input object some one else provided
                % what happens to any other varargin pairs? 
            elseif isa(varargin{1}, 'img.Photometry') % giving a photometry object is the prefered method! 
                
                obj.this_input = obj.makeInputVars;
                
                list = obj.this_input.list_added_properties;

                for ii = 1:length(list)
                    if isprop(varargin{1}, list{ii})
                        obj.this_input.(list{ii}) = varargin{1}.(list{ii}); 
                    end
                end
                
            else % just make an input object and parse the varargin pairs
                obj.this_input = obj.makeInputVars;
                obj.this_input.scan_vars(varargin{:}); 
            end
            
            % make sure we got fluxes, timestamps and all the aux data
            if isempty(obj.this_input.timestamps) || isempty(obj.this_input.fluxes) || ...
                    isempty(obj.this_input.errors) || isempty(obj.this_input.areas) ||...
                    isempty(obj.this_input.backgrounds) || isempty(obj.this_input.variances) || ...
                    isempty(obj.this_input.offsets_x) || isempty(obj.this_input.offsets_y) || ...
                    isempty(obj.this_input.centroids_x) || isempty(obj.this_input.centroids_y) || ...
                    isempty(obj.this_input.widths) || isempty(obj.this_input.bad_pixels) || ...
                    isempty(obj.this_input.flags)

                error('Some photometric products are empty!'); 

            end
            
            if length(varargin)>1 && isa(varargin{2}, 'img.ModelPSF') % can also pick up the FWHM from the ModelPSF
                if obj.pars.use_fwhm_per_star
                    obj.checker.fwhm = varargin{2}.fwhm_per_star'; 
                    if ~isempty(obj.star_indices)
                        obj.checker.fwhm = obj.checker.fwhm(1, obj.star_indices); % keep only the chosen stars
                    end
                else
                    obj.checker.fwhm = varargin{2}.fwhm; 
                end
                obj.fwhm_log = vertcat(obj.fwhm_log, varargin{2}.fwhm); 
                obj.juldates_log = vertcat(obj.juldates_log, obj.head.get_juldates(mean(obj.this_input.timestamps))); % keep a log of the julian dates of the middle of each batch as well
            end
            
            if ~isempty(obj.head) && ~isempty(obj.this_input.pars_struct)
                obj.head.PHOT_PARS = obj.this_input.pars_struct; % update the header with the photometric parameters used in the analysis
            end
            
            if ndims(obj.this_input.fluxes)>3
                error('This class cannot handle more than 3D fluxes!'); 
            end
            
            if ndims(obj.this_input.widths)>3
                error('This class cannot handle more than 3D auxiliary measurements!'); 
            end
            
            obj.clear; % get rid of the data from last batch
            
            obj.frame_counter = obj.frame_counter + size(obj.this_input.fluxes,1); % how many frames we processed
            
            obj.storeExtraFluxes; % take the extra fluxes before we get rid of them
                
            if obj.pars.use_threshold && obj.is_done_burn % we need to keep only good stars / aperture at the input level
                
                if isempty(obj.star_indices) % only do this once, at the first batch after the burn-in is done
                    obj.calcGoodStarsAndApertures; % this produces the star indices and finds the best aperture
                    obj.parseExtraFluxes; % find the indices for extra fluxes we want to save
                    obj.dumpBadStarsAndAperturesFromBuffer; % this reduces the flux and aux buffers to only the good stars and the best aperture
                end
                
                obj.this_input.fluxes = obj.this_input.fluxes(:,obj.star_indices,obj.aperture_index); % get rid of unwanted stars and apertures
                
                obj.this_input.cutouts = obj.this_input.cutouts(:,:,:,obj.star_indices); % get rid of unwanted stars and apertures
                
                for ii = 1:length(obj.aux_names)
                    new_data = obj.this_input.(obj.aux_names{ii});
                    obj.this_input.(obj.aux_names{ii}) = new_data(:,obj.star_indices,obj.aperture_index); % get rid of unwanted stars and apertures
                end
                
            elseif ~isempty(obj.aperture_index) % if not using the burn in and threshold, we must define the aperture index manually after calling reset()
                
                obj.this_input.fluxes = obj.this_input.fluxes(:,:,obj.aperture_index); 
                
                for ii = 1:length(obj.aux_names)
                    new_data = obj.this_input.(obj.aux_names{ii});
                    obj.this_input.(obj.aux_names{ii}) = new_data(:,:,obj.aperture_index); 
                end
                
            end % we need to keep only good stars / aperture at the input level
            
            if obj.pars.use_remove_cosmic_rays
                obj.removeCosmicRays;
            end
            
            obj.removeFluxNans; % we may as well interpolate the NaNs right here...
            
            obj.removeLinearFit; % put the detrended fluxes in this_input
            
            obj.storeBuffers; % store the new fluxes and timestamps
            
            obj.storeCutouts; % store the cutouts for the extended region only
            
            obj.storeAuxiliary; % store the new aux data
            
            obj.calcSubBuffers; % cut the background, extended and search regions from the relevant buffers
            
            if obj.is_done_burn
                obj.checker.input(obj); % feed the data into the quality checker only after the burn-in is done
            end
            
        end
        
        function calcSubBuffers(obj) % cut the background, extended and search regions from the relevant buffers
            
            % the extended batch region reaches from the end of the flux buffer back "length_extended"
            extended_start_idx = max(1, size(obj.flux_buffer,1)-obj.pars.length_extended+1);
            extended_start_idx_aux = max(1, size(obj.aux_buffer,1)-obj.pars.length_extended+1);
            obj.extended_flux = obj.flux_buffer(extended_start_idx:end,:,:); 
            obj.extended_detrend = obj.detrend_buffer(extended_start_idx:end,:,:); 
            obj.extended_aux = obj.aux_buffer(extended_start_idx_aux:end,:,:,:);
            obj.extended_timestamps = obj.timestamps_buffer(extended_start_idx:end); 
            obj.extended_juldates = obj.juldates_buffer(extended_start_idx:end); 
            obj.extended_filenames = obj.filename_buffer(extended_start_idx:end); 
            obj.extended_frame_num= obj.frame_num_buffer(extended_start_idx:end); 
            
            % the search region is defined in the middle of the extended batch
            margins = floor((obj.pars.length_extended - obj.pars.length_search)/2); 
            
            obj.search_start_idx = margins + 1; 
            obj.search_end_idx = size(obj.extended_flux,1) - margins; 
            obj.search_flux = obj.extended_flux(obj.search_start_idx:obj.search_end_idx,:,:); % cut the search region out of the extended batch
            obj.search_detrend = obj.extended_detrend(obj.search_start_idx:obj.search_end_idx,:,:); % cut the search region out of the extended batch
            
            obj.search_timestamps = obj.extended_timestamps(obj.search_start_idx:obj.search_end_idx);
            obj.search_juldates = obj.extended_juldates(obj.search_start_idx:obj.search_end_idx);
            obj.search_filenames = obj.extended_filenames(obj.search_start_idx:obj.search_end_idx); 
            obj.search_frame_num = obj.extended_frame_num(obj.search_start_idx:obj.search_end_idx); 
            
            obj.search_aux = obj.extended_aux(obj.search_start_idx:obj.search_end_idx,:,:,:); % cut the search region out of the extended batch
            
            % the background region goes back from the start of extended region up to "length_background" before that
            background_start_idx = max(1, extended_start_idx-obj.pars.length_background);
            background_end_idx = extended_start_idx-1;
            obj.background_flux = obj.flux_buffer(background_start_idx:background_end_idx,:,:);
            obj.background_detrend = obj.detrend_buffer(background_start_idx:background_end_idx,:,:);

            obj.background_timestamps = obj.timestamps_buffer(background_start_idx:background_end_idx);
            obj.background_juldates = obj.juldates_buffer(background_start_idx:background_end_idx);
            obj.background_filenames = obj.filename_buffer(background_start_idx:background_end_idx);
            obj.background_frame_num = obj.frame_num_buffer(background_start_idx:background_end_idx);
            
            obj.background_aux = obj.aux_buffer(background_start_idx:background_end_idx,:,:,:); % the size of the aux buffer is just big enough to get the background+extended regions
            
        end
        
        function calcGoodStarsAndApertures(obj) % produces the star indices and finds the best aperture
            
            f = obj.flux_buffer; 
            b = permute(obj.aux_buffer(:,:,obj.aux_indices.backgrounds,:),[1,2,4,3]); % We permute to move dim4 to dim3, the aperture number
            a = permute(obj.aux_buffer(:,:,obj.aux_indices.areas,:),[1,2,4,3]); % We permute to move dim4 to dim3, the aperture number
            
            % permute puts star number in dim1 and aperture number in dim2
            S = permute(nanmean(f - a.*b), [2,3,1]); % the mean flux is taken minus the background
            N = permute(nanmedian(util.series.binning(obj.detrend_buffer, obj.pars.length_extended ,'func', 'std')), [2,3,1]); % for noise calculations we don't subtract background, instead we remove the trends
            % This is the same as calculating the RMS of each extended region (the median ignores outlier batches)
            
            obj.star_snr = S./N; % dim1 is star index, dim2 is aperture index
            
            passed = obj.star_snr >= obj.pars.threshold; 
            
            num_passed = sum(passed,1); % count how many stars are above threshold in each aperture number
            
            if obj.pars.use_reject_gaussian_photometry
                not_gaussian = cellfun(@isempty, regexp(obj.head.PHOT_PARS.types, 'gaussian.*')); % which aperture indices refer to gaussian photometry
                num_passed = num_passed.*not_gaussian; % get rid of these apertures
            end
            
            if obj.pars.use_reject_aperture_photometry
                not_aperture = cellfun(@isempty, regexp(obj.head.PHOT_PARS.types, 'aperture.*')); % which aperture indices refer to aperture photometry
                num_passed = num_passed.*not_aperture; % get rid of these apertures
            end
            
            [~,obj.aperture_index] = max(num_passed, [], 2); % pick the best aperture
            
            passed = passed(:,obj.aperture_index); % keep only the stars on the column of the best aperture
            
            obj.bad_ratios = false(size(S,1),1); 
            
            % compare the flux signal for the picked aperture and the smallest aperture, get rid of stars with large ratios
            if obj.aperture_index>1
                r = S(:,obj.aperture_index)./S(:,1);
                obj.bad_ratios = obj.bad_ratios | (r-nanmedian(r))./nanstd(r) > 3;
            end
            
            % compare the flux signal for the picked aperture and the biggest aperture, get rid of stars with large ratios
            if obj.aperture_index<size(S,2)
                r = S(:,end)./S(:,obj.aperture_index);
                obj.bad_ratios = obj.bad_ratios | (r-nanmedian(r))./nanstd(r) > 3;
            end
            
            passed = passed & ~obj.bad_ratios; 
            
            obj.star_snr = obj.star_snr(:,obj.aperture_index); % keep the S/N from the chosen aperture only
                        
            obj.star_indices = find(passed)'; 
            
            % also get the fit of stellar S/N to size
            s = obj.star_snr;
            R = obj.star_sizes; 
            R(s<0) = [];
            s(s<0) = []; 
            
            fr = util.fit.polyfit(s, R, 'order', 2, 'sigma', 3); 
            
            obj.size_snr_coeffs = fr.coeffs; 
            
        end
        
        function dumpBadStarsAndAperturesFromBuffer(obj) % reduces the flux and aux buffers to only the good stars and the best aperture
            
            obj.clear;
            obj.flux_buffer = obj.flux_buffer(:,obj.star_indices,obj.aperture_index); 
            obj.detrend_buffer = obj.detrend_buffer(:,obj.star_indices,obj.aperture_index); 
            obj.aux_buffer = obj.aux_buffer(:,obj.star_indices,:,obj.aperture_index); 
            obj.cutouts = obj.cutouts(:,:,:,obj.star_indices); 
            
            obj.extended_fluxes_extra = obj.extended_fluxes_extra(:,obj.star_indices,obj.extra_fluxes_indices); 
            
            if size(obj.checker.fwhm, 2)>1
                obj.checker.fwhm = obj.checker.fwhm(1, obj.star_indices); % keep only the chosen stars
            end
            
        end
        
        function parseExtraFluxes(obj)
            
            obj.extra_fluxes_indices = []; 
            
            for ii = 1:length(obj.pars.extra_flux_types)
                
                s = obj.pars.extra_flux_types{ii}; 
                
                if ischar(s)

                    idx = regexp(s, '\d+');

                    if isempty(idx) % no number given, use the same number as the main flux
                        
                        % get the currently used aperture size
                        ap_size = util.text.extract_numbers(obj.head.PHOT_PARS.types);
                        ap_size = ap_size{obj.aperture_index}; 
                        
                        s = strip(s); % remove whitespace
                        s = sprintf('%s %4.2f', s, ap_size);
                        
                        ap_index = find(strcmp(obj.head.PHOT_PARS.types, s)); 
                        
                        obj.extra_fluxes_indices(end+1,1) = ap_index;
                        
                    elseif idx==1 % number is given, without any text, use this as an index
                        obj.extra_fluxes_indices(end+1,1) = str2double(s); % just assign it as an index
                    else
                        ap_index = find(strcmp(obj.head.PHOT_PARS.types, s)); 
                        obj.extra_fluxes_indices(end+1,1) = ap_index;
                    end

                elseif isnumeric(s)
                    obj.extra_fluxes_indices(end+1,1) = s; % just assign it as an index
                end
                
            end
            
        end
        
        function removeCosmicRays(obj) % gets rid of isolated spikes in the flux
            
            low_thresh = 3.5; % set the limit for the neighbors of cosmic rays to be some low value like 3.5 sigma
                
            f = obj.this_input.fluxes;

            noise = util.stat.rstd(f); 
            average = nanmedian(f); 

            f2 = (f-average)./noise; % normalized flux

            idx = f2>obj.pars.cosmic_ray_threshold; % mark regions where there is a cosmic ray
            list = find(idx); % linear indices in a vector

            for ii = 1:length(list) % go over the list and check these are individual detections

                [frame, star, ap] = ind2sub(size(idx), list(ii));

                if frame==1 % first frame of this input
                    if f2(frame+1,star,ap)>low_thresh % the very next frame is also above a few sigma, so this doesn't look like a cosmic ray
                        idx(frame, star,ap) = false; % this is no longer seen as a cosmic ray index
                    end
                elseif frame==size(f,1) % last frame of this input
                    if f2(frame-1,star,ap)>low_thresh % the previous frame is also above a few sigma, so this doesn't look like a cosmic ray
                        idx(frame, star,ap) = false; % this is no longer seen as a cosmic ray index
                    end
                else
                    if f2(frame+1,star,ap)>low_thresh && f2(frame-1,star,ap)>low_thresh % both of the two nearest neighbors is also above a few sigma, so this doesn't look like a cosmic ray
                        idx(frame, star,ap) = false; % this is no longer seen as a cosmic ray index
                    end
                end

            end % for ii in list

            obj.this_input.fluxes(idx) = NaN; % replace the cosmic rays frames with NaNs 
            
        end
        
        function removeFluxNans(obj)
            
            obj.this_input.fluxes = fillmissing(obj.this_input.fluxes, 'spline'); 
            
        end
        
        function removeLinearFit(obj)
            
            f = obj.this_input.fluxes; 

            obj.this_input.detrend = util.series.detrend(f, 'iterations', 2); 
            
%             for ii = 1:size(f,3)
%                 
%                 fit_result = util.fit.polyfit(1:size(f,1), f(:,:,ii), 'order', 1, 'double', 1, 'sigma', 3, 'iterations', 2); 
% 
%                 f2 = f(:,:,ii) - [fit_result.ym]; % remove the fit
% 
%                 obj.this_input.detrend(:,:,ii) = f2; 
% 
%             end
            
        end
        
        function data_out = appendData(obj, data_existing, data_new, buffer_size) % append buffer data inside a try/catch to protect against temporary memory shortages
            
            if nargin<4 || isempty(buffer_size)
                buffer_size = obj.pars.length_psd;
            end
            
            try
            
                data_out = vertcat(data_existing, data_new); 

                if size(data_out,1)>buffer_size
                    data_out = data_out(end-buffer_size+1:end,:,:,:,:); 
                end

            catch ME
                
                if strcmp(ME.identifier, 'MATLAB:nomem')

                    util.text.date_printf('Out of memory while appending data! Trying again...'); 
                    
                    pause(30); % let the memory from other runs free up
                    
                    data_out = vertcat(data_existing, data_new); 

                    if size(data_out,1)>buffer_size
                        data_out = data_out(end-buffer_size+1:end,:,:,:,:); 
                    end

                else
                    rethrow(ME); 
                end
                
            end
            
        end
        
        function storeBuffers(obj)
            
            obj.flux_buffer = obj.appendData(obj.flux_buffer, single(obj.this_input.fluxes));
            obj.detrend_buffer = obj.appendData(obj.detrend_buffer, single(obj.this_input.detrend)); 
            obj.timestamps = vertcat(obj.timestamps, single(obj.this_input.timestamps)); % this tracks timestamps for the entire run
            obj.timestamps_buffer = obj.appendData(obj.timestamps_buffer, single(obj.this_input.timestamps));
            obj.juldates_buffer = obj.appendData(obj.juldates_buffer, obj.this_input.juldates); 
            obj.filename_buffer = obj.appendData(obj.filename_buffer, repmat({obj.this_input.filename}, [length(obj.this_input.timestamps), 1])); 
            obj.frame_num_buffer = obj.appendData(obj.frame_num_buffer, single(1:length(obj.this_input.timestamps))'); % assume the frame numbers are just run continuously from 1->number of frames in batch
            
            obj.extended_average_offsets = obj.appendData(obj.extended_average_offsets, obj.this_input.offset_averages, obj.pars.length_extended);
            
%             if size(obj.flux_buffer,1)>obj.pars.length_psd % only save the recent data
%                 obj.flux_buffer = obj.flux_buffer(end-obj.pars.length_psd+1:end,:,:); 
%                 obj.detrend_buffer = obj.detrend_buffer(end-obj.pars.length_psd+1:end,:,:); 
%                 obj.timestamps_buffer = obj.timestamps_buffer(end-obj.pars.length_psd+1:end); 
%                 obj.juldates_buffer = obj.juldates_buffer(end-obj.pars.length_psd+1:end); 
%                 obj.filename_buffer = obj.filename_buffer(end-obj.pars.length_psd+1:end); 
%                 obj.frame_num_buffer = obj.frame_num_buffer(end-obj.pars.length_psd+1:end); 
%             end
         
        end
        
        function storeCutouts(obj)
            
            try
            obj.cutouts = cat(3, obj.cutouts, obj.this_input.cutouts);
            catch ME
                if strcmp(ME.identifier, 'MATLAB:nomem')
                    util.text.date_printf('Out of memory while storing cutouts! Trying again...'); 
                    pause(30); 
                    obj.cutouts = cat(3, obj.cutouts, obj.this_input.cutouts);
                else    
                    rethrow(ME);
                end
                
            end
            
            try 
             
                if size(obj.cutouts,3)>obj.pars.length_extended
                    obj.cutouts = obj.cutouts(:,:,end-obj.pars.length_extended+1:end,:); 
                end
            
            catch ME
                if strcmp(ME.identifier, 'MATLAB:nomem')
                    util.text.date_printf('Out of memory while reducing size of cutouts! Trying again...'); 
                    pause(30); 
                    if size(obj.cutouts,3)>obj.pars.length_extended
                        obj.cutouts = obj.cutouts(:,:,end-obj.pars.length_extended+1:end,:); 
                    end
                else
                    rethrow(ME); 
                end
            end
            
        end
        
        function storeAuxiliary(obj)
            
            list = obj.aux_names;
            new_aux = NaN(size(obj.this_input.errors,1), size(obj.this_input.widths,2), length(list), size(obj.this_input.widths,3), 'like', obj.this_input.widths); % preallocate
            
            for ii = 1:length(list)
                new_aux(:,:,ii,:) = permute(obj.this_input.(list{ii}), [1,2,4,3]); % allow multiple apertures (3D aux matrices from photometry) 
            end
            
            obj.aux_buffer = obj.appendData(obj.aux_buffer, new_aux); 
            
%             obj.aux_buffer = vertcat(obj.aux_buffer, new_aux); % add the new auxiliary data
% 
%             if size(obj.aux_buffer,1)>obj.pars.length_psd % make sure to only save the recent data
%                 obj.aux_buffer = obj.aux_buffer(end-(obj.pars.length_psd)+1:end,:,:,:); 
%             end
            
        end
        
        function storeExtraFluxes(obj)
            
            star_idx = obj.star_indices;
            if isempty(star_idx)
                star_idx = 1:size(obj.flux_buffer,2); 
            end
            
            ap_idx = obj.extra_fluxes_indices;
            if isempty(ap_idx)
                ap_idx = 1:size(obj.flux_buffer,3);
            end
            
            if isempty(obj.extended_fluxes_extra)
                obj.extended_fluxes_extra = obj.this_input.fluxes(:, star_idx, ap_idx); 
            else
                obj.extended_fluxes_extra = vertcat(obj.extended_fluxes_extra, obj.this_input.fluxes(:, star_idx, ap_idx)); 
            end
            
            if size(obj.extended_fluxes_extra,1)>obj.pars.length_extended
                obj.extended_fluxes_extra = obj.extended_fluxes_extra(end-obj.pars.length_extended+1:end,:,:); 
            end
            
        end
        
        function saveHours(obj, bad_batch) % store the star hours from this batch in the checker.hours object. The bad_batch argument tells the StarHours to disqualify the entire batch (default is 0)
            
            if nargin<3 || isempty(bad_batch)
                bad_batch = 0;
            end
            
            obj.checker.hours.input(obj.checker, bad_batch); 
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function plot(obj, varargin) % show the data fromt the flux/aux buffers for a subset of stars. 
            
            input = util.text.InputVars;
            input.input_var('stars', 1:5); % which stars to plot
            input.input_var('aux', []); % which aux indices to plot (leave empty to not plot them at all)
            input.input_var('axes', [], 'axis'); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            next_state = input.axes.NextPlot;
            
            plot(input.axes, obj.timestamps_buffer, obj.flux_buffer(:,input.stars,1), '-'); 
            
            input.axes.NextPlot = 'add'; 
            
            input.axes.ColorOrderIndex = 1;
            plot(input.axes, obj.background_timestamps, obj.background_flux(:,input.stars,1), 's'); 
            
            input.axes.ColorOrderIndex = 1;
            plot(input.axes, obj.extended_timestamps, obj.extended_flux(:,input.stars,1), 'v'); 
            
            input.axes.ColorOrderIndex = 1;
            plot(input.axes, obj.search_timestamps, obj.search_flux(:,input.stars,1), '^'); 
            
            if ~isempty(input.aux)
            
                if ischar(input.aux)
                    input.aux = obj.aux_indices.(input.aux); % turn a string into an index
                end
                
                input.axes.ColorOrderIndex = 1;
                plot(input.axes, obj.background_timestamps, obj.background_aux(:,input.stars,input.aux), '.'); 

                input.axes.ColorOrderIndex = 1;
                plot(input.axes, obj.extended_timestamps, obj.extended_aux(:,input.stars,input.aux), 'o'); 

                input.axes.ColorOrderIndex = 1;
                plot(input.axes, obj.search_timestamps, obj.search_aux(:,input.stars,input.aux), 'p'); 
            
            end
            
            input.axes.NextPlot = next_state;
            
        end
        
    end    
    
end

