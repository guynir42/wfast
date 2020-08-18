classdef DataStore < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        checker@trig.QualityChecker;
        this_input@util.text.InputVars; 
        
    end
    
    properties % inputs/outputs
        
        timestamps; % keep a journal of the timestamps for the entire run
        
        flux_buffer; % biggest flux buffer used for PSD (we cut from it smaller buffers)
        aux_buffer; % derive all other aux buffers from this (by default this buffer is just equal to the background_aux buffer, but we may change that in the future)
        timestamps_buffer; % timestamps for the flux buffer
        juldates_buffer; % julian dates for each flux measurement
        filename_buffer = {}; % full path + filename for each measurement in the buffer
        frame_num_buffer = []; % frame inside each file
        
        background_flux; % cut out the background flux for calculating the variance outside the search region
        background_aux; % cut out the background aux for calculating the variance outside the search region
        background_timestamps; % timestamps for the background region
        background_juldates; % for each timestamp calculate the julian date
        background_filenames = {}; % a cell array the same length as the timestamps, with the filename of the source of each measurement
        background_frame_num = []; 
        
        extended_flux; % cutout of the flux from the flux_buffer extended around the search region
        extended_aux; % cutout of the aux from the aux_buffer extended around the search region 
        extended_timestamps; % timestamps for the extended batch region
        extended_juldates; % for each timestamp calculate the julian date
        extended_filenames = {}; % a cell array the same length as the timestamps, with the filename of the source of each measurement
        extended_frame_num = [];  % frame inside each file
        
        search_start_idx; % starting index for the search region out of the EXTENDED BATCH!
        search_end_idx;  % end index for the search region out of the EXTENDED BATCH!
        
        search_flux; % the flux in the search region
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
        
    end
    
    properties % switches/controls
        
        % number of frames to keep for each type of calculation
        length_burn_in = 5e3; % this mininal number of frames is needed to start event finding! 
        length_psd = 2e4; % this many flux measurements used for calculating the PSD
        length_background = 2000; % this many flux/aux measurements used for calculating the variance of the background (outside the search region)
        length_extended = 200; % extended region around the search region, used for overlap, filter edges, etc. 
        length_search = 100; % the area searched in each batch (ideally equal to each batch so there's no need to search the same region twice)
        
        use_threshold = true; % use the minimal S/N to keep out bad stars and not even store them
        threshold = 5; % stars with S/N lower than this are disqualified after the burn-in period!
        func_rms = @nanstd; % which function should be used for calculating the noise in S/N 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = DataStore(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.DataStore')
                if obj.debug_bit>1, fprintf('DataStore copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('DataStore constructor v%4.2f\n', obj.version); end
            
                obj.checker = trig.QualityChecker;
                
                obj.reset;

            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.checker.reset;
            
            obj.timestamps = [];
            obj.flux_buffer = [];
            obj.aux_buffer = [];
            obj.timestamps_buffer = [];
            obj.juldates_buffer = [];
            obj.filename_buffer = {};
            obj.frame_num_buffer = []; 
            obj.frame_counter = 0; 
        
            obj.calcAuxIndices;
            
            obj.star_indices = []; 
            obj.star_snr = [];
            
            obj.clear;
            
        end
        
        function clear(obj)
        
            obj.background_flux = [];
            obj.background_aux = [];
            obj.background_timestamps = [];
            obj.background_juldates = []; 
            obj.background_filenames = {};
            obj.background_frame_num = [];
            
            obj.extended_flux = [];
            obj.extended_aux = [];
            obj.extended_timestamps = [];
            obj.extended_juldates = [];
            obj.extended_filenames = {}; 
            obj.extended_frame_num = [];
            
            obj.search_start_idx = [];
            obj.search_end_idx = [];
            
            obj.search_flux = [];
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
            
            val = obj.frame_counter>obj.length_burn_in;
            
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
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function val = makeInputVars(obj)
            
            val = util.text.InputVars;
            val.use_ordered_numeric = 1; % can just give the numeric values in the correct order, without naming each one
            
            % these inputs are used for explicitely giving the photometric products
            val.input_var('timestamps', []); 
            val.input_var('fluxes', []);
            
            for ii = 1:length(obj.aux_names)
                if ~strcmp(obj.aux_names{ii}, 'fluxes')
                    val.input_var(obj.aux_names{ii}); 
                end
            end
            
            val.input_var('filename', ''); 
            val.input_var('juldates', [], 'julian_dates'); 
            
        end
        
        function loadFile(obj, filename)
            
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
                obj.this_input = util.oop.full_copy(varargin{1});
                % what happens to any other varargin pairs?
            elseif isa(varargin{1}, 'img.Photometry')
                obj.this_input = obj.makeInputVars;
                
                list = obj.this_input.list_added_properties;

                for ii = 1:length(list)
                    obj.this_input.(list{ii}) = varargin{1}.(list{ii}); 
                end
                    
            else
                obj.this_input = obj.makeInputVars;
                obj.this_input.scan_vars(varargin{:}); 
            end
            
            if isempty(obj.this_input.timestamps) || isempty(obj.this_input.fluxes) || ...
                    isempty(obj.this_input.errors) || isempty(obj.this_input.areas) ||...
                    isempty(obj.this_input.backgrounds) || isempty(obj.this_input.variances) || ...
                    isempty(obj.this_input.offsets_x) || isempty(obj.this_input.offsets_y) || ...
                    isempty(obj.this_input.centroids_x) || isempty(obj.this_input.centroids_y) || ...
                    isempty(obj.this_input.widths) || isempty(obj.this_input.bad_pixels) || ...
                    isempty(obj.this_input.flags)

                error('Some photometric products are empty!'); 

            end
            
            if ndims(obj.this_input.fluxes)>3
                error('This class cannot handle more than 3D fluxes!'); 
            end
            
            if ndims(obj.this_input.widths)>2
                error('This class cannot handle more than 2D auxiliary measurements!'); 
            end
                
            obj.clear;
            
            obj.frame_counter = obj.frame_counter + size(obj.this_input.fluxes,1); 
            
            if obj.use_threshold && obj.is_done_burn % we need to keep only good stars at the input level
                
                if isempty(obj.star_indices)
                    obj.calcGoodStars; % this produces the star indices
                    obj.dumpBadStarsFromBuffer; % this reduces the flux and aux buffers to only the good stars
                end
                
                obj.this_input.fluxes = obj.this_input.fluxes(:,obj.star_indices,:); 
                
                for ii = 1:length(obj.aux_names)
                    new_data = obj.this_input.(obj.aux_names{ii});
                    obj.this_input.(obj.aux_names{ii}) = new_data(:,obj.star_indices); 
                end
                
            end
            
            % store the new fluxes and timestamps
            obj.flux_buffer = vertcat(obj.flux_buffer, obj.this_input.fluxes);
            obj.timestamps = vertcat(obj.timestamps, obj.this_input.timestamps); % this tracks timestamps for the entire run
            obj.timestamps_buffer = vertcat(obj.timestamps_buffer, obj.this_input.timestamps);
            obj.juldates_buffer = vertcat(obj.juldates_buffer, obj.this_input.juldates); 
            obj.filename_buffer = vertcat(obj.filename_buffer, repmat({obj.this_input.filename}, [length(obj.this_input.timestamps), 1])); 
            obj.frame_num_buffer = vertcat(obj.frame_num_buffer, (1:length(obj.this_input.timestamps))'); % assume the frame numbers are just run continuously from 1->number of frames in batch
            
            if size(obj.flux_buffer,1)>obj.length_psd
                obj.flux_buffer = obj.flux_buffer(end-obj.length_psd+1:end,:,:); 
                obj.timestamps_buffer = obj.timestamps_buffer(end-obj.length_psd+1:end); 
                obj.juldates_buffer = obj.juldates_buffer(end-obj.length_psd+1:end); 
                obj.filename_buffer = obj.filename_buffer(end-obj.length_psd+1:end); 
                obj.frame_num_buffer = obj.frame_num_buffer(end-obj.length_psd+1:end); 
            end
            
            % store the new aux data
            list = obj.aux_names;
            new_aux = NaN(size(obj.this_input.errors,1), size(obj.this_input.errors,2), length(list), 'like', obj.this_input.errors); % preallocate
            
            for ii = 1:length(list)
                new_aux(:,:,ii) = obj.this_input.(list{ii}); 
            end
            
            obj.aux_buffer = vertcat(obj.aux_buffer, new_aux);

            if size(obj.aux_buffer,1)>obj.length_background+obj.length_extended % the aux buffer is smaller than flux buffer: it is only long enough for background+extended region
                obj.aux_buffer = obj.aux_buffer(end-(obj.length_background+obj.length_extended)+1:end,:,:); 
            end

            obj.calcSubBuffers;
            
            if obj.is_done_burn
                obj.checker.input(obj); % feed the data into the quality checker
            end
            
        end
        
        function calcSubBuffers(obj)
            
            % the extended batch region reaches from the end of the flux buffer back "length_extended"
            extended_start_idx = max(1, size(obj.flux_buffer,1)-obj.length_extended+1);
            extended_start_idx_aux = max(1, size(obj.aux_buffer,1)-obj.length_extended+1);
            obj.extended_flux = obj.flux_buffer(extended_start_idx:end,:,:); 
            obj.extended_aux = obj.aux_buffer(extended_start_idx_aux:end,:,:);
            obj.extended_timestamps = obj.timestamps_buffer(extended_start_idx:end); 
            obj.extended_juldates = obj.juldates_buffer(extended_start_idx:end); 
            obj.extended_filenames = obj.filename_buffer(extended_start_idx:end); 
            obj.extended_frame_num= obj.frame_num_buffer(extended_start_idx:end); 
            
            % the search region is defined in the middle of the extended batch
            margins = floor((obj.length_extended - obj.length_search)/2); 
            
            obj.search_start_idx = margins + 1; 
            obj.search_end_idx = size(obj.extended_flux,1) - margins; 
            obj.search_flux = obj.extended_flux(obj.search_start_idx:obj.search_end_idx,:,:); % cut the search region out of the extended batch
            
            obj.search_timestamps = obj.extended_timestamps(obj.search_start_idx:obj.search_end_idx);
            obj.search_juldates = obj.extended_juldates(obj.search_start_idx:obj.search_end_idx);
            obj.search_filenames = obj.extended_filenames(obj.search_start_idx:obj.search_end_idx); 
            obj.search_frame_num = obj.extended_frame_num(obj.search_start_idx:obj.search_end_idx); 
            
            obj.search_aux = obj.extended_aux(obj.search_start_idx:obj.search_end_idx,:,:); % cut the search region out of the extended batch
            
            % the background region goes back from the start of extended region up to "length_background" before that
            background_start_idx = max(1, extended_start_idx-obj.length_background);
            background_end_idx = extended_start_idx-1;
            obj.background_flux = obj.flux_buffer(background_start_idx:background_end_idx,:,:);

            obj.background_timestamps = obj.timestamps_buffer(background_start_idx:background_end_idx);
            obj.background_juldates = obj.juldates_buffer(background_start_idx:background_end_idx);
            obj.background_filenames = obj.filename_buffer(background_start_idx:background_end_idx);
            obj.background_frame_num = obj.frame_num_buffer(background_start_idx:background_end_idx);
            
            background_end_idx_aux = extended_start_idx_aux-1;
            obj.background_aux = obj.aux_buffer(1:background_end_idx_aux,:,:); % the size of the aux buffer is just big enough to get the background+extended regions
            
        end
        
        function calcGoodStars(obj)
            
            S = nanmean(obj.flux_buffer); % signal
            N = obj.func_rms(obj.flux_buffer); % noise
            
            obj.star_snr = S./N;
            obj.star_indices =  find(obj.star_snr >= obj.threshold); 
            
        end
        
        function dumpBadStarsFromBuffer(obj)
            
            obj.clear;
            obj.flux_buffer = obj.flux_buffer(:,obj.star_indices,:); 
            obj.aux_buffer = obj.aux_buffer(:,obj.star_indices,:); 
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function plot(obj, varargin)
            
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

