classdef Event < handle

    properties(Transient=true)
        
        dates; 
        
    end
    
    properties % objects
        
        head@head.Header; % this doesn't have to be filled, but is useful to attach when saving events...
        
        phot_pars; % a struct with some housekeeping about how the photometry was done
        
        fit_pars; % fit to the occultation parameters
        
        sim_pars; % if simulated, this should contain the fit parameters...
        
        gaia_data; % an object with all rows from the GAIA match (a cat.data table row)
        
    end
    
    properties(Dependent=true)
        
        cutouts;
        positions;
        stack;
        batch_index;
        filename;
        
    end
    
    properties % inputs/outputs
        
        timestamps; % all timestamps for the 2-batch window
        flux_raw_all; % all the fluxes in the 2-batch window *
        flux_filtered; % flux of the best star and best kernel, after filtering and normalizing
        flux_detrended; % flux of the best star after removing slow changes (trends)
        std_flux; % noise level of the detrended flux
        flux_fitted; % best fit flux model (to be filled later, adding a occult.Parameters object)
        
        % NOTE: these can be absolute thresholds (positive) or relative to threshold (negative). 
        time_range_thresh; % all times around the peak (for same star/kernel) above this threshold are also part of the event
        kern_range_thresh; % any kernels in the time range above this threshold are also included in the event
        star_range_thresh; % any stars in the time range above this threshold are also included in the event
        
        time_index; % index in the 2-batch time window where the maximum value occured
        kern_index; % index of the best kernel
        star_index; % index of the best star 
        
        time_indices; % what continuous range of timestamp indices passed the time_range_thresh
        kern_indices; % what assorted kernel indices passed the kern_range_thresh inside the time range for one star
        star_indices; % what assorted star indices passed the kern_range_thresh inside the time range for one kernel
        
        peak_timestamp; % timestamp for time_index
        best_kernel; % full kernel that goes with kern_index

        duration;
        star_snr;
        magnitude;
        RA;
        Dec;
        
        which_batch = 'first'; % can be "first" or "second", depending on where is the peak in the 2-batch window
        
        serial; % which event in the run (including unkept events)
        
        snr; % trigger value (after adjusting for "previous_std")
        threshold; % trigger threshold used (in S/N units)
        previous_std; % average rms of previous batches
        is_positive; % is the filter-response positive (typically for occultations)
        
        was_psd_corrected;
        was_var_buffered;
        
        notes = ''; % why this event was not kept
        keep = 1; % only keep the "real" events that are not duplicates
        
    end
    
    properties % switches/controls
        
        corr_flux = 'raw'; 
        max_corr;
        max_num_nans;
        
        debug_bit = 1;
        
    end
    
    properties(Hidden=true)
       
        % many many options for being not-real 
        is_duplicate = 0; % has a better event in the other batch
        is_offset = 0; % the star's position offset is too big to be a real event
        is_global = 0; % this event is triggered on many stars at the same time, suggesting a non-astronomical event
        is_cosmic_ray = 0; % cosmic ray
        is_artefact = 0; % some other sort of image artefact
        is_bad_batch = 0; % this batch is noisy and has lots of events
        is_bad_star = 0; % this star is constantly triggering (black listed)
        is_corr_offsets = 0; % if there is strong correlation to star offsets in the cutout 
        is_corr_bg = 0; % if there is strong anti-correlation to background
        is_corr_area = 0; % if there is strong correlation to aperture area
        is_corr_width = 0; % if there is a strong correlation to PSF width
        is_nan_flux = 0; % this event has too many NaN values in the raw flux around the peak
        is_nan_offsets = 0; % this event has too many NaN values in the offsets of x or y around the peak
        is_simulated = 0; % we've added this event on purpose
        is_forced = 0; % we've forced the system to trigger on a subtrheshold event
        % ...
        
        % additional housekeeping info
        used_background_sub; % did the photometery object use local (e.g., annulus) b/g subtraction?
        time_step;
        
        % The following are additional metadata collected by Photometry class
        % "at_peak" means values for all stars at the time of event peak (row vector)
        % "at_star" means values for this star on all timestamps (column vector)
        % They intersect at star_index and time_index 
        % "time_average" means for each star, averaged over time (row vector)
        % "star_average" means for each timestamp, averaged between stars (column vector)
        areas_at_peak;
        areas_at_star;
        areas_time_average;
        areas_star_average;
        
        backgrounds_at_peak;
        backgrounds_at_star;
        backgrounds_time_average;
        backgrounds_star_average;
        
        variances_at_peak;
        variances_at_star;
        variances_time_average;
        variances_star_average;
        
        offsets_x_at_peak;
        offsets_x_at_star;
        offsets_x_time_average;
        offsets_x_star_average;
        
        offsets_y_at_peak;
        offsets_y_at_star;
        offsets_y_time_average;
        offsets_y_star_average;
        
        widths_at_peak;
        widths_at_star;
        widths_time_average;
        widths_star_average;
        
        bad_pixels_at_peak;
        bad_pixels_at_star;
        bad_pixels_time_average;
        bad_pixels_star_average;
        
        corr_x; % correlation with offsets_x
        corr_y; % correlation with offsets_y
        corr_r; % correlation with offset size 
        corr_b; % correlation with background
        corr_a; % correlation with area
        corr_w; % correlation with PSF width
        
        covar_x; % covariance with offsets_x
        covar_y; % covariance with offsets_y
        covar_r; % covariance with offset size 
        covar_b; % covariance with background
        covar_a; % covariance with area
        covar_w; % covariance with PSF width
        
        aperture; % aperture radius used by photometry (empty if not using aperture photometry)
        gauss_sigma; % gaussian width used by photometry (empty if not using PSF photometry)
        
        % these are chosen based on "which_batch"
        cutouts_first;
        cutouts_second;
        positions_first;
        positions_second;
        stack_first;
        stack_second;
        batch_index_first;
        batch_index_second;
        filename_first;
        filename_second;
        
        t_end;
        t_end_stamp;
        
        version = 1.05;
        
    end
    
    methods % constructor
        
        function obj = Event(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.Event')
                if obj.debug_bit>1, fprintf('Event copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('Event constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function clearImages(obj)
            
            obj.cutouts_first = [];
            obj.cutouts_second = [];
            obj.positions_first = [];
            obj.positions_second = [];
            obj.stack_first = [];
            obj.stack_second = [];
            
        end
        
    end
    
    methods % getters
        
        function val = get.dates(obj)
            
            if isempty(obj.dates) && ~isempty(obj.t_end) && ~isempty(obj.t_end_stamp)
                val = util.text.str2time(obj.t_end, obj.t_end_stamp, obj.timestamps);
            else
                val = obj.dates;
            end
            
        end
        
        function val = get.cutouts(obj)
            
            val = obj.(['cutouts_' obj.which_batch]);
            
        end
        
        function val = get.positions(obj)
            
            val = obj.(['positions_' obj.which_batch]);
            
        end
        
        function val = get.stack(obj)
            
            val = obj.(['stack_' obj.which_batch]);
            
        end
        
        function val = get.batch_index(obj)
            
            val = obj.(['batch_index_' obj.which_batch]);
            
        end
        
        function val = get.filename(obj)
            
            val = obj.(['filename_' obj.which_batch]);
            
        end
        
        function val = peak_datetime(obj)
            
            if isempty(obj.peak_timestamp) && ~isempty(obj.t_end) && ~isempty(obj.t_end_stamp)
                val = util.text.str2time(obj.t_end, obj.t_end_stamp, obj.peak_timestamp);
            else
                val = [];
            end
            
        end
        
        function val = kernel_timestamps(obj)
            
            dt = median(diff(obj.timestamps));
            
            val = (-floor(length(obj.best_kernel)/2):floor(length(obj.best_kernel)/2)).*dt;
            
        end
        
        function val = frame_index(obj) % which frame the peak occured (inside the batch)
            
            val = obj.time_index;
            
            if util.text.cs(obj.which_batch, 'second')
                val = val - size(obj.cutouts_first,3);
            end
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function input(obj, table_row, fluxes, timestamps, kernels, stds) % to be depricated! 
            
            box = table_row{1, 'BoundingBox'};
            idx = table_row{1, 'VoxelIdxList'}{1};
            val = table_row{1, 'VoxelValues'}{1};
            
            if nargin>3 && ~isempty(timestamps)
                obj.timestamps = timestamps; 
            end
            
            obj.range_time_idx = ceil(box(2)):floor(box(2))+box(5);
            obj.range_kernels_idx = ceil(box(1)):floor(box(1))+box(4); 
            obj.range_stars = ceil(box(3)):floor(box(3))+box(6); 
            
            [obj.snr, max_idx] = max(val); 
            [y,x,z] = ind2sub(size(fluxes), idx(max_idx));
            obj.peak_timestamp = timestamps(y);
            
            if nargin>4 && ~isempty(kernels)
                obj.best_kernel = kernels(:,x);
            end
            
            if nargin>5 && ~isempty(stds)
                obj.stds = squeeze(stds);
            end
            
            obj.star_index = z;
            
            obj.flux_filtered = fluxes(:, x, z);
            
            obj.which_batch = 'first';
            obj.frame_index = y;
                        
        end
        
        function val = is_same(obj, other)
            
            val = 0;
            
            for ii = 1:length(other)
                
                if obj.star_index==other(ii).star_index
                    t_self = obj.timestamps(obj.time_indices);
                    t_other = other(ii).timestamps(other(ii).time_indices);
                    if ~isempty(intersect(t_self, t_other))
                        val = 1;
                        return;
                    end
                    
                end
                
            end
            
        end
        
        function self_check(obj)
            
            % check if there are too many NaNs in the raw flux during the event
            NN = nnz(isnan(obj.flux_raw_all(obj.time_indices, obj.star_index)));
            if NN>=obj.max_num_nans
                obj.keep = 0;
                obj.is_nan_flux = 1;
                obj.addNote(sprintf('flux with %d NaNs', NN));
            end
            
            % check if there are too many NaNs in the offsets
            Nidx = isnan(obj.offsets_x_at_star(obj.time_indices)) | isnan(obj.offsets_y_at_star(obj.time_indices)); 
            if nnz(Nidx)>=obj.max_num_nans % && ~all(Nidx) % if all frames have NaN offsets (i.e., forced photometry), it is ok. If some of them do, it is bad shape
                obj.keep = 0;
                obj.is_nan_offsets = 1;
                obj.addNote(sprintf('offsets with %d NaNs', nnz(Nidx)));
            end
            
            % check for correlation with offsets_x (negative or positive)
            [obj.corr_x, obj.covar_x] = obj.correlation(obj.offsets_x_at_star);
            if abs(obj.corr_x)>obj.max_corr
                obj.keep = 0;
                obj.is_corr_offsets = 1;
                obj.addNote(sprintf('correlation->offsets x= %f', obj.corr_x));
            end
            
            % check for correlation with offsets_y (negative or positive)
            [obj.corr_y, obj.covar_y] = obj.correlation(obj.offsets_y_at_star);
            if abs(obj.corr_y)>obj.max_corr
                obj.keep = 0;
                obj.is_corr_offsets = 1;
                obj.addNote(sprintf('correlation->offsets y= %f', obj.corr_y));
            end
            
            x = obj.offsets_x_at_star;
            y = obj.offsets_y_at_star;
            r = sqrt(x.^2+y.^2); 
            
            % check for correlation with size of offsets (x and y combined)
            [obj.corr_r, obj.covar_r] = obj.correlation(r);
            if abs(obj.corr_r)>obj.max_corr
                obj.keep = 0;
                obj.is_corr_offsets = 1;
                obj.addNote(sprintf('correlation->offset size= %f', obj.corr_r));
            end
                        
            % check for negative correlation with the background
            [obj.corr_b, obj.covar_b] = obj.correlation(obj.backgrounds_at_star);
            if obj.corr_b<-obj.max_corr
                obj.keep = 0;
                obj.is_corr_bg = 1;
                obj.addNote(sprintf('correlation->background= %f', obj.corr_b));
            end
            
            % check for correlation with area
            [obj.corr_a, obj.covar_a] = obj.correlation(obj.areas_at_star);
            if obj.corr_a>obj.max_corr
                obj.keep = 0;
                obj.is_corr_area = 1;
                obj.addNote(sprintf('correlation->area= %f', obj.corr_a));
            end
            
            % check for correlation with PSF width
            [obj.corr_w, obj.covar_w] = obj.correlation(obj.widths_at_star);
            if abs(obj.corr_w)>obj.max_corr
                obj.keep = 0;
                obj.is_corr_width = 1;
                obj.addNote(sprintf('correlation->width= %f', obj.corr_w));
            end
            
            % add any other checks you can think about
            % ...
            
        end
        
        function [corr, covar] = correlation(obj, vector, flux)
            
            if nargin<3 || isempty(flux)
                flux = obj.flux_detrended;
            end
            
            if ischar(flux)
                
                if util.text.cs(flux, 'detrended')
                    flux = obj.flux_detrended;
                elseif util.text.cs(flux, 'raw')
                    flux = obj.flux_raw_all(:,obj.star_index); 
                else
                    error('Unknown flux option "%s". Try "detrended" or "raw"', flux); 
                end

            end
            
            f = flux;
            t = obj.time_indices;
            
            f_temp = f;
            f_temp(t) = NaN;
            Mf = nanmean(f_temp); % get the mean outside the event
            Sf = nanstd(f_temp); % get the std outside the event
            f = (f - Mf)./Sf;
            
            v_temp = vector;
            v_temp(t) = NaN;
%             Mv = nanmean(v_temp); % get the mean outside the event
            Mv = nanmean(vector(t)); % get the mean inside the event!
            Sv = nanstd(v_temp); % get the std outside the event!
            vector = (vector - Mv)./Sv;
            
            denom = sqrt(nansum(f(t).^2).*nansum(vector(t).^2));
            if denom==0
                corr = 0;
            else
                corr = nansum(f(t).*vector(t))./denom; % ignore the STD normalization and just get the correlation coefficient (the normalized projection)
            end
            
            covar = nansum(f(t).*vector(t)); % get the covariance between the two vectors in the event region, normalized by the STD outside the event region
            
        end
        
    end
        
    methods % utilities
        
        function addNote(obj, str)
            
            if isempty(obj.notes)
                obj.notes = str;
            else
                obj.notes = [obj.notes ', ' str];
            end
            
        end
        
        function new_obj = reduce_memory_copy(obj)
            
            if util.text.cs(obj.which_batch, 'first')
                other_batch = 'second';
            else
                other_batch = 'first';
            end
            
            % make shallow, temporary copies of arrays
            stack_temp = obj.(['stack_' obj.which_batch]);
            cutouts_temp = obj.(['cutouts_' obj.which_batch]);
            pos_temp = obj.(['positions_', obj.which_batch]);
            
            stack_temp2 = obj.(['stack_' other_batch]);
            cutouts_temp2 = obj.(['cutouts_' other_batch]);
            pos_temp2 = obj.(['positions_', other_batch]);
            
            obj.clearImages; % first clear the images to make copying fast
            new_obj = util.oop.full_copy(obj);
            
            % new object gets only the relevant arrays (50% storage) and
            % only gets shallow copies of the data (quick copy)
            new_obj.(['stack_' obj.which_batch]) = stack_temp;
            new_obj.(['cutouts_' obj.which_batch]) = cutouts_temp;
            new_obj.(['positions_', obj.which_batch]) = pos_temp;
            
            % recover the temporary arrays so input event returns to normal
            obj.(['stack_' obj.which_batch]) = stack_temp;
            obj.(['cutouts_' obj.which_batch]) = cutouts_temp;
            obj.(['positions_', obj.which_batch]) = pos_temp;
            
            obj.(['stack_' other_batch]) = stack_temp2;
            obj.(['cutouts_' other_batch]) = cutouts_temp2;
            obj.(['positions_', other_batch]) = pos_temp2;
            
        end
        
        function reload_memory(obj) % loads cutouts/positions/stack from file
            
            obj.cutouts_first = h5read(obj.filename_first, '/cutouts'); 
            obj.positions_first = h5read(obj.filename_first, '/positions'); 
            obj.stack_first = h5read(obj.filename_first, '/stack'); 

            obj.cutouts_second = h5read(obj.filename_second, '/cutouts'); 
            obj.positions_second = h5read(obj.filename_second, '/positions'); 
            obj.stack_second = h5read(obj.filename_second, '/stack'); 
            
        end
        
        function val = getRunName(obj)
            
            val = ''; 

            if ~isempty(obj.head) && ~isempty(obj.head.OBJECT)
                val = obj.head.OBJECT;
            elseif ~isempty(obj.filename)
                [~, val] = fileparts(fileparts(obj.filename)); 
            end

        end
        
        function save(obj, filename)
            
            if nargin<2 || isempty(filename)
                
                filename = 'event'; 

                run_name = obj.getRunName;

                if ~isempty(run_name)
                    filename = [filename '_' run_name];
                end

            end
            
            event = obj;
            
            save(filename, 'event'); 
            
        end
        
        function saveDialog(obj, filename)
            
            if nargin<2 || isempty(filename)
                
                filename = 'event'; 

                run_name = ''; 

                run_name = obj.getRunName;

                if ~isempty(run_name)
                    filename = [filename '_' run_name];
                end

            end
            
            [filepath,filename,ext] = fileparts(filename);
            
            if ~isempty(ext)
                filename = [filename ext];
            else
                filename = [filename '.mat']; 
            end
            
            if isempty(filepath)
                filepath = fullfile(getenv('DATA'), 'WFAST/saved/events/'); 
            end
            
            [filename, filepath] = uiputfile(fullfile(filepath, filename)); 
            
            if ~isequal(filename, 0)
                obj.save(fullfile(filepath, filename)); 
            end
            
        end
        
        function str_out = print_gaia_data(obj)
            
            if ~isempty(obj.gaia_data)
                Mag = obj.gaia_data.Mag_BP;
                Teff = obj.gaia_data.Teff;
            else
                Mag = [];
                Teff = [];
            end
            
            if ~isempty(obj.head)
                RA = obj.head.RA;
                Dec = obj.head.DEC;
            else
                RA = '';
                Dec = '';
            end

            str = sprintf('run: %s | RA= %s | Dec= %s | Mag BP= %4.2f | Teff= % 5d', obj.getRunName, RA, Dec, Mag, round(Teff)); 
                        
            if nargout>0
                str_out = str;
            else
                disp(str);
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('parent', []);
            input.input_var('font_size', 18);
            input.scan_vars(varargin{:});
            
            if isempty(input.parent)
                input.parent = gcf;
            end
            
            delete(input.parent.Children);
            
            ax1 = axes('Parent', input.parent, 'Position', [0.1 0.25 0.35 0.5]);
            obj.showFlux('ax', ax1);
            
            ax2 = axes('Parent', input.parent, 'Position', [0.48 0.3 0.20 0.4]);
            obj.showCutouts('ax', ax2);
            
            ax3 = axes('Parent', input.parent, 'Position', [0.7 0.3 0.25 0.4]);
            obj.showStack('ax', ax3);
            
        end
        
        function showFlux(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('font_size', 12);
            input.scan_vars(varargin{:});
            
            if isempty(input.ax), input.ax = gca; end
            
            input.ax.NextPlot = 'replace';
            h1 = plot(input.ax, obj.timestamps, obj.flux_filtered, 'LineWidth', 2);
            h1.DisplayName = 'Filtered LC';
            
            input.ax.NextPlot = 'add';
            range = obj.time_indices;
            if ~isempty(range)
                h2 = plot(input.ax, obj.timestamps(range), obj.flux_filtered(range), 'LineWidth', 2);
                h2.DisplayName = 'trigger region';
            end
            
            f = obj.flux_detrended;
            s = obj.std_flux;
            
            input.ax.ColorOrderIndex = 3;
            h3 = plot(input.ax, obj.timestamps, f./s, '-');
            h3.DisplayName = 'Raw LC';
            
            sign = 1;
            if obj.is_positive==0
                sign = -1;
            end

            h4 = plot(input.ax, obj.kernel_timestamps+obj.peak_timestamp, obj.best_kernel*5*sign, ':');
            h4.DisplayName = 'best filter';
            
            h5 = plot(input.ax, obj.timestamps, obj.backgrounds_at_star, '--'); 
            h5.DisplayName = 'background';
            
            h6 = plot(input.ax, obj.timestamps, obj.widths_at_star, 'p', 'MarkerSize', 3);
            h6.DisplayName = 'PSF width';
            
            h6 = plot(input.ax, obj.timestamps, obj.offsets_x_at_star, 'x', 'MarkerSize', 4);
            h6.DisplayName = 'offset x';
            
            h8 = plot(input.ax, obj.timestamps, obj.offsets_y_at_star, '+', 'MarkerSize', 4);
            h8.DisplayName = 'offset y';
            
%             h8 = bar(input.ax, obj.timestamps, obj.bad_pixels_at_star-mean(obj.bad_pixels_at_star)-5, 'BaseValue', -5, 'FaceAlpha', 0.5);
%             h8.DisplayName = 'relative bad pixels';
            
            xlabel(input.ax, 'timestamp (seconds)');
            ylabel(input.ax, 'flux S/N');
            
%             if strcmp(obj.which_batch, 'first')
%                 lh = legend(input.ax, 'Location', 'SouthEast');
%             else
%                 lh = legend(input.ax, 'Location', 'SouthWest');
%             end

            lh = legend(input.ax, 'Location', 'NorthOutside', 'Orientation', 'Vertical');
            
            lh.FontSize = input.font_size;
            lh.NumColumns = 3;

            title(input.ax, obj.notes, 'FontSize', input.font_size);
            
            util.plot.inner_title(sprintf('id: %d | star: %d | batches: %d-%d | event S/N= %4.2f | star S/N= %4.2f', ...
                obj.serial, obj.star_index, obj.batch_index_first, obj.batch_index_second, obj.snr, obj.star_snr),...
                'ax', input.ax, 'Position', 'NorthWest', 'FontSize', input.font_size);
            
            mx = max(max(abs(obj.flux_filtered)), max(abs(f./s))); 
            
            input.ax.YLim = [-1 1].*1.2.*mx;
            
%             input.ax.YLim(2) = max(abs(input.ax.YLim));
%             input.ax.YLim(1) = -max(abs(input.ax.YLim))-2;
            
            input.ax.XLim = [obj.timestamps(1) obj.timestamps(end)];

            input.ax.NextPlot = 'replace';
            
        end
        
        function showCutouts(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('parent', []); 
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('position', []);
            input.input_var('number', 9);
            input.input_var('bias', []);
            input.input_var('dynamic_range', []);
            input.input_var('font_size', 10);
            input.scan_vars(varargin{:});
            
            if ndims(obj.cutouts)>=4
                cutouts = cat(3, obj.cutouts_first(:,:,:,obj.star_index), obj.cutouts_second(:,:,:,obj.star_index));
            else
                cutouts = cat(3, obj.cutouts_first, obj.cutouts_second);
            end
            
            if ~isempty(cutouts)
                
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
                    
                end
                
                idx = obj.time_index; 
                
                idx_start = idx - floor(input.number/2); 
                idx_end = idx + ceil(input.number/2) - 1;
                
                if idx_start<1
                    idx_end = idx_end + 1 -idx_start;
                    idx_start = 1;
                end
                
                if idx_end>size(cutouts,3)
                    idx_start = idx_start + size(cutouts,3) - idx_end;
                    idx_end = size(cutouts,3);
                end
                
                rad = [];
                if ~isempty(obj.offsets_x_at_star) && ~isempty(obj.offsets_x_at_star)

                    cen = floor([size(cutouts,2), size(cutouts,1)]/2)+1;
                    cen = cen + [obj.offsets_x_at_star(idx_start:idx_end) obj.offsets_y_at_star(idx_start:idx_end)];
                    
%                     if ~isempty(obj.gauss_sigma)
%                         rad = obj.gauss_sigma;
%                         str = sprintf('sigma= %4.2f', rad(end)); 
%                         col = 'magenta';
%                     elseif ~isempty(obj.aperture)
                    if ~isempty(obj.aperture)
                        rad = obj.aperture;
                        str = sprintf('ap= %4.2f', rad(end)); 
                        col = 'green';
                    end

                end

                Nrows = ceil(sqrt(input.number));
                Ncols = Nrows;

                if isempty(input.bias) || isempty(input.dynamic_range)
                    dyn = util.img.autodyn(cutouts);
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

                    imagesc(ax{ii}, cutouts(:,:,idx_start+ii-1));
                    ax{ii}.CLim = dyn;
                    axis(ax{ii}, 'image');
                    ax{ii}.XTick = [];
                    ax{ii}.YTick = [];
                    
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
            
        end
        
        function showStack(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('range', 100);
            input.scan_vars(varargin{:});
            
            if isempty(input.ax), input.ax = gca; end
            
            if ~isempty(obj.stack)
                
                util.plot.show(obj.stack, 'fancy', 1, 'autodyn', 1, 'ax', input.ax); 
                title(input.ax, '');
                
                if ~isempty(obj.cutouts) && ~isempty(obj.positions)
                    
                    S = size(obj.cutouts);
                    S = S(1:2); 
                    
                    P = obj.positions(obj.star_index, :); 
                    
                    C = P - floor(S/2) - 0.5; % corner of the rectangle
                    
                    rectangle('Position', [C S], 'Parent', input.ax); 
                    
                    input.ax.XLim = [P(1)-input.range P(1)+input.range];
                    input.ax.YLim = [P(2)-input.range P(2)+input.range];
                    
                    if input.ax.XLim(1)<1
                        input.ax.XLim(1) = 1;
                    elseif input.ax.XLim(2)>size(obj.stack,2)
                        input.ax.XLim(2) = size(obj.stack,2);
                    end
                    
                    if input.ax.YLim(1)<1
                        input.ax.YLim(1) = 1;
                    elseif input.ax.YLim(2)>size(obj.stack,1)
                        input.ax.YLim(2) = size(obj.stack,1);
                    end
                    
                end
                
            end
            
        end
        
    end    
    
end

