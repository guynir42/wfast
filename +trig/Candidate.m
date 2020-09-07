classdef Candidate < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        head; % header information
        kern_props; % struct with the data about the best kernel
        star_props; % table with one row from astrometry/GAIA
        sim_pars = struct; % struct with the simulation parameters (for simulated events only)

        finder_pars;
        store_pars;
        checker_pars; 
        
    end
    
    properties % inputs/outputs
        
        serial; % each candidates's id in the run
        
        analysis_time = ''; % timestring when this event was identified by analysis (UTC)
        
        notes = {}; % anything added to the event, like failed cuts
        
        cutouts; % cutouts of the selected star, for the duration of the extended region
        stack; % stack of the images in which this event is detected (optional)
        kernel; % the lightcurve of the matched-filter kernel used 
        
        time_index; % out of the extended region
        kern_index; % from the full filter bank
        star_index; % from the stars that passed the initial burn-in (not from the subset that survived the pre-filter
        star_index_global; % from all the stars in this run (compare this to the catalog, cutouts, etc)
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
        flux_corrected; % fluxes for all stars, corrected by PSD or by removing linear fit
        flux_filtered; % flux after matched-filtering, for the peak star and kernel
        psd; % power spectral density for this star only
        freq_psd; % frequency values for the PSD
        
        auxiliary; % values of auxiliary measurements (e.g., background, offsets) for all stars and all frames in the extended region
        aux_names; % cell array with the names of each auxiliary
        aux_indices; % struct with a field for each auxiliary name, containing the index in the big matrix
        
        relative_dx; % the offsets_x for this star, minus the flux weighted mean offsets_x of all stars
        relative_dy; % the offsets_x for this star, minus the flux weighted mean offsets_x of all stars
        
        cut_matrix; % the different cuts from the QualityChecker (for this star, on the extended region)
        cut_names; % a cell array with the names of the cuts, e.g., 'shakes', 'offset_size', 'corr_x_25',...
        cut_indices; % a struct with each cut name as a field containing the index of that cut in the matrix's 2nd dim
        cut_thresh; % which threshold is used for each cut (a row vector)
        cut_two_sided; % true when the cut value should be compared in abs() to the threshold (a row vector)
        
        cut_string = {}; % cell string discribing the cuts that this event hit
        cut_hits = []; % which cuts coincide with the event peak
        cut_value = []; % what is the value that triggered the cut (maximum within range of the peak)
        
        batch_number;
        
        snr;
        
    end
    
    properties % switches/controls
        
        kept = 1; % by default any new event is considered real (kept)
        
        classification = ''; 
        
        display_cuts = []; % user can add more cuts to display even if the event did not happen to cross that cut
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
        
        use_show_secrets = true; % this should be turned off at some point
        
        run_name; 
        run_date;
        
        is_positive; 
        
        is_simulated = 0; % by default events are not simulated

        flux_buffer; % flux buffer going back as far as possible, for this star only
        timestamp_buffer; % timestamps for the above
        filtered_flux_past_values; % normalization data for this star and this kernel only. Can be either raw filtered_flux values for the last "backround" period, or the variance values from the var_buf
        
        flux_raw_all; % the raw flux (over the extended region) for all stars
        flux_corrected_all; % the corrected flux (over the exteneded region) for all stars
        auxiliary_all; % the auxiliary data (for the extended region) for all stars. 
        
        search_start_idx;
        search_end_idx;
        
        kern_extra; % any other kernels that passed the lower threshold for kernels
        star_extra; % any other stars that passed the lower threshold for stars
        
        version = 1.00;
        
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
    
    methods % reset/clear
        
    end
    
    methods % getters
        
        function classes = getListOfClasses(obj)
            
            classes = {'cosmic ray', 'satellite', 'flare', 'occultation', 'mystery', 'edge effect', 'artefact'};
            
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
        
        function val = star_snr(obj)
            
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
            
            str = sprintf('id: %d | star: %d | frame= %d | time: %d-%ds | event S/N= %4.2f | star S/N= %4.2f', ...
                    obj.serial, obj.star_index, obj.time_index, round(obj.time_start), round(obj.time_end), obj.snr, obj.star_snr);
            
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
                str = sprintf('R= %4.2f | r= %4.2f | b= %4.2f | v= %4.1f | invert= %d', ...
                    obj.kern_props.R, obj.kern_props.r, obj.kern_props.b, obj.kern_props.v, ~obj.is_positive);
            end 
            
        end
        
        function str = printStellarProps(obj)
            
            if isempty(obj.star_props)
                str = '';
            else
                str = sprintf('mag= %4.2f | temp= %dK', obj.star_props.Mag_BP, round(obj.star_props.Teff)); 
            end
            
        end
        
        function val = getStellarSize(obj)
            
            if ~isempty(obj.star_props) && ismember('Mag_BP', obj.star_props.Properties.VariableNames) && ~isnan(obj.star_props.Mag_BP) && ...
                    ismember('Teff', obj.star_props.Properties.VariableNames) && ~isnan(obj.star_props.Teff) % catalog entry exists and has non NaN mag/temperature
                
                if ~ismember('BolMag', obj.star_props.Properties.VariableNames) || isnan(obj.star_props.BolMag) || ...
                    ~ismember('BolTemp', obj.star_props.Properties.VariableNames) || ~isnan(obj.star_props.BolTemp) % no bolometric data, need to calculate it now
                
                    [C,T] = util.ast.bol_corr(obj.star_props.Teff, obj.star_props.Mag_BP, obj.star_props.Mag_RP, 'BP', 'RP', 'GAIA', 0);
                
                    bol_mag = obj.star_props.Mag_BP + C; 
                
                    obj.star_props.BolMag = bol_mag;
                    obj.star_props.BolTemp = T; 
                end
                
                val = util.ast.stellar_size(obj.star_props.BolMag, obj.star_props.BolTemp, 'units', 'FSU', 'dist', 40/206000); % distance is given in parsec (40AU)
                
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
            
            str = sprintf('Run: "%s" | ecl lat= %4.1f deg \ncoords: %s %s | V= %d', obj.run_name, ECL, RA, Dec, round(v)); 
            
        end
        
        function str = printSimulationData(obj)
            
            if ~obj.is_simulated
                str = '';
            else
                str = sprintf('SIM: R= %4.2f | r= %4.2f | b= %4.2f | v= %4.1f', ...
                    obj.sim_pars.R, obj.sim_pars.r, obj.sim_pars.b, obj.sim_pars.v); 
            end
            
        end
        
    end
    
    methods % setters
        
        function addNote(obj, str)
            
            if isempty(obj.notes)
                obj.notes{1} = str;
            else
                obj.notes{end+1} = str;
            end
            
        end
        
        function setRunNameDate(obj)
            
            path = fileparts(obj.getFilename); 
            
            if ~isempty(path)
                [path, obj.run_name] = fileparts(path); 
                [path, obj.run_date] = fileparts(path); 
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
        
        function val = correlateSlopes(obj, timescale)
            
            if nargin<2 || isempty(timescale)
                timescale = 10;
            end
            
            t = obj.timestamps;
            x = obj.auxiliary(:,obj.aux_indices.offsets_x);
            y = obj.auxiliary(:,obj.aux_indices.offsets_y);
            
            cx = util.series.correlation(x,t, timescale).*sqrt(timescale);
            cy = util.series.correlation(y,t, timescale).*sqrt(timescale);
            
            val = (abs(cx)+abs(cy));
            
        end
        
        function val = getSlopes(obj, timescale)
            
            if nargin<2 || isempty(timescale)
                timescale = 10;
            end
            
            t = obj.timestamps;
            x = obj.auxiliary(:,obj.aux_indices.offsets_x);
            y = obj.auxiliary(:,obj.aux_indices.offsets_y);
            
            x = x-nanmedian(x); 
            y = y-nanmedian(y); 
            
            k = (-timescale/2:timescale/2)'; 
%             k = (1:timescale)'; 
            k = k./sqrt(sum(k.^2)); 
            
            val = [x,y];
            val = abs(filter2(k,val)); 
%             val = sqrt(sum(val.^2,2)); 
            val = [x,y,val]; 
                        
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj_vec, varargin)
            
            input = util.text.InputVars;
            input.input_var('index', []); 
            input.input_var('kept', [], 'show_kept'); 
            input.input_var('cuts', []); 
            input.input_var('parent', []);
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
            
            if isempty(input.parent.UserData)
                input.parent.UserData = struct('number', length(obj_vec), 'index', input.index, ...
                    'show_kept', input.kept, 'show_cuts', input.cuts); % add other state variables here
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
            
            obj = obj_vec(idx); 
            
            delete(input.parent.Children);
            
            margin_left = 0.07;
            
            ax1 = axes('Parent', input.parent, 'Position', [margin_left 0.55 0.55 0.35]);
            obj.showRawFlux('ax', ax1, 'on_top', 1, 'equal', 0, 'title', 0);
            
            ax2 = axes('Parent', input.parent, 'Position', [margin_left 0.2 0.55 0.35]);
            obj.showFilteredFlux('ax', ax2, 'title', 0, 'cuts', input.cuts); 
            
            ax3 = axes('Parent', input.parent, 'Position', [0.68 0.2 0.3 0.5]);
            obj.showCutouts('ax', ax3);
            
            %%%%%%%%%%%%%%%%%%%%%% panel info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            info_panel = uipanel(input.parent, 'Units', 'Normalized', 'Position', [margin_left 0.9 1-margin_left 0.1], 'title', 'info'); 
            
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
                'Units', 'Normalized', 'Position', [0.05 0 0.9 1], 'FontSize', 14, 'HorizontalAlignment', 'Left'); 
            
            %%%%%%%%%%%%%%%%%%%%%% panel controls %%%%%%%%%%%%%%%%%%%%%%%%%
            
            control_panel = uipanel(input.parent, 'Units', 'Normalized', 'Position', [margin_left 0.0 1-margin_left*2 0.1], 'title', 'controls'); 
            
            button = uicontrol(control_panel, 'Style', 'pushbutton', 'String', '-', 'FontSize', 20, ...
                'Units', 'Normalized', 'Position', [0.0 0.1 0.1 0.9], ...
                'Callback', @obj_vec.callback_prev, 'UserData', input.parent, ...
                'Tooltip', 'go back to previous candidate'); 
            
            button = uicontrol(control_panel, 'Style', 'edit', 'String', sprintf('index= %d / %d', input.parent.UserData.index, length(obj_vec)), 'FontSize', 16, ...
                'Units', 'Normalized', 'Position', [0.1 0.1 0.2 0.9], ...
                'Callback', @obj_vec.callback_index, 'UserData', input.parent, ...
                'Tooltip', 'jump to any candidate index'); 
            
            button = uicontrol(control_panel, 'Style', 'pushbutton', 'String', '+', 'FontSize', 20, ...
                'Units', 'Normalized', 'Position', [0.3 0.1 0.1 0.9], ...
                'Callback', @obj_vec.callback_next, 'UserData', input.parent, ...
                'Tooltip', 'go forward to next candidate'); 
            
            if input.parent.UserData.show_kept
                kept_string = 'kept';
            else
                kept_string = 'all'; 
            end
            
            button = uicontrol(control_panel, 'Style', 'pushbutton', 'String', kept_string, 'FontSize', 20, ...
                'Units', 'Normalized', 'Position', [0.43 0.1 0.1 0.9], ...
                'Callback', @obj_vec.callback_kept, 'UserData', input.parent); 
            
            if input.parent.UserData.show_kept
                button.TooltipString = 'Showing only the kept candidates';
            else
                button.TooltipString = 'Showing all candidates';
            end
            
            if input.cuts
                str = 'with cuts';
            else
                str = 'no cuts';
            end
            
            button = uicontrol(control_panel, 'Style', 'pushbutton', 'String', str, ...
                'Units', 'Normalized', 'Position', [0.55 0.1 0.1 0.9], 'FontSize', 12, ...
                'Callback', @obj_vec.callback_cuts, 'UserData', input.parent);
            
            if input.cuts
                button.TooltipString = 'show cuts on top of the filtered flux';
            else
                button.TooltipString = 'show only the filtered flux without cuts'; 
            end
            
            button = uicontrol(control_panel, 'Style', 'edit', 'String', ['Cuts= ' util.text.print_vec([obj.cut_hits, obj.display_cuts], ', ')], 'FontSize', 18, ...
                'Units', 'Normalized', 'Position', [0.67 0.1 0.2 0.9], ...
                'Callback', @obj_vec.callback_parse_cuts, 'UserData', input.parent);
            
            button.TooltipString = sprintf('Showing the cuts for: %s', strjoin(obj.cut_names([obj.cut_hits, obj.display_cuts]), ', '));  
            
            button = uicontrol(control_panel, 'Style', 'pushbutton', 'String', 'Classify', 'FontSize', 18, ...
                'Units', 'Normalized', 'Position', [0.89 0.1 0.1 0.9], ...
                'Callback', @obj_vec.callback_classify, 'UserData', input.parent); 
            
            button.TooltipString = 'choose a classification for this candidate'; 
            
%             ax4 = axes('Parent', input.parent, 'Position', [0.7 0.3 0.25 0.4]);
%             obj.showStack('ax', ax4);

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
            
            h1 = plot(input.ax, x, obj.flux_raw, '-', 'LineWidth', 2, 'Color', [0.929 0.694 0.125]);
            h1.DisplayName = 'raw flux';
            
            input.ax.NextPlot = 'add';
            h2 = plot(input.ax, xk, obj.flux_mean*obj.kernel_lightcurve(~obj.is_positive), ':', 'LineWidth', 2, 'Color', input.ax.Colormap(4,:));
            h2.DisplayName = 'best kernel';
            
            input.ax.NextPlot = 'replace';
            
            input.ax.FontSize = input.font_size;
            input.ax.YAxis(1).Color = [0 0 0];
            
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
            
            input.ax.YLim = [mn-0.25.*abs(mn) mx+0.25.*abs(mx)]; 
            
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
            
            str = 'Classify this candidate:';
            
            classes = obj.getListOfClasses;
            
            rep = listdlg('ListString', classes, 'PromptString', str, 'ListSize',[150,250],...
                'InitialValue',length(classes), 'SelectionMode','single', 'Name','classification');
            
            if ~isempty(rep)
                obj.classification = classes{rep}; 
            end
            
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
            
            obj.show('index', idx); 
            
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
            
            obj.show('index', idx);
            
        end
        
        function callback_next(obj, hndl, ~)
            
            num = hndl.UserData.UserData.number;
            idx = hndl.UserData.UserData.index;
            use_kept = hndl.UserData.UserData.show_kept;
            
            if use_kept
                N = num+1;
            else
                N = 1;
            end
            
            for ii = 1:N
                
                idx = idx + 1; 
                if idx>num
                    idx = 1;
                end
                
                if use_kept && obj(idx).kept % when we stumble on a kept event
                    break;
                end
                
            end
            
            obj.show('index', idx); 
            
        end
        
        function callback_kept(obj, hndl, ~)
            
            use_kept = ~hndl.UserData.UserData.show_kept;
            
            hndl.UserData.UserData.show_kept = use_kept;
            
            hndl.Value = use_kept;
            
            if use_kept
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
            
            obj.show;
            
        end
        
        function callback_classify(obj, hndl, ~)
            
            idx = hndl.UserData.UserData.index;
            
            obj(idx).popupClassifier;
            
            obj.show;
            
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
            
            obj.show; 
            
        end
        
    end
    
end

