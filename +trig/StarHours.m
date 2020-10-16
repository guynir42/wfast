classdef StarHours < handle
% Container for storing useful (and also disqualified) star-time at different
% values of S/N and for each star separately. 
% 
% The only way to use this object is to give it a QualityChecker object as
% the first argument to the input() function. The second optional input to  
% this function, bad_batch=0, can be set to true to tell StarHours that the
% incoming batch is completely disqualified. This can happen if e.g., there
% are too many event candidates in a single batch (and it is black-listed). 
% 
% This object does not have many user-defined parameters. They are:
% -snr_bin_min: minimal S/N values to store in the histogram. Default 5. 
% -snr_bin_max: maximal S/N values to store in the histogram. Default 40. 
% -snr_bin_width: the bin size of the histogram. Default 0.5. 
%
% Once these are defined, the object simply sorts the star-time it receives 
% from the QualityChecker for each batch, placing the amount of time (in 
% seconds) into the right bin based on the star index and S/N value at that 
% time. That means each star can contribute, during the course of a run, to
% several different S/N bins. 
% The reason to keep each star individually is twofold:
% (a) we can debug certain runs, finding out where star-time was lost. 
% (b) if we black list any stars at the end of the run, we must remove
%     their contribution to the star-time histograms. 
% After the run is done and stars are removed, we can safely sum the data
% into 1D histograms for different S/N bins only. This is done when these
% histograms are transferred to the RunSummary object. 
% 
% The outputs from this object are:
% -runtime: the total time, in seconds, that was processed (not including 
%           the burn-in period. 
% -histogram: a 2D histogram of the useful star-time for each S/N and each 
%             star index. This is the main data product. 
% -histogram_with_losses: the same as histogram, but without disqualifying 
%                         any time for any reason. This represents the total
%                         time that was scanned. 
% -losses_exclusive: the amount of time, in each S/N bin and for each star
%                    index, for each data quality cut type. This is a 3D 
%                    matrix where the last index can be matched to each cut
%                    type using "cut_indices" (see example in QualityChecker). 
%                    The "exclusive" losses mean times that were disqualified
%                    only by a single cut. This helps answer "how much time
%                    would we gain if we did not place this cut". 
% -losses_inclusive: The same as the above, only including time lost to each 
%                    cut regardless of whether other cuts also affected that
%                    star at that time. It answers the question: "how much 
%                    time would we lose if we only place this cut."
% -losses_bad_stars: The amount of time lost because stars were black listed
%                    and removed from the run. This only happens if running
%                    event-finding and too many candidates show up in one star. 
% -snr_bin_edges: the edges of the bins are stored for plotting, etc. 
% -star_bin_edges: same as above. This will just be 1:number of stars. 
% 
% NOTE: all histograms are saved in SECONDS, so star-hours can be recovered
%       by dividing the values by 3600 (as done in the plotting tools). 
%
% Plotting: use the showHistogram() function to display the results. 
% Some optional parameters that can be used with this function are:
% -type: specificy which cut you want to view losses for, as a number or as 
%        a name of the cut. Default is [], which shows just the good star
%        hours instead of the losses. 
% -exclusive: when true, show exclusive losses. Only works when specifying 
%             the "type" option. Default is false. 
% -min_snr and max_snr: put limits on what to plot in the S/N range. The 
%                       default is full bin range (snr_bin_min to snr_bin_max). 
% -sum: if false, will show the values for each star individually (in 2D). 
%       if true, will show the summed result for all stars. Default false. 
% -hours: translate the star-seconds to star-hours. Default true. 
% -axes: which axes to plot into. Default is gca(). 
% -font_size: for the axes and labels. Default is 18. 
% -log: show the results in log-scale. Default is true. 
% -color: give a specific color to the bars (for multiple plots). Default 
%         is [], which just picks the figure's next color. 
% -alpha: the transparency of the bars. Default is 0.7. 
% 
% Another diagnostic tool is the printReport() function. 
% It prints a summary of the star hours and the losses to each cut. 
% Optional arguments:
% -units: choose "seconds" or "hours" (default) or "percent" (or "%"). 
%         This shows the losses in different units, or as a fraction of the 
%         total star-time. The total good times are always shown as hours. 
% -format: choose "text" (default) or "latex". The default just prints a 
%          table in regular ascii. The "latex" option formats the table as 
%          a latex table that can be copied into a tabular environment. 


    properties(Transient=true)
        
    end
    
    properties % objects
        
        head@head.Header;
        
    end
    
    properties % inputs/outputs
        
        snr_bin_edges; % this is calculated from the parameters below
        star_bin_edges; % this is calculated after we know how many stars we've got
        
        runtime = 0; % total runtime since last reset()
        histogram; % the number of useful seconds accumulated for each star and each S/N bin (this is saved after subtracting losses)
        histogram_with_losses; % the number of seconds accumulated without excluding anything
        losses_exclusive; % same as histogram, only counting the losses due to each cut (exclusive means times where ONLY this cut was responsible)
        losses_inclusive; % same as histogram, only counting the losses due to each cut (inclusive means times where this cut was ALSO responsible)
        losses_bad_stars; % number of seconds lost when disqualifying stars from the black list 
        
        cut_names = {}; % get these from the QualityChecker
        cut_indices = {}; % get these from the QualityChecker
        
        % these are temporary values used each batch:
        timestamps; % timestamps for the extended region
        flux; % flux for the extended region
        idx_start; % start of the search region
        idx_end; % end of the search region
        star_snr; % get the S/N for each star based on the extended region
        time_step; % average dt (0.04s)
        batch_length; % how many seconds are in the batch, total
        flags; % a matrix with all the regions flagged for each cut
        
    end
    
    properties % switches/controls
        
        % the following properties are used to build the S/N bin edges
        snr_bin_min = 2; % stars with lower S/N are not even tested for events (here S/N is calculated on the raw fluxes)
        snr_bin_max = 40; % biggest value we expect to see in the star-hour histogram
        snr_bin_width = 0.5; % for the star-hour histogram
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        star_snr_dist; % output of histcounts2(). Should be 1 in the right star and S/N bin, and zeros elsewhere. 
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = StarHours(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.StarHours')
                if obj.debug_bit>1, fprintf('StarHours copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('StarHours constructor v%4.2f\n', obj.version); end
            
                obj.reset;
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.runtime = 0;
            obj.histogram = [];
            obj.histogram_with_losses = [];
            obj.losses_exclusive = [];
            obj.losses_inclusive = [];
            obj.losses_bad_stars = [];
            obj.snr_bin_edges = [];
            obj.star_bin_edges = []; 
            
            obj.clear; 
            
        end
        
        function clear(obj)
            
            obj.timestamps = [];
            obj.flux = [];
            obj.idx_start = [];
            obj.idx_end = [];
            obj.star_snr = [];
            obj.time_step = []; 
            obj.batch_length = [];
            obj.flags = []; 
            obj.star_snr_dist = [];
            
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function input(obj, checker, bad_batch)
            
            if nargin<3 || isempty(bad_batch)
                bad_batch = 0;
            end
            
            obj.clear;
            
            % get the data from checker
            obj.timestamps = checker.extended_timestamps;
            obj.flux = checker.extended_flux - ...
                checker.extended_aux(:,:,checker.aux_indices.backgrounds).*checker.extended_aux(:,:,checker.aux_indices.areas);
            obj.idx_start = checker.search_start_idx; 
            obj.idx_end = checker.search_end_idx; 
        
            obj.flags = checker.cut_flag_matrix(obj.idx_start:obj.idx_end,:,:); % only look at the search region! 
            
            obj.cut_names = checker.cut_names; 
            obj.cut_indices = checker.cut_indices;
        
            % calculate some stuff based on that
            obj.star_snr = nanmean(obj.flux,1)./nanstd(checker.extended_flux); % S/N on the extended region
            obj.time_step = median(diff(obj.timestamps)); 
            obj.batch_length = (obj.idx_end - obj.idx_start + 1).*obj.time_step; 
            
            if isempty(obj.star_bin_edges) || isempty(obj.snr_bin_edges) || isempty(obj.histogram)
                obj.makeHistograms;
            end
            
            % fill up the cumulative stuff for the whole run:
            obj.runtime = obj.runtime + obj.batch_length; 
            
            obj.star_snr_dist = histcounts2(obj.star_snr, obj.star_bin_edges(1:end-1), obj.snr_bin_edges, obj.star_bin_edges-0.5); % resulting histogram is stars on the X axis and S/N on the Y axis
            
            if isempty(obj.histogram_with_losses)
                obj.histogram_with_losses = obj.star_snr_dist.*obj.batch_length; % add the amount of time for each star in the right bin
            else
                obj.histogram_with_losses = obj.histogram_with_losses + obj.star_snr_dist.*obj.batch_length; % add the amount of time for each star in the right bin
            end
            if ~bad_batch % do not include anything in histogram or losses when the batch is bad
                
                bad_times = checker.bad_times; % this has "true" in every frame/star that is disqualified for anything
                lost_time = sum(bad_times(obj.idx_start:obj.idx_end,:)).*obj.time_step; 

                obj.histogram = obj.histogram + obj.star_snr_dist.*(obj.batch_length - lost_time); % add the amount of time for each star in the right bin

    %             flag_also = false(size(bad_times,1), size(bad_times,2), length(obj.cut_names)); % logical "true" in each frame/star where this cut is also or only included
    %             flag_only = false(size(bad_times,1), size(bad_times,2), length(obj.cut_names)); % logical "true" in each frame/star where this cut is only included (with no others!)

                flag_only = false(size(obj.flags)); % logical "true" in each frame/star where this cut is only included (with no others!)

                for ii = 1:length(obj.cut_names)

                    num_xors = sum(xor(obj.flags(:,:,ii), obj.flags),3); % each frame/star counts the number of places in other cuts that were not "true" at the same time as the current cut
                    flag_only(:,:,ii) = num_xors==length(obj.cut_names)-1; % for each frame/star if all other cuts have zero here, then the xor would be "true" for all other cuts except the current cut

                    obj.losses_inclusive(:,:,ii) = single(obj.losses_inclusive(:,:,ii) + obj.star_snr_dist.*sum(obj.flags(:,:,ii)).*obj.time_step); % add the amount of time (per star, per S/N) for this cut (even if this time is shared with other cuts)
                    obj.losses_exclusive(:,:,ii) = single(obj.losses_exclusive(:,:,ii) + obj.star_snr_dist.*sum(flag_only(:,:,ii)).*obj.time_step); % add the amount of time (per star, per S/N) for this cut (even if this time is shared with other cuts)

                end

            end
            
        end
        
        function makeHistograms(obj)
            
            obj.star_bin_edges = 1:size(obj.flux,2)+1;
            obj.snr_bin_edges = obj.snr_bin_min:obj.snr_bin_width:obj.snr_bin_max;
            
            obj.histogram = zeros(length(obj.snr_bin_edges)-1, length(obj.star_bin_edges)-1, 'single'); 
            obj.losses_exclusive = zeros(length(obj.snr_bin_edges)-1, length(obj.star_bin_edges)-1, length(obj.cut_names), 'single');
            obj.losses_inclusive = zeros(length(obj.snr_bin_edges)-1, length(obj.star_bin_edges)-1, length(obj.cut_names), 'single'); 
            obj.losses_bad_stars = zeros(length(obj.snr_bin_edges)-1, length(obj.star_bin_edges)-1, 'single'); 
            
        end
        
        function removeStars(obj, idx)
            
            obj.losses_bad_stars(:,idx) = obj.histogram(:,idx); % store the number of star seconds that were lost when we remove bad stars from the histograms
            
            % maybe replace NaN with zero? 
            obj.histogram(:,idx) = NaN;
            obj.losses_exclusive(:,idx,:) = NaN;
            obj.losses_inclusive(:,idx,:) = NaN;
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function viewer(obj, varargin)
            
            import util.stat.sum2; 
            
            parent = []; 
            for ii = 1:2:length(varargin)
                if util.text.cs(varargin{ii}, 'parent') && length(varargin)>ii
                    parent = varargin{ii+1};
                end
            end
            
            if isempty(parent) || ~isvalid(parent) % we did not get a legal figure/panel
                parent = gcf;
            end
            
            if isempty(parent.UserData) || ~isa(parent.UserData, 'util.text.InputVars')
                
                input = util.text.InputVars;
                input.input_var('good', true); % whether to show the good times
                input.input_var('types', []); % choose multiple cut types (use the numbers!)
                input.input_var('exclusive', false); % if false, show inclusive losses
                input.input_var('max_snr', obj.snr_bin_max); 
                input.input_var('rebin', []); % to be added later
                input.input_var('sum', true); % we still need to figure out how to display non-summed multiple histograms...
                input.input_var('hours', true);             
                input.input_var('log', true, 'use_log'); 
                input.input_var('alpha', 0.7); 
                input.input_var('parent', []); 
                
                parent.UserData = input; % save for later
                
            else
                input = parent.UserData; % grab an existing input object
            end
            
            input.scan_vars(varargin{:}); 
            input.parent = parent; 
            
            % done parsing inputs, all is saved in the parent's UserData. 
            
            delete(parent.Children); % get rid of existing panels/plots
            
            %%%%%%%%%%%%%%% histogram display %%%%%%%%%%%%%%%%%%%%
            
            hist_width = 0.8; 
            
            panel_image = uipanel(parent, 'title', '', ...
                'Units', 'Normalized', 'Position', [0 0.1 hist_width 0.8]); 
            
            ax = axes('Parent', panel_image); 
            
            if input.good
                obj.showHistogram('type', [], 'max_snr', input.max_snr, ...
                    'sum', input.sum, 'log', input.log, 'hours', input.hours, ...
                    'alpha', input.alpha); 
            end
            
            ax.NextPlot = 'add'; 
            
            types = unique(input.types); % sorted list, removing redundant entries
            
            for ii = 1:length(types)
                
                obj.showHistogram('type', types(ii), 'max_snr', input.max_snr, ...
                    'sum', input.sum, 'log', input.log, 'hours', input.hours, ...
                    'alpha', input.alpha, 'exclusive', input.exclusive); 
                
            end
            
            ax.NextPlot = 'replace'; 
            
            %%%%%%%%%%%%%%% display controls %%%%%%%%%%%%%%%%%%%%
            
            panel_display = uipanel(parent, 'title', 'display', ...
                'Units', 'Normalized', 'Position', [0 0 hist_width 0.1]); 
            
            margin = 0.03; 
            width = 0.17;
            bottom = 0.1; 
            height = 0.8; 
            
            color_on = [0.1 0.3 1]; 
            
            button.Position = [-margin/2 0 0 0]; % left edge
            
            button = uicontrol(panel_display, 'Style', 'pushbutton', 'String', 'log scale', ...
                'Units', 'Normalized', 'Position', [margin + button.Position(1)+button.Position(3), bottom, width, height], ...
                'Callback', @obj.callback_log, 'FontSize', 14, 'UserData', input); 
            
            if input.log, button.ForegroundColor = color_on; end
            
            button = uicontrol(panel_display, 'Style', 'pushbutton', 'String', 'sum', ...
                'Units', 'Normalized', 'Position', [margin + button.Position(1)+button.Position(3), bottom, width, height], ...
                'Callback', @obj.callback_sum, 'FontSize', 14, 'UserData', input); 
            
            if input.sum, button.ForegroundColor = color_on; end
            
            button = uicontrol(panel_display, 'Style', 'pushbutton', 'String', 'seconds', ...
                'Units', 'Normalized', 'Position', [margin + button.Position(1)+button.Position(3), bottom, width, height], ...
                'Callback', @obj.callback_hours, 'FontSize', 14, 'UserData', input); 
            
            if input.hours, button.String = 'hours'; end
            
            button = uicontrol(panel_display, 'Style', 'edit', 'String', sprintf('max S/N= %4.1f', input.max_snr), ...
                'Units', 'Normalized', 'Position', [margin + button.Position(1)+button.Position(3), bottom, width, height], ...
                'Callback', @obj.callback_snr, 'FontSize', 14, 'UserData', input); 
            
            button = uicontrol(panel_display, 'Style', 'edit', 'String', sprintf('alpha= %4.2f', input.alpha), ...
                'Units', 'Normalized', 'Position', [margin + button.Position(1)+button.Position(3), bottom, width, height], ...
                'Callback', @obj.callback_alpha, 'FontSize', 14, 'UserData', input); 
            
            %%%%%%%%%%%%%%% info display %%%%%%%%%%%%%%%%%%%%
            
            panel_info = uipanel(parent, 'title', 'info', ...
                'Units', 'Normalized', 'Position', [0 0.9 hist_width 0.1]); 
            
            str = sprintf('runtime= %d hours | star-hours= %d/%d', round(obj.runtime/3600), ...
                round(sum2(obj.histogram)/3600), round(sum2(obj.histogram_with_losses)/3600)); % any other info to show? 
            
            uicontrol(panel_info, 'Style', 'pushbutton', 'String', str, 'FontSize', 14, ...
                'Units', 'Normalized', 'Position', [margin/2 bottom 1-margin height]); 
            
            %%%%%%%%%%%%%%% choose types %%%%%%%%%%%%%%%%%%%%
            
            panel_types = uipanel(parent, 'title', 'types', ...
                'Units', 'Normalized', 'Position', [hist_width 0 1-hist_width 1]); 
            
            N = length(obj.cut_names); % number of type buttons we need
            
            total_seconds = sum2(obj.histogram_with_losses); 
            good_seconds = sum2(obj.histogram); 
            
            margin = 0.1; % left/right margin
            width = 1-margin;
            gap = 0.01; 
            height_total = 1/(2+N); % number of buttons
            height = height_total - gap; 
            
            font_size = 12;
            
            button = uicontrol(panel_types, 'Style', 'pushbutton', ...
                'String', sprintf('%-18s [%d%%]', 'good times', round(100*good_seconds/total_seconds)), ...
                'Units', 'Normalized', 'Position', [margin/2, 1-(gap/2+height), width, height], ...
                'Callback', @obj.callback_good, 'FontSize', font_size, 'UserData', input); 
            
            if input.good, button.ForegroundColor = color_on; end
            
            for ii = 1:N
                
                name = obj.cut_names{ii};
                if length(name)>13
                    name = [name(1:13) '...']; 
                end
                
                if input.exclusive
                    cut_seconds = sum2(obj.losses_exclusive(:,:,ii)); 
                else
                    cut_seconds = sum2(obj.losses_inclusive(:,:,ii)); 
                end
                
                button = uicontrol(panel_types, 'Style', 'pushbutton', ...
                    'String', sprintf('%-18s [%d%%]', name, round(100*cut_seconds/total_seconds)), ...
                    'Units', 'Normalized', 'Position', [margin/2, 1+gap/2-(gap+height)*(ii+1), width, height], ...
                    'Callback', @obj.callback_types, 'FontSize', font_size, 'UserData', {input, ii});
                
                if ismember(ii, input.types), button.ForegroundColor = color_on; end
            
                
            end
            
            button = uicontrol(panel_types, 'Style', 'pushbutton', 'String', 'inclusive', ...
                'Units', 'Normalized', 'Position', [margin/2, gap/2, width, height], ...
                'Callback', @obj.callback_exclusive, 'FontSize', font_size, 'UserData', input); 
            
            if input.exclusive, button.String = 'exclusive'; end
            
            
        end
        
        function callback_log(obj, hndl, ~)
            
            input = hndl.UserData; 
            
            input.log = ~input.log;
            
            obj.viewer('parent', input.parent); 
            
        end
        
        function callback_sum(obj, hndl, ~)
            
            input = hndl.UserData; 
            
            input.sum = ~input.sum;
            
            obj.viewer('parent', input.parent); 
            
        end
        
        function callback_hours(obj, hndl, ~)
            
            input = hndl.UserData; 
            
            input.hours = ~input.hours;
            
            obj.viewer('parent', input.parent); 
            
        end
        
        function callback_snr(obj, hndl, ~)
            
            input = hndl.UserData; 
            value = util.text.extract_numbers(hndl.String); 
            value = value{1};
            
            input.max_snr = value; 
            
            obj.viewer('parent', input.parent); 
            
        end
        
        function callback_alpha(obj, hndl, ~)
            
            input = hndl.UserData; 
            value = util.text.extract_numbers(hndl.String); 
            value = value{1};
            
            input.alpha = value; 
            
            obj.viewer('parent', input.parent); 
            
        end
        
        function callback_good(obj, hndl, ~)
            
            input = hndl.UserData; 
            
            input.good = ~input.good;
            
            obj.viewer('parent', input.parent); 
            
        end
        
        function callback_exclusive(obj, hndl, ~)
            
            input = hndl.UserData; 
            
            input.exclusive = ~input.exclusive;
            
            obj.viewer('parent', input.parent); 
            
        end
        
        function callback_types(obj, hndl, ~)
            
            input = hndl.UserData{1}; 
            type = hndl.UserData{2}; 
            
            if ismember(type, input.types)
                input.types(input.types==type) = []; % remove from the list
            else
                input.types = [input.types type]; % add to list
            end
            
            obj.viewer('parent', input.parent); 
            
        end
        
        function showHistogram(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('type', []); % choose cut type or leave empty to get the good hours
            input.input_var('exclusive', false); % if false, show inclusive losses
            input.input_var('max_snr', obj.snr_bin_max); 
            input.input_var('min_snr', obj.snr_bin_min); 
            input.input_var('rebin', []); % to be added later
            input.input_var('sum', false);
            input.input_var('hours', true); 
            input.input_var('axes', [], 'axis'); 
            input.input_var('font_size', 18); 
            input.input_var('log', true, 'use_log'); 
            input.input_var('color', []); 
            input.input_var('alpha', 0.7); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            if strcmp(input.axes.NextPlot, 'replace')
                reset(input.axes); 
            end
            
            convert = 1;
            
            if input.hours
                convert = 3600; 
            end
                
            if isempty(input.type)
                
                values = obj.histogram;
                str = sprintf('Effective star-time (%4.2f ', util.stat.sum2(values)/convert); 
                
            else
                if ischar(input.type)
                    input.type = obj.cut_indices.(input.type); 
                end

                if input.exclusive
                    values = obj.losses_exclusive(:,:,input.type); 
                    exc = '(exclusively)';
                    
                else
                    values = obj.losses_inclusive(:,:,input.type); 
                    exc = '(inclusively)'; 
                end
                
                str = sprintf('Time lost %s to "%s" (%4.2f ', exc, strrep(obj.cut_names{input.type}, '_', ' '), util.stat.sum2(values)/convert);
                
            end
            
            units = 'seconds)'; 
            
            if input.hours
                values = values./3600; % transform to hours instead of seconds!
                units = 'hours)'; 
            end
            
            if input.sum
                
                values = nansum(values,2);
                h = bar(input.axes, obj.snr_bin_edges(1:end-1), values, 'FaceAlpha',input.alpha); 
                
                if ~isempty(input.color)
                    h.FaceColor = input.color;
                end
                
                mx = nanmax(values); 
                
                if mx<=0
                    mx = 1;
                else
                    mx = mx.*1.1;
                end
                
                if strcmp(input.axes.NextPlot, 'replace')
                    input.axes.YLim = [0 mx]; 
                end
                
                if input.log
                    input.axes.YScale = 'log'; 
                end
                
                h.DisplayName = [str ' ' units];
            
                legend(input.axes, 'Location', 'NorthEast'); 
            
            else
                
                cla(input.axes); 
                
                imagesc(input.axes, 'XData', obj.snr_bin_edges(1:end-1)+obj.snr_bin_width/2, 'CData', values');
                ylabel(input.axes, 'Star index');                 
                colorbar(input.axes); 
                input.axes.YLim = [0 length(obj.star_bin_edges)-1]; 
                
                input.axes.YScale = 'linear'; 
                
                if input.log
                    input.axes.ColorScale = 'log'; 
                end
                
                legend(input.axes, 'off'); 
              
                util.plot.inner_title([str ' ' units], 'ax', input.axes, 'Position', 'NorthEast', 'Color', 'white'); 
            
            end
            
            input.axes.XLim = [input.min_snr, input.max_snr-obj.snr_bin_width];
            xlabel(input.axes, 'S/N'); 
            
            input.axes.FontSize = input.font_size;
            axis(input.axes, 'normal'); 
            axis(input.axes, 'xy'); 
            
        end
        
        function printReport(obj, varargin)
            
            import util.text.cs;
            
            input = util.text.InputVars;
            input.input_var('units', 'hours'); % can also choose "hours" or "percent"
            input.input_var('format', 'text'); % can aso choose "latex"
            input.scan_vars(varargin{:}); 
            
            if cs(input.units, 'seconds')
                unit_str = 's'; 
            elseif cs(input.units, 'hours')
                unit_str = 'h'; 
            elseif cs(input.units, 'percent', '%')
                unit_str = '%'; 
            else
                error('Unknown option "%s" to "units" argument. Try using "seconds" or "hours" or "percent"...', input.units); 
            end
            
            if cs(input.format, 'text')
                sep = '|';
                line = '--------------------------+--------------+------------- \n';
                ending = newline;
                code = @(str) str;
            elseif cs(input.format, 'latex')
                sep = '&';
                line = ['\\hline' newline]; 
                ending = ['\\' newline]; 
                code = @(str) sprintf('\\code{%s}', str); 
            else
                error('Unknown "format" option "%s". Use "text" or "latex" instead...', input.format); 
            end
            
            fprintf('%-25s %s inclusive[%s] %s exclusive[%s] %s', 'Cut name', sep, unit_str, sep, unit_str, ending); 
            fprintf(line); 
            
            useful = util.stat.sum2(obj.histogram); 
            total = util.stat.sum2(obj.histogram_with_losses); 
            
            % get the total star hours, not including losses
            
            for ii = 1:length(obj.cut_names)
                
                inc = util.stat.sum2(obj.losses_inclusive(:,:,ii));
                if unit_str=='h'
                    inc = inc/3600;
                elseif unit_str=='%'
                    inc = inc/total.*100; 
                end
                
                exc = util.stat.sum2(obj.losses_exclusive(:,:,ii));
                if unit_str=='h'
                    exc = exc/3600;
                elseif unit_str=='%'
                    exc = exc/total.*100; 
                end
                
                fprintf('%-25s %s %12.2f %s %12.2f %s', code(obj.cut_names{ii}), sep, inc, sep, exc, ending);
                
                
            end
            
            stars = util.stat.sum2(obj.losses_bad_stars); 
            
            if unit_str=='h'
                stars= stars/3600;
            elseif unit_str=='%'
                stars = stars/total.*100; 
            end
            
            fprintf('%-25s %s %12.2f %s         ---  %s', code('Bad stars'), sep, stars, sep, ending); 
            
            fprintf(line); 
            
            fprintf('%-25s %s %11.1fh %s %11.1fh %s', 'Good / total time', sep, useful/3600, sep, total/3600, ending); % always show it as hours!  

            
        end
        
    end
    
end







