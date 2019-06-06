classdef Event < handle

    properties(Transient=true)
        
        dates; 
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        range_time_idx;
        range_kernels_idx;
        range_stars;
        
        peak_timestamp;
        best_kernel;
        star_idx;
        snr;
        
        timestamps;
        flux_filtered;
        flux_fitted;
        flux_raw_all;
        stds;
        
        which_batch = 'first'; % can be "first" or "second", depending on where is the peak
        which_frame = []; 
        
        
    end
    
    properties % switches/controls
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        cutouts;
        positions;
        stack;
        batch_index;
        filename;
        
    end
    
    properties(Hidden=true)
        
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
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Event(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.Event')
                if obj.debug_bit, fprintf('Event copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            elseif length(varargin)>=1 && istable(varargin{1})
                obj.input(varargin{:});
            else
                if obj.debug_bit, fprintf('Event constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
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
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function input(obj, table_row, fluxes, timestamps, kernels, stds)
            
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
            
            obj.star_idx = z;
            
            obj.flux_filtered = fluxes(:, x, z);
            
            obj.which_batch = 'first';
            obj.which_frame = y;
            if y>floor(size(fluxes,1)/2)
                obj.which_batch = 'second'; % the peak is in the second half of the two-batch fluxes vector
                obj.which_frame = y - floor(size(fluxes,1)/2); 
            end
            
        end
        
        function val = is_same(obj, other)
            
            val = 0;
            
            for ii = 1:length(other)
                
                if obj.star_idx==other(ii).star_idx
                    t_self = obj.timestamps(obj.range_time_idx);
                    t_other = other(ii).timestamps(other(ii).range_time_idx);
                    if ~isempty(intersect(t_self, t_other))
                        val = 1;
                        return;
                    end
                    
                end
                
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
        end
        
        function showFlux(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('font_size', 14);
            input.scan_vars(varargin{:});
            
            if isempty(input.ax), input.ax = gca; end
            
            input.ax.NextPlot = 'replace';
            h1 = plot(input.ax, obj.timestamps, obj.flux_filtered);
            h1.DisplayName = 'Filtered LC';
            
            input.ax.NextPlot = 'add';
            range = obj.range_time_idx;
            h2 = plot(input.ax, obj.timestamps(range), obj.flux_filtered(range));
            h2.DisplayName = 'trigger region';
            f = obj.flux_raw_all(:,obj.star_idx);
            m = mean(f);
            if ~isempty(obj.stds)
                s = obj.stds(obj.star_idx);
            else
                s = std(obj.flux_raw_all);
                s = s(obj.star_idx);
            end
            
            h3 = plot(input.ax, obj.timestamps, (f-m)./s, '-');
            h3.DisplayName = 'Raw LC';
            
            h4 = plot(input.ax, obj.kernel_timestamps+obj.peak_timestamp, obj.best_kernel*5, ':');
            h4.DisplayName = 'best filter';
            
            xlabel(input.ax, 'timestamp (seconds)');
            ylabel(input.ax, 'flux S/N');   
            
            if strcmp(obj.which_batch, 'first')
                lh = legend(input.ax, 'Location', 'NorthEast');
            else
                lh = legend(input.ax, 'Location', 'SouthWest');
            end
            
            lh.FontSize = input.font_size;
            
            util.plot.inner_title(sprintf('star: %d | batches: %d-%d | S/N= %4.2f | \\sigma= %4.2f', ...
                obj.star_idx, obj.batch_index_first, obj.batch_index_second, obj.snr, obj.stds(obj.star_idx)),...
                'ax', input.ax, 'Position', 'NorthWest', 'FontSize', input.font_size);
            
            input.ax.YLim(1) = -input.ax.YLim(2);
%             input.ax.YLim = input.ax.YLim*1.2;
            input.ax.XLim = [obj.timestamps(1) obj.timestamps(end)];

            input.ax.NextPlot = 'replace';
            
        end
        
        function showStack(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.scan_vars(varargin{:});
            
            if isempty(input.ax), input.ax = gca; end
            
            if ~isempty(obj.stack)
                
                util.plot.show(obj.stack, 'fancy', 0, 'autodyn', 1, 'ax', input.ax); 
                
                if ~isempty(obj.cutouts) && ~isempty(obj.positions)
                    
                    S = size(obj.cutouts);
                    S = S(1:2); 
                    
                    P = obj.positions(obj.star_idx, :); 
                    
                    C = P - floor(S/2) - 0.5; % corner of the rectangle
                    
                    rectangle('Position', [C S], 'Parent', input.ax); 
                    
                end
                
            end
            
        end
        
        function showCutouts(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('number', 9);
            input.scan_vars(varargin{:});
            
            if isempty(input.ax) 
                fig = gcf;
                delete(fig.Children);
                input.ax = axes('Parent', fig);
            end
            
            cutouts = cat(3, obj.cutouts_first(:,:,:,obj.star_idx), obj.cutouts_second(:,:,:,obj.star_idx));
            
            if ~isempty(obj.cutouts)
                
                idx = obj.which_frame; 
                if strcmp(obj.which_batch, 'second')
                    idx = idx + size(obj.cutouts,3);
                end
                
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
                
            end
            
%             fprintf('idx= %d | idx_start= %d | idx_end= %d\n', idx, idx_start, idx_end);
            
            input.ax.Units = 'pixels';
            pos = input.ax.Position;
            if pos(3)>pos(4)
                pos(3) = pos(4);
            elseif pos(4)>pos(3)
                pos(4) = pos(3);
            end
            parent = input.ax.Parent;
            panel = uipanel('Units', 'pixels', 'Position', pos, 'Parent', parent);
            delete(input.ax);
            
            Nrows = ceil(sqrt(input.number));
            Ncols = Nrows;
            
            for ii = 1:input.number
                
                x = mod(ii-1, Nrows);
                y = floor((ii-1)/Nrows);
                
                ax{ii} = axes('Position', [x/Ncols y/Nrows 1/Ncols 1/Nrows], 'Parent', panel);
                
                util.plot.show(cutouts(:,:,idx_start+ii-1), 'fancy', 0, 'autodyn', 0);
                
                if idx_start+ii-1==idx
                    util.plot.inner_title([num2str(idx_start+ii-1) ' (peak)'], 'Position', 'NorthWest');
                else
                    util.plot.inner_title(num2str(idx_start+ii-1), 'Position', 'NorthWest');
                end
                
            end

        end
        
    end    
    
end

