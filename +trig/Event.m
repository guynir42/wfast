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
        
        which_frame = 'first'; % can be 1 or 2, depending on where is the peak
        
        
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
            
            val = obj.(['cutouts_' obj.which_frame]);
            
        end
        
        function val = get.positions(obj)
            
            val = obj.(['positions_' obj.which_frame]);
            
        end
        
        function val = get.stack(obj)
            
            val = obj.(['stack_' obj.which_frame]);
            
        end
        
        function val = get.batch_index(obj)
            
            val = obj.(['batch_index_' obj.which_frame]);
            
        end
        
        function val = get.filename(obj)
            
            val = obj.(['filename_' obj.which_frame]);
            
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
                obj.stds = stds;
            end
            
            obj.star_idx = z;
            
            obj.flux_filtered = fluxes(:, x, z);
            
            if y>floor(size(fluxes,1)/2)
                obj.which_frame = 'second'; % the peak is in the second half of the two-batch fluxes vector
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
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
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
            s = obj.stds(obj.star_idx);
            
            h3 = plot(input.ax, obj.timestamps, (f-m)./s, '-');
            h3.DisplayName = 'Raw LC';
            
            h4 = plot(input.ax, obj.kernel_timestamps+obj.peak_timestamp, obj.best_kernel*5, ':');
            h4.DisplayName = 'best filter';
            
            xlabel(input.ax, 'timestamp (seconds)');
            ylabel(input.ax, 'flux S/N');   
            legend(input.ax);
            
            
            
        end
        
    end    
    
end

