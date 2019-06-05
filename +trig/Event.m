classdef Event < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        range_timestamps;
        range_kernels;
        range_stars;
        
        peak_timestamp;
        best_kernel;
        star_idx;
        snr;
        
        timestamps;
        flux_filtered;
        flux_fitted;
        flux_raw_all;
        
    end
    
    properties % switches/controls
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
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
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function input(obj, table_row, fluxes, timestamps, kernels)
            
            box = table_row{1, 'BoundingBox'};
            idx = table_row{1, 'VoxelIdxList'}{1};
            val = table_row{1, 'VoxelValues'}{1};
            
            if nargin>3 && ~isempty(timestamps)
                obj.timestamps = timestamps; 
            end
            
            obj.range_timestamps = ceil(box(2)):floor(box(2))+box(5);
            obj.range_kernels = ceil(box(1)):floor(box(1))+box(4); 
            obj.range_stars = ceil(box(3)):floor(box(3))+box(6); 
            
            [obj.snr, max_idx] = max(val); 
            [y,x,z] = ind2sub(size(fluxes), idx(max_idx));
            obj.peak_timestamp = timestamps(y);
            
            if nargin>4 && ~isempty(kernels)
                obj.best_kernel = kernels(:,x);
            end
            
            obj.star_idx = z;
            
            obj.flux_filtered = fluxes(:, x, z);
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.scan_vars(varargin{:});
            
            if isempty(input.ax), input.ax = gca; end
            
            input.ax.NextPlot = 'replace';
            plot(input.ax, obj.timestamps, obj.flux_filtered);
            
            input.ax.NextPlot = 'add';
            range = obj.range_timestamps;
            plot(input.ax, obj.timestamps(range), obj.flux_filtered(range));
            
        end
        
    end    
    
end

