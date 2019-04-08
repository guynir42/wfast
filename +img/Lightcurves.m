classdef Lightcurves < handle

    properties(Transient=true)
        
        gui;
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        frame_index = 0;
        num_frames = 0;
        
    end
    
    properties % switches/controls
        
        use_subtract_backgrounds = 0;
        
        show_flux_cal = 'raw'; % can choose "raw" or "cal" or "both".
        show_what = 'flux'; % can choose "flux", "weight", "offsets", "widths" or "backgrounds"
        use_show_full = 0; % if turned on, will show entire run, if off will show just what was already recorded. 
        show_num_stars = 10; % up to this number of stars are shown.
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        timestamps;
        fluxes;
        fluxes_cal;
        weights;
        offsets_x;
        offsets_y;
        widths;
        backgrounds;
        
    end
    
    properties(Hidden=true)
        
        fluxes_full;
        timestamps_full;
        weights_full;
        offsets_x_full;
        offsets_y_full;
        widths_full;
        backgrounds_full;
        
        show_cal_list = {'raw', 'cal', 'both'};
        sow_what_list = {'flux', 'weight', 'offsets', 'widths', 'backgrounds'};
        
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = Lightcurves(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'Lightcurves')
                if obj.debug_bit, fprintf('Lightcurves copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Lightcurves constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.fluxes_full = [];
            obj.weights_full = [];
            obj.offsets_x_full = [];
            obj.offsets_y_full = [];
            obj.widths_full = [];
            obj.backgrounds_full = [];
            
            obj.frame_index = 0;
            
        end
        
    end
    
    methods % getters
        
        function val = get.fluxes(obj)
            
            val = obj.fluxes_full(1:obj.frame_index,:);
            
            if all(isnan(val))
                val = [];
            end
            
        end
        
        function val = get.fluxes_cal(obj)
            
            if isempty(obj.fluxes)
                val = [];
            else
            
                f = obj.fluxes;
                
                if obj.use_subtract_backgrounds
                    f = f - obj.backgrounds.*obj.weights;
                end
                
                f_frames_average = median(f,2); 
                f_frames_average_norm = f_frames_average./median(f_frames_average,1); 
                val = f./f_frames_average_norm;
            
            end
        end
        
        function val = get.timestamps(obj)
            
            val = obj.timestamps_full(1:obj.frame_index,:);
            
            if all(isnan(val))
                val = 1:obj.num_frames;
            end
            
        end
        
        function val = get.weights(obj)
            
            val = obj.weights_full(1:obj.frame_index,:);
            
            if all(isnan(val))
                val = [];
            end
            
        end
        
        function val = get.offsets_x(obj)
            
            val = obj.offsets_x_full(1:obj.frame_index,:);
            
            if all(isnan(val))
                val = [];
            end
            
        end
        
        function val = get.offsets_y(obj)
            
            val = obj.offsets_y_full(1:obj.frame_index,:);
            
            if all(isnan(val))
                val = [];
            end
            
        end
        
        function val = get.widths(obj)
            
            val = obj.widths_full(1:obj.frame_index,:);
            
            if all(isnan(val))
                val = [];
            end
            
        end
        
        function val = get.backgrounds(obj)
            
            val = obj.backgrounds_full(1:obj.frame_index,:);
            
            if all(isnan(val))
                val = [];
            end
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function startup(obj, num_frames, num_stars)
            
            obj.fluxes_full = NaN(num_frames, num_stars);
            obj.fluxes_cal_full = NaN(num_frames, num_stars);
            obj.weights_full = NaN(num_frames, num_stars);
            obj.offsets_x_full = NaN(num_frames, num_stars);
            obj.offsets_y_full = NaN(num_frames, num_stars);
            obj.widths_full = NaN(num_frames, num_stars);
            obj.backgrounds_full = NaN(num_frames, num_stars);
            
            obj.frame_index = 0;
            
        end
        
        function input(obj, varargin)
            
            if isscalar(varargin) && isa(varargin{1}, 'util.text.InputVars')
                input = varargin{1};
            else
                input = util.text.InputVars;
                input.use_ordered_numeric = 1;
                input.input_var('fluxes', [], 'fluxes_raw');
                input.input_var('timestamps', [], 'times');
                input.input_var('weights', []);
                input.input_var('offsets_x', []);
                input.input_var('offsets_y', []);
                input.input_var('widths', []);
                input.input_var('backgrounds', []);
                input.scan_vars(varargin{:});
            end
            
            N = size(input.fluxes,1);
            
            obj.fluxes_full(obj.frame_index+1:obj.frame_index+N,:) = input.fluxes;
            
            if ~isempty(input.timestamps)
                obj.timestamps_full(obj.frame_index+1:obj.frame_index+N,:) = input.timestamps;
            end
            
            if ~isempty(input.weights)
                obj.weights_full(obj.frame_index+1:obj.frame_index+N,:) = input.weights;
            end
            
            if ~isempty(input.offsets_x)
                obj.offsets_x_full(obj.frame_index+1:obj.frame_index+N,:) = input.offsets_x;
            end
            
            if ~isempty(input.offsets_y)
                obj.offsets_y_full(obj.frame_index+1:obj.frame_index+N,:) = input.offsets_y;
            end
            
            if ~isempty(input.widths)
                obj.widths_full(obj.frame_index+1:obj.frame_index+N,:) = input.widths;
            end
            
            if ~isempty(input.backgrounds)
                obj.backgrounds_full(obj.frame_index+1:obj.frame_index+N,:) = input.backgrounds;
            end
            
            obj.frame_index = obj.frame_index + N;
            
        end
        
        function getData(obj, photometry, type)
            
            if nargin<3 || isempty(type)
                type = '';
            end
            
            list = {'fluxes', 'timestamps', 'weights', 'offsets_x', 'offsets_y', 'widths', 'backgrounds'};
            
            if ~isempty(type)
                list2 = strcat(list, type);
            else
                list2 = list;
            end
            
            input = util.text.InputVars;
            
            for ii = 1:length(list)
                input.input_var(list{ii}, photometry.(list2{ii}));
            end
            
            obj.input(input);
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
            import util.text.cs;
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.scan_vars(varargin{:});
            
            if isempty(input.ax)
                if ~isempty(obj.gui) && obj.gui.check
                    input.ax = obj.gui.image_axes;
                else
                    input.ax = gca;
                end
            end
            
            cla(input.ax);
            
            if cs(obj.show_what, 'fluxes')
                
                if cs(obj.show_cal, 'raw')
                    addPlots(input.ax, obj.fluxes);
                elseif cs(obj.show_cal, 'cal')
                    addPlots(input.ax, obj.fluxes_cal);
                elseif cs(obj.show_cal, 'both')
                    obj.plotData(input.ax, obj.fluxes, ':');
                    obj.plotData(input.ax, obj.fluxes_cal, '-');
                end
                
            elseif cs(obj.show_what, 'weights')
                obj.addPlots(input.ax, obj.weights);
            elseif cs(obj.show_what, 'offsets')
                obj.addPlots(input.ax, obj.offsets_x, '-');
                obj.addPlots(input.ax, obj.offsets_y, ':');
            elseif cs(obj.show_what, 'widths')                
                obj.addPlots(input.ax, obj.widths);
            elseif cs(obj.show_what, 'backgrounds')
                obj.addPlots(input.ax, obj.backgrounds);
            else
                error('Unknown data to show "%s", use "fluxes" or "offset" etc...', obj.show_what);
            end
            
            hold(input.ax, 'off');
            
        end
        
        function addPlots(obj, ax, data, line_str)
            
            if nargin<4 || isempty(line_str)
                line_str = '-';
            end
            
            ax.NextPlot = 'add';
            ax.ColorOrderIndex = 1;
            
            for ii = 1:obj.show_num_stars
                plot(ax, obj.timestamps, data, line_str);
            end
            
        end
        
    end    
    
end
