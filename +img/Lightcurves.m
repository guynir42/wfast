classdef Lightcurves < handle

    properties(Transient=true)
        
        gui@img.gui.LightGUI;
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        frame_index = 1;
        num_frames = 0;
        total_frames;
        
    end
    
    properties % switches/controls
        
        use_subtract_backgrounds = 0;
        
        show_what = 'flux'; % can choose "flux", "weight", "offsets", "widths" or "backgrounds"
        show_flux_cal = 'raw'; % can choose "raw" or "cal" or "both".
        show_num_stars = 10; % up to this number of stars are shown.
        use_smooth = 1;
        smooth_interval = 10;
        
        use_double_up = 1; % choose if you want to expand the data storage by factor of 2 each time when space runs out... 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        timestamps;
        fluxes;
        fluxes_cal;
        weights;
        offsets_x;
        offsets_y;
        centroids_x;
        centroids_y;
        widths;
        backgrounds;
        
    end
    
    properties(Hidden=true)
        
        fluxes_full;
        timestamps_full;
        weights_full;
        offsets_x_full;
        offsets_y_full;
        centroids_x_full;
        centroids_y_full;
        widths_full;
        backgrounds_full;
        
        show_what_list = {'flux', 'weight', 'offsets', 'centroids', 'widths', 'backgrounds'};
        show_cal_list = {'raw', 'cal', 'both'};
        
        version = 1.02;
        
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
            obj.timestamps_full = [];
            obj.weights_full = [];
            obj.offsets_x_full = [];
            obj.offsets_y_full = [];
            obj.widths_full = [];
            obj.backgrounds_full = [];
            
            obj.frame_index = 1;
            
        end
        
    end
    
    methods % getters
        
        function val = get.fluxes(obj)
            
            val = obj.fluxes_full;
            
            if isempty(val)
                return;
            end
            
            val = val(1:obj.frame_index-1,:);
            
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
                
                f_frames_average = median(f,2, 'omitnan'); 
                f_frames_average_norm = f_frames_average./median(f_frames_average,1, 'omitnan'); 
                val = f./f_frames_average_norm;
            
            end
        end
        
        function val = get.timestamps(obj)
            
            val = obj.timestamps_full;
            
            if isempty(val)
                return;
            end
            
            val = val(1:obj.frame_index-1);
            
            if all(isnan(val))
                val = 1:obj.num_frames;
            end
            
        end
        
        function val = get.weights(obj)
            
            val = obj.weights_full;
            
            if isempty(val)
                return;
            end
            
            val = val(1:obj.frame_index-1,:);
            
            if all(isnan(val))
                val = [];
            end
            
        end
        
        function val = get.offsets_x(obj)
            
            val = obj.offsets_x_full;
            
            if isempty(val)
                return;
            end
            
            val = val(1:obj.frame_index-1,:);
            
            if all(isnan(val))
                val = [];
            end
            
        end
        
        function val = get.offsets_y(obj)
            
            val = obj.offsets_y_full;
            
            if isempty(val)
                return;
            end
            
            val = val(1:obj.frame_index-1,:);
            
            if all(isnan(val))
                val = [];
            end
            
        end
        
        function val = get.centroids_x(obj)
            
            val = obj.centroids_x_full;
            
            if isempty(val)
                return;
            end
            
            val = val(1:obj.frame_index-1,:);
            
            if all(isnan(val))
                val = [];
            end
            
        end
        
        function val = get.centroids_y(obj)
            
            val = obj.centroids_y_full;
            
            if isempty(val)
                return;
            end
            
            val = val(1:obj.frame_index-1,:);
            
            if all(isnan(val))
                val = [];
            end
            
        end
        
        function val = get.widths(obj)
            
            val = obj.widths_full;
            
            if isempty(val)
                return;
            end
            
            val = val(1:obj.frame_index-1,:);
            
            if all(isnan(val))
                val = [];
            end
            
        end
        
        function val = get.backgrounds(obj)
            
            val = obj.backgrounds_full;
            
            if isempty(val)
                return;
            end
            
            val = val(1:obj.frame_index-1,:);
            
            if all(isnan(val))
                val = [];
            end
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function startup(obj, total_frames) % do we need this??
            
            obj.total_frames = total_frames;
            
            obj.frame_index = 1;
            
        end
        
        function input(obj, varargin)
            
            import util.vec.insert_matrix;
            
            if isscalar(varargin) && isa(varargin{1}, 'util.text.InputVars')
                input = varargin{1};
            else
                input = util.text.InputVars;
                input.use_ordered_numeric = 1;
                input.input_var('fluxes', [], 'fluxes_raw');
                input.input_var('timestamps', [], 'times');
                input.input_var('weights', []);
                input.input_var('centroids_x', []);
                input.input_var('centroids_y', []);
                input.input_var('offsets_x', []);
                input.input_var('offsets_y', []);
                input.input_var('widths', []);
                input.input_var('backgrounds', []);
                input.scan_vars(varargin{:});
            end
            
            N = size(input.fluxes,1);
            
            obj.fluxes_full = insert_matrix(obj.fluxes_full, input.fluxes, [obj.frame_index,1], NaN, obj.use_double_up);
            obj.timestamps_full = insert_matrix(obj.timestamps_full, input.timestamps, [obj.frame_index,1], NaN, obj.use_double_up);
            obj.weights_full = insert_matrix(obj.weights_full, input.weights, [obj.frame_index,1], NaN, obj.use_double_up);
            obj.offsets_x_full = insert_matrix(obj.offsets_x_full, input.offsets_x, [obj.frame_index,1], NaN, obj.use_double_up);
            obj.offsets_y_full = insert_matrix(obj.offsets_y, input.offsets_y, [obj.frame_index,1], NaN, obj.use_double_up);
            obj.centroids_x_full = insert_matrix(obj.centroids_x_full, input.centroids_x, [obj.frame_index,1], NaN, obj.use_double_up);
            obj.centroids_y_full = insert_matrix(obj.centroids_y_full, input.centroids_y, [obj.frame_index,1], NaN, obj.use_double_up);
            obj.widths_full = insert_matrix(obj.widths_full, input.widths, [obj.frame_index,1], NaN, obj.use_double_up);
            obj.backgrounds_full = insert_matrix(obj.backgrounds_full, input.backgrounds, [obj.frame_index,1], NaN, obj.use_double_up);
            
            obj.frame_index = obj.frame_index + N;
            
        end
        
        function getData(obj, photometry, type)
            
            if nargin<3 || isempty(type)
                type = '';
            end
            
            list = {'fluxes', 'weights', 'offsets_x', 'offsets_y', 'centroids_x', 'centroids_y', 'widths', 'backgrounds'};
            
            if ~isempty(type)
                list2 = strcat(list, ['_' type]);
            else
                list2 = list;
            end
            
            input = util.text.InputVars;
            
            for ii = 1:length(list)
                input.input_var(list{ii}, photometry.(list2{ii}));
            end
            
            input.input_var('timestamps', photometry.timestamps);
            
            obj.input(input);
            
        end
        
        function diskDump(obj, filename)
            
            timestamps = obj.timestamps;
            fluxes = obj.fluxes;
            widths = obj.widths;
            weights = obj.weights;
            centroids_x = obj.centroids_x;
            centroids_y = obj.centroids_y;
            
            save(filename, 'timestamps', 'fluxes', 'widths', 'weights', 'centroids_x', 'centroids_y');
            
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
                    input.ax = obj.gui.axes_image;
                else
                    input.ax = gca;
                end
            end
            
            delete(allchild(input.ax));
            delete(findobj(input.ax.Parent, 'type', 'legend'));
            
            if cs(obj.show_what, 'fluxes')
                
                if cs(obj.show_flux_cal, 'raw')
                    obj.addPlots(input.ax, obj.fluxes);
                elseif cs(obj.show_flux_cal, 'cal')
                    obj.addPlots(input.ax, obj.fluxes_cal);
                elseif cs(obj.show_flux_cal, 'both')
                    obj.addPlots(input.ax, obj.fluxes, ':');
                    obj.addPlots(input.ax, obj.fluxes_cal, '-');
                    legend(input.ax, 'flux (uncalibrated)', 'flux (self-calibrated)', 'location', 'NorthEast');                    
                end
                
                ylabel(input.ax, 'flux (counts)');
            
            elseif cs(obj.show_what, 'weights')
                obj.addPlots(input.ax, obj.weights);
                ylabel(input.ax, 'weights (pixels in aperture)');
            elseif cs(obj.show_what, 'offsets')
                obj.addPlots(input.ax, obj.offsets_x, '-');
                obj.addPlots(input.ax, obj.offsets_y, ':');
                ylabel(input.ax, 'offset inside cutouts (pixels)');
            elseif cs(obj.show_what, 'centroids')
                obj.addPlots(input.ax, obj.centroids_x, '-');
                obj.addPlots(input.ax, obj.centroids_y, ':');
                ylabel(input.ax, 'centroid positions in image (pixels)');
            elseif cs(obj.show_what, 'widths')
                obj.addPlots(input.ax, obj.widths);
                ylabel(input.ax, 'widths using 2nd moment (pixels)');
            elseif cs(obj.show_what, 'backgrounds')
                obj.addPlots(input.ax, obj.backgrounds);
                ylabel(input.ax, 'background (counts/pixel)');
            else
                error('Unknown data to show "%s", use "fluxes" or "offset" etc...', obj.show_what);
            end
            
            hold(input.ax, 'off');
            
            xlabel(input.ax, 'timestamps (seconds)');
            
        end
        
        function addPlots(obj, ax, data, line_str)
            
            if nargin<4 || isempty(line_str)
                line_str = '-';
            end
            
            if isempty(data)
                return;
            end
            
            if obj.use_smooth
                data_smoothed = util.img.conv_f(ones(obj.smooth_interval,1)./obj.smooth_interval, data, 'crop', 'same'); 
            end
            
            ax.NextPlot = 'add';
            ax.ColorOrderIndex = 1;
            
            for ii = 1:obj.show_num_stars
                
                if obj.use_smooth
                    h = scatter(ax, obj.timestamps, data(:,ii), 1, 'o', 'MarkerEdgeAlpha', 0.3);
                    if ~isempty(h), h.HandleVisibility = 'off'; end
                else
                    h = plot(ax, obj.timestamps, data(:,ii), line_str, 'LineWidth', 1);
                    if ii>1 && ~isempty(h), h.HandleVisibility = 'off'; end
                end
                
                if ~isempty(h)
                    h.UserData = ii;
                    h.ButtonDownFcn = @obj.callback_cutout_number;
                end
                
                if obj.use_smooth
                    step = obj.smooth_interval;                    
                    h = plot(ax, obj.timestamps(step:1:end-step), data_smoothed(step:1:end-step, ii), line_str, 'LineWidth', 3, 'Color', h.CData);
                    if ii>1 && ~isempty(h), h.HandleVisibility = 'off'; end
                    if ~isempty(h)
                        h.UserData = ii;
                        h.ButtonDownFcn = @obj.callback_cutout_number;
                    end
                end
                
            end
            
        end
        
        function callback_cutout_number(obj, hndl, ~)
            
            if obj.gui.check
                obj.gui.button_cut_number.String = ['cut: ' num2str(hndl.UserData)];
            else
                disp(['Cutout number is ' num2str(hndl.UserData)]);
            end
            
        end
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = img.gui.LightGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
    end
    
end
