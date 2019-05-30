classdef AutoFocus < handle

    properties(Transient=true)
        
        fig; % for displaying the focus curve
        ax; 
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        pos; % focuser position at each iteration
        widths; % one or multiple measurements per focus position
        weights; % weight of each sampling point (proportional to star flux)
        weights_reduced; % sum all exposures, not including NaNs and negative values (and NaNs in min_positions)
        xy_pos; % position of sampling points
        xy_pos_reduced; % positions of good points
        
        % need to update these at some point
        x_max = 2160;
        y_max = 2560;
        
        fit_results = {};
        min_positions; % minimum of curves for each star
        min_positions_reduced; % replace NaNs with zeros
        
        surface_coeffs;
        
        found_pos; % put your focuser to this position
        found_tip; % put your focuser to this tip value
        found_tilt; % put your focuser to this tilt value
        
    end
    
    properties % switches/controls
        
        use_fit_tip_tilt = 0;
       
        step = 0.01;
        range = 0.1;
        
        angle = 0; % between tip axis and pixel y axis (degrees)
        spider_diameter = 100; % in cm
        pixel_size = 6.5; % in microns
        num_pixels = 2000; % across the sensor (roughly)
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = AutoFocus(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.focus.AutoFocus')
                if obj.debug_bit, fprintf('AutoFocus copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('AutoFocus constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.pos = [];
            obj.widths = [];
            obj.weights = [];
            obj.weights_reduced = [];
            obj.xy_pos = [];
            obj.xy_pos_reduced = [];
            
            obj.fit_results = {};
            obj.min_positions = [];
            obj.min_positions_reduced = [];
            obj.surface_coeffs = [];
            
            obj.found_pos = [];
            obj.found_tip = [];
            obj.found_tilt = [];
            
            obj.clear;
            
        end
        
        function clear(obj)
            
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function startup(obj, focuser_pos, xy_positions, num_stars)
            
            obj.pos = focuser_pos + (-obj.range:obj.step:obj.range);
            obj.widths = NaN(num_stars, length(obj.pos));
            obj.weights = NaN(num_stars, length(obj.pos));
            obj.min_positions = NaN(num_stars,1);
            obj.xy_pos = xy_positions;
            
        end
        
        function input(obj, idx, position, widths, fluxes)
            
            obj.pos(idx) = position;
            obj.widths(:,idx) = util.vec.tocolumn(widths);
            obj.weights(:,idx) = util.vec.tocolumn(fluxes);
            
        end
        
        function calculate(obj)
            
            obj.fitCurves;
            if obj.use_fit_tip_tilt
                obj.fitSurface;
                obj.findPosTipTilt;
            else
                obj.findPosOnly;
            end
            
        end
        
        function fitCurves(obj)
            
            obj.fit_results = {};
            
            for ii = 1:size(obj.widths, 1) % number of stars
                
                x = obj.pos';
                y = obj.widths(ii,:)';
                
                % remove the NaN values
                x(isnan(y)) = [];
                y(isnan(y)) = [];
                
                if length(y)<5 % not enough points to do the fit
                    continue;
                end
                
                try % fit failed for some reason
                    obj.fit_results{ii} = fit(x, y, 'poly2');
                catch ME
                    warning(ME.getReport);
                    continue;
                end
                
                if ~obj.checkFit(x, y, obj.fit_results{ii}) % fit is not good enough, or doesn't have a minimum
                    continue;
                end
                
                obj.min_positions(ii) = -obj.fit_results{ii}.p2./(2*obj.fit_results{ii}.p1);
                
            end
            
            obj.weights_reduced = sum(obj.weights,2);
            
            bad_points = isnan(obj.min_positions) | obj.weights_reduced<=0 | isnan(obj.weights_reduced);
            
            obj.min_positions_reduced = obj.min_positions(~bad_points);
            obj.weights_reduced = obj.weights_reduced(~bad_points);
            obj.weights_reduced = obj.weights_reduced./mean(obj.weights_reduced);
            obj.xy_pos_reduced = obj.xy_pos(~bad_points,:);
            
        end
        
        function val = checkFit(obj, x, y, fr)
            
            if fr.p1<0 % no minimum
                val = 0;
                return;
            end
            
            y_model = feval(fr, x);
            chi2 = sum((y_model-y).^2); 
            
            if chi2>length(x)*2 % what would be a reasonable check on chi2 if the errors are not normalized??
                val = 0;
                return;
            end
            
            val = 1;
            
        end
        
        function fitSurface(obj)
            
            % add check that there is anything to fit (maybe before the
            % call to this function?)
            
            m = size(obj.xy_pos_reduced,1); % number of measurements
            
            B = obj.min_positions_reduced; % measured best position for each location
            w = obj.weights_reduced;  % total flux of each position 
            
            % rotate to the direction of the actuators... 
            xy_rot = obj.xy_pos_reduced*[cosd(obj.angle), -sind(obj.angle); sind(obj.angle), cosd(obj.angle)]; 
            
            A = [ones(m,1) xy_rot(:,1) xy_rot(:,2)]; % design matrix! 

            obj.surface_coeffs = lscov(A,B,w); % coeffs are: piston, x (tilt) and y (tip)
            
        end
        
        function findPosOnly(obj)
            
            obj.found_pos = mean(obj.min_positions, 'omitnan');
        
        end
        
        function findPosTipTilt(obj)
            
            obj.found_pos = obj.surface_coeffs(1);
            obj.found_pos = mean(obj.min_positions, 'omitnan');
            if obj.debug_bit, fprintf('BEST POS: mean= %f | surface piston term= %f\n', mean(obj.min_positions, 'omitnan'), obj.surface_coeffs(1)); end
            
            obj.found_tilt = obj.surface_coeffs(2).*obj.spider_diameter.*1e4./obj.pixel_size;
            obj.found_tip = obj.surface_coeffs(3).*obj.spider_diameter.*1e4./obj.pixel_size;
            
            if obj.debug_bit, fprintf('BEST TIP= %f | BEST tilt= %f\n', obj.found_tip, obj.found_tilt); end
            
        end
            
    end
    
    methods % plotting tools / GUI
        
        function plot(obj, varargin)
            
            if isempty(obj.fig) || ~isvalid(obj.fig)
                obj.fig = figure;
            end
            
            if isempty(obj.ax) || ~isvalid(obj.ax)
                obj.ax = axes('Parent', obj.fig);
            end
            
            cla(obj.ax);
            hold(obj.ax, 'on');
            
            N = min(size(obj.widths,2), 10); 
            
            for ii = 1:N
                
                h(ii) = plot(obj.ax, obj.pos, obj.widths(ii,:));
                
                if length(obj.fit_results)>=ii
                    plot(obj.ax, obj.pos, feval(obj.fit_results{ii}, obj.pos), ':', 'Color', h(ii).Color);
                end
                
            end
            
            hold(obj.ax, 'off');
            
            xlabel(obj.ax, 'focuser position (mm)');
            ylabel(obj.ax, 'width (second moment)');
            
%             legend(h, {'center', 'upper left', 'lower left', 'upper right', 'lower right'}, 'Parent', obj.fig);
            
        end
        
        function show(obj, varargin)
           
            if isempty(obj.fig) || ~isvalid(obj.fig)
                obj.fig = figure;
            end 
            
            if isempty(obj.ax) || ~isvalid(obj.ax)
                obj.ax = axes('Parent', obj.fig);
            end
            
            [X,Y] = meshgrid(1:obj.x_max, 1:obj.y_max);
            
            util.plot.show(obj.surface_coeffs(1) + obj.surface_coeffs(2).*Y + obj.surface_coeffs(3).*X);
            hold(obj.ax, 'on');
            scatter(obj.xy_pos_reduced(:,1), obj.xy_pos_reduced(:,2), obj.weights_reduced*10, obj.min_positions_reduced, 'filled'); 
%             axis(obj.ax, 'image');
%             colorbar(obj.ax); 
%             obj.ax.XLim = [1,obj.x_max];
%             obj.ax.YLim = [1,obj.y_max];
            
            hold(obj.ax, 'off');

        end
        
    end
    
end

