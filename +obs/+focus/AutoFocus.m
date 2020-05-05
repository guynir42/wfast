classdef AutoFocus < handle

    properties(Transient=true)
        
        gui@obs.focus.gui.AutoGUI;
        aux_figure;
        cam; % optional link back to camera object
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        pos; % focuser position at each iteration
        widths; % one or multiple measurements per focus position
        weights; % weight of each sampling point (proportional to star flux)
        xy_pos; % position of sampling points
        xy_pos_reduced; % positions of good points
        
        % need to update these at some point
        x_max = 2160;
        y_max = 2560;
        
        pos_range_vector;
        
        fit_results = [];
        min_positions; % minimum of curves for each star
        min_weights; % the weight to give each star result
        min_positions_reduced; % replace NaNs with zeros
        min_weights_reduced; % average weights not including NaNs and negative values (and NaNs in min_positions)
        
        surface_coeffs; % results of 2D fit to surface (const, x coeff, y coeff)
        
        found_pos; % put your focuser to this position
        found_tip; % put your focuser to this tip value
        found_tilt; % put your focuser to this tilt value
        
    end
    
    properties % switches/controls
        
        expT = 0.03;
        frame_rate = NaN;
        batch_size = 5;
        
        use_loop_back = 0;
        
        cut_size = 25;
        aperture = 9; 
        annulus = [12 0]; 
        gaussian = 5;
        index_flux = 1; % the power index we set for the input flux when calculating the weight for each point
        
        use_fit_curves = 0; % use fit to 2nd order polynomial instead of the minimal position
        use_fit_tip_tilt = 0; % tell the camera/acquisition object to use the tip/tilt results
        
        step = 0.01; % step size for scanning the focus positions (mm)
        range = 0.15; % range on either side of the best focus position (mm)
        
        angle = -60; % between tip axis and pixel y axis (degrees)
        spider_diameter = 92.4; % in cm
        pixel_size = 12; % in microns
        
        num_plots = 12;
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        phot_struct; % save the results of the photometry made to get the widths/flux of each star
        
        version = 1.02;
        
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
            obj.min_weights_reduced = [];
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
        
        function val = bad_indices(obj)
            
            val = isnan(obj.min_weights) | obj.min_weights<=0 | isnan(obj.min_positions);
            
            if obj.use_fit_curves && ~isempty(obj.fit_results)
                coeffs = [obj.fit_results.coeffs];
                val = val | coeffs(3,:)'<0; % add fits for parabolas without a minimum
            end
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function val = getPosScanValues(obj, current_pos)
            
            val = current_pos + (-obj.range:obj.step:obj.range);
            
            obj.pos_range_vector = val;
            
        end

        function input(obj, idx, position, widths, fluxes, xy_pos)
            
            if nargin<5 || isempty(fluxes)
                fluxes = 1;
            end
            
            if nargin>=6 && ~isempty(xy_pos)
                obj.xy_pos = xy_pos;
            end
            
            obj.pos(idx) = position;
                
            widths = nanmean(widths, 1); 
            weights = nanmean(fluxes.^obj.index_flux, 1); 
            weights = weights./nanmean(weights); 
            
            obj.widths(:,idx) = util.vec.tocolumn(widths);
            obj.weights(:,idx) = util.vec.tocolumn(weights);
            
            obj.fit_results = [];
            
        end
        
        function calculate(obj)
            
            if isempty(obj.pos) || isempty(obj.widths)
                return;
            end
            
            if obj.use_fit_curves

                if isempty(obj.fit_results)
                    obj.fit_results = util.fit.polyfit(obj.pos', obj.widths', 'errors', 1./obj.weights, 'order', 2); 
                end
                
                coeffs = [obj.fit_results.coeffs];
                a = coeffs(3,:); 
                b = coeffs(2,:); 
                c = coeffs(1,:); 

                x_min = -b./2./a; % minimal position for each parabola
                y_min = c-b.^2./a./4; % depth of each parabola

                obj.min_positions = util.vec.tocolumn(x_min); 
                obj.min_weights = nanmedian(obj.weights,2)./util.vec.tocolumn(y_min); % the average flux is one measure of the goodness of that star but also the smallness of the minimal width! 

            else
                
                [mn, idx] = nanmin(obj.widths, [], 2); % find the minimal value for each star
                
                obj.min_positions = util.vec.tocolumn(obj.pos(idx)); 
                
                obj.min_weights = nanmedian(obj.weights,2)./mn; % the average flux is one measure of the goodness of that star but also the smallness of the minimal width! 
                
            end
            
            obj.min_weights = obj.min_weights./nanmean(obj.min_weights); % normalize the weights
            
            % get positions and weights after removing the bad measurements
            obj.min_positions_reduced = obj.min_positions(~obj.bad_indices);
            obj.min_weights_reduced = obj.min_weights(~obj.bad_indices); 
            obj.xy_pos_reduced = obj.xy_pos(~obj.bad_indices,:); 
            
        end
        
        function fitCurves(obj) % to be deprecated

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
            
            obj.min_weights_reduced = sum(obj.weights,2);
            
            bad_points = isnan(obj.min_positions) | obj.min_weights_reduced<=0 | isnan(obj.min_weights_reduced);
            
            obj.min_positions_reduced = obj.min_positions(~bad_points);
            obj.min_weights_reduced = obj.min_weights_reduced(~bad_points);
            obj.min_weights_reduced = obj.min_weights_reduced./mean(obj.min_weights_reduced);
            obj.xy_pos_reduced = obj.xy_pos(~bad_points,:);
            
        end
        
        function val = checkFit(obj, x, y, fr) %  to be deprecated
            
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
            
            % add check that there is anything to fit (maybe before the call to this function?)
            
            m = size(obj.xy_pos_reduced,1); % number of measurements
            
            B = obj.min_positions_reduced; % measured best position for each location
            w = obj.min_weights_reduced;  % total flux of each position 
            
            % rotate to the direction of the actuators...
            xy = obj.xy_pos_reduced - nanmean(obj.xy_pos_reduced, 1); 
            xy_rot = xy*[cosd(obj.angle), -sind(obj.angle); sind(obj.angle), cosd(obj.angle)]; 
            
            A = [ones(m,1) xy_rot(:,1) xy_rot(:,2)]; % design matrix! 

            obj.surface_coeffs = lscov(A,B,w); % coeffs are: piston, x (tilt) and y (tip)
                        
        end
        
        function findPosOnly(obj) % to be deprecated
            
            [~,idx] = nanmin(nanmean(obj.widths, 1));
            obj.found_pos = obj.pos(idx);
                
%             obj.found_pos = mean(obj.min_positions, 'omitnan');
        
        end
        
        function findPosTipTilt(obj)
            
%             obj.found_pos = obj.surface_coeffs(1);
            obj.found_pos = nansum(obj.min_positions.*obj.min_weights)./nansum(obj.min_weights);
%             if obj.debug_bit, fprintf('BEST POS: mean= %f | surface piston term= %f\n', nanmean(obj.min_positions), obj.surface_coeffs(1)); end
            
            obj.found_tip = obj.surface_coeffs(2).*obj.spider_diameter.*1e4./obj.pixel_size; % in this case tip is X slope
            obj.found_tilt = obj.surface_coeffs(3).*obj.spider_diameter.*1e4./obj.pixel_size; % in this case tip is Y slope
            
%             if obj.debug_bit, fprintf('BEST TIP= %f | BEST tilt= %f\n', obj.found_tip, obj.found_tilt); end
            
        end
        
        function str = printout(obj)
            
            str = sprintf('pos= %f | tip= %f | tilt= %f', obj.found_pos, obj.found_tip, obj.found_tilt);
            
            if ~isempty(obj.surface_coeffs)
                str = sprintf('%s | coeffs= %s', str, util.text.print_vec(obj.surface_coeffs, '  ')); 
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function plot(obj, varargin)
            
            if isempty(obj.pos) || isempty(obj.widths)
                return;
            end
            
            input = util.text.InputVars;
            input.input_var('font_size', 20); 
            input.input_var('ax', [], 'axes', 'axis'); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.ax)
                if ~isempty(obj.gui) && obj.gui.check
                    input.ax = obj.gui.axes_image;
                else
                    input.ax = gca;
                end
            end
            
            cla(input.ax);
            hold(input.ax, 'on');
            
            N = min(size(obj.widths,1), obj.num_plots); 

            mn = util.stat.min2(obj.widths(1:N,:));
            mx = util.stat.max2(obj.widths(1:N,:));
            
            for ii = 1:N

                input.ax.ColorOrderIndex = ii;

                h = plot(input.ax, obj.pos, obj.widths(ii,:), '-', 'LineWidth', 1.5); % show the raw data
%                 h.DisplayName = sprintf('pos= %4.3f | weight= %5.2f', obj.min_positions(ii), obj.min_weights(ii)); 

                if ~isempty(obj.min_positions)
                    plot(input.ax, obj.min_positions(ii), mn-0.05, 'v', 'MarkerSize', sqrt(abs(obj.min_weights(ii)))*5, 'Color', h.Color); 
                end

                if obj.use_fit_curves && ~isempty(obj.fit_results)
                    plot(input.ax, obj.pos, obj.fit_results(ii).ym, ':', 'Color', h.Color, 'LineWidth', 1); % show the 2nd order polynomial fit
                end
                
            end
            
            
            plot(input.ax, obj.pos, nanmean(obj.widths,1), 'LineWidth', 3, 'Color', 'k');
%             plot(input.ax, obj.pos, nanmean(obj.widths), 'LineWidth', 3, 'Color', 'k', 'DisplayName', 'average');
            
%             hl = legend(input.ax, 'Location', 'North', 'NumColumns', 3); 
%             hl.FontSize = 12;
            
            if ~isempty(obj.found_pos)
                plot(input.ax, obj.found_pos, input.ax.YLim, '--g'); 
            end

            hold(input.ax, 'off');
            
            xlabel(input.ax, 'focuser position (mm)');
            ylabel(input.ax, 'width (second moment)');
            
            if ~isempty(obj.pos_range_vector)
                input.ax.XLim = [min(obj.pos_range_vector), max(obj.pos_range_vector)];
            end
            
            box(input.ax, 'on'); 
            
            input.ax.YLim = [mn-0.1 mx+0.1];
            
            input.ax.FontSize = input.font_size;
            
        end
        
        function show(obj, varargin)
           
            obj.fitSurface;
            
            if isempty(obj.aux_figure) || ~isvalid(obj.aux_figure)
                obj.aux_figure = figure;
            else
                clf(obj.aux_figure); 
                figure(obj.aux_figure);
            end 
            
            ax = axes('Parent', obj.aux_figure);
            
            [X,Y] = meshgrid(-obj.x_max/2:obj.x_max/2, -obj.y_max/2:obj.y_max/2);
            
            % rotate to the direction of the actuators...
            xy = obj.xy_pos_reduced - nanmean(obj.xy_pos_reduced, 1); 
            xy_rot = xy*[cosd(obj.angle), -sind(obj.angle); sind(obj.angle), cosd(obj.angle)]; 
            
            util.plot.show(obj.surface_coeffs(1) + obj.surface_coeffs(2).*X + obj.surface_coeffs(3).*Y, 'ax', ax, 'fancy', 'off', ...
                'xvalues', min(xy_rot(:,1)):max(xy_rot(:,1)), 'yvalues', min(xy_rot(:,2)):max(xy_rot(:,2))); 
%                 'xvalues', -obj.x_max/2:obj.x_max/2, 'yvalues', -obj.y_max/2:obj.y_max/2);
            
            hold(ax, 'on');
            
            
            scatter(ax, xy_rot(:,1), xy_rot(:,2), obj.min_weights_reduced*10, obj.min_positions_reduced, 'filled');

            axis(ax, 'image');
            colorbar(ax); 
            
            hold(ax, 'off');

        end
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = obs.focus.gui.AutoGUI(obj); 
            end
            
            obj.gui.make;
            
        end
        
    end
    
end

