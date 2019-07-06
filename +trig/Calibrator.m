classdef Calibrator < handle

    properties(Transient=true)
        
        gui; 
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        timestamps;
        fluxes;
        errors;
        fluxes_no_outliers;
        fluxes_detrended;
        stds_detrended
        fluxes_global_cal;
        fluxes_global_cal_no_outliers;
        
        a_pars;
        b_pars;
        
        idx_outliers;
        
    end
    
    properties % switches/controls
        
        iterations = 2; % number of iterations of fit+outlier rejection
        num_sigma = 5; % for outlier rejection
        
        display_number = 5;
        display_what = 'raw'; % can also choose "cal"
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.02;
        
    end
    
    methods % constructor
        
        function obj = Calibrator(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.Calibrator')
                if obj.debug_bit, fprintf('Calibrator copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Calibrator constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            obj.timestamps = [];
            obj.fluxes = [];
            obj.errors = [];
            obj.fluxes_no_outliers = [];
            obj.fluxes_detrended = [];
            obj.stds_detrended = [];
            obj.fluxes_global_cal = [];
            obj.fluxes_global_cal_no_outliers = [];
            
            obj.a_pars = [];
            obj.b_pars = [];
            
            obj.idx_outliers = [];
            
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function input(obj, fluxes, varargin)
            
            if nargin<2
                help('trig.Calibrator.input'); return;
            end
            
            if isempty(fluxes)
                return;
            end
            
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('errors', []);
            input.input_var('timestamps', []);
            input.scan_vars(varargin{:});
            
            obj.clear;
            
            obj.fluxes = fluxes;
            obj.fluxes_no_outliers = fluxes;
            obj.errors = input.errors;
            
            if ~isempty(input.timestamps)
                obj.timestamps = input.timestamps;
            else
                obj.timestamps = (1:size(obj.fluxes,1))'/25;
            end
            
            obj.cleanup_fluxes;
            
            obj.global_calibration;
            
        end
        
        function cleanup_fluxes(obj)
            
            for ii = 1:obj.iterations
                
                obj.fit;
                obj.outlier_removal;
                
            end
            
        end
        
        function fit(obj)
            
            tt = tic;
            
            % some useful shorthands
            s = @(x) sum(x,1, 'omitnan');
            f = obj.fluxes_no_outliers;
            t = obj.timestamps;
            n = size(f,1);
            % need to adapt this to using the error estimates
            
            denom = n.*s(t.^2)-s(t).^2;
            
            obj.b_pars = (s(f).*s(t.^2)-s(t).*s(t.*f))./denom;
            obj.a_pars = (n.*s(t.*f)-s(t).*s(f))./denom;
            
            if obj.debug_bit>1, fprintf('runtime "fit": %f seconds\n', toc(tt)); end
            
        end
        
        function outlier_removal(obj)
            
            t = tic;
            
            obj.fluxes_detrended = obj.fluxes_no_outliers - obj.b_pars - obj.a_pars.*obj.timestamps;
            obj.fluxes_detrended = obj.fluxes_detrended - mean(obj.fluxes_detrended, 'omitnan');
            
            obj.stds_detrended = std(obj.fluxes_detrended, [], 1, 'omitnan');
            
            idx = abs(obj.fluxes_detrended./obj.stds_detrended)>obj.num_sigma;
            
            if isempty(obj.idx_outliers)
                obj.idx_outliers = idx;
            else
                obj.idx_outliers = obj.idx_outliers | idx;
            end
            
            obj.fluxes_no_outliers = obj.fluxes;
            obj.fluxes_no_outliers(obj.idx_outliers) = NaN;
            
            if obj.debug_bit>1, fprintf('runtime "outlier_removal": %f seconds\n', toc(t)); end
            
        end
        
        function global_calibration(obj)
            
            t = tic;
            
            F = obj.fluxes_no_outliers;
            V = var(obj.fluxes_no_outliers, [], 'omitnan');
            
            S = sum(F, 2, 'omitnan'); % sum of fluxe for each frame
%             S = sum(F.^2./V,2,'omitnan'); % sum of fluxe for each frame
            
            S_average = mean(S);
            
            T = S./S_average; % transparency for each image (relative to mean image)
            
            obj.fluxes_global_cal = obj.fluxes./T;
            
            obj.fluxes_global_cal_no_outliers = obj.fluxes_global_cal;
            obj.fluxes_global_cal_no_outliers(obj.idx_outliers) = NaN;
            
            if obj.debug_bit>1, fprintf('runtime "global_calibration": %f seconds\n', toc(t)); end

        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('number', obj.display_number);
            input.input_var('outliers', true);
            input.input_var('display_what', obj.display_what);
            input.scan_vars(varargin{:});
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            input.ax.NextPlot = 'replace'; 
            
            for ii = 1:input.number
                    
                if util.text.cs(input.display_what, 'raw')

                    str = sprintf('Star %d: f = %6.2f+%6.4ft', ii, obj.b_pars(ii), obj.a_pars(ii));

                    h = plot(input.ax, obj.timestamps, obj.fluxes_reduced(:,ii), '.');
                    h.DisplayName = str;
                    h.ButtonDownFcn = @obj.callback_touch;
                    h.HandleVisibility = 'off';

                    input.ax.NextPlot = 'add'; 

                    h2 = plot(input.ax, obj.timestamps, obj.a_pars(ii).*obj.timestamps + obj.b_pars(ii), '-', 'Color', h.Color);
                    h2.DisplayName = str;
                    h2.ButtonDownFcn = @obj.callback_touch;

                    if input.outliers
                        h3 = plot(input.ax, obj.timestamps(obj.idx_outliers(:,ii)), obj.fluxes(obj.idx_outliers(:,ii),ii), 'x', 'Color', h.Color);
                        if ~isempty(h3)
                            h3.DisplayName = str;
                            h3.ButtonDownFcn = @obj.callback_touch;
                            h3.HandleVisibility = 'off';
                        end
                    end
                    
                else
                    
                    cal_mean = mean(obj.fluxes_global_cal_no_outliers, 'omitnan');
                    cal_var = var(obj.fluxes_global_cal_no_outliers, [], 'omitnan');
                    
                    str = sprintf('Star %d: mean flux = %6.2f | var: %6.2f', ii, cal_mean(ii), cal_var(ii));

                    h = plot(input.ax, obj.timestamps, obj.fluxes_global_cal(:,ii), '-');
                    h.DisplayName = str;
                    h.ButtonDownFcn = @obj.callback_touch;
                    
                    input.ax.NextPlot = 'add'; 
                    
                end
   
            end
            
            input.ax.NextPlot = 'replace';
            xlabel(input.ax, 'timestamp (seconds)');
            ylabel(input.ax, 'flux (counts)');
            legend(input.ax);
            
        end
        
        function callback_touch(obj, hndl, event)
            
            str = sprintf('%s (%4.2f, %6.1f)', hndl.DisplayName, event.IntersectionPoint(1),event.IntersectionPoint(2));
            
            if obj.debug_bit, disp(['callback: touch. ' str]); end
            
            % add other GUI changes here...
            
        end
        
    end    
    
end

