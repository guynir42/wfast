classdef Calibrator < handle

    properties(Transient=true)
        
        gui; 
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        timestamps;
        fluxes;
        fluxes_reduced;
        fluxes_subtracted;
        fluxes_no_outliers;
        fluxes_calibrated;
        fluxes_cal_no_outliers;
        
        a_pars;
        b_pars;
        
        idx_outliers;
        sub_var;
        out_var;
        out_mean;
        cal_mean;
        cal_var;
        
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
       
        version = 1.00;
        
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
            obj.fluxes_reduced = [];
            obj.fluxes_subtracted = [];
            obj.fluxes_no_outliers = [];
            obj.fluxes_calibrated = [];
            obj.fluxes_cal_no_outliers = [];
            
            obj.a_pars = [];
            obj.b_pars = [];
            
            obj.idx_outliers = [];
            obj.sub_var = [];
            obj.out_mean = [];
            obj.out_var = [];
            obj.cal_mean = [];
            obj.cal_var = [];
            
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
            input.input_var('timestamps', []);
            input.scan_vars(varargin{:});
            
            obj.clear;
            
            obj.fluxes = fluxes;
            obj.fluxes_reduced = obj.fluxes;
            
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
            
            obj.fluxes_no_outliers = obj.fluxes;
            obj.fluxes_no_outliers(obj.idx_outliers) = NaN;
            obj.out_mean = mean(obj.fluxes_no_outliers, 1, 'omitnan');
            obj.out_var = var(obj.fluxes_no_outliers, [], 1, 'omitnan');
            
        end
        
        function fit(obj)
            
            % some useful shorthands
            s = @(x) sum(x,1, 'omitnan');
            f = obj.fluxes_reduced;
            t = obj.timestamps;
            n = size(f,1);
            
            denom = n.*s(t.^2)-s(t).^2;
            
            obj.b_pars = (s(f).*s(t.^2)-s(t).*s(t.*f))./denom;
            obj.a_pars = (n.*s(t.*f)-s(t).*s(f))./denom;
            
        end
        
        function outlier_removal(obj)
            
            obj.fluxes_subtracted = obj.fluxes_reduced - obj.b_pars - obj.a_pars.*obj.timestamps;
            
            obj.sub_var = var(obj.fluxes_subtracted, [], 1, 'omitnan');
            
            idx = abs(obj.fluxes_subtracted./sqrt(obj.sub_var))>obj.num_sigma;
            
            if isempty(obj.idx_outliers)
                obj.idx_outliers = idx;
            else
                obj.idx_outliers = obj.idx_outliers | idx;
            end
            
            obj.fluxes_reduced(obj.idx_outliers) = NaN;
            
        end
        
        function global_calibration(obj)
            
            F = obj.fluxes_no_outliers;
            V = obj.out_var;
            
            S = sum(F, 2, 'omitnan'); % sum of fluxe for each frame
%             S = sum(F.^2./V,2,'omitnan'); % sum of fluxe for each frame
            
            S_average = mean(S);
            
            T = S./S_average; % transparency for each image (relative to mean image)
            
            obj.fluxes_calibrated = obj.fluxes./T;
            
            obj.fluxes_cal_no_outliers = obj.fluxes_calibrated;
            obj.fluxes_cal_no_outliers(obj.idx_outliers) = NaN;
            
            obj.cal_mean = mean(obj.fluxes_cal_no_outliers, 1, 'omitnan');
            obj.cal_var = var(obj.fluxes_cal_no_outliers, [], 1, 'omitnan');
            
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
                    
                    str = sprintf('Star %d: mean flux = %6.2f | var: %6.2f', ii, obj.cal_mean(ii), obj.cal_var(ii));

                    h = plot(input.ax, obj.timestamps, obj.fluxes_calibrated(:,ii), '-');
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

