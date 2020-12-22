classdef KuiperBeltModel < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        index_power_law = -3.8;
        index_lower = -4.0;
        index_upper = -3.6;
        
        start_radius = 0.25; % km
        
        normalization = 1.1e7; % number of objects above "start_radius"
        norm_lower = 0.4e7; % lower limit on normalization
        norm_upper = 2.6e7; % upper limit on normalization
        
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
        
        function obj = KuiperBeltModel(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.KuiperBeltModel')
                if obj.debug_bit>1, fprintf('KuiperBeltModel copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('KuiperBeltModel constructor v%4.2f\n', obj.version); end
            
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
        
        function [N, N_l, N_u] = numDensityIntervals(obj, r_edges, power_law_range) % how many objects per square degree in each radius bin
            
            if nargin<3 || isempty(power_law_range)
                power_law_range = 1;
            end
            
            q = abs(obj.index_power_law); 
            
            if power_law_range
                ql = abs(obj.index_lower);
                qu = abs(obj.index_upper); 
            else
                ql = q;
                qu = q;
            end
            
            N = (r_edges(1:end-1)./obj.start_radius).^(1-q) - (r_edges(2:end)./obj.start_radius).^(1-q);
            N = obj.normalization.*N';
            
            N_l = (r_edges(1:end-1)./obj.start_radius).^(1-qu) - (r_edges(2:end)./obj.start_radius).^(1-qu);
            N_l = obj.norm_lower.*N_l';
            
            N_u = (r_edges(1:end-1)./obj.start_radius).^(1-ql) - (r_edges(2:end)./obj.start_radius).^(1-ql);
            N_u = obj.norm_upper.*N_u';
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('r_edges', obj.start_radius:0.1:3); 
            input.input_var('power_law_range', true); 
            input.input_var('log', true, 'logarithm'); 
            input.input_var('axes', [], 'axis'); 
            input.input_var('font_size', 18); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            [N, N_l, N_u] = obj.numDensityIntervals(input.r_edges, input.power_law_range); 
            
            r = util.vec.tocolumn(input.r_edges);
            r = r(1:end-1) + diff(r)/2; 
            
            [h_line, h_shade] = util.plot.shaded(r, N, [N-N_l, N_u-N], 'axes', input.axes); 
            
            h_line.DisplayName = 'KBO model: N(>r)=N_0 r^{-q}';
            
            power_of_ten = floor(log10(obj.normalization)); 
            N_0 = obj.normalization./10.^power_of_ten;
            N_p = obj.norm_upper./10.^power_of_ten - N_0; 
            N_m = N_0 - obj.norm_lower./10.^power_of_ten; 
            
            if input.power_law_range
                h_shade.DisplayName = sprintf('N_0= %3.1f_{-%3.1f}^{+%3.1f}x10^{%d} | q= %3.1f_{-%3.1f}^{+%3.1f}', N_0, N_m, N_p, power_of_ten, ...
                    abs(obj.index_power_law), abs(obj.index_power_law) - abs(obj.index_lower), abs(obj.index_upper) - abs(obj.index_power_law)); 
            else
                h_shade.DisplayName = sprintf('%g < N_0 < %g', obj.norm_lower, obj.norm_upper); 
            end
            
            if input.log
                input.axes.YScale = 'log';
            else
                input.axes.YScale = 'linear';
            end
            
            xlabel(input.axes, 'Occulter radius [km]'); 
            ylabel(input.axes, 'Number density [deg^{-2}]'); 
            
            input.axes.FontSize = input.font_size; 
            
        end
        
    end    
    
end

