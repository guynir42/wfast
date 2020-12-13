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
            
            q = abs(obj.index_power_law); 
            
            if input.power_law_range
                ql = abs(obj.index_lower);
                qu = abs(obj.index_upper); 
            else
                ql = q;
                qu = q;
            end
            
            N = (input.r_edges(1:end-1)./obj.start_radius).^(1-q) - (input.r_edges(2:end)./obj.start_radius).^(1-q);
            N = obj.normalization.*N';
            
            N_l = (input.r_edges(1:end-1)./obj.start_radius).^(1-qu) - (input.r_edges(2:end)./obj.start_radius).^(1-qu);
            N_l = obj.norm_lower.*N_l';
            
            N_u = (input.r_edges(1:end-1)./obj.start_radius).^(1-ql) - (input.r_edges(2:end)./obj.start_radius).^(1-ql);
            N_u = obj.norm_upper.*N_u';
                        
            r = util.vec.tocolumn(input.r_edges);
            r = r(1:end-1) + diff(r)/2; 
            
            util.plot.shaded(r, N, [N-N_l, N_u-N], 'axes', input.axes); 
            
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

