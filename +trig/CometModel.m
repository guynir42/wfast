classdef CometModel < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        name = 'KBOs'; % by default, create a KBO model
        
        index_power_law = 3.8;
        index_lower = 3.6;
        index_upper = 4.0;
        
        start_radius = 0.25; % km
        
        normalization = 1.1e7; % number of objects above "start_radius"
        norm_lower = 0.4e7; % lower limit on normalization
        norm_upper = 2.6e7; % upper limit on normalization
        
        distance = 40; % AU
        
        
    end
    
    properties % switches/controls
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        wavelength = 550; 
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = CometModel(varargin)
            
            import util.text.cs;
            
            if ~isempty(varargin) && isa(varargin{1}, 'trig.CometModel')
                if obj.debug_bit>1, fprintf('CometModel copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            elseif ~isempty(varargin) && ischar(varargin{1})
                
                if cs(varargin{1}, 'kbos', 'kuiper belt model')
                    % leave all defaults
                elseif cs(varargin{1}, 'hills cloud', 'inner oort cloud')
                    obj.name = 'Hills'; 
                    obj.start_radius = 5; 
                    obj.normalization = 1e12/(4*180^2/pi); % total number of comets, over entire sky
                    obj.norm_lower = obj.normalization/10; 
                    obj.norm_upper = obj.normalization*10; 
                    obj.distance = 3000;
                elseif cs(varargin{1}, 'oort cloud')
                    obj.name = 'Oort'; 
                    obj.start_radius = 5; 
                    obj.normalization = 1e12/(4*180^2/pi); % total number of comets, over entire sky
                    obj.norm_lower = obj.normalization/10; 
                    obj.norm_upper = obj.normalization*10; 
                    obj.distance = 10000; 
                else
                    error('Unknown comet model: "%s". Use "KBOs" or "Hills" or "Oort"', varargin{1}); 
                end
                
            else
                if obj.debug_bit>1, fprintf('CometModel constructor v%4.2f\n', obj.version); end
            
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
        
        function [N, N_l, N_u] = numDensityCumulative(obj, r_edges, power_law_range) % how many objects per square degree larger or equal to that radius
            
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
            
            N = (r_edges(1:end-1)./obj.start_radius).^(1-q);
            N = obj.normalization.*N';
            
            N_l = (r_edges(1:end-1)./obj.start_radius).^(1-qu);
            N_l = obj.norm_lower.*N_l';
            
            N_u = (r_edges(1:end-1)./obj.start_radius).^(1-ql);
            N_u = obj.norm_upper.*N_u';
            
        end
        
        function [N, N_l, N_u] = numDensityIntervals(obj, r_edges, power_law_range) % how many objects per square degree in each radius bin
            
            if nargin<3 || isempty(power_law_range)
                power_law_range = 1;
            end
            
            q = obj.index_power_law; 
            
            if power_law_range
                ql = obj.index_lower;
                qu = obj.index_upper; 
            else
                ql = q;
                qu = q;
            end
            
            N = (r_edges(1:end-1)./obj.start_radius).^(1-q) - (r_edges(2:end)./obj.start_radius).^(1-q);
            N = obj.normalization.*N';
            
            N_l = (r_edges(1:end-1)./obj.start_radius).^(1-ql) - (r_edges(2:end)./obj.start_radius).^(1-ql);
            N_l = obj.norm_lower.*N_l';
            
            N_u = (r_edges(1:end-1)./obj.start_radius).^(1-qu) - (r_edges(2:end)./obj.start_radius).^(1-qu);
            N_u = obj.norm_upper.*N_u';
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('r_edges', []); 
            input.input_var('power_law_range', true); 
            input.input_var('log', true, 'logarithm'); 
            input.input_var('axes', [], 'axis'); 
            input.input_var('font_size', 18); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            if isempty(input.r_edges)
                fsu2km = sqrt(obj.wavelength.*1e-12.*obj.distance.*1.496e+8/2); 
                r_km = (0.3:0.1:3)*fsu2km;
                input.r_edges = r_km;
            end
            
            r = util.vec.tocolumn(input.r_edges);
            r = r(1:end-1) + diff(r)/2; 
            
            [N, N_l, N_u] = obj.numDensityCumulative(input.r_edges, input.power_law_range); 
            
            [h_line, h_shade] = util.plot.shaded(r, N, [N-N_l, N_u-N], 'axes', input.axes); 
            
            h_line.DisplayName = 'KBO model: N(>r)=N_0 r^{-q}';
            
            power_of_ten = floor(log10(obj.normalization)); 
            N_0 = obj.normalization./10.^power_of_ten;
            N_p = obj.norm_upper./10.^power_of_ten - N_0; 
            N_m = N_0 - obj.norm_lower./10.^power_of_ten; 
            
            if input.power_law_range
                h_shade.DisplayName = sprintf('N_0= %3.1f_{-%3.1f}^{+%3.1f}x10^{%d} | q= %3.1f_{-%3.1f}^{+%3.1f}', N_0, N_m, N_p, power_of_ten, ...
                    obj.index_power_law, obj.index_power_law - obj.index_lower, obj.index_upper - obj.index_power_law); 
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

