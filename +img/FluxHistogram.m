classdef FluxHistogram < handle

    properties % inputs/outputs
        
        batch_number = 1; % index of current batch
        juldates; % julian date of each batch
        airmass; % airmass of each batch
        flux_frac; % 2D matrix, dim 1 is batch/juldate/airmass, dim2 is number of epochs in each fractional-flux deviation interval
        edges_linear; % edges of the flux histogram
        flux_log; % 2D matrix, dim 1 is batch/juldate/airmass, dim2 is number of epochs in each log-flux deviation interval
        edges_log; % edges of the log-flux histogram
        
    end
    
    properties % switches/controls
        
        min_flux = 1000; % must have this much (batch average) flux to be counted
        
        binning_factors = [1, 3, 5]; 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        bins_linear;
        bins_log;
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = FluxHistogram(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'img.FluxHistogram')
                if obj.debug_bit>1, fprintf('FluxHistogram copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('FluxHistogram constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.batch_number = 1;
            obj.juldates = [];
            obj.airmass = [];
            obj.flux_frac = [];
            obj.edges_linear = [];
            obj.flux_log = [];
            obj.edges_log = [];
            
        end
        
    end
    
    methods % getters
        
        function val = getLinearEdges(obj)
            
            val = -1:0.01:1;
            
        end
        
        function val = getLogEdges(obj)
            
            val = -0.1:0.001:0.1;
            
        end
        
        function val = get.bins_linear(obj)
            
            val = obj.edges_linear(1:end-1)+diff(obj.edges_linear);
            
        end
        
        function val = get.bins_log(obj)
            
            val = obj.edges_log(1:end-1)+diff(obj.edges_log);
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function startup(obj, num_batches)
            
            if nargin > 1 && ~isempty(num_batches) % preallocate
            
                if isempty(obj.edges_linear)
                    obj.edges_linear = obj.getLinearEdges;
                end

                obj.flux_frac = NaN(num_batches, length(obj.edges_linear) - 1, length(obj.binning_factors));

                if isempty(obj.edges_log)
                    obj.edges_log = obj.getLogEdges;
                end

                obj.flux_log = NaN(num_batches, length(obj.edges_log) - 1, length(obj.binning_factors));

            end
            
        end
        
        function finishup(obj)
            
            obj.flux_frac = obj.flux_frac(1:obj.batch_number - 1,:,:); 
            obj.flux_log = obj.flux_log(1:obj.batch_number - 1,:,:); 
            
        end
        
        function input(obj, flux_corrected, juldates, airmass)
            
            F = flux_corrected;
            J = nanmean(juldates);
            A = nanmean(airmass); 
            
            mean_F = nanmean(F, 1); 
            idx = mean_F > obj.min_flux;
            
            F = F(:,idx); 
            mean_F = mean_F(idx); 
            
            if isempty(obj.edges_linear)
                obj.edges_linear = obj.getLinearEdges;
            end
            
            if isempty(obj.edges_log)
                obj.edges_log = obj.getLogEdges;
            end
            
            for ii = 1:length(obj.binning_factors)
                    
                FF = F; % filtered flux
                
                b = obj.binning_factors(ii);
                if b > 1
                    FF = filter2(ones(b,1)/b, F); 
                end
                
                flux_dev = (FF - mean_F)./mean_F; % flux deviation
                
                N = histcounts(flux_dev(:), obj.edges_linear);

                obj.flux_frac(obj.batch_number, :, ii) = N; 

                flux_log = log10(FF./mean_F);
                N = histcounts(flux_log(:), obj.edges_log);

                obj.flux_log(obj.batch_number, :, ii) = N; 

                obj.juldates(obj.batch_number, 1) = J;
                obj.airmass(obj.batch_number, 1) = A;
            
            end
            
            obj.batch_number = obj.batch_number + 1;
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

