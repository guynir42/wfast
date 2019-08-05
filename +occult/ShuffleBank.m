classdef ShuffleBank < handle

    properties(Transient=true)
        
        gen@occult.CurveGenerator;
        prog@util.sys.ProgressBar;
        
    end
    
    properties % objects
        
        
    end
    
    properties % inputs/outputs
        
        % outputs:
        time_axis; % a 1D time axis for all LCs
        kernels; % a 2D map, with kernels in 1st dimension, and parameters in 2nd
        pars; % a struct vector, one struct per lightcurve, with R,r,b,v values
        
        % input/output data
        fluxes; 
        fluxes_filtered;
        stds;
        timestamps;
        
        snrs_tested; 
        
    end
    
    properties % switches/controls
        
        R_range = [0 0.5]; % star radius, FSU
        r_range = [0.3 2]; % occulter radius, FSU
        b_range = [0 2]; % impact parameter, FSU
        v_range = [5 30]; % crossing velocity (projected), FSU/second
        % lets assume t=0 for all!
        
        W = 4; % time window, seconds
        T = 30; % integration time, ms
        f = 25; % frame rate, Hz
        
        number = 1e6;
        threshold = 0.95; % if other kernel loses information below this limit, it needs its own filter
        test_number = 5000; 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = ShuffleBank(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'occult.ShuffleBank')
                if obj.debug_bit, fprintf('ShuffleBank copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('ShuffleBank constructor v%4.2f\n', obj.version); end
            
                obj.gen = occult.CurveGenerator;
                obj.prog = util.sys.ProgressBar;
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function clear(obj)
            
            obj.fluxes = [];
            obj.fluxes_filtered = [];
            obj.stds = [];
            obj.timestamps = [];
        
        end
        
    end
    
    methods % getters
        
        function val = get.gen(obj)
            
            if isempty(obj.gen)
                obj.gen = occult.CurveGenerator;
            end
            
            val = obj.gen;
            
        end
        
        function val = get.prog(obj)
            
            if isempty(obj.prog)
                obj.prog = util.sys.ProgressBar;
            end
            
            val = obj.prog;
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function makeKernels(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('number', obj.number, 'N');
            input.input_var('threshold', obj.threshold);
            input.scan_vars(varargin{:});
            
            obj.kernels = [];
            obj.pars = struct([]);
            
            obj.prog.start(input.number);
            
            counter = 0; 
            
            for ii = 1:input.number
                
                obj.gen.R = obj.R_range(1) + rand.*(obj.R_range(2)-obj.R_range(1));
                obj.gen.r = obj.r_range(1) + rand.*(obj.r_range(2)-obj.r_range(1));
                obj.gen.b = obj.b_range(1) + rand.*(obj.b_range(2)-obj.b_range(1));
                obj.gen.v = obj.v_range(1) + rand.*(obj.v_range(2)-obj.v_range(1));
                obj.gen.getLightCurves;
                
                if isempty(obj.kernels)
                    obj.kernels = single(obj.gen.lc.flux - 1);
                    obj.pars = struct('R', obj.gen.R, 'r', obj.gen.r, 'b', obj.gen.b, 'v', obj.gen.v, 'norm', sqrt(sum(obj.kernels(:,1).^2)));
                    obj.kernels = obj.kernels./obj.pars.norm;
                    obj.time_axis = obj.gen.lc.time;
                else
                    snr = occult.compareKernels(single(obj.gen.lc.flux-1), obj.kernels);
                    counter = counter + 1;
                    
                    if max(snr)<obj.threshold
                        if obj.debug_bit>1, fprintf('Adding new kernel with R= %4.2f | r= %4.2f | b= %4.2f | v= %5.3f\n', obj.gen.R, obj.gen.r, obj.gen.b, obj.gen.v); end
                        obj.kernels = horzcat(obj.kernels, single(obj.gen.lc.flux - 1));
                        obj.pars = horzcat(obj.pars, struct('R', obj.gen.R, 'r', obj.gen.r, 'b', obj.gen.b, 'v', obj.gen.v, 'norm', sqrt(sum(obj.kernels(:,end).^2)))); 
                        obj.kernels(:,end) = obj.kernels(:,end)./obj.pars(end).norm;
                        counter = 0; % popcorn method
                    end
                end
                
                if counter>obj.test_number % popcorn method
                    break;
                end
                
                obj.prog.showif(ii);
                
            end
            
            fprintf('Finished shuffling filter kernels. Iterations= %d | pop count= %d | num filters= %d\n', ii, counter, size(obj.kernels,2));
            
        end
        
        function runTest(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('number', obj.test_number, 'N');
            input.input_var('threshold', obj.threshold);
            input.scan_vars(varargin{:});
            
            if isempty(obj.kernels)
                error('Need to first fill the filter kernels before running a test. Try using "makeKernels"...');
            end
            
            obj.snrs_tested = NaN(input.number,1);
            
            obj.prog.start(input.number);
            
            for ii = 1:input.number
                
                obj.gen.R = obj.R_range(1) + rand.*(obj.R_range(2)-obj.R_range(1));
                obj.gen.r = obj.r_range(1) + rand.*(obj.r_range(2)-obj.r_range(1));
                obj.gen.b = obj.b_range(1) + rand.*(obj.b_range(2)-obj.b_range(1));
                obj.gen.v = obj.v_range(1) + rand.*(obj.v_range(2)-obj.v_range(1));
                obj.gen.getLightCurves;
                
                snr = occult.compareKernels(single(obj.gen.lc.flux-1), obj.kernels);
                
                obj.snrs_tested(ii) = max(snr);
                
                if max(snr)<obj.threshold
                    if obj.debug_bit, fprintf('Kernel test failed with S/N= %4.2f | R= %4.2f | r= %4.2f | b= %4.2f | v= %5.3f\n', max(snr), obj.gen.R, obj.gen.r, obj.gen.b, obj.gen.v); end
                end

                obj.prog.showif(ii);
                
            end
            
        end
        
        function input(obj, varargin)
            
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('flux', [], 'fluxes');
            input.input_var('stds', []);
            input.input_var('times', [], 'timestamps');
            input.scan_vars(varargin{:});
            
            if isempty(input.flux)
                error('Must input "flux" to filter...');
            end
            
            if size(input.flux,2)>1
                obj.fluxes = permute(input.flux, [1,3,2]);
            else
                obj.fluxes = input.flux;
            end
            
            obj.fluxes = fillmissing(obj.fluxes, 'spline'); 
            
            if size(input.stds,2)>1
                obj.stds = permute(input.stds, [1,3,2]);
            else
                obj.stds = input.stds;
            end
            
            obj.timestamps = input.times;
            
            obj.fluxes_filtered = util.vec.convolution(obj.kernels, obj.fluxes)./obj.stds;
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

