classdef ShuffleBank < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        gen@occult.CurveGenerator;
        prog@util.sys.ProgressBar;
        
    end
    
    properties % inputs/outputs
        
        R_range = [0 0.5]; % star radius, FSU
        r_range = [0.3 2]; % occulter radius, FSU
        b_range = [0 2]; % impact parameter, FSU
        v_range = [5 30]; % crossing velocity (projected), FSU/second
        % lets assume t=0 for all!
        
        W = 4; % time window, seconds
        T = 30; % integration time, ms
        f = 25; % frame rate, Hz
        
        % outputs:
        timestamps; % a 1D time axis for all LCs
        bank; % a 2D map, with lightcurves in 1st dimension, and parameters in 2nd
        pars; % a struct vector, one struct per lightcurve, with R,r,b,v values
        
        snrs_tested; 
        
    end
    
    properties % switches/controls
        
        number = 1000;
        threshold = 0.95; % if other kernel loses information below this limit, it needs its own filter
        
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
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function run(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('number', obj.number, 'N');
            input.input_var('threshold', obj.threshold);
            input.scan_vars(varargin{:});
            
            obj.bank = [];
            obj.pars = struct([]);
            
            obj.prog.start(input.number);
            
            for ii = 1:input.number
                
                obj.gen.R = obj.R_range(1) + rand.*(obj.R_range(2)-obj.R_range(1));
                obj.gen.r = obj.r_range(1) + rand.*(obj.r_range(2)-obj.r_range(1));
                obj.gen.b = obj.b_range(1) + rand.*(obj.b_range(2)-obj.b_range(1));
                obj.gen.v = obj.v_range(1) + rand.*(obj.v_range(2)-obj.v_range(1));
                obj.gen.getLightCurves;
                
                if isempty(obj.bank)
                    obj.bank = obj.gen.lc.flux;
                    obj.timestamps = obj.gen.lc.time;
                else
                    snr = occult.compareKernels(obj.gen.lc.flux, obj.bank);
                    if max(snr)<obj.threshold
                        if obj.debug_bit>1, fprintf('Adding new kernel with R= %4.2f | r= %4.2f | b= %4.2f | v= %5.3f\n', obj.gen.R, obj.gen.r, obj.gen.b, obj.gen.v); end
                        obj.bank = horzcat(obj.bank, obj.gen.lc.flux);
                        obj.pars = horzcat(obj.pars, struct('R', obj.gen.R, 'r', obj.gen.r, 'b', obj.gen.b, 'v', obj.gen.v)); 
                    end
                end
                
                obj.prog.showif(ii);
                
            end
            
        end
        
        function test(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('number', obj.number, 'N');
            input.input_var('threshold', obj.threshold);
            input.scan_vars(varargin{:});
            
            if isempty(obj.bank)
                error('Need to first fill the filter bank before running a test. Try using "run"...');
            end
            
            obj.snrs_tested = NaN(input.number,1);
            
            obj.prog.start(input.number);
            
            for ii = 1:input.number
                
                obj.gen.R = obj.R_range(1) + rand.*(obj.R_range(2)-obj.R_range(1));
                obj.gen.r = obj.r_range(1) + rand.*(obj.r_range(2)-obj.r_range(1));
                obj.gen.b = obj.b_range(1) + rand.*(obj.b_range(2)-obj.b_range(1));
                obj.gen.v = obj.v_range(1) + rand.*(obj.v_range(2)-obj.v_range(1));
                obj.gen.getLightCurves;
                
                snr = occult.compareKernels(obj.gen.lc.flux, obj.bank);
                
                obj.snrs_tested(ii) = max(snr);
                
                if max(snr)<obj.threshold
                    if obj.debug_bit, fprintf('Kernel test failed with S/N= %4.2f | R= %4.2f | r= %4.2f | b= %4.2f | v= %5.3f\n', max(snr), obj.gen.R, obj.gen.r, obj.gen.b, obj.gen.v); end
                end

                obj.prog.showif(ii);
                
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

