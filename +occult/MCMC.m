classdef MCMC < handle

    properties(Transient=true)
        
        gui;
        
    end
    
    properties % objects
        
        gen@occult.CurveGenerator;
        
        points@occult.Parameters; % a 2D matrix, each column contains one chain of points in parameter space
        
        prog@util.sys.ProgressBar;
        
    end
    
    properties % inputs/outputs
        
        input_flux; % the input lightcurve to be matched
        input_errors; % the rms error on the given flux
        
        runtime = 0;
        
    end
    
    properties % switches/controls
        
        num_chains = 1;
        num_steps = 1000; 
        num_burned = 100; 
        
        step_size = 0.01;
        
        default_error = 0.1; % if input lightcurves don't have any errors per sample, use this instead... 
        
        par_list = {'r', 'b', 'v'}; % these parameters are chosen randomly each step. The rest are taken from the generator's parameters
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = MCMC(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'occult.,MCMC')
                if obj.debug_bit, fprintf('MCMC copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('MCMC constructor v%4.2f\n', obj.version); end
            
                obj.gen = occult.CurveGenerator;
                obj.prog = util.sys.ProgressBar;
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.points = occult.Parameters.empty;
            
            obj.runtime = 0;
            
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function startup(obj)
            
            if length(obj.gen.r)>1
                error('Must setup generator with all scalar parameters!');
            end
            
            obj.points(obj.num_steps, obj.num_chains) = occult.Parameters;
            
            obj.runtime = 0;
            
        end
        
        function initial_point(obj)
            
            for ii = 1:length(obj.par_list)
                
                if strcmp(obj.par_list{ii}, 'r')
                    obj.gen.r = (obj.gen.max_r-obj.gen.min_r)*rand + obj.gen.min_r;
                elseif strcmp(obj.par_list{ii}, 'r2')
                    obj.gen.r2 = (obj.gen.max_r2-obj.gen.min_r2)*rand + obj.gen.min_r2;
                elseif strcmp(obj.par_list{ii}, 'd')
                    obj.gen.d = (obj.gen.max_d-obj.gen.min_d)*rand + obj.gen.min_d;
                elseif strcmp(obj.par_list{ii}, 'th')
                    obj.gen.th = (obj.gen.max_th-obj.gen.min_th)*rand + obj.gen.min_th;
                elseif strcmp(obj.par_list{ii}, 'R')
                    obj.gen.R = (obj.gen.max_R-obj.gen.min_R)*rand + obj.gen.min_R;
                elseif strcmp(obj.par_list{ii}, 'b')
                    obj.gen.b = (obj.gen.max_b-obj.gen.min_b)*rand + obj.gen.min_b;
                elseif strcmp(obj.par_list{ii}, 'v')
                    obj.gen.v = (obj.gen.max_v-obj.gen.min_v)*rand + obj.gen.min_v;
                elseif strcmp(obj.par_list{ii}, 't')
                    obj.gen.t = (obj.gen.max_t-obj.gen.min_t)*rand + obj.gen.min_t;
                end
                
            end
                
        end
        
        function step(obj)
            
            for ii = 1:length(obj.par_list)
                
                if strcmp(obj.par_list{ii}, 'r')
                    
                    obj.gen.r = normrnd(obj.gen.r, obj.step_size);
                    if obj.gen.r<obj.gen.min_r
                        obj.gen.r = obj.gen.min_r;
                    elseif obj.gen.r>obj.gen.max_r
                        obj.gen.r = obj.gen.max_r;
                    end
                    
                elseif strcmp(obj.par_list{ii}, 'r2')
                    
                    obj.gen.r2 = normrnd(obj.gen.r2, obj.step_size);
                    if obj.gen.r2<obj.gen.min_r2
                        obj.gen.r2 = obj.gen.min_r2;
                    elseif obj.gen.r2>obj.gen.max_r2
                        obj.gen.r2 = obj.gen.max_r2;
                    end
                    
                elseif strcmp(obj.par_list{ii}, 'd')
                    
                    obj.gen.d = normrnd(obj.gen.d, obj.step_size);
                    if obj.gen.d<obj.gen.min_d
                        obj.gen.d = obj.gen.min_d;
                    elseif obj.gen.d>obj.gen.max_d
                        obj.gen.d = obj.gen.max_d;
                    end
                    
                elseif strcmp(obj.par_list{ii}, 'th')
                    
                    obj.gen.th = normrnd(obj.gen.th, 1); % note angles step is taken as 1 deg for simplicity
                    if obj.gen.th<obj.gen.min_th
                        obj.gen.th = obj.gen.min_th;
                    elseif obj.gen.th>obj.gen.max_th
                        obj.gen.th = obj.gen.max_th;
                    end
                    
                elseif strcmp(obj.par_list{ii}, 'R')
                    
                    obj.gen.R = normrnd(obj.gen.R, obj.step_size);
                    if obj.gen.R<obj.gen.min_R
                        obj.gen.R = obj.gen.min_R;
                    elseif obj.gen.R>obj.gen.max_R
                        obj.gen.R = obj.gen.max_R;
                    end
                    
                elseif strcmp(obj.par_list{ii}, 'b')
                    
                    obj.gen.b = normrnd(obj.gen.b, obj.step_size);
                    if obj.gen.b<obj.gen.min_b
                        obj.gen.b = obj.gen.min_b;
                    elseif obj.gen.b>obj.gen.max_b
                        obj.gen.b = obj.gen.max_b;
                    end
                    
                elseif strcmp(obj.par_list{ii}, 'v')
                    
                    obj.gen.v = normrnd(obj.gen.v, 1); % note the velocity step is 1 FSU/s for simplicity
                    if obj.gen.v<obj.gen.min_v
                        obj.gen.v = obj.gen.min_v;
                    elseif obj.gen.v>obj.gen.max_v
                        obj.gen.v = obj.gen.max_v;
                    end
                    
                elseif strcmp(obj.par_list{ii}, 't')
                    
                    obj.gen.t = normrnd(obj.gen.t, 1); % note the time offset step is 1 ms for simplicity
                    if obj.gen.t<obj.gen.min_t
                        obj.gen.t = obj.gen.min_t;
                    elseif obj.gen.t>obj.gen.max_t
                        obj.gen.t = obj.gen.max_t;
                    end
                    
                end
                
            end
            
        end
        
        function calcLikelihood(obj)
            
            if ~isempty(obj.input_errors)
                e = obj.input_errors;
            else
                e = obj.default_error;
            end
            
            f1 = obj.gen.lc.flux;
            f2 = obj.input_flux;
            
            chi2 = nansum(((f1-f2)./e).^2); 
            
            obj.gen.lc.pars.chi2 = chi2;
            obj.gen.lc.pars.likelihood = chi2cdf(chi2, length(f1)-length(obj.par_list), 'upper'); % do we need to subtract these parameters from the dof??
            
        end
        
        function run(obj)
            
            obj.reset;
            obj.startup;
            
            t = tic;
            
            obj.prog.start(obj.num_chains.*obj.num_steps);
            
            for ii = 1:obj.num_chains
                
                obj.initial_point;
                obj.gen.getLightCurves;                
                obj.calcLikelihood; % get the first point likelihood/chi2
                
                obj.points(1,ii).copy_from(obj.gen.lc.pars); % automatically accept first point... 
                
                for jj = 2:obj.num_steps
                    
                    obj.prog.showif((ii-1)*obj.num_steps + jj);
                    
                    obj.step; % move the generator to a nearby random point
                    
                    obj.gen.getLightCurves;
                    obj.calcLikelihood;
                    
                    ratio = obj.gen.lc.pars.likelihood./obj.points(jj-1,ii).likelihood;
                    
                    if rand<=ratio % if ratio is big, we are likely to take the new point
                        obj.points(jj,ii).copy_from(obj.gen.lc.pars);
                    else % point is rejected, repeat the previous point
                        obj.points(jj,ii).copy_from(obj.points(jj-1,ii));
                    end
                    
                end
                
            end
            
            obj.runtime = toc(t);
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function plot(obj, varargin)
            
            
            
        end
        
    end    
    
end

