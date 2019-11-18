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
        
        num_successes = 0;
        num_failures = 0;
        
        runtime = 0;
        
    end
    
    properties % switches/controls
        
        num_chains = 1;
        num_steps = 1000; 
        num_burned = 100; 
        
        step_size = 0.1;
        
        default_error = 0.1; % if input lightcurves don't have any errors per sample, use this instead... 
        
        par_list = {'r', 'b', 'v'}; % these parameters are chosen randomly each step. The rest are taken from the generator's parameters
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)

        version = 1.01;
        
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

            obj.num_successes = 0;
            obj.num_failures = 0;
            
            obj.runtime = 0;
            
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function useSimulatedInput(obj, varargin)
            
            obj.gen.getLightCurves;
            obj.input_flux = obj.gen.lc.flux;
            
            
        end
        
        function startup(obj)
            
            if length(obj.gen.r)>1
                error('Must setup generator with all scalar parameters!');
            end
            
            obj.points(obj.num_steps, obj.num_chains) = occult.Parameters;
            
            obj.runtime = 0;
            
        end
        
        function gotoInitialPoint(obj)
            
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
        
        function takeStep(obj)
            
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
            obj.gen.lc.pars.likelihood = chi2cdf(chi2, nnz(~isnan(f1))-length(obj.par_list), 'upper'); % do we need to subtract these parameters from the dof??
            
        end
        
        function run(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('plot', true);
            input.input_var('ax', [], 'axes', 'axis');
            input.scan_vars(varargin{:});
            
            obj.reset;
            obj.startup;
            
            if input.plot
                
                if isempty(input.ax)
                    input.ax = gca;
                end
                
            end
            
            t = tic;
            
            obj.prog.start(obj.num_chains.*obj.num_steps);
            
            for ii = 1:obj.num_chains
                
                obj.gotoInitialPoint;
                obj.gen.getLightCurves;
                obj.calcLikelihood; % get the first point likelihood/chi2
                
                obj.points(1,ii).copy_from(obj.gen.lc.pars); % automatically accept first point... 
                
                for jj = 2:obj.num_steps
                    
                    obj.prog.showif((ii-1)*obj.num_steps + jj);
                    
                    obj.takeStep; % move the generator to a nearby random point
                    
                    obj.gen.getLightCurves;
                    obj.calcLikelihood;
                    
                    ratio = obj.gen.lc.pars.likelihood./obj.points(jj-1,ii).likelihood;
                    
                    if rand<=ratio % if ratio is big, we are likely to take the new point
                        obj.num_successes = obj.num_successes + 1;
                        obj.points(jj,ii).copy_from(obj.gen.lc.pars);
                        obj.points(jj,ii).counts = 1;
                    else % point is rejected, repeat the previous point
                        obj.num_failures = obj.num_failures + 1;
                        obj.points(jj,ii).copy_from(obj.points(jj-1,ii));
                        obj.points(jj,ii).counts = obj.points(jj-1,ii).counts + 1;
                    end
                    
                end
                
            end
            
            obj.runtime = toc(t);
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function plot(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('pars', []); 
            input.scan_vars(varargin{:});
                        
            if isempty(input.ax)
                input.ax = gca;
            end
            
            if isempty(input.pars)
                input.pars = obj.par_list(1:min(3,length(obj.par_list)));
            end
            
            if ischar(input.pars)
                input.pars = {input.pars};
            end
            
            p1 = [obj.points.(input.pars{1})];
            p1 = reshape(p1, [obj.num_steps, obj.num_chains]);
            s = [obj.points.counts];
            s = reshape(s, [obj.num_steps, obj.num_chains]);
            
            if length(input.pars)==1
                histogram(input.ax, p1);
                
            elseif length(input.pars)==2
                
                p2 = [obj.points.(input.pars{2})];
                p2 = reshape(p2, [obj.num_steps, obj.num_chains]);
                p2 = p2(obj.num_burned:end,:);
                
                plot(input.ax, p1, p2, '.');
                
                ylabel(input.ax, input.pars{2});
                input.ax.YLim = obj.gen.([input.pars{2} '_range']);
                
            elseif length(input.pars)==3
                
                p2 = [obj.points.(input.pars{2})];
                p2 = reshape(p2, [obj.num_steps, obj.num_chains]);
                
                p3 = [obj.points.(input.pars{3})];
                p3 = reshape(p3, [obj.num_steps, obj.num_chains]);
                
                holding_pattern = input.ax.NextPlot;
                
                scatter3(input.ax, p1(obj.num_burned:end,1), p2(obj.num_burned:end,1), p3(obj.num_burned:end,1), 20+s(obj.num_burned:end,1), '.');
                
                input.ax.NextPlot = 'add';
                
                for ii = 2:size(p1,2)
                    scatter3(input.ax, p1(obj.num_burned:end,ii), p2(obj.num_burned:end,ii), p3(obj.num_burned:end,ii), 20+s(obj.num_burned:end,ii), '.');
                end
                
                input.ax.NextPlot = holding_pattern;
                
                ylabel(input.ax, input.pars{2});
                zlabel(input.ax, input.pars{3});
                
                input.ax.YLim = obj.gen.([input.pars{2} '_range']);
                input.ax.ZLim = obj.gen.([input.pars{3} '_range']);
                
            end
        
            xlabel(input.ax, input.pars{1});
            input.ax.FontSize = 24;
            
            input.ax.XLim = obj.gen.([input.pars{1} '_range']);
            
        end
        
    end    
    
end

