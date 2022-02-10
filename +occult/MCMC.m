classdef MCMC < handle

    properties(Transient=true)
        
        gui;
        gen@occult.CurveGenerator; % used to create lightcurves to fit to the measured
        bank@occult.ShuffleBank; % for finding an initial guess
        
    end
    
    properties % objects
        
        init_point@occult.Parameters; % keep track of the initial point 
        true_point@occult.Parameters; % if the data is simulated and the true value is known
        best_point@occult.Parameters;

        points@occult.Parameters; % a chain of trial points
        results@table; % translate the MCMC results from a matrix of points to a table
        prog@util.sys.ProgressBar; % print out the progress

    end
    
    properties % inputs/outputs
        
        input_flux; % the input lightcurve to be matched
        input_errors; % the rms error on the given flux
        
        input_R; % best estimate for the star's size
        input_v; % best estimate for the fields' velocity
        
        counter = 0; 
        num_successes = 0; % number of steps that succeeded
        
        posterior_prob; % a matrix of the probability distribution found by summing the points
        post_pars; % a cell array of the parameters used when calculating the posterior (the other parameters are integrated)
        post_x_axis; % axis points for the first parameter
        post_y_axis; % axis points for the second parameter
        post_options; % an InputVars object with the options used to generate the posterior
        
        title_strings; % a string translating the shorthands like 'r' to a title like 'Occulter radius [FSU]'
        
    end
    
    properties % switches/controls
        
        num_chains = 1; % can run multiple chains at the same time
        max_num_chains = 10; % using automated methods to choose starting points
        num_steps = 10000; % total number of steps (including burn-in)
        num_burned = 1000; % number of steps to burn at the begining of each chain
        
        step_sizes = [0.1 0.1 3 0.1]; % in order of parameters: step size for each parameter
        circ_bounds = [0 0 0 0]; % in order of parameters: which par gets a circular boundary condition
        
        par_list = {'r', 'b', 'v', 'R'}; % these parameters are chosen randomly each step. The rest are taken from the generator's parameters
        
%         use_bank = true; % use a filter bank (needs to be loaded manually!) to find an initial position (otherwise chose random starting point)
        
        initialization_method = 'search'; 

        use_priors = false; % a general switch to turn on/off the use of prior functions        
        prior_log_functions = {}; % can input a cell array with a different function per parameter (leave empty for uniform prior). min/max values are taken from generator
        
        plot_every = 100; % when plotting, how many steps go by between plots (set to 0 for no plotting)
        show_num_chains = 4; % the maximum number of chains to show on the plot
        
        show_chain_pars = {'r', 'b', 'v'}; % which parameters are shown, by default, when plotting the chain
        show_posterior_pars = {'r', 'v'}; % which parameters are shown, by default, when plotting posteriors
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)

        default_plot_every = 1;
        default_show_num_chains = 4;
        default_show_chain_pars = {'r', 'b', 'v'}; 
        default_show_posterior_pars = {'r', 'v'}; 
        
        velocity_disp = 3.5; % dispersion of velocities (e.g., in the Kuiper belt ~3.5 FSU/s)
        
        brake_bit = 1; % set to 0 when running, set to 1 to stop run
        
        % save when we run sampleGridLikelihood
        log_likelihood_map; 
        map_r;
        map_b;
        map_v;
        map_R; 
        
        version = 1.06;
        
    end
    
    methods % constructor
        
        function obj = MCMC(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'occult.,MCMC')
                if obj.debug_bit>1, fprintf('MCMC copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('MCMC constructor v%4.2f\n', obj.version); end
            
                obj.gen = occult.CurveGenerator;
                
                obj.init_point = occult.Parameters; 
                
                obj.best_point = occult.Parameters; 
                
                obj.prog = util.sys.ProgressBar;
                
                obj.results = table;
                
                obj.generateTitles;
                
            end
            
        end
        
        function generateTitles(obj)
            
            obj.title_strings.R = 'Stellar radius [FSU]';
            obj.title_strings.r = 'Occulter radius [FSU]';
            obj.title_strings.b = 'Impact parameter [FSU]'; 
            obj.title_strings.v = 'Velocity [FSU/s]';
            obj.title_strings.t = 'Mid-occultation time [s]'; 
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.points = occult.Parameters.empty; 
            obj.results = table.empty; 
            
            obj.init_point = occult.Parameters;
            obj.best_point = occult.Parameters;
            
            obj.counter = 0; 
            obj.num_successes = 0;
            
            if ~isempty(obj.gui) && obj.gui.check
                cla(obj.gui.axes_posterior); 
            end
            
            obj.gen.reset;
            
            if ~isempty(obj.true_point)
                obj.gen.lc.pars.copy_from(obj.true_point);
                obj.gen.getLightCurves;
            end
            
        end
        
        function resetLogLikelihoodMap(obj)
            
            obj.log_likelihood_map = []; 
            obj.map_r = [];
            obj.map_b = [];
            obj.map_v = [];
            obj.map_R = []; 
            
        end
        
    end
    
    methods % getters
        
        function val = get.gen(obj)
            
            if isempty(obj.gen)
                obj.gen = occult.CurveGenerator;
            end
            
            val = obj.gen;
            
        end
        
        function val = getSummary(obj)
            
            val = ''; % to be continued
            
        end
        
        function val = getProgress(obj)
            
            val = sprintf('chain: %d / %d points | %s', obj.counter, obj.num_steps, obj.prog.show); 
            
        end
        
        function T = getTable(obj, varargin) % depricated, replaced by pointsToTable
            
            input = util.text.InputVars;
            input.input_var('good', false); % only get good chains and skip points in the burn sample
            input.scan_vars(varargin{:});
            
            p = obj.points;
            
            for ii = 1:size(p,2)
                [p(:,ii).chain] = deal(ii); % make sure each point knows which chain it came from
            end
            
            pass = obj.findGoodChains;
            
            if input.good
                p = p(obj.num_burned+1:end, pass);
            end
            
            R = [p.R]';
            r = [p.r]';
            b = [p.b]';
            v = [p.v]';
            chi2 = [p.chi2]';
            ndof = [p.ndof]';
            logl = [p.logl]';
            likelihood = [p.likelihood]';
            counts = [p.counts]';
            weight = [p.weight]'; 
            chain = [p.chain]'; 
            
            T = table(R,r,b,v,chi2,ndof,logl,likelihood,counts,weight,chain); 
            
            if ~input.good % add burn and good columns to the table
                burn = [p.burn]';
                good = [p.good]'; 
                T2 = table(burn, good); 
                T = [T,T2]; 
            end
            
        end
        
        function val = getResult(obj, name, skip_burn, chain_num)
            
            if nargin<3 || isempty(skip_burn)
                skip_burn = true;
            end
            
            if nargin<4 || isempty(chain_num)
                chain_num = [];
            end
            
            if ~isempty(obj.results)
                
                val = obj.results;
                
                if skip_burn
                    val = val(val.burn==0,:);
                end
                
                if ~isempty(chain_num)
                    val = val(ismember(val.chain, chain_num),:);
                end
                
                val = val.(name); 
                
            elseif ~isempty(obj.points)
                
                val = obj.points;
                
                if skip_burn
                    val = val(obj.num_burned+1:end,:); 
                end
                
                if ~isempty(chain_num)
                    val = val(:,chain_num);
                end
                
                val = [val.(name)]';
                
            else
                val = [];
            end
            
        end
        
        function val = getNumChains(obj)
            
            if ~isempty(obj.points)
                val = size(obj.points,2);
            elseif ~isempty(obj.results)
                val = length(unique(obj.results.chain));
            else
                val = 0;
            end
            
        end
        
    end
    
    methods % setters
        
        function set.input_flux(obj, val)
            
            if ~isequal(obj.input_flux, val)
                
                obj.input_flux = val;
                
                obj.gen.getLightCurves; 
                
                obj.input_flux = util.img.pad2size(obj.input_flux, [size(obj.gen.lc.flux,1),1], NaN); 
                
                if ~isempty(obj.gui) && obj.gui.check
                    cla(obj.gui.axes_lightcurve); 
                end
                
                obj.resetLogLikelihoodMap;
                
            end
            
        end
        
        function set.input_R(obj, val)
            
            if ~isequal(obj.input_R, val)
                
                obj.input_R = val;
                
                obj.resetLogLikelihoodMap; 
                
            end
            
        end
        
        function set.true_point(obj, val)
            
            if isempty(val)
                obj.true_point = occult.Parameters.empty;
                obj.gen.reset; 
            else
                obj.true_point = val;                
                obj.gen.reset;
                obj.gen.lc.pars.copy_from(val); % also update the generator so it shows the true LC if plotting before a run starts
                obj.gen.getLightCurves;
            end
            
        end
        
    end
    
    methods % default setups
        
        function setupQuickScan(obj)
            
            obj.use_priors = 0;
            obj.initialization_method = 'random'; 
            obj.gen.v_range = [0.2 30]; 
            obj.par_list = {'r'  'b'  'v'};
%             obj.step_sizes = [0.5 0.5 3];
            obj.circ_bounds = [0 0 0];
%             obj.num_steps = 10000;
%             obj.num_burned = 1000; 

        end
        
        function setupDeepScan(obj)
            
            obj.par_list = {'r'  'b'  'v', 'R'};
            obj.step_sizes = [0.25 0.25 3 0.1];
            
            obj.initialization_method = 'random'; 
            obj.circ_bounds = [0 0 0 0];
            obj.gen.R_range = [0 10]; 
            
            obj.use_priors = 1; 
            obj.prior_log_functions = {'', '', '', @(R) -(R-obj.input_R).^2 ./ (2.*(obj.input_R.*0.1).^2) }; 
%             obj.num_steps = 10000;
%             obj.num_burned = 1000; 

        end
        
        function setupNarrowVelocityPrior(obj) % TODO: allow user input to change the range or input physical parameters
        
            obj.gen.v_range = obj.input_v - [5.15 3.65]; % apply the range of KBO motion between circular and parabolic orbits (in FSU/s)
            if obj.gen.v_range(1)<1, obj.gen.v_range(1) = 1; end % make sure the slowest velocity is not too slow
            
            obj.step_sizes(3) = 0.1; 
            
        end
        
        function setupWideVelocityPrior(obj)
            
            obj.gen.v_range = [1 40]; 
            obj.step_sizes(3) = 1;
            
        end
                
    end
    
    methods % calculations
        
        function loadTemplateBank(obj, varargin)
            
            input = util.text.InputVars; 
            input.input_var('distance', 40, 'distance_au');
            input.input_var('frame_rate', 25); 
            input.scan_vars(varargin{:}); 
            
            f = sprintf('templates_%dAU_%dHz.mat', floor(input.distance), floor(input.frame_rate)); 
            
            f = fullfile(getenv('DATA'), 'WFAST/occultations/', f); 
            
            if exist(f, 'file')
                if obj.debug_bit 
                    fprintf('Loading templates from %s\n', f); 
                end
                
                load(f); 
                
                obj.bank = bank;
                
                % maybe copy the ranges into this object's main generator? 
                
            else
                error('Cannot find file: %s', f); 
            end
            
        end
        
        function getRandomLightcurve(obj)
            
            obj.gen.lc.pars.reset; 
            obj.gen.randomLC; 
            
        end
        
        function useSimulatedInput(obj, varargin)
            
            obj.gen.getLightCurves;
            obj.input_flux = obj.gen.lc.flux;
            
            obj.true_point = occult.Parameters;
            obj.true_point.copy_from(obj.gen.lc.pars); 
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.updateLightcurves; 
            end
            
        end
        
        function useSimulatedNoisyInput(obj, varargin)
            
            obj.gen.getLightCurves;
            obj.gen.generateNoise;
            obj.input_flux = obj.gen.lc.flux_noisy;
            
            obj.true_point = occult.Parameters;
            obj.true_point.copy_from(obj.gen.lc.pars); 
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.updateLightcurves; 
            end
            
        end
        
        function startup(obj)
            
%             if length(obj.gen.r)>1
%                 error('Must setup generator with all scalar parameters!');
%             end
            
%             obj.points(obj.num_steps,obj.num_chains) = occult.Parameters; 
            
            if ~isempty(obj.gui) && obj.gui.check
                cla(obj.gui.axes_chain);
                cla(obj.gui.axes_lightcurve);
            end
            
            obj.init_point(1,obj.num_chains) = occult.Parameters; 
            
            obj.num_successes = zeros(1, obj.num_chains); 
            
        end
        
        function finishup(obj)
            
            obj.points = obj.points(1:obj.counter, :); 
            obj.prog.finish;
            
            obj.brake_bit = 1; 
            
            if ~isempty(obj.gui) && obj.gui.check 
                
                obj.gui.update;
                
                if obj.counter>obj.num_burned
                    obj.showPosterior; 
                end
            end
            
            
            
        end
        
        function gotoBestTemplateInitialPoint(obj)
            
            if isempty(obj.bank)
                error('Must load a template bank first!'); 
            end
            
            f = obj.input_flux - 1;
            
            ff = util.vec.convolution(obj.bank.kernels, f, 'cross', 1); % filtered flux
            
            [~,idx] = util.stat.max2(abs(ff)); % find the peak in the filtered flux
            
            pars = obj.bank.pars(idx(2)); % grab a struct with some parameters
            
            for ii = 1:obj.num_chains
                
                obj.init_point(ii).copy_from(obj.gen.lc.pars); % start by getting all parameters from generator

                for jj = 1:length(obj.par_list)

                    name = obj.par_list{jj};

                    if isfield(pars, name)                    
                        obj.init_point(ii).(name) = pars.(name); 
                    end

                end

                if ismember('t', obj.par_list)
                    obj.init_point(ii).t = 0; % the template banks all have t=0 
                end

            end
            
            obj.gen.lc.pars.copy_from(obj.init_point); 
            
        end
        
        function gotoRandomInitialPoint(obj)
            
            for ii = 1:obj.num_chains

                for jj = 1:length(obj.par_list)

                    name = obj.par_list{jj}; 
                    range = obj.gen.([name '_range']); 

                    S = substruct('.', name, '()', {ii});
                    obj.gen = subsasgn(obj.gen, S, diff(range).*rand + range(1)); 

                end
                
                obj.gen.lc.pars.chi2(ii) = NaN;
                obj.gen.lc.pars.ndof(ii) = NaN; 
                obj.gen.lc.pars.logl(ii) = NaN; 
                obj.gen.lc.pars.likelihood(ii) = 0; 
    
            end 
            
            obj.init_point = occult.Parameters.empty;
            obj.init_point(1,obj.num_chains) = occult.Parameters;
            obj.init_point.copy_from(obj.gen.lc.pars); 
            
        end
                
        function [log_likelihood, r, b, v, R] = sampleGridLikelihood(obj, varargin)
            
            if isempty(obj.input_flux)
                error('Cannot run MCMC without an input flux! Use useSimulatedInput() or useSimulatedNoisyInput() or input a lightcurve manually.'); 
            end
            
            input = util.text.InputVars;
            input.input_var('star_size', obj.input_R, 'stellar size');
            input.input_var('radius_step', 0.25); 
            input.input_var('impact_step', 0.25); 
            input.input_var('velocity_step', 2.5); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.star_size)
                input.star_size = 1; % in case input_R is not given
            end
            
            r = obj.gen.r_range(1):input.radius_step:obj.gen.r_range(2); 
            b = obj.gen.b_range(1):input.impact_step:obj.gen.b_range(2); 
            v = obj.gen.v_range(1):input.velocity_step:obj.gen.v_range(2); 
            R = input.star_size; 
            
            M = ones(length(b), length(v)); % matrix to fit the sizes of b and v
            
            B = b'.*M;
            V = v.*M;
            
            obj.gen.reset;
            obj.gen.R = input.star_size;
            obj.gen.t = 0;
            obj.gen.b = B(:); % linearize to get all options
            obj.gen.v = V(:); % linearize to get all options
            
            log_likelihood = NaN(length(r), length(b), length(v)); 
            
            for ii = 1:length(r)
                
                obj.gen.r = r(ii); 
                
                LL = obj.calcLikelihood; 
                
                log_likelihood(ii,:,:) = reshape(LL, [length(b), length(v)]); 
                
            end
            
            obj.log_likelihood_map = log_likelihood;
            obj.map_r = r;
            obj.map_b = b;
            obj.map_v = v;
            obj.map_R = R; 
            
        end
        
        function points = searchGoodPoints(obj, varargin) % find good starting points
            
            if ~isempty(obj.input_errors)
                e = obj.input_errors;
            else
                e = 1./obj.gen.snr;
            end
            
            chi2 = nansum(((obj.input_flux - 1)/e).^2);
            ndof = nnz(~isnan(obj.input_flux)) - 3; % assume parameters are r,b,v, only
            
%             base_likelihood = chi2cdf(chi2, ndof, 'upper'); 
            logl_base = -0.5*chi2;

            if isempty(obj.log_likelihood_map)
                [LL, r, b, v, R] = obj.sampleGridLikelihood(varargin{:});
            else
                LL = obj.log_likelihood_map; 
                r = obj.map_r;
                b = obj.map_b;
                v = obj.map_v;
                R = obj.map_R; 
            end
            
            idx = LL>base_likelihood+log(2);
            
            stats = regionprops3(idx);
                        
            points = occult.Parameters.empty;
            
            N = min(obj.max_num_chains, height(stats)); 
            
            for ii = 1:N
                
                pos = round(stats{ii, 'Centroid'}); % the order is x,y,z so b,r,v
                
                points(ii) = occult.Parameters;
                points(ii).r = r(pos(2)); 
                points(ii).b = b(pos(1)); 
                points(ii).v = v(pos(3)); 
                points(ii).R = R; 
                
            end
            
        end
        
        function takeStep(obj)
            
            for jj = 1:length(obj.par_list)

                name = obj.par_list{jj}; 
                range = obj.gen.([name '_range']); 
                obj.gen.(name) = normrnd(obj.gen.(name), obj.step_sizes(jj).*ones(1,obj.num_chains)); 
                
                for ii = 1:obj.num_chains

                    S = substruct('.', name, '()', {ii});

                    if obj.circ_bounds(jj)

                        for ii = 1:1000
                            
                            if subsref(obj.gen, S)>range(2)
                                delta = subsref(obj.gen, S) - range(2);
                                obj.gen = subsasgn(obj.gen, S, range(1) + delta); % circle back the rest of the way from the lower bound
                            elseif subsref(obj.gen, S)<range(1)
                                delta = range(1) - subsref(obj.gen, S);
                                obj.gen = subsasgn(obj.gen, S, range(2) - delta); % circle back the rest of the way from the upper bound
                            else
                                break; % not out of bounds! 
                            end
                            
                        end

                    else % reflect back
                        
                        for ii = 1:1000
                            
                            if subsref(obj.gen, S)>range(2)
                                delta = subsref(obj.gen, S) - range(2);
                                obj.gen = subsasgn(obj.gen, S, range(2) - delta); % reflect the next point off the upper edge of the range
                            elseif subsref(obj.gen, S)<range(1)
                                delta = range(1) - subsref(obj.gen, S);
                                obj.gen = subsasgn(obj.gen, S, range(1) + delta); % leave the point at the lower edge of the range
                            else
                                break; % not out of bounds! 
                            end
                            
                        end
                    end 

                end % for ii (num_chains)

            end % for jj (par list)

        end
        
        function LL = calcLikelihood(obj)
            
            obj.gen.getLightCurves; 
            
            if ~isempty(obj.input_errors)
                e = obj.input_errors;
            else
                e = 1./obj.gen.snr;
            end

            f1 = obj.gen.lc.flux;
            f2 = obj.input_flux;

            chi2 = nansum(((f1-f2)./e).^2, 1); 

            obj.gen.lc.pars.chi2 = chi2;
            obj.gen.lc.pars.ndof = sum(~isnan(f1),1)-length(obj.par_list);
            
            for ii = 1:length(chi2)

%                 L = chi2cdf(chi2(ii),, 'upper'); % do we need to subtract these parameters from the ndof??
                LL = -0.5.*chi2(ii); % moved to log-likelihood for better stability
                
                if obj.use_priors
                    for jj = 1:length(obj.prior_log_functions)
                        if ~isempty(obj.prior_log_functions{jj})
                            values = obj.gen.(obj.par_list{jj});
                            % multiply by the value of the function at the given point in parameter space
                            LL = LL + obj.prior_log_functions{jj}(values(ii)); % assume prior functions give log as well 
                        end
                    end
                end
                
                obj.gen.lc.pars.logl(ii) = LL; 
                obj.gen.lc.pars.likelihood(ii) = exp(LL);

                % likelihood is zero for parameters that are out of bounds
                for jj = 1:length(obj.par_list)

                    S = substruct('.', obj.par_list{jj}, '()', {ii});

                    value = subsref(obj.gen, S);
                    range = obj.gen.([obj.par_list{jj} '_range']); 

                    if value<range(1) || value>range(2)
                        obj.gen.lc.pars.logl(ii) = -Inf; 
                        obj.gen.lc.pars.likelihood(ii) = 0;
                        break; % no need to check other parameters
                    end

                end
                
            end
            
            if nargout>0
                LL = obj.gen.lc.pars.logl;
            end

        end
        
        function run(obj, varargin)
            
            import util.text.cs;
            
            input = util.text.InputVars;
            input.input_var('plot', false);
            input.input_var('ax', [], 'axes', 'axis');
            input.scan_vars(varargin{:});
            
            if ~ismember('plot', input.list_scan_properties) && ~isempty(obj.gui) && obj.gui.check
                input.plot = true; 
            end
            
            if isempty(obj.input_flux)
                error('Cannot run MCMC without an input flux! Use useSimulatedInput() or useSimulatedNoisyInput() or input a lightcurve manually.'); 
            end
            
            obj.reset;
            obj.startup;
            
            on_cleanup = onCleanup(@() obj.finishup); % make sure to wrap up: truncate the chain, shut down the progress bar
            
            if input.plot
                
                if isempty(input.ax)
                    if ~isempty(obj.gui) && obj.gui.check
                        input.ax = obj.gui.axes_chain; 
                    else
                        input.ax = gca;
                    end
                end
                
            end
            
            obj.prog.start(obj.num_steps);

            if cs(obj.initialization_method, 'search')
                obj.init_point = obj.searchGoodPoints;
                if isempty(obj.init_point)
                    error('Could not find any good starting points...'); 
                end
                obj.num_chains = length(obj.init_point); 
            elseif cs(obj.initialization_method, 'bank')
                obj.gotoBestTemplateInitialPoint; 
            elseif cs(obj.initialization_method, 'random')
                obj.gotoRandomInitialPoint;
            else
                error('Unknown initialization option "%s". Use "search", "bank" or "random" instead...', obj.initialization_method); 
            end
            
            obj.gen.reset;
            obj.gen.lc.pars.copy_from(obj.init_point); 
            
            obj.calcLikelihood; % get the first point likelihood/chi2

            obj.points(obj.num_steps, obj.num_chains) = occult.Parameters; 

            obj.points(1,:).copy_from(obj.gen.lc.pars); % automatically accept first point... 
            
            obj.brake_bit = 0;
            
            obj.counter = 1; 
            
            for jj = 2:obj.num_steps

                if obj.brake_bit, break; end
                
                obj.prog.showif(jj);

                obj.takeStep; % move the generator to a nearby random point

                obj.calcLikelihood;
                
                if obj.debug_bit>2
                    obj.gen.lc.pars.printout;
                end
                
                obj.points(jj,:).copy_from(obj.gen.lc.pars); % copy all parameters from generator into different points on the chains
                
                for ii = 1:obj.num_chains
                    
%                     ratio = obj.gen.lc.pars.likelihood(ii)./obj.points(jj-1,ii).likelihood; % compare the new potential point with the last point

                    log_ratio = obj.gen.lc.pars.logl(ii) - obj.points(jj-1,ii).logl;
                    
                    if isnan(log_ratio) || log(rand)<log_ratio % if ratio is big, we are likely to take the new point
                        if length(obj.num_successes)<ii
                            obj.num_successes(ii) = 1;
                        else
                            obj.num_successes(ii) = obj.num_successes(ii) + 1;
                        end
                        
                        obj.points(jj,ii).counts = 1;
                        
                        if obj.counter>obj.num_burned && ~(obj.best_point.logl>obj.points(jj,ii).logl) % keep track of the best fit point
                            obj.best_point.copy_from(obj.points(jj,ii)); 
                        end

                    else % point is rejected, repeat the previous point
                        obj.points(jj,ii).copy_from(obj.points(jj-1,ii));
                        obj.points(jj,ii).counts = obj.points(jj-1,ii).counts + 1;
                        obj.points(jj,ii).weight = 0; % don't count these points! 
                    end

                end
                
                obj.gen.lc.pars.copy_from(obj.points(jj,:)); % get the new point positions back into the generator in case we decided to go back to the previous point 
                
                obj.counter = obj.counter + 1;
                
                if input.plot && obj.plot_every>0 && (obj.plot_every==1 || mod(obj.counter, obj.plot_every)==1)
                    
                    obj.plot('burn', 1, 'ax', input.ax, 'sparse', obj.plot_every); 
                    
                    if ~isempty(obj.gui) && obj.gui.check
                        obj.gui.updateLightcurves; 
                    end
                    
                    drawnow;
                    
                end
                                
                if ~isempty(obj.gui) && obj.gui.check
                    obj.gui.update;
                end
                
            end
                        
        end
        
        function [pass, likelihood, acceptance] = calcChainProps(obj, varargin)

            input = util.text.InputVars;
            input.input_var('success_ratio', 0.05); % minimal number of successes (relative to total steps) for a chain to be included
            input.input_var('likelihood', []); % minimal mean log-likelihood needed for a chain to be included
            input.scan_vars(varargin{:}); 
            
            if isempty(input.likelihood)
                input.likelihood = -Inf; 
            end
            
            for ii = 1:obj.getNumChains
                
                acceptance(ii) = obj.num_successes(ii)./obj.counter;
                likelihood(ii) = nanmean(obj.getResult('logl', true, ii));
                
                if ~isempty(input.success_ratio)
                    pass(ii) = acceptance(ii) >= input.success_ratio && ... 
                        likelihood(ii) >= input.likelihood;
                end
                
            end

        end
        
        function [pass, mean_chi2] = findGoodChains(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('sigma', 3); % how many MADs above the noise should outliers be
            input.scan_vars(varargin{:}); 
            
            N = obj.getNumChains;
            
            mean_chi2 = zeros(1, N);
            
            for ii = 1:N
%                 mean_chi2(ii) = nanmean([obj.points(:,ii).chi2]); 
                mean_chi2(ii) = nanmean(obj.getResults('chi2', true, ii)); 
            end
            
            if N > 4
                scat = mad(mean_chi2,1); % median deviation skips outliers
                aver = median(mean_chi2); % median average skips outliers
                pass = mean_chi2 - aver < scat * input.sigma; % low chi2 are also accepted, but high outliers are removed
            else
                pass = true(1,N); 
            end
            
        end
        
        function [N, x, y] = calcPosterior(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('pars', obj.show_posterior_pars); 
            input.input_var('success_ratio', 0.05); % minimal number of successes (relative to total steps) for a chain to be included
            input.input_var('likelihood', []); % minimal mean likelihood needed for a chain to be included
            input.input_var('weights', false, 'use_weights');
            input.input_var('burn', false); 
            input.input_var('oversample', 5); 
            input.scan_vars(varargin{:}); 
            
            if obj.counter==0
                N = 0;
                x = [];
                y = [];
                return;
            end
            
            chain_num = find(obj.calcChainProps('success', input.success_ratio, 'likelihood', input.likelihood));
            
            x = obj.getResult(input.pars{1}, ~input.burn, chain_num); 
            y = obj.getResult(input.pars{2}, ~input.burn, chain_num); 
            w = obj.getResult('weight', ~input.burn, chain_num); 
            
%             x = [obj.points(idx1, idx2).(input.pars{1})];
%             y = [obj.points(idx1, idx2).(input.pars{2})];
%             w = logical([obj.points(idx1, idx2).weight]);
            
            if input.weights
                x = x(w);
                y = y(w); 
            end
            
            idx_x = find(strcmp(input.pars{1}, obj.par_list), 1, 'first');
            dx = obj.step_sizes(idx_x)/input.oversample; 

            idx_y = find(strcmp(input.pars{2}, obj.par_list), 1, 'first');
            dy = obj.step_sizes(idx_y)/input.oversample; 
            
            if isempty(x)
                N = [];
                x = [];
                y = [];
            else
                [N, Ex, Ey] = histcounts2(x, y, 'BinWidth', [dx, dy], 'Normalization','probability');
                N = N';
                x = repmat(Ex(1:end-1)+dx/2, [size(N,1),1]);                 
                y = repmat(Ey(1:end-1)'+dy/2, [1, size(N,2)]); 
            end
            
            obj.posterior_prob = N;
            obj.post_pars = input.pars;
            obj.post_x_axis = x;
            obj.post_y_axis = y; 
            obj.post_options = input; 
            
        end
        
        function [med, lower, upper] = getMedianAndBounds(obj, par, percentile, skip_burn)
            
            if nargin<3 || isempty(percentile)
                percentile = 68; 
            end
            
            if percentile>1
                percentile = percentile/100;
            end
            
            if nargin<4 || isempty(skip_burn)
                skip_burn = 1;
            end
            
%             p = obj.points;
%             if skip_burn
%                 p = p(1:obj.num_burned,:); 
%             end
%             
%             values = [p(:).(par)];
            
            values = obj.getResult(par, skip_burn); 
            N = length(values); 
            
            % for even N these will be the same
            mid_idx1 = floor((N+1)/2); 
            mid_idx2 = ceil(N/2); 
            
            lower_frac = (1-percentile)/2; 
            lower_idx = floor(lower_frac.*N); 
            
            upper_frac = (1+percentile)/2; 
            upper_idx = ceil(upper_frac.*N); 
            
            v_sort = sort(values); 
            
            med = (v_sort(mid_idx1) + v_sort(mid_idx2))/2; % median value is the average of middle two elements, unless N is odd, so they are the same element
            lower = v_sort(lower_idx); 
            upper = v_sort(upper_idx); 
            
        end
        
        function [center, lower, upper] = getSmallestBounds(obj, par, percentile, skip_burn)
            
            if nargin<3 || isempty(percentile)
                percentile = 68; 
            end
            
            if percentile>1
                percentile = percentile/100;
            end
            
            if nargin<4 || isempty(skip_burn)
                skip_burn = 1;
            end
            
%             p = obj.points;
%             if skip_burn
%                 p = p(1:obj.num_burned,:); 
%             end
%             
%             values = [p(:).(par)];
            values = obj.getResult(par, skip_burn); 
            
            [lower, upper, center] = util.stat.dist_bounds(values, 'fraction', percentile); 
            
%             N = length(values); 
%             
%             sorted_values = sort(values); 
%             
%             func = @(lb) obj.find_upper_bound(sorted_values, lb, percentile) - lb; % loss function: distance between lower and upper bounds
%             
%             lower_bound_guess = sorted_values(floor((1-percentile)/2.*N));
%             
%             lower = fminsearch(func, lower_bound_guess); 
%             upper = obj.find_upper_bound(sorted_values, lower, percentile); 
%             center = (upper+lower)/2; 
                        
        end
        
        function upper_bound = find_upper_bound(~, sorted_values, lower_bound, percentile)
        
            N = length(sorted_values);
            idx_low = find(sorted_values>lower_bound, 1, 'first'); 
            idx_high = ceil(idx_low+percentile*N); 
            if idx_high<=N
                upper_bound = sorted_values(idx_high); 
            else
                upper_bound = Inf;
            end
            
%             fprintf('low: %4.2f | high: %4.2f\n', lower_bound, upper_bound); 
            
        end
        
        function flagBadPoints(obj)
            
            pass = obj.findGoodChains; 
            
            [obj.points(1:obj.num_burned, :).burn] = deal(true);
            [obj.points(1:obj.num_burned, :).good] = deal(false);
            [obj.points(:, ~pass).good] = deal(false);
            
        end
        
        function pointsToTable(obj)
            
            if ~isempty(obj.points)
                
                obj.flagBadPoints;
                
                names = [obj.par_list, {'chi2', 'ndof', 'logl', ...
                    'counts', 'weight', 'chain', 'burn', 'good'}];
                
                obj.results = cell2table(cell(numel(obj.points), length(names)), ...
                    'VariableNames', names);
                
                for ii = 1:size(obj.points,2)
                    [obj.points(:,ii).chain] = deal(ii); % make sure each point knows which chain it came from
                end
                
                [obj.points(1:obj.num_burned,:).burn] = deal(1); 
                
                for ii = 1:length(names)
                    obj.results.(names{ii}) = [obj.points.(names{ii})]'; 
                end
                
                obj.points = occult.Parameters.empty;
                
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function plot(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('pars', obj.show_chain_pars); 
            input.input_var('burn', false); 
            input.input_var('sparse', []); % show only one point every N points
            input.input_var('full_titles', false, 'titles'); 
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('font_size', 18); 
            input.input_var('hold', false); 
            input.scan_vars(varargin{:});
            
            if obj.counter==0
                return; 
            end
            
            if isempty(input.ax)
                
                if ~isempty(obj.gui) && obj.gui.check
                    input.ax = obj.gui.axes_chain; 
                else
                    input.ax = gca;
                end
                
            end
            
            hold_state = input.ax.NextPlot; 
            
            if input.hold
                input.ax.NextPlot = 'add'; 
            end
            
            if isempty(input.pars)
                input.pars = obj.par_list(1:min(3,length(obj.par_list))); % just show the first 3 pars
            end
            
            if ischar(input.pars)
                input.pars = {input.pars};
            end
            
            if isempty(input.sparse)
                step = 1;
            else
                step = input.sparse; 
            end
            
            sz = ones(obj.counter, 1); 
%             colors = {'b', 'r', 'g', 'm'};
            shapes = {'o', 'p', 's', 'v'}; 
            
            N = min(obj.counter, obj.num_burned); 
            if input.burn
                sz(1:N) = 5;
            else
                sz(1:N) = NaN;
            end
            
            sz(N+1:end) = 20; 
            
            if length(input.pars)==1 % just show the distribution of this one parameter I guess... 
                
                histogram(input.ax, p1);
                
            elseif length(input.pars)==2
                
                error('need to fix this option to be more like 3 pars...'); 
                
                p2 = [obj.points.(input.pars{2})];
                                
                scatter(input.ax, p1(obj.num_burned:end), p2(obj.num_burned:end), 20, counts(obj.num_burned:end), 'o', 'filled');
                
                if input.burn
                    hold(input.ax, 'on'); 
                    scatter(input.ax, p1(1:obj.num_burned), p2(1:obj.num_burned), 10, counts(1:obj.num_burned), 'o');
                end

                % the X limits and labels are added in the end
                ylabel(input.ax, input.pars{2});
                input.ax.YLim = obj.gen.([input.pars{2} '_range']);

                hcb = colorbar(input.ax); 
                ylabel(hcb, 'number of repeats'); 
                
            elseif length(input.pars)==3
                
                h = findobj(input.ax, 'Type', 'Scatter'); 
                
                N = min(obj.num_chains, obj.show_num_chains); 
                
                sz = sz(1:step:obj.counter);
                
                for ii = 1:N
                    % TODO: refactor to use the table not just the points
                    p1 = [obj.points(1:step:obj.counter,ii).(input.pars{1})]';
                    p2 = [obj.points(1:step:obj.counter,ii).(input.pars{2})]';
                    p3 = [obj.points(1:step:obj.counter,ii).(input.pars{3})]';
                    lkl = [obj.points(1:obj.counter,ii).likelihood]'; % used only for calculating the mean likelihood of each chain
                    
                    if isempty(h) || length(h)~=N || ~isvalid(h(ii))
                        N_colors = size(input.ax.ColorOrder,1);
                        col_idx = mod(ii-1, N_colors) + 1; 
                        sha_idx = floor((ii-1)/N_colors) + 1; 
                        h(ii) = scatter3(input.ax, p1, p2, p3, sz, input.ax.ColorOrder(col_idx,:), shapes{sha_idx}, 'filled');
                        h(ii).DisplayName = sprintf('Chain %d (acc=%3.1f%%, lkl=%4.2g)', ...
                            ii, 100*obj.num_successes(ii)./obj.counter, nanmean(lkl(obj.num_burned+1:end,1))); 
                        input.ax.NextPlot = 'add';                 
                    else
                        h(N-ii+1).XData = p1; 
                        h(N-ii+1).YData = p2;
                        h(N-ii+1).ZData = p3;
                        h(N-ii+1).SizeData = sz;
                        h(N-ii+1).DisplayName = sprintf('Chain %d (acc=%3.1f%%, lkl=%4.2g)', ...
                            ii, 100*obj.num_successes(ii)./obj.counter, nanmean(lkl(obj.num_burned+1:end,1),1)); 
                        input.ax.NextPlot = 'add';  
                    end

                end
                
                input.ax.NextPlot = 'replace'; 
                
                % the X limits and labels are added in the end
                if input.full_titles
                    ylabel(input.ax, obj.title_strings.(input.pars{2}));
                else
                    ylabel(input.ax, input.pars{2});
                end
                
                input.ax.YLim = obj.gen.([input.pars{2} '_range']);
                
                if input.full_titles
                    zlabel(input.ax, obj.title_strings.(input.pars{3}));
                else
                    zlabel(input.ax, input.pars{3});
                end
                
                input.ax.ZLim = obj.gen.([input.pars{3} '_range']);
                
%                 input.ax.CLim = [0 1]; 
%                 hcb = colorbar(input.ax); 
%                 ylabel(hcb, 'likelihood'); 
                
            end
        
            if input.full_titles
                xlabel(input.ax, obj.title_strings.(input.pars{1}));
            else
                xlabel(input.ax, input.pars{1});
            end
            
            input.ax.FontSize = input.font_size;
            
            input.ax.XLim = obj.gen.([input.pars{1} '_range']);
            
            input.ax.NextPlot = hold_state;
            
            for ii = 1:length(h)
                h(ii).ButtonDownFcn = @obj.callback_click_point; 
            end
            
            legend(input.ax, 'Location', 'NorthEastOutside'); 
            
        end
        
        function callback_click_point(obj, ~, ev)
            
            if obj.debug_bit, fprintf('clicked point at '); end
            
            p = obj.points(:); 
            
            deltas = zeros(length(p), length(obj.par_list)); 

            for ii = 1:length(obj.par_list)
                deltas(:,ii) = ([p.(obj.par_list{ii})] - ev.IntersectionPoint(ii)).^2;
            end

            deltas = sum(deltas,2); 
            
            [~, idx] = min(deltas); 
            
            clicked_point = p(idx);
            
            chain_idx = floor((idx-1)/size(obj.points,1)) + 1;
            
            obj.gen.lc.pars.reset; 
            obj.gen.lc.pars.copy_from(clicked_point); 

            if obj.debug_bit, clicked_point.printout; end

            if ~isempty(obj.gui) && obj.gui.check
                cla(obj.gui.axes_lightcurve); 
                obj.plotLightcurves('ax', obj.gui.axes_lightcurve, 'color', chain_idx); 
                obj.gui.update;
            end
            
        end
        
        function plotLightcurves(obj, varargin)
            
            import util.text.print_vec; 
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('color', []); 
            input.input_var('font_size', 16); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.ax)
                
                if ~isempty(obj.gui) && obj.gui.check
                    input.ax = obj.gui.axes_lightcurve; 
                else
                    input.ax = gca;
                end
                
            end
            
            obj.gen.getLightCurves;
            t = obj.gen.lc.time;
            f = obj.gen.lc.flux; 
            
            h = findobj(input.ax, 'Type', 'Line');
            
            if ~isempty(obj.points) % there are active chains
                N = min(obj.num_chains, obj.show_num_chains); 
            elseif ~isempty(obj.true_point) % no chains, but we can show the true point
                N = 1;
            else % no chains and no true point, will show nothing
                N = 0;
            end
            
            if isempty(h) || ~isvalid(h(1)) || N~=length(h)-1
                
                delete(input.ax.Children); 
                f_input = util.img.pad2size(obj.input_flux, [size(f,1),1], 1); 
                plot(input.ax, t, f_input, 'dk', 'LineWidth', 3, 'DisplayName', 'Input flux'); 
                
                if ~isempty(obj.points) || ~isempty(obj.true_point)
                    h = obj.gen.lc.plot('ax', input.ax, 'hold', true, 'legend', true, 'noise', false, 'FontSize', input.font_size, 'number', N);
                end
                
                N = length(h); % make sure N is smaller if there aren't that many plots... 
                
                if ~isempty(input.color) % use this color for the first non-input-flux plot

                    if ischar(input.color) || length(input.color)==3
                        h.Color = input.color;
                    elseif isscalar(input.color)
                        h.Color = input.ax.ColorOrder(input.color,:);
                    end

                end

            else
                for ii = 1:N
                    h(N-ii+1).XData = t; 
                    h(N-ii+1).YData = f(:,ii); 
                    h(N-ii+1).DisplayName = sprintf('r= %4.2f | b= %4.2f | v= %4.2f', obj.gen.r(ii), obj.gen.b(ii), obj.gen.v(ii)); 
                end
            end
            
            delete(findobj(input.ax, 'Type', 'Text')); 
            
            title(input.ax, ''); 
            input.ax.YLim = [-0.2 1.4]; 
            
            input.ax.FontSize = input.font_size; 
            
            chi2 = obj.gen.lc.pars.chi2;
            
            ndof = nnz(~isnan(obj.input_flux))-length(obj.par_list);
            
            if ~isempty(obj.points)
                text(input.ax, t(10), 0.3, sprintf('chi2 (dof=%d) = %s', ...
                    ndof, print_vec(round(chi2(1:N)), ', ')), 'FontSize', input.font_size); 

                text(input.ax, t(10), 0.1, sprintf('chi2/dof = %s', ...
                    print_vec(round(chi2(1:N)./ndof,1), ', ')), 'FontSize', input.font_size); 

                text(input.ax, t(10), -0.1, sprintf('likelihood= %s', ...
                    sprintf('10^{%4.1f} ',obj.gen.lc.pars.logl(1:N)./log(10))),...
                    'FontSize', input.font_size); 
            end
            
        end
        
        function showPosterior(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('pars', obj.show_posterior_pars); 
            input.input_var('burn', false); 
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('font_size', 18); 
            input.input_var('contour', true); 
            input.input_var('inner_titles', true); 
            input.input_var('full_titles', false, 'titles'); 
            input.input_var('bounds', true); 
            input.input_var('true_point', true); 
            input.input_var('best_point', true); 
            input.input_var('horizontal', false); 
            input.input_var('hold', false); 
            input.input_var('legend', false); 
            input.scan_vars(varargin{:});
            
            if obj.counter==0
                return;
            end
            
            if ischar(input.pars)
                input.pars = {input.pars}; 
            end
            
            if isempty(input.ax)
                
                if ~isempty(obj.gui) && obj.gui.check
                    input.ax = obj.gui.axes_posterior; 
                else
                    input.ax = gca;
                end
                
            end
            
            if input.hold
                hold(input.ax, 'on'); 
            end
            
            if isempty(input.pars)
                error('Must specify which parameters to show!'); 
            elseif isscalar(input.pars)
                
                [med, lower, upper] = obj.getSmallestBounds(input.pars{1}, 68, ~input.burn);
            
                if ~input.horizontal % vertical option
                    
%                     h = histogram(input.ax, [obj.points(1:obj.num_burned,:).(input.pars{1})]); 
                    h = histogram(input.ax, obj.getResult(input.pars{1}, true)); 
                    
                    x_title = input.pars{1};
                    if input.full_titles
                        x_title = obj.title_strings.(x_title);                    
                    end

                    y_title = 'Number of points'; 

                    input.ax.YLim = input.ax.YLim;
                    
                    hold(input.ax, 'on'); 

                    if input.bounds
                        h_bounds = area(input.ax, [lower upper], [1 1].*input.ax.YLim(2), ...
                            'EdgeColor', 'none', 'FaceColor', 'm', 'FaceAlpha', 0.2); 
                        h_middle = plot(input.ax, med.*[1 1], input.ax.YLim, ':m', 'LineWidth', 2); 
                    end
                    
                    if ~isempty(obj.true_point) && input.true_point
                        h_true = plot(input.ax, obj.true_point.(input.pars{1}).*[1 1], input.ax.YLim, '-g', 'LineWidth', 2); 
                        h_true.DisplayName = sprintf('True: %s= %4.1f', ...
                            input.pars{1}, obj.true_point.(input.pars{1}));
                    end

                    if ~isempty(obj.best_point) && input.best_point
                        h_best = plot(input.ax, obj.best_point.(input.pars{1}).*[1 1], input.ax.YLim, '--r', 'LineWidth', 2); 
                        h_best.DisplayName = sprintf('Best: %s= %4.1f', ...
                            input.pars{1}, obj.best_point.(input.pars{1}));
                    end
                    
                    input.ax.XLim = obj.gen.([input.pars{1} '_range']); 
                    
                    input.ax.Children = circshift(input.ax.Children, 1); 
                    
                    h_bounds.HandleVisibility = 'off';
                    h_middle.HandleVisibility = 'off';
                    
                else % horizontal option
                    
%                     h = histogram(input.ax, [obj.points(1:obj.num_burned,:).(input.pars{1})], 'Orientation', 'Horizontal');
                    h = histogram(input.ax, obj.getResult(input.pars{1}, true), 'Orientation', 'Horizontal');
                
                    y_title = input.pars{1};
                    if input.full_titles
                        y_title = obj.title_strings.(y_title);                    
                    end

                    x_title = 'Number of points'; 

                    input.ax.XLim = input.ax.XLim;
                    
                    hold(input.ax, 'on'); 

                    if input.bounds
                        h_bounds = area(input.ax, input.ax.XLim, upper.*[1 1], lower, ...
                            'EdgeColor', 'None', 'FaceColor', 'm', 'FaceAlpha', 0.2); 
                        h_middle = plot(input.ax, input.ax.XLim, med.*[1 1], ':m', 'LineWidth', 2); 
                    end
                    
                    if ~isempty(obj.true_point) && input.true_point
                        h_true = plot(input.ax, input.ax.XLim, obj.true_point.(input.pars{1}).*[1 1], '-g', 'LineWidth', 2); 
                        h_true.DisplayName = sprintf('True: %s= %4.1f', ...
                            input.pars{1}, obj.true_point.(input.pars{1}));
                    end

                    if ~isempty(obj.best_point) && input.best_point
                        h_best = plot(input.ax, input.ax.XLim, obj.best_point.(input.pars{1}).*[1 1], '--r', 'LineWidth', 2); 
                        h_best.DisplayName = sprintf('Best: %s= %4.1f', ...
                            input.pars{1}, obj.best_point.(input.pars{1}));
                    end

                    input.ax.YLim = obj.gen.([input.pars{1} '_range']); 
                    
                    input.ax.Children = circshift(input.ax.Children, 1); 
                    
                    h_bounds.HandleVisibility = 'off';
                    h_middle.HandleVisibility = 'off';
                    
                end
                
                h.FaceAlpha = 1;
                h.FaceColor = [1 1 1].*0.7; 
                
                if input.legend
                    
                    h.DisplayName = sprintf('%s= %4.2f_{-%4.1f}^{+%4.1f}', ...
                        input.pars{1}, med, med-lower, upper-med); 
                    if input.horizontal
                        if med<mean(input.ax.YLim)
                            l_pos = 'NorthEast'; 
                        else
                            l_pos = 'SouthEast'; 
                        end
                    else
                        if med<mean(input.ax.XLim)
                            l_pos = 'NorthEast'; 
                        else
                            l_pos = 'NorthWest'; 
                        end
                    end
                end
                
            elseif length(input.pars)==2
                
                [N,x,y] = obj.calcPosterior(varargin{:}); 
                if isempty(N)
                    disp('None of the chains passed the posterior calculation conditions'); 
                    return;
                end
                
                x_title = input.pars{1};
                y_title = input.pars{2};
                if input.full_titles
                    x_title = obj.title_strings.(x_title);
                    y_title = obj.title_strings.(y_title);
                end

                if input.contour
                    [M,c] = contour(input.ax, x, y, N, '-k'); 
                else
                    imagesc(input.ax, x(1,:), y(:,1), N); 
                    axis(input.ax, 'xy'); 
                    colormap(flip(gray)); 
                    input.ax.CLim = [0 prctile(N(:), 98)]; 
                end

                
                input.ax.XLim = obj.gen.([input.pars{1} '_range']); 
                input.ax.YLim = obj.gen.([input.pars{2} '_range']); 
                
                hold(input.ax, 'on'); 

                if ~isempty(obj.true_point) && input.true_point
                    h_true = plot(input.ax, obj.true_point.(input.pars{1}), obj.true_point.(input.pars{2}), 'og', 'MarkerSize', 10, 'MarkerFaceColor', 'g'); 
                    h_true.DisplayName = sprintf('True point: %s= %4.2f | %s= %4.2f', ...
                            input.pars{1}, obj.true_point.(input.pars{1}), ...
                            input.pars{2}, obj.true_point.(input.pars{2}));
                end

                if ~isempty(obj.best_point) && input.best_point
                    h_best = plot(input.ax, obj.best_point.(input.pars{1}), obj.best_point.(input.pars{2}), '+r', 'MarkerSize', 20); 
                    h_best2 = plot(input.ax, obj.best_point.(input.pars{1}), obj.best_point.(input.pars{2}), 'or', 'MarkerSize', 20); 
                    h_best.DisplayName = sprintf('True point: %s= %4.2f | %s= %4.2f', ...
                            input.pars{1}, obj.true_point.(input.pars{1}), ...
                            input.pars{2}, obj.true_point.(input.pars{2}));
                    h_best2.HandleVisibility = 'off'; 
                end
                
            else % more than 2 pars
                error('Must specify one or two parameters!'); 
            end

            if input.inner_titles
                util.plot.inner_title(x_title, 'Position', 'Bottom', 'ax', input.ax); 
                util.plot.inner_title(y_title, 'Position', 'Left', 'ax', input.ax); 
                util.plot.inner_title('Posterior', 'Position', 'Top', 'ax', input.ax); 
            else
                xlabel(input.ax, x_title); 
                ylabel(input.ax, y_title); 
            end

            input.ax.FontSize = input.font_size; 
            
            if input.hold
                hold(input.ax, 'on'); 
            else
                hold(input.ax, 'off'); 
            end
            
            if input.legend
                legend(input.ax, 'Location', l_pos); 
            end
            
        end
        
        function showResults(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('best', false); 
            input.input_var('parent', []); 
            input.input_var('font_size', 18); 
            input.scan_vars(varargin{:});
            
            if isempty(input.parent)
               input.parent = gcf;
            end
            
            if obj.counter==0
                disp('Cannot show results without any chains!'); 
                return;
            end
            
            delete(input.parent.Children); 
            
            %%%%%%%%%%%%%%%% DENSITY %%%%%%%%%%%%%%%%%%%%%%
            
            ax_density = axes('Parent', input.parent, 'Position', [0.1 0.15, 0.4 0.4]);             
            obj.showPosterior('ax', ax_density, 'font_size', input.font_size, 'inner', 0, ...
                'title', 1, 'contour', 0, 'best', input.best, varargin{:});
            P = ax_density.Position;
            
            if ax_density.XLim(2)==ax_density.XTick(end)
                ax_density.XTick(end) = [];
            end
            
            if ax_density.YLim(2)==ax_density.YTick(end)
                ax_density.YTick(end) = [];
            end
            
            low_y = ax_density.YTick(2); 
            low_x = ax_density.XTick(1); 
            high_x = ax_density.XLim(2)*0.95; 
            
            N = sum(obj.calcChainProps);
            
            text(ax_density, high_x, low_y, sprintf('%d points\n %d chains', obj.num_steps*N, N), ...
                'FontSize', input.font_size, 'HorizontalAlignment', 'right', 'FontWeight', 'bold'); 
            
            
            %%%%%%%%%%%%%%%% VERTICAL %%%%%%%%%%%%%%%%%%%%%%
            
            ax_vertical = axes('Parent', input.parent, 'Position', [P(1) P(2)+P(4) P(3) P(4)]); 
            obj.showPosterior('ax', ax_vertical, 'font_size', input.font_size, ...
                'inner', 0, 'pars', 'r', 'Legend', 1, 'best', input.best, varargin{:}); 
            ax_vertical.XTick = []; 
            xlabel(ax_vertical, ''); 

            idx = find(ax_vertical.YTick==0); 
            if ~isempty(idx)
                ax_vertical.YTick(idx) = [];
            end
            
            %%%%%%%%%%%%%%%% HORIZONTAL %%%%%%%%%%%%%%%%%%%%%%
            
            ax_horizontal = axes('Parent', input.parent, 'Position', [P(1)+P(3) P(2) P(3) P(4)]); 
            obj.showPosterior('ax', ax_horizontal, 'font_size', input.font_size, 'inner', 0, ...
                'pars', 'v', 'horizontal', 1, 'Legend', 1, 'best', input.best, varargin{:}); 
            ax_horizontal.YTick = []; 
            ylabel(ax_horizontal, ''); 
            
            idx = find(ax_horizontal.XTick==0); 
            if ~isempty(idx)
                ax_horizontal.XTick(idx) = [];
            end
            
            %%%%%%%%%%%%%%%% LIGHTCURVE %%%%%%%%%%%%%%%%%%%%%%
            
            margin_x = 0.025;
            margin_y = 0.1;
            ax_lc = axes('Parent', input.parent, 'Position', [P(1)+P(3)+margin_x P(2)+P(4)+margin_y P(3)-margin_x P(4)-margin_y]); 
            obj.gen.reset;
            obj.gen.lc.pars.copy_from([obj.best_point, obj.true_point]); 
            obj.gen.getLightCurves; 
            t = obj.gen.lc.time;
            f = obj.gen.lc.flux; 
                        
            h = plot(ax_lc, t, util.img.pad2size(obj.input_flux, size(t), 1), 'kd', 'MarkerSize', 8, 'LineWidth', 2); 
            h.DisplayName = 'input flux';
            
            hold(ax_lc, 'on'); 
            
            if size(f,2)>1
                h_true = plot(ax_lc, t, f(:,2), '-g', 'LineWidth', 2); 
                h_true.DisplayName = sprintf('r=%4.1f | v=%4.1f', obj.true_point.r, obj.true_point.v); 
            end
            
            if input.best
                h_best = plot(ax_lc, t, f(:,1), '--r', 'LineWidth', 2); 
                h_best.DisplayName = sprintf('r=%4.1f | v=%4.1f', obj.best_point.r, obj.best_point.v); 
            end
            
            hold(ax_lc, 'off'); 
            
            ax_lc.FontSize = input.font_size;
            
            legend(ax_lc, 'Location', 'SouthEast'); 
            
            ax_lc.YAxisLocation='Right';
            ax_lc.XAxisLocation='Bottom';
%             xtickformat(ax_lc, '%ds'); 
            
            ax_lc.XLim = [-1 1.5].*1; 

            xlabel(ax_lc, 'time [s]'); 
            ylabel(ax_lc, 'Flux [normalized]'); 
            
        end
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = occult.gui.MCMC_GUI(obj); 
            end
            
            obj.gui.make;
            
        end
        
    end    
    
end

