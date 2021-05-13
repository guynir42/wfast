classdef MCMC < handle

    properties(Transient=true)
        
        gui;
        
    end
    
    properties % objects
        
        gen@occult.CurveGenerator; % used to create lightcurves to fit to the measured
        
        bank@occult.ShuffleBank; % for finding an initial guess
        
        init_point@occult.Parameters; % keep track of the initial point 
        true_point@occult.Parameters; % if the data is simulated and the true value is known
        best_point@occult.Parameters;

        points@occult.Parameters; % a chain of trial points
        
        prog@util.sys.ProgressBar; % print out the progress

    end
    
    properties % inputs/outputs
        
        input_flux; % the input lightcurve to be matched
        input_errors; % the rms error on the given flux
        
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
        
        step_sizes = [0.25, 0.25, 3]; % in order of parameters: step size for each parameter
        circ_bounds = [0 0 0]; % in order of parameters: which par gets a circular boundary condition
        
        par_list = {'r', 'b', 'v'}; % these parameters are chosen randomly each step. The rest are taken from the generator's parameters
        
%         use_bank = true; % use a filter bank (needs to be loaded manually!) to find an initial position (otherwise chose random starting point)
        
        initialization_method = 'search'; 

        use_priors = false; % a general switch to turn on/off the use of prior functions        
        prior_functions = {}; % can input a cell array with a different function per parameter (leave empty for uniform prior). min/max values are taken from generator
        
        plot_every = 1; % when plotting, how many steps go by between plots (set to 0 for no plotting)
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
        
        brake_bit = 1; % set to 0 when running, set to 1 to stop run
        
        % save when we run sampleGridLikelihood
        likelihood_map; 
        map_r;
        map_b;
        map_v;
        map_R; 
        
        version = 1.05;
        
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
                
%                 obj.true_point = occult.Parameters; 
                obj.best_point = occult.Parameters; 
                
                obj.prog = util.sys.ProgressBar;
                
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
            
            obj.init_point = occult.Parameters;
            obj.best_point = occult.Parameters;
%             obj.pick_point = occult.Parameters;
            
            obj.counter = 0; 
            obj.num_successes = 0;
            
            if ~isempty(obj.gui) && obj.gui.check
                cla(obj.gui.axes_posterior); 
            end
            
        end
        
    end
    
    methods % getters
        
        function val = getSummary(obj)
            
            val = ''; % to be continued
            
        end
        
        function val = getProgress(obj)
            
            val = sprintf('chain: %d / %d points | %s', obj.counter, obj.num_steps, obj.prog.show); 
            
        end
        
    end
    
    methods % setters
        
        function set.input_flux(obj, val)
            
            if ~isequal(obj.input_flux, val)
                
                obj.input_flux = val;
                
                obj.true_point = occult.Parameters.empty; 
                
                if ~isempty(obj.gui) && obj.gui.check
                    cla(obj.gui.axes_lightcurve); 
                end
                
                obj.likelihood_map = []; 
                obj.map_r = [];
                obj.map_b = [];
                obj.map_v = [];
                obj.map_R = []; 
                
            end
            
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
                obj.gen.lc.pars.likelihood(ii) = 0; 
    
            end 
            
            obj.init_point = occult.Parameters.empty;
            obj.init_point(1,obj.num_chains) = occult.Parameters;
            obj.init_point.copy_from(obj.gen.lc.pars); 
            
        end
                
        function [likelihood, r, b, v, R] = sampleGridLikelihood(obj, varargin)
            
            if isempty(obj.input_flux)
                error('Cannot run MCMC without an input flux! Use useSimulatedInput() or useSimulatedNoisyInput() or input a lightcurve manually.'); 
            end
            
            input = util.text.InputVars;
            input.input_var('star_size', 1, 'stellar size');
            input.input_var('radius_step', 0.25); 
            input.input_var('impact_step', 0.25); 
            input.input_var('velocity_step', 2.5); 
            input.scan_vars(varargin{:}); 
            
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
            
            likelihood = NaN(length(r), length(b), length(v)); 
            
            for ii = 1:length(r)
                
                obj.gen.r = r(ii); 
                
                L = obj.calcLikelihood; 
                
                likelihood(ii,:,:) = reshape(L, [length(b), length(v)]); 
                
            end
            
            obj.likelihood_map = likelihood;
            obj.map_r = r;
            obj.map_b = b;
            obj.map_v = v;
            obj.map_R = R; 
            
        end
        
        function points = searchGoodPoints(obj, varargin)
            
            if ~isempty(obj.input_errors)
                e = obj.input_errors;
            else
                e = 1./obj.gen.snr;
            end
            
            chi2 = nansum(((obj.input_flux - 1)/e).^2);
            dof = length(obj.input_flux) - 3; % assume parameters are r,b,v, only
            
            base_likelihood = chi2cdf(chi2, dof, 'upper'); 
            
            if isempty(obj.likelihood_map)
                [L, r, b, v, R] = obj.sampleGridLikelihood(varargin{:});
            else
                L = obj.likelihood_map; 
                r = obj.map_r;
                b = obj.map_b;
                v = obj.map_v;
                R = obj.map_R; 
            end
            
            idx = L>2*base_likelihood;
            
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

                        if subsref(obj.gen, S)>range(2)
                            delta = subsref(obj.gen, S) - range(2);
                            obj.gen = subsasgn(obj.gen, S, range(1) + delta); % circle back the rest of the way from the lower bound
                        elseif subsref(obj.gen, S)<range(1)
                            delta = range(1) - subsref(obj.gen, S);
                            obj.gen = subsasgn(obj.gen, S, range(2) - delta); % circle back the rest of the way from the upper bound
                        end

                    else % reflect back
                        
                        if subsref(obj.gen, S)>range(2)
                            delta = subsref(obj.gen, S) - range(2);
                            obj.gen = subsasgn(obj.gen, S, range(2) - delta); % reflect the next point off the upper edge of the range
                        elseif subsref(obj.gen, S)<range(1)
                            delta = range(1) - subsref(obj.gen, S);
                            obj.gen = subsasgn(obj.gen, S, range(1) + delta); % leave the point at the lower edge of the range
                        end

                    end 

                end % for ii (num_chains)

            end % for jj (par list)

        end
        
        function L = calcLikelihood(obj)
            
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
            
            for ii = 1:length(chi2)

                L = chi2cdf(chi2(ii), nnz(~isnan(f1(:,ii)))-length(obj.par_list), 'upper'); % do we need to subtract these parameters from the dof??
                
                if obj.use_priors
                    for jj = 1:length(obj.prior_functions)
                        if ~isempty(obj.prior_functions{jj})
                            values = obj.gen.(obj.par_list{jj});
                            % multiply by the value of the function at the given point in parameter space
                            L = L.*obj.prior_functions{jj}(values(ii)); 
                        end
                    end
                end
                
                obj.gen.lc.pars.likelihood(ii) = L;

                % likelihood is zero for parameters that are out of bounds
                for jj = 1:length(obj.par_list)

                    S = substruct('.', obj.par_list{jj}, '()', {ii});

                    value = subsref(obj.gen, S);
                    range = obj.gen.([obj.par_list{jj} '_range']); 

                    if value<range(1) || value>range(2)
                        obj.gen.lc.pars.likelihood(ii) = 0;
                        break; % no need to check other parameters
                    end

                end
                
            end
            
            if nargout>0
                L = obj.gen.lc.pars.likelihood;
            end

        end
        
        function run(obj, varargin)
            
            import util.text.cs;
            
            input = util.text.InputVars;
            input.input_var('plot', false);
            input.input_var('ax', [], 'axes', 'axis');
            input.scan_vars(varargin{:});
            
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

                    log_ratio = log(obj.gen.lc.pars.likelihood(ii)) - log(obj.points(jj-1,ii).likelihood);
                    
                    if isnan(log_ratio) || log(rand)<log_ratio % if ratio is big, we are likely to take the new point
                        obj.num_successes(ii) = obj.num_successes(ii) + 1;
                        obj.points(jj,ii).counts = 1;
                        
                        if obj.best_point.likelihood<obj.points(jj,ii).likelihood % keep track of the best fit point
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
                    
                    obj.plot('burn', 1, 'ax', input.ax); 
                    
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
        
        function [N, x, y] = calcPosterior(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('pars', obj.show_posterior_pars); 
            input.input_var('success_ratio', 0.2); % minimal number of successes (relative to total steps) for a chain to be included
            input.input_var('likelihood', 0.01); % minimal mean likelihood needed for a chain to be included
            input.input_var('weights', false, 'use_weights');
            input.input_var('burn', false); 
            input.scan_vars(varargin{:}); 
            
            if input.burn
                idx1 = 1:obj.counter; 
            else
                idx1 = obj.num_burned:obj.counter; 
            end
            
            for ii = 1:size(obj.points,2)
                
                if ~isempty(input.success_ratio)
                    idx2(ii) = obj.num_successes(ii)./obj.counter > input.success_ratio && ... 
                        nanmean([obj.points(idx1,ii).likelihood])> input.likelihood;
                end
                
            end
            
            x = [obj.points(idx1, idx2).(input.pars{1})];
            y = [obj.points(idx1, idx2).(input.pars{2})];
            w = logical([obj.points(idx1, idx2).weight]);
            
            if input.weights
                x = x(w);
                y = y(w); 
            end
            
            idx_x = find(strcmp(input.pars{1}, obj.par_list), 1, 'first');
            dx = obj.step_sizes(idx_x); 

            idx_y = find(strcmp(input.pars{2}, obj.par_list), 1, 'first');
            dy = obj.step_sizes(idx_y); 

            [N, Ex, Ey] = histcounts2(x, y, 'BinWidth', [dx, dy], 'Normalization','probability');
            N = N';
            x = repmat(Ex(1:end-1)+dx/2, [size(N,1),1]);                 
            y = repmat(Ey(1:end-1)'+dy/2, [1, size(N,2)]); 

            obj.posterior_prob = N;
            obj.post_pars = input.pars;
            obj.post_x_axis = x;
            obj.post_y_axis = y; 
            obj.post_options = input; 
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function plot(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('pars', obj.show_chain_pars); 
            input.input_var('burn', false); 
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('font_size', 18); 
            input.input_var('hold', false); 
            input.scan_vars(varargin{:});
            
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
                
                for ii = 1:N

                    p1 = [obj.points(1:obj.counter,ii).(input.pars{1})]';
                    p2 = [obj.points(1:obj.counter,ii).(input.pars{2})]';
                    p3 = [obj.points(1:obj.counter,ii).(input.pars{3})]';
                    lkl = [obj.points(1:obj.counter,ii).likelihood]';
                    
                    if isempty(h) || length(h)~=N || ~isvalid(h(ii))
                        N_colors = size(input.ax.ColorOrder,1);
                        col_idx = mod(ii-1, N_colors) + 1; 
                        sha_idx = floor((ii-1)/N_colors) + 1; 
                        h(ii) = scatter3(input.ax, p1, p2, p3, sz, input.ax.ColorOrder(col_idx,:), shapes{sha_idx}, 'filled');
                        h(ii).DisplayName = sprintf('Chain %d (acc=%3.1f%%, lkl=%4.2f)', ...
                            ii, 100*obj.num_successes(ii)./obj.counter, nanmean(lkl(obj.num_burned+1:end,1))); 
                        input.ax.NextPlot = 'add';                 
                    else
                        h(N-ii+1).XData = p1; 
                        h(N-ii+1).YData = p2;
                        h(N-ii+1).ZData = p3;
                        h(N-ii+1).SizeData = sz;
                        h(N-ii+1).DisplayName = sprintf('Chain %d (acc=%3.1f%%, lkl=%4.2f)', ...
                            ii, 100*obj.num_successes(ii)./obj.counter, nanmean(lkl(obj.num_burned+1:end,1),1)); 
                        input.ax.NextPlot = 'add';  
                    end

                end
                
                input.ax.NextPlot = 'replace'; 
                
                % the X limits and labels are added in the end
%                 ylabel(input.ax, obj.title_strings.(input.pars{2}));
                ylabel(input.ax, input.pars{2});
                input.ax.YLim = obj.gen.([input.pars{2} '_range']);
                
%                 zlabel(input.ax, obj.title_strings.(input.pars{3}));
                zlabel(input.ax, input.pars{3});
                input.ax.ZLim = obj.gen.([input.pars{3} '_range']);
                
%                 input.ax.CLim = [0 1]; 
%                 hcb = colorbar(input.ax); 
%                 ylabel(hcb, 'likelihood'); 
                
            end
        
%             xlabel(input.ax, obj.title_strings.(input.pars{1}));
            xlabel(input.ax, input.pars{1});
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
            
            N = min(obj.num_chains, obj.show_num_chains); 
            
            if isempty(h) || ~isvalid(h(1)) || N~=length(h)-1
                
                delete(input.ax.Children); 
                
                plot(input.ax, t, obj.input_flux, 'dk', 'LineWidth', 3, 'DisplayName', 'Input flux'); 
                
                h = obj.gen.lc.plot('ax', input.ax, 'hold', true, 'legend', true, 'noise', false, 'FontSize', input.font_size, 'number', N);
                
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
            
            chi2 = obj.gen.lc.pars.chi2;
            
            dof = nnz(~isnan(obj.gen.lc.flux(:,1)))-length(obj.par_list);
            
            text(input.ax, t(10), 0.3, sprintf('chi2 (dof=%d) = %s', ...
                dof, print_vec(round(chi2(1:N)), ', ')), 'FontSize', input.font_size); 
            
            text(input.ax, t(10), 0.1, sprintf('chi2/dof = %s', ...
                print_vec(round(chi2(1:N)./dof,1), ', ')), 'FontSize', input.font_size); 
            
            text(input.ax, t(10), -0.1, sprintf('likelihood= %s', ...
                print_vec(obj.gen.lc.pars.likelihood(1:N))), 'FontSize', input.font_size); 
            
        end
        
        function showPosterior(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('pars', obj.show_posterior_pars); 
            input.input_var('burn', false); 
            input.input_var('ax', [], 'axes', 'axis');
            input.input_var('font_size', 18); 
            input.input_var('hold', false); 
            input.scan_vars(varargin{:});
            
            if isempty(input.ax)
                
                if ~isempty(obj.gui) && obj.gui.check
                    input.ax = obj.gui.axes_posterior; 
                else
                    input.ax = gca;
                end
                
            end
            
            if isempty(input.pars)
                error('Must specify which parameters to show!'); 
            elseif isscalar(input.pars)
                histogram(input.ax, [obj.points(1:obj.num_burned).(input.pars{1})]); 
            elseif length(input.pars)==2
                
                [N,x,y] = obj.calcPosterior(varargin{:}); 

                [M,c] = contour(input.ax, x, y, N, '-k'); 
                
            else % more than 2 pars
                error('Must specify one or two parameters!'); 
            end
            
            util.plot.inner_title(input.pars{1}, 'Position', 'Bottom', 'ax', input.ax); 
            util.plot.inner_title(input.pars{2}, 'Position', 'Left', 'ax', input.ax); 
            util.plot.inner_title('Posterior', 'Position', 'Top', 'ax', input.ax); 

            input.ax.FontSize = input.font_size; 

            if ~isempty(obj.true_point)
                
                hold(input.ax, 'on'); 
                
                if ~isempty(obj.best_point)
                    plot(input.ax, obj.best_point.(input.pars{1}), obj.best_point.(input.pars{2}), 'og', 'MarkerSize', 10, 'MarkerFaceColor', 'g'); 
                end
                
                plot(input.ax, obj.true_point.(input.pars{1}), obj.true_point.(input.pars{2}), '+r', 'MarkerSize', 20); 
                plot(input.ax, obj.true_point.(input.pars{1}), obj.true_point.(input.pars{2}), 'or', 'MarkerSize', 20); 
                
                hold(input.ax, 'off'); 
                
            end
            
        end
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = occult.gui.MCMC_GUI(obj); 
            end
            
            obj.gui.make;
            
        end
        
    end    
    
end

