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
        num_failures = 0; % number of steps that failed
        
        
    end
    
    properties % switches/controls
        
        num_steps = 1000; % total number of steps (including burn-in)
        num_burned = 100; % number of steps to burn at the begining of each chain
        
        step_sizes = [0.1, 0.1, 1]; % in order of parameters: step size for each parameter
        circ_bounds = [0 1 0]; % in order of parameters: which par gets a circular boundary condition
        
        par_list = {'r', 'b', 'v'}; % these parameters are chosen randomly each step. The rest are taken from the generator's parameters
        
        use_bank = true; % use a filter bank (needs to be loaded manually!) to find an initial position (otherwise chose random starting point)
        
        use_priors = false; % a general switch to turn on/off the use of prior functions        
        prior_functions = {}; % can input a cell array with a different function per parameter (leave empty for uniform prior). min/max values are taken from generator
        
        plot_every = 1; % when plotting, how many steps go by between plots (set to 0 for no plotting)
        
        show_chain_pars = {'r', 'b', 'b'}; % which parameters are shown, by default, when plotting the chain
        show_posterior_pars = {'r', 'v'}; % which parameters are shown, by default, when plotting posteriors
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)

        brake_bit = 1; % set to 0 when running, set to 1 to stop run
        
        version = 1.03;
        
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
                
            end
            
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
            obj.num_failures = 0;
            
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
            obj.input_flux = obj.gen.lc.flux_noisy;
            
            obj.true_point = occult.Parameters;
            obj.true_point.copy_from(obj.gen.lc.pars); 
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.updateLightcurves; 
            end
            
        end
        
        function startup(obj)
            
            if length(obj.gen.r)>1
                error('Must setup generator with all scalar parameters!');
            end
            
            obj.points(obj.num_steps,1) = occult.Parameters; 
            
            if ~isempty(obj.gui) && obj.gui.check
                cla(obj.gui.axes_chain);
                cla(obj.gui.axes_lightcurve);
            end
            
        end
        
        function finishup(obj)
            
            obj.points = obj.points(1:obj.counter); 
            obj.prog.finish;
            
            obj.brake_bit = 1; 
            
        end
        
        function gotoBestTemplateInitialPoint(obj)
            
            if isempty(obj.bank)
                error('Must load a template bank first!'); 
            end
            
            f = obj.input_flux - 1;
            
            ff = util.vec.convolution(obj.bank.kernels, f, 'cross', 1); % filtered flux
            
            [~,idx] = util.stat.max2(abs(ff)); % find the peak in the filtered flux
            
            pars = obj.bank.pars(idx(2)); % grab a struct with some parameters
            
            obj.init_point.copy_from(obj.gen.lc.pars); % start by getting all parameters from generator
            
            for ii = 1:length(obj.par_list)
                
                name = obj.par_list{ii};
                
                if isfield(pars, name)                    
                    obj.init_point.(name) = pars.(name); 
                end
                
            end
            
            if ismember('t', obj.par_list)
                obj.init_point.t = 0; % the template banks all have t=0 
            end
            
            obj.gen.lc.pars.copy_from(obj.init_point); 
            
        end
        
        function gotoRandomInitialPoint(obj)
            
            for ii = 1:length(obj.par_list)
                
                name = obj.par_list{ii}; 
                range = obj.gen.([name '_range']); 
                
                obj.gen.(name) = diff(range).*rand + range(1); 
                
            end
            
            obj.init_point.copy_from(obj.gen.lc.pars); 
            
        end
        
        function takeStep(obj)
            
            for ii = 1:length(obj.par_list)
                
                name = obj.par_list{ii}; 
                range = obj.gen.([name '_range']); 
                obj.gen.(name) = normrnd(obj.gen.(name), obj.step_sizes(ii)); 
                
                if obj.circ_bounds(ii)
                    
                    if obj.gen.(name)>range(2)
                        delta = obj.gen.(name) - range(2);
                        obj.gen.(name) = range(1) + delta; % circle back the rest of the way from the lower bound
                    elseif obj.gen.(name)<range(1)
                        delta = range(1) - obj.gen.(name);
                        obj.gen.(name) = range(2) - delta; % circle back the rest of the way from the upper bound
                    end
                    
                end
                
            end
            
        end
        
        function calcLikelihood(obj)
            
            for ii = 1:length(obj.par_list)
                
                value = obj.gen.lc.pars.(obj.par_list{ii});
                range = obj.gen.([obj.par_list{ii} '_range']);                
                
                if value<range(1) || value>range(2)
                    obj.gen.lc.pars.likelihood = 0; 
                    return;
                end
                
            end
            
            obj.gen.getLightCurves; 
            
            if ~isempty(obj.input_errors)
                e = obj.input_errors;
            else
                e = 1./obj.gen.snr;
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

            if obj.use_bank && ~isempty(obj.bank)
                obj.gotoBestTemplateInitialPoint; 
            else
                obj.gotoRandomInitialPoint;
            end
            
            obj.calcLikelihood; % get the first point likelihood/chi2

            obj.points(1).copy_from(obj.gen.lc.pars); % automatically accept first point... 
            
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
                
                ratio = obj.gen.lc.pars.likelihood./obj.points(jj-1).likelihood; % compare the new potential point with the last point

                if rand<=ratio % if ratio is big, we are likely to take the new point
                    obj.num_successes = obj.num_successes + 1;
                    obj.points(jj).copy_from(obj.gen.lc.pars);
                    obj.points(jj).counts = 1;
                    
                    if obj.best_point.likelihood<obj.points(jj).likelihood % keep track of the best fit point
                        obj.best_point.copy_from(obj.points(jj)); 
                    end
                    
                else % point is rejected, repeat the previous point
                    obj.num_failures = obj.num_failures + 1;
                    obj.points(jj).copy_from(obj.points(jj-1));
                    obj.points(jj).counts = obj.points(jj-1).counts + 1;
                end

                obj.gen.lc.pars.copy_from(obj.points(jj)); 
                
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
            
            lkl = [obj.points.likelihood]';
            
            p1 = [obj.points.(input.pars{1})]';
            p1 = p1(1:size(lkl,1),:); 
            
            sz = ones(size(lkl)); 
            
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
                
                p2 = [obj.points.(input.pars{2})]';
                p2 = p2(1:size(lkl,1),:); 
                
                p3 = [obj.points.(input.pars{3})]';
                p3 = p3(1:size(lkl,1),:); 
                
                h = findobj(input.ax, 'Type', 'Scatter'); 
                
                if isempty(h) || ~isvalid(h)
                    h = scatter3(input.ax, p1, p2, p3, sz, lkl, 'o', 'filled');
                else
                    h.XData = p1; 
                    h.YData = p2;
                    h.ZData = p3;
                    h.SizeData = sz;
                    h.CData = lkl;
                end

                % the X limits and labels are added in the end
                ylabel(input.ax, input.pars{2});
                input.ax.YLim = obj.gen.([input.pars{2} '_range']);
                
                zlabel(input.ax, input.pars{3});
                input.ax.ZLim = obj.gen.([input.pars{3} '_range']);
                
                input.ax.CLim = [0 1]; 
                hcb = colorbar(input.ax); 
                ylabel(hcb, 'likelihood'); 
                
            end
        
            xlabel(input.ax, input.pars{1});
            input.ax.FontSize = input.font_size;
            
            input.ax.XLim = obj.gen.([input.pars{1} '_range']);
            
            input.ax.NextPlot = hold_state;
            
            if ~isempty(h)
                h.ButtonDownFcn = @obj.callback_click_point; 
            end
            
        end
        
        function callback_click_point(obj, hndl, ev)
            
            if obj.debug_bit, fprintf('clicked point at '); end
            
            for ii = 1:length(obj.par_list)
                
                obj.gen.lc.pars.(obj.par_list{ii}) = ev.IntersectionPoint(ii);
                
                if obj.debug_bit, fprintf('%s= %4.2f ', obj.par_list{ii}, ev.IntersectionPoint(ii)); end
                
            end
            
            obj.calcLikelihood; 
            
            if obj.debug_bit, fprintf('lkl= %4.2f\n', obj.gen.lc.pars.likelihood); end
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gen.lc.plot('ax', obj.gui.axes_lightcurve); 
            end
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.updateLightcurves; 
                obj.gui.update;
            end
            
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
                x = [obj.points(obj.num_burned+1:end).(input.pars{1})];
                y = [obj.points(obj.num_burned+1:end).(input.pars{2})];
                
                idx_x = find(strcmp(input.pars{1}, obj.par_list), 1, 'first');
                dx = obj.step_sizes(idx_x); 
                
                idx_y = find(strcmp(input.pars{2}, obj.par_list), 1, 'first');
                dy = obj.step_sizes(idx_y); 
                
                [N, Ex, Ey] = histcounts2(x, y, 'BinWidth', [dx, dy], 'Normalization','probability');
                N = N';
                x = repmat(Ex(1:end-1)+dx/2, [size(N,1),1]);                 
                y = repmat(Ey(1:end-1)'+dy/2, [1, size(N,2)]); 
                
                [M,c] = contour(input.ax, x, y, N, '-k'); 
                
            else % more than 2 pars
                error('Must specify one or two parameters!'); 
            end
            
            xlabel(input.ax, input.pars{1}); 
            ylabel(input.ax, input.pars{2}); 

            input.ax.FontSize = input.font_size; 

            if ~isempty(obj.true_point)
                
                hold(input.ax, 'on'); 
                
                plot(input.ax, obj.true_point.(input.pars{1}), obj.true_point.(input.pars{2}), '+r', 'MarkerSize', 20); 
                plot(input.ax, obj.true_point.(input.pars{1}), obj.true_point.(input.pars{2}), 'or', 'MarkerSize', 20); 
                
                hold(input.ax, 'off'); 
                
            end
            
        end
        
        function plotLightcurves(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
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
            
            if isempty(h) || ~isvalid(h(1))
                delete(input.ax.Children); 
                plot(input.ax, t, obj.input_flux, 'dk', 'LineWidth', 3, 'DisplayName', 'Input flux'); 
                
                obj.gen.lc.plot('ax', input.ax, 'hold', true, 'legend', true, 'noise', false, 'FontSize', input.font_size);
                
            else
                h = h(1); % recover only the first plot (which is the model fit)
                h.XData = t; 
                h.YData = f; 
            end
            
            delete(findobj(input.ax, 'Type', 'Text')); 
            
            title(input.ax, ''); 
            input.ax.YLim = [0 1.2]; 
            
            chi2 = obj.gen.lc.pars.chi2;
            dof = nnz(~isnan(obj.gen.lc.flux))-length(obj.par_list);
            
            text(input.ax, t(10), 0.6, sprintf('chi2 / dof = %4.2f / %d = %4.2f', ...
                chi2, dof, chi2/dof), 'FontSize', input.font_size); 
            
            text(input.ax, t(10), 0.4, sprintf('likelihood= %g', ...
                obj.gen.lc.pars.likelihood), 'FontSize', input.font_size); 
            
        end
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = occult.gui.MCMC_GUI(obj); 
            end
            
            obj.gui.make;
            
        end
        
    end    
    
end

