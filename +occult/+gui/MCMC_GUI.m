classdef MCMC_GUI < handle
    
    properties 
        
        owner@occult.MCMC; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        menus = {};
        panels = {};
        
        latest_error = '';
        latest_warning = '';
        
        font_size = 12;
        big_font_size = 16;
        edit_font_size = 11;
        small_font_size = 10;
        
        color_on = [0 0.3 1];
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        menu_options;
        
        panel_controls;
        panel_simulations;
        panel_display; 
        
        panel_progress;
        
        panel_chain;
        axes_chain;        
        axes_posterior;
        
        panel_lightcurve; 
        axes_lightcurve;
        
        panel_close;
        button_close;
        
    end
    
    properties (Hidden=true)
        
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = MCMC_GUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'occult.MCMC')
                
                if obj.debug_bit>1, fprintf('MCMC_GUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input an occult.MCMC to constructor of MCMC_GUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            import util.plot.MenuItem;
            
            obj.buttons = {};
            obj.menus = {};
            obj.panels = {}; 
            
            obj.fig = util.plot.FigHandler('MCMC GUI');
            obj.fig.clear;
            obj.fig.bottom = 5;
            obj.fig.height = 18;
            obj.fig.width = 36;
            obj.fig.center;
%             obj.fig.maximize;
            
            
            %%%%%%%%%%%%%%%%%%%%%%% MENUS %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % MenuItem(parent, text, type, variable, tooltip, separator)
            % obj.addButton(name, text, type, variable, tooltip, separator)
            % menu types: menu, toggle, push, input, input_text, info, custom
            
            % NOTE: to enable tooltips, run obj.menu_XXX.assignJavaObjectsTopLevel 
            %       on each top-level menu object, after adding all items.
            
            
            obj.menu_options = MenuItem(obj, '&Options', 'menu'); 
            obj.menu_options.addButton('menu_setup', '&Setup', 'menu'); 
            obj.menu_options.menu_setup.addButton('button_gen', '&Generator GUI', 'push', 'gen'); 
            obj.menu_options.menu_setup.addButton('button_sim', '&Use simulated', 'push', 'useSimulatedInput'); 
            obj.menu_options.menu_setup.addButton('button_noisy', 'Use &Noisy sim', 'push', 'useSimulatedNoisyInput'); 
            
            
            %%%%%%%%%%%%%%%%%%% LEFT SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N = 12; % number of buttons on left side
            
            pos = N;
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            % Add buttons using obj.addButton(button_name, var_name='', type='', str1='', str2='', font_size='', split=1, color_on=[], color_off=[], tooltip)
            
            num_buttons = 7;
            pos = pos-num_buttons;
            obj.panel_controls = GraphicPanel(obj.owner, [0 pos/N 0.2 num_buttons/N], 'controls');
            obj.panel_controls.addButton('button_run', 'run', 'custom', 'RUN', '', '', [], '', '', 'start a new MCMC run'); 
            obj.panel_controls.addButton('button_num_chains', 'num_chains', 'input', '', '', '', [], '', '', 'number of chains to run in parallel'); 
            obj.panel_controls.addButton('button_num_steps', 'num_steps', 'input', '', '', '', [], '', '', 'number of steps in the chain'); 
            obj.panel_controls.addButton('button_num_burn', 'num_burned', 'input', '', '', '', [], '', '', 'number of steps to burn in the beginning of the chain');
            obj.panel_controls.addButton('button_step_sizes', 'step_sizes', 'input', '', '', '', [], '', '', 'the step size for each parameter'); 
            obj.panel_controls.addButton('button_circ_bounds', 'circ_bounds', 'input', '', '', '', [], '', '', 'determine which parameter gets a circular boundary condition'); 
            obj.panel_controls.addButton('button_priors', 'use_bank', 'toggle', 'random start', 'using bank', '', 0.5, obj.color_on, '', 'use a random starting position or a template bank best kernel'); 
            obj.panel_controls.addButton('button_priors', 'use_priors', 'toggle', 'ignoring priors', 'using priors', '', 0.5, obj.color_on, '', 'apply the prior functions on some of the parameters'); 
            obj.panel_controls.number = num_buttons;
            
            obj.panel_controls.make;
            obj.panel_controls.button_run.Callback = @obj.callback_run; 
            
            %%%%%%%%%%% panel simulations %%%%%%%%%%%%%%%%%%
            
            num_buttons = 2;
            pos = pos-num_buttons;
            obj.panel_simulations = GraphicPanel(obj.owner, [0 pos/N 0.2 num_buttons/N], 'controls');
            obj.panel_simulations.addButton('button_gen', 'gen', 'push', 'Generator GUI', '', '', 0.8, '', '', 'show the generator GUI'); 
            obj.panel_simulations.addButton('button_use_lc', 'useSimulatedNoisyInput', 'push', 'input', '', '', 0.2, '', '', 'use the current generator LC as input'); 
            obj.panel_simulations.addButton('button_random', '', 'custom', 'random seed', '', '', [], '', '', 'pick random pars for generator'); 
            obj.panel_simulations.make;
            
            obj.panel_simulations.button_random.Callback = @obj.callback_random; 
            
            %%%%%%%%%%% panel display %%%%%%%%%%%%%%%%%%
            
            num_buttons = 2; 
            pos = pos - num_buttons;
            obj.panel_display = GraphicPanel(obj.owner, [0 pos/N 0.2 num_buttons/N], 'display'); 
            obj.panel_display.addButton('button_every', 'plot_every', 'input', 'plot every= ', '', '', [], '', '', 'update the plot every N steps'); 
            obj.panel_display.addButton('button_chains', 'show_num_chains', 'input', 'show= ', ' chains', '', [], '', '', 'how many (max) number of chains to display'); 
            
            obj.panel_display.number = num_buttons;
            obj.panel_display.make;
            
            %%%%%%%%%%% panel info %%%%%%%%%%%%%%%%%%
            
            obj.panel_progress = GraphicPanel(obj.owner, [0.2 (N-1)/N 0.8 1/N], 'progress'); 
            obj.panel_progress.addButton('button_info', 'getProgress', 'info'); 
            obj.panel_progress.make;
            
            %%%%%%%%%%% panel chain %%%%%%%%%%%%%%%%%%
            
            obj.panel_chain = GraphicPanel(obj.owner, [0.2 (N-7)/N 0.8 6/N], 'chain'); 
            obj.panel_chain.make;
            
            obj.axes_chain = axes('Parent', obj.panel_chain.panel, 'Position', [0.15 0.3 0.7 0.7]); 
            
            obj.axes_posterior = axes('Parent', obj.panel_chain.panel, 'Position', [0.8 0.1 0.18 0.4]); 
            
            %%%%%%%%%%% panel lightcurve %%%%%%%%%%%%%%%%%%
            
            obj.panel_lightcurve = GraphicPanel(obj.owner, [0.2 0 0.8 3/N], 'lightcurve'); 
            obj.panel_lightcurve.make; 
            
            obj.axes_lightcurve = axes('Parent', obj.panel_lightcurve.panel); 
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Position', [0 0 0.2 1/N]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE GUI');
            obj.button_close.Callback = @obj.callback_close;
            
            obj.update;
            
        end
        
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            for ii = 1:length(obj.menus)
                obj.menus{ii}.update;
            end
            
            if obj.owner.brake_bit
                obj.panel_controls.button_run.String = 'RUN'; 
            else
                obj.panel_controls.button_run.String = 'STOP'; 
            end
            
        end
        
        function updateLightcurves(obj)
            
            obj.owner.plotLightcurves('ax', obj.axes_lightcurve);
            
        end
        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_run(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: run'); end
            
            if obj.owner.brake_bit
                obj.owner.run('plot', 1, 'ax', obj.axes_chain); 
            else
                obj.owner.brake_bit = 1; 
            end
            
            obj.update; 
            
        end
        
        function callback_random(obj, ~, ~)
            
            if obj.debug_bit>1, disp('callback: random'); end
            
            obj.owner.getRandomLightcurve; 
            obj.owner.useSimulatedNoisyInput;             
            obj.updateLightcurves; 
            
        end
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end








