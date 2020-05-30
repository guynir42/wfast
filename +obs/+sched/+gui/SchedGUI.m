classdef SchedGUI < handle
    
    properties 
        
        owner@obs.sched.Scheduler; % link back to containg object

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
        
        panel_current; 
        panel_controls;        
        panel_display;
        panel_sim;
        panel_contrast;
    
        panel_report;
        
        panel_filename;
        
        panel_close;
        button_close;
        
        panel_image;
        button_reset_axes;
        axes_image;
    
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = SchedGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'obs.sched.Scheduler')
                
                if obj.debug_bit>1, fprintf('SchedGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input an obs.sched.Scheduler to constructor of SchedGUI!');
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
            
            obj.fig = util.plot.FigHandler('Scheduler');
            obj.fig.clear;
            obj.fig.height = 20;
            obj.fig.width = 36;
            movegui(obj.fig.fig, 'center');
            
            obj.fig.bottom = 3;
            
            %%%%%%%%%%%%%%%%%%%%%%% MENUS %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % MenuItem(parent, text, type, variable, tooltip, separator)
            % obj.addButton(name, text, type, variable, tooltip, separator)
            % menu types: menu, toggle, push, input, input_text, info, custom
            
            
            
            N = 12; % number of buttons on left side
            
            %%%%%%%%%%%%%%%%%%% LEFT SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            pos = N;
            
            %%%%%%%%%%% panel current %%%%%%%%%%%%%%%
            
            num_buttons = 2;
            pos = pos-num_buttons;
            obj.panel_current = GraphicPanel(obj.owner, [0 pos/N 0.2 num_buttons/N], 'current state', 1); % last input is for vertical (default)
            obj.panel_current.number = num_buttons;
            obj.panel_current.addButton('button_ra', 'current_RA', 'info', 'RA= ', '', 'edit', 0.5, '', '', 'right ascention last object to be observed');
            obj.panel_current.addButton('button_dec', 'current_Dec', 'info', 'Dec= ', '', 'edit', 0.5, '', '', 'declination last object to be observed');
            obj.panel_current.addButton('button_side', 'current_side', 'info', 'side= ', '', '', 0.5, '', '', 'side (hemisphere) of last object to be observed');
            obj.panel_current.addButton('button_wind', 'wind_string', 'info', 'wind= ', '', 'edit', 0.5, '', '', 'wind state right now (from Manager or from simulation)');
            obj.panel_current.margin = [0.02 0.05]; 
            obj.panel_current.make;
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            % Add buttons using obj.addButton(button_name, var_name='', type='', str1='', str2='', font_size='', split=1, color_on=[], color_off=[], tooltip)
            
            num_buttons = 3;
            pos = pos-num_buttons;
            obj.panel_controls = GraphicPanel(obj.owner, [0 pos/N 0.2 num_buttons/N], 'controls', 1); % last input is for vertical (default)
            obj.panel_controls.number = num_buttons;
            obj.panel_controls.addButton('button_stay', 'use_stay_on_side', 'toggle', 'any side', 'stay on side', 'edit', 0.5, obj.color_on, '', 'choose if the telescope should stay on the same side all night'); 
            obj.panel_controls.addButton('input_wind', 'max_wind', 'input', 'max wind= ', 'km/h', 'edit', 0.5, '', '', 'maximum wind speed to observe towards the west'); 
            obj.panel_controls.addButton('button_constraints', 'ephem.constraints', 'push', 'constraints', '', '', 0.5, '', '', 'choose the default constraints to apply BEFORE reading the target list'); 
            obj.panel_controls.margin = [0.02 0.05];
            obj.panel_controls.make;
            
            %%%%%%%%%%% panel display %%%%%%%%%%%%%%%
            
            num_buttons = 3;
            pos = pos-num_buttons;
            obj.panel_display = GraphicPanel(obj.owner, [0 pos/N 0.2 num_buttons/N], 'display', 1); % last input is for vertical (default)
            obj.panel_display.number = num_buttons;
            obj.panel_display.addButton('button_ecliptic', 'map.show_ecliptic', 'toggle', 'ecliptic', 'ecliptic', '', 0.5, 'red', '', 'show the ecliptic altitude overlay'); 
            obj.panel_display.addButton('input_altitude', 'map.alt_limit', 'input', 'alt_lim= ', '', '', 0.5, '', '', 'alt limit for display on the map only!'); 
            obj.panel_display.margin = [0.02 0.05];
            obj.panel_display.make;
            
            
            %%%%%%%%%%% panel simulations %%%%%%%%%%%%%%%
            
            num_buttons = 3;
            pos = pos-num_buttons;
            obj.panel_sim = GraphicPanel(obj.owner, [0 pos/N 0.2 num_buttons/N], 'simulations', 1); % last input is for vertical (default)
            obj.panel_sim.number = num_buttons;
            obj.panel_sim.addButton('input_step', 'sim_time_step', 'input', 'step= ', 'min', '', 0.5, '', '', 'time interval for checking out new targets'); 
            obj.panel_sim.addButton('input_side', 'sim_starting_side', 'input_text', 'start= ', '', '', 0.5, '', '', 'side on which the telescope begins the simulation'); 
            % add wind start/end time here if you want...
            obj.panel_sim.addButton('input_pause', 'sim_plot_pause', 'input', 'pause= ', 's', '', 0.5, '', '', 'additional pause time between simulation steps, for slowing down the plotting'); 
            obj.panel_sim.addButton('button_resume', 'run_simulation', 'push', 'continue', '', '', 0.5, '', '', 'continue the simulation from where it was stopped'); 
            obj.panel_sim.addButton('button_start', '', 'custom', 'Start simulation', '', '', 1, '', '', 'start or stop the simulation'); 
            obj.panel_sim.margin = [0.02 0.05];
            obj.panel_sim.make;
            
            obj.panel_sim.button_start.Callback = @obj.callback_start_stop;
            
            %%%%%%%%%%% panel contrast %%%%%%%%%%%%%%%
            
%             pos = pos - 5; 
%             obj.panel_contrast = util.plot.ContrastLimits(obj.axes_image, obj.fig.fig, [0 pos/N 0.2 5/N], 1); % last input is for vertical (default)
%             obj.panel_contrast.font_size = obj.font_size;
%             obj.panel_contrast.big_font_size = obj.big_font_size;
%             obj.panel_contrast.small_font_size = obj.small_font_size;
%             obj.panel_contrast.edit_font_size = obj.edit_font_size;
            
            
            %%%%%%%%%%% panel report %%%%%%%%%%%%%%%%%
            
            obj.panel_report = GraphicPanel(obj.owner, [0.2 (N-1)/N 0.8 1/N], 'report', 1); % last input is for vertical (default)
            obj.panel_report.addButton('button_report', 'report', 'info'); 
            obj.panel_report.margin = [0.02 0.01];
            obj.panel_report.make;
            
            %%%%%%%%%%% panel filename %%%%%%%%%%%%%%%
            
            obj.panel_filename = GraphicPanel(obj.owner, [0.2 0 0.8 1/N], 'report', 1);  % last input is for vertical (default)
            obj.panel_filename.addButton('input_filename', 'filename', 'input', '', '', '', 0.6, '', '', 'give the filename to load the target list from'); 
            obj.panel_filename.addButton('button_read', 'readFile', 'push', 'readFile', '', '', 0.2, '', '', 'read the targets from current file'); 
            obj.panel_filename.addButton('button_browse', '', 'custom', 'browse', '', '', 0.2, '', '', 'choose a new file for the target list'); 
            obj.panel_filename.margin = [0.01 0.02];
            obj.panel_filename.make;
            
            obj.panel_filename.button_browse.Callback = @obj.callback_browse;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 1/N 0.8 (N-2)/N]); 
                        
            obj.makeAxes;
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.85 0.95 0.15 0.05], obj.owner, '', 'custom', 'new axes', '');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            obj.button_reset_axes.Tooltip = 'Create a new image axis, zoomed out and with default contrast limits'; 
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Position', [0 0 0.2 1/N]);
            obj.button_close = GraphicButton(obj.panel_close, [0.05 0.05 0.9 0.9], obj.owner, '', 'custom', 'CLOSE GUI');
            obj.button_close.Callback = @obj.callback_close;
            
            obj.update;
            
        end
            
        function makeAxes(obj, ~, ~)
            
            delete(obj.axes_image);
            
            obj.axes_image = axes('Parent', obj.panel_image);
            
            obj.panel_contrast.ax = obj.axes_image;
%             colorbar(obj.axes_image);
            axis(obj.axes_image, 'image');
            
            obj.panel_contrast.ax = obj.axes_image;
            
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
                obj.panel_sim.button_start.String = 'Start simulation';
            else
                obj.panel_sim.button_start.String = 'Stop simulation';
            end
            
            obj.owner.show;
            
%             obj.panel_contrast.update;
            
            drawnow;
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_start_stop(obj, ~, ~)
            
            if obj.debug_bit>1, disp('callback: start/stop'); end
            
            if obj.owner.brake_bit
                obj.owner.run_simulation(1); % this first argument is used to start a new simulation 
            else
                obj.owner.brake_bit = 1;
            end
            
            obj.update;
            
        end
   
        function callback_browse(obj, ~, ~)
            
            if obj.debug_bit>1, disp('callback: browse'); end
            
            [filename, path] = uigetfile('*.txt', 'Choose a new target list', fullfile(getenv('DATA'), 'WFAST/target_lists/target_list.txt')); 
            
            if ischar(filename)
                obj.owner.filename = fullfile(path, filename);
                obj.owner.readFile;
            end
            
            obj.update;
            
        end
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end