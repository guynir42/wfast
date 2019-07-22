classdef FinderGUI < handle
    
    properties 
        
        owner@trig.Finder; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 12;
        big_font_size = 16;
        edit_font_size = 10;
        small_font_size = 8;
        
        color_on = [0 0.3 1];
        
        brake_bit = 1;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_threshold;
        
        panel_rejection;
        
        panel_memory;
        
        panel_statistics;
        
        panel_close;
        button_close;
        
        panel_obs_pars;
        
        panel_phot_pars;
        
        panel_event_pars;
        
        panel_image;
        
        panel_chooser;
        
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = FinderGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'trig.Finder')
                
                if obj.debug_bit, fprintf('FinderGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input a trig.Finder to constructor of FinderGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            
            obj.fig = util.plot.FigHandler('event finder');
            obj.fig.clear;
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 25;
            movegui(obj.fig.fig, 'center');
            
            N_left = 15; % number of buttons on left side
            
            %%%%%%%%%%%%%%%%%%% LEFT SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            pos = N_left;
            
            %%%%%%%%%%% panel threshold %%%%%%%%%%%%%%%
            
            % Add buttons using obj.addButton(button_name, var_name='', type='', str1='', str2='', font_size='', split=1, color_on=[], color_off=[], tooltip)
            
            num_buttons = 5;
            pos = pos-num_buttons;
            obj.panel_threshold = GraphicPanel(obj.owner, [0 pos/N_left 0.2 num_buttons/N_left], 'threshold', 1); % last input is for vertical (default)
            obj.panel_threshold.addButton('input_star_snr', 'min_star_snr', 'input', 'star S/N= ', '', '', 1, [], [], 'Stars with lower S/N thant this are not tested for events');
            obj.panel_threshold.addButton('input_threshold', 'threshold', 'input', 'threshold= ', '', '', 1, [], [], 'Threshold for individual events');
            obj.panel_threshold.addButton('input_time_thresh', 'time_range_thresh', 'input', 'time thresh= ', '', '', 1, [], [], 'Area of lightcurve around peak that is still considered part of the event');
            obj.panel_threshold.addButton('input_kern_thresh', 'kern_range_thresh', 'input', 'kern thresh= ', '', '', 1, [], [], 'Other kerenls around best one that are above this threshold are counted in the same event');
            obj.panel_threshold.addButton('input_star_thresh', 'star_range_thresh', 'input', 'star thresh= ', '', '', 1, [], [], 'Other stars around best one that are above this threshold are counted in the same event');
            obj.panel_threshold.number = num_buttons;
            obj.panel_threshold.margin = [0.05 0.02];
            obj.panel_threshold.make;
            
            %%%%%%%%%%% panel rejection %%%%%%%%%%%%%%%
            
            num_buttons = 5;
            pos = pos-num_buttons;
            obj.panel_rejection = GraphicPanel(obj.owner, [0 pos/N_left 0.2 num_buttons/N_left], 'rejection parameters', 1); % last input is for vertical (default)
            obj.panel_rejection.addButton('input_max_events', 'max_events', 'input', 'max events= ', '', '', 1, [], [], 'How many events can be found in one time-window');
            obj.panel_rejection.addButton('input_max_stars', 'max_stars', 'input', 'max stars= ', '', '', 1, [], [], 'How may stars can trigger at the same time');
            obj.panel_rejection.addButton('input_max_frames', 'max_frames', 'input', 'max frames= ', '', '', 1, [], [], 'How many frames can the event last, before being discarded because it is too long');
            obj.panel_rejection.addButton('input_max_nans', 'max_num_nans', 'input', 'max NaNs= ', '', '', 1, [], [], 'How many NaN values in the flux/centroid values can we accept in the event before it is disqualified');
            obj.panel_rejection.addButton('input_max_corr', 'max_corr', 'input', 'max corr= ', '', '', 1, [], [], 'Highest correlation between flux and background/aperture area before disqualifying events');
            obj.panel_rejection.number = num_buttons;
            obj.panel_rejection.margin = [0.05 0.02];
            obj.panel_rejection.make;
             
            %%%%%%%%%%% panel statistics %%%%%%%%%%%%%%%
            
            num_buttons = 3;
            pos = pos-num_buttons;
            obj.panel_statistics = GraphicPanel(obj.owner, [0 pos/N_left 0.2 num_buttons/N_left], 'statistics', 1); % last input is for vertical (default)
            obj.panel_statistics.addButton('button_num_batches', 'total_batches', 'info', 'Nbtch= ', '', 'small', 1/3, [], [], 'Total number of batches processed');
            obj.panel_statistics.addButton('button_num_events', 'num_events', 'info', 'Nev= ', '', 'small', 1/3, [], [], 'Total number of found events');
            obj.panel_statistics.addButton('button_num_kept', 'num_kept', 'info', 'Nkept= ', '', 'small', 1/3, [], [], 'Number of kept events');
            obj.panel_statistics.addButton('button_star_hours', 'star_hours_total', 'info', 'star hours= ', '', 'small', 2/3, [], [], 'How may star hours were scanned, total');
            obj.panel_statistics.addButton('button_star_hours_lost', 'star_hours_lost', 'info', 'lost= ', '', 'small', 1/3, [], [], 'How may star hours were lost due to overlap with events');
            obj.panel_statistics.addButton('button_make_histogram', 'histogram', 'push', 'show hist', '', '', 1, [], [], 'Plot a histogram of all S/N results');
            obj.panel_statistics.number = num_buttons;
            
            obj.panel_statistics.make;
            
            %%%%%%%%%%% panel memory %%%%%%%%%%%%%%%
            
            num_buttons = 1;
            pos = pos-num_buttons;
            obj.panel_memory = GraphicPanel(obj.owner, [0 pos/N_left 0.2 num_buttons/N_left], 'memory control', 1); % last input is for vertical (default)
            obj.panel_memory.addButton('button_clear', 'clearOneEventMemory', 'push', 'clear mem ', '', '', 0.5, [], [], 'Remove images and raw fluxes for current event');
            obj.panel_memory.addButton('button_load', 'loadOneEventMemory', 'push', 'load mem ', '', '', 0.5, [], [], 'Load from file images and raw fluxes for current event');
            obj.panel_memory.number = num_buttons;
            
            obj.panel_memory.make;
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            num_buttons = 1; pos = pos - num_buttons;            
            obj.panel_close = uipanel('Position', [0 pos 0.2 num_buttons/N_left]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            %%%%%%%%%%%%%%%%%%%%%% MIDDLE %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_middle = 15; pos = N_middle;
            
            %%%%%%%%%%% panel obs pars %%%%%%%%%%%%%%%%%%
            
            num_buttons = 1; pos = pos - 1;
            obj.panel_obs_pars = GraphicPanel(obj.owner, [0.2 pos/N_middle 0.8 num_buttons/N_middle], 'observational parameters', 0); % last input is for horizontal
            obj.panel_obs_pars.addButton('button_pars', 'obs_pars_str', 'info', ' ', '', 'small');
            obj.panel_obs_pars.number = num_buttons;
            obj.panel_obs_pars.make;
            
            %%%%%%%%%%% panel phot pars %%%%%%%%%%%%%%%%%%
            
            num_buttons = 1; pos = pos - 1;
            obj.panel_phot_pars = GraphicPanel(obj.owner, [0.2 pos/N_middle 0.8 num_buttons/N_middle], 'photmetric parameters', 0); % last input is for horizontal
            obj.panel_phot_pars.addButton('button_pars', 'phot_pars_str', 'info', ' ', '', 'small');
            obj.panel_phot_pars.number = num_buttons;
            obj.panel_phot_pars.make;
            
            %%%%%%%%%%% panel event pars %%%%%%%%%%%%%%%%%%
            
            num_buttons = 1; pos = pos - 1;
            obj.panel_phot_pars = GraphicPanel(obj.owner, [0.2 pos/N_middle 0.8 num_buttons/N_middle], 'event parameters', 0); % last input is for horizontal
            obj.panel_phot_pars.addButton('button_pars', 'event_pars_str', 'info', ' ', '', 'small');
            obj.panel_phot_pars.number = num_buttons;
            obj.panel_phot_pars.make;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            num_buttons = 11; pos = pos - num_buttons;
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 pos/N_middle 0.8 num_buttons/N_middle]);
            
            %%%%%%%%%%% panel chooser %%%%%%%%%%%%%%%%%%
            
            obj.panel_chooser = GraphicPanel(obj.owner, [0.2 0/N_middle 0.8 1/N_middle], '', 0); % last input is for horizontal
            obj.panel_chooser.addButton('button_play', '', 'custom', 'PLAY', '', '', [], [], [], 'Show all events in a loop');
            obj.panel_chooser.addButton('button_prev', 'display_prev_event', 'push', 'PREV', '', '', [], [], [], 'Show previous event');
            obj.panel_chooser.addButton('input_index', 'display_event_idx', 'input', 'idx= ', '', '', [], [], [], 'Index of currently displayed event');
            obj.panel_chooser.addButton('button_next', 'display_next_event', 'push', 'NEXT', '', '', [], [], [], 'Show next event');
            obj.panel_chooser.addButton('button_kept', 'use_display_kept_events', 'toggle', 'ALL', 'KEPT');
            obj.panel_chooser.margin = [0.02 0.1];
            obj.panel_chooser.number = 5;
            obj.panel_chooser.make;
            
            obj.panel_chooser.button_play.Callback = @obj.callback_play;
            
            obj.update;
            
            obj.owner.show;
            
        end
        
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
                        
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
        function callback_play(obj, ~, ~)
            
            if obj.debug_bit, disp('callback: play/stop'); end
            
            if obj.owner.use_display_kept_events
                if obj.owner.num_kept<2
                    obj.owner.show;
                    drawnow;
                    obj.update;
                    return;
                end
            else
                if obj.owner.num_events<2
                    obj.owner.show;
                    drawnow;
                    obj.update;
                    return;
                end
            end
            
            if obj.brake_bit
                for ii = 1:10000
                    
                    if obj.brake_bit, break; end
                    if ~obj.check, break; end
                    obj.owner.display_next_event;
                    obj.owner.show; 
                    drawnow;
                    
                end
            else
                obj.brake_bit = 1;
            end
            
        end
        
    end
    
end