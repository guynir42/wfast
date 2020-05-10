classdef EphemGUI < handle
    
    properties 
        
        owner@head.Ephemeris; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        menus = {};
        
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
        
        panel_inputs;
        
        panel_outputs;
            
        panel_close;
        button_close;
        
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = EphemGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'head.Ephemeris')
                
                if obj.debug_bit>1, fprintf('EphemGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                                
            else
                error('Input a head.Ephemeris object to constructor of EphemGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.MenuItem;
            
            obj.buttons = {};
            obj.menus = {};
            
            obj.fig = util.plot.FigHandler('ephemeris');
            obj.fig.clear;
            obj.fig.bottom = 5;
            obj.fig.height = 12;
            obj.fig.width = 8;
            movegui(obj.fig.fig, 'center');
            obj.fig.left = obj.fig.left +5;
            
            %%%%%%%%%%%%%%%%%%%%%%% MENUS %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % MenuItem(parent, text, type, variable, tooltip, separator)
            % obj.addButton(name, text, type, variable, tooltip, separator)
            % menu types: menu, toggle, push, input, input_text, info, custom
            
            %%%%%%%%%%%%%%%%%%% BUTTONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N = 10; 
            pos = N;
            
            % Add buttons using obj.addButton(button_name, var_name='', type='', str1='', str2='', font_size='', split=1, color_on=[], color_off=[], tooltip)
            
            %%%%%%%%%%%%%%%%%%%% panel inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            num_buttons = 5;
            pos = pos-num_buttons;
            obj.panel_inputs = GraphicPanel(obj.owner, [0 pos/N 1 num_buttons/N], 'inputs', 1); % last input is for vertical (default)
            obj.panel_inputs.addButton('button_constraints', 'constraints', 'push', 'constraints', '', '', 0.5, '', '', 'choose constraints for automatically choosing fields'); 
            obj.panel_inputs.addButton('button_field_picker', '', 'custom', 'pick field', '', '', 0.5, '', '', 'choose a field from a list of default field types'); 
            obj.panel_inputs.addButton('input_name', 'name', 'input_text', 'name= ', '', '', 0.8, '', '', 'enter a name for the object'); 
            obj.panel_inputs.addButton('button_resolve', 'resolve', 'push', 'resolve', '', '', 0.2, '', '', 'resolve the current name using default fields or using SIMBAD'); 
            obj.panel_inputs.addButton('input_time', 'time_str', 'input_text', 'time= ', '', '', 0.8, '', '', 'input the time for this object (full date or time of day in hours)'); 
            obj.panel_inputs.addButton('button_update', 'update', 'push', 'now', '', '', 0.2, '', '', 'update object time to current time'); 
            obj.panel_inputs.addButton('input_ra', 'RA', 'input_text', 'RA= ', '', '', 0.5, '', '', 'input the right ascention of the object (hours)'); 
            obj.panel_inputs.addButton('input_dec', 'Dec', 'input_text', 'Dec= ', '', '', 0.5, '', '', 'input the declination of the object'); 
            obj.panel_inputs.number = num_buttons;
            
            obj.panel_inputs.make;
            
            obj.panel_inputs.button_field_picker.Callback = @obj.callback_field_picker;
            
            %%%%%%%%%%%%%%%%%%%% panel outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            num_buttons = 4;
            pos = pos-num_buttons;
            obj.panel_outputs = GraphicPanel(obj.owner, [0 pos/N 1 num_buttons/N], 'outputs', 1); % last input is for vertical (default)
            
            obj.panel_outputs.addButton('button_lst', 'LST', 'info', 'LST= ', '', '', 0.5, '', '', 'local sidereal time'); 
            obj.panel_outputs.addButton('button_ha', 'HA', 'info', 'HA= ', '', '', 0.5, '', '', 'hour angle'); 
            
            obj.panel_outputs.addButton('button_alt', 'Alt_deg', 'info', 'ALT= ', '', '', 0.5, '', '', 'altitude above horizon (degrees)'); 
            obj.panel_outputs.addButton('button_az', 'Az_deg', 'info', 'AZ= ', '', '', 0.5, '', '', 'azimuth angle (degrees)'); 
            
            obj.panel_outputs.addButton('button_ecl', 'ECL_LAT', 'info', 'ecl= ', '', '', 0.5, '', '', 'ecliptic latitude'); 
            obj.panel_outputs.addButton('button_gal', 'GAL_LAT', 'info', 'gal= ', '', '', 0.5, '', '', 'galactic latitude'); 
            
            obj.panel_outputs.addButton('button_moon_dist', 'moon.Dist', 'info', 'moon dist.= ', '', '', 0.5, '', '', 'moon angular distance (degrees)');
            obj.panel_outputs.addButton('button_moon_ill', 'moon.IllF', 'info', 'moon ill.= ', '', '', 0.5, '', '', 'moon illumination fraction'); 
            
            obj.panel_outputs.addButton('button_airmass', 'AIRMASS', 'info', 'A.M.= ', '', '', 0.5, '', '', 'airmass'); 
            
            obj.panel_outputs.make;
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Position', [0 0 1 1/N]);
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
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_field_picker(obj, ~, ~)
            
            if obj.debug_bit>1, disp('callback: field picker'); end
            
            rep = questdlg('Choose a keyword field:', 'keyword field', 'ecliptic', 'galactic', 'moon', 'ecliptic'); 
            
            if isempty(rep)
                return;
            else
                obj.owner.input(rep);
            end
            
            obj.update;
            
        end
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end