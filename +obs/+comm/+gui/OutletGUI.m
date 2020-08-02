classdef OutletGUI < handle
    
    properties 
        
        owner@obs.comm.OutletControl; % link back to containg object

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
        
        panel_controls;
        
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = OutletGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'obs.comm.OutletControl')
                
                if obj.debug_bit>1, fprintf('OutletGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input an OutletControl object to constructor of OutletGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            
            obj.buttons = {};
            obj.panels = {}; 
            
            obj.fig = util.plot.FigHandler('Outlets');
            obj.fig.clear;
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 10;
            obj.fig.center;
            
            obj.panel_controls = GraphicPanel(obj.owner, [0 0 1 1], '', 1); % last input is for vertical (default)
            
            obj.panel_controls.addButton('button_cam_pc', 'cam_pc', 'toggle', 'cam-PC is off', 'cam-PC is on', '', 0.7, '', '', 'control the power to the camera control computer'); 
            obj.panel_controls.addButton('button_cam_pc_cycle', 'cycle_cam_pc', 'push', 'cycle cam-PC', '', '', 0.3, '', '', 'power cycle (15s) the camera control computer'); 
            
            obj.panel_controls.addButton('button_dome_pc', 'dome_pc', 'toggle', 'dome-PC is off', 'dome-PC is on', '', 0.7, '', '', 'control the power to the observatory control computer'); 
            obj.panel_controls.addButton('button_dome_pc_cycle', 'cycle_dome_pc', 'push', 'cycle dome-PC', '', '', 0.3, '', '', 'power cycle (15s) the observatory control computer'); 
            
            obj.panel_controls.addButton('button_router', 'router', 'toggle', 'router is off', 'router is on', '', 0.7, '', '', 'control the power to the router'); 
            obj.panel_controls.addButton('button_router_cycle', 'cycle_router', 'push', 'cycle router', '', '', 0.3, '', '', 'power cycle (15s) the router'); 
            
            obj.panel_controls.addButton('button_lights', 'lights', 'toggle', 'lights are off', 'lights are on', '', 0.7, '', '', 'control the power to the LED lights'); 
            obj.panel_controls.addButton('button_lights_cycle', 'cycle_lights', 'push', 'cycle lights', '', '', 0.3, '', '', 'power cycle (15s) the LED lights'); 
            
            obj.panel_controls.addButton('button_balor', 'balor', 'toggle', 'Balor camera is off', 'Balor camera is on', '', 0.7, '', '', 'control the power to the Balor camera'); 
            obj.panel_controls.addButton('button_balor_cycle', 'cycle_balor', 'push', 'cycle Balor', '', '', 0.3, '', '', 'power cycle (15s) the Balor camera'); 
            
            obj.panel_controls.addButton('button_arduinos', 'arduinos', 'toggle', 'Arduinos are off', 'Arduinos are on', '', 0.7, '', '', 'control the power to the Arduinos'); 
            obj.panel_controls.addButton('button_arduinos_cycle', 'cycle_arduinos', 'push', 'cycle Arduinos', '', '', 0.3, '', '', 'power cycle (15s) the Arduinos'); 
            
            obj.panel_controls.addButton('button_dome', 'dome', 'toggle', 'Dome is off', 'Dome is on', '', 0.7, '', '', 'control the power to the dome'); 
            obj.panel_controls.addButton('button_dome_cycle', 'cycle_dome', 'push', 'cycle dome', '', '', 0.3, '', '', 'power cycle (15s) the dome'); 
            
            obj.panel_controls.addButton('button_mount', 'mount', 'toggle', 'Mount is off', 'Mount is on', '', 0.7, '', '', 'control the power to the mount'); 
            obj.panel_controls.addButton('button_mount_cycle', 'cycle_mount', 'push', 'cycle mount', '', '', 0.3, '', '', 'power cycle (15s) the mount'); 
            
            obj.panel_controls.addButton('button_close', '', 'custom', 'CLOSE GUI'); 
            
            obj.panel_controls.make;
            
            obj.panel_controls.button_close.Callback = @obj.callback_close;
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            obj.update;
            
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
           
            c = ~isempty(obj) && ~isempty(obj.panel_controls) && isvalid(obj.panel_controls.panel);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end