classdef LightGUI < handle
    
    properties 
        
        owner@img.Lightcurves; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 16;
        edit_font_size = 14;
        small_font_size = 12;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_methods;
        
        panel_control;
        
        panel_display;
        
        panel_process;
        
        panel_info;
        
        panel_close;
        button_close;
        
        panel_image;
        button_reset_axes;
        button_cut_number;
        axes_image;
    
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = LightGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'img.Lightcurves')
                
                if obj.debug_bit, fprintf('LightGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input an img.Lightcurves to constructor of LightGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            
            obj.fig = util.plot.FigHandler('lightcurves');
            obj.fig.clear;
            
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 25;
            movegui(obj.fig.fig, 'center');
            
            %%%%%%%%%%% panel methods %%%%%%%%%%%%%%%
            
            obj.panel_methods = GraphicPanel(obj.owner, [0 0.9 0.2 0.1], 'methods'); 
            obj.panel_methods.addButton('button_signal_method', 'signal_method', 'info', ' ', '', '', 0.5);
            obj.panel_methods.addButton('button_background_method', 'background_method', 'info', ' ', '', '', 0.5);
            obj.panel_methods.make;
            
            %%%%%%%%%%% panel control %%%%%%%%%%%%%%%
            
            obj.panel_control = GraphicPanel(obj.owner, [0 0.6 0.2 0.3], 'control'); 
            obj.panel_control.number = 4;
            obj.panel_control.addButton('button_double_up', 'use_double_up', 'toggle', 'double up', 'double up', '', [], 'red');
            obj.panel_control.make;
            
            %%%%%%%%%%% panel display %%%%%%%%%%%%%%%%
            
            obj.panel_display = GraphicPanel(obj.owner, [0 0.3 0.2 0.3], 'display');
            obj.panel_display.number = 4;
            obj.panel_display.addButton('button_show_what', 'show_what', 'picker', 'flux'); 
            obj.panel_display.addButton('button_flux_cal', 'show_flux_cal', 'picker', 'raw'); 
            obj.panel_display.addButton('input_num_stars', 'show_num_stars', 'input', 'Nstars= ');
            obj.panel_display.addButton('button_smooth', 'use_smooth', 'toggle', 'smooth', 'smooth', '', 0.5, 'red');
            obj.panel_display.addButton('input_smooth', 'smooth_interval', 'input', ' ', '', '', 0.5);
            
            obj.panel_display.make;
            
            obj.panel_display.button_show_what.String = obj.owner.show_what_list;
            obj.panel_display.button_show_what.Callback = @obj.callback_show_what;
            
            obj.panel_display.button_flux_cal.String = obj.owner.show_cal_list;
            obj.panel_display.button_flux_cal.Callback = @obj.callback_flux_cal;
            
            %%%%%%%%%%% panel process %%%%%%%%%%%%%%%%
            
            obj.panel_process = GraphicPanel(obj.owner, [0 0.1 0.2 0.2], 'process');
            obj.panel_process.number = 2;
            obj.panel_process.addButton('button_background', 'use_subtract_backgrounds', 'toggle', 'sub b/g', 'sub b/g', '', [], 'red');
            
            obj.panel_process.make;
            
            %%%%%%%%%%% panel info %%%%%%%%%%%%%%%%%%%
            
            obj.panel_info = GraphicPanel(obj.owner, [0.2 0.9 0.8 0.1], 'info', 0);
            obj.panel_info.number = 4;
            
            obj.panel_info.make;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 0.0 0.8 0.9]);
                        
            obj.makeAxes;
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.9 0.95 0.1 0.05], obj.owner, '', 'custom', 'reset');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            
            obj.button_cut_number = GraphicButton(obj.panel_image, [0.0 0.95 0.1 0.05], obj.owner, '', 'custom', '');
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
                        
            obj.panel_close = uipanel('Position', [0 0 0.2 0.1]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            obj.update;
            
        end
            
        function makeAxes(obj, ~, ~)
            
            delete(obj.axes_image);
            
            obj.axes_image = axes('Parent', obj.panel_image);
            
            obj.update;
            
        end
                
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            obj.owner.show;
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_show_what(obj, hndl, ~)
            
            if obj.debug_bit, disp('callback: show_what'); end
            
            obj.owner.show_what = hndl.String{hndl.Value};

            obj.update;
            
        end
        
        function callback_flux_cal(obj, hndl, ~)
            
            if obj.debug_bit, disp('callback: flux_cal'); end
            
            obj.owner.show_flux_cal = hndl.String{hndl.Value};

            obj.update;
            
        end
        
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end