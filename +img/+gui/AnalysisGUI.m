classdef AnalysisGUI < handle
    
    properties 
        
        owner@img.Analysis; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 16;
        edit_font_size = 14;
        small_font_size = 12;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_controls;
        
        panel_contrast;
    
        panel_close;
        button_close;
    
        panel_objects;
            
        panel_progress
        
        panel_image;
        button_reset_axes;
        axes_image;
    
        panel_run;
        
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = AnalysisGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'img.Analysis')
                
                if obj.debug_bit, fprintf('AnalysisGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input an img.Analysis to constructor of AnalysisGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            
            obj.fig = util.plot.FigHandler('analysis');
            obj.fig.clear;
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 25;
            movegui(obj.fig.fig, 'center');
            
            %%%%%%%%%%%%%%%%%%% LEFT SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_left = 10; pos = N_left;
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            % Add buttons using obj.addButton(button_name, var_name='', type='', str1='', str2='', font_size='', split=1, color_on=[], color_off=[])
            
            N = 4; pos = pos - N;
            
            obj.panel_controls = GraphicPanel(obj.owner, [0 pos/N_left 0.2 N/N_left], 'controls', 1); % last input is for vertical (default)
            obj.panel_controls.number = N;
            obj.panel_controls.addButton('button_bg_stack', 'use_background_stack', 'toggle', 'b/g stack', 'b/g stack', 'small', 0.5, 'red');
            obj.panel_controls.addButton('button_bg_cutouts', 'use_background_cutouts', 'toggle', 'b/g cutouts', 'b/g cutouts', 'small', 0.5, 'red');
            obj.panel_controls.make;
            
            %%%%%%%%%%% panel contrast %%%%%%%%%%%%%%%
            
            N = 5; pos = pos - N;
            
            obj.panel_contrast = util.plot.ContrastLimits(obj.axes_image, obj.fig.fig, [0 pos/N_left 0.2 5/N_left], 1); % last input is for vertical (default)
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            N = 1; pos = pos - N;
            
            obj.panel_close = uipanel('Position', [0 pos 0.2 N/N_left]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            %%%%%%%%%%%%%%%%%%% RIGHT SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_right = 10; pos = N_right;
            
            %%%%%%%%%%% panel objects %%%%%%%%%%%%%%%%
            
            N = 10; pos = pos - N;
            
            obj.panel_objects = GraphicPanel(obj.owner, [0.8 pos/N_right 0.2 N/N_right], 'objects'); 
            obj.panel_objects.number = N;
            obj.panel_objects.addButton('button_parameters', 'pars', 'push', 'Parameters');
            obj.panel_objects.addButton('button_reader', 'reader', 'push', 'Reader');
            obj.panel_objects.addButton('button_calibration', 'cal', 'push', 'Calibration');
            obj.panel_objects.addButton('button_clipper', 'clip', 'push', 'Clipper');
            obj.panel_objects.addButton('button_background', 'back', 'push', 'Background');
            obj.panel_objects.addButton('button_photometry', 'phot', 'push', 'Photometry');
            obj.panel_objects.addButton('button_light_basic', 'light_basic', 'push', 'basic', '', '', 1/3);
            obj.panel_objects.addButton('button_light_ap', 'light_ap', 'push', 'aper', '', '', 1/3);
            obj.panel_objects.addButton('button_light_gauss', 'light_gauss', 'push', 'gauss', '', '', 1/3);
            obj.panel_objects.addButton('button_finder', 'finder', 'push', 'Finder');
            
            obj.panel_objects.make;
            
            %%%%%%%%%%%%%%%%%%%%%% MIDDLE %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_middle = 10; pos = N_middle;
            
            %%%%%%%%%%% panel progres %%%%%%%%%%%%%%%%
            
            N = 1; pos = pos - N;
            obj.panel_progress = GraphicPanel(obj.owner, [0.2 pos/N_middle 0.6 N/N_middle]);
            obj.panel_progress.addButton('button_progress', '', 'custom', ' ', '', 'small');
            obj.panel_progress.make;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            N = 8; pos = pos - N;
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 pos/N_middle 0.6 N/N_middle]);
                        
            obj.makeAxes;
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.9 0.95 0.1 0.05], obj.owner, '', 'custom','reset');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            
            
            %%%%%%%%%%% panel run/stop %%%%%%%%%%%%%%%
            
            N = 1; pos = pos - N;
            
            obj.panel_run = GraphicPanel(obj.owner, [0.2 pos/N_middle 0.6 N/N_middle]);
            
            obj.panel_run.addButton('button_run', '', 'custom', 'RUN');
            obj.panel_run.make;
            obj.panel_run.button_run.Callback = @obj.callback_run;
            
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
            
            obj.panel_progress.button_progress.String = obj.owner.prog.show;
            
            if obj.owner.brake_bit
                obj.panel_run.button_run.String = 'RUN';
                
            else
                obj.panel_run.button_run.String = 'STOP';
            end
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
    
        function callback_run(obj, ~, ~)
            
            if obj.owner.brake_bit
                if obj.debug_bit, disp('callback: run'); end            
                obj.owner.run;
            else
                if obj.debug_bit, disp('callback: stop'); end            
                obj.owner.brake_bit = 1;
            end
            
            obj.update;
            
        end
             
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end