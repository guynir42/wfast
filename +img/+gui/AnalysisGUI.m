classdef AnalysisGUI < handle
    
    properties 
        
        owner@img.Analysis; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 12;
        big_font_size = 16;
        edit_font_size = 10;
        small_font_size = 8;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_controls;
        
        panel_contrast;
    
        panel_close;
        button_close;
    
        panel_objects;
            
        panel_progress
        
        panel_info;
        
        panel_image;
        button_reset_axes;
        button_batch_counter;
        input_num_rect;
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
            
            N_left = 20; pos = N_left;
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            % Add buttons using obj.addButton(button_name, var_name='', type='', str1='', str2='', font_size='', split=1, color_on=[], color_off=[])
            
            N = 13; pos = pos - N;
            
            obj.panel_controls = GraphicPanel(obj.owner, [0 pos/N_left 0.2 N/N_left], 'controls', 1); % last input is for vertical (default)
            obj.panel_controls.number = N;
            obj.panel_controls.addButton('button_reset', 'reset', 'push', 'RESET', '', '', 0.5, '', '', 'Start a new run by reseting all events and lightcurves');
            obj.panel_controls.addButton('input_num_batches', 'num_batches', 'input', 'Nbatch= ', '', 'small', 0.5, '', '', 'Maximum batches, limited by user input or by number of files in reader'); 
            obj.panel_controls.addButton('button_bg_stack', 'use_background_stack', 'toggle', 'b/g stack off', 'b/g stack on', 'small', 0.5, 'red', '', 'Subtract background from stack images');
            obj.panel_controls.addButton('button_bg_cutouts', 'use_background_cutouts', 'toggle', 'b/g cutouts off', 'b/g cutouts on', 'small', 0.5, 'red', '', 'Subtract background from cutouts');
            obj.panel_controls.margin = [0.03 0.01];
            obj.panel_controls.make;
            
            %%%%%%%%%%% panel contrast %%%%%%%%%%%%%%%
            
            N = 5; pos = pos - N;
            
            obj.panel_contrast = util.plot.ContrastLimits(obj.axes_image, obj.fig.fig, [0 pos/N_left 0.2 5/N_left], 1, [0.01 0.01]); % 4th input is for vertical (default), 5th is margin
            obj.panel_contrast.font_size = obj.font_size;
            obj.panel_contrast.big_font_size = obj.big_font_size;
            obj.panel_contrast.small_font_size = obj.small_font_size;
            obj.panel_contrast.edit_font_size = obj.edit_font_size;
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            N = 2; pos = pos - N;
            
            obj.panel_close = uipanel('Position', [0 pos 0.2 N/N_left]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1]+[0.1 0.1 -0.2 -0.2], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            %%%%%%%%%%%%%%%%%%% RIGHT SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_right = 20; pos = N_right;
            
            %%%%%%%%%%% panel objects %%%%%%%%%%%%%%%%
            
            N = 20; pos = pos - N;
            
            obj.panel_objects = GraphicPanel(obj.owner, [0.8 pos/N_right 0.2 N/N_right], 'objects'); 
            obj.panel_objects.number = N;
            obj.panel_objects.addButton('button_parameters', 'pars', 'push', 'Parameters');
            obj.panel_objects.addButton('button_reader', 'reader', 'push', 'Reader');
            obj.panel_objects.addButton('button_calibration', 'cal', 'push', 'Calibration');
            obj.panel_objects.addButton('button_clipper', 'clip', 'push', 'Clipper');
            obj.panel_objects.addButton('button_background', 'back', 'push', 'Background');
            obj.panel_objects.addButton('button_photometry', 'phot', 'push', 'Photometry');
%             obj.panel_objects.addButton('button_light_basic', 'light_basic', 'push', 'basic', '', '', 1/3);
            obj.panel_objects.addButton('button_light_ap', 'light_ap', 'push', 'lightcurves', '', '', 1);
%             obj.panel_objects.addButton('button_light_gauss', 'light_gauss', 'push', 'gauss', '', '', 1/3);
            obj.panel_objects.addButton('button_finder', 'finder', 'push', 'Finder');
            obj.panel_objects.margin = [0.1 0.005];
            obj.panel_objects.make;
            
            %%%%%%%%%%%%%%%%%%%%%% MIDDLE %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_middle = 10; pos = N_middle;
            
            %%%%%%%%%%% panel progres %%%%%%%%%%%%%%%%
            
            N = 1; pos = pos - N;
            obj.panel_progress = GraphicPanel(obj.owner, [0.2 pos/N_middle 0.6 N/N_middle]);
            obj.panel_progress.addButton('button_progress', '', 'custom', ' ', '', 'edit');
            obj.panel_progress.margin = [0.0 0.0];
            obj.panel_progress.make;
            
            %%%%%%%%%%% panel info %%%%%%%%%%%%%%%%%%%
            
            N = 1; pos = pos - N;
            obj.panel_info = GraphicPanel(obj.owner, [0.2 pos/N_middle 0.6 N/N_middle], 0); 
            obj.panel_info.addButton('button_fwhm', 'FWHM', 'info', 'FWHM= ', 'pix', '', 0.5); 
            obj.panel_info.addButton('button_seeing', 'seeing', 'info', 'seeing= ', '"', '', 0.5);
            obj.panel_info.margin = 0.1;
            obj.panel_info.make;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            N = 7; pos = pos - N;
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 pos/N_middle 0.6 N/N_middle]);
                        
            obj.makeAxes;
            
            obj.button_batch_counter = GraphicButton(obj.panel_image, [0.0 0.0 0.15 0.05], obj.owner, 'batch_counter', 'info', 'N= ');
            obj.button_batch_counter.Tooltip = 'How many batches already finished';
            
            obj.input_num_rect = GraphicButton(obj.panel_image, [0.0 0.95 0.15 0.05], obj.owner, 'display_num_rect_stars', 'input', ' ', ' rect');
            obj.input_num_rect.Tooltip = 'How many rectangles (maximum) to show on screen';
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.85 0.95 0.15 0.05], obj.owner, '', 'custom', 'new axes', '');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            obj.button_reset_axes.Tooltip = 'Create a new image axis, zoomed out and with default contrast limits'; 
            
            
            %%%%%%%%%%% panel run/stop %%%%%%%%%%%%%%%
            
            N = 1; pos = pos - N;
            
            obj.panel_run = GraphicPanel(obj.owner, [0.2 pos/N_middle 0.6 N/N_middle]);
            
            obj.panel_run.addButton('button_run', '', 'custom', 'RUN', '', 'big');
            obj.panel_run.margin = [0.02 0.1];
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
                if obj.owner.batch_counter==0
                    obj.panel_run.button_run.String = 'START NEW RUN';
                else
                    obj.panel_run.button_run.String = 'CONTINUE';
                end
            else
                obj.panel_run.button_run.String = 'STOP';
            end
            
            obj.panel_contrast.update;
            
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