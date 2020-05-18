classdef AndorGUI < handle
    
    properties 
        
        owner@obs.cam.Andor; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        big_font_size = 16;
        font_size = 13;
        edit_font_size = 12;
        small_font_size = 10;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_controls;
        
        panel_contrast;
    
        panel_objects;
        
        panel_aux;
        
        panel_info;
        
        panel_start;
        
        panel_close;
        button_close;
        
        panel_image;
        button_reset_axes;
        button_flip;
        button_number;
        button_clipping;
        axes_image;
    
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = AndorGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'obs.cam.Andor')
                
                if obj.debug_bit>1, fprintf('AndorGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                                
            else
                error('Input an obs.cam.Andor to constructor of AndorGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            
            obj.fig = util.plot.FigHandler('Andor camera');
            obj.fig.clear;
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 25;
            movegui(obj.fig.fig, 'center');
            
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            obj.panel_controls = GraphicPanel(obj.owner, [0 7/15 0.2 8/15], 'controls');
            obj.panel_controls.number = 8;
            obj.panel_controls.addButton('button_record', 'record', 'push', 'RECORD', '', 'big', 0.6);
            obj.panel_controls.addButton('button_mode', 'cycleModes', 'push', 'science', '', 'small', 0.4);
            obj.panel_controls.addButton('button_preview', 'preview', 'push', 'PREVIEW', '', '', 0.8);
            obj.panel_controls.addButton('input_time_deep', 'expT_deep', 'input', ' ', 's','', 0.2);
            obj.panel_controls.addButton('input_num_batches', 'num_batches', 'input', 'batches= ', '','small',0.5);
            obj.panel_controls.addButton('input_batch_size', 'batch_size', 'input', 'frames= ', '', 'small', 0.5);
            obj.panel_controls.addButton('input_exp_time', 'expT', 'input', 'T= ', 's', '', 0.5);
            obj.panel_controls.addButton('input_frame_rate', 'frame_rate', 'input', 'f= ', 'Hz', '', 0.5);
            obj.panel_controls.addButton('button_frame_rate_measured', 'frame_rate_measured', 'info', 'f= ', 'Hz', '', 1);
            obj.panel_controls.addButton('button_zoom', 'use_roi', 'toggle', 'ZOOM', 'UNZOOM', 'small', 0.3);
            obj.panel_controls.addButton('input_im_size', 'im_size', 'input', ' ', '', 'edit', 0.7);
            obj.panel_controls.addButton('button_pick_center', '', 'custom', 'center', '', 'small', 0.3);
            obj.panel_controls.addButton('input_center_region', 'center_region', 'input', 'center= ', '', 'edit', 0.7);
            obj.panel_controls.addButton('button_autofocus', 'autofocus', 'push', 'Auto focus'); 
            obj.panel_controls.make;
            obj.panel_controls.button_pick_center.Callback = @obj.callback_pick_center;
            
            obj.panel_controls.button_record.Tooltip = 'Start recording full frame images. Takes <number of batches>X<number frames> total';
            obj.panel_controls.button_mode.Tooltip = 'Clcik to cycle through record modes: science, dark, flat';
            obj.panel_controls.button_preview.Tooltip = 'Take a single image with a longer exposure time. Does not save the image';
            obj.panel_controls.input_time_deep.Tooltip = 'Control the time of the PREVIEW exposure.';
            obj.panel_controls.input_num_batches.Tooltip = ['Number of files to save when using RECORD. Each file contains ' num2str(obj.owner.batch_size) ' images'];
            obj.panel_controls.input_batch_size.Tooltip = 'Number of individual frames in each file when using RECORD.';
            obj.panel_controls.input_exp_time.Tooltip = 'Exposure time for RECORD and for LIVE.';
            obj.panel_controls.input_frame_rate.Tooltip = 'Nominal frame rate the camera tries to maintain when in RECORD mode only. Set NaN to take images as fast as possible';
            obj.panel_controls.button_frame_rate_measured.Tooltip = 'Measured frame rate (calculated from the last few batches)';
            obj.panel_controls.button_zoom.Tooltip = ['Toggle zoom (ROI) mode on and off. Currently the camera is in ' obj.owner.getZoomStr ' mode'];
            obj.panel_controls.input_im_size.Tooltip = 'Set the height and width of the image (enables zoom mode)';
            obj.panel_controls.button_pick_center.Tooltip = 'Pick the center point of the ROI region by clicking the full frame image';
            obj.panel_controls.input_center_region.Tooltip = 'Manually pick the center of the ROI region';
            obj.panel_controls.button_autofocus.Tooltip = 'Start a focusing run. ';
            
            
            %%%%%%%%%%% panel contrast %%%%%%%%%%%%%%%
            
            obj.panel_contrast = ContrastLimits([], obj.fig.fig, [0 2/15 0.2 5/15]);
            
            %%%%%%%%%%% panel objects %%%%%%%%%%%%%%%%
            
            obj.panel_objects = GraphicPanel(obj.owner, [0.8 7/15 0.2 8/15], 'object'); 
            obj.panel_objects.number = 5;
            obj.panel_objects.addButton('button_buffers', 'buffers', 'push', 'Buffers');
            obj.panel_objects.addButton('button_header', 'head', 'push', 'Header');
            obj.panel_objects.addButton('button_focuser', 'focuser', 'push', 'Focuser');
            
            obj.panel_objects.make;
            
            %%%%%%%%%%% panel aux %%%%%%%%%%%%%%%%%%%%
            
            obj.panel_aux = GraphicPanel(obj.owner, [0.8 0/15 0.2 7/15], 'auxiliary');
            obj.panel_aux.number = 5;
            
            obj.panel_aux.make;
            
            %%%%%%%%%%% panel info %%%%%%%%%%%%%%%%%%%
            
            obj.panel_info = GraphicPanel(obj.owner, [0.2 14/15 0.6 1/15], 'hardware information');
            obj.panel_info.addButton('button_info', '', 'custom');
            obj.panel_info.make;
            obj.panel_info.button_info.font_size = 'small';
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 2/15 0.6 12/15]);
                        
            obj.makeAxes;
            obj.button_flip = GraphicButton(obj.panel_image, [0.0 0.95 0.15 0.05], obj.owner, 'use_show_flipped', 'toggle', 'no flip', 'flip');
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.85 0.95 0.15 0.05], obj.owner, '', 'custom','new axes');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            obj.button_number = GraphicButton(obj.panel_image, [0.0 0.00 0.15 0.1], obj.owner, 'batch_counter', 'info', 'N= ');
            obj.button_clipping = GraphicButton(obj.panel_image, [0.85 0.00 0.15 0.1], obj.owner, '', 'custom', 'ok');
           
            obj.button_number.Tooltip = 'How many batches already taken in this run';
            obj.button_flip.Tooltip = 'Flip the image 180 degrees (for viewing beyond the meridian)';
            obj.button_reset_axes.Tooltip = 'Generate a new image axes with contrast and zoom initialized';            
            obj.button_clipping.Tooltip = 'Display warning if image is saturated (currently disabled)';
            
            %%%%%%%%%%% panel start %%%%%%%%%%%%%%%%%%
            
            obj.panel_start = GraphicPanel(obj.owner, [0.2 0/15 0.6 2/15]);
            obj.panel_start.addButton('button_start', '', 'custom', 'START');
            obj.panel_start.make;
            obj.panel_start.button_start.Callback = @obj.callback_start_stop;
            
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
                        
            obj.panel_close = uipanel('Position', [0 0/15 0.2 2/15]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            obj.update;
            
        end
            
        function makeAxes(obj, ~, ~)
            
            delete(obj.axes_image);
            
            obj.axes_image = axes('Parent', obj.panel_image);
            
            obj.panel_contrast.ax = obj.axes_image;
            colorbar(obj.axes_image);
            axis(obj.axes_image, 'image');
            xlabel(obj.axes_image, 'South');
            ylabel(obj.axes_image, 'East');
            
            obj.panel_contrast.ax = obj.axes_image;
            
        end
                
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            obj.panel_contrast.font_size = obj.font_size;
            obj.panel_contrast.edit_font_size = obj.edit_font_size;
            obj.panel_contrast.update;
            
            obj.panel_controls.input_num_batches.Tooltip = ['Number of files to save when using RECORD. Each file contains ' num2str(obj.owner.batch_size) ' images'];
            obj.panel_controls.button_zoom.Tooltip = ['Toggle zoom (ROI) mode on and off. Currently the camera is in ' obj.owner.getZoomStr ' mode'];
            
            if obj.owner.brake_bit
                
                obj.panel_info.button_info.String = obj.owner.printout_HW;
            
                obj.panel_start.button_start.String = 'START VIDEO PREVIEW';
                
            else
                obj.panel_start.button_start.String = 'STOP';
            end
            
            obj.panel_controls.button_mode.String = obj.owner.mode;
            
            % choose a color for control button
            if util.text.cs(obj.owner.mode, 'science')
                obj.panel_controls.button_mode.BackgroundColor = util.plot.GraphicButton.defaultColor;
            elseif util.text.cs(obj.owner.mode, 'dark')
                obj.panel_controls.button_mode.BackgroundColor = [0.3 0.7 0.95]; % custom light blue
            elseif util.text.cs(obj.owner.mode, 'flat')
                obj.panel_controls.button_mode.BackgroundColor = 'green';
            else
                
            end
            
            % set image panel background color when recording...
            if obj.owner.brake_bit==0 && obj.owner.use_save % when recording! 
                if util.text.cs(obj.owner.mode, 'science')
                    chosen_color = 'red';
                elseif util.text.cs(obj.owner.mode, 'dark')
                    chosen_color = [0.3 0.7 0.95]; % custom light blue
                elseif util.text.cs(obj.owner.mode, 'flat')
                    chosen_color = 'green';
                else
                    chosen_color = util.plot.GraphicButton.defaultColor;
                end
                obj.panel_image.BackgroundColor = chosen_color;
            else
                obj.panel_image.BackgroundColor = util.plot.GraphicButton.defaultColor;
            end
            
            % grey out controls when taking images... 
            for ii = 1:length(obj.panel_controls.buttons)
                if obj.owner.brake_bit
                    obj.panel_controls.buttons(ii).Enable = 'on';
                else
                    obj.panel_controls.buttons(ii).Enable = 'off';
                end
            end
            
            obj.panel_controls.button_frame_rate_measured.Enable = 'on';
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_pick_center(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: pick_center'); end
            
            h = findobj(obj.axes_image, 'type','image');
            
            if isempty(h) || size(h.CData,1)~=obj.owner.max_height || size(h.CData,2)~=obj.owner.max_width
                return; % don't show a rectangle unless it is on the full frame image... 
            end
            
            obj.owner.showROI;
            
            [x,y] = getpts(obj.axes_image);
            
            if ~isempty(x) && ~isempty(y)
                
                obj.owner.center_region = [y(end) x(end)];
                
                obj.owner.showROI;
                
                obj.owner.use_roi = 1;
            
            end
            
            obj.update;
            
        end
        
        function callback_start_stop(obj, ~, ~)
            
            if obj.debug_bit>1, disp('callback: start/stop'); end
            
            if obj.owner.brake_bit
                obj.owner.live;
            else
                obj.owner.stop;
            end
            
            obj.update;
            
        end
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit>1, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end