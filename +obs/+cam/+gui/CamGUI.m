classdef CamGUI < handle
    
    properties 
        
        owner@obs.cam.CameraControl; % link back to containg object

        fig@util.plot.FigHandler;
        
        buttons = {};
        
        font_size = 16;
        edit_font_size = 14;
        small_font_size = 12;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        grey_list = {}; % which buttons should be disabled (greyed-out) when camera is running...
        
        panel_controls;
        panel_contrast;
        panel_inputs;
        
        panel_loupe_image;
        axes_loupe;
        panel_loupe_control;
        panel_focus;
        panel_objects;
        
        panel_run;
        panel_image;
        button_flip;
        button_reset_axes;
        button_number;
        button_clipping;
        axes_image;
    
        panel_close;
        button_close;
        
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = CamGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'obs.cam.CameraControl')
                
                if obj.debug_bit, fprintf('CamGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                                
            else
                error('Input an obs.cam.CameraControl to constructor of CamGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            
            if isempty(obj.fig)
                obj.fig = util.plot.FigHandler('Camera Control');
            end
            
            obj.fig.reset;
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 25;
            movegui(obj.fig.fig, 'center');
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            obj.panel_controls = GraphicPanel(obj.owner, [0 10/15 0.2 5/15], 'controls');
            obj.panel_controls.number = 5;
            obj.panel_controls.addButton('button_record', 'record', 'push', 'RECORD', '', '', 0.6);
            obj.panel_controls.addButton('button_mode', 'cycleModes', 'push', 'stars', '', '', 0.2);            
            obj.panel_controls.addButton('input_num_frames', 'num_frames_record', 'input', ' ', '','',0.2);
            obj.panel_controls.addButton('button_preview', 'preview', 'push', 'PREVIEW', '', '', 0.8);
            obj.panel_controls.addButton('input_preview_time', 'preview_time', 'input', ' ', 's','', 0.2);
            obj.panel_controls.addButton('button_roi_zoom', 'zoom', 'push', 'ZOOM', '', '', 0.7);
            obj.panel_controls.addButton('button_roi_cycle', 'cycleROIsize', 'push','512', '', 'small', 0.3);
            obj.panel_controls.addButton('button_roi_info', 'getROItext', 'info', ' ', '', '');
            obj.panel_controls.addButton('button_temperature', 'getTemperatureHW', 'info', 'temp= ', '', 'small', 0.5);
            obj.panel_controls.addButton('button_frame_rate', 'mean_frame_rate', 'info', 'f= ', 'Hz', 'small', 0.5);
            obj.panel_controls.make;
            
            obj.panel_controls.button_record.Callback = @obj.callback_record;
%             obj.panel_controls.button_roi_info.Callback = obj.panel_controls.button_roi_zoom.Callback;
            
            obj.grey_list = {obj.panel_controls.button_record, obj.panel_controls.button_mode, obj.panel_controls.input_num_frames, ...
                obj.panel_controls.button_preview, obj.panel_controls.input_preview_time, obj.panel_controls.button_roi_zoom, obj.panel_controls.button_temperature};
            
            %%%%%%%%%%% panel contrast %%%%%%%%%%%%%%%
            
            obj.panel_contrast = ContrastLimits(obj.axes_image, [], [0 5/15 0.2 5/15]);
            
            %%%%%%%%%%% panel inputs %%%%%%%%%%%%%%%%%
            
            obj.panel_inputs = GraphicPanel(obj.owner, [0 2/15 0.2 3/15], 'inputs');
            obj.panel_inputs.number = 3;
            obj.panel_inputs.addButton('input_batch_size', 'batch_size', 'input');
            obj.panel_inputs.addButton('input_exp_time', 'expT', 'input');
            obj.panel_inputs.make;
            
            obj.grey_list = [obj.grey_list, {obj.panel_inputs.input_batch_size, obj.panel_inputs.input_exp_time}];
            
            %%%%%%%%%%% panel loupe image %%%%%%%%%%%%
            
            obj.panel_loupe_image = uipanel('Title', 'loupe magnifier', 'Position', [0.8 10/15 0.2 5/15]);
            obj.axes_loupe = axes('Parent', obj.panel_loupe_image);
            
            %%%%%%%%%%% panel loupe control %%%%%%%%%%           
            
            obj.panel_loupe_control = GraphicPanel(obj.owner.loupe, [0.8 7/15 0.2 3/15]);
            obj.panel_loupe_control.addButton('button_mode', '', 'picker', 'zoom','','',0.5);
            obj.panel_loupe_control.addButton('button_autodyn', 'use_autodyn', 'toggle', 'autodyn', 'autodyn','',0.5);
            obj.panel_loupe_control.addButton('input_size', 'rect_size', 'input', 'size= ', '','',0.5);
            obj.panel_loupe_control.addButton('button_draw_rect', 'use_show_rect', 'toggle', 'draw', 'draw', '', 0.5);
            obj.panel_loupe_control.addButton('button_pick', 'choosePoint', 'push', 'pick', '','', 0.5);
            obj.panel_loupe_control.addButton('button_autofind', 'use_autofind', 'toggle', 'autofind', 'autofind', '', 0.5);
            obj.panel_loupe_control.make;
            obj.panel_loupe_control.button_autodyn.color_on = 'red';
            obj.panel_loupe_control.button_draw_rect.color_on = 'red';
            obj.panel_loupe_control.button_autofind.color_on = 'red';
            
            obj.panel_loupe_control.button_mode.Callback = @obj.callback_loupe_mode;            
            
            %%%%%%%%%%% panel focus %%%%%%%%%%%%%%%%%%
            
            obj.makeFocusPanel;
            
            %%%%%%%%%%% panel objects %%%%%%%%%%%%%%%%
            
            obj.panel_objects = GraphicPanel(obj.owner, [0.8 0/15 0.2 5/15], 'objects');
            obj.panel_objects.addButton('button_camera', 'cam', 'push', 'Cam HW GUI');
            obj.panel_objects.addButton('button_buffers', 'buffers', 'push', 'Buffers');
            obj.panel_objects.addButton('button_pars', 'pars', 'push', 'Parameters');
            obj.panel_objects.number = 5;
            obj.panel_objects.make;
            
            %%%%%%%%%%% panel run %%%%%%%%%%%%%%%%%%%%
            
            obj.panel_run = GraphicPanel(obj.owner, [0.2 0/15 0.6 2/15]);
            obj.panel_run.addButton('button_run', '', 'custom', 'RUN');            
            obj.panel_run.make;
            obj.panel_run.button_run.Callback = @obj.callback_run;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 2/15 0.6 13/15]);
                        
            obj.makeAxes;
            
            obj.button_flip = GraphicButton(obj.panel_image, [0.05 0.95 0.15 0.05], obj.owner, 'use_flip_image', 'toggle', 'no flip', 'flip');
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.85 0.95 0.15 0.05], obj.owner, '', 'custom','reset');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            
            obj.button_number = GraphicButton(obj.panel_image, [0.0 0.00 0.15 0.1], obj.owner, 'batch_counter', 'info', 'N= ');
            obj.button_clipping = GraphicButton(obj.panel_image, [0.85 0.00 0.15 0.1], obj.owner, '', 'custom', 'ok');
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
                        
            obj.panel_close = uipanel('Position', [0 0 0.2 2/15]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            obj.update;
            
        end
            
        function makeFocusPanel(obj)

            obj.panel_focus = util.plot.GraphicPanel(obj.owner.focuser, [0.8 5/15 0.2 2/15], 'focus', 1, 'gui_cam');
            obj.panel_focus.number = 2;
            
            if ~isempty(obj.owner.focuser) % && ~isempty(obj.owner.focuser.hndl) && obj.owner.focuser.hndl.Connected
                obj.panel_focus.addButton('button_down', 'down', 'push', 'DOWN', '', 'small', 0.3);
                obj.panel_focus.addButton('input_step', 'step', 'input', 'step= ', '', '', 0.4);
                obj.panel_focus.addButton('button_up', 'up', 'push', 'UP', '', 'small', 0.3);
                obj.panel_focus.addButton('button_connect', 'connect', 'push', 'connect', '', 'small', 0.3);
                obj.panel_focus.addButton('input_pos', 'pos', 'input', 'pos= ', '', '', 0.4);
                obj.panel_focus.addButton('button_af', '', 'custom', 'AF', '', 'small', 0.3);
            
            else
                disp('Cannot connect to focuser!');
            end
            
            obj.panel_focus.make;
            
        end
        
        function makeAxes(obj, ~, ~)
            
            delete(obj.axes_image);
            
            obj.axes_image = axes('Parent', obj.panel_image);
            imagesc(obj.axes_image, zeros(512));
            obj.panel_contrast.ax = obj.axes_image;
            colorbar(obj.axes_image);
            axis(obj.axes_image, 'image');
            xlabel(obj.axes_image, 'South');
            ylabel(obj.axes_image, 'East');
            title('');
                        
            obj.panel_contrast.ax = obj.axes_image;
                        
        end
                
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            obj.panel_controls.button_roi_cycle.String = ['->' num2str(obj.owner.current_ROI_size)];
            if obj.owner.checkIsZoomed
                obj.panel_controls.button_roi_zoom.String = 'UNZOOM';
            else
                obj.panel_controls.button_roi_zoom.String = 'ZOOM';
            end
            
            blue = [0.1 0.4 0.95]; % custom light blue
            
            if util.text.cs(obj.owner.mode, 'dark')
                obj.panel_controls.button_mode.BackgroundColor = blue;
            elseif util.text.cs(obj.owner.mode, 'flat')
                obj.panel_controls.button_mode.BackgroundColor = 'green';
            else
                obj.panel_controls.button_mode.BackgroundColor = util.plot.GraphicButton.defaultColor;
            end
            
            if obj.owner.brake_bit
                obj.panel_run.button_run.String = 'RUN';
            else
                obj.panel_run.button_run.String = 'STOP';
            end
            
            obj.panel_controls.button_mode.String = obj.owner.mode;
            
            if obj.owner.is_clipping
                obj.button_clipping.control.BackgroundColor = 'yellow';
                obj.button_clipping.String = 'Clipping!';
            else
                obj.button_clipping.control.BackgroundColor = util.plot.GraphicButton.defaultColor;
                obj.button_clipping.String = 'ready';
                if obj.owner.now_live
                    obj.button_clipping.String = 'live';
                elseif obj.owner.now_recording
                    obj.button_clipping.String = 'recording';
                end
            end
            
            if obj.owner.now_recording                
                
                obj.panel_image.BackgroundColor = 'red';
                if util.text.cs(obj.owner.mode, 'dark')
                    obj.panel_image.BackgroundColor = blue;
                elseif util.text.cs(obj.owner.mode, 'flat')
                    obj.panel_image.BackgroundColor = 'green';
                end
                
            else
                obj.panel_image.BackgroundColor = util.plot.GraphicButton.defaultColor;
            end
            
            if ~obj.owner.brake_bit
                for ii = 1:length(obj.grey_list)
                    obj.grey_list{ii}.Enable = 'off';
                end
            else
                for ii = 1:length(obj.grey_list)
                    obj.grey_list{ii}.Enable = 'on';
                end
            end
            
        end
        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_loupe_mode(obj, hndl, ~)
            
            
            
        end
        
        function callback_record(obj, ~, ~)
            
            N = util.text.extract_numbers(obj.panel_controls.input_num_frames.String);
            N = N{1};
            
            if obj.debug_bit, disp(['callback: record. N= ' num2str(N)]); end
            
            obj.owner.record('num_batches', N);
            
        end
        
        function callback_run(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: run/stop'); end
            
            if obj.owner.brake_bit
                obj.owner.live;
            else
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