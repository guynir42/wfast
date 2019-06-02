classdef AcqGUI < handle
    
    properties 
        
        owner@img.Acquisition; % link back to containg object

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
    
        panel_objects;
        
        panel_save;
        
        panel_run;
        button_run;
        
        panel_info;
        panel_image;
        button_reset_axes;
        button_batch_num;
        button_show_what;
        button_flip;
        axes_image;
    
        panel_close;
        button_close;
        
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = AcqGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'img.Acquisition')
                
                if obj.debug_bit, fprintf('AcqGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input an Acquisition object to constructor of AcqGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            
            obj.fig = util.plot.FigHandler('Acquisition');
            
            obj.fig.bottom = 5;
            obj.fig.height = 20;
            obj.fig.width = 36;
            movegui(obj.fig.fig, 'center');
            obj.fig.clear;
            
            N_left = 15;
            N_right = 15;
            W_left = 0.25;
            W_right = 0.15;
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%
            
            pos = N_left;
            N = 9;
            pos = pos - N;
            obj.panel_controls = GraphicPanel(obj.owner, [0 pos/N_left W_left N/N_left], 'controls');
            obj.panel_controls.number = N;
            obj.panel_controls.addButton('button_source_choose', 'chooseSource', 'push', 'choose src', '', '', 2/3);
            obj.panel_controls.addButton('button_source_gui', 'src', 'push', 'GUI', '', '', 1/3);
            
            obj.panel_controls.addButton('button_reset', 'reset', 'push', 'RESET', '', '', 1/3);
            obj.panel_controls.addButton('input_name', 'run_name', 'input_text', 'name= ', '', '', 2/3);
            
            obj.panel_controls.addButton('button_num_files', 'num_files', 'info', 'max= ', '', '', 1/3);
            obj.panel_controls.addButton('input_num_batches', 'num_batches', 'input', 'Nbatches= ', '', '', 1/3);
            obj.panel_controls.addButton('input_batch_size', 'batch_size', 'input', 'Nframes= ', '', '', 1/3);
            
            obj.panel_controls.addButton('button_stars_found', 'num_stars_found', 'info', ' ', '', '', 1/3);
            obj.panel_controls.addButton('input_num_stars', 'num_stars', 'input', 'Nstars= ', '', '', 1/3);
            obj.panel_controls.addButton('input_cut_size', 'cut_size', 'input', 'size= ', '', '', 1/3);
            
            obj.panel_controls.addButton('input_edges', 'avoid_edges', 'input', 'edge= ', '', '', 1/3);
            obj.panel_controls.addButton('input_num_bgs', 'num_backgrounds', 'input', 'Nbgs= ', '', '', 1/3);
            obj.panel_controls.addButton('input_cut_size_bg', 'cut_size_bg', 'input', 'size= ', '', '', 1/3);
            
            obj.panel_controls.addButton('input_expT', 'expT', 'input', 'T= ', 's', '', 0.5);
            obj.panel_controls.addButton('input_frame_rate', 'frame_rate', 'input', 'f= ', 'Hz', '', 0.5);
            
            obj.panel_controls.addButton('button_mextractor', 'use_mextractor', 'toggle', 'mextractor', 'mextractor', '', 0.5, 'red');
            obj.panel_controls.addButton('button_arbitrary_pos', 'use_arbitrary_pos', 'toggle', 'arbitrary', 'arbitrary', '', 0.5, 'red');
            
            obj.panel_controls.addButton('button_adjust', 'use_adjust_cutouts', 'toggle', 'adjust', 'adjust', '', 0.5, 'red');
            obj.panel_controls.addButton('button_background', 'use_background', 'toggle', 'b/g', 'b/g', '', 0.5, 'red');
            
            obj.panel_controls.addButton('button_preview', 'preview', 'push', 'PREVIEW', '', '', 0.5);
            obj.panel_controls.addButton('button_live', 'live', 'push', 'LIVE', '', '', 0.5);
            obj.panel_controls.make;
            obj.panel_controls.button_num_files.Callback = @obj.callback_num_files;
            
            %%%%%%%%%%% panel contrast %%%%%%%%%%%%%%%
            
            N = 5;
            pos = pos - N;
            obj.panel_contrast = util.plot.ContrastLimits(obj.axes_image, obj.fig.fig, [0 pos/N_left W_left N/N_left]);
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Position', [0 0 W_left 1/N_left]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RIGHT SIDE %%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%% panel objects %%%%%%%%%%%%%%%%
            
            pos = N_right;
            N = 12;
            pos = pos - N;
            obj.panel_objects = GraphicPanel(obj.owner, [1-W_right pos/N_right W_right N/N_right], 'objects');
            obj.panel_objects.number = N;
            obj.panel_objects.addButton('button_pars', 'pars', 'push', 'Parameters GUI');
            obj.panel_objects.addButton('button_buffers', 'buffers', 'push', 'Buffers GUI');
            obj.panel_objects.addButton('button_reader', 'reader', 'push', 'Reader GUI');
            obj.panel_objects.addButton('button_calibration', 'cal', 'push', 'Calibration GUI');
            obj.panel_objects.addButton('button_background', 'back', 'push', 'Background GUI');
            obj.panel_objects.addButton('button_clipper', 'clip', 'push', 'Clipper GUI');
            obj.panel_objects.addButton('button_clipper_bg', 'clip_bg', 'push', 'Clipper BG GUI');
            obj.panel_objects.addButton('button_photometry', 'phot', 'push', 'Photometry GUI');
            obj.panel_objects.addButton('button_phot_stack', 'phot_stack', 'push', 'stack Phot GUI');
            obj.panel_objects.addButton('button_lightcurves', 'lightcurves', 'push', 'Lightcurves GUI');
            obj.panel_objects.addButton('button_deflator', 'deflator', 'push', 'Deflator GUI');
            obj.panel_objects.make;
            
            %%%%%%%%%%% panel save %%%%%%%%%%%%%%%%%%%
            
            N = 3;
            pos = pos - N;
            obj.panel_save = GraphicPanel(obj.owner, [1-W_right pos/N_right W_right N/N_right], 'save');
            obj.panel_save.number = N;
            obj.panel_save.addButton('button_save', 'use_save', 'toggle', 'save', 'save', '', 1, 'red');
            obj.panel_save.addButton('button_trig', 'use_triggered_save', 'toggle', 'trig_save', 'trig_save', '', 1, 'red');
            
            obj.panel_save.make;
            
            %%%%%%%%%%% panel run %%%%%%%%%%%%%%%%%%
            
            obj.panel_run = uipanel('Title','', 'Position', [W_left, 0, 1-W_left-W_right, 0.1]);
            obj.button_run = GraphicButton(obj.panel_run, [0 0 1 1], obj.owner, '', 'custom', 'RUN');
            obj.button_run.Callback = @obj.callback_run;
            
            %%%%%%%%%%% panel info %%%%%%%%%%%%%%%%%%
            
            obj.panel_info = GraphicPanel(obj.owner, [W_left 0.9 1-W_left-W_right 0.1], 'info', 0);
            obj.panel_info.addButton('button_frame_rate', 'frame_rate_average', 'info', 'f= ', ' Hz'); 
            obj.panel_info.addButton('button_width', 'average_width', 'info', 'width= ', ' px', 'small', 0.5);
            obj.panel_info.addButton('button_seeing', 'average_width', 'custom', 'seeing= ', '"', 'small', 0.5);
            obj.panel_info.addButton('button_offset_y', 'average_offsets', 'custom', 'dy= ', '', 'small', 0.5);
            obj.panel_info.addButton('button_offset_x', 'average_offsets', 'custom', 'dx= ', '', 'small', 0.5);
            obj.panel_info.addButton('button_flux', 'average_flux', 'info', 'flux= ', '', 'small', 0.5);
            obj.panel_info.addButton('button_bg', 'average_background', 'info', 'b/g= ', '', 'small', 0.5);
            obj.panel_info.addButton('button_temperature', 'sensor_temperature', 'info', 's.temp= '); 
            
            obj.panel_info.make;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [W_left 0.1 1-W_left-W_right 0.8]);
                        
            obj.makeAxes;
            
            obj.button_batch_num = GraphicButton(obj.panel_image, [0 0.95 0.1 0.05], obj.owner, 'batch_counter', 'info', 'N= ', '', 'small');
            
            obj.button_show_what = GraphicButton(obj.panel_image, [0 0.00 0.15 0.05], obj.owner, 'show_what', 'picker', 'full', '', 'small');
            obj.button_show_what.Callback = @obj.callback_show_what;
            obj.button_show_what.String = obj.owner.show_what_list;
            
            obj.button_flip = GraphicButton(obj.panel_image, [0.9 0.00 0.1 0.05], obj.owner, 'use_flip', 'toggle', 'flip', '', 'small');
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.9 0.95 0.1 0.05], obj.owner, '', 'custom','reset');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            
            obj.update;
            
        end
            
        function makeAxes(obj, ~, ~)
            
            delete(obj.axes_image);
            
            obj.axes_image = axes('Parent', obj.panel_image);
            
            obj.panel_contrast.ax = obj.axes_image;
%             colorbar(obj.axes_image);
            axis(obj.axes_image, 'image');
            
            obj.panel_contrast.ax = obj.axes_image;
            
            if ~isempty(obj.owner.images)
                obj.owner.show;
            end
            
        end
                
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            if isa(obj.owner.src, 'file.Reader')
                obj.panel_controls.button_source_choose.String = 'src: Reader';
            elseif isa(obj.owner.src, 'img.Simulator')
                obj.panel_controls.button_source_choose.String = 'src: Simulator';
            elseif isa(obj.owner.src, 'obs.cam.CameraControl') || isa(obj.owner.src, 'obs.cam.Andor')
                obj.panel_controls.button_source_choose.String = 'src: Camera';
            end
            
            obj.panel_info.button_seeing.String = sprintf('seeing= %4.2f', obj.owner.average_width.*obj.owner.pars.plate_scale.*2.355);
            
%             obj.panel_info.button_offsets.String = sprintf('dx,dy= %4.2f,%4.2f', obj.owner.average_offsets(2), obj.owner.average_offsets(1));
            if length(obj.owner.average_offsets)==2
                obj.panel_info.button_offset_x.String = sprintf('dx= %4.2f', obj.owner.average_offsets(2));
                obj.panel_info.button_offset_y.String = sprintf('dy= %4.2f', obj.owner.average_offsets(1));
            else
                obj.panel_info.button_offset_x.String = 'dx= ';
                obj.panel_info.button_offset_y.String = 'dy= ';
            end
            
            if obj.owner.brake_bit
                if obj.owner.start_index>1
                    obj.button_run.String = 'CONTINUE';
                else
                    obj.button_run.String = 'RUN';
                end
%                 obj.button_run.Enable = 'on';
            else
                obj.button_run.String = 'STOP';
%                 obj.button_run.Enable = 'off';
            end
            
            for ii = 1:length(obj.button_show_what.String)
                if util.text.cs(obj.owner.show_what, obj.button_show_what.String{ii})
                    obj.button_show_what.Value = ii;
                    break;
                end
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
        
        function callback_num_files(obj, ~, ~)
            
            
            if obj.debug_bit, disp('callback: num_files'); end
            
            if ~isinf(obj.owner.num_files)
                obj.owner.num_batches = obj.owner.num_files;
            end
            
            obj.update;
        
        end
        
        function callback_show_what(obj, hndl, ~)
            
            if obj.debug_bit, disp('callback: show_what'); end
            
%             obj.owner.show_what_index = hndl.Value;
%             obj.owner.show_what = obj.owner.show_what_list{obj.owner.show_what_index};

            obj.owner.show_what = hndl.String{hndl.Value};

            obj.owner.show;
            obj.panel_contrast.autodyn;
            
            obj.update;
            
        end
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end