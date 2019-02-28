classdef ReaderGUI < handle
    
    properties 
        
        owner@file.Reader; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 16;
        edit_font_size = 14;
        small_font_size = 12;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        % left side
        panel_files;
        panel_control;
        panel_contrast;
        
        % right side
        panel_limits;
        panel_tbd;
        
        % top/bottom panels
        panel_info;
        panel_stop;
        button_stop;
        
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
       
        function obj = ReaderGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'file.Reader')
                
                if obj.debug_bit, fprintf('ReaderGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input a file.Reader to constructor of ReaderGUI!');
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
                obj.fig = util.plot.FigHandler('reader');
            end
            
            obj.fig.reset;
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 25;
            movegui(obj.fig.fig, 'center');
            
            %%%%%%%%%%% panel files %%%%%%%%%%%%%%%%%%

            obj.panel_files = GraphicPanel(obj.owner, [0 9/12 0.2 3/12], 'files');
            obj.panel_files.addButton('button_browse', 'browseDir', 'push', 'browse');
            obj.panel_files.addButton('button_glob', '', 'custom', 'glob: ', '','',0.3);
            obj.panel_files.addButton('input_glob', 'glob_string', 'input_text', ' ','','',0.7);
            obj.panel_files.addButton('button_num_files', 'num_files', 'info', 'N= ');
            obj.panel_files.make;
            
            %%%%%%%%%%% panel control %%%%%%%%%%%%%%%%
            
            obj.panel_control = GraphicPanel(obj.owner, [0 6/12 0.2 3/12], 'control');
            obj.panel_control.addButton('button_reset', 'reset', 'push', 'RESET');
            obj.panel_control.addButton('button_quick', 'quickScan', 'push', 'quick scan');
            obj.panel_control.addButton('button_batch', 'batch', 'push', 'single batch');
            obj.panel_control.make;
            
            %%%%%%%%%%% panel contrast %%%%%%%%%%%%%%%

            obj.panel_contrast = ContrastLimits(obj.axes_image, obj.fig.fig, [0 1/12 0.2 5/12], 1);
            
            %%%%%%%%%%% panel limits %%%%%%%%%%%%%%%%%

            obj.panel_limits = GraphicPanel(obj.owner, [0.8 3/12 0.2 9/12], 'limits');
            obj.panel_limits.addButton('input_num_batches', 'num_batches', 'input', 'num_batch= ');
            obj.panel_limits.addButton('input_files_per_batch', 'num_files_per_batch', 'input', 'files/batch= ');
            obj.panel_limits.addButton('input_file1', 'file_index_start', 'input', 'file1= ', '','',0.5);
            obj.panel_limits.addButton('input_file2', 'file_index_finish', 'input', 'file2= ', '','',0.5);
            obj.panel_limits.addButton('button_wrap', 'use_wrap_around', 'toggle', 'wrap around');
            obj.panel_limits.addButton('input_frames_per_batch', 'num_frames_per_batch', 'input', 'frames/batch= ');
            obj.panel_limits.addButton('input_frame1', 'frame_index_start', 'input', 'frame1= ');
            obj.panel_limits.addButton('input_frame2', 'frame_index_finish', 'input', 'frame2= ');
            obj.panel_limits.addButton('input_left', 'AOI_left', 'input', 'left= ', '', '', 0.5);
            obj.panel_limits.addButton('input_width', 'AOI_width', 'input', 'width= ', '', '', 0.5);
            obj.panel_limits.addButton('input_top', 'AOI_top', 'input', 'top= ', '', '', 0.5);
            obj.panel_limits.addButton('input_height', 'AOI_height', 'input', 'height= ', '', '', 0.5);
            obj.panel_limits.addButton('button_reset_AOI', 'resetAOI', 'push', 'reset AOI');
            obj.panel_limits.make;
            
            %%%%%%%%%%% panel limits %%%%%%%%%%%%%%%%%
            
            obj.panel_tbd = GraphicPanel(obj.owner, [0.8 0/12 0.2 3/12], '');
            obj.panel_tbd.number = 3;
            obj.panel_tbd.make;
            
            %%%%%%%%%%% panel info %%%%%%%%%%%%%%%%%%%

            obj.panel_info = GraphicPanel(obj.owner, [0.2 0.9 0.6 0.1], 'info');
            obj.panel_info.addButton('button_folder', 'current_dir', 'info', 'dir: ','', 'small');
            obj.panel_info.addButton('button_filename', 'shortname', 'info', 'file: ', '', 'small');
            obj.panel_info.make;            
            
            %%%%%%%%%%% panel stop %%%%%%%%%%%%%%%%%%%

            obj.panel_stop = uipanel('Position', [0.2 0 0.6 0.1]);
            obj.button_stop = GraphicButton(obj.panel_stop, [0 0 1 1], obj.owner, '', 'custom', 'RUN');
            obj.button_stop.Callback = @obj.callback_stop;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 0.1 0.6 0.8]);
                        
            obj.makeAxes;
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.9 0.95 0.1 0.05], obj.owner, '', 'custom','reset');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Position', [0 0 0.2 1/12]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            obj.update;
            
        end
            
        function makeAxes(obj, ~, ~)
            
            delete(obj.axes_image);
            
            obj.axes_image = axes('Parent', obj.panel_image);
            
            obj.panel_contrast.ax = obj.axes_image;
            
            axis(obj.axes_image, 'image');
            
            obj.panel_contrast.ax = obj.axes_image;
            
            im = findobj(obj.axes_image, 'type', 'image');
            if isempty(im)
                im = imagesc(zeros(512));
                axis image;
            end
            
            colorbar(obj.axes_image);
            
            if ~isempty(obj.owner.images)
                im.CData = obj.owner.images(:,:,1);
            elseif ~isempty(obj.owner.psfs)
                im.CData = obj.owner.psfs(:,:,1);
            end
                        
        end
                
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            if obj.owner.brake_bit
                obj.button_stop.String = 'LOOP';
            else
                obj.button_stop.String = 'STOP';
            end
            
            this_num = min(obj.owner.num_files, obj.owner.this_file_index);
            obj.panel_files.button_num_files.String = util.text.print_vec([this_num, obj.owner.num_files], '/');
            
            if ~isempty(obj.owner.temp_num_files_per_batch) && obj.owner.temp_num_files_per_batch==1 && ...
                    ~isempty(obj.owner.temp_frame_index_start) && obj.owner.temp_frame_index_start==1 && ...
                    ~isempty(obj.owner.temp_frame_index_finish) && obj.owner.temp_frame_index_finish==1
                
                obj.panel_control.button_quick.control.BackgroundColor = [0.7 0.3 0.1];
                
            else
                obj.panel_control.button_quick.control.BackgroundColor = [0.94 0.94 0.94];
            end

            if isempty(obj.owner.this_filename)
                [~,a,b] = fileparts(obj.owner.prev_filename);
                obj.panel_info.button_filename.String = ['file: ' a b];
            end
            
        end
                        
        function c = check(obj)
           
            c =  ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_stop(obj, ~, ~)
            
            if obj.debug_bit, disp(['callback: stop. brake_bit= ' num2str(obj.owner.brake_bit)]); end
            
            if obj.owner.brake_bit
                obj.owner.brake_bit = 0;
                obj.owner.loop;
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