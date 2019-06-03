classdef DeflatorGUI < handle

    properties 
        
        def@file.Deflator; % link back to containg object
        fig@util.plot.FigHandler;
        buf_gui@file.gui.BufferGUI;
        
        panel_filename;
        button_filename;
                
        panel_objects;
        button_reader;
        button_buffers;
        button_src_dir;
        button_out_dir;
        button_out_dir_backup;
            
        panel_controls;
        button_backup;
        button_auto_delete;
        button_num_files;
        
        panel_run_stop;
        button_run_stop; 

        image_panel;
        image_axes; 
        
        font_size = 18;        
        
        debug_bit = 1;
        
    end
    
    properties (Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
       
        function obj = DeflatorGUI(parent)
            
            % later add other options like copy constructor
            
            if obj.debug_bit, fprintf('DeflatorGUI constructor v%4.2f\n', obj.version); end
            
            assert(isa(parent, 'file.Deflator'), 'Input a file.Deflator to constructor of DeflatorGUI!');
            
            obj.def = parent;
            
        end
        
    end
    
    methods % make/update gui
        
        function makeGUI(obj)
           
            obj.fig = util.plot.FigHandler('deflator gui');
            
            %%%%%%%%%%%%%%%%%%%%%%%% objects panel %%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_objects = uipanel('Title', 'objects', 'Position', [0 0 0.4 0.5], 'Units','Normalized');
            
            N = 5;
            
            obj.button_out_dir_backup = uicontrol(obj.panel_objects, 'Style', 'pushbutton', 'String', 'backup dir',...
                'Units', 'Normalized', 'Position', [0 0/N 1 1/N], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_out_dir_backup);
            
            obj.button_out_dir = uicontrol(obj.panel_objects, 'Style', 'pushbutton', 'String', 'output dir',...
                'Units', 'Normalized', 'Position', [0 1/N 1 1/N], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_out_dir);
            
            obj.button_src_dir = uicontrol(obj.panel_objects, 'Style', 'pushbutton', 'String', 'source dir',...
                'Units', 'Normalized', 'Position', [0 2/N 1 1/N], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_src_dir);
                        
            obj.button_buffers = uicontrol(obj.panel_objects, 'Style', 'pushbutton', 'String', 'Buffers',...
                'Units', 'Normalized', 'Position', [0 3/N 1 1/N], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_buffers);
                        
            obj.button_reader = uicontrol(obj.panel_objects, 'Style', 'pushbutton', 'String', 'Reader',...
                'Units', 'Normalized', 'Position', [0 4/N 1 1/N], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_reader);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%% controls panel %%%%%%%%%%%%%%%%%%%%%
                        
            obj.panel_controls = uipanel('Title', 'controls', 'Position', [0 0.5 0.4 0.5], 'Units','Normalized');
                        
            N = 5;
            
            obj.button_backup = uicontrol(obj.panel_controls, 'Style', 'pushbutton', 'String', 'backup ON',...
                'Units', 'Normalized', 'Position', [0 0/N 1 1/N], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_backup);
            
            obj.button_auto_delete = uicontrol(obj.panel_controls, 'Style', 'pushbutton', 'String', 'auto delete',...
                'Units', 'Normalized', 'Position', [0 1/N 1 1/N], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_auto_delete, 'Enable', 'off');
            
            obj.button_num_files = uicontrol(obj.panel_controls, 'Style', 'pushbutton', 'String', '',...
                'Units', 'Normalized', 'Position', [0 4/N 1 1/N], 'FontSize', obj.font_size);
                        
            
            %%%%%%%%%%%%%%%%%%%%%%%% file display panel %%%%%%%%%%%%%%%%%%%
                          
            obj.panel_filename = uipanel('Title', 'current file', 'Position', [0.4 0.9 0.6 0.1], 'Units','Normalized');
                      
            obj.button_filename = uicontrol(obj.panel_filename, 'Style', 'pushbutton', 'String', 'no file loaded yet',...
                'Units', 'Normalized', 'Position', [0 0 1 1], 'FontSize', obj.font_size);
            
            %%%%%%%%%%%%%%%%%%%%%%%% run/stop panel %%%%%%%%%%%%%%%%%%%%%%%
            
            obj.panel_run_stop = uipanel('Title', '', 'Position', [0.4 0 0.6 0.1], 'Units','Normalized');
                        
            obj.button_run_stop = uicontrol(obj.panel_run_stop, 'Style', 'pushbutton', 'String', 'RUN',...
                'Units', 'Normalized', 'Position', [0 0 1 1], 'FontSize', 26,...
                'Callback', @obj.callback_run_stop);            
            
            %%%%%%%%%%%%%%%%%%%%%%%%% image panel %%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.image_panel = uipanel('Title', 'image', 'Position', [0.4 0.1 0.6 0.8], 'Units', 'Normalized');
            
            obj.image_axes = axes('Parent', obj.image_panel);
            
            obj.updateGUI;
            
        end
        
        function updateGUI(obj)
           
            if isempty(obj.fig) || ~obj.fig.isvalid || isempty(obj.image_panel) || ~isvalid(obj.image_panel)
                return; % short circuit this function if there's no GUI
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%% objects panel %%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.button_src_dir.String = ['src:' obj.def.src_dir.pwd];
            obj.button_out_dir.String = ['out:' obj.def.out_dir.pwd];
            obj.button_out_dir_backup.String = ['bkp:' obj.def.out_dir_backup.pwd];
           
            %%%%%%%%%%%%%%%%%%%%%%%%%% controls panel %%%%%%%%%%%%%%%%%%%%%
            
            if obj.def.use_auto_backup
                obj.button_backup.String = 'backup ON';
            else
                obj.button_backup.String = 'backup OFF';
            end
                        
            if obj.def.use_auto_delete
                obj.button_auto_delete.String = 'delete ON';
                obj.button_auto_delete.BackgroundColor = 'red';
            else            
                obj.button_auto_delete.String = 'delete OFF';
                obj.button_auto_delete.BackgroundColor = [0.7 0.7 0.7];
            end
            
            if isempty(obj.def.reader.filenames)
                obj.button_num_files.String = 'file: ';
            else
                N = obj.def.reader.this_file_index;
                if N<=length(obj.def.reader.filenames)
                    obj.button_num_files.String = ['file: ' num2str(N) ' / ' num2str(length(obj.def.reader.filenames)) '      '];
                else
                    obj.button_num_files.String = ['file: ' num2str(length(obj.def.reader.filenames)) ' / ' num2str(length(obj.def.reader.filenames)) ' done!'];
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%% file display panel %%%%%%%%%%%%%%%%%%%
            
            if isempty(obj.def.reader.this_filename)
                obj.button_filename.String = 'No file loaded yet!';
            else
                [~, filename, ext] = fileparts(obj.def.reader.this_filename);
                filename = [filename ext];
                obj.button_filename.String = filename;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%% run/stop panel %%%%%%%%%%%%%%%%%%%%%%%
            
            if obj.def.brake_bit
                obj.button_run_stop.String = 'RUN';
            else
                obj.button_run_stop.String = 'STOP';
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%% image panel %%%%%%%%%%%%%%%%%%%%%%%%%
            
            if isempty(obj.image_axes) || ~isvalid(obj.image_axes)
                obj.image_axes = axes('Parent', obj.image_panel); 
            end
            
            if ~isempty(obj.def.reader.images)
                util.plot.show(obj.def.reader.images(:,:,1,1), 'ax', obj.image_axes, 'autodyn', 'on');
            end
            
            drawnow;
            
        end
        
    end
    
    methods % callbacks
    
          %%%%%%%%%%%%%%%%%%%%%%%%%% objects panel %%%%%%%%%%%%%%%%%%%%%%%%
          
          function callback_reader(obj, ~, ~)
              
              obj.def.reader.makeGUI;
              
              obj.updateGUI;
              
          end
          
          function callback_buffers(obj, ~, ~)
             
%               if isempty(obj.def.buffers.gui)
%                   obj.def.buffers.makeGUI;
%                   obj.def.buffers.gui.debug_bit = 0;
%               end
              
              obj.def.buffers.makeGUI;
              
              obj.updateGUI; 
              
          end
          
          function callback_src_dir(obj, ~, ~)
             
              obj.def.src_dir.browse;
              
              obj.updateGUI;
                            
          end
          
          function callback_out_dir(obj, ~, ~)
             
              obj.def.out_dir.browse;
              
              obj.updateGUI;
                            
          end
          
          function callback_out_dir_backup(obj, ~, ~)
             
              obj.def.out_dir_backup.browse;
              
              obj.updateGUI;
                            
          end
        
          %%%%%%%%%%%%%%%%%%%%%%%%%%%% controls panel %%%%%%%%%%%%%%%%%%%%%
          
          function callback_backup(obj, ~, ~)
             
              obj.def.use_auto_backup = ~obj.def.use_auto_backup;
              
              obj.updateGUI;
              
          end
          
          function callback_auto_delete(obj, ~, ~)
             
              obj.def.use_auto_delete = ~obj.def.use_auto_delete;
              
              obj.updateGUI;
              
          end
            
          %%%%%%%%%%%%%%%%%%%%%%%%%% run/stop panel %%%%%%%%%%%%%%%%%%%%%%%
          
          function callback_run_stop(obj, ~, ~)
             
              if obj.debug_bit, disp('callback: run/stop'); end
              
              if obj.def.brake_bit
                  
                  obj.def.run;
                  
              else
                  
                  obj.def.brake_bit = 1;
                  
              end
              
              obj.updateGUI;
              
          end
          
          
    end
        
end

