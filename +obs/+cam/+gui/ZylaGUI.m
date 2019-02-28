classdef ZylaGUI <  handle
   
    properties
       
        cam@obs.ZylaControl;
        fig@util.FigHandler;        
             
    end
    
    properties % graphic objects
    
        panel_counters;
        button_accumulate;
        button_acquiring;
                
        panel_AOI;
        button_AOI_size;
        button_AOI_stride;
        button_zoom;
        button_AOI_center;
        
        panel_cooling;
        button_fan_speed; % input
        button_sensor_cooling; % input
        button_sensor_temp; 
        
        panel_sensor;
        button_pixel_size;
        button_sensor_size;
        
        panel_data_flow;
        button_baseline;
        button_bit_depth;
        button_pixel_encoding;
        button_transfer_rate;
        button_gain_control; % input
        button_readout_time;

        panel_parameters;
        button_exp_time; % input
        button_frame_count; % input
        button_frame_rate; % input
                
        panel_shutter;
        button_overlap; % input
        button_cycle_mode; % input
        button_trigger_mode; % input
        button_shutter_mode; % input 
        button_global_clear; % input 
        
        panel_metadata;
        button_metadata_enable; % input 
%         button_metadata_timestamps; % input
        
        panel_processing;
        button_noise_filter; % input
        button_blemish; % input
        
        panel_additional;
%         button_busy;
        
        panel_info;
        button_serial;
        button_version;
        button_interface;
        
        panel_close;
        button_close;
        
    end
    
    properties % switches and controls
       
        font_size = 12;
        
    end
    
    methods % constructor
       
        function obj = ZylaGUI(other)
            
            if nargin>0 && ~isempty(other) && isa(other, 'obs.ZylaControl')
                obj.cam = other;
            else
                error('need a obs.ZylaControl object as input');
            end
            
        end
        
    end
    
    methods % make gui
        
        function makeGUI(obj)
        
            obj.fig = util.FigHandler(obj.fig);
            obj.fig.name = 'Zyla GUI';
            obj.fig.reset;
        
            N = 7; % number of rows
            width = 0.5; % always 2 columns... 
            height = 1/N;
            bottom = 1-height;
            left = 0;
            
%             obj.panel_counters = uipanel('Title', 'counters', 'Position', [left, bottom, width, height]);
%             
%             obj.button_accumulate = uicontrol(obj.panel_counters, 'style', 'pushbutton',...
%                 'Units', 'Normalized', 'Position', [0, 0, 0.8, 1], 'FontSize', obj.font_size);
%             
%             obj.button_acquiring = uicontrol(obj.panel_counters, 'style', 'pushbutton',...
%                 'Units', 'Normalized', 'Position', [0.8, 0, 0.2, 1], 'FontSize', obj.font_size);
            
            %%%%%%%%%%%%%%%%%%%%%%% AOI panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            bottom = bottom - height;
            obj.panel_AOI = uipanel('Title', 'AOI', 'Position', [left, bottom, width, 2*height]);
            
            obj.button_AOI_size = uicontrol(obj.panel_AOI, 'style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [0, 0.5, 0.4, 0.5], 'FontSize', obj.font_size, ...
                'Callback', @obj.update, 'Enable', 'off');
            
            obj.button_AOI_stride = uicontrol(obj.panel_AOI, 'style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [0.4, 0.5, 0.4, 0.5], 'FontSize', obj.font_size,...
                'Callback', @obj.update, 'Enable', 'off');
            
            obj.button_zoom = uicontrol(obj.panel_AOI, 'style', 'toggle',...
                'Units', 'Normalized', 'Position', [0.8, 0.5, 0.2, 0.5], 'FontSize', obj.font_size, ...
                'Callback', @obj.callback_zoom);
            
            obj.button_AOI_center = uicontrol(obj.panel_AOI, 'style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [0, 0, 1, 0.5], 'FontSize', obj.font_size, ...
                'Callback', @obj.update, 'Enable', 'off');
            
            %%%%%%%%%%%%%%%%%%%%%%% cooling panel %%%%%%%%%%%%%%%%%%%%%%%%%
            bottom = bottom - height;  
            obj.panel_cooling = uipanel('Title', 'cooling', 'Position', [left, bottom, width, height]);
                    
            obj.button_fan_speed = uicontrol(obj.panel_cooling, 'style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [0/3, 0, 1/3, 1], 'FontSize', obj.font_size, ...
                'Callback', @obj.callback_fan_speed);
            
            obj.button_sensor_cooling = uicontrol(obj.panel_cooling, 'style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [1/3, 0, 1/3, 1], 'FontSize', obj.font_size, ...
                'Callback', @obj.callback_sensor_cooling);
            
            obj.button_sensor_temp = uicontrol(obj.panel_cooling, 'style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [2/3, 0, 1/3, 1], 'FontSize', obj.font_size,...
                'Callback', @obj.update, 'Enable', 'off');
            
            %%%%%%%%%%%%%%%%%%%%%%% sensor panel %%%%%%%%%%%%%%%%%%%%%%%%%%
            bottom = bottom - height;  
            obj.panel_sensor = uipanel('Title', 'sensor', 'Position', [left, bottom, width, height]);
            
            obj.button_sensor_size = uicontrol(obj.panel_sensor, 'style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [0/2, 0, 1/2, 1], 'FontSize', obj.font_size,...
                'Callback', @obj.update, 'Enable', 'off');
            
            obj.button_pixel_size = uicontrol(obj.panel_sensor, 'style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [1/2, 0, 1/2, 1], 'FontSize', obj.font_size,...
                'Callback', @obj.update, 'Enable', 'off');
                        
            %%%%%%%%%%%%%%%%%%%%%%% data flow panel %%%%%%%%%%%%%%%%%%%%%%%
            bottom = bottom - 2*height;  
            obj.panel_data_flow = uipanel('Title', 'data flow', 'Position', [left, bottom, width, 2*height]);
                                
            obj.button_baseline = uicontrol(obj.panel_data_flow, 'style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [0/3, 0.5, 1/3, 0.5], 'FontSize', obj.font_size,...
                'Callback', @obj.update, 'Enable', 'off');
            
            obj.button_bit_depth = uicontrol(obj.panel_data_flow, 'style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [1/3, 0.5, 1/3, 0.5], 'FontSize', obj.font_size,...
                'Callback', @obj.update, 'Enable', 'off');
            
            obj.button_pixel_encoding = uicontrol(obj.panel_data_flow, 'style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [2/3, 0.5, 1/3, 0.5], 'FontSize', obj.font_size,...
                'Callback', @obj.update, 'Enable', 'off');
            
            obj.button_transfer_rate = uicontrol(obj.panel_data_flow, 'style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [0/3, 0, 1/3, 0.5], 'FontSize', obj.font_size,...
                'Callback', @obj.update, 'Enable', 'off');
            
            obj.button_readout_time = uicontrol(obj.panel_data_flow, 'style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [1/3, 0, 1/3, 0.5], 'FontSize', obj.font_size,...
                'Callback', @obj.update, 'Enable', 'off');
            
            obj.button_gain_control = uicontrol(obj.panel_data_flow, 'style', 'popupmenu',...
                'String', {'16bit low-noise & high well', '11bit low noise', '11bit high well'},...
                'Units', 'Normalized', 'Position', [2/3, 0, 1/3, 0.5], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_gain_control);
            
            %%%%%%%%%%%%%%%%%%%%%%% parameters panel %%%%%%%%%%%%%%%%%%%%%%
            bottom = bottom - height;  
            obj.panel_parameters = uipanel('Title', 'parameters', 'Position', [left, bottom, width, height]);
                        
            obj.button_exp_time = uicontrol(obj.panel_parameters, 'Style', 'edit',...
                'Units', 'Normalized', 'Position', [0/3, 0, 1/3, 1], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_exp_time);
            
            obj.button_frame_rate = uicontrol(obj.panel_parameters, 'Style', 'edit',...
                'Units', 'Normalized', 'Position', [1/3, 0, 1/3, 1], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_frame_rate);
            
            obj.button_frame_count = uicontrol(obj.panel_parameters, 'Style', 'edit',...
                'Units', 'Normalized', 'Position', [2/3, 0, 1/3, 1], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_frame_count);
            
            %%%%%%%%%%%%%%%%%%%%%%% shutter panel %%%%%%%%%%%%%%%%%%%%%%%%%
            left = 0.5;
            bottom = 1 - 2*height;  
            obj.panel_shutter = uipanel('Title', 'shutter', 'Position', [left, bottom, width, 2*height]);
            
            obj.button_cycle_mode = uicontrol(obj.panel_shutter, 'Style', 'popupmenu',...
                'Units', 'Normalized', 'Position', [0/3, 0.5, 1/3, 0.5], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_cycle_mode, 'String', {'fixed', 'continuous'});
                        
            obj.button_trigger_mode = uicontrol(obj.panel_shutter, 'Style', 'popupmenu',...
                'Units', 'Normalized', 'Position', [1/3, 0.5, 1/3, 0.5], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_trigger_mode, 'String', {'internal', 'software', 'external', 'external start', 'external exposure'});
            
            obj.button_shutter_mode = uicontrol(obj.panel_shutter, 'Style', 'popupmenu',...
                'Units', 'Normalized', 'Position', [2/3, 0.5, 1/3, 0.5], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_shutter_mode, 'String', {'global', 'rolling'});
            
            obj.button_overlap = uicontrol(obj.panel_shutter, 'Style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [0/3, 0, 1/3, 0.5], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_overlap);
            
            obj.button_global_clear = uicontrol(obj.panel_shutter, 'Style', 'push',...
                'Units', 'Normalized', 'Position', [1/3, 0, 1/3, 0.5], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_global_clear);
                        
            %%%%%%%%%%%%%%%%%%%%%%% metadata panel %%%%%%%%%%%%%%%%%%%%%%%%
            
            bottom = bottom - height;  
            obj.panel_metadata = uipanel('Title', 'metadata', 'Position', [left, bottom, width, height]);
            
            obj.button_metadata_enable = uicontrol(obj.panel_metadata, 'Style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [0/2, 0, 1/2, 1], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_metadata_enable);
            
%             obj.button_metadata_timestamps = uicontrol(obj.panel_metadata, 'Style', 'pushbutton',...
%                 'Units', 'Normalized', 'Position', [1/2, 0, 1/2, 1], 'FontSize', obj.font_size,...
%                 'Callback', @obj.callback_metadata_timestamps);
                        
            %%%%%%%%%%%%%%%%%%%%%%% additional panel %%%%%%%%%%%%%%%%%%%%%%
            bottom = bottom - height;  
            obj.panel_additional = uipanel('Title', 'additional', 'Position', [left, bottom, width, height]);
                              
            %%%%%%%%%%%%%%%%%%%%%%% processing panel %%%%%%%%%%%%%%%%%%%%%%
            bottom = bottom - height;  
            obj.panel_processing = uipanel('Title', 'processing', 'Position', [left, bottom, width, height]);
                          
            obj.button_noise_filter = uicontrol(obj.panel_processing, 'Style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [0/2, 0, 1/2, 1], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_noise_filter);
            
            obj.button_blemish = uicontrol(obj.panel_processing, 'Style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [1/2, 0, 1/2, 1], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_blemish);
            
            %%%%%%%%%%%%%%%%%%%%%%% info panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            bottom = bottom - height;  
            obj.panel_info = uipanel('Title', 'info', 'Position', [left, bottom, width, height]);
                    
            [rc, val] = AT_GetString(obj.cam.hndl, 'SerialNumber', 255); AT_CheckWarning(rc);
            
            obj.button_serial = uicontrol(obj.panel_info, 'style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [0/3, 0, 1/3, 1], 'FontSize', obj.font_size,...
                'Callback', @obj.update, 'String', ['S/N: ' val], 'Enable', 'off');
            
            [rc, val] = AT_GetString(1, 'SoftwareVersion', 255); AT_CheckWarning(rc);
            
            obj.button_version = uicontrol(obj.panel_info, 'style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [1/3, 0, 1/3, 1], 'FontSize', obj.font_size,...
                'Callback', @obj.update, 'String', ['version: ' val], 'Enable', 'off');
            
            [rc, val] = AT_GetString(obj.cam.hndl, 'InterfaceType', 255); AT_CheckWarning(rc);
            
            obj.button_interface = uicontrol(obj.panel_info, 'style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [2/3, 0, 1/3, 1], 'FontSize', obj.font_size,...
                'Callback', @obj.update, 'String', ['interface: ' val], 'Enable', 'off');
            
            %%%%%%%%%%%%%%%%%%%%%%% close panel %%%%%%%%%%%%%%%%%%%%%%%%%%%
            bottom = bottom - height;  
            obj.panel_close = uipanel('Title', 'close', 'Position', [left, bottom, width, height]);
            
            obj.button_close = uicontrol(obj.panel_close, 'Style', 'pushbutton',...
                'Units', 'Normalized', 'Position', [0, 0, 1, 1], 'FontSize', obj.font_size,...
                'Callback', @obj.callback_close, 'String', 'CLOSE');

            
            obj.update;
            
        end
            
    end
    
    methods % update gui / check gui
        
        function c = check(obj)
           
            c = ~isempty(obj.fig) && obj.fig.check && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
        function update(obj, ~, ~)
           
            import util.text.cs;
            import util.text.f2s;
            
            % counters panel
%             [rc, val] = AT_GetInt(obj.cam.hndl, 'AccumulateCount'); AT_CheckWarning(rc);
%             obj.button_accumulate.String = ['accumulate count= ' num2str(val)];
            
%             [rc, val] = AT_GetInt(obj.cam.hndl, 'CameraAcquiring'); AT_CheckWarning(rc);
%             obj.button_acquiring.String = ['acquiring= ' num2str(val)];
            
            
            %%%%%%%%%%%%%%%%%%%%%%% AOI panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            [rc, width] = AT_GetInt(obj.cam.hndl, 'AOIWidth'); AT_CheckWarning(rc);
            [rc, height] = AT_GetInt(obj.cam.hndl, 'AOIHeight'); AT_CheckWarning(rc);
            obj.button_AOI_size.String = [num2str(width) 'x' num2str(height)];
            
            [rc, left] = AT_GetInt(obj.cam.hndl, 'AOILeft'); AT_CheckWarning(rc);
            [rc, top] = AT_GetInt(obj.cam.hndl, 'AOITop'); AT_CheckWarning(rc);
            obj.button_AOI_center.String = ['center: (' num2str(floor(left+width/2)) ',' num2str(floor(top+height/2)) ')'];
                        
            [rc, val] = AT_GetInt(obj.cam.hndl, 'AOIStride'); AT_CheckWarning(rc);
            obj.button_AOI_stride.String = ['AOI stride: ' num2str(val)];
            
            if width<obj.cam.maxWidth || height<obj.cam.maxHeight
                obj.button_zoom.String = 'unzoom';
            else
                obj.button_zoom.String = 'zoom';
            end
                                
            %%%%%%%%%%%%%%%%%%%%%%% cooling panel %%%%%%%%%%%%%%%%%%%%%%%%%
            
            [rc, val] = AT_GetEnumIndex(obj.cam.hndl, 'FanSpeed'); AT_CheckWarning(rc);
            if val
                obj.button_fan_speed.String = 'fan on';
            else
                obj.button_fan_speed.String = 'fan off';
            end
            
            [rc, val] = AT_GetBool(obj.cam.hndl, 'SensorCooling'); AT_CheckWarning(rc);
            if val
                obj.button_sensor_cooling.String = 'cooling on';
            else
                obj.button_sensor_cooling.String = 'cooling off';
            end
            
            [rc, val] = AT_GetFloat(obj.cam.hndl, 'SensorTemperature'); AT_CheckWarning(rc);
            obj.button_sensor_temp.String = [num2str(val) 'C'];
            
            %%%%%%%%%%%%%%%%%%%%%%% sensor panel %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            [rc, s_height] = AT_GetInt(obj.cam.hndl, 'SensorHeight'); AT_CheckWarning(rc);            
            [rc, s_width] = AT_GetInt(obj.cam.hndl, 'SensorWidth'); AT_CheckWarning(rc);
            [rc, p_height] = AT_GetFloat(obj.cam.hndl, 'PixelHeight'); AT_CheckWarning(rc);            
            [rc, p_width] = AT_GetFloat(obj.cam.hndl, 'PixelWidth'); AT_CheckWarning(rc);
            
            obj.button_sensor_size.String = ['size: ' num2str(p_width*s_width*1e-3) 'x' num2str(p_height*s_height*1e-3) 'mm'];
            
            obj.button_pixel_size.String = ['pixel size: ' num2str(p_width) 'x' num2str(p_height) 'um'];
            
            %%%%%%%%%%%%%%%%%%%%%%% data flow panel %%%%%%%%%%%%%%%%%%%%%%%
                        
            [rc, val] = AT_GetInt(obj.cam.hndl, 'Baseline'); AT_CheckWarning(rc); 
            obj.button_baseline.String = ['baseline: ' num2str(val)];
            
            [rc, val] = AT_GetEnumIndex(obj.cam.hndl, 'BitDepth'); AT_CheckWarning(rc); 
            [rc, val] = AT_GetEnumStringByIndex(obj.cam.hndl, 'BitDepth', val, 255); AT_CheckWarning(rc); 
            obj.button_bit_depth.String = ['depth: ' val];
            
            [rc, val] = AT_GetEnumIndex(obj.cam.hndl, 'PixelEncoding'); AT_CheckWarning(rc); 
            [rc, val] = AT_GetEnumStringByIndex(obj.cam.hndl, 'PixelEncoding', val, 255); AT_CheckWarning(rc); 
            obj.button_pixel_encoding.String = ['encoding: ' val];
                                    
            [rc, val] = AT_GetFloat(obj.cam.hndl, 'MaxInterfaceTransferRate'); AT_CheckWarning(rc); 
            obj.button_transfer_rate.String = ['rate: ' f2s(val) 'MB/s'];
            
            [rc, val] = AT_GetFloat(obj.cam.hndl, 'ReadoutTime'); AT_CheckWarning(rc); 
            obj.button_readout_time.String = ['read time: ' f2s(val)];
            
            [rc, val] = AT_GetEnumIndex(obj.cam.hndl, 'SimplePreAmpGainControl'); AT_CheckWarning(rc);
            obj.button_gain_control.Value = 3-val;
%             [rc, str] = AT_GetEnumStringByIndex(obj.cam.hndl, 'SimplePreAmpGainControl', val, 255); AT_CheckWarning(rc);

            
            %%%%%%%%%%%%%%%%%%%%%%% parameters panel %%%%%%%%%%%%%%%%%%%%%%
            obj.button_exp_time.String = ['expT= ' num2str(obj.cam.getExpTime)];
            obj.button_frame_rate.String = ['frame rate= ' num2str(obj.cam.getFrameRate)];
            [rc, val] = AT_GetInt(obj.cam.hndl, 'FrameCount'); AT_CheckWarning(rc);
            obj.button_frame_count.String = ['frame count= ' num2str(val)];

            %%%%%%%%%%%%%%%%%%%%%%% shutter panel %%%%%%%%%%%%%%%%%%%%%%%%%

            [rc, val] = AT_GetEnumIndex(obj.cam.hndl, 'CycleMode'); AT_CheckWarning(rc);
            [rc, str] = AT_GetEnumStringByIndex(obj.cam.hndl, 'CycleMode', val, 255); AT_CheckWarning(rc);
            for ii = 1:length(obj.button_cycle_mode.String)
                if cs(str, obj.button_cycle_mode.String{ii})
                    obj.button_cycle_mode.Value = ii;
                    break;
                end
            end
            
            [rc, val] = AT_GetEnumIndex(obj.cam.hndl, 'TriggerMode'); AT_CheckWarning(rc);
            [rc, str] = AT_GetEnumStringByIndex(obj.cam.hndl, 'TriggerMode', val, 255); AT_CheckWarning(rc);
            for ii = 1:length(obj.button_trigger_mode.String)
                if cs(str, obj.button_trigger_mode.String{ii})
                    obj.button_trigger_mode.Value = ii;
                    
                end
            end
            
            [rc, val] = AT_GetEnumIndex(obj.cam.hndl, 'ElectronicShutteringMode'); AT_CheckWarning(rc);
            [rc, str] = AT_GetEnumStringByIndex(obj.cam.hndl, 'ElectronicShutteringMode', val, 255); AT_CheckWarning(rc);
            for ii = 1:length(obj.button_shutter_mode.String)
                if cs(str, obj.button_shutter_mode.String{ii})
                    obj.button_shutter_mode.Value = ii;
                    break;
                end
            end
            
            [rc, val] = AT_GetBool(obj.cam.hndl, 'Overlap'); AT_CheckWarning(rc);
            if val
                obj.button_overlap.String = 'overlap';
            else                
                obj.button_overlap.String = 'no overlap';
            end
            
            %%%%%%%%%%%%%%%%%%%%%%% metadata panel %%%%%%%%%%%%%%%%%%%%%%%%
            
            [rc, val] = AT_GetBool(obj.cam.hndl, 'MetadataEnable'); AT_CheckWarning(rc);
            if val
                obj.button_metadata_enable.String = 'metadata enabled';
            else
                obj.button_metadata_enable.String = 'metadata disabled';
            end
            
%             if meta_enabled
%                 
%                 [rc, val] = AT_GetBool(obj.cam.hndl, 'MetadataTimestamps'); AT_CheckWarning(rc);
%                 if val
%                     obj.button_metadata_timestamps.String = 'timestamps enabled';
%                 else
%                     obj.button_metadata_timestamps.String = 'timestamps disabled';
%                 end
%                 
%             end
                        
            %%%%%%%%%%%%%%%%%%%%%%% processing panel %%%%%%%%%%%%%%%%%%%%%%
                        
            [rc, val] = AT_GetBool(obj.cam.hndl, 'SpuriousNoiseFilter'); AT_CheckWarning(rc);
            if val
                obj.button_noise_filter.String = 'noise filter enabled';
            else
                obj.button_noise_filter.String = 'noise filter disabled';
            end 
            
            [rc, val] = AT_GetBool(obj.cam.hndl, 'StaticBlemishCorrection'); AT_CheckWarning(rc);
            if val
                obj.button_blemish.String = 'blemish correction enabled';
            else
                obj.button_blemish.String = 'blemish correction disabled';
            end
            
            %%%%%%%%%%%%%%%%%%%%%%% info panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % updated once inside the makeGUI function...
            
            
        end
        
    end
    
    methods % callbacks
        
        function callback_zoom(obj, hndl, ~)
           
             val = ~hndl.Value;
                         
             if ~val
                obj.cam.zoom;
            else
                obj.cam.full_frame;
            end
            
            obj.button_zoom.Value = hndl.Value;
            
            obj.update;
            
        end
                
        function callback_fan_speed(obj, ~, ~)
                                
            [rc, val] = AT_GetEnumIndex(obj.cam.hndl, 'FanSpeed'); AT_CheckWarning(rc);
            
            if val
                rc = AT_SetEnumIndex(obj.cam.hndl, 'FanSpeed', 0); AT_CheckWarning(rc);
            else                
                rc = AT_SetEnumIndex(obj.cam.hndl, 'FanSpeed', 1); AT_CheckWarning(rc);
            end
            obj.update;
            
        end
                
        function callback_sensor_cooling(obj, ~, ~)
                
            [rc, val] = AT_GetBool(obj.cam.hndl, 'SensorCooling'); AT_CheckWarning(rc);
            
            if val
                rc = AT_SetBool(obj.cam.hndl, 'SensorCooling', 0); AT_CheckWarning(rc);
            else
                rc = AT_SetBool(obj.cam.hndl, 'SensorCooling', 1); AT_CheckWarning(rc);
            end
            
            obj.update;
            
        end
        
        function callback_gain_control(obj, hndl, ~)
            
            rc = AT_SetEnumIndex(obj.cam.hndl, 'SimplePreAmpGainControl', 3-hndl.Value); AT_CheckWarning(rc);
            
            obj.update;
            
        end
        
        function callback_exp_time(obj, hndl, ~)
            
            val = util.text.extract_numbers(hndl.String);
            val = val{1};
            
            if ~isempty(val)
                obj.cam.setExpTime(val);
            end
            
            obj.update;
            
        end
        
        function callback_frame_rate(obj, hndl, ~)
            
            val = util.text.extract_numbers(hndl.String);
            val = val{1};
            
            if ~isempty(val)
                obj.cam.setFrameRate(val);
            end
            
            obj.update;
            
        end
        
        function callback_frame_count(obj, hndl, ~)
            
            val = util.text.extract_numbers(hndl.String);
            val = val{1};
            
            if ~isempty(val)
                rc = AT_SetInt(obj.cam.hndl, 'FrameCount', val); AT_CheckWarning(rc);
            end
            
            obj.update;
            
        end
        
        function callback_cycle_mode(obj, hndl, ~)
            
            if hndl.Value==1
                rc = AT_SetEnumString(obj.cam.hndl, 'CycleMode', 'Fixed'); AT_CheckWarning(rc);
            elseif hndl.Value==2
                rc = AT_SetEnumString(obj.cam.hndl, 'CycleMode', 'Continuous'); AT_CheckWarning(rc);
            end
            
            obj.update;
            
        end   
        
        function callback_trigger_mode(obj, hndl, ~)
            
            if hndl.Value==1
                rc = AT_SetEnumString(obj.cam.hndl, 'TriggerMode', 'Internal'); AT_CheckWarning(rc);
            elseif hndl.Value==2
                rc = AT_SetEnumString(obj.cam.hndl, 'TriggerMode', 'Software'); AT_CheckWarning(rc);
            elseif hndl.Value==3
                rc = AT_SetEnumString(obj.cam.hndl, 'TriggerMode', 'External'); AT_CheckWarning(rc);
            elseif hndl.Value==4
                rc = AT_SetEnumString(obj.cam.hndl, 'TriggerMode', 'External Start'); AT_CheckWarning(rc);
            elseif hndl.Value==5
                rc = AT_SetEnumString(obj.cam.hndl, 'TriggerMode', 'External Exposure'); AT_CheckWarning(rc);
            end
            
            obj.update;
            
        end
        
        function callback_shutter_mode(obj, hndl, ~)
            
            if hndl.Value==1
                rc = AT_SetEnumString(obj.cam.hndl, 'ElectronicShutteringMode', 'Global'); AT_CheckWarning(rc);
            elseif hndl.Value==2
                rc = AT_SetEnumString(obj.cam.hndl, 'ElectronicShutteringMode', 'Rolling'); AT_CheckWarning(rc);
            end
            
            obj.update;
            
        end  
        
        function callback_global_clear(obj, ~, ~)
            
            [rc, val] = AT_GetBool(obj.cam.hndl, 'RollingShutterGlobalClear'); AT_CheckWarning(rc);
            
            if val
                AT_SetBool(obj.cam.hndl, 'RollingShutterGlobalClear', 0); AT_CheckWarning(rc);
            else                
                AT_SetBool(obj.cam.hndl, 'RollingShutterGlobalClear', 1); AT_CheckWarning(rc);
            end
            
            obj.update;
            
        end
        
        function callback_overlap(obj, ~, ~)
            
            [rc, val] = AT_GetBool(obj.cam.hndl, 'Overlap'); AT_CheckWarning(rc);
            
            if val
                AT_SetBool(obj.cam.hndl, 'Overlap', 0); AT_CheckWarning(rc);
            else                
                AT_SetBool(obj.cam.hndl, 'Overlap', 1); AT_CheckWarning(rc);
            end
            
            obj.update;
            
        end
                
        function callback_metadata_enable(obj, ~, ~)

            [rc, val] = AT_GetBool(obj.cam.hndl, 'MetadataEnable'); AT_CheckWarning(rc);
            
            if val
                rc = AT_SetBool(obj.cam.hndl, 'MetadataEnable', 0); AT_CheckWarning(rc);
            else
                rc = AT_SetBool(obj.cam.hndl, 'MetadataEnable', 1); AT_CheckWarning(rc);
            end
            
            obj.update;
            
        end  
                        
        function callback_noise_filter(obj, ~, ~)

            [rc, val] = AT_GetBool(obj.cam.hndl, 'SpuriousNoiseFilter'); AT_CheckWarning(rc);
            
            if val
                rc = AT_SetBool(obj.cam.hndl, 'SpuriousNoiseFilter', 0); AT_CheckWarning(rc);
            else
                rc = AT_SetBool(obj.cam.hndl, 'SpuriousNoiseFilter', 1); AT_CheckWarning(rc);
            end
            
            obj.update;
            
        end  
                        
        function callback_blemish(obj, ~, ~)

            [rc, val] = AT_GetBool(obj.cam.hndl, 'StaticBlemishCorrection'); AT_CheckWarning(rc);
            
            if val
                rc = AT_SetBool(obj.cam.hndl, 'StaticBlemishCorrection', 0); AT_CheckWarning(rc);
            else
                rc = AT_SetBool(obj.cam.hndl, 'StaticBlemishCorrection', 1); AT_CheckWarning(rc);
            end
            
            obj.update;
            
        end  
        
        function callback_close(obj, ~, ~)
           
            delete(obj.fig.fig);
            
        end
        
%         function callback_metadata_timestamps(obj, ~, ~)
% 
%             [rc, val] = AT_GetBool(obj.cam.hndl, 'MetadataTimestamps'); AT_CheckWarning(rc);
%             
%             if val
%                 rc = AT_SetBool(obj.cam.hndl, 'MetadataTimestamps', 0); AT_CheckWarning(rc);
%             else
%                 rc = AT_SetBool(obj.cam.hndl, 'MetadataTimestamps', 1); AT_CheckWarning(rc);
%             end
%             
%             obj.update;
%             
%         end
                
    end
        
end