classdef BufferGUI < handle
   
    properties % handles and stuff
       
        buffers@file.BufferWheel;
        
        fig_handler;
        
        debug_bit = 1;
        
    end
    
    properties(Hidden=true)
       
        version = 1.02;
        
    end
        
    methods % constructor & master functions
        
        function obj = BufferGUI(other)
            if nargin>0 && isa(other, 'file.gui.BufferGUI')                
                if obj.debug_bit>1, fprintf('BufferGUI copy-constructor (nothing to do really...) v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(other);
            elseif nargin>0 && isa(other, 'file.BufferWheel')
                if obj.debug_bit>1, fprintf('BufferGUI BufferWheel constructor v%4.2f\n', obj.version); end
                obj.buffers = other;
            else
                if obj.debug_bit, fprintf('BufferGUI constructor v%4.2f\n', obj.version); end
            end
            
        end
        
        function make(obj)
            
            if isempty(obj.fig_handler)
                obj.fig_handler = util.plot.FigHandler;
            else
                obj.fig_handler.reset;
            end
            
            obj.fig_handler.name = 'Buffers GUI';
            
            obj.fig_handler.width = 8;
            obj.fig_handler.height= 10;
            obj.fig_handler.bottom = 4;
            obj.makeControlPanel;
            obj.makeDisplayPanel;
            obj.makeIndicatorPanel;
            obj.makeClosePanel;
            
            obj.updateGUI;
            
        end
        
        function update(obj,~,~)
            
            obj.updateGUI;
            
        end
        
        function updateGUI(obj)
            
            if ~obj.check
                return;
            end
            
            if obj.debug_bit, disp('updating buffer GUI'); end
            
            obj.updateControlPanel;
            obj.updateDisplayPanel;
            obj.updateIndicatorPanel;
            
        end
        
        function c = check(obj)
            
            c = ~isempty(obj) && ~isempty(obj.fig_handler) && obj.fig_handler.check && ~isempty(obj.panel_controls) && isvalid(obj.panel_controls);
            
        end
        
        function set.buffers(obj, val)
            
            obj.buffers = val;
           
            if obj.check
                obj.makeIndicatorPanel;
            end
            
            obj.updateGUI;
            
        end
        
    end
            
    properties % basic controls
        
        panel_controls;
        
        button_async_write;
        button_mex_write;
        button_pars_write;
        button_file_type;
        input_deflate;
        input_product_type;
        input_debug_bit;
        
    end
    
    methods % control panel create/update/callbacks
        
        function makeControlPanel(obj)
            
            obj.panel_controls = uipanel('Title', 'controls', 'Position', [0 0.3 0.9 0.7]);
            
            Ntot = 7;
            N = Ntot;
            
            N = N-1;
            obj.button_async_write = uicontrol(obj.panel_controls, 'Style', 'toggle', 'String', 'async',...
                'Units', 'Normalized', 'Position', [0 N/Ntot 1 1/Ntot], 'FontSize', 10,...
                'Callback', @obj.callback_async_write);
            
            N = N-1;
            obj.button_mex_write = uicontrol(obj.panel_controls, 'Style', 'toggle', 'String', 'mex',...
                'Units', 'Normalized', 'Position', [0 N/Ntot 1 1/Ntot], 'FontSize', 10,...
                'Callback', @obj.callback_mex_write);
            
            N = N-1;
            obj.button_pars_write = uicontrol(obj.panel_controls, 'Style', 'toggle', 'String', 'pars',...
                'Units', 'Normalized', 'Position', [0 N/Ntot 1 1/Ntot], 'FontSize', 10,...
                'Callback', @obj.callback_pars_write);
            
            
            N = N-1;
            obj.button_file_type = uicontrol(obj.panel_controls, 'Style', 'pushbutton', 'String', 'deflate',...
                'Units', 'Normalized', 'Position', [0 N/Ntot 1 1/Ntot], 'FontSize', 10,...
                'Callback', @obj.callback_file_type);
                        
            N = N-1;
            obj.input_deflate = uicontrol(obj.panel_controls, 'Style', 'edit', 'String', 'deflate',...
                'Units', 'Normalized', 'Position', [0 N/Ntot 1 1/Ntot], 'FontSize', 10,...
                'Callback', @obj.callback_deflate);
            
            N = N-1;
            obj.input_product_type = uicontrol(obj.panel_controls, 'Style', 'edit', 'String', 'numImages',...
                'Units', 'Normalized', 'Position', [0 N/Ntot 1 1/Ntot], 'FontSize', 10,...
                'Callback', @obj.callback_product_type);
            
            N = N-1;
            obj.input_debug_bit = uicontrol(obj.panel_controls, 'Style', 'edit', 'String', 'debug_bit',...
                'Units', 'Normalized', 'Position', [0 N/Ntot 1 1/Ntot], 'FontSize', 10,...
                'Callback', @obj.callback_debug_bit);
            
        end
        
        function updateControlPanel(obj)
            
            if obj.buffers(1).use_async
                obj.button_async_write.Value = 0;
                obj.button_async_write.String = 'async';
            else                
                obj.button_async_write.Value = 1;
                obj.button_async_write.String = 'no async';
            end
                        
            if obj.buffers(1).use_mex
                obj.button_mex_write.Value = 0;
                obj.button_mex_write.String = 'mex';
            else                
                obj.button_mex_write.Value = 1;
                obj.button_mex_write.String = 'no mex';
            end  
            
            if obj.buffers(1).use_write_header
                obj.button_pars_write.Value = 0;
                obj.button_pars_write.String = 'write pars';
            else                
                obj.button_pars_write.Value = 1;
                obj.button_pars_write.String = 'no write pars';
            end
            
            obj.button_file_type.String = ['type: ' obj.buffers(1).file_type];
            
            obj.input_deflate.String = ['deflate= ' num2str(obj.buffers(1).use_deflate)];
            obj.input_product_type.String = obj.buffers(1).product_type;
            obj.input_debug_bit.String = ['debug_bit= ' num2str(obj.buffers(1).debug_bit)];
            
            
        end
        
        function callback_async_write(obj, hndl, ~)
            
            val = ~hndl.Value;
            
            if obj.debug_bit>1, disp(['callback: async write: ' num2str(val)]); end
            
            for ii = 1:length(obj.buffers)
                obj.buffers(ii).use_async = val;
            end
            
            obj.updateGUI;
            
        end
        
        function callback_mex_write(obj, hndl, ~)
                       
            val = ~hndl.Value;
            
            if obj.debug_bit>1, disp(['callback: mex write: ' num2str(val)]); end
            
            for ii = 1:length(obj.buffers)
                obj.buffers(ii).use_mex = val;
            end
            
            obj.updateGUI;
            
        end
        
        function callback_pars_write(obj, hndl, ~)
                       
            val = ~hndl.Value;
            
            if obj.debug_bit>1, disp(['callback: pars write: ' num2str(val)]); end
            
            for ii = 1:length(obj.buffers)
                obj.buffers(ii).use_write_pars = val;
            end
            
            obj.updateGUI;
            
        end
        
        function callback_file_type(obj, ~, ~)
            
            if obj.debug_bit>1, disp('callback: file_type'); end
            
            for ii = 1:length(obj.buffers)
%                 obj.buffers(ii).cycle_file_type;
            end
            
            obj.updateGUI;
                                    
        end
        
        function callback_deflate(obj, hndl, ~)
                        
            val = util.text.extract_numbers(hndl.String);
            val = val{1};
            
            if isempty(val)
                val = 0; % default value!
            end
            
            if val<0 || val>9
                error(['cannot set a deflate value of ' val ' choose a number between 0 and 9...']);
            end
            
            if obj.debug_bit>1, disp(['callback: deflate: ' num2str(val)]); end
            
            for ii = 1:length(obj.buffers)
                obj.buffers(ii).use_deflate = val;
            end
            
            obj.updateGUI;
            
        end
        
        function callback_product_type(obj, hndl, ~)
                                    
            val = hndl.String;
            
            if obj.debug_bit>1, disp(['callback: product type: ' val]); end
            
            for ii = 1:length(obj.buffers)
                obj.buffers(ii).product_type = val;
            end
            
            obj.updateGUI;
            
        end
        
        function callback_debug_bit(obj, hndl, ~)
            
            val = util.text.extract_numbers(hndl.String);
            val = val{1};
            
            if isempty(val)
                val = 0; % default value!
            end
                        
            if obj.debug_bit>1, disp(['callback: debug_bit: ' num2str(val)]); end
            
            for ii = 1:length(obj.buffers)
                obj.buffers(ii).debug_bit = val;
            end
            
            obj.updateGUI;
            
        end
        
    end
    
    properties % display statistics
        
        panel_display;
        
        button_save_time;
        button_image_size;
        
    end
    
    methods % display panel create/update
       
        function makeDisplayPanel(obj)
            
            obj.panel_display = uipanel('Title', 'display', 'Position', [0 0.1 0.9 0.2]);
            
            Ntot = 2;
            N = 0;
            
            obj.button_save_time = uicontrol(obj.panel_display, 'Style', 'pushbutton', ...
                'Units', 'Normalized', 'Position', [N/Ntot, 0, 1/Ntot, 1], 'FontSize', 10, ...
                'Callback', @obj.update);
            
            N = N+1;
            
            obj.button_image_size = uicontrol(obj.panel_display, 'Style', 'pushbutton', ...
                'Units', 'Normalized', 'Position', [N/Ntot, 0, 1/Ntot, 1], 'FontSize', 10,...
                'Callback', @obj.update);
                        
            
        end
        
        function updateDisplayPanel(obj)
            
            N = length(obj.buffers);
            N_used = 0;
            t_measured = zeros(N,1);
            
            for ii = 1:N
                                
                t= obj.buffers(ii).getMeanSaveTime;
                
                if t>0
                    t_measured(ii) = t;
                    N_used = N_used + 1;
                end
                
            end
            
            obj.button_save_time.String = ['write time= ' num2str(sum(t_measured)./N_used,3) 's'];
            
            obj.button_image_size.String = ['im.size= ' util.text.print_vec(obj.buffers(1).im_size, 'x')];
            
        end
        
    end
   
    properties % writing indicators panel
       
        indicator_panel;
%         indicator_lights@matlab.ui.control.UIControl; % a vector of buttons
        indicator_lights;
    
    end
    
    methods % indicator make/update/callback
        
        function makeIndicatorPanel(obj)
           
            obj.indicator_panel = uipanel('Title', 'busy', 'Position', [0.9 0 0.1 1]);
            
            N = length(obj.buffers);
                        
            for ii = N:-1:1
            
                obj.indicator_lights(ii) = uicontrol(obj.indicator_panel, 'Style', 'Pushbutton', ...
                    'Units', 'Normalized', 'Position', [0, (ii-1)/N, 1, 1/N], 'Callback', @obj.update);
            
            end
            
            obj.updateIndicatorPanel;
                
            
        end
        
        function updateIndicatorPanel(obj)
           
            for ii = 1:length(obj.buffers)
               
                if obj.buffers(ii).isWriting
                    set(obj.indicator_lights(ii), 'BackgroundColor', 'Red');
                else
                    set(obj.indicator_lights(ii), 'BackgroundColor', 'Green');
                end
            end
            
        end
        
    end
    
    properties % close panel
        
        panel_close;
        button_close;
        
    end
    
    methods % close panel make/callback 
    
        function makeClosePanel(obj)
            
            obj.panel_close = uipanel('Position', [0 0 0.9 0.1]);
            
            obj.button_close = uicontrol(obj.panel_close, 'Style', 'pushbutton', 'String', 'CLOSE',...
                'Units', 'Normalized', 'Position', [0 0 1 1], 'FontSize', 14, ...
                'Callback', @obj.callback_close);
            
        end
        
        function callback_close(obj, ~, ~)
            
            if obj.debug_bit>1, disp('callback: close!'); end
            
            delete(obj.fig_handler.fig);
            
        end
        
    end
        
end