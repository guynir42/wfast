classdef OutletControl < handle
    
    properties(Transient=true)
        
        gui@obs.comm.gui.OutletGUI; 
        figure_dialog; 
        button_text;
        button_confirm;
        button_cancel;
        
    end
    
    properties % switches/controls
        
        sockets = {'cam_pc', 'dome_pc', 'router', 'lights', 'balor', 'arduinos', 'dome', 'mount'}; 
        
        use_demo_mode = 0; % when true, it would not send any on/off/cycle commands to the switcher
        
        status = false(8); 
        
        debug_bit = 1; 
        
    end
        
    properties (Dependent=true)
        
        cam_pc;
        dome_pc;
        router;
        lights;
        balor;
        arduinos;
        dome;
        mount;
        
    end
    
    properties (Hidden=true) % connection details like IP, user, password
        
        username = 'admin';
        password = 'kbos';
        ip_address = '192.168.1.104:8000'; 
        
        cam_pc_;
        dome_pc_;
        router_;
        lights_;
        balor_;
        arduinos_;
        dome_;
        mount_;
        
        version = 1.00; 
        
    end
    
    methods % constructor
        
        function obj = OutletControl(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.comm.OutletControl')
                if obj.debug_bit>1, fprintf('OutletControl copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('OutletControl constructor v%4.2f\n', obj.version); end
            
                obj.update;
                
            end
            
        end
        
    end
    
    methods % getters
        
        function val = get.cam_pc(obj)
            
            val = obj.cam_pc_; 
            
        end
        
        function val = get.dome_pc(obj)
            
            val = obj.dome_pc_; 
            
        end
        
        function val = get.router(obj)
            
            val = obj.router_; 
            
        end
        
        function val = get.lights(obj)
            
            val = obj.lights_; 
            
        end
        
        function val = get.balor(obj)
            
            val = obj.balor_; 
            
        end
        
        function val = get.arduinos(obj)
            
            val = obj.arduinos_; 
            
        end
        
        function val = get.dome(obj)
            
            val = obj.dome_; 
            
        end
        
        function val = get.mount(obj)
            
            val = obj.mount_; 
            
        end
        
    end
    
    methods % setters
        
        function set_value(obj, dev, val, dialog)
            
            if nargin<4 || isempty(dialog)
                dialog = 1;
            end
            
            if util.text.parse_bool(val)
                str = 'ON'; 
            else
                str = 'OFF';
            end
            
            if dialog
                obj.popupDialog(dev, str); % create a pop-up dialog that the user must confirm
            else
                obj.send_string(dev, str); % send the command without further warning
            end
            
        end
        
        function cycle(obj, dev, dialog)
            
            if nargin<3 || isempty(dialog)
                dialog = 1;
            end
            
            if dialog
                obj.popupDialog(dev, 'cycle'); % create a pop-up dialog that the user must confirm
            else
                obj.send_string(dev, 'cycle'); % send the command without further warning
            end
            
        end
        
        function set.cam_pc(obj, val)
            
            val = util.text.parse_bool(val); 
            dialog = ~val; % only create a confirmation dialog when turning off
            obj.set_value('cam_pc', val, dialog);
            
        end
        
        function cycle_cam_pc(obj)
            
            obj.cycle('cam_pc', 1); 
            
        end
        
        function set.dome_pc(obj, val)

            val = util.text.parse_bool(val); 
            dialog = ~val; % only create a confirmation dialog when turning off
            obj.set_value('dome_pc', val, dialog); 
            
        end
        
        function cycle_dome_pc(obj)
            
            obj.cycle('dome_pc', 1); 
            
        end
        
        function set.router(obj, val)
            
            val = util.text.parse_bool(val); 
            dialog = ~val; % only create a confirmation dialog when turning off
            obj.set_value('router', val, dialog); 
            
        end
        
        function cycle_router(obj)
            
            obj.cycle('router', 1); 
            
        end
        
        function set.lights(obj, val)
            
            val = util.text.parse_bool(val); 
            dialog = 0; % never ask for confirmation to turn this on/off
            obj.set_value('lights', val, dialog); 
            
        end
        
        function cycle_lights(obj)
            
            obj.cycle('lights', 0); 
            
        end
        
        function set.balor(obj, val)
            
            val = util.text.parse_bool(val); 
            dialog = ~val; % only create a confirmation dialog when turning off
            obj.set_value('balor', val, dialog); 
            
        end
        
        function cycle_balor(obj)
            
            obj.cycle('balor', 1); 
            
        end
        
        function set.arduinos(obj, val)
            
            val = util.text.parse_bool(val); 
            dialog = ~val; % only create a confirmation dialog when turning off
            obj.set_value('arduinos', val, dialog); 
            
        end
        
        function cycle_arduinos(obj)
            
            obj.cycle('arduinos', 1); 
            
        end
        
        function set.dome(obj, val)
            
            val = util.text.parse_bool(val); 
            dialog = ~val; % only create a confirmation dialog when turning off
            obj.set_value('dome', val, dialog); 
            
        end
        
        function cycle_dome(obj)
            
            obj.cycle('dome', 1); 
            
        end
        
        function set.mount(obj, val)
            
            val = util.text.parse_bool(val); 
            dialog = ~val; % only create a confirmation dialog when turning off
            obj.set_value('mount', val, dialog); 
            
        end
        
        function cycle_mount(obj)
            
            obj.cycle('mount', 1); 
            
        end
        
    end
        
    methods % commands 
        
        function update(obj)
            
            [rc, rv] = obj.command('status', 0); % second argument says that we do not use demo-mode
            
            if rc
                error('Could not connect to power switcher... rc= %d', rc); 
            end
            
            [~, idx] = regexp(rv, '<div id="state">'); 
            
            obj.status = flip(hexToBinaryVector(rv(idx+1:idx+2)));
            
            if length(obj.status)==length(obj.sockets)
                for ii = 1:length(obj.sockets)

                    if ~isempty(obj.sockets{ii})
                        obj.([obj.sockets{ii} '_']) = obj.status(ii); 
                    end

                end
            end
            
        end
        
    end
    
    methods (Hidden=true) % internal functions that have no protections against accidental power-down of the wrong device
        
        function rc = send_string(obj, dev, str)
            
            idx = find(strcmpi(dev, obj.sockets));
            
            if isempty(idx)
                error('Coult not find device "%s" in socket list. Try "%s" or "%s" etc... ', dev, obj.sockets{1}, obj.sockets{2}); 
            end
            
            if util.text.cs(str, 'cycle')
                str = 'CCL';
            end
            
            str = upper(str); 
            
            rc = obj.command(sprintf('outlet?%d=%s', idx, str), obj.use_demo_mode);
            
            obj.update; 
            
        end
        
        function [rc, rv] = command(obj, str, demo)
            
            if nargin<3 || isempty(demo)
                demo = 1;
            end
            
            command_string = sprintf('curl -u %s:%s http://%s/%s', obj.username, obj.password, obj.ip_address, str);
            
            if demo==0
                [rc, rv] = system(command_string); 
            else
                rc = 0;
                rv = '';
                disp(command_string); 
            end
            
        end
        
    end
    
    methods % printouts and GUI
    
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = obs.comm.gui.OutletGUI(obj); 
            end
            
            obj.gui.make; 
            
        end
        
        function popupDialog(obj, device, action)
            
            import util.text.cs;
            
            if ~ismember(device, obj.sockets)
                error('Unknown device "%s". Try "%s" or "%s" etc... ', obj.sockets{1}, obj.sockets{2}); 
            end
            
            if isempty(obj.figure_dialog) || ~isa(obj.figure_dialog, 'matlab.ui.Figure') || ~isvalid(obj.figure_dialog)
                obj.figure_dialog = figure('Name', 'confirm action', 'Position', [1000 600 400 300], 'KeyPressFcn', @obj.callback_key_press);
            else
                figure(obj.figure_dialog);
            end
            
            if util.text.cs(action, 'on', 'off')
                action_string = ['turn ' action]; 
            elseif util.text.cs(action, 'cycle')
                action_string = ['power ' action]; 
            else
                error('Unknown action "%s". Use "on" or "off" or "cycle". ', action); 
            end
            
            clf(obj.figure_dialog); 
            
            obj.button_text = uicontrol('Style', 'text', 'String' , sprintf('Are you sure you want to %s the %s?', action_string, upper(device)), ...
                'Units', 'Normalized', 'Position', [0 0.5 1 0.4], 'Parent', obj.figure_dialog, 'FontSize', 24, ...
                'KeyPressFcn', @obj.callback_key_press); 
            
            obj.button_confirm = uicontrol('Style', 'pushbutton', 'string', 'CONFIRM', 'FontSize', 16, ...
                'Units', 'Normalized', 'Position', [0.1 0.1 0.4 0.3], 'ForegroundColor', 'red', ...
                'Callback', @obj.callback_send_string, 'Parent', obj.figure_dialog, ...
                'KeyPressFcn', @obj.callback_key_press); 
            
            obj.button_confirm.UserData = {device, action};
            
            obj.button_cancel = uicontrol('Style', 'pushbutton', 'string', 'Abort', 'FontSize', 22, ...
                'Units', 'Normalized', 'Position', [0.6 0.1 0.3 0.3], ...
                'Callback', @obj.callback_close_prompt, 'Parent', obj.figure_dialog, ...
                'KeyPressFcn', @obj.callback_key_press); 
            
            uicontrol(obj.button_cancel); 
            
        end
        
        function callback_send_string(obj, hndl, ~)
            
            dev = hndl.UserData{1}; 
            str = hndl.UserData{2}; 
            
            if obj.debug_bit>1, fprintf('Sending command "%s" for device "%s"\n', str, dev); end
            
            obj.send_string(dev, str); 
            
            obj.callback_close_prompt;
            
        end
        
        function callback_close_prompt(obj, ~, ~)
            
            delete(obj.figure_dialog);
            obj.figure_dialog = []; 
            
        end
        
        function callback_key_press(obj, hndl, event)
            
            if strcmp(event.Key, 'leftarrow')
                if isvalid(obj.button_confirm), uicontrol(obj.button_confirm); end
            elseif strcmp(event.Key, 'rightarrow')
                if isvalid(obj.button_cancel), uicontrol(obj.button_cancel); end
            elseif strcmp(event.Key, 'return')
                if isvalid(hndl) && isequal(hndl, obj.button_confirm)
                    obj.callback_send_string(hndl); 
                elseif isvalid(hndl) && isequal(hndl, obj.button_cancel)
                    obj.callback_close_prompt; 
                end
            end
                
        end
        
    end
    
end

