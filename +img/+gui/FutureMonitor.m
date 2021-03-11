classdef FutureMonitor < handle

    properties(Transient=true)
        
        fig; 
        aux_fig;
        timer; 
        
    end
    
    properties % objects
        
        owner;
        
    end
    
    properties % inputs/outputs
        
        menu;
        buttons;
        delete_buttons; 
        
    end
    
    properties % switches/controls
        
        period = 5; 
        
        font_size = 16; 
        color_default = 0.94.*[1 1 1];
        color_running = [0 0.9 0.2];
        color_error = [0.9 0.2 0.1]; 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
        
        key_status_shift = 0;
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = FutureMonitor(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'img.gui.FutureMonitor')
                if obj.debug_bit>1, fprintf('FutureMonitor copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            elseif ~isempty(varargin) && isa(varargin{1}, 'img.Analysis')
                if obj.debug_bit>1, fprintf('FutureMonitor analysis constructor v%4.2f\n', obj.version); end
                obj.owner = varargin{1};
            else
                if obj.debug_bit>1, fprintf('FutureMonitor constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % make and update GUI
        
        function show(obj)
            
            if isempty(obj.fig) || ~isvalid(obj.fig)
                obj.fig = figure('Name', 'Future monitor',...
                    'Units', 'centimeters', 'Position', [5, 5, 30, 20]); 
                
                set(obj.fig, 'WindowKeyPressFcn', @obj.callback_key_press);
                set(obj.fig, 'KeyReleaseFcn', @obj.callback_key_released);
            
            else
                obj.fig = figure(obj.fig); 
            end
            
            obj.update; 
            
            obj.setup_timer; 
            
        end
        
        function update(obj)
            
            if isempty(obj.fig) || ~isvalid(obj.fig)
                return;
            end
            
            delete(obj.fig.Children); 
            
            f = obj.owner.futures;
            d = obj.owner.futures_dir; 
            b = obj.owner.futures_batches; 
            
            N = max(10, length(f)); 
            
            for ii = 1:length(f)
                
                if ~isempty(f{ii}) && isa(f{ii}, 'parallel.Future') && isvalid(f{ii})
                
                    is_running = strcmp(f{ii}.State, 'running'); 
                    is_error = ~isempty(f{ii}.Error); 
                    
                    color = obj.color_default;
                    if is_running, color = obj.color_running; end
                    if is_error, color = obj.color_error; end

                    data.index = ii;
                    data.future = f{ii}; 
                    data.dir = d{ii};
                    data.num_batches = b{ii}; 
                    
                    obj.buttons{ii} = uicontrol(obj.fig, 'Style', 'Pushbutton', ...
                        'Units', 'Normalized', 'Position', [0 (N-ii)/N 0.8 1/N], ...
                        'String', obj.print_single(ii, f{ii}, d{ii}, b{ii}), ...
                        'Callback', @obj.callback_button, 'UserData', data, ...
                        'BackgroundColor', color, 'FontSize', obj.font_size);

                    obj.delete_buttons{ii} = uicontrol(obj.fig, 'Style', 'Pushbutton', ...
                        'Units', 'Normalized', 'Position', [0.8 (N-ii)/N 0.2 1/N], ...
                        'String', 'delete future', ...
                        'Callback', @obj.callback_delete, 'UserData', data, ...
                        'BackgroundColor', color, 'FontSize', obj.font_size);

                    if is_running, obj.delete_buttons{ii}.String = 'cancel future'; end

                end                
                
            end
            
        end
        
        function val = print_single(obj, idx, future, directory, num_batches)
            
            finish_datetime = future.FinishDateTime;
            if isempty(finish_datetime)
                finish_datetime = datetime('now', 'TimeZone', 'Local');
            end
            
            start_datetime = future.StartDateTime; 
            
            val = sprintf('% 2d) %s (%d batches) --- runtime: %s ', ...
                idx, directory, num_batches, char(finish_datetime-start_datetime)); 
            
        end
        
    end
    
    methods % callbacks from buttons
        
        function callback_key_press(obj, hndl, event)
            
            if any(contains(event.Modifier,'shift'))
                
                obj.key_status_shift = 1;
                
            end
            
        end
        
        
        
        function callback_key_released(obj, hndl, event)
            
            obj.key_status_shift = 0;
            
        end
        
        function callback_button(obj, hndl, event)
            
            f = hndl.UserData.future; % get the future object
            if isempty(f) || ~isa(f, 'parallel.Future') || ~isvalid(f)
                return;
            end
            
            is_running = strcmp(f.State, 'running'); 
            is_error = ~isempty(f.Error); 
            
            if obj.key_status_shift && ispc
                d = fullfile(getenv('DATA'), 'WFAST', hndl.UserData.dir(1:4), hndl.UserData.dir);
                if exist(d, 'dir')
                    winopen(d);
                end
                obj.key_status_shift = 0; 
            elseif is_running
                disp(f.Diary); 
            elseif is_error
                delete(obj.aux_fig); 
                obj.aux_fig = figure('Name', 'Error', 'Position', [300, 200, 1000, 600]); 
                uicontrol(obj.aux_fig, 'Style', 'Text', 'String', f.Error.getReport('extended','hyperlinks','off'), ...
                    'Units', 'Normalized', 'Position', [0.02 0.02 0.96 0.96], 'FontSize', 12, ...
                    'ButtonDownFcn', @obj.callback_remove_aux_fig); 
            else
                disp(f.Diary); 
            end
            
        end
        
        function callback_delete(obj, hndl, ~)
            
            f = hndl.UserData.future; % get the future object
            if isempty(f) || ~isa(f, 'parallel.Future') || ~isvalid(f)
                return;
            end
            
            idx = hndl.UserData.index;
            
            is_running = strcmp(f.State, 'running'); 
            is_error = ~isempty(f.Error); 
            
            if is_running
                rep = questdlg('Stop this run?', sprintf('Stop run on worker %d?', idx), 'Yes', 'No', 'No'); 
                if strcmp(rep, 'Yes')
                    cancel(f);                     
                end
            else
                rep = questdlg('Delete this future?', sprintf('Delete future number %d?', idx), 'Yes', 'No', 'No'); 
                if strcmp(rep, 'Yes')                    
                    cancel(f); 
                    obj.owner.futures{idx} = [];
                end
            end
            
            obj.update;
            
        end
        
        function callback_self_destruct(obj, hndl, ~)
            
            delete(hndl); 
            
        end
        
        function callback_remove_aux_fig(obj, ~, ~)
            
            delete(obj.aux_fig); 
            
        end
        
    end
    
    methods % timer related
        
        function callback_timer(obj, ~, ~)
            
            obj.update;
            
        end
        
        function setup_timer(obj, ~, ~)
            
            if ~isempty(obj.timer) && isa(obj.timer, 'timer') && isvalid(obj.timer)
                if strcmp(obj.timer.Running, 'on')
                    stop(obj.timer);
                    delete(obj.timer);
                    obj.timer = [];
                end
            end
            
            delete(timerfind('name', 'future_monitor-timer'));
            
            obj.timer = timer('BusyMode', 'drop', 'ExecutionMode', 'fixedRate', 'Name', 'future_monitor-timer', ...
                'Period', obj.period, 'StartDelay', obj.period, 'TimerFcn', @obj.callback_timer, 'ErrorFcn', @obj.setup_timer);
            
            start(obj.timer);
            
        end
        
        function stop_timer(obj, ~, ~)
            
            if ~isempty(obj.timer) && isa(obj.timer, 'timer') && isvalid(obj.timer)
                stop(obj.timer); 
            end
            
        end
        
    end
    
end

