classdef ScopeAssistant < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        hndl; % serial object
        telescope; % link back to telescope object (e.g., for sending stop signal)
        data@util.vec.CircularBuffer;
        
    end
    
    properties % inputs/outputs
        
        dataCol = {'X', 'Y', 'Z', 'dist'};
        
        distance; % in cm (from ultrasonic sensor)
        
        period; % what interval we use for getting feedback from the sensor (zero for no automatic updates)
        
        % acc_vec = acc_vec_data/gain + bias
        gain = [256 256 256]; % x,y,z 
        bias = [0 0 0]; % x,y,z
        down = [0 0 -1]; % define what is the down direction
        
        calibration_data;
        
    end
    
    properties % switches/controls
        
        port_name = 'COM20'; % used only for USB (not used now)
        bluetooth_name = 'HC-06'; 
        bluetooth_id = ''; 
        
        use_check_alt = 1;
        alt_limit = 10;
        
        default_period = 0.1;
        
        status = 0;
        time;
        jd;
        reply = '';
        acc_vec_raw;
        
        log_message_sent = 0;
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        acc_vec;
        angle;
        ALT;
        
    end
    
    properties(Hidden=true)
       
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = ScopeAssistant(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.sens.ScopeAssistant')
                if obj.debug_bit, fprintf('ScopeAssistant copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('ScopeAssistant constructor v%4.2f\n', obj.version); end
                
                obj.connect(varargin{:});
                obj.data = util.vec.CircularBuffer(100);
                obj.reset;
                
            end
            
        end
        
        function connect(obj, varargin)
            
%             obj.connectUSB(varargin{:});

            obj.disconnect;

            disp('connecting arduino bluetooth');

            obj.connectBluetooth(varargin{:});
            
            pause(0.1);
            
            obj.setupTimer(obj.default_period);
            
        end
        
        function connectBluetooth(obj, name, id)
            
            if nargin>1 && ~isempty(name)
                obj.bluetoot_name = name;
            end
            
            if nargin>2 && ~isempty(id)
                obj.bluetooth_id = id;
            end
            
            % first, make sure to close existing connections...
            try
                fclose(obj.hndl);
                delete(obj.hndl);
            end
            
            if isempty(obj.bluetooth_name)
                error('Must supply a name for a bluetooth device (e.g., HC-06)');
            end
            
            % must be paired to the bluetooth device! 
            obj.hndl = Bluetooth(obj.bluetooth_name, 1); % second argument is channel==1
            
            if isempty(obj.bluetooth_id) || isnumeric(obj.bluetooth_id)
                
                in = instrhwinfo('bluetooth', obj.bluetooth_name);
                
                if isnumeric(obj.bluetooth_id) && obj.bluetooth_id>0
                    idx = obj.bluetooth_id;
                else
                    idx = 1;
                end
                
                if idx>length(in)
                    error('Cannot open device number %d, there are only %d devices named "%s".', idx, length(in), obj.bluetooth_name);
                end
                
                if isempty(in(idx).RemoteID)
                    error('Device not found. Make sure to pair the device!');
                end
                
                if obj.debug_bit
                    fprintf('Found a bluetooth device with ID: %s\n', obj.bluetooth_id);
                end
                
                obj.bluetooth_id = in(idx).RemoteID(9:end);
                
            end
            
            obj.hndl.RemoteID = obj.bluetooth_id; 
            
            fopen(obj.hndl);
            
            pause(0.1);
            
            obj.update;
            
        end
        
        function connectUSB(obj, port_name)
            
            if nargin>1 && ~isempty(port_name)
                obj.port_name = port_name;
            end
            
            % first, make sure to close existing connections...
            try
                fclose(obj.hndl);
                delete(obj.hndl);
            end
            
            if isempty(obj.port_name)
                error('Must supply a port name serial device (e.g., COM1)');
            end
            
            % must be paired to the bluetooth device! 
            obj.hndl = serial(obj.port_name);
            
            fopen(obj.hndl);
            
        end
        
        function disconnect(obj)
            
            if ~isempty(obj.hndl)
                fclose(obj.hndl);
                delete(obj.hndl);
                obj.hndl = [];
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset_calibration(obj)
            
            obj.calibration_data = [];
            obj.gain = [];
            obj.bias = [];
            
        end
        
        function reset(obj)
            
            obj.data.reset;
            
        end
        
    end
    
    methods % getters
        
        function val = is_connected(obj)
            val = ~isempty(obj.hndl) && isvalid(obj.hndl) && strcmp(obj.hndl.Status, 'open');
        end
        
        function val = get.acc_vec(obj)
            
            if isempty(obj.gain) || isempty(obj.bias)
                val = obj.acc_vec_raw;
            else
                val = (obj.acc_vec_raw - obj.bias)./obj.gain;
            end
            
        end
        
        function val = acc_vec_mean(obj)
            
            val = mean(obj.data.data(:,2:4));
            
        end
        
        function val = acc_vec_std(obj)
            
            val = std(obj.data.data(:,2:4));
            
        end
        
        function val = distance_mean(obj)
            
            val = mean(obj.data.data(:,5));
            
        end
        
        function val = get.angle(obj)
            
            if isempty(obj.acc_vec) || isempty(obj.down)
                val = [];
            else
                val = acosd(sum(obj.down.*obj.acc_vec)./sqrt(sum(obj.down.^2).*sum(obj.acc_vec.^2)));
            end
            
        end
        
        function val = get.ALT(obj)
            
            val = 90 - obj.angle;
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % commands
        
        function update(obj)
            
            obj.status = 0;
            
            if ~obj.is_connected
                error('Device is closed, use fopen or connect function');
            end
            
            obj.hndl.BytesAvailableFcn = @obj.read_data;
            
            obj.use_check_alt = 1;
            
            fprintf(obj.hndl, 'measure;');
            
        end
        
        function setupTimer(obj, period)
            
            if ~strcmp(obj.hndl.Status, 'open')
                error('Device is closed, use fopen or connect function');
            end
            
            fprintf(obj.hndl, 'timer, %f;', period);
            
            obj.update;
            
        end
        
        function read_data(obj, ~, ~)
            
            obj.reply = strip(fgetl(obj.hndl)); % text reply

%             reply = str2double(regexp(obj.reply,'-?\d*','Match'));
            numeric_reply = str2double(split(obj.reply, ','))';
            
            if isnan(numeric_reply), disp(obj.reply); end
            
            if length(numeric_reply)<5, return; end
            
            obj.acc_vec_raw = numeric_reply(1:3);
            obj.distance = numeric_reply(4);
            obj.period = numeric_reply(5);
            obj.time = datetime('now', 'timezone', 'UTC');
            obj.jd = juliandate(obj.time);
            obj.status = 1;
            
            obj.data.input([obj.jd obj.acc_vec obj.distance]);
            
            if obj.debug_bit>1
                disp(['ALT= ' num2str(obj.ALT)]);
            end
            
            if ~isempty(obj.telescope) && obj.telescope.use_accelerometer && obj.use_check_alt && obj.ALT<obj.alt_limit
                
                try 
                    
                    obj.telescope.stop;

                    if obj.log_message_sent==0
                        obj.telescope.log.input(['Arduino stopped telescope at angle ALT: ' num2str(obj.ALT) ' degrees...']);
                        fprintf('%s: Arduino sending stop signal to telescope!\n', obj.telescope.log.report(1:8));
                        obj.log_message_sent = 1;
                    end

                catch ME
                    warning(ME.getReport);
                end
                
            end
            
            if obj.ALT>30
                obj.log_message_sent = 0;
            end
            
        end
        
        function measureCalibrationdata(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('plotting', 0);
            input.input_var('time', 10, 'duration');
            input.input_var('N', 1000, 'number', 'measurements');
            input.scan_vars(varargin{:});
            
            obj.reset_calibration;
            
            for ii = 1:input.N
                
                obj.update;
                
                obj.calibration_data(ii,:) = obj.acc_vec;
                
                pause(input.time./input.N);
                
                if input.plotting
                    obj.showCalibration;
                    drawnow;
                end
                
            end
            
        end
        
        function sum_err_square = runCalibration(obj, varargin)
            
            data = obj.calibration_data;
            
            B_initial = mean(data);
            G_initial = mean(sqrt(sum((data-B_initial).^2,2)))*[1 1 1]; % assume equal gain at first

            b_initial = [B_initial, G_initial];

            func = @(b) obs.sens.Accelerometer.calcErrorsGainBias(data, b(1:3), b(4:6), B_initial, G_initial);

            b_final = fminsearch(func, b_initial);

            obj.bias = b_final(1:3);
            obj.gain = b_final(4:6);

            data_cal = (data-obj.bias)./obj.gain;
            
            S = sum( sum(data_cal.^2,2) - 1, 1);

            if obj.debug_bit
                
                fprintf('Calibration complete. GAIN= %f %f %f | BIAS= %f %f %f\n', ...
                    obj.gain(1), obj.gain(2), obj.gain(3), obj.bias(1), obj.bias(2), obj.bias(3));
                
                fprintf('---> Total summed errors= %f\n', S);
                
            end
            
            if nargout>0
                sum_err_square = S;
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function showCalibration(obj, varargin)
            
            X = obj.calibration_data(:,1);
            Y = obj.calibration_data(:,2);
            Z = obj.calibration_data(:,3);
            
            scatter3(X,Y,Z,'.b');
            
            if ~isempty(obj.gain) && ~isempty(obj.bias) % also show the calibrated results
                hold on;
                scatter3(X./obj.gain(1)+obj.bias(1), Y./obj.gain(2)+obj.bias(2), Z./obj.gain(3)+obj.bias(3), '.r');
                hold off;
            end
            
        end
        
        function liveView(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('N', 1000, 'number', 'measurements');
            input.input_var('interval', 0.1);
            input.input_var('ax', [], 'axes', 'axis');
            input.scan_vars(varargin{:});
            
            obj.reset;
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            for ii = 1:input.N
                
                obj.update; 
                
                pause(input.interval);
                
                color = linspace(0,1,obj.data.N);
                color = [color', zeros(size(color,2),2)];
                
                if ~isvalid(input.ax)
                    return;
                end
                
                if ~isempty(obj.data.data)
                    scatter3(input.ax, obj.data.data(:,2), obj.data.data(:,3), obj.data.data(:,4), 3, color(1:size(obj.data.data,1),:));
                end
                
                input.ax.XLim = [-1 1]*1.2;
                input.ax.YLim = [-1 1]*1.2;
                input.ax.ZLim = [-1 1]*1.2;
                
                title(input.ax, ['Acc (x,y,z)= (' num2str(obj.acc_vec,2) ') | d= ' num2str(obj.distance) 'cm'], 'FontSize', 14);
                drawnow;
                
            end
            
        end
        
        function livePlot(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('N', 1000, 'number', 'measurements');
            input.input_var('interval', 0.1);
            input.input_var('ax', [], 'axes', 'axis');
            input.scan_vars(varargin{:});
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            start_jd = juliandate(datetime('now', 'timezone', 'utc'));
            
            for ii = 1:input.N
                
                if obj.period==0
                    obj.update; 
                end
                
                pause(input.interval);
                
                if ~isvalid(input.ax)
                    return;
                end
                
                d = vertcat(obj.data.data_ordered);
                
                if ~isempty(obj.data.data)
                    plot(input.ax, (d(:,1)-start_jd)*24*3600, d(:,2:4));
                end
                
                drawnow;
                
            end
            
        end
        
    end
   
    methods (Static=true)
       
        function sum_err_square = calcErrorsGainBias(data, Biases, Gains, B_initial, G_initial) % calculate the calibration fit

            new_data = (data-Biases)./Gains;

            errors = sqrt(sum((new_data).^2,2)) - 1;

            bayes_term = sum(cosh(Biases./B_initial).*cosh(Gains./G_initial)); % this keeps the parameters from getting too big

            errors = errors.*bayes_term.*1e-5;

            sum_err_square = sum(errors.^2,1);

        end 
        
    end
    
end

