classdef ScopeAssistant < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        hndl; % serial object
        telescope; % link back to telescope object (e.g., for sending stop signal)
        data@util.vec.CircularBuffer;
        reco@obs.comm.Reconnect;
        
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
        
        port_name = 'COM25'; % used only for USB (not used now)
%         bluetooth_name = 'HC-06'; 
        bluetooth_name = 'btspp://0021130386D9';
        bluetooth_id = ''; 
        
        use_check_alt = 1;
        alt_limit = 10;
        
        default_period = 0; % replaced 0.1, trying to work without internal timer
        timeout = 5;
        
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
       
        version = 1.02;
        
    end
    
    methods % constructor
        
        function obj = ScopeAssistant(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.sens.ScopeAssistant')
                if obj.debug_bit>1, fprintf('ScopeAssistant copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('ScopeAssistant constructor v%4.2f\n', obj.version); end
                
                obj.data = util.vec.CircularBuffer(100);
                obj.reco = obs.comm.Reconnect;
                
                try 
                    obj.connect(varargin{:});
                catch ME
                    warning(ME.getReport)
                end
                
                obj.reset;
                
            end
            
        end
        
        function connect(obj, varargin)
            
%             obj.connectUSB(varargin{:});

            if obj.reco.should

                try
                
%                     obj.disconnect; % we don't have to disconnect since we could still find this bluetooth object using instrfindall 

                    t = datetime('now', 'TimeZone', 'UTC');
                    fprintf('%s: connecting arduino bluetooth\n', t);

                    obj.connectBluetooth(varargin{:});

                    pause(0.1);

                    obj.setupTimer(obj.default_period);

                    obj.reco.inputSuccess;
                
                catch ME
                    obj.reco.inputFailure(ME.getReport);
                    rethrow(ME); 
                end
                
            end
            
        end
        
        function connectBluetooth(obj, name, id)
            
            if nargin>1 && ~isempty(name)
                obj.bluetooth_name = name;
            end
            
            if nargin>2 && ~isempty(id)
                obj.bluetooth_id = id;
            end
            
            if isempty(obj.bluetooth_name)
                error('Must supply a name for a bluetooth device (e.g., HC-06)');
            end
            
            % first try to find if there is an orphan Blluetooth object with the righ ID
            inst=instrfindall; 
            
            ind = find(contains(inst.Name, obj.bluetooth_name), 1, 'first');
            
            if isempty(ind)
                
                if ~isempty(obj.hndl)
                    fclose(obj.hndl);
                    delete(obj.hndl);
                end
                
                % must be paired to the bluetooth device! 
                obj.hndl = Bluetooth(obj.bluetooth_name, 1); % second argument is channel==1
                obj.hndl.Timeout = obj.timeout;

                try
                    fopen(obj.hndl);
                catch
                    try % try, try again
                        fopen(obj.hndl);
                    catch 
                        t = datetime('now', 'TimeZone', 'UTC');
                        fprintf('%s: Cannot open bluetooth to ScopeAssistant!\n', t); 
                        delete(obj.hndl); 
                        obj.hndl = [];
                        return; 
                    end
                end

                pause(0.1);

            else
                obj.hndl = inst(ind);
            end
            
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
            
            if ~isempty(obj) && ~isempty(obj.hndl)  
                
                if isvalid(obj.hndl)
                    fclose(obj.hndl);
                end
                
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
            
            if isempty(obj.gain) || isempty(obj.bias) || isempty(obj.acc_vec_raw)
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
        
        function ok = update(obj)
            
            if ~obj.is_connected
                obj.status = 0;
                ok = 0; 
                return;
            end
            
            obj.hndl.BytesAvailableFcn = @obj.read_data;
            
            ok = obj.send('measure'); 
            
            if ok
                obj.use_check_alt = 1;
            end
            
        end
        
        function setupTimer(obj, period)
            
            if nargin<2 || isempty(period)
                period = obj.default_period;
            end
            
            if isempty(obj.hndl) || ~isvalid(obj.hndl) || ~strcmp(obj.hndl.Status, 'open')
                error('Device is closed, use fopen or connect function');
            end
            
            fprintf(obj.hndl, 'timer, %f;', period);
            
            obj.update;
            
        end
        
        function ok = send(obj, str)
            
            if nargin<2 || isempty(str)
                return;
            end
            
            if ~strcmp(obj.hndl.Status, 'open')
                warning('Device is closed, use fopen or connect function');
                ok = 0;
                return;
            end
            
            try 
                fprintf(obj.hndl, str);
                ok = 1;
            catch ME
%                 disp('Problem writing to Arduino'); 
                
                try % try again
                    fprintf(obj.hndl, str);
                    ok = 1;
                catch ME
%                     disp('Failed second attempt to write'); 
                    ok = 0;
                    obj.status = 0;
                end
                
            end
            
            
        end
        
        function read_data(obj, ~, ~)
            
            obj.reply = strip(fgetl(obj.hndl)); % text reply

%             t = datetime('now', 'TimeZone', 'UTC'); 
%             fprintf('%s Arduino reply: %s\n', t, obj.reply); 
            
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

