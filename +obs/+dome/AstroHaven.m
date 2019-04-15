% AstroHaven class  
% Package: +obs/+dome/
% Description: A class for controlling the AstroHaven dome.
%              The class opens a serial port to the instrument
% Tested : Matlab R2018a
%     By : Guy Nir                    Dec 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: %
% Reliable: ?
%--------------------------------------------------------------------------

classdef AstroHaven < handle
    
    properties % objects
        
        hndl; % serial port object
        
        acc1@obs.sens.Accelerometer; % accelerometer for shutter1
        acc2@obs.sens.Accelerometer; % accelerometer for shutter2
        
        log@util.sys.Logger;
        
    end
    
    properties % definitions, switches, conrols
        
        status = 0;
        id = 'dome';
        
        port_name = 'COM13'; % change this later
        
        use_accelerometers = 0;
        
        reply = '';
        
        debug_bit = 1;
        
    end
    
    properties (Dependent = true)
        
        is_closed;
        shutter1;
        shutter1_deg;
        shutter2;
        shutter2_deg;
        
    end
    
    properties (Hidden = true)
        
        fail_reply_counter = 0;
        max_fail_reply = 3;
        fail_connect_counter = 0;
        max_fail_connect = 3;
        
        acc_name = 'HC-06';
        acc_id1 = '';
        acc_id2 = '';
        
        % vector of acceleration "g" for calibration of accelerometers
        cal_g_acc_vec1;
        cal_g_acc_vec2;
        
        % open and close timing for opening angle without accelerometers
        open_time1 = 0;
        close_time1 = 0;
        open_time2 = 0;
        close_time2 = 0;
        
        % calibrate close/open
        cal_open_time1 = 25;
        cal_close_time1 = 25;
        cal_open_time2 = 25;
        cal_close_time2 = 25;
        
        timeout = 10; % how many seconds to wait before returning the control 
        loop_res = 10; % how many times to call "open" or "close" when in a loop
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = AstroHaven(varargin)
            
            if isempty(varargin)
            
                if obj.debug_bit, fprintf('Boltwood default constructor v%4.2f\n', obj.version); end
                
                obj.log = util.sys.Logger('AstroHaven_dome');
                
            end
            
            
            obj.connect;
            
        end
        
        function delete(obj)
            
            obj.disconnect;
            
        end
        
    end
    
    methods % resetters
        
        function reset(obj)
            
            obj.fail_reply_counter = 0;
            obj.fail_reply_counter = 0;
            
        end
        
    end
    
    methods % getters
        
        function val = get.is_closed(obj)
            
            val = strcmp(obj.reply, '0');
            
        end
        
        function val = get.shutter1(obj)
            
            % 0: all closed, 1: shutter 2 open, 2: shutter 1 open 3: both open
            
            if any(strcmp(obj.reply, {'2', '3'}))
                val = 'open';
            elseif any(strcmp(obj.reply, {'0', '1'}))
                val = 'closed';
            else
                val = 'error';
            end
            
        end
        
        function val = get.shutter1_deg(obj)
            
            if strcmp(obj.shutter1, 'closed')
                val = 90;
            else
                val = obj.calcAngle(1);
            end
            
        end
        
        function val = get.shutter2(obj)
            
            % 0: all closed, 1: shutter 2 open, 2: shutter 1 open 3: both open
            
            if any(strcmp(obj.reply, {'1', '3'}))
                val = 'open';
            elseif any(strcmp(obj.reply, {'0', '2'}))
                val = 'closed';
            else
                val = 'error';
            end
            
        end
        
        function val = get.shutter2_deg(obj)
                                    
            if strcmp(obj.shutter2, 'closed')
                val = 90;
            else
                val = obj.calcAngle(2);
            end
            
        end
        
    end
    
    methods % commands
                
        function emergencyClose(obj)
            
            obj.log.input('Emergency close!');
            
            try 
                % obj.counter = 0;
                obj.send('C');
                obj.closeBoth(100);
                obj.update;
            catch ME
                obj.log.error(ME.getReport);
                warning(ME.getReport);
            end
            
        end
        
        function openBoth(obj, number)
            
            if nargin<2 || isempty(number)
                number = 1;
            end
            
            obj.log.input(['Open both. N= ' num2str(number)]);
            
            try
            
                obj.update;
                
                for ii = 1:number

                    t = tic;

                    % obj.counter = 0;
                    
                    obj.send('a');
                    
                    pause(0.1); 
                    
                    obj.send('b');

                    pause(0.1);

                    obj.open_time1 = obj.open_time1 + toc(t);
                    obj.open_time2 = obj.open_time2 + toc(t);

                end
                
                obj.update;

            catch ME
                obj.log.error(ME.getReport);
                warning(ME.getReport);
            end
            
        end
        
        function closeBoth(obj, number)
            
            obj.log.input(['Close both. N= ' num2str(number)]);
            
            if nargin<2 || isempty(number)
                number = 1;
            end
            
            try
                
                obj.update;
                
                for ii = 1:number

                    t = tic;
                    
                    % obj.counter = 0;
                    obj.send('A');
                    obj.send('B');

                    pause(0.2);

                    obj.close_time1 = obj.close_time1 + toc(t);
                    obj.close_time2 = obj.close_time2 + toc(t);

                end
                
                obj.update;
                
            catch ME
                obj.log.error(ME.getReport);
                warning(ME.getReport);
            end
            
        end
                
        function open1(obj, number)
            
            if nargin<2 || isempty(number)
                number = 1;
            end
            
            obj.log.input(['Open shutter 1. N= ' num2str(number)]);
            
            try

                command = 'a';

                obj.update;
                
                for ii = 1:number

                    t = tic;

                    % obj.counter = 0; % this prevents an infinite loop of send calling itself...
                    obj.send(command);
                    
                    pause(0.2);
                    
                    obj.open_time1 = obj.open_time1 + toc(t);

                end

                obj.update;
                    
            catch ME
                obj.log.error(ME.getReport);
                warning(ME.getReport);
            end
            
        end
        
        function open2(obj, number)
            
            if nargin<2 || isempty(number)
                number = 1;
            end
            
            obj.log.input(['Open shutter 2. N= ' num2str(number)]);

            try
                
                command = 'b';

                obj.update;
                
                for ii = 1:number

                    t = tic;

                    % obj.counter = 0; % this prevents an infinite loop of send calling itself...
                    obj.send(command);

                    pause(0.2);
                    
                    obj.open_time2 = obj.open_time2 + toc(t);

                end
                
                obj.update;
                
            catch ME
                obj.log.error(ME.getReport);
                warning(ME.getReport);
            end
            
        end
        
        function close1(obj, number)
            
            if nargin<2 || isempty(number)
                number = 1;
            end
            
            obj.log.input(['Close shutter 1. N= ' num2str(number)]);
            
            try
                
                obj.update;
                
                command = 'A';

                for ii = 1:number

                    t = tic;

                    % obj.counter = 0; % this prevents an infinite loop of send calling itself...
                    obj.send(command);
    %                 fprintf(obj.hndl, command);

                    pause(0.2);

                    obj.close_time1 = obj.close_time1 + toc(t);

                end
                
                obj.update;

            catch ME
                obj.log.error(ME.getReport);
                warning(ME.getReport);
            end
            
        end
        
        function close2(obj, number)
                      
            if nargin<2 || isempty(number)
                number = 1;
            end
            
            obj.log.input(['Close shutter 2. N= ' num2str(number)]);

            try

                command = 'B';

                obj.update;
                
                for ii = 1:number

                    t = tic;

                    % obj.counter = 0; % this prevents an infinite loop of send calling itself...
                    obj.send(command);

                    pause(0.2);

                    obj.close_time2 = obj.close_time2 + toc(t);

                end
                
                obj.update;

            catch ME
                obj.log.error(ME.getReport);
                warning(ME.getReport);
            end
            
        end
        
        function open1Full(obj)
            
        end
        
        function close1Full(obj)
            
        end
        
        function open2Full(obj)
            
        end
        
        function close2Full(obj)
            
        end
        
    end
    
    methods % internal functions (to be made hidden/private)
        
        function connect(obj)
            
            str = sprintf('connecting to dome! attempt %d', obj.fail_connect_counter);
            
            if obj.debug_bit>1, disp(str); end
            
            obj.log.input(str);
            
            if obj.fail_connect_counter>=3
                disp(['Cannot connect to dome, giving up after ' num2str(obj.fail_connect_counter) ' tries']);
                obj.hndl = [];
                err_msg = sprintf('Cannot attempt additional connections to dome (limited to %d attempts)', obj.max_fail_connect);
                obj.log.error(err_msg);
                error(err_msg);
            end
            
            try
                
                if ~isempty(obj.hndl) && isvalid(obj.hndl)
                    obj.disconnect;
                end
                
                pause(0.1);
                
                obj.hndl = serial(obj.port_name);

                obj.hndl.Terminator = '';

                for ii = 1:3
                
                    if strcmp(obj.hndl.Status, 'open'), break; end
                    
                    try 
                        fopen(obj.hndl);
                    catch ME
                        
                        if strcmp(ME.identifier, 'MATLAB:serial:fopen:opfailed')
                            if obj.debug_bit>1, fprintf('failed to open serial port, attempt %d\n', ii); end 
                        else 
                            rethrow(ME);
                        end
                        
                    end
                    
                    pause(0.1);
                    
                end % loop over trying to connect
                   
                if ~strcmp(obj.hndl.Status, 'open')
                    
                    if obj.debug_bit>1, disp('Giving up on opening serial port...'); end
                    
                    obj.fail_connect_counter = obj.fail_connect_counter + 1;
                    if ~isempty(obj.hndl)
                        delete(obj.hndl);
                        obj.hndl = [];
                    end
                    
                else
                    
                    if obj.debug_bit>1, disp('Successful reconnect!'); end
                    
                    obj.fail_connect_counter = 0;
                    
                end
                    
                obj.update;

                if obj.use_accelerometers
                    obj.connectAccelerometers;
                end

            catch ME
                obj.log.error(ME.getReport);
                warning(ME.getReport);
            end
            
        end
        
        function disconnect(obj)
            
            if obj.debug_bit>1, disp('disconnecting from dome!'); end
            
            obj.log.input('Disconnecting from dome');
            
            try 
                if ~isempty(obj.hndl)
                    fclose(obj.hndl);
                    delete(obj.hndl);
                end
                obj.hndl = [];
            catch ME
                obj.log.error(ME.getReport);
                warning(ME.getReport);
            end
            
        end
        
        function connectAccelerometers(obj)
            
            obj.log.input('Connecting to accelerometers.');
            
            try
                
                obj.acc1 = obs.sens.Accelerometer(obj.acc_name, obj.acc_id1);
                
            catch ME
                
                obj.log.error(ME.getReport);
                obj.acc1 = obs.sens.Accelerometer.empty;
                warning(ME.getWarning);
                
            end
            
            try
                
                obj.acc2 = obs.sens.Accelerometer(obj.acc_name, obj.acc_id2);
                
            catch ME
                
                obj.log.error(ME.getReport);
                obj.acc2 = obs.sens.Accelerometer.empty;
                warning(ME.getWarning);
                
            end
            
        end

        function send(obj, command)
            
            if obj.debug_bit>1, disp(['sending command ' command]); end
            
            for ii = 1:3
                
                try 
                    if ~isempty(obj.hndl)
                        fprintf(obj.hndl, command);
                    end
                    return;
                catch ME

                    if strcmp(ME.identifier, 'MATLAB:serial:fprintf:opfailed')
                        if obj.debug_bit>1, disp('Failed to send command, trying again...'); end
                    else
                        rethrow(ME);
                    end

                end
                
                pause(0.05);

            end
            
            if obj.debug_bit>1, disp('Gave up on sending this command... trying to reconnect!'); end
            
            pause(0.1);
            obj.connect;
            
        end
        
        function update(obj)
            
            if obj.debug_bit>1, fprintf('updating data...\n'); end
            
            if isempty(obj.hndl)
                obj.status = 0;
                return;
            else
                flushinput(obj.hndl);
            end
            
            obj.getReply;
            
            if isempty(obj.reply)
                obj.status = 0;
            else
                obj.status = 1;
            end

        end
        
        function getReply(obj, ~, ~)
            
            num=0;
            
            try
                [obj.reply, num] = fread(obj.hndl, 1);
                obj.reply = char(obj.reply);
            catch ME
                disp('--read failed--');
                obj.reply = '';
                
            end
            
            if num==0 || isempty(obj.reply)
                if obj.debug_bit>1, disp(['Couldnt read reply from serial port... reconnecting (fail_reply_counter= ' num2str(obj.fail_reply_counter) ')']); end
                warning('off', 'MATLAB:serial:fread:unsuccessfulRead');

                if obj.fail_reply_counter>obj.max_fail_reply % too many attempts to reconnect
                    error('Could not get reply from dome...');
                else
                    obj.fail_reply_counter = obj.fail_reply_counter + 1;
                end
                
                obj.connect;
                
            end
            
            if strcmp(obj.reply, '0') % reset all time estimates when closed
                obj.open_time1 = 0;
                obj.close_time1 = 0;
                obj.open_time2 = 0;
                obj.close_time2 = 0; 
            end
            
            obj.fail_reply_counter = 0; % after successful reply we can reset the counter
            
            if obj.debug_bit>1, fprintf('                reply: %s\n', obj.reply); end
            
        end
        
        function val = calcAngle(obj, shutter)

            if nargin<2 || isempty(shutter)
                shutter = 1;
            end
            
            if shutter==1
                
                if obj.use_accelerometers && ~isempty(obj.acc1) && ~isempty(obj.cal_g_acc_vec1) % use accelerometer
                    
                    val = asind(sum(obj.acc1.acc_vec.*obj.cal_g_acc_vec1)./sqrt(sum(obj.acc1.acc_vec.^2).*sum(obj.acc1.acc_vec.^2))); % use dot procuct to calculate the angle, use asind because dome angle is 90-theta
                    
                elseif ~isempty(obj.cal_open_time1) && ~isempty(obj.cal_open_time1)
                    
                    fractional_open = obj.open_time1./obj.cal_open_time1;
                    fractional_close = obj.close_time1./obj.cal_close_time1;
                    
                    current_fraction = 1-fractional_open+fractional_close;
                    current_fraction = max(current_fraction, 0);
                    current_fraction = min(current_fraction, 1);
                    
                    val = current_fraction*90;
                    
                else
                
                    val = [];
                    
                end
                
            elseif shutter==2
                
                if obj.use_accelerometers && ~isempty(obj.acc2) && ~isempty(obj.cal_g_acc_vec2) % use accelerometer
                    
                    val = asind(sum(obj.acc2.acc_vec.*obj.cal_g_acc_vec2)./sqrt(sum(obj.acc2.acc_vec.^2).*sum(obj.acc2.acc_vec.^2))); % use dot procuct to calculate the angle, use asind because dome angle is 90-theta
                    
                elseif ~isempty(obj.cal_open_time2) && ~isempty(obj.cal_open_time2)
                    
                    fractional_open = obj.open_time2./obj.cal_open_time2;
                    fractional_close = obj.close_time2./obj.cal_close_time2;
                    
                    current_fraction = 1-fractional_open+fractional_close;
                    current_fraction = max(current_fraction, 0);
                    current_fraction = min(current_fraction, 1);
                    
                    val = current_fraction*90;
                    
                else
                
                    val = [];
                    
                end
                
            end
            
        end
        
    end
    
end
