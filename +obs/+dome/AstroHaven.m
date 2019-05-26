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
    
    properties(Transient=true)
        
        gui;
        
    end
    
    properties % objects
        
        hndl; % serial port object
        
        acc1@obs.sens.Accelerometer; % accelerometer for shutter1
        acc2@obs.sens.Accelerometer; % accelerometer for shutter2
        
        log@util.sys.Logger;
        
    end
    
    properties % definitions, switches, conrols
        
        status = 0;
        id = 'dome';
        
        port_name = 'COM14'; % change this later
        
        use_accelerometers = 0;
        
        max_fail_reply = 3;
        max_fail_connect = 3;
        
        reply = '';
        
        number1 = 20; 
        number2 = 20;
        number_both = 20;
        
        brake_bit = 1;
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
        
        default_number1;
        default_number2;
        default_number_both;
        
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = AstroHaven(varargin)
            
            if isempty(varargin)
            
                if obj.debug_bit, fprintf('Boltwood default constructor v%4.2f\n', obj.version); end
                
                obj.log = util.sys.Logger('AstroHaven_dome', obj);
                
                util.oop.save_defaults(obj);
                
            end
            
            obj.connect;
            obj.log.heartbeat(600);
            
        end
        
        function delete(obj)
            
            obj.disconnect;
            
        end
        
    end
    
    methods % resetters
        
        function reset(obj)
            
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
                if ~isempty(obj.number_both) && obj.number_both>0
                    number = obj.number_both;
                else
                    number = 1;
                end
            end
            
            obj.log.input(['Open both. N= ' num2str(number)]);
            
            try
            
                obj.update;
                
                t = tic;
                
                reply = obj.command('ab', number); % can add a "max_duration" argument to allow the loop to stop after so long
                
                if isempty(reply)
                    % what to do if a command didn't succeed? should this be an error?
                else
                    obj.open_time1 = obj.open_time1 + toc(t);
                    obj.open_time2 = obj.open_time2 + toc(t);
                end
                
                obj.update;

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function closeBoth(obj, number)
            
            if nargin<2 || isempty(number)
                if ~isempty(obj.number_both) && obj.number_both>0
                    number = obj.number_both;
                else
                    number = 1;
                end
            end
            
            obj.log.input(['Close both. N= ' num2str(number)]);
            
            try
            
                obj.update;
                
                t = tic;
                
                reply = obj.command('AB', number); % can add a "max_duration" argument to allow the loop to stop after so long
                
                if isempty(reply)
                    % what to do if a command didn't succeed? should this be an error?
                else
                    obj.close_time1 = obj.close_time1 + toc(t);
                    obj.close_time2 = obj.close_time2 + toc(t);
                end
                
                obj.update;

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
                
        function open1(obj, number)
            
            if nargin<2 || isempty(number)
                if ~isempty(obj.number1) && obj.number1>0
                    number = obj.number1;
                else
                    number = 1;
                end
            end
            
            obj.log.input(['Open shutter1. N= ' num2str(number)]);
            
            try
            
                obj.update;
                
                t = tic;
                
                reply = obj.command('a', number); % can add a "max_duration" argument to allow the loop to stop after so long
                
                if isempty(reply)
                    % what to do if a command didn't succeed? should this be an error?
                else
                    obj.open_time1 = obj.open_time1 + toc(t);
                end
                
                obj.update;

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function open2(obj, number)
            
            if nargin<2 || isempty(number)
                if ~isempty(obj.number2) && obj.number2>0
                    number = obj.number2;
                else
                    number = 1;
                end
            end
            
            obj.log.input(['Open shutter2. N= ' num2str(number)]);
            
            try
            
                obj.update;
                
                t = tic;
                
                reply = obj.command('b', number); % can add a "max_duration" argument to allow the loop to stop after so long
                
                if isempty(reply)
                    % what to do if a command didn't succeed? should this be an error?
                else
                    obj.open_time2 = obj.open_time2 + toc(t);
                end
                
                obj.update;

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function close1(obj, number)
            
            if nargin<2 || isempty(number)
                if ~isempty(obj.number1) && obj.number1>0
                    number = obj.number1;
                else
                    number = 1;
                end
            end
            
            obj.log.input(['Close shutter1. N= ' num2str(number)]);
            
            try
            
                obj.update;
                
                t = tic;
                
                reply = obj.command('A', number); % can add a "max_duration" argument to allow the loop to stop after so long
                
                if isempty(reply)
                    % what to do if a command didn't succeed? should this be an error?
                else
                    obj.close_time1 = obj.close_time1 + toc(t);
                end
                obj.update;

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function close2(obj, number)
             
            if nargin<2 || isempty(number)
                if ~isempty(obj.number2) && obj.number2>0
                    number = obj.number2;
                else
                    number = 1;
                end
            end
            
            obj.log.input(['Close shutter2. N= ' num2str(number)]);
            
            try
            
                obj.update;
                
                t = tic;
                
                reply = obj.command('B', number); % can add a "max_duration" argument to allow the loop to stop after so long
                
                if isempty(reply)
                    % what to do if a command didn't succeed? should this be an error?
                else
                    obj.close_time2 = obj.close_time2 + toc(t);
                end
                
                obj.update;

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function open1Full(obj)
            
            obj.open1(1000);
            
        end
        
        function close1Full(obj)
            
            obj.close1(1000);
            
        end
        
        function open2Full(obj)
            
            obj.open2(1000);
            
        end
        
        function close2Full(obj)
            
            obj.close2(1000);
            
        end
        
        function openBothFull(obj)
            
            obj.openBoth(1000);
            
        end
        
        function closeBothFull(obj)
            
            obj.closeBoth(1000);
            
        end
        
    end
    
    methods % internal functions (to be made hidden/private)
        
        function connect(obj)
            
            try
            
                for ii = 1:obj.max_fail_connect

                    str = sprintf('connecting to dome! attempt %d', ii);

                    if obj.debug_bit>1, disp(str); end

                    obj.log.input(str);

                    if ~isempty(obj.hndl) && isvalid(obj.hndl)
                        obj.disconnect;
                    end

                    pause(0.1);

                    obj.hndl = serial(obj.port_name);

                    obj.hndl.Terminator = '';

                    try 
                        fopen(obj.hndl);
                    catch ME

                        if strcmp(ME.identifier, 'MATLAB:serial:fopen:opfailed') % if this is the regular connection error, just report it and try again
                            str = sprintf('failed to open serial port, attempt %d\n', ii); 
                            if obj.debug_bit>1, disp(str); end 
                            obj.log.input(str);
                        else 
                            rethrow(ME); % if this is some other error, throw it up the line
                        end

                    end

                    if strcmp(obj.hndl.Status, 'open'), break; end % if we succeed, no need to continue with the loop

                end

                if strcmp(obj.hndl.Status, 'open') % if the loop ended with success

                    if obj.debug_bit>1, disp('Successful reconnect!'); end

                    obj.update;

                    if obj.use_accelerometers
                        obj.connectAccelerometers;
                    end

                else % failed to connect after so many tries

                    if obj.debug_bit>1, disp('Giving up on opening serial port...'); end
                    obj.log.input('Failed to connect to serial port :(');

                    if ~isempty(obj.hndl)
                        obj.disconnect
                    end

                end
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
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
                rethrow(ME);
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

        function reply = command(obj, command_vector, number)
            
            if nargin<3 || isempty(number)
                number = 1;
            end
            
            if isempty(obj.hndl) || ~isvalid(obj.hndl) || ~strcmp(obj.hndl.Status, 'open')
                obj.connect;
            end
            
            obj.log.input(['Sending commands ' command_vector ' for ' num2str(number) ' times.']);
            
            try 

                on_cleanup = onCleanup(@obj.stop);
                obj.brake_bit = 0;
                
                for ii = 1:number
                    
                    list_replies = '';
                    
                    for jj = 1:length(command_vector)

                        for kk = 1:obj.max_fail_reply
                            
                            if obj.brake_bit
                                return;
                            end

                            if ~isempty(obj.gui)
                                obj.gui.update;
                            end
                            
                            reply = obj.send(command_vector(jj)); 
                            
                            pause(0.05); 
                            
                            if ~isempty(obj.reply)
                                list_replies = [list_replies reply]; % keep track of what the dome returned
                                break;
                            end
                        
                        end % for kk (attempts to send)
                        
                        if isempty(reply)
                            error('Failed to send command %s after %d attempts!', command_vector(jj), kk);
                        end
                        
                    end % for jj (command list)
                    
                    % NOTE: when e.g., closing shutter 1, the command is A.
                    %       The reply when it is successfully closed is X, 
                    %       which is 23 bigger (in ASCII) than A. 
                    %       Thus we can identify when a command can be removed.
                    idx = list_replies-23==command_vector; % which command has gotten confirmation that dome is closed/open
                    
                    command_vector(idx) = [];
                    
                    if isempty(command_vector)
                        return;
                    end
                    
                end % for ii (number)

            catch ME                
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function reply = send(obj, command)
            
            reply = ''; % empty reply means failed to send
            
            try 
                if ~isempty(obj.hndl)
                    
                    flushinput(obj.hndl);
                    fprintf(obj.hndl, command);
                    reply = obj.getReply;
                    if obj.debug_bit>4, fprintf('sent: %s | reply: %s\n', command, reply); end
                end
                
            catch ME

                if strcmp(ME.identifier, 'MATLAB:serial:fprintf:opfailed')
                    if obj.debug_bit>1, disp('Failed to send command...'); end
                    reply = '';
                elseif strcmp(ME.identifier, 'MATLAB:serial:flushinput:opfailed')
                    if obj.debug_bit>1, disp('Failed to flush input...'); end
                    reply = '';
                else
                    rethrow(ME); % any other error is reported up the line
                end

            end

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
            
            if ~isempty(obj.gui)
                obj.gui.update;
            end

        end
        
        function reply = getReply(obj, ~, ~)
            
            warning('off', 'MATLAB:serial:fread:unsuccessfulRead');
            
            num=0;
            reply = '';
            
            try
                [reply, num] = fread(obj.hndl, 1);
                reply = char(reply);
            catch ME
                if strcmp(ME.identifier, 'MATLAB:serial:fread:opfailed')
                    if obj.debug_bit>1, disp('Failed to read command...'); end
                    obj.reply = '';
                else
                    rethrow(ME); % any other error is reported up the line
                end    
            end
            
            if num==0 || isempty(reply)
                obj.reply = '';
                if obj.debug_bit>1, disp('Failed to read command...'); end
                return;
            else
                obj.reply = reply;
            end
            
            if strcmp(obj.reply, '0') % reset all time estimates when closed
                obj.open_time1 = 0;
                obj.close_time1 = 0;
                obj.open_time2 = 0;
                obj.close_time2 = 0; 
            end
            
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
        
        function stop(obj)
            
            obj.brake_bit = 1;
            
            if ~isempty(obj.gui)
                obj.gui.update;
            end
            
        end
        
    end
   
    methods % plotting / GUI
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = obs.dome.gui.AstroHavenGUI(obj);
            end
            
            obj.gui.make;
        end
        
    end
    
end
