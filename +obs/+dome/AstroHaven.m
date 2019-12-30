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
        
        owner@obs.Manager;
        
        acc_W@obs.sens.Accelerometer; % accelerometer for west shutter
        acc_E@obs.sens.Accelerometer; % accelerometer for east shutter
        
        sync@obs.comm.PcSync;
        
        log@util.sys.Logger;
        
    end
    
    properties % definitions, switches, conrols
        
        status = 0;
        id = 'dome';
        
        port_name = 'COM29'; % change this later
        
        use_accelerometers = 0;
        
        max_fail_reply = 3;
        max_fail_connect = 3;
        
        reply = '';
        
        number_west = 50; 
        number_east = 50;
        number_both = 50;
        
        brake_bit = 1;
        debug_bit = 1;
        
    end
    
    properties (Dependent = true)
        
        is_closed;
        shutter_west;
        shutter_west_deg;
        shutter_east;
        shutter_east_deg;
        
    end
    
    properties (Hidden = true)
        
        acc_name = 'HC-06';
        acc_id_west = '';
        acc_id_east = '';
        
        % vector of acceleration "g" for calibration of accelerometers
        cal_g_acc_vec_west;
        cal_g_acc_vec_east;
        
        % open and close timing for opening angle without accelerometers
        open_time_west = 0;
        close_time_west = 0;
        open_time_east = 0;
        close_time_east = 0;
        
        % calibrate close/open
        cal_open_time_west = 25;
        cal_close_time_west = 25;
        cal_open_time_east = 25;
        cal_close_time_east = 25;
        
        timeout = 10; % how many seconds to wait before returning the control 
        loop_res = 10; % how many times to call "open" or "close" when in a loop
        
        default_number_west;
        default_number_east;
        default_number_both;
        
        version = 1.02;
        
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
        
        function val = get.shutter_west(obj)
            
            % 0: all closed, 1: east shutter open, 2: west shutter open 3: both open
            % a: west shutter closing, A: west shutter opening
            % b: east shutter closing, B: east shutter opening
            
            if any(strcmp(obj.reply, {'2', '3'}))
                val = 'open';
            elseif any(strcmp(obj.reply, {'0', '1'}))
                val = 'closed';
            elseif any(strcmp(obj.reply, {'a', 'A'}))
                val = 'West move';
            elseif any(strcmp(obj.reply, {'b', 'B'}))
                val = 'East move';
            else
                val = 'error';
            end
            
        end
        
        function val = get.shutter_west_deg(obj)
            
            if strcmp(obj.shutter_west, 'closed')
                val = 90;
            else
                val = obj.calcAngle('West');
            end
            
        end
        
        function val = get.shutter_east(obj)
            
            % 0: all closed, 1: east shutter open, 2: west shutter open 3: both open
            % a: west shutter closing, A: west shutter opening
            % b: east shutter closing, B: east shutter opening
            
            if any(strcmp(obj.reply, {'1', '3'}))
                val = 'open';
            elseif any(strcmp(obj.reply, {'0', '2'}))
                val = 'closed';
            elseif any(strcmp(obj.reply, {'a', 'A'}))
                val = 'West move';
            elseif any(strcmp(obj.reply, {'b', 'B'}))
                val = 'East move';
            else
                val = 'error';
            end
            
        end
        
        function val = get.shutter_east_deg(obj)
                                    
            if strcmp(obj.shutter_east, 'closed')
                val = 90;
            else
                val = obj.calcAngle('East');
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
                    obj.open_time_west = obj.open_time_west + toc(t);
                    obj.open_time_east = obj.open_time_east + toc(t);
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
                    obj.close_time_west = obj.close_time_west + toc(t);
                    obj.close_time_east = obj.close_time_east + toc(t);
                end
                
                obj.update;

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
                
        function openWest(obj, number)
            
            if nargin<2 || isempty(number)
                if ~isempty(obj.number_west) && obj.number_west>0
                    number = obj.number_west;
                else
                    number = 1;
                end
            end
            
            obj.log.input(['Open shutter West. N= ' num2str(number)]);
            
            try
            
                obj.update;
                
                t = tic;
                
                reply = obj.command('a', number); % can add a "max_duration" argument to allow the loop to stop after so long
                
                if isempty(reply)
                    % what to do if a command didn't succeed? should this be an error?
                else
                    obj.open_time_west = obj.open_time_west + toc(t);
                end
                
                obj.update;

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function openEast(obj, number)
            
            if nargin<2 || isempty(number)
                if ~isempty(obj.number_east) && obj.number_east>0
                    number = obj.number_east;
                else
                    number = 1;
                end
            end
            
            obj.log.input(['Open shutter East. N= ' num2str(number)]);
            
            try
            
                obj.update;
                
                t = tic;
                
                reply = obj.command('b', number); % can add a "max_duration" argument to allow the loop to stop after so long
                
                if isempty(reply)
                    % what to do if a command didn't succeed? should this be an error?
                else
                    obj.open_time_east = obj.open_time_east + toc(t);
                end
                
                obj.update;

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function closeWest(obj, number)
            
            if nargin<2 || isempty(number)
                if ~isempty(obj.number_west) && obj.number_west>0
                    number = obj.number_west;
                else
                    number = 1;
                end
            end
            
            obj.log.input(['Close shutter West. N= ' num2str(number)]);
            
            try
            
                obj.update;
                
                t = tic;
                
                reply = obj.command('A', number); % can add a "max_duration" argument to allow the loop to stop after so long
                
                if isempty(reply)
                    % what to do if a command didn't succeed? should this be an error?
                else
                    obj.close_time_west = obj.close_time_west + toc(t);
                end
                obj.update;

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function closeEast(obj, number)
             
            if nargin<2 || isempty(number)
                if ~isempty(obj.number_east) && obj.number_east>0
                    number = obj.number_east;
                else
                    number = 1;
                end
            end
            
            obj.log.input(['Close shutter East. N= ' num2str(number)]);
            
            try
            
                obj.update;
                
                t = tic;
                
                reply = obj.command('B', number); % can add a "max_duration" argument to allow the loop to stop after so long
                
                if isempty(reply)
                    % what to do if a command didn't succeed? should this be an error?
                else
                    obj.close_time_east = obj.close_time_east + toc(t);
                end
                
                obj.update;

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function openWestFull(obj)
            
            obj.openWest(1000);
            
        end
        
        function closeWestFull(obj)
            
            obj.closeWest(1000);
            
        end
        
        function openEastFull(obj)
            
            obj.openEast(1000);
            
        end
        
        function closeEastFull(obj)
            
            obj.closeEast(1000);
            
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

                    delete(instrfind('Name', ['Serial-', obj.port_name])); % get rid of any leftover serial connections to this object's hndl
                    
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
                
                obj.send('R');
                
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
                
                obj.acc_W = obs.sens.Accelerometer(obj.acc_name, obj.acc_id_west);
                
            catch ME
                
                obj.log.error(ME.getReport);
                obj.acc_E = obs.sens.Accelerometer.empty;
                warning(ME.getWarning);
                
            end
            
            try
                
                obj.acc_E = obs.sens.Accelerometer(obj.acc_name, obj.acc_id_east);
                
            catch ME
                
                obj.log.error(ME.getReport);
                obj.acc_E = obs.sens.Accelerometer.empty;
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
                            
                            if ~isempty(obj.owner) && ~isempty(obj.owner.gui) && obj.owner.gui.check
                                obj.owner.gui.updateStopButton;
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
                obj.open_time_west = 0;
                obj.close_time_west = 0;
                obj.open_time_east = 0;
                obj.close_time_east = 0; 
            end
            
        end
        
        function val = calcAngle(obj, shutter)

            import util.text.cs;
            
            if nargin<2 || isempty(shutter)
                shutter = 'West';
            end
            
            if cs(shutter, 'West')
                
                if obj.use_accelerometers && ~isempty(obj.acc_W) && ~isempty(obj.cal_g_acc_vec_west) % use accelerometer
                    
                    val = asind(sum(obj.acc_W.acc_vec.*obj.cal_g_acc_vec_west)./sqrt(sum(obj.acc_W.acc_vec.^2).*sum(obj.acc_W.acc_vec.^2))); % use dot procuct to calculate the angle, use asind because dome angle is 90-theta
                    
                elseif ~isempty(obj.cal_open_time_west) && ~isempty(obj.cal_open_time_west)
                    
                    fractional_open = obj.open_time_west./obj.cal_open_time_west;
                    fractional_close = obj.close_time_west./obj.cal_close_time_west;
                    
                    current_fraction = 1-fractional_open+fractional_close;
                    current_fraction = max(current_fraction, 0);
                    current_fraction = min(current_fraction, 1);
                    
                    val = current_fraction*90;
                    
                else
                
                    val = [];
                    
                end
                
            elseif cs(shutter, 'East')
                
                if obj.use_accelerometers && ~isempty(obj.acc_E) && ~isempty(obj.cal_g_acc_vec_east) % use accelerometer
                    
                    val = asind(sum(obj.acc_E.acc_vec.*obj.cal_g_acc_vec_east)./sqrt(sum(obj.acc_E.acc_vec.^2).*sum(obj.acc_E.acc_vec.^2))); % use dot procuct to calculate the angle, use asind because dome angle is 90-theta
                    
                elseif ~isempty(obj.cal_open_time_east) && ~isempty(obj.cal_open_time_east)
                    
                    fractional_open = obj.open_time_east./obj.cal_open_time_east;
                    fractional_close = obj.close_time_east./obj.cal_close_time_east;
                    
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
