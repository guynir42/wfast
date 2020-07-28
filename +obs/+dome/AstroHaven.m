classdef AstroHaven < handle
% Control the AstroHaven clamshell dome. 
% The interface to the dome is via serial port. 
% This is fairly convinient as long as the COM port does not change. 
% Since this is a clamshell there are very few controls available. 
% You can open each shutter using openEast(number) or openWest(number)
% or openBoth(number). Same for closeXXX(number). 
% The number is for how many keystrokes are sent to the serial controller, 
% and correspond to degrees (very roughly). 
% The default number for each side or both is saved in the object properties, 
% so these commands can be called without arguments (e.g., from the GUI). 
%
% There are limit sensors telling the controller when each shutter is fully
% closed or fully open. If it is in the middle, there is no sensor to tell
% how closed it is, so we guess this using the time it takes to open and 
% close, and then interpolate from that. This is not very accurate so take
% the shutter position angle with a grain of salt. 
%
% The dome uses the same "brake_bit" strategy as other devices, so setting
% that to 1 stops the dome, e.g., using the GUI stop button. 
% Also, clicking Ctrl+C in the command line breaks out of the loop and stops
% the dome from opening/closing. 
%
% In the future we will implement a timed open/close loop to slowly track 
% the telescope from East to West. 
% Also we plan to install accelerometers to know the real angle of each 
% shutter, instead of trying to guess it. 


    properties(Transient=true)
        
        gui;
        
    end
    
    properties % objects
        
        hndl; % serial port object
        
        owner@obs.Manager; % link back to the top level manager object
        
        acc_W@obs.sens.Accelerometer; % accelerometer for west shutter
        acc_E@obs.sens.Accelerometer; % accelerometer for east shutter
        
        cam_pc@obs.comm.PcSync; % communicate with the camera-PC
        
        log@util.sys.Logger; % write reports to text file
        
        timer; % for dome tracking
        
    end
    
    properties % definitions, switches, conrols
        
        status = 0;
        id = 'dome';
        
        port_name = 'COM13'; % change this later
        
        use_accelerometers = 0; % for future installation
        
        use_tracking = 0; % when on, the dome timer will slowly open the west shutter and close the east shutter 
        track_rate = 25; % how many "steps" to open west / close east to take every 30 minutes
        
        max_fail_reply = 3; % how many sent messages to try before failing the command
        max_fail_connect = 3; % how many reconnect attempts before failing to connect
        
        reply = ''; % reply (single letter) from the serial controller
        
        number_west = 50; % default value to move West shutter
        number_east = 50; % default value to move East shutter
        number_both = 50; % default value to move both shutters
        
        brake_bit = 1;
        debug_bit = 1;
        
    end
    
    properties (Dependent = true)
        
        is_closed; % if both shutter report being closed
        shutter_west; % West shutter situation (string)
        shutter_west_deg; % West shutter position angle (degrees)
        shutter_east; % East shutter situation (string)
        shutter_east_deg; % East shutter position angle (degrees)
        
    end
    
    properties (Hidden = true)
        
        dome_radius = 275; % cm
        dome_angle = 30; % degrees (relative to East-West direction, clockwise) 
        tel_offset = 50; % cm (offset of telescope from center of dome)
        angle_spare = 20; % additional altitude to spare from the calculated shutter position
        
        % accelerometer names (for future)
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
        
        default_number_west;
        default_number_east;
        default_number_both;
        
        version = 1.04;
        
    end
    
    methods % constructor
        
        function obj = AstroHaven(varargin)
            
            if isempty(varargin)
            
                if obj.debug_bit>1, fprintf('Boltwood default constructor v%4.2f\n', obj.version); end
                
                obj.log = util.sys.Logger('AstroHaven_dome', obj);
                
                util.oop.save_defaults(obj);
                
            end
            
            obj.connect;
            obj.log.heartbeat(600);
            
        end
        
        function delete(obj) % make sure serial port is deleted when clearing this object
            
            obj.disconnect;
            
        end
        
        function connect(obj) % connect to serial port
            
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
                obj.log.error(ME);
                rethrow(ME);
            end

        end
        
        function disconnect(obj) % close, delete and clear the serial object "hndl"
            
            if obj.debug_bit>1, disp('disconnecting from dome!'); end
            
            obj.log.input('Disconnecting from dome');
            
            try 
                if ~isempty(obj.hndl)
                    fclose(obj.hndl);
                    delete(obj.hndl);
                end
                obj.hndl = [];
            catch ME
                obj.log.error(ME);
                rethrow(ME);
            end
            
        end
        
        function connectAccelerometers(obj) % for future installation
            
            obj.log.input('Connecting to accelerometers.');
            
            try
                
                obj.acc_W = obs.sens.Accelerometer(obj.acc_name, obj.acc_id_west);
                
            catch ME
                
                obj.log.error(ME);
                obj.acc_E = obs.sens.Accelerometer.empty;
                warning(ME.getWarning);
                
            end
            
            try
                
                obj.acc_E = obs.sens.Accelerometer(obj.acc_name, obj.acc_id_east);
                
            catch ME
                
                obj.log.error(ME);
                obj.acc_E = obs.sens.Accelerometer.empty;
                warning(ME.getWarning);
                
            end
            
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
%             elseif strcmp(obj.shutter_west, 'open full')
%                 val = 0; 
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
        
        function val = get.use_tracking(obj)
            
            if obj.use_tracking && ~isempty(obj.timer) && isa(obj.timer, 'timer') && isvalid(obj.timer) && strcmp(obj.timer.Running, 'on')
                val = 1;
            else
                val = 0;
            end
            
        end
        
    end
    
    methods % setters
        
        function set.use_tracking(obj, val)
            
            if val
            
                if isempty(obj.timer) || ~isa(obj.timer, 'timer') || ~isvalid(obj.timer) || strcmp(obj.timer.Running, 'off')
                    obj.setup_timer;
                end
            
            end
            
            obj.use_tracking = val;
                
        end
        
    end
    
    methods % commands to move or stop shutters
        
        function adjustDomeMaxOpen(obj, side) % open the dome as far as you can
            
            if obj.is_closed
                error('cannot adjust the dome when it is closed!'); 
%                 return;
            end
            
            if util.text.cs(side, 'East')
                obj.openEastFull;
                obj.closeWestFull;
                obj.openWest(100); % make sure the top shutter is off
            elseif util.text.cs(side, 'West')
                obj.openBothFull;
            else
                error('Unknown option "%s" to "side" parameter. Use "East" or "West". ', side); 
            end
            
        end
        
        function adjustDomeNew(obj, side, Az, Alt)
        % Usage: obj.adjustDome(side, Az, Alt)
        % Automatically choose the correct position for the dome shutters, 
        % based on the prefered side and telescope Az/Alt. 
        % Inputs: -side: 'East' or 'West', depending on declination and 
        %                proximity to meridian. 
        %          -Az: telescope azimuth (degrees). 
        %          -Alt: telescope altitude (degrees). 
        % 
        % Will open/close the shutters and start the tracking if needed.
        %
        % NOTE: will not do anything if dome is closed. 
        %
        % To adjust the dome position use the different dome parameters:
        % dome_radius, dome_angle, tel_offset, angle_spare
        % Also we depend on the shutter angle estimates, 
        
            if nargin<4
                help('obs.dome.AstroHaven.adjustDome'); 
                return;
            end
            
            if obj.is_closed
                error('cannot adjust the dome when it is closed!'); 
%                 return;
            end
            
            if isempty(HA_deg)
                error('Must supply a valid hour angle (HA)!'); 
            end
            
            if isempty(Dec_deg) || Dec_deg<-35
                error('Unable to view targets with %d declinations!', Dec_dec);
            end
            
            if util.text.cs(side, 'East') % east side is still controlled by heuristics
                
                obj.openEastFull;
                obj.closeWestFull;
                
                % estimate for how much we need to open West shutter when 
                % reaching meridian (at Dec+30 need 100, at Dec+90 need 200, 
                % at -30 don't need to open at all). 
                west_num = 5/3*(Dec_deg+30); 
                if west_num>0
                    obj.openWest(west_num);
                end
                
                obj.use_tracking = 0; % no need for tracking on East side
                
            elseif util.text.cs(side, 'West') % new geometry driven formula
                
                
                obj.openWestFull; % get a reference position to move from 
                obj.closeEastFull; % check if we need to open this at all...
                
                % estimate for how much we need to open East shutter when
                % starting observations on meridian (at Dec+30 need 20, 
                % at Dec 0 need 80, at Dec -30 need almost full open). 
                east_num = -2*(Dec_deg-40); 
                if east_num>0
                    obj.openEast(east_num); % east side still uses heuristic formula
                end
                
                angle_west = Az - 90 + obj.dome_angle; % the angle of telescope from center of West shutter. 
                
                angle_alt = obj.offset_angle(Alt); % get the altitude angle the dome should have to match the telescope altitude angle
                angle_alt = angle_alt - obj.angle_spare; % get some spare distance
                
                obj.use_tracking = 1; % I don't see any reason not to track on the West side
                
                % do we need to set the tracking rate in case someone changed it?
                
            else
                error('Unknown option "%s" to "side" parameter. Use "East" or "West". ', side); 
            end
            
        end
        
        function adjustDome(obj, side, HA_deg, Dec_deg)
        % Usage: obj.adjustDome(side, HA_deg, Dec_deg)
        % Automatically choose the correct position for the dome shutters, 
        % based on the prefered side and target coordinates. 
        % Inputs: -side: 'East' or 'West', depending on declination and 
        %                proximity to meridian. 
        %          -HA_deg: target hour angle in degrees. 
        %          -Dec_deg: target declination in degrees. 
        % 
        % Will open/close the shutters and start the tracking if needed.
        %
        % NOTE: will not do anything if dome is closed. 
        %
        % To get the correct parameters you need to know the dome geometry,
        % relative to the mount+OTA, and also know the conversion btw the 
        % open/close number and the angle it opens (assuming it is linear, 
        % which it almost certainly isn't). 
        % In practice, we just measured this empirically and used
        % parameters that gave reasonable results at all angles. 
        
            if nargin<4
                help('obs.dome.AstroHaven.adjustDomeOld'); 
                return;
            end
            
            if obj.is_closed
                error('cannot adjust the dome when it is closed!'); 
%                 return;
            end
            
            if isempty(HA_deg)
                error('Must supply a valid hour angle (HA)!'); 
            end
            
            if isempty(Dec_deg) || Dec_deg<-35
                error('Unable to view targets with %d declinations!', Dec_dec);
            end
            
            if util.text.cs(side, 'East')
                
                obj.openEastFull;
                obj.closeWestFull;
                
                % estimate for how much we need to open West shutter when 
                % reaching meridian (at Dec+30 need 100, at Dec+90 need 200, 
                % at -30 don't need to open at all). 
                west_num = 5/3*(Dec_deg+30); 
                if west_num>0
                    obj.openWest(west_num);
                end
                
                obj.use_tracking = 0; % no need for tracking on East side
                
            elseif util.text.cs(side, 'West')
                
                obj.openWestFull; % get a reference position to move from 
                obj.closeEastFull; % check if we need to open this at all...
                
                % estimate for how much we need to open East shutter when
                % starting observations on meridian (at Dec+30 need 20, 
                % at Dec 0 need 80, at Dec -30 need almost full open). 
                east_num = -2*(Dec_deg-40); 
                if east_num>0
                    obj.openEast(east_num); 
                end
                
                % need to figure out how much to close the West shutter, 
                % assuming the dome will also track (slowly open this side
                % and close the East side during the observations). 
                
                angle_to_horizon = 90 - 20 - HA_deg; % degrees
                west_num = 160 - 1.2*(Dec_deg+30) + angle_to_horizon*2;
                if west_num>0
                    obj.closeWest(west_num);
                end
                
                obj.use_tracking = 1; % I don't see any reason not to track on the West side
                
                % do we need to set the tracking rate in case someone changed it?
                
            else
                error('Unknown option "%s" to "side" parameter. Use "East" or "West". ', side); 
            end
            
        end
        
        function emergencyClose(obj) % do everything you can to close dome (including sending the uncancelable 'C' command) 
            
            obj.log.input('Emergency close!');
            
            try 
                % obj.counter = 0;
                obj.use_tracking = 0;
                obj.send('C');
                obj.closeBoth(100);
                obj.update;
            catch ME
                obj.log.error(ME);
                warning(ME.getReport);
            end
            
        end
        
        function openBoth(obj, number) % open both shutters a certain amount
            
            if nargin<2 || isempty(number)
                if ~isempty(obj.number_both) && obj.number_both>0
                    number = obj.number_both;
                else
                    number = 1;
                end
            end
            
            number = ceil(number); % avoid fractional numbers
            
            obj.log.input(['Open both. N= ' num2str(number)]);
            
            try
            
                obj.update;
                
                obj.use_tracking = 0;
                
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
                obj.log.error(ME);
                rethrow(ME);
            end
            
        end
        
        function closeBoth(obj, number) % close both shutters a certain amount
            
            if nargin<2 || isempty(number)
                if ~isempty(obj.number_both) && obj.number_both>0
                    number = obj.number_both;
                else
                    number = 1;
                end
            end
            
            number = ceil(number); % avoid fractional numbers
            
            obj.log.input(['Close both. N= ' num2str(number)]);
            
            try
            
                obj.update;
                
                obj.use_tracking = 0;
                
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
                obj.log.error(ME);
                rethrow(ME);
            end
            
        end
            
        function openBothFull(obj) % send command to both shutters until they are fully open
            
            obj.openBoth(1000);
            
        end
        
        function closeBothFull(obj) % send command to both shutters until they are fully open
            
            obj.closeBoth(1000);
            
        end
            
        function openWest(obj, number) % open West shutter a certain amount
            
            if nargin<2 || isempty(number)
                if ~isempty(obj.number_west) && obj.number_west>0
                    number = obj.number_west;
                else
                    number = 1;
                end
            end
            
            number = ceil(number); % avoid fractional numbers
            
            obj.log.input(['Open shutter West. N= ' num2str(number)]);
            
            try
            
                obj.update;
                
                obj.use_tracking = 0;
                
                t = tic;
                
                reply = obj.command('a', number); % can add a "max_duration" argument to allow the loop to stop after so long
                
                if isempty(reply)
                    % what to do if a command didn't succeed? should this be an error?
                else
                    obj.open_time_west = obj.open_time_west + toc(t);
                end
                
                obj.update;

            catch ME
                obj.log.error(ME);
                rethrow(ME);
            end
            
        end
        
        function closeWest(obj, number) % close West shutter a certain amount
            
            if nargin<2 || isempty(number)
                if ~isempty(obj.number_west) && obj.number_west>0
                    number = obj.number_west;
                else
                    number = 1;
                end
            end
            
            number = ceil(number); % avoid fractional numbers
            
            obj.log.input(['Close shutter West. N= ' num2str(number)]);
            
            try
            
                obj.update;
                
                obj.use_tracking = 0;
                
                t = tic;
                
                reply = obj.command('A', number); % can add a "max_duration" argument to allow the loop to stop after so long
                
                if isempty(reply)
                    % what to do if a command didn't succeed? should this be an error?
                else
                    obj.close_time_west = obj.close_time_west + toc(t);
                end
                obj.update;

            catch ME
                obj.log.error(ME);
                rethrow(ME);
            end
            
        end
        
        function openWestFull(obj) % send command to West shutter until it is fully open
            
            obj.openWest(1000);
            
        end
        
        function closeWestFull(obj) % send command to West shutter until it is fully closed
            
            obj.closeWest(1000);
            
        end
        
        function openEast(obj, number) % open East shutter a certain amount
            
            if nargin<2 || isempty(number)
                if ~isempty(obj.number_east) && obj.number_east>0
                    number = obj.number_east;
                else
                    number = 1;
                end
            end
            
            number = ceil(number); % avoid fractional numbers
            
            obj.log.input(['Open shutter East. N= ' num2str(number)]);
            
            try
            
                obj.update;
                
                obj.use_tracking = 0;
                
                t = tic;
                
                reply = obj.command('b', number); % can add a "max_duration" argument to allow the loop to stop after so long
                
                if isempty(reply)
                    % what to do if a command didn't succeed? should this be an error?
                else
                    obj.open_time_east = obj.open_time_east + toc(t);
                end
                
                obj.update;

            catch ME
                obj.log.error(ME);
                rethrow(ME);
            end
            
        end
        
        function closeEast(obj, number) % close East shutter a certain amount
             
            if nargin<2 || isempty(number)
                if ~isempty(obj.number_east) && obj.number_east>0
                    number = obj.number_east;
                else
                    number = 1;
                end
            end
            
            number = ceil(number); % avoid fractional numbers
            
            obj.log.input(['Close shutter East. N= ' num2str(number)]);
            
            try
            
                obj.update;
                
                obj.use_tracking = 0;
                
                t = tic;
                
                reply = obj.command('B', number); % can add a "max_duration" argument to allow the loop to stop after so long
                
                if isempty(reply)
                    % what to do if a command didn't succeed? should this be an error?
                else
                    obj.close_time_east = obj.close_time_east + toc(t);
                end
                
                obj.update;

            catch ME
                obj.log.error(ME);
                rethrow(ME);
            end
            
        end
        
        function openEastFull(obj) % send command to East shutter until it is fully open
            
            obj.openEast(1000);
            
        end
        
        function closeEastFull(obj) % send command to East shutter until it is fully closed
            
            obj.closeEast(1000);
            
        end
        
        function stop(obj) % stop the motion of the shutters
            
            obj.brake_bit = 1;
            
            obj.use_tracking = 0;
            
            if ~isempty(obj.gui)
                obj.gui.update;
            end
            
        end
        
        function setup_timer(obj, ~, ~)
            
            if ~isempty(obj.timer) && isa(obj.timer, 'timer') && isvalid(obj.timer)
                
                if strcmp(obj.timer.Running, 'on')
                    stop(obj.timer);
                end
                
                delete(obj.timer);
                obj.timer = [];
                
            end
            
            delete(timerfind('name', 'dome-timer'));
            
            obj.timer = timer('BusyMode', 'queue', 'ExecutionMode', 'fixedRate', 'Name', 'dome-timer', ...
                'Period', 30*60, 'StartDelay', 30*60, 'TimerFcn', @obj.callback_timer, 'ErrorFcn', @obj.setup_timer);
            
            start(obj.timer);
            
        end
        
        function callback_timer(obj, ~, ~)
            
            if obj.use_tracking
                
                obj.update;

                obj.log.input('Dome is tracking (open West / close East)');

                try 
                   
                    % first move the West shutter down 
                    date_str = util.text.time2str(datetime('now', 'TimeZone', 'UTC')); 
                    if obj.debug_bit, fprintf('%s: Opening West shutter by %d steps\n', date_str, obj.track_rate); end
                    
                    t = tic;
                
                    reply = obj.command('a', obj.track_rate); 

                    if isempty(reply)
                        % what to do if a command didn't succeed? should this be an error?
                    else
                        obj.open_time_west = obj.open_time_west + toc(t);
                    end

                    % now close the East shutter
                    date_str = util.text.time2str(datetime('now', 'TimeZone', 'UTC'));
                    if obj.debug_bit, fprintf('%s: Closing East shutter by %d steps\n', date_str, obj.track_rate); end
                    
                    t = tic;

                    reply = obj.command('B', obj.track_rate); 

                    if isempty(reply)
                        % what to do if a command didn't succeed? should this be an error?
                    else
                        obj.close_time_east = obj.close_time_east + toc(t);
                    end

                    obj.use_tracking = 1;
                    
                    obj.update;

                    if ~isempty(obj.gui)
                        obj.gui.update;
                    end
                    
                catch ME
                    obj.use_tracking = 0;
                    obj.log.error(ME); 
                    rethrow(ME); 
                end
                
            end
            
        end
        
    end
    
    methods(Hidden=true) % internal functions
        
        function val = offset_angle(obj, base_angle, offset, radius) % calculate the dome angle, given the telescope angle, telescope offset, and dome radius
            
            if nargin<3 || isempty(offset)
                offset = obj.tel_offset;
            end
            
            if nargin<4 || isempty(radius)
                radius = obj.dome_radius; 
            end
            
            if base_angle>=90
                val = 90;
            elseif base_angle<0
                val = 0;
            else
                func = @(b) abs(tand(base_angle)-(radius*sind(b))/(offset+radius*cosd(b))); 
                val = fminsearch(func, base_angle);    
            end
            
        end
        
        function reply = command(obj, command_vector, number) % generic interface to move shutters ("command_vector" is a string sent "number" of times to serial port)
            
            if nargin<3 || isempty(number)
                number = 1;
            end
            
            number = ceil(number); % avoid fractional numbers
            
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
                obj.log.error(ME);
                rethrow(ME);
            end
            
        end
        
        function reply = send(obj, command) % send a single string to the serial port
            
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
        
        function update(obj) % communicate with hardware to make sure it is still connected
            
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
        
        function reply = getReply(obj, ~, ~) % read the serial port
            
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
        
        function val = calcAngle(obj, shutter) % calculate the shutter angle by interpolating the time it takes to open/close

            import util.text.cs;
            
            if nargin<2 || isempty(shutter)
                shutter = 'West';
            end
            
            if cs(shutter, 'West')
                
                if obj.use_accelerometers && ~isempty(obj.acc_W) && ~isempty(obj.cal_g_acc_vec_west) % use accelerometer
                    
                    val = asind(sum(obj.acc_W.acc_vec.*obj.cal_g_acc_vec_west)./sqrt(sum(obj.acc_W.acc_vec.^2).*sum(obj.acc_W.acc_vec.^2))); % use dot procuct to calculate the angle, use asind because dome angle is 90-theta
                    
                elseif ~isempty(obj.open_time_west) && ~isempty(obj.cal_open_time_west) && ...
                        ~isempty(obj.close_time_west) && ~isempty(obj.cal_close_time_west)
                    
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
