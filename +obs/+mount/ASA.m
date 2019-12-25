classdef (CaseInsensitiveProperties, TruncatedProperties) ASA < handle
% Control ASA mount. 
% This class works mainly by giving an "object" and then letting the mount
% slew to that object. It also has shortcut fields to change various 
% parameters in the mount object itself (which is in "hndl"). 
%
% To give a target use inputTarget(star_name) or inputTarget(RA,DEC). 
% Make sure the target parmeters are correctly translated to the object 
% coordinates objRA and objDEC. Also make sure the objALT is big enough. 
% To make the telescope go to target just use "slew" without arguments. 
% 
% NOTE: All parameters of the mount / target are given in degrees or in 
%       sexagesimal strings (HH:MM:SS for RA and DD:MM:SS for DEC). 
%
% NOTE: this class performs some tests on the target object before the slew
%       and during the movement, and will stop the motion in case there's 
%       trouble. In addition the arduino accelerometer can be enabled to 
%       make constant checks for altitude and stop the motion of the telescope
%       if needed (even during manual slew). Make sure use_acceleromter=1. 
%       If there is a problem communicating with arduino then set it to 0. 
%
% Additional commands will be implemented later on. 
% 
% 
% PLEASE READ THE COMMENTS ON PROPERTIES FOR MORE DETAIL!

    properties(Transient=true)
        
        gui@obs.mount.gui.ASAGUI;
        
        audio@util.sys.AudioControl;
        
        prev_objects = {}; 
        
    end
    
    properties % objects
        
        hndl; % ASCOM object
        
        owner@obs.Manager;
        
        object@head.Ephemeris; % all position and time calculations done using Ephemeris class
        
        ard@obs.sens.ScopeAssistant; % connect to accelerometer and ultrasonic sensor
        
        sync@obs.comm.PcSync; % to be given from Manager to pass data on to camera computer
        
        timer; % timer object
        
        guiding_history@util.vec.CircularBuffer;
        
        log@util.sys.Logger; % keep track of all commands and errors
        
    end
    
    properties % inputs/outputs
        
        status = 0; % if connected and responsive, set this to 1
        
    end
    
    properties % switches/controls
        
        use_guiding = 1;
        use_accelerometer = 1; % make constant checks for altitude outside of the mounts own sensors
        use_ultrasonic = 0; % make constant checks that there is nothing in front of the telescope
        
        use_motor_toggle = 0; % when slewing is done, turn motor off then on again
        
        move_rate = 1; % manual slew rate in deg/sec
        
        step_arcsec = 5; % not used yet
        
        brake_bit = 1; % when slewing, set this to 0. Interrupts may set it back to 1 to stop the slew. 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        limit_alt; % degrees above horizon where telescope is not allowed to go
        limit_flip; % degrees beyond meridian where telescope can go before needing to flip
        
        % these are parameters of the target object
        % they are kept in the Ephemeris class
        objName;
        objRA_deg; 
        objDEC_deg;
        
        objRA; % string in HH:MM:SS format
        objDEC; % string in DD:MM:SS format
        
        % these are read-only values for inspection
        objHA; % string in HH:MM:SS format
        objALT; % degrees
        objAZ; % degrees 
        
        % these are read out directly from the mount        
        % these are read-only values for inspection
        telRA_deg; 
        telDEC_deg;
        
        telRA; % string in HH:MM:SS format
        telDEC; % string in DD:MM:SS format
        
        telHA; % string in HH:MM:SS format
        telALT; % degrees
        telAZ; % degrees
        
        LST; % string in HH:MM:SS format
        
        % can we get motor status from hndl??
        tracking; % 1 when tracking is on
        
        rate; % do we need this??
        
        rate_RA;
        rate_DE;
        
        pier_side;
        
    end
    
    properties(Hidden=true)
       
        was_tracking = 0;
        default_move_rate;
        
        num_prev_objects = 10; % how many prev_objects do we keep?
        
        version = 1.03;
        
    end
    
    methods % constructor and connect commands
        
        function obj = ASA(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.mount.ASA')
                if obj.debug_bit, fprintf('ASA copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('ASA constructor v%4.2f\n', obj.version); end
            
            end
            
            obj.guiding_history = util.vec.CircularBuffer;
            
            obj.log = util.sys.Logger('ASA_mount', obj);
            
            try
            
                obj.object = head.Ephemeris;

                obj.connect;

                obj.update;

            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function delete(obj)
           
            obj.disconnect;
            
        end
        
        function connect(obj)
            
            obj.log.input('Connecting to mount');
            
            try 
                
                if ~isempty(obj.hndl)
                    obj.disconnect;
                    pause(0.1);
                end
                
                obj.loadServer;
                
                pause(5);
                
                obj.hndl = actxserver('AstrooptikServer.Telescope');
            
                if ~obj.checkServer
                    return;
                end
                
                obj.hndl.SiteLatitude  = obj.object.latitude;
                obj.hndl.SiteLongitude = obj.object.longitude;
                
                obj.hndl.Connected = 1;
            
                obj.hndl.MotorOn;
                
                obj.update;
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function connectArduino(obj)

            obj.log.input('Connecting Arduino.');
            
            if isempty(obj.ard) 

                try
                    
                    if isempty(obj.ard)
                        obj.ard = obs.sens.ScopeAssistant;
                    end
                    
                catch ME
%                     obj.use_accelerometer = 0;
                    warning(ME.getReport);
                end

            end
            
            if ~isempty(obj.ard)
                
                try
                    
                    if isempty(obj.ard.telescope)
                        obj.ard.telescope = obj;
                    end
                    obj.ard.connect;

                    obj.ard.update;
                    
                catch ME
%                     obj.use_accelerometer = 0;
                    warning(ME.getReport);
                end
                
            end
            
        end
        
        function val = checkServer(obj)
            
            val = 0;
            
            [~, str] = system('tasklist /fi "imagename eq AstroOptikServer.exe"'); % check 
                
            str = strip(str);
            idx = strfind(str, 'AstroOptikServer.exe');
            
            if ~isempty(idx)
                
                num = util.text.extract_numbers(str(idx:end)); % found the server, now check the memory usage
                
                if num{1}(3)>30
                    val = 1;
                end
                
            end
            
        end
        
        function loadServer(obj)
            
%             system('C:\Program Files (x86)\Autoslew\AstroOptikServer.exe &'); % call the system command to load the server outside of matlab 
%             system('D:\matlab\wfast\+obs\+mount\launch_server.bat &'); % call the batch file to load the server then exit the cmd
            system(fullfile(getenv('WFAST'), '+obs\+mount\launch_server.bat &')); % call the batch file to load the server then exit the cmd

            tic;
            
            for ii = 1:100
                
                if obj.checkServer
                    return;
                else
                    pause(0.05);
                end
                
            end
            
            obj.killServer;
            error('Timeout while waiting for AstroOptikServer.exe to load for %f seconds!', toc); 
            
        end
        
        function disconnect(obj)
            
            obj.log.input('Disconnecting from mount');
            
            try
                obj.killServer;
                delete(obj.hndl);
                obj.hndl = [];
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
            
        end
        
        function killServer(obj)
            
            [ret, str] = system('taskkill /fi "imagename eq AstroOptikServer.exe" /f'); % force killing of the server.
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.object.reset;
            
        end
        
        function resetPrevObjects(obj)
            
            obj.prev_objects = {};
            
        end
        
    end
    
    methods % getters
        
        function val = get.limit_alt(obj)
            
            val = obj.hndl.AltitudeLimit;
            
        end
        
        function val = get.limit_flip(obj)
            
            val = obj.hndl.MeridianFlipMaxAngle;
            
        end
        
        function val = get.objName(obj)
            
            if isempty(obj.sync) || isempty(obj.sync.outgoing) || ~isfield(obj.sync.outgoing, 'OBJECT')
                val = '';
            else
                val = obj.sync.outgoing.OBJECT;
            end
            
        end
        
        function val = get.objRA_deg(obj)

            val = obj.object.RA_deg;
            
        end
        
        function val = get.objDEC_deg(obj)
            
            val = obj.object.Dec_deg;
            
        end
        
        function val = get.DE_target_deg(obj)
            
            val = obj.object.Dec_deg;
            
        end
        
        function val = get.objRA(obj)
            
            val = obj.object.RA; 
            
        end
        
        function val = get.objDEC(obj)
            
            val = obj.object.DE;
            
        end
        
        function val = get.DE_target(obj)
            
            val = obj.object.DE;
            
        end
        
        function val = get.objHA(obj)
            
            val = obj.object.HA;
            
        end
        
        function val = get.objALT(obj)
            
            val = obj.object.Alt_deg;
            
        end
        
        function val = get.objAZ(obj)
            
            val = obj.object.Az_deg;
            
        end
        
        function val = get.telRA_deg(obj)
            
            try 
                val = obj.hndl.RightAscension*15; % convert hours to degrees!
            catch
                val = [];
            end
            
        end
        
        function val = get.telDEC_deg(obj)
            
            try
                val = obj.hndl.Declination;
            catch
                val = [];
            end
            
        end
        
        function val = get.telDE_deg(obj)
            
            try
                val = obj.hndl.Declination;
            catch
                val = [];
            end
            
        end
        
        function val = get.telRA(obj)
            
            val = head.Ephemeris.deg2hour(obj.telRA_deg);
            
        end
        
        function val = get.telDEC(obj)
            
            val = head.Ephemeris.deg2sex(obj.telDEC_deg);
            
        end
        
        function val = get.telDE(obj)
            
            val = head.Ephemeris.deg2sex(obj.telDEC_deg);
            
        end
        
        function val = get.telALT(obj)
            
            try
                val = obj.hndl.Altitude;
            catch
                val = [];
            end
            
        end
        
        function val = get.telAZ(obj)
            
            try 
                val = obj.hndl.Azimuth;
            catch 
                val = [];
            end
            
        end
        
        function val = get.telHA(obj)
            
            val = head.Ephemeris.deg2hour(obj.LST_deg - obj.telRA_deg);
            
        end
        
        function val = LST_deg(obj)
            
            try 
                val = obj.hndl.SiderealTime*15; % in degrees...
            catch 
                val = [];
            end
            
        end
        
        function val = get.LST(obj)
            
            val = head.Ephemeris.deg2hour(obj.LST_deg);
            
        end
        
        function val = get.tracking(obj)
            
            try 
                val = obj.hndl.Tracking;
            catch 
                val = [];
            end
        
        end
        
        function val = latitutde(obj)
            
%             val = obj.hndl.SiteLatitude;
            val = obj.object.latitude;
            
        end
        
        function val = longitude(obj)
           
%             val = obj.hndl.SiteLongitude;
            val = obj.object.longitude;

        end
        
        function val = get.rate_RA(obj)
            
            val = obj.hndl.RightAscensionRate;
            
        end
        
        function val = get.rate_DE(obj)
            
            val = obj.hndl.DeclinationRate;
            
        end
        
        function val = get.pier_side(obj)
            
            val = obj.hndl.SideOfPier;
            
        end
        
        function val = obj_pier_side(obj)
            
            if isempty(obj.objRA_deg) || isempty(obj.objDEC_deg)
                val = '';
            else
                val = obj.hndl.DestinationSideOfPier(obj.objRA_deg/15, obj.objDEC_deg);
            end
            
        end
        
        function val = tel_time_to_limit(obj)
            
            val = NaN; % need to implement this! 
            
        end
        
        function val = obj_time_to_limit(obj)
            
            val = NaN; % need to implement this! 
            
        end
        
    end
    
    methods % setters
        
        function set.objName(obj, val)
            
            if ~isempty(obj.sync)
                obj.sync.outgoing.OBJECT = val;
            end
            
        end
        
        function set.objRA(obj, val)
            
            obj.object.RA = val; % note this is given in HOURS!
            obj.object.update;
            obj.objName = '';
            
            try
                obj.hndl.TargetRightAscension = obj.object.RA_deg/15; % also update the telescope's target field...
            catch ME
%                 obj.log.error(ME.getReport);
%                 rethrow(ME);
            end
            
        end
        
        function set.objRA_deg(obj, val)
            
            obj.object.RA_deg = val;
            obj.object.update;
            obj.objName = '';
            
            try
                obj.hndl.TargetRightAscension = obj.object.RA_deg/15; % also update the telescope's target field...
            catch ME
%                 obj.log.error(ME.getReport);
%                 rethrow(ME);
            end
            
        end
        
        function set.objDEC(obj, val)
            
            obj.object.Dec = val;
            obj.object.update;
            obj.objName = '';
            
            try
                obj.hndl.TargetDeclination = obj.object.Dec_deg;
            end
            
        end
        
        function set.objDEC_deg(obj, val)
            
            obj.object.Dec_deg = val;
            obj.object.update;
            obj.objName = '';
            
            try
                obj.hndl.TargetDeclination = obj.object.Dec_deg;
            end
            
        end
        
        function set.DE_target_deg(obj, val)
            
            obj.object.Dec_deg = val;
            obj.object.update;
            obj.objName = '';
            
            try
                obj.hndl.TargetDeclination = obj.object.Dec_deg;
            end
            
        end
        
        function set.tracking(obj, val)
            
            try 
                
                obj.was_tracking = val;
                
                if obj.hndl.Tracking~=val

                    obj.hndl.Tracking = val;

                    res = 0.01;
                    N = 100;
                    tic;
                    for ii = 1:N

                        if obj.hndl.Tracking==val
                            return;
                        end

                        pause(res);

                    end

                    error('Timeout after %f seconds waiting for mount to set tracking to %d.', toc, val);

                end
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function set.rate_RA(obj, val)
            
            obj.hndl.RightAscensionRate = val;
            
        end
        
        function set.rate_DE(obj, val)
            
            obj.hndl.DeclinationRate = val;
            
        end
        
    end
    
    methods % calculations / commands
        
        function inputTarget(obj, varargin)
            
            if isempty(varargin) && isempty(obj.objName)
                error('Must supply an object name or coordinates (or fill objName field)');
            elseif isempty(varargin)
                varargin{1} = obj.objName;
            end
            
            if length(varargin)>=2
                obj.objName = '';
            end
            
            obj.object.input(varargin{:}); % Ephemeris now uses Eran's name resolver
            
        end
        
        function parseTargetString(obj, str)
            
            C = strip(strsplit(str, {'{','}',','}));
            
            C = C(~cellfun(@isempty,C));
            
            if length(C)==2
                obj.objName = '';
                obj.objRA = C{1};
                obj.objDEC = C{2};
            elseif length(C)==3
                obj.objRA = C{2};
                obj.objDEC = C{3};
                obj.objName = C{1};
            else
                error('Wrong number of extracted strings (%d).', length(C));
            end
            
            
        end
        
        function val = check_before_slew(obj)
            
            val = 0;
            
            obj.object.update;
            
            if obj.object.Alt_deg<obj.limit_alt
                disp(['Target alt (' num2str(obj.object.Alt_deg) ') is below alt limit (' num2str(obj.limit_alt) ')']);
                return;
            end
            
            if obj.use_accelerometer
                if obj.ard.status==0
                    
                    obj.ard.update;
                    
                    pause(0.1);
                    
                    if isempty(obj.ard) || obj.ard.status==0 || obj.checkArduinoTime==0
                        obj.connectArduino;
                    end
                    
                end
                
                if obj.ard.status==0
                    error('Cannot slew without a responsive ScopeAssistant (arduino)');
                end
                
            end
            
            try
                if isempty(obj.audio)
                    obj.audio = util.sys.AudioControl;
                end
                
                obj.audio.playBlackAlert;
                
            catch ME
                warning(ME.getReport);
            end
            
            val = 1;
            
        end
        
        function val = check_while_moving(obj)
            
            val = 0;
            
            try 
                
                if obj.telALT<obj.limit_alt
                    return;
                end

                val = 1;
            
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function val = check_need_flip(obj)
            
            if strcmp(obj.obj_pier_side, 'pierUnknown')
                error('Mount cannot find PierSide for this object!');
            end
            
            val = ~strcmp(obj.pier_side, obj.obj_pier_side);
            
        end
        
        function slew(obj, varargin)
            
            if isempty(obj.object.RA) || isempty(obj.object.DE) 
                error('Please provide a target with viable RA/DE');
            end
            
            obj.log.input(sprintf('Slewing to target. RA= %s | DE= %s | ALT= %4.2f', obj.object.RA, obj.object.DE, obj.object.ALT_deg));
            
            try 
                
                if ~obj.check_before_slew
                    error('Prechecks failed, aborting slew');
                end
                
                if obj.check_need_flip
                    
                    if obj.telDE_deg<30 || obj.telDE_deg>60
                        obj.slewWithoutPrechecks(obj.hndl.RightAscension, 45); % do a preslew to dec +70 so we can make the flip! 
                    end
                    
                    
                    % other things to do before a flip...?
                    
                end
                
                try % keep a history of all targets
                    str = [obj.objName ' {' obj.objRA ',' obj.objDEC '}'];
                    if isempty(obj.prev_objects) || all(~strcmp(str, obj.prev_objects))
                        obj.prev_objects = [str obj.prev_objects];
                        if length(obj.prev_objects)>obj.num_prev_objects
                            obj.prev_objects = obj.prev_objects(1:obj.num_prev_objects);
                        end
                    end
                catch ME
                    rethrow(ME);
                end
                
                ra_hours_Jnow = obj.object.RA_deg_now./15; % convert to hours! 
                dec_deg_Jnow = obj.object.Dec_deg_now;
                
                obj.slewWithoutPrechecks(ra_hours_Jnow, dec_deg_Jnow);
                
                if obj.use_motor_toggle
                    
                    obj.hndl.MotorOff; 
                    
                    pause(0.01);
                    
                    obj.hndl.MotorOn;
                    
                    obj.slewWithoutPrechecks(ra_hours_Jnow, dec_deg_Jnow);
                    
                end
                
                obj.tracking = 1;
                
                if ~isempty(obj.sync)
                    obj.sync.outgoing.stop_camera = 0; % need to tell cam-pc to start working! 
                    obj.updateCamera;
                    obj.sync.update;
                end
                
                obj.setup_timer;
                
                if ~isempty(obj.sync)
                    obj.sync.incoming = struct;
                end
                
                obj.guiding_history.reset;
                obj.guiding_history.N = 100; 
                obj.resetRate;
                
                obj.audio.stop;
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function slewWithoutPrechecks(obj, ra_hours_Jnow, dec_deg_Jnow)
            
            obj.brake_bit = 0;
            
            obj.hndl.SlewToCoordinatesAsync(ra_hours_Jnow, dec_deg_Jnow);

            for ii = 1:100000

                pause(0.01);

                if obj.brake_bit
                    obj.stop;
                    break;
                end

                if ~obj.check_while_moving
                    obj.emergency_stop;
                    break;
                end

                if obj.hndl.Slewing==0
                    break;
                end

            end

        end
        
        function setup_timer(obj, ~, ~)
            
            if ~isempty(obj.timer) && isa(obj.timer, 'timer') && isvalid(obj.timer)
                if strcmp(obj.timer.Running, 'on')
                    stop(obj.timer);
                    delete(obj.timer);
                    obj.timer = [];
                end
            end
            
            delete(timerfind('name', 'mount-timer'));
            
            obj.timer = timer('BusyMode', 'queue', 'ExecutionMode', 'fixedRate', 'Name', 'mount-timer', ...
                'Period', 2, 'StartDelay', 1, 'TimerFcn', @obj.callback_timer, 'ErrorFcn', @obj.setup_timer);
            
            start(obj.timer);
            
        end
        
        function callback_timer(obj, ~, ~) % update sensors and GUI
            
            if isempty(obj) || isempty(obj.hndl)
                return;
            end
            
            try 
            
                if obj.use_guiding && obj.tracking % && obj.sync.status
                    
                    direction = 1;
                    if strcmp(obj.pier_side, 'pierEast')
                        direction  = -1;
                    end
                    
                    if ~isempty(obj.sync.incoming) && isfield(obj.sync.incoming, 'RA_rate_delta') && ~isempty(obj.sync.incoming.RA_rate_delta)
                        obj.rate_RA = obj.rate_RA + direction*obj.sync.incoming.RA_rate_delta;
                        obj.sync.incoming.RA_rate_delta = 0; % must zero this out, so if we lose connection we don't keep adding these deltas
                        obj.sync.outgoing.RA_rate = obj.rate_RA;
                    end
                    
                    if ~isempty(obj.sync.incoming) && isfield(obj.sync.incoming, 'DE_rate_delta') && ~isempty(obj.sync.incoming.DE_rate_delta)
                        obj.rate_DE = obj.rate_DE + direction*obj.sync.incoming.DE_rate_delta;
                        obj.sync.incoming.DE_rate_delta = 0; % must zero this out, so if we lose connection we don't keep adding these deltas
                        obj.sync.outgoing.DE_rate = obj.rate_DE;
                    end
                    
                else
                    % what to do here? reconnect or leave that to t1?
                end
                
                if ~isempty(obj.rate_RA) && ~isempty(obj.rate_DE) 
                    obj.guiding_history.input([obj.rate_RA, obj.rate_DE]);
                end
                
                if ~isempty(obj.gui) && obj.gui.check
                    obj.gui.update;
                end
                
                if ~isempty(obj.owner) && ~isempty(obj.owner.gui) && obj.owner.gui.check
                    obj.owner.gui.updateStopButton;
                end
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function resetRate(obj)
            
            obj.rate_RA = 0;
            obj.rate_DE = 0;
                
        end
        
        function adjustPosition(obj, RA_deg, DE_deg)
            
            error('This doesnt work, dont use it');
            
            obj.log.input(sprintf('Adjusting position by RA: %f deg | DE: %f deg', RA_deg, DE_deg));
            
            try

                total_time = 1;
                time_res = 0.1;
                N = ceil(total_time./time_res);
                total_time = N.*time_res; % if the division was not integer

                if obj.prechecks==0
                    obj.stop;
                    return;
                end

                obj.brake_bit = 0;
                on_cleanup = onCleanup(@obj.stop);
            
                % convert the total adjustment to the rate in arcsec/sec
                obj.hndl.RightAscensionRate = DE_deg.*3600/total_time;
                obj.hndl.DeclinationRate = DE_deg.*3600/total_time;
            
                t = tic;
                
                for ii = 1:N
                
                    if obj.brake_bit || obj.checks==0
                        obj.emergency_stop;
                        return;
                    end
                    
                    if toc(t)>total_time
                        obj.stop;
                        return;
                    end
                    
                    pause(time_res);
                    
                end
            
            catch ME
                obj.stop;
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
                
        end
        
        function engineeringSlew(obj,Alt,Az)
            
            try 

                if Alt<obj.limit_alt
                    error('Input altitude is below limit of %f degress', obj.limit_alt);
                end
                
                obj.brake_bit = 0;
                
                obj.tracking = 0;
                obj.hndl.SlewToAltAzAsync(Alt, Az);
                
                for ii = 1:100000
                    
                    pause(0.01);
                    
                    if obj.brake_bit
                        obj.stop;
                        break;
                    end
                    
                    if ~obj.check_while_moving
                        obj.emergency_stop;
                        break;
                    end
                    
                    if obj.hndl.Slewing==0
                        break;
                    end
                    
                end
            
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
            
        end
        
        function syncToTarget(obj)
            
%             obj.hndl.SyncToTarget;
            obj.hndl.SyncToCoordinates;
            
        end
        
        function park(obj)
            
            obj.hndl.Park;
            obj.tracking = 0;
            
        end
        
        function zenith_west(obj)
            
        end
        
        function zenith_east(obj)
            
        end
        
        function update(obj)
            
            obj.object.update;
            
            try
                obj.hndl.Connected;
            catch 
                obj.status = 0; 
                return;
            end
            
            % check server is running
            [~, str] = system('tasklist /fi "imagename eq AstroOptikServer.exe"'); % check 
                
            str = strip(str);
            idx = strfind(str, 'AstroOptikServer.exe');
            
            if isempty(idx)
                obj.status = 0;
                return;
            end
            
            if obj.telALT<-20 % this occurs when software is disconnected from mount (e.g., on power out)
                obj.status = 0;
                return;
            end
            
            obj.setup_timer;
            
            try % logger hearteat
                if ~obj.log.check_heartbeat
                    obj.log.heartbeat(300, obj); 
                end
            catch ME
                warning(ME.getReport);
            end
            
            try % arduino
               
                if obj.use_accelerometer
                    if isempty(obj.ard) || obj.ard.status==0 || obj.checkArduinoTime==0
                        obj.connectArduino;
                    end
                end
                
            catch ME
                warning(ME.getReport);
            end
            
            % add additional tests?
            
            obj.status = 1;
            
        end
        
        function stop(obj)
            
%             obj.log.input('stopping telescope');
            
            try 
                
                obj.brake_bit = 1;
            
                obj.hndl.AbortSlew;
                obj.hndl.tracking = 0;
                obj.resetRate;
                
                if ~isempty(obj.sync)

                    if obj.sync.outgoing.stop_camera==0
                        obj.log.input('Telescope stopped, sending camera stop command');
                        disp(obj.log.report);
                    end

                    obj.sync.outgoing.stop_camera = 1;
                    
                end
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function emergency_stop(obj)
            
            obj.log.input('stopping telescope');
            
            disp('stopping telescope!');
            
            obj.stop;
            
        end
        
        function updateCamera(obj)
            
            obj.sync.outgoing.RA = obj.objRA;
            obj.sync.outgoing.DEC = obj.objDEC;            
            obj.sync.outgoing.RA_DEG = obj.objRA_deg;
            obj.sync.outgoing.DEC_DEG = obj.objDEC_deg;
            obj.sync.outgoing.TELRA = obj.telRA;
            obj.sync.outgoing.TELDEC = obj.telDEC;
            obj.sync.outgoing.TELRA_DEG = obj.telRA_deg;
            obj.sync.outgoing.TELDEC_DEG = obj.telDEC_deg;
            
            if isempty(obj.tracking) || obj.tracking==0 || (~isempty(obj.objRA_deg) && abs(obj.objRA_deg-obj.telRA_deg)>1) % if mount stops tracking or is 4 time-minutes away from target RA, stop the camera (e.g., when reaching limit)
                obj.sync.outgoing.stop_camera = 1;
            else
%                 obj.sync.outgoing.stop_camera = 0;
            end
            
        end
        
        function val = checkArduinoTime(obj)
            
            val = 1;
            
            if ~isempty(obj.log.time)
            
                dt = seconds(obj.log.time- obj.ard.time);
            
                if dt>600 % arduino has not updated in a few minutes
                    val = 0;
                    return;
                end
                
            end
            
        end
        
        function val = printout(obj)
        
            val = sprintf('Status= %d, LST= %s, RA= %s, Dec= %s, HA= %s, ALT= %4.2f, %s', ...
                obj.status, obj.LST, obj.telRA, obj.telDE, obj.telHA, obj.telALT, obj.pier_side); 
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = obs.mount.gui.ASAGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
        function plot_rate(obj, ax)
            
            if ~isempty(obj.guiding_history.data)
                hold(ax, 'off');
                plot(ax, obj.guiding_history.data_ordered(:,1), 'x');
                hold(ax, 'on');
                plot(ax, obj.guiding_history.data_ordered(:,2), '+');
                ylabel(ax, 'rates');
                hold(ax, 'off');
            end
            
        end
        
    end    
    
    properties(Transient=true, Hidden=true, Dependent=true) % for backward compatibility
        
        DE_target_deg;
        DE_target;
        telDE_deg;
        telDE;
        
    end
    
end

