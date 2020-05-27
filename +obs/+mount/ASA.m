classdef (CaseInsensitiveProperties, TruncatedProperties) ASA < handle
% Control ASA mount. 
% This class works mainly by giving an "object" and then letting the mount
% slew to that object. It also has shortcut fields to change various 
% parameters in the ASCOM mount object itself (which is in "hndl"). 
%
% To give a target use inputTarget(star_name) or inputTarget(RA,DEC). 
% Make sure the target parmeters are correctly translated to the object 
% coordinates objRA and objDec. Also make sure the objALT is high enough. 
% To make the telescope go to target just use "slew" without arguments. 
% 
% NOTE: All parameters of the mount / target are given in hours/degrees as
%       numeric values or sexagesimal strings (HH:MM:SS or DD:MM:SS). 
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
        
        audio@util.sys.AudioControl; % to play warning sound when slewing
        
        prev_objects = {}; % list of latest objects 
        
    end
    
    properties % objects
        
        hndl; % ASCOM object
        
        owner@obs.Manager; % link back to the Manager object
        
        object@head.Ephemeris; % the object/target/field that we want to point at. 
        
        ard@obs.sens.ScopeAssistant; % connect to accelerometer and ultrasonic sensor
        
        cam_pc@obs.comm.PcSync; % communications object to camera PC, to be given from Manager to pass data on to camera computer
        
        timer; % timer object
        
        guiding_history@util.vec.CircularBuffer; % latest values from camera guiding system
        correct_history@util.vec.CircularBuffer; % latest corrections applied when guiding
        
        reco@obs.comm.Reconnect; % controls reconnection attempts, and prevents recurring failed attempts
        
        log@util.sys.Logger; % writes all commands and errors to text file
        
    end
    
    properties % inputs/outputs
        
        status = 0; % if connected and responsive, set this to 1
        
    end
    
    properties % switches/controls
        
        use_guiding = 1; % apply the corrections from cam_pc
        use_integral_guiding = 1; % use integration over previous guiding data to adjust rates
        use_pid_guiding = 0; % use more complicated 
        integral_number = 50; 
        pid_P = 0.3;
        pid_I = 0.6;
        pid_D = 0.1;
        
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
        objDec_deg;
        
        objRA; % string in HH:MM:SS format
        objDec; % string in DD:MM:SS format
        
        % these are read-only values for inspection
        objHA; % string in HH:MM:SS format
        objALT; % degrees
        objAZ; % degrees 
        
        % these are read out directly from the mount        
        % these are read-only values for inspection
        telRA_deg; 
        telDec_deg;
        
        telRA; % string in HH:MM:SS format
        telDec; % string in DD:MM:SS format
        
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
       
        % the default alt-az is parking 2 (zenith west)
        target_altitude = 90;
        target_azimuth = 270;
        default_target_altitude = 90;
        default_target_azimuth = 270; 
        
        was_tracking = 0; % when pressing the NSEW buttons, to temporarily turn off tracking (and maybe bring it back on)
        default_move_rate; % move rate for the NSEW buttons (in degrees/second)

        prev_dRA = 0; % tracking rate used before
        prev_dDE = 0; % tracking rate used before

        num_prev_objects = 10; % how many prev_objects do we keep?
        
        version = 1.04;
        
    end
    
    methods % constructor and connect commands
        
        function obj = ASA(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.mount.ASA')
                if obj.debug_bit>1, fprintf('ASA copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('ASA constructor v%4.2f\n', obj.version); end
            
            end
            
            obj.guiding_history = util.vec.CircularBuffer;
            obj.correct_history = util.vec.CircularBuffer;
            
            obj.reco =obs.comm.Reconnect;
            
            obj.log = util.sys.Logger('ASA_mount', obj);
            
            try
            
                obj.object = head.Ephemeris;

                try
                    obj.connect;
                catch ME
                    warning(ME.getReport); 
                end

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
            
            if obj.reco.should

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
                    
                    if obj.status==0
                        error('Mount connected but status is 0...'); 
                    end
                    
                    obj.reco.inputSuccess;

                catch ME
                    obj.reco.inputFailure(ME.getReport);
                    obj.log.error(ME.getReport);
                    rethrow(ME);
                end

            end
            
        end
        
        function connectArduino(obj)

            obj.log.input('Connecting Arduino.');
            
            if isempty(obj.ard) 

                try
                    
                    obj.ard = obs.sens.ScopeAssistant;
                    
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
%                     obj.use_accelerometer = 0; % this prevents endless loops of reconnect but also doesn't allow multiple reconnects
                    warning(ME.getReport);
                end
                
            end
            
        end
        
        function val = checkServer(obj) % check if the Autoslew application is running
            
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
        
        function loadServer(obj) % restart the Autoslew application
            
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
        
        function disconnect(obj) % remove the server and clear the hndl
            
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
        
        function killServer(obj) % kill the applications using the system "taskkill" command
            
            [ret, str] = system('taskkill /fi "imagename eq AstroOptikServer.exe" /f'); % force killing of the server.
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.object.reset;
            
        end
        
        function resetPrevObjects(obj) % clear recent target list
            
            obj.prev_objects = {};
            
        end
        
    end
    
    methods % getters
        
        function val = get.limit_alt(obj)
            
            try
                val = obj.hndl.AltitudeLimit;
            catch
                val = [];
            end
        end
        
        function val = get.limit_flip(obj)
            
            try
                val = obj.hndl.MeridianFlipMaxAngle;
            catch
                val = [];
            end
            
        end
        
        function val = get.objName(obj)
            
%             if isempty(obj.cam_pc) || isempty(obj.cam_pc.outgoing) || ~isfield(obj.cam_pc.outgoing, 'OBJECT')
%                 val = '';
%             else
%                 val = obj.cam_pc.outgoing.OBJECT;
%             end

            val = obj.object.name;

            
        end
        
        function val = get.objRA_deg(obj)

            val = obj.object.RA_deg;
            
        end
        
        function val = get.objDec_deg(obj)
            
            val = obj.object.Dec_deg;
            
        end
        
        function val = get.DE_target_deg(obj)
            
            val = obj.object.Dec_deg;
            
        end
        
        function val = get.objRA(obj)
            
            val = obj.object.RA; 
            
        end
        
        function val = get.objDec(obj)
            
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
        
        function val = get.telDec_deg(obj)
            
            try
                val = obj.hndl.Declination;
            catch
                val = [];
            end
            
        end
        
        function val = get.telRA(obj)
            
            val = head.Ephemeris.deg2hour(obj.telRA_deg);
            
        end
        
        function val = get.telDec(obj)
            
            val = head.Ephemeris.deg2sex(obj.telDec_deg);
            
        end
        
        function val = get.telDE(obj)
            
            val = head.Ephemeris.deg2sex(obj.telDec_deg);
            
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
            
            ha = obj.LST_deg - obj.telRA_deg;
            
            if ha>180
                ha = ha -360;
            end
            
            val = head.Ephemeris.deg2hour(ha);
            
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
        
        function val = is_slewing(obj)
            
            try
                val = obj.hndl.Slewing;
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
            
            try 
                val = obj.hndl.RightAscensionRate;
            catch
                val = [];
            end
            
        end
        
        function val = get.rate_DE(obj)
            
            try
                val = obj.hndl.DeclinationRate;
            catch
                val = [];
            end
                
        end
        
        function val = get.pier_side(obj)
            
            try
                val = obj.hndl.SideOfPier;
            catch ME
                val = [];
            end
            
        end
        
        function val = telHemisphere(obj)
           
            if isempty(obj.pier_side)
                val = [];
            elseif strcmp(obj.pier_side, 'pierEast')
                val = 'West';
            elseif strcmp(obj.pier_side, 'pierWest')
                val = 'East';
            end
            
        end
        
        function val = obj_pier_side(obj)
            
            if isempty(obj.objRA_deg) || isempty(obj.objDec_deg)
                val = '';
            else
                try
                    val = obj.hndl.DestinationSideOfPier(obj.objRA_deg/15, obj.objDec_deg);
                catch 
                    val = '';
                end
            end
            
        end
        
        function val = objHemisphere(obj)
            
            side = obj.obj_pier_side;
            
            if isempty(side)
                val = '';
            elseif strcmp(side, 'pierEast')
                val = 'West'; 
            elseif strcmp(side, 'pierWest')
                val = 'East';
            else
                val = 'unobservable'; 
            end
            
        end
        
        function val = tel_time_to_limit(obj)
            
            val = NaN; % need to implement this! 
            
        end
        
        function val = obj_time_to_limit(obj)
            
            val = obj.object.calcObsTimeMinutes;
            
        end
        
    end
    
    methods % setters
        
        function set.objName(obj, val)
            
            obj.object.name = val; 
            
            if ~isempty(obj.cam_pc)
                obj.cam_pc.outgoing.OBJECT = val;
            end
            
        end
        
        function set.objRA(obj, val)
            
            obj.object.RA = val; % note this is given in HOURS!
            obj.object.update;
            obj.object.name = '';
            obj.object.keyword = '';
            
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
            obj.object.name = '';
            obj.object.keyword = '';
            
            try
                obj.hndl.TargetRightAscension = obj.object.RA_deg/15; % also update the telescope's target field...
            catch ME
%                 obj.log.error(ME.getReport);
%                 rethrow(ME);
            end
            
        end
        
        function set.objDec(obj, val)
            
            obj.object.Dec = val;
            obj.object.update;
            obj.objName = '';
            
            try
                obj.hndl.TargetDeclination = obj.object.Dec_deg;
            end
            
        end
        
        function set.objDec_deg(obj, val)
            
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
            
            if isempty(obj.hndl) || obj.status==0
                return;
            end
            
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
            
            if ~isempty(obj.hndl) && obj.status
                obj.hndl.RightAscensionRate = val; % arcsec per second
            end
            
        end
        
        function set.rate_DE(obj, val)
            
            if ~isempty(obj.hndl) && obj.status
                obj.hndl.DeclinationRate = val; % arcsec per second
            end
            
        end
        
    end
    
    methods % parsing targets
        
        function inputTarget(obj, varargin) % pass the arguments to the Ephemeris object (give RA and Dec or object name)
            
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
        
        function C_out = scanTargetList(obj, filename) % read list of targets from text file
            
            if ~exist(filename, 'file')
                error('Cannot find filename "%s".', filename);
            end
            
            f = fopen(filename, 'r'); 
            clean_up_file = onCleanup(@() fclose(f));
            
            C = {};
            
            for ii = 1:1e3 % arbitrary loop length with break inside
                
                L = fgetl(f); 
                
                if isnumeric(L)
                    break;
                end
                
                [name, RA, Dec ] = obj.parseTargetString(L); 
                
                if isempty(RA) || isempty(Dec)
                    str = name;
                elseif isempty(name)
                    str = sprintf('{%s, %s}', RA, Dec);
                else
                    str = sprintf('%s {%s, %s}', name, RA, Dec);
                end
                
                C{end+1} = str;
                
            end
            
            if nargout==0
                
                if length(C)>obj.num_prev_objects
                    obj.num_prev_objects = length(C);
                end
                    
                for ii = length(C):-1:1
                    obj.addTargetList(C{ii}); 
                end
               
            else
                C_out = C;
            end
            
        end
        
        function [name, RA, Dec] = parseTargetString(obj, str) % get a single line of text and turn it into a target for the list
            
            e = head.Ephemeris; % use this to format the RA/Dec inputs
            
            C = strip(strsplit(str, {'{','}',','}));
            
            C = C(~cellfun(@isempty,C));
            
            for ii = 1:min(2,length(C))
                
                if contains(C{ii}, '+')
                    new_cells = strip(strsplit(C{ii}, '+'));
                    new_cells{2} = ['+' new_cells{2}];
                elseif contains(C{ii}, '-')
                    new_cells = strip(strsplit(C{ii}, '-'));
                    new_cells{2} = ['-' new_cells{2}];
                else
                    new_cells = C{ii}; 
                end
                
                if ii<length(C)
                    C = [C(1:ii-1) new_cells C(ii+1:end)];
                else
                    C = [C(1:ii-1) new_cells];
                end
                
            end
            
            if length(C)==1
                name = C{1};
                RA = '';
                Dec = '';
            elseif length(C)==2
                name = '';
                RA = C{1};
                Dec = C{2};
            elseif length(C)==3
                name = C{1};
                RA = C{2};
                Dec = C{3};
            end
            
            if ~isempty(RA)
                e.RA = RA;
                RA = e.RA; 
            end
            
            if ~isempty(Dec)
                e.Dec = Dec;
                Dec = e.Dec; 
            end
            
            if isempty(RA) || isempty(Dec)
                e.input(name, []);
                if ~isnan(e.RA_deg) && ~isnan(e.Dec_deg)
                    RA = e.RA;
                    Dec = e.Dec;
                end
            end
            
            if nargout==0
                obj.objRA = RA;
                obj.objDec = Dec;
                obj.objName = name;
            end
            
        end
        
        function addTargetList(obj, str) % add an object to the list of recent targets
            
            if nargin<2 || isempty(str)
                str = [obj.objName ' {' obj.objRA ', ' obj.objDec '}'];
            end
            
            if isempty(obj.prev_objects) || all(~strcmp(str, obj.prev_objects)) % check this string is not on the list already... 
                obj.prev_objects = [str obj.prev_objects]; % add new items up front
                if length(obj.prev_objects)>obj.num_prev_objects % clip the end of the list if needed
                    obj.prev_objects = obj.prev_objects(1:obj.num_prev_objects);
                end
            end
            
        end
        
    end
    
    methods % calculations / commands
        
        function val = check_before_slew(obj) % check object is inside limits etc.
            
            val = 0;
            
            obj.object.update; % the Ephemeris object can calculate the ALT and whatever other parameters
            
            if obj.object.Alt_deg<obj.limit_alt
                disp(['Target alt (' num2str(obj.object.Alt_deg) ') is below alt limit (' num2str(obj.limit_alt) ')']);
                return;
            end
            
            if obj.use_accelerometer
                
                if ~isempty(obj.ard) && obj.ard.status==0
                    
                    ok = obj.ard.update;
                    
                    pause(0.1);
                    
                    if isempty(obj.ard) || obj.ard.status==0 || obj.checkArduinoTime==0 || ok==0
                        obj.connectArduino;
                    end
                    
                end
                
                if isempty(obj.ard) || obj.ard.status==0
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
        
        function val = check_while_moving(obj) % altitude check to be done while slewing
            
            val = 0;
            
            try 
                
                if obj.telALT<obj.limit_alt
                    fprintf('Alt limit crossed!'); 
                    return;
                end

                % I think it is better to rely on the timer for these checks
%                 ok = obj.ard.update;
%                 
%                 if obj.ard.ALT<obj.limit_alt || ok==0
%                     return;
%                 end
                
                val = 1;
            
            catch ME
                warning(ME.getReport);
            end
            
        end
        
        function val = check_after_slew(obj)
            
            val = 1; 
            
            % the RA/Dec condition for a successfull slew is 10 arcmin,
            % because there is some bug in turning coordinates from J2000
            % to Jnow (not sure if the problem is in the mount software or
            % in Eran's coco function. 
            if abs(obj.telRA_deg-obj.objRA_deg)*3600>600
                val = 0; 
                return;
            end
            
            if abs(obj.telDec_deg-obj.objDec_deg)*3600>600
                val = 0;
                return;
            end
            
            pause(1); 
            
            v = obj.getVibrationStrength; 
            
            if sqrt(v)*3600>10
                val = 0;
                return;
            end
            
        end
        
        function val = check_need_flip(obj) % check if target is on the other side of the pier
            
            if strcmp(obj.obj_pier_side, 'pierUnknown')
                error('Mount cannot find PierSide for this object!');
            end
            
            val = ~strcmp(obj.pier_side, obj.obj_pier_side);
            
        end
        
        function slew(obj, varargin) % begin slewing after doing pre-checks
        % Usage: slew(obj, varargin)
        % Begin slewing to the current coordinates in "object". 
        % Will preform some prechecks (e.g., altitude of target) and then 
        % start to slew, making additional checks during the motion. 
        %
        % To stop the slew, use the mount.stop() method, or set the value
        % of "brake_bit" to 1 (e.g., by pressing the GUI stop button). 
        %
        % OPTIONAL ARGUMENTS:
        %   -RA: input a different right ascention to the "object". Must be
        %        in numeric hours or sexagesimal string of hours. 
        %   -Dec: input a different declination to the "object". Must be in 
        %         numeric degrees or sexagesimal degrees. 
        %   -skip prechecks: do not perform the prechecks and continue directly
        %                    to the slewing. Not recommended (default false). 
        %   -history: use this to add the current target to the list of 
        %             recent targets (default true). 
        
            input = util.text.InputVars;
            input.input_var('ra', [], 'right ascention'); 
            input.input_var('dec', [], 'declination'); 
            input.input_var('skip', false, 'skip prechecks');
            input.input_var('history', true, 'previous_targets'); 
            input.scan_vars(varargin{:}); 
            
            if ~isempty(input.ra)
                obj.object.RA = input.ra;
            end
            
            if ~isempty(input.dec)
                obj.object.Dec = input.dec;
            end
            
            if isempty(obj.object.RA) || isempty(obj.object.Dec) 
                error('Please provide a target with viable RA/DE');
            end
            
            obj.log.input(sprintf('Slewing to target. RA= %s | DE= %s | ALT= %4.2f', obj.object.RA, obj.object.Dec, obj.object.ALT_deg));
            
            try 
                
                if ~obj.check_before_slew
                    error('Prechecks failed, aborting slew');
                end
                
                % translate the object coordinates to current epoch
                ra_hours_Jnow = obj.object.RA_deg_now./15; % convert to hours! 
                dec_deg_Jnow = obj.object.Dec_deg_now;
                
                if obj.check_need_flip % do pre-slews to manage to do the flip without getting stuck! 
                    
                    if obj.telDec_deg<30 || obj.telDec_deg>60
                        
                        obj.slewWithoutPrechecks(obj.hndl.RightAscension, 45); % do a preslew to dec +70 so we can make the flip! 
                        
                        try 
                            if strcmp(obj.obj_pier_side, 'pierEast') % object is on the WEST SIDE!
                                obj.slewWithoutPrechecks(mod(obj.object.LST_deg/15+4,24), 45); % do a second preslew to the same Dec with RA that is easier to flip from
                            elseif strcmp(obj.obj_pier_side, 'pierWest') % object is on the EAST SIDE!
                                obj.slewWithoutPrechecks(mod(obj.object.LST_deg/15-4,24), 45); % do a second preslew to the same Dec with RA that is easier to flip from
                            end
                        catch ME
                            disp('Could not complete the second pre-slew:'); 
                            warning(ME.getReport); 
                        end
                        
                    end
                    
                    
                    % other things to do before a flip...?
                    
                end
                
                if input.history
                    try % keep a history of all targets
                        str = [obj.objName ' {' obj.objRA ', ' obj.objDec '}'];
                        obj.addTargetList(str);
                    catch ME
                        rethrow(ME);
                    end
                end
                
                obj.slewWithoutPrechecks(ra_hours_Jnow, dec_deg_Jnow);
                
                if obj.use_motor_toggle % this does not work! 
                    
                    obj.hndl.MotorOff; 
                    
                    pause(0.01);
                    
                    obj.hndl.MotorOn;
                    
                    obj.slewWithoutPrechecks(ra_hours_Jnow, dec_deg_Jnow);
                    
                end
                
                obj.tracking = 1;
                                
                pause(3); 
                
                % check that we've reached the right position and not vibrating 
                if ~obj.check_after_slew
                    
                    obj.log.input('Slew post-check failed, trying to slew again...');
                    disp(obj.log.report); 
                    
                    obj.slewWithoutPrechecks(ra_hours_Jnow, dec_deg_Jnow);
                    
                    obj.tracking = 1;
                
                    pause(1); 
                    
                    if ~obj.check_after_slew % if this fails, slew aside a little and come back
                        
                        obj.log.input('Slew post-check failed again, trying to slew to different coordinate and return...'); 
                        disp(obj.log.report); 

                        obj.slewWithoutPrechecks(ra_hours_Jnow+0.1, dec_deg_Jnow+1);
                        
                        pause(3); 
                        
                        obj.slewWithoutPrechecks(ra_hours_Jnow, dec_deg_Jnow);
                        
                        if ~obj.check_after_slew
                            obj.log.error('Slewing post-checks failed after all attempts to re-slew. '); 
                            error(obj.log.report); 
                        end
                        
                    end
                    
                end
                
                if ~isempty(obj.cam_pc)
                    obj.cam_pc.outgoing.stop_camera = 0; % need to tell cam-pc to start working! 
                    obj.updateCamera;
                    obj.cam_pc.update;
                end
                
                obj.setup_timer;
                
                if ~isempty(obj.cam_pc)
                    obj.cam_pc.incoming = struct;
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
        
        function slewWithoutPrechecks(obj, ra_hours_Jnow, dec_deg_Jnow) % internal command to move, to be called after doing pre-checks
            
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
        
        function [total_var, RA_var, DE_var] = getVibrationStrength(obj, varargin) % variance in RA/DE in units of degrees (run this only when tracking...)
            
            input = util.text.InputVars;
            input.input_var('N', 50, 'number', 'iterations'); 
            input.input_var('delay', 0.01, 'pause'); 
            input.scan_vars(varargin{:}); 
            
            RA = nan(input.N, 1); 
            DE = nan(input.N, 1); 
            
            for ii = 1:input.N
                
                RA(ii) = obj.telRA_deg; 
                DE(ii) = obj.telDec_deg;
                
                pause(input.delay); 
                
            end
            
            RA_var = nanvar(RA); 
            DE_var = nanvar(DE); 
            
            total_var = RA_var + DE_var;
            
        end
        
        function setup_timer(obj, ~, ~) % start the timer that checks Arduino and updates GUI
            
            if ~isempty(obj.timer) && isa(obj.timer, 'timer') && isvalid(obj.timer)
                if strcmp(obj.timer.Running, 'on')
                    stop(obj.timer);
                end
                
                delete(obj.timer);
                obj.timer = [];
                
            end
            
            delete(timerfind('name', 'mount-timer'));
            
            obj.timer = timer('BusyMode', 'queue', 'ExecutionMode', 'fixedRate', 'Name', 'mount-timer', ...
                'Period', 1, 'StartDelay', 1, 'TimerFcn', @obj.callback_timer, 'ErrorFcn', @obj.setup_timer);
            
            start(obj.timer);
            
        end
        
        function callback_timer(obj, ~, ~) % update sensors and GUI
            
            if isempty(obj) || isempty(obj.hndl)
                return;
            end
            
            try 
            
                obj.object.update;
                
                if obj.use_guiding && ~isempty(obj.tracking) && obj.tracking % && obj.cam_pc.status
                    
                    direction = 1;
                    if strcmp(obj.pier_side, 'pierEast')
                        direction  = -1;
                    end
                    
                    if ~isempty(obj.cam_pc.incoming) && isfield(obj.cam_pc.incoming, 'RA_rate_delta') && ~isempty(obj.cam_pc.incoming.RA_rate_delta)...
                            && isfield(obj.cam_pc.incoming, 'DE_rate_delta') && ~isempty(obj.cam_pc.incoming.DE_rate_delta)
                        
                        dRA = obj.cam_pc.incoming.RA_rate_delta;
                        if isempty(dRA) || isnan(dRA), dRA = 0; end
                        
                        dDE = obj.cam_pc.incoming.DE_rate_delta;
                        if isempty(dDE) || isnan(dDE), dDE = 0; end
                        
                        if ~isequal(dRA, obj.prev_dRA) && ~isequal(dDE, obj.prev_dDE)
                            
                            obj.prev_dRA = dRA;
                            obj.prev_dDE = dDE;
%                             fprintf('dRA= %6.4f | dDE= %6.4f\n', dRA, dDE); 
                            
                            obj.correct_history.input([dRA, dDE]); 
                            
                            % rates in arcsec per second
                            if obj.use_integral_guiding
                                obj.rate_RA = direction*obj.correct_history.mean(1); 
                                obj.rate_DE = direction*obj.correct_history.mean(2); 
                            elseif obj.use_pid_guiding
                                obj.rate_RA = obj.rate_RA + direction*obj.getPID_RA; 
                                obj.rate_DE = obj.rate_DE + direction*obj.getPID_DE;
                            else
                                obj.rate_RA = obj.rate_RA + direction*dRA;
                                obj.rate_DE = obj.rate_DE + direction*dDE;
                            end
                            
                            
                        end
                        
                    end
                    
%                     if ~isempty(obj.cam_pc.incoming) && isfield(obj.cam_pc.incoming, 'RA_rate_delta') && ~isempty(obj.cam_pc.incoming.RA_rate_delta)
%                         dRA = obj.cam_pc.incoming.RA_rate_delta;
%                         if isempty(dRA) || isnan(dRA), dRA = 0; end
%                         
%                         if ~isequal(dRA, obj.prev_dRA)
%                             obj.rate_RA = obj.rate_RA + direction*dRA;
%                             obj.prev_dRA = dRA;
%                             fprintf('dRA= %6.4f ', dRA); 
%                         end
%                         
% %                         obj.cam_pc.incoming.RA_rate_delta = 0; % must zero this out, so if we lose connection we don't keep adding these deltas
%                         
%                         obj.cam_pc.outgoing.RA_rate = obj.rate_RA;
%                         
%                     end
%                     
%                     if ~isempty(obj.cam_pc.incoming) && isfield(obj.cam_pc.incoming, 'DE_rate_delta') && ~isempty(obj.cam_pc.incoming.DE_rate_delta)
%                         dDE = obj.cam_pc.incoming.DE_rate_delta;
%                         if isempty(dDE) || isnan(dDE), dDE = 0; end
%                         
%                         if ~isequal(dDE, obj.prev_dDE)
%                             obj.rate_DE = obj.rate_DE + direction*dDE;
%                             obj.prev_dDE = dDE;
%                             fprintf('| dDE= %6.4f\n', dDE); 
%                         end
%                         
% %                         obj.cam_pc.incoming.DE_rate_delta = 0; % must zero this out, so if we lose connection we don't keep adding these deltas
%                         obj.cam_pc.outgoing.DE_rate = obj.rate_DE;
%                         
%                     end
                    
%                     fprintf('dRA= %6.4f | dDe= %6.4f\n', dRA, dDE); 
                    
                else
                    % what to do here? reconnect or leave that to t1?
                end
                
                if ~isempty(obj.rate_RA) && ~isempty(obj.rate_DE) 
                    obj.guiding_history.input([obj.rate_RA, obj.rate_DE]);
                end
                
                if obj.use_accelerometer
                    
                    if ~isempty(obj.ard)
                        ok = obj.ard.update;
                    end
                    
                    if isempty(obj.ard) || ok==0
                        obj.connectArduino;
                    end
                    
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
        
        function val = getPID_RA(obj) % Proportional-Integral-Differential control of RA rate
            
            P = obj.correct_history.data(obj.correct_history.idx,1); 
            I = nansum(obj.correct_history.data(:,1)); 
            
            idx = obj.correct_history.idx;
            idx = idx - 1;
            if idx<1
                idx = obj.correct_history.N;
            end
            
            if idx>0 && idx<obj.correct_history.N
                D = P-obj.correct_history.data(idx,1); 
            else
                D = 0;
            end
            
            K_total = obj.pid_P + obj.pid_I + obj.pid_D;
            
            K_P = obj.pid_P/K_total;
            K_I = obj.pid_I/K_total;
            K_D = obj.pid_D/K_total;
            
            val = K_P.*P + K_I.*I + K_D.*D;
            
        end
        
        function val = getPID_DE(obj) % Proportional-Integral-Differential control of Dec rate
            
            P = obj.correct_history.data(obj.correct_history.idx,2); 
            I = nansum(obj.correct_history.data(:,2)); 
            
            idx = obj.correct_history.idx;
            idx = idx - 1;
            if idx<1
                idx = obj.correct_history.N;
            end
            
            if idx>0 && idx<obj.correct_history.N
                D = P-obj.correct_history.data(idx,2); 
            else
                D = 0;
            end
            
            K_total = obj.pid_P + obj.pid_I + obj.pid_D;
            
            K_P = obj.pid_P/K_total;
            K_I = obj.pid_I/K_total;
            K_D = obj.pid_D/K_total;
            
            val = K_P.*P + K_I.*I + K_D.*D;
            
        end
        
        function resetRate(obj) % clear guiding history 
            
            obj.rate_RA = 0;
            obj.rate_DE = 0;
            
            obj.prev_dRA = 0;
            obj.prev_dDE = 0;
            
            obj.cam_pc.incoming.RA_rate_delta = 0;
            obj.cam_pc.incoming.DE_rate_delta = 0;
            
            obj.guiding_history.reset;
            obj.correct_history.reset(obj.integral_number);
            
        end
        
        function adjustPosition(obj, RA_deg, DE_deg) % not working yet
            
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
        
        function engineeringSlew(obj,Alt,Az) % alternative way to slew to Alt/Az
            
            if nargin<2 || isempty(Alt)
                Alt = obj.target_altitude;
            end
            
            if nargin<3 || isempty(Az)
                Az = obj.target_azimuth;
            end
            
            try 

                if Alt<obj.limit_alt
                    error('Input altitude is below limit of %f degress', obj.limit_alt);
                end
                
                obj.brake_bit = 0;
                
                obj.hndl.SlewToAltAzAsync(Az, Alt); % note that in ASCOM, SlewToAltAz expects Az and then Alt !!!
                
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
            
                obj.tracking = 0;

            catch ME
                obj.log.error(ME.getReport);
                obj.tracking = 0;
                rethrow(ME);
            end
            
            
        end
        
        function sync(obj) % this is still not working! 
            
            obj.hndl.SyncToCoordinates(obj.hndl.TargetRightAscension, obj.hndl.TargetDeclination);
            
        end
        
        function park(obj) % go to parking position 
            
            obj.hndl.Park;
            obj.tracking = 0;
            
        end
        
        function zenith_west(obj) % not yet implemented
            
        end
        
        function zenith_east(obj) % not yet implemented
            
        end
        
        function update(obj) % check connection to hardware and set status to 0 or 1
            
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
            
            obj.setup_timer; % the timer also updates GUI and arduino
            
            try % logger hearteat
                if ~obj.log.check_heartbeat
                    obj.log.heartbeat(300, obj); 
                end
            catch ME
                warning(ME.getReport);
            end
            
%             try % arduino
%                
%                 if obj.use_accelerometer
%                     
%                     obj.ard.update;
%                     
%                     if isempty(obj.ard) || obj.ard.status==0 || obj.checkArduinoTime==0
%                         obj.connectArduino;
%                     end
%                     
%                 end
%                 
%             catch ME
%                 warning(ME.getReport);
%             end
            
            % add additional tests?
            
            obj.status = 1;
            
        end
        
        function stop(obj) % stop moving the telescope, reset guiding, stop tracking, during normal operations (does not log this)
            
%             obj.log.input('stopping telescope');
            
            try 
                
                obj.brake_bit = 1;
                
                if ~isempty(obj.hndl) && obj.status
                    obj.hndl.AbortSlew;
                    obj.hndl.tracking = 0;
                end
                
                obj.resetRate;
                
%                 if ~isempty(obj.cam_pc)
% 
%                     if isfield(obj.cam_pc.outgoing, 'stop_camera') && obj.cam_pc.outgoing.stop_camera==0
%                         obj.log.input('Telescope stopped, sending camera stop command');
%                         disp(obj.log.report);
%                     end
% 
%                     obj.cam_pc.outgoing.stop_camera = 1;
%                     
%                 end
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function emergency_stop(obj) % stops slew, tracking and guiding and reports it on log
            
            obj.log.input('stopping telescope');
            
            disp('stopping telescope!');
            
            obj.stop;
            
        end
        
        function updateCamera(obj) % send details on object coordinates to cam_pc
            
            obj.cam_pc.outgoing.RA = obj.objRA;
            obj.cam_pc.outgoing.DEC = obj.objDec;            
            obj.cam_pc.outgoing.RA_DEG = obj.objRA_deg;
            obj.cam_pc.outgoing.DEC_DEG = obj.objDec_deg;
            obj.cam_pc.outgoing.TELRA = obj.telRA;
            obj.cam_pc.outgoing.TELDEC = obj.telDec;
            obj.cam_pc.outgoing.TELRA_DEG = obj.telRA_deg;
            obj.cam_pc.outgoing.TELDEC_DEG = obj.telDec_deg;
            
%             if isempty(obj.tracking) || obj.tracking==0 || (~isempty(obj.objRA_deg) && abs(obj.objRA_deg-obj.telRA_deg)>1) % if mount stops tracking or is 4 time-minutes away from target RA, stop the camera (e.g., when reaching limit)
%                 
%                 if ~isempty(obj.cam_pc) && isfield(obj.cam_pc.outgoing, 'stop_camera') && obj.cam_pc.outgoing.stop_camera==0
%                     obj.cam_pc.outgoing.stop_camera = 1;
%                     obj.log.input('Telescope not tracking, stopping camera');
%                     disp(obj.log.report); 
%                 end
%                 
%             else
%                 obj.cam_pc.outgoing.stop_camera = 0;
%             end
            
        end
        
        function stress_test(obj, varargin) % slew the telescope multiple times to random pointings to test for problems
        % Usage: stress_test(obj, varargin) 
        % Will choose random points in the sky, that are observable, using 
        % the Ephemeris.random() method to choose coordinates. 
        % It will slew and pause for a few seconds, then slew again. 
        % This continues for several targets, testing that the telescope
        % is working properly and does not run into any limits. 
        %
        % OPTIONAL ARGUMENTS:
        %   -number: how many random targets to choose. Default 10. 
        %   -pause: what length of delay between slews. Default 10 seconds. 
        % 
        % In addition to these, the varargin is also passed on to the 
        % head.Ephemeris.random() function, which looks at:
        %   -alt_limit: minimal altitude for targets. Default 25 degrees. 
        %   -south: minimal declination in degrees. Default -20. 
        %   -side: choose "east" or "west" or "south". Default "both". 
        %   -meridian: distance (on both sides) from meridian. Default 5 deg. 
        % 
        
            input = util.text.InputVars;
            input.input_var('number', 10, 'iterations');
            input.input_var('pause', 10, 'delay_time'); 
            input.scan_vars(varargin{:});
            
            for ii = 1:input.number
                
                obj.object.random(varargin);
                
                obj.slew('history', false); 
                
                pause(input.pause);
                
            end
            
        end
        
        function val = checkArduinoTime(obj) % check if arduino has been in contact in the last few minutes
            
            val = 1;
            
            if ~isempty(obj.log.time)
            
                dt = seconds(obj.log.time - obj.ard.time);
            
                if dt>600 % arduino has not updated in a few minutes
                    val = 0;
                    return;
                end
                
            end
            
        end
        
        function val = printout(obj) % quick summary of mount status
        
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
        
        function plot_rate(obj, ax) % show the recent guiding rates on a graph
            
            if ~isempty(obj.guiding_history.data)
                hold(ax, 'off');
                plot(ax, obj.guiding_history.data_ordered(:,1)*15, 'x');
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
        telDE;
        
    end
    
end

