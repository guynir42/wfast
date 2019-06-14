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
        
    end
    
    properties % objects
        
        hndl; % ASCOM object
        
        object@head.Ephemeris; % all position and time calculations done using Ephemeris class
        
        ard@obs.sens.ScopeAssistant; % connect to accelerometer and ultrasonic sensor
        
        log@util.sys.Logger; % keep track of all commands and errors
        
    end
    
    properties % inputs/outputs
        
        status = 0; % if connected and responsive, set this to 1
        
    end
    
    properties % switches/controls
        
        limit_alt = 15; % degrees above horizon where telescope is not allowed to go
        
        use_accelerometer = 0; % make constant checks for altitude outside of the mounts own sensors
        use_ultrasonic = 0; % make constant checks that there is nothing in front of the telescope
        
        step_arcsec = 5; % not used yet
        
        brake_bit = 1; % when slewing, set this to 0. Interrupts may set it back to 1 to stop the slew. 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        % these are parameters of the target object
        % they are kept in the Ephemeris class
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
        
    end
    
    properties(Hidden=true)
       
        version = 1.02;
        
    end
    
    methods % constructor
        
        function obj = ASA(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.mount.ASA')
                if obj.debug_bit, fprintf('ASA copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('ASA constructor v%4.2f\n', obj.version); end
            
            end
            
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
                
                obj.hndl = actxserver('AstrooptikServer.Telescope');
            
                obj.hndl.SiteLatitude  = obj.object.latitude;
                obj.hndl.SiteLongitude = obj.object.longitude;
                
                obj.hndl.Connected = 1;
            
                obj.hndl.MotorOn;
                
                if obj.use_accelerometer && isempty(obj.ard) 
                    
                    try
                        obj.ard = obs.sens.ScopeAssistant;
                        obj.ard.telescope = obj;
                    catch ME
                        obj.use_accelerometer = 0;
                        warning(ME.getReport);
                    end

                end
                
                if obj.use_accelerometer
                    obj.ard.connect;
                end
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function loadServer(obj)
            
%             system('C:\Program Files (x86)\Autoslew\AstroOptikServer.exe &'); % call the system command to load the server outside of matlab 
            system('D:\matlab\wfast\+obs\+mount\launch_server.bat &'); % call the batch file to load the server then exit the cmd
            
            tic;
            
            for ii = 1:100
                
                [~, str] = system('tasklist /fi "imagename eq AstroOptikServer.exe"'); % check 
                
                str = strip(str);
                idx = strfind(str, 'AstroOptikServer.exe');
                
                if isempty(idx)
                    continue; % can't find the server on the tasklist
                else
                    num = util.text.extract_numbers(str(idx:end)); % found the server, now check the memory usage
                    if num{1}(3)>30
                        return; % server has loaded and is connected to telescope
                    end
                    % if this last condition fails, repeatedly, it means
                    % the server is stuck on a dialog... 
                end
                
                pause(0.05);
            
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
        
    end
    
    methods % getters
        
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
        
    end
    
    methods % setters
        
        function set.objRA(obj, val)
            
            obj.object.RA = val; % note this is given in HOURS!
            obj.object.update;
            
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
            
            try
                obj.hndl.TargetDeclination = obj.object.Dec_deg;
            end
            
        end
        
        function set.objDEC_deg(obj, val)
            
            obj.object.Dec_deg = val;
            obj.object.update;
            
            try
                obj.hndl.TargetDeclination = obj.object.Dec_deg;
            end
            
        end
        
        function set.DE_target_deg(obj, val)
            
            obj.object.Dec_deg = val;
            obj.object.update;
            
            try
                obj.hndl.TargetDeclination = obj.object.Dec_deg;
            end
            
        end
        
        function set.tracking(obj, val)
            
            try 
                obj.hndl.Tracking = val;
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
    end
    
    methods % calculations / commands
        
        function inputTarget(obj, varargin)
            
            obj.object.input(varargin{:}); % Ephemeris now uses Eran's name resolver
            
        end
        
        function val = check_before_slew(obj)
            
            val = 0;
            
            obj.object.update;
            
            if obj.object.Alt_deg<obj.limit_alt
                disp(['Target alt (' num2str(obj.object.Alt_deg) ') is below alt limit (' num2str(obj.limit_alt) ')']);
                return;
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
        
        function slew(obj, varargin)
            
            if isempty(obj.object.RA) || isempty(obj.object.DE) 
                error('Please provide a target with viable RA/DE');
            end
            
            obj.log.input(sprintf('Slewing to target. RA= %s | DE= %s | ALT= %4.2f', obj.object.RA, obj.object.DE, obj.object.ALT_deg));
            
            try 
                
                if ~obj.check_before_slew
                    error('Prechecks failed, aborting slew');
                end
                
                ra_hours_Jnow = obj.object.RA_deg_now./15; % convert to hours! 
                dec_deg_Jnow = obj.object.Dec_deg_now;
                
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
            
                obj.tracking = 1;
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
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
        
        function park(obj)
            
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
            
            try
               
                if isempty(obj.ard) && obj.use_accelerometer
                    obj.ard = obs.sens.ScopeAssistant;
                    obj.ard.telescope = obj;
                end
                
                if ~isempty(obj.ard) && obj.use_accelerometer
                    obj.ard.update;
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
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
        function emergency_stop(obj)
            
            obj.log.input('stopping telescope');
            
            try 
                
                obj.brake_bit = 1;
            
                obj.hndl.AbortSlew;
                
            catch ME
                obj.log.error(ME.getReport);
                rethrow(ME);
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
    properties(Transient=true, Hidden=true, Dependent=true) % for backward compatibility
        
        DE_target_deg;
        DE_target;
        telDE_deg;
        telDE;
        
    end
    
end

