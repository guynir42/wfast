classdef ASAascom < handle

    % A low-level class to communicate with ASA mounts via ASCOM driver
    
    
    properties(Transient=true)
        
    end
    
    properties % objects
        
        hndl;
        log@util.sys.Logger;
        %target@head.Ephemeris;
        
    end
    
    properties % inputs/outputs
        
        targetRA
        targetDec
        targetHA
        targetAz
        targetAlt
        
        obsLon          = 34.9; % [deg]
        obsLat          = 30.5; % [deg]
        
        limitMinAlt     = 15;
        flipOption      = 'preferWest';
        flipStayOnSide  = [];
        
        % park - there is a method called park
                
        whileMoveAsync  = true;
        
        LogError        = true;
%         LogErrorFile    = 'Mount_ASADDM160_LogErr_%s.txt';  % %s is [YYYYMMDD]
        LogCmd          = true;
%         LogCmdFile      = 'Mount_ASADDM160_LogCmd_%s.txt';  % %s is [YYYYMMDD]
        
    end
    
    properties % switches/controls
        
        step = 1; % degree... 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        status
        
        telRA     % telescope current RA [deg]
        telDec    % telescope current Dec [deg]
        telAz    
        telAlt
        telHA
        telLST
        telTime
        
        
        
        tracking  % tracking true/false
        rate      % tracking rate
        
  
        
    end
    
    properties(Hidden=true)
       
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = ASAascom(varargin)
            % ASAascom constructor
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.mount.ASAascom')
                if obj.debug_bit
                    fprintf('ASA copy-constructor v%4.2f\n', obj.version);
                end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit
                    fprintf('ASA constructor v%4.2f\n', obj.version);
                end
            
            end
            
            obj.connect;
            
            obj.update;
            
            % set default tracking rate
            
            % set default limits
            
            % set default observatory coordinates
            
        end
        
        function connect(obj)
            % Connect to the ASA mount using the ASCOM driver
            
            obj.log.input('Connecting to mount.');
            
            try 
                
                obj.hndl = actxserver('AstrooptikServer.Telescope');
            
                obj.hndl.SiteLatitude  = obj.obsLat;
                obj.hndl.SiteLongitude = obj.obsLon;
                
                obj.hndl.Connected = 1;
            
                obj.hndl.MotorOn;
               
            catch ME
                % should we do anything else?
                warning(ME.getReport);
                
                % write log Error
                log_error(obj,'Error while trying to connect to ASA mount');
                
            end
            
        end
        
        function disconnect(obj)
            
        end
        
    end

    
    
    methods % reset/clear
        
        function Ans=motorOn(obj)
            % set MotorOn, return false if error
            % 
            
            do_start_cmd(obj);
            
            try
                obj.hndl.MotorOn;
                Ans = true;
            catch
                Ans=false;
                % Log error
                log_error(obj,'MotorOn command failed');
                
            end
            
        end
        
        function motorOff(obj)
            % set MotorOff
            
            do_start_cmd(obj);
            
            try
                obj.hndl.MotorOff;
                Ans = true;
            catch
                Ans=false;
                % Log error
                log_error(obj,'MotorOff command failed');
            end
            
        end
            
        
        function reset(obj)
            % ???
                        
        end
        
        function update(obj)
            % ???
            
            
        end
        
        
    end
    
    methods % do start/end
        function do_start_cmd(obj)
            % This command is running in the begining of each command
            
            % Log command in log file
            if (obj.LogCmd)
                % Log command into log file
                
                DB = dbstack; % get from dbstack the name of the caller function
                
                obj.log.input(DB(2).name); % write the name of the calling function, outside of the "do_start_cmd"
                % LOG FFU
            end
            
        end
        
        function log_error(obj,String)
            % Log string into error file
            
            if (obj.LogError)
                obj.log.error(String); 
            end
            
        end
        
        function verify_equatorial_coo(obj)
            % convert coordinate types and verify validity/observability
            % Description: If coordinates are invalid than the targetRA/Dec
            %              will be populated with NaN.
            
            RAD           = 180./pi;
            SEC_IN_DAY    = 86400;
            MAX_LST_DIFF  = 10;     % [sec]
            
            obj.do_start_cmd;
            
            RA  = obj.targetRA;
            Dec = obj.targetDec;
            
            if isempty(RA)
                % RA is not provided
                obj.log_error('RA is empty');
                RA = NaN;
            else
                if ischar(RA)
                    % RA is sexagesimal / convert to deg
                    RA = celestial.coo.convertdms(RA,'SH','d');
                else
                    % RA is something else
                    % FFU
                    RA = NaN;
                end
            end
            if isempty(Dec)
                % Dec is not provided
                obj.log_error('Dec is empty');
                Dec = NaN;
            else
                if ischar(Dec)
                    % Dec is sexagesimal / convert to deg
                    Dec = celestial.coo.convertdms(Dec,'SD','d');
                else
                    % Dec is something else
                    % FFU
                    Dec = NaN;
                end
            end
            
            
            % verify coordinates validity
            if RA>=0 && RA<=360
                % RA range is valid
            else
                % invalid range
                obj.log_error(sprintf('RA is not in valid range: %f',RA));
                RA = NaN;
            end
            if Dec>=-90 && RA<=90
                % Dec range is valid
            else
                % invalid range
                obj.log_error(sprintf('Dec is not in valid range: %f',Dec));
                Dec = NaN;
            end
            
            % Check if RA/Dec is observable
            JD            = celestial.time.julday;
            LST_computer  = celestial.time.lst(JD,obj.obsLon./RAD,'m');  % [frac of day]
            LST_telescope = obj.telLST;   % [frac of day]
            % compare telescope and computer LST
            if abs(LST_computer - LST_telescope)>(MAX_LST_DIFF./SEC_IN_DAY)
                % There is a big offset between computer LST and telescope
                % LST
                % set targetRA/Dec to NaN and log error
                RA  = NaN;
                Dec = NaN;
                obj.log_error(sprintf('Big offset between computer and telescope LST : %f seconds',LST_computer - LST_telescope));
            else
                % Check that target is above the horizon
                HorizCoo = celestial.coo.horiz_coo([RA Dec]./RAD, JD, [obj.obsLon obj.obsLat]./RAD, 'h');
                Az       = HorizCoo(1).*RAD;  % [deg]
                Alt      = HorizCoo(2).*RAD;  % [deg]
                
                if (Alt<obj.limitMinAlt)
                    % requested Alt is below telescope limit
                    % set targetRA/Dec to NaN and log error
                    obj.log_error(sprintf('Requested RA/Dec are below altitude limit : %f %f %f',RA,Dec,Alt));
                else
                    % any additional verification come here
                    
                    %--- Coordinates verified ---
                end
            end
            
            % populate the target RA/Dec in deg
            obj.targetRA  = RA;
            obj.targetDec = Dec;
            
            
        end
        
        function do_before_move(obj)
            % Aggregation of commands to perform before telescope motion
            
            obj.do_start_cmd;
            
            % check that telescope is not parking, if so Unpark
            IsParking = obj.hndl.AtPark;
            if IsParking
                % unPark
                try
                    obj.hndl.Unpark;
                catch
                    % Log error
                    obj.log_error('Unpark command while slewCoo failed');
                end
            end
            
        end
        
        function while_slewing(obj)
            % Upper aggregation of commands to perform while telescope is moving
            % Description: This is calling the do_while_slewing set of
            %              commands, wither using timer or a while loop.
            
            if obj.whileMoveAsync
                % Aynchronous mode (using timer)
                
                T          = timer;
                T.Period   = 1;
                T.ExecuationMode = 'fixedDelay';
                T.TimerFcn = @obj.timer_slewing;
                
                
            else
                % Synchrounous mode (while loop)
                
                while obj.isSlewing
                    % do as long as telescope is slewing
                   
                    obj.do_while_slewing;
                    
                end
                
                obj.do_after_slewing;
                
            end
            
        end
        
        function timer_slewing(obj,T)
            % Callback function for slewing timer
            % Input  : - Mount object
            %          - Timer object
            
            if obj.isSlewing
                obj.do_while_slewing;
            else
                stop(T);
                clear(T);
                
                obj.do_after_slewing;
            end
            
        end
        
        function do_while_slewing(obj)
             % Lower aggregation of commands to perform while telescope is moving
             % Description: Here are the actual commands to perform while
             %              the telescope is moving.
             
             
             % 1. calculate current altitude and verify it is above the
             % horizon
             Alt = obj.telAlt;
             if (Alt<obj.limitMinAlt)
                 % Problem detected - telescope is below horizon limit
                 obj.AbortSlew
                 
                 obj.log_error(sprintf('Slew aborted because telescope altitude %f is below limit',Alt));
             end
             
             % Check the status of the accelerometer
             
             
        end
        
        
        function do_after_slewing(obj)
            % Commands to do after slewing
            
            
        end
        
    end
    
    methods % getters
        
         function val=get.status(obj)
            % get mount status
            % Package: obs.mount
            % Description: Return true if mount is responsive
            
            obj.do_start_cmd;
            
            try
                val = obj.hndl.RightAscension;
            catch
                val = NaN; 
            end
            
            if isnan(Val)
                % telescope status is not ok
                val = false;
            else
                % telescope is responding assume status is ok
                val = true;
            end
            
            
        end
        
         function val = get.telRA(obj)
            % get RA of telescope in deg
            % Is this J2000 or equinox of date?
            % Output : - Mount R.A. J2000(?) in deg.
            %            Return NaN if error.
            
            obj.do_start_cmd;
            
            
            try
                val = obj.hndl.RightAscension;
            catch
                val = NaN;
                
                if (obj.LogError)
                    % Log Error
                    
                    % Log FFU
                end
                
            end
            
            % Assuming the output is in hours - convert to deg
            val = val.*15;
            
         end
        
         function val = get.telDec(obj)
            % get Dec of telescope in deg
            
            % copy from RA: Declination
            
            
         end
        
         function val = get.telAz(obj)
            % get Az of telescope in deg
            
            % copy from RA: Azimuth
             
            
         end
         
         function val = get.telAlt(obj)
            % get Alt of telescope in deg
            
            % copy from RA: Altitude
            
         end
         
         function val = get.telHA(obj)
            % get HA of telescope in deg
            
            val = 
            
         end
         
         function val = get.telLST(obj)
            % get LST of telescope in fraction of day
            
            val = 
            
         end
         
         function val = get.telTime(obj)
            % get Time of telescope in JD
            
            val = 
            
         end
         
         function val = get.tracking(obj)
            % get tracking on/off
            
            val = 
            
         end
        
         function val = get.rate(obj)
            % get tracking rate
            
            val = 
            
         end
         
         function val = get.obsLon(obj)
            % get obs Longitude in [deg]
            
            val = 
            
         end
         
         function val = get.obsLat(obj)
            % get obs Latitude in [deg]
            
            val = 
            
         end
         
         function val = get.targetRA(obj)
             % get tragetRA 
            
             if isempty(targetRA)
             
         end
         
         
    end
    
    methods % setters
        function Ans=set.tracking(obj, val)
            % set tracking [true/false]
            % Return true if sucssful
            
            
        end
        
        function Ans=set.rate(obj, val)
            % set tracking rate speed: <val ["/s]>, ['sidereal'], 'lunar',
            % 'solar'
            % Return true if sucssful
            
        end
        
    end
    
    methods % critical commands
        
        
        function Ans=abortSlew(obj)
            % stop/abort telescope motion
            % Package: obs.mount
            % Description: Abort telescope slew. This command is effective
            %              only after call to SlewToTargetAsync, SlewToCoordinatesAsync,
            %              SlewToAltAzAsync, MoveAxis.
            %              Tracking is returned to its pre-slew state.
            %              See p.2 in The ITtelescopeV3 Methods
            
            obj.do_start_cmd;
            
            try
                obj.hndl.AbortSlew;
                Ans = true;
            catch
                Ans = false;
                % Log Error
                obj.log_error('AbortSlew command failed');
            end
            
            
        end
        
        function Ans=ShutDown
            % Shut down the mount
            
            % ShutDown
            
            
        end
        
    end
    
    methods (Static)
        function Long=ang_in_range(Long,Period)
            % Convert an angle to the 0 to 360 range
            % Description: Convert an angle to the range 0 to 360.
            % Input  : - Matrix of angles.
            %          - Period of angles. Default is 360.
            % Output : - Angle in the allowed range.

            Long = (Long./Period - floor(Long./Period)).*Period;

            
        end
        
    end
    
    methods % telescope motion
        
        function [Val,Res]=isTracking(obj)
            % Return true if telescope is tracking
            % Output  : - True/false for telescope tracking
            %           - Error/sucess message
            
            
            try
                Val = obj.hndl.Tracking;
                Res = 'valid response';
            catch
                Val = NaN;
                Res = 'Exception while calling isTracking';
                
                % log error
                log_error(obj,sprintf('isTracking command failed'));
            end
                
        end
        
        function Ans=slewCoo(obj,RA,Dec)
            % Slew to J2000.0 RA/Dec [deg] return immidetly [unless]
            
            do_start_cmd(obj);
            
            if (nargin>=3)
                % user supplied coordinates in arguments
                obj.targetRA  = RA;
                obj.targetDec = Dec;
            else
                % get Coorddinates from targetRA, targetDec
                RA  = obj.targetRA;
                Dec = obj.targetDec;
            end
            
            % convert equatorial coordinates in targetRA/Dec to deg
            % and verify validity and observability of coordinates
            
            obj.verify_equatorial_coo;
            RA  = obj.targetRA;    % deg
            Dec = ibj.targetDec;   % deg
            
            
            % check that telscope is not AtPark, if so Unpark
            obj.do_before_move;
            
            
            RA_hour = RA./15;
            
            if isnan(RA) || isnan(Dec)
                % RA/Dec are NaN
                % don't move telescope
                
            else
                % move telescope
                try
                    % Ans= obj.hndl.DestinationSideOfPier(RA_hour,Dec);

                    obj.hndl.SlewToCoordinatesAsync(RA_hour,Dec);  % input: hour,deg
                    Ans = true;
                catch
                    Ans = false;

                    log_error(obj,'SlewToCoordinatesAsync command failed');
                end
            end
            
            obj.while_slewing;
            
        end
        
        function Ans=slewAzAlt(obj,Az,Alt)
            % Slew to Az, Alt [deg]
            % Return false if operation is not allowed
            
            % obj.hndl.SlewToAltAzAsync
        end
        
        function Ans=slewHADec(obj,HA,Dec)
            % Slew to HA, Dec [deg]
            % Return false if operation is not allowed
            
        end
        
        function Ans=moveTel(obj,DRA,DDec)
            % move telescope in RA and Dec [deg]
            % Input  : - Mount object
            %          - Delta RA to move [deg of arc]
            %          - Delta Dec to move [deg of arc]
            % Output : - false if command move command failed.
            
            do_start_cmd;
            
            CurRA  = obj.RA .* 15;  % deg
            CurDec = obj.Dec;       % deg;
            
            RA  = CurRA + DRA./cosd(CurDec);
            Dec = CurDec + DDec;
            
            % make sure in [0 360] range
            RA  = mod(RA,360);
            
            Ans = true;
            if (Dec>90)
                Dec = 90;
                warning('Requested Dec is out of range - set to Dec=90');
            end
            if (Dec<-90)
                Dec = -90;
                warning('Requested Dec is out of range - set to Dec=-90');
            end
            
            
            Ans = obj.slewCoo(RA,Dec)
            
        end
        
        function Ans=park(obj,ParkPos)
            % Goto praking position
            % Need to define ParkPos according to cabailities / TBD
            
            obj.do_start_cmd;
            
            try
                % shut trcking off
                
                % SlewToAltAzAsync
                
                
                
                
                obj.hndl.park;
                Ans = true;
            catch
                Ans = false;
                % Log Error
                log_error(obj,'park command failed');
            end
            
            % verify that AtPark is true
            AtPark = obj.hndl.AtPark;
            if (~AtPark)
                % error
                log_error(obj,'park command return, but telescope is not AtPark');
                Ans = false;
            end
            
            % Verify that telescope is at park position by reading
            % its Az/Alt
            try
                Az  = obj.hndl.Azimuth;
                Alt = obj.hndl.Altitude;
                
                % FFU
                
            catch
                % Log Error
                log_error('Error while trying to read Az/Alt to verify parking position');
            end
            
            
            
        end
        
        function [Ans]=isParking(obj)
            % True if telescope is parking
            % Output : - True if telescope parking, NaN if unknown.
            
            obj.do_start_cmd;
            
            try
                Ans = obj.hndl.AtPark;
            catch
                Ans = NaN;
                % log error
                log_error('isParking command failed');
            end
            
        end
        
        function Ans=setPark(obj)
            % set the Telescope parking position to its current position
            
        end
        
        
        function Ans=isSlewing(obj)
            % Retrun true if telescope is slewing, NaN if error.
           
            obj.do_start_cmd;
            
            try
                Ans = obj.hndl.Slewing;
            catch
                Ans = NaN;
                % log error
                log_error('isSlewing command failed');
            end
            
        end
        
        function Ans=sideOfPier(obj,RA,Dec)
            % Predict side of pier for German equatorial mount
            % Input  : - ASAascom object
            %          - RA [deg]
            %          - Dec [deg]
            % Output : - Side of pier ?, NaN if error.
            
            do_start_cmd(obj);
            
            RA_hour = RA./15;
            
            try
                Ans = obj.hndl.DestinationSideOfPier(RA_hour,Dec);  % input hours,deg
            
                % Check answer validity
                % FFU?????
                
            catch
                Ans = NaN;
                % Log File
                log_error(obj,'sideOfPier command failed');
            end
            
            
        end
        
        function [Ans,Res]=findHome(obj)
            % find telescope home (synchronous)
            % Output  : - true/false
            %           - Error/sucess Message
            
            do_start_cmd;
            
            IsPark = obj.isParking;
            if isnan(IsPark)
                % can't find home - isParking unknown
                Ans = false;
                Res = 'isParking returne NaN - cant find home';
            else
                if (IsPark)
                    % Can't find home while telescope is parking
                    Ans = false;
                    Res = 'unpark telescope - cant find home while telescope is parking';
                else
                    try
                        obj.hndl.FindHome;
                        % synchronous
                        Ans = true;
                        Res = 'sucess';
                    catch
                        Ans = false;
                        Res = 'FindHome exception';
                        
                        % Log Error
                        log_error(obj,'findHome command failed');
                    end
                end
            end
        end
        
    end
    
end