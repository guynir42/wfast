classdef ASAascom < handle

    % A low-level class to communicate with ASA mounts via ASCOM driver
    
    
    properties(Transient=true)
        
    end
    
    properties % objects
        
        hndl;
        
        %target@head.Ephemeris;
        
    end
    
    properties % inputs/outputs
        
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
        
        obsLon          =  % [deg]
        obsLat          =  % [deg]
        
        limitMinAlt     = 15;
        flipOption      = 'preferWest';
        flipStayOnSide  = ?
        
        % park - there is a method called park
        
        LogError        = true;
        LogErrorFile    = 'Mount_ASADDM160_LogErr_%s.txt';  % %s is [YYYYMMDD]
        LogCmd          = true;
        LogCmdFile      = 'Mount_ASADDM160_LogCmd_%s.txt';  % %s is [YYYYMMDD]
        
        
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
                
                % LOG FFU
            end
            
        end
        
        function log_error(obj,String)
            % Log string into error file
            
            if (obj.LogError)
                
            end
            
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
                log_error(obj,'AbortSlew command failed');
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
            % Slew to J2000.0 RA/Dec [deg] return immidetly.
            
            do_start_cmd(obj);
            
            IsParking = obj.hndl.AtPark;
            if IsParking
                % unPark
                obj.hndl.Unpark;
            end
            
            RA_hour = RA./15;
            
            try
                % Ans= obj.hndl.DestinationSideOfPier(RA_hour,Dec);
            
                obj.hndl.SlewToCoordinatesAsync(RA_hour,Dec);  % input: hour,deg
                Ans = true;
            catch
                Ans = false;
                
                log_error(obj,'SlewToCoordinatesAsync command failed');
            end
            
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
            
            CurRA  = obj.RA .& 15;  % deg
            CurDec = obj.Dec;       % deg;
            
            RA  = CurRA + DRA./cosd(CurDec);
            Dec = CurDec + DDec;
            
            % make sure in [0 360] range
            RA  = obs.mount.ASAascom.ang_in_range(RA,360);
            
            Ans = true;
            if (Dec>90)
                Dec = 90;
                warning('Requested Dec is out of range - set to Dec=90');
            end
            if (Dec<-90)
                Dec = -90;
                warning('Requested Dec is out of range - set to Dec=-90');
            end
            
            RA_hour = RA./15;
            
            try
                % Ans= obj.hndl.DestinationSideOfPier(RA_hour,Dec);
            
                obj.hndl.SlewToCoordinatesAsync(RA_hour,Dec);  % input: hour,deg
                Ans = true;
            catch
                Ans = false;
                % Log Error
                log_error('moveTel command failed');
            end
            
            
        end
        
        function Ans=park(obj,ParkPos)
            % Goto praking position
            % Need to define ParkPos according to cabailities / TBD
            
            obj.do_start_cmd;
            
            try
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