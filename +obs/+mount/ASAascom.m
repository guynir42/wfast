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
        
        park % ?
        
        
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
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.mount.ASA')
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
            
                obj.hndl.SiteLatitude  = obj.target.obsLat;
                obj.hndl.SiteLongitude = obj.target.obsLon;
                
                obj.hndl.Connected = 1;
            
                obj.hndl.MotorOn;
               
            catch ME
                % should we do anything else?
                warning(ME.getReport);
            end
            
        end
        
        function disconnect(obj)
            
        end
        
    end

    
    
    methods % reset/clear
        
        function reset(obj)
            % ???
                        
        end
        
        function update(obj)
            % ???
            
            
        end
        
        
    end
    
    
    methods % getters
         function val = get.telRA(obj)
            % get RA of telescope in deg
            % Is this J2000 or equinox of date?
            % Output : - Mount R.A. J2000(?) in deg.
            %            Return NaN if error.
            
            try
                val = obj.hndl.RightAscension;
            catch
                val = NaN;
                
                % log error
                
            end
            
            % Assuming the output is in hours - convert to deg
            val = val.*15;
            
         end
        
         function val = get.telDec(obj)
            % get Dec of telescope in deg
            
            val = 
            
         end
        
         function val = get.telAz(obj)
            % get Az of telescope in deg
            
            val = 
            
         end
         
         function val = get.telAlt(obj)
            % get Alt of telescope in deg
            
            val = 
            
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
            % Return tru if sucssful
            
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
            
            obj.hndl.AbortSlew;
            
            
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
                
                % log file
            end
                
        end
        
        function slewCoo(obj,RA,Dec)
            % Slew to J2000.0 RA/Dec [deg] return immidetly.
            
            RA_hour = RA./15;
            obj.hndl.SlewToCoordinatesAsync(RA_hour,Dec);  % input: hour,deg
            
            
            
        end
        
        function Ans=slewAzAlt(obj,Az,Alt)
            % Slew to Az, Alt [deg]
            % Return false if operation is not allowed
        end
        
        function Ans=slewHADec(obj,HA,Dec)
            % Slew to HA, Dec [deg]
            % Return false if operation is not allowed
            
        end
        
        function Ans=park(obj,ParkPos)
            % Goto praking position
            % Need to define ParkPos according to cabailities / TBD
            
            
            
        end
        
        function Ans=isParking(obj)
            % True if telescope is parking
            % Output : - True if telescope parking, NaN if error.
            
            try
                Ans = obj.hndl.AtPark;
            catch
                Ans = NaN;
                
                % log error
            end
            
        end
        
        function Ans=setPark(obj)
            % set the Telescope parking position to its current position
            
        end
        
        
        function Ans=isSlewing(obj)
            % Retrun true if telescope is slewing, false if on target
            
        end
        
        function Ans=sideOfPier(obj,RA,Dec)
            % Predict side of pier for German equatorial mount
            % Input  : - ASAascom object
            %          - RA [deg]
            %          - Dec [deg]
            % Output : - Side of pier ?
            
            RA_hour = RA./15;
            Ans = obj.hndl.DestinationSideOfPier(RA_hour,Dec);  % input hours,deg
            
            % Check answer validity
            
            
        end
        
        function [Ans,Res]=findHome(obj)
            % find telescope home (synchronous)
            % Output  : - true/false
            %           - Error/sucess Message
            
            
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
                    end
                end
            end
        end
        
    end
    
end