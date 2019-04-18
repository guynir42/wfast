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
            
            val = 
            
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
        function Ans=stop(obj)
            % stop telescope motion
            % Return true if stop sucessfully
            
        end
        
        
    end
    
    methods % telescope motion
        function Ans=slewCoo(obj,RA,Dec)
            % Slew to J2000.0 RA/Dec [deg]
            % Return false if operation is not allowed
            
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
        
        function Ans=isSlewing(obj)
            % Retrun true if telescope is slewing, false if on target
            
        end
    end
    
end