classdef ASA < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        hndl;
        
    end
    
    properties % inputs/outputs
        
    end
    
    properties % switches/controls
        
        step = 1; % degree... 
        
        debug_bit = 1;
        
        RA_target; % in degrees        
        DE_target; % in degrees
        
    end
    
    properties(Dependent=true)
        
        RA_target_str;
        DE_target_str;
        
        RA_str;
        DE_str;
        LST_str;
        
        RA;
        DE;
        LST;
        HA;
        
        ALT;
        AZ;
        
        HA_target;
        ALT_target;
        AZ_target;
        
        % can we get motor status from hndl??
        tracking;
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = ASA(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.mount.ASA')
                if obj.debug_bit, fprintf('ASA copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('ASA constructor v%4.2f\n', obj.version); end
            
            end
            
            obj.connect;
            
        end
        
        function connect(obj)
            
            try 
                obj.hndl = actxserver('AstrooptikServer.Telescope');
            
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
        
    end
    
    methods % getters
        
        function val = get.RA_target_str(obj)
            
            val = head.Ephemeris.rad2ra(deg2rad(obj.RA_target));
            
        end
        
        function val = get.DE_target_str(obj)
            
            val = head.Ephemeris.rad2dec(deg2rad(obj.DE_target));
            
        end
        
        function val = get.RA(obj)
            
            val = obj.hndl.RightAscension/24*360;
            
        end
        
        function val = get.RA_str(obj)
            
            val = head.Ephemeris.rad2ra(deg2rad(obj.RA));
            
        end
        
        function val = get.DE(obj)
            
            val = obj.hndl.Declination;
            
        end
        
        function val = get.DE_str(obj)
            
            val = head.Ephemeris.rad2dec(deg2rad(obj.DE));
            
        end
        
        function val = get.ALT(obj)
            
            val = obj.hndl.Altitude;
            
        end
        
        function val = get.AZ(obj)
            
            val = obj.hndl.Azimuth;
            
        end
        
        function val = get.LST(obj)
            
            val = obj.hndl.SiderealTime/24*360; % in degrees...
            
        end
        
        function val = get.HA(obj)
            
            val = obj.LST - obj.RA;
            
        end
        
        function val = get.LST_str(obj)
            
            val = head.Ephemeris.rad2ra(deg2rad(obj.LST));
            
        end
        
        function val = get.HA_target(obj)
            
            val = obj.LST - obj.RA_target;
            
        end
        
        function val = get.ALT_target(obj)
            
            
            
        end
        
        function val = get.AZ_target(obj)
            
        end
        
        function val = latitutde(obj)
            
            val = obj.hndl.SiteLatitude;
            
        end
        
        function val = longitude(obj)
           
            val = obj.hndl.SiteLongitude;
            
        end
        
    end
    
    methods % setters
        
        function set.RA_target(obj, val)
            
            obj.hndl.TargetRightAscension = val*24/360;
            
            obj.RA_target = val; % should we make this dependent?
            
        end
        
        function set.DE_target(obj, val)
            
            obj.hndl.TargetDeclination = val;
            
            obj.DE_target = val; % should we make this dependent?
            
        end
        
    end
    
    methods % calculations
        
        function input(obj, varargin)
            
            if ~isempty(which('celestial.coo.coo_resolver', 'function'))
                [obj.RA_target, obj.DE_target] = celestial.coo.coo_resolver(varargin{:}, 'OutUnits', 'deg');
            else
                error('no name resolver has been found... try adding MAAT to the path.');
            end
            
        end
        
        function goto(obj, varargin)
            
            
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
    methods(Static=true)
        
        
        
    end
    
end

