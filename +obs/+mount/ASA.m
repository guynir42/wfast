classdef ASA < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        hndl;
        
        target@head.Ephemeris;
        
    end
    
    properties % inputs/outputs
        
    end
    
    properties % switches/controls
        
        step = 1; % degree... 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        RA_target; % in degrees        
        DE_target; % in degrees
        
        RA_target_str;
        DE_target_str;
        
        HA_target;
        ALT_target;
        AZ_target;
        
        RA;
        DE;        
        
        RA_str;
        DE_str;
        
        HA;
        ALT;
        AZ;
        
        LST;
        LST_str;
        
        % can we get motor status from hndl??
        tracking;
        
        rate;
        
    end
    
    properties(Hidden=true)
       
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = ASA(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'obs.mount.ASA')
                if obj.debug_bit, fprintf('ASA copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('ASA constructor v%4.2f\n', obj.version); end
            
            end
            
            obj.target = head.Ephemeris;
            
            obj.connect;
            
            obj.update;
            
        end
        
        function connect(obj)
            
            try 
                
                obj.hndl = actxserver('AstrooptikServer.Telescope');
            
                obj.hndl.SiteLatitude = obj.target.latitude;
                obj.hndl.SiteLongitude = obj.target.longitude;
                
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
            
            obj.target.reset;
            
        end
        
    end
    
    methods % getters
        
        function val = get.RA_target(obj)
            
            val = rad2deg(obj.target.RA_rad);
            
        end
        
        function val = get.DE_target(obj)
            
            val = rad2deg(obj.target.DE_rad);
            
        end
        
        function val = get.RA_target_str(obj)
            
            val = obj.target.RA; 
            
        end
        
        function val = get.DE_target_str(obj)
            
            val = obj.target.DE;
            
        end
        
        function val = get.HA_target(obj)
            
            val = obj.target.HA;
            
        end
        
        function val = get.ALT_target(obj)
            
            val = obj.target.ALT;
            
        end
        
        function val = get.AZ_target(obj)
            
            val = obj.target.AZ;
            
        end
        
        function val = get.RA(obj)
            
            val = obj.hndl.RightAscension/24*360;
            
        end
        
        function val = get.DE(obj)
            
            val = obj.hndl.Declination;
            
        end
        
        function val = get.RA_str(obj)
            
            val = head.Ephemeris.rad2ra(deg2rad(obj.RA));
            
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
        
        function val = get.HA(obj)
            
            val = obj.LST - obj.RA;
            
        end
        
        function val = get.LST(obj)
            
            val = obj.hndl.SiderealTime/24*360; % in degrees...
            
        end
        
        function val = get.LST_str(obj)
            
            val = head.Ephemeris.rad2ra(deg2rad(obj.LST));
            
        end
        
        function val = latitutde(obj)
            
%             val = obj.hndl.SiteLatitude;
            val = obj.target.latitude;
            
        end
        
        function val = longitude(obj)
           
%             val = obj.hndl.SiteLongitude;
            val = obj.target.longitude;

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
        
        function update(obj)
            
            obj.target.update;
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
    methods(Static=true)
        
        
        
    end
    
end

