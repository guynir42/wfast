classdef (CaseInsensitiveProperties, TruncatedProperties) Ephemeris < handle
% This class keeps track of a specific time and sky coordinate. 
% Typically each observation would have a specific Ephemeris. 
% It uses these three data (RA,DEC and time) to calculate all sorts of 
% additional information such as LST, altitude, moon/sun positions, 
% galactic and ecliptic coordinates and so on. 
%
% It can accept new targets using the input() function, that can get either
% a pair of arguments (RA/DEC) or a name of a star which is resolved using
% Eran's interface to SIMBAD. In this case give an empty DEC. Also, you can
% sepcify the time of the observation as the third argument, or leave the 
% time that the object was set to before. Use "now" to update to current time. 
% IMPORTANT: the RA must be input as hours, in numeric value or sexagesimal. 
%
% If you use other ways to update the time/coordinates of the object, make 
% sure to call update() or updateSecondaryCoords() to get the Jnow, ecliptic, 
% galactic, and sun/moon positions updated to the new data. The difference
% is that update() also moves the time to current time. 
%
% More useful functions include timeTravelHours() to move forward (or back)
% in time to check observational conditions then. Use calcObsTimeMinutes()
% to figure out how many minutes can the field be observed until reaching
% either the meridien or the altitude limit (default is 25 degrees). 
%
% It also has a few static functions that convert strings to decimal degrees, 
% using hour2deg for RA-type hour string, and sex2deg for DEC type sexagesimal 
% degree string, converting them to numeric degrees (and the inverse functions). 
%
% 
    
    properties % objects
        
        time@datetime; % time of the observation (could be past or future, usually sets to current time)
        moon; % structure with some details on the moon position etc
        sun; % structure with some details on the sun position etc
        
    end
    
    properties % switches/controls
        
        % observational constraints
        alt_limit = 25; % minimal angle above horizon for being observable (degrees)
        min_duration = 0.5; % minimal time until reaching meridian or horizon (hours)
        
        debug_bit = 0;
        
    end
    
    properties % inputs/outputs
        
        RA_deg; % numeric value in degrees
        Dec_deg; % numeric value in degrees
        
    end
    
    properties(Dependent=true)
        
        RA; % sexagesimal hour string
        Dec; % sexagesimal degree string
        
        STARTTIME; % observation start
        JD; % observation start in julian date
        MJD; % modified julian date

        LST; % local sidereal time (sexagesimal hour string)
        LST_deg; % in numeric degrees
        
        HA; % hour angle (sexagesimal hour string)
        HA_deg; % in numeric degrees
        
        Alt_deg; % altitude above horizon in numeric degrees
        Az_deg; % azimut 
        
        AIRMASS; 
        
        ECL_long; % ecliptic longitude
        ECL_lat; % ecliptic latitude
        
        GAL_long; % galactic longitude
        GAL_lat; % galactic latitude
        
    end
    
    properties(Hidden=true)
        
        RA_deg_now; % in current epoch (rather than in J2000) 
        Dec_deg_now; % in current epoch (rather than in J2000) 
        
        % these are updated when calling updateSecondaryCoords()
        ecliptic_lambda; 
        ecliptic_beta;
        galactic_longitude;
        galactic_latitude;
        
        sidereal_time_type = 'mean'; % can also choose "apparent"
        
        % measured to an accuracy of about 10m with "GPS data" app on my LG G6 phone, with an open dome, at 2020 Feb 16
        elevation = 876; % meters above sea level
        latitude = 30.596806899999997; % degrees North
        longitude = 34.762054900000003; % degrees East
        
%         latitude = 30.59678; % from my phone inside the dome (using whatsmylocation.net)
%         longitude = 34.76202; % from my phone inside the dome (using whatsmylocation.net)
        
%         latitude = 30.5968322; % from my phone (outside the dome)
%         longitude = 34.7619663; % from my phone (outside the dome)
        
%         latitude = 30.59583333333; % WISE observatory (from wikipedia)
%         longitude = 34.76333333333; % WISE observatory (from wikipedia)
        
%         latitude = 30.597296; % updated from google maps
%         longitude = -34.761897; % updated from google maps
        
%         latitude = 31.907867; % of the observatory
%         longitude = 34.811363; % of the observatory
                
        default_fields; % save a few default places to goto when no explicite targets are given
        star_density; % an estimate for the number of stars (right now this is only available for default fields)
        
        version = 1.04;
        
    end
    
    methods % constructor
        
        function obj = Ephemeris(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'head.Ephemeris')
                if obj.debug_bit, fprintf('Ephemeris copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Ephemeris constructor v%4.2f\n', obj.version); end
                
                obj.update;
                
                for ii = 1:length(varargin)
                    
                    if isa(varargin{ii}, 'datetime')
                        obj.time = varargin{ii};
                    end
                    
                end
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj) % remove all coordinate data, update to current time
            
            obj.RA_deg = [];
            obj.Dec_deg = [];
            obj.moon = [];
            obj.sun = [];
            obj.RA_deg_now = [];
            obj.Dec_deg_now = [];

            % these are updated when calling updateSecondaryCoords()
            obj.ecliptic_lambda = []; 
            obj.ecliptic_beta = [];
            obj.galactic_longitude = [];
            obj.galactic_latitude = [];
            
            obj.update;
            
        end
        
    end
    
    methods % getters
        
        function val = get.RA(obj)
            
            val = obj.deg2hour(obj.RA_deg);
            
        end
        
        function val = get.Dec(obj)
            
            val = obj.deg2sex(obj.Dec_deg);
            
        end
        
        function val = get.STARTTIME(obj)
            
            val = util.text.time2str(obj.time);
            
        end
        
        function val = get.JD(obj)

            if isempty(obj.STARTTIME)
                val = [];
            else
                val = juliandate(obj.time); 
%                 val = juliandate(obj.STARTTIME, 'yyyy-mm-ddThh:MM:ss');
            end
             
        end
        
        function val = get.MJD(obj)
            
            if isempty(obj.STARTTIME)
                val = [];
            else
                val = obj.JD - 2400000.5;
            end
        end
        
        function val = get.LST_deg(obj) % stolen the math from Eran's lst function 
            
            import util.text.cs;
            
            % convert JD to integer day + fraction of day
            TJD = floor(obj.JD - 0.5) + 0.5;
            DayFrac = obj.JD - TJD;

            T = (TJD - 2451545.0)./36525.0;

            GMST0UT = 24110.54841 + 8640184.812866.*T + 0.093104.*T.*T - 6.2e-6.*T.*T.*T;

            % convert to fraction of day in range [0 1)
            GMST0UT = GMST0UT./86400.0;

            GMST0UT = GMST0UT - floor(GMST0UT);
            val = GMST0UT + 1.0027379093.*DayFrac + abs(obj.longitude)./360; % note this will fail for targets on the West longitude! 
            val = val - floor(val);
    
            if cs(obj.sidereal_time_type, 'mean')
                
            elseif cs(obj.sidereal_time_type, 'apparent')
                % calculate nutation
                nut_long = obj.nutation(obj.JD);
                Obl    = obj.obliquity(obj.JD);
                EquationOfEquinox = 3600.*nut_long.*cosd(Obl)./15; % convert to time seconds
                val = val + EquationOfEquinox./86400; % convert to fractions of day
            else
                error('Unknown sidereal_time_type. Use "mean" or "apparent"'); 
            end

            val = val.*360; % convert to degrees
            
        end
        
        function val = get.LST(obj) % in sexagesimal format (HH:MM:SS)
            
            val = obj.deg2hour(obj.LST_deg);
            
        end
        
        function val = get.HA_deg(obj)
            
            val = obj.LST_deg - obj.RA_deg;
            
            if val>180
                val = val - 360;
            end

        end
        
        function val = get.HA(obj)
            
            s = sign(obj.HA_deg);
            val = obj.deg2hour(abs(obj.HA_deg));
            if s<0
                val = ['-' val]; 
            end
            
        end
        
        function val = get.Alt_deg(obj)
            
%             val = obj.ha2alt(obj.HA_deg, obj.DEC_deg, obj.latitude);
            
            % reference: http://www.jgiesen.de/elevaz/basics/
            val = asind(sind(obj.Dec_deg).*sind(obj.latitude) + cosd(obj.Dec_deg).*cosd(obj.latitude).*cosd(obj.HA_deg));
            
        end
        
        function val = get.Az_deg(obj)
            
%             val = obj.ha2az(obj.HA_deg, obj.DEC_deg, deg2rad(obj.latitude));
            
            SinAlt = sind(obj.DEC_deg).*sind(obj.latitude) + cosd(obj.DEC_deg).*cosd(obj.HA_deg).*cosd(obj.latitude);
            CosAlt = sqrt(1-SinAlt.*SinAlt);

            SinAz  = (-cosd(obj.DEC_deg).*sind(obj.HA_deg))./CosAlt;
            CosAz  = (sind(obj.DEC_deg).*cosd(obj.latitude) - cosd(obj.DEC_deg).*cosd(obj.HA_deg).*sind(obj.latitude))./CosAlt;

            val = atan2d(SinAz, CosAz);
            
        end
        
        function val = get.AIRMASS(obj)
           
            alt = obj.Alt_deg;
            
            if alt<3% unreliable values near the horizon
                val = []; 
            else
                
                SecZ = 1./cosd(90-alt);

                val = SecZ - 0.0018167.*(SecZ - 1) - 0.002875.*(SecZ - 1).*(SecZ - 1)...
                    - 0.0008083.*(SecZ - 1).*(SecZ - 1).*(SecZ - 1); % Hardie's polynomial formula

            end
            
        end
        
        function val = RA_now(obj)
            
            val = obj.deg2hour(obj.RA_deg_now);
            
        end
        
        function val = Dec_now(obj)
            
            val = obj.deg2sex(obj.Dec_deg_now);
            
        end
        
        function val = get.ECL_long(obj)
            
            val = obj.ecliptic_lambda;
            
        end
        
        function val = get.ECL_lat(obj)
            
            val = obj.ecliptic_beta;
            
        end
        
        function val = get.GAL_long(obj)
            
            val = obj.galactic_longitude;
            
        end
        
        function val = get.GAL_lat(obj)
            
            val = obj.galactic_latitude;
            
        end
        
    end
    
    methods % setters
        
        function set.RA(obj, val) % expect RA to be given in hours, not degrees! 

            if isempty(val)
                obj.RA_deg = [];
            elseif isnumeric(val) && isscalar(val)
                obj.RA_deg = val*15;
            elseif isnumeric(val) && length(val)==3
                val = val*15;
                obj.RA_deg = val(1) + val(2)/60 + val(3)/3600;
            elseif ischar(val)
                obj.RA_deg = obj.hour2deg(val);
            else
                error('Input a string RA in HH:MM:SS format or a scalar hour or a 3-vector [H,M,S]!');
            end
            
            obj.updateSecondaryCoords;
            
        end 
        
        function set.Dec(obj, val)
            
            if isempty(val)
                obj.Dec_deg = [];
            elseif isnumeric(val) && isscalar(val)
                obj.Dec_deg = val;
            elseif isnumeric(val) && length(val)==3
                obj.Dec_deg = sign(val(1))*(abs(val(1)) + val(2)/60 + val(3)/3600);
            elseif ischar(val)
                obj.Dec_deg = obj.sex2deg(val);
            else
                error('Input a string DEC in DD:MM:SS format or a scalar in degrees or a 3-vector [d,m,s]!');
            end
            
            obj.updateSecondaryCoords;
            
        end
        
    end
    
    methods % update moon/sun and coordinates
        
        function update(obj) % set to current time and update everything else
            
            obj.time = datetime('now', 'timezone', 'UTC');
            
            obj.updateSecondaryCoords;
            
        end
        
        function updateSecondaryCoords(obj) % update the moon/sun and secondary coordinates 
            
            obj.updateMoon;
            obj.updateSun;
            obj.updateEquatorialNow;
            obj.updateEcliptic;
            obj.updateGalactic; 
            
        end
        
        function updateSun(obj) % make a struct with sun properties
            
            RAD = pi./180;
            
            if ~isempty(which('celestial.SolarSys.get_sun'))
                obj.sun = celestial.SolarSys.get_sun(obj.JD, [obj.longitude, obj.latitude].*pi./180);
                obj.sun.RA = obj.sun.RA./RAD; % convert to degrees
                obj.sun.Dec = obj.sun.Dec./RAD; % convert to degrees
                obj.sun.Az = obj.sun.Az./RAD; % convert to degrees
                obj.sun.Alt = obj.sun.Alt./RAD; % convert to degrees
                obj.sun.dAzdt = obj.sun.dAzdt./RAD; % convert to degrees
                obj.sun.dAltdt = obj.sun.dAltdt./RAD; % convert to degrees
                
            else
                obj.sun = [];
            end
            
        end
        
        function updateMoon(obj) % make a struct with moon properties
            
            RAD = pi./180;
            
            if ~isempty(which('celestial.SolarSys.get_moon'))
                obj.moon = celestial.SolarSys.get_moon(obj.JD, [obj.longitude, obj.latitude].*RAD);
                obj.moon.RA = obj.moon.RA./RAD; % convert to degrees
                obj.moon.Dec = obj.moon.Dec./RAD; % convert to degrees
                obj.moon.Az = obj.moon.Az./RAD; % convert to degrees
                obj.moon.Alt = obj.moon.Alt./RAD; % convert to degrees
                obj.moon.Phase = obj.moon.Phase./pi; % convert to fraction
                obj.moon.Dist = obj.getMoonDistance; 
            else
                obj.moon = [];
            end
            
        end
        
        function updateEquatorialNow(obj) % convert J2000 to current equinox
           
            if ~isempty(which('celestial.coo.coco')) && ~isempty(obj.RA_deg) && ~isempty(obj.Dec_deg)
                
                JY = convert.time(obj.JD,'JD','J');
                CurrentEquinox = sprintf('J%8.3f',JY);
                out_coord = celestial.coo.coco([obj.RA_deg, obj.Dec_deg], 'J2000', CurrentEquinox, 'd', 'd');

                obj.RA_deg_now = out_coord(1);
                obj.Dec_deg_now = out_coord(2);
                
            else
                obj.RA_deg_now = [];
                obj.Dec_deg_now = [];
            end
            
        end
        
        function updateEcliptic(obj) % calculate ecliptic coordinates
            
            if ~isempty(which('celestial.coo.coco')) && ~isempty(obj.RA_deg) && ~isempty(obj.Dec_deg)
                out_coord = celestial.coo.coco([obj.RA_deg, obj.Dec_deg], 'J2000', 'e', 'd', 'd');
                obj.ecliptic_lambda = out_coord(1);
                obj.ecliptic_beta = out_coord(2);
            else
                obj.ecliptic_lambda = [];
                obj.ecliptic_beta = [];
            end
            
        end
        
        function updateGalactic(obj) % calculate galactic coordinates
            
            if ~isempty(which('celestial.coo.coco')) && ~isempty(obj.RA_deg) && ~isempty(obj.Dec_deg)
                out_coord = celestial.coo.coco([obj.RA_deg, obj.Dec_deg], 'J2000', 'g', 'd', 'd');
                obj.galactic_longitude = out_coord(1);
                obj.galactic_latitude = out_coord(2);
            else
                obj.ecliptic_lambda = [];
                obj.ecliptic_beta = [];
            end
            
        end
        
    end
    
    methods % calculations
        
        function input(obj, RA, DEC, time) % give RA/DEC or RA=<star name>, DEC=[] (optional 3rd argument is time, can use "now")
        % Usage: input(obj, RA, DEC, <time>) 
        % Give the coordinates in one of two ways:
        % (1) give the RA and Dec directly. This can be in numeric or string
        %     format. In both formats RA is given in HOURS!
        % (2) give the name of the object and leave DEC empty. 
        %     Will use Eran's name resolver via SIMBAD. 
        %
        % An optional third parameter is the time of observation. This can 
        % be one of several options: 
        % (1) string "update" or "now" simply updates to current time. 
        % (2) datetime object. Must use timezone UTC! 
        % (3) string in YYYY-MM-DDThh:mm:ss.sss format. 
        % (4) numeric value, interpreted as julian date. 
        % 
        % In either case, the secondary coordinates are updated to the new
        % RA/Dec and time. 
            
            import util.text.cs;
        
            if nargin==1, help('head.Ephemeris.input'); return; end
            
            if nargin<3
                DEC = [];
            end
            
            if ischar(RA) && isempty(DEC)
                name = RA; 
                if cs(name, 'ecliptic', 'kbos')
                    obj.gotoDefaultField('ecliptic'); 
                    RA = obj.RA;
                    DEC = obj.Dec;
                elseif cs(name, 'galactic')
                    obj.gotoDefaultField('galactic'); 
                    RA = obj.RA;
                    DEC = obj.Dec;
                elseif cs(name, 'moon')
                    RA = 0;
                    DEC = 0;
                elseif ~isempty(which('celestial.coo.coo_resolver', 'function'))
                    [RA, DEC] = celestial.coo.coo_resolver(name, 'OutUnits', 'deg', 'NameServer', @VO.name.server_simbad);
                    RA = RA/15;
                else
                    error('no name resolver has been found... try adding MAAT to the path.');
                end 
                
            end
            
            obj.RA = RA;
            obj.Dec = DEC;
            
            if nargin>3 && ~isempty(time)
                
                if isa(time, 'datetime')
                    obj.time = time;
                    obj.updateSecondaryCoords;
                elseif ischar(time) && util.text.cs(time, 'now', 'update')
                    obj.update;
                elseif ischar(time)
                    obj.time = util.text.str2time(time);
                    obj.updateSecondaryCoords;
                elseif isnumeric(time) % juldate??
                    obj.time = datetime(time, 'ConvertFrom', 'juliandate', 'TimeZone', 'UTC');
                    obj.updateSecondaryCoords;
                end
                
            end
            
            if cs(name, 'moon')
                obj.RA = obj.moon.RA/15;
                obj.Dec = obj.moon.Dec;
            end
            
            
        end
        
        function timeTravelHours(obj, H) % move a few hours to the future (positive) or past (negative)
        % Usage: timeTravelHours(obj, H)
        % Push the time of this object by H hours. Positive H is to the future, 
        % negative H is to the past. RA and Dec remain the same. 
        % All other coordinates are adjusted accordingly. 
        % 
        % Use this to check what observing conditions would be for a future 
        % observation, or check what they were for a past one. 
        
            if nargin==1, help('head.Ephemeris.timeTravelHours'); return; end
        
            obj.time = obj.time + hours(H);
            
            obj.updateSecondaryCoords;
            
        end
        
        function val = calcObsTimeMinutes(obj, alt_limit, resolution_minutes) % how much time is left until object reaches meridien or alt-limit (minutes)
        % Usage: val = calcObsTimeMinutes(obj, alt_limit=25, resolution_minutes=5) 
        % 
        % Push the same coordinates forward in time until reaching the 
        % meridien or until going under the "alt_limit" (default is 25 degrees). 
        % Will return the time rounded to "resolution_minutes" (default 5). 
        % The output is in minutes! 
        % For objects below the alt_limit at the stored time, returns NaN. 
        
            if nargin<2 || isempty(alt_limit)
                alt_limit = obj.alt_limit;
            end
            
            if nargin<3 || isempty(resolution_minutes)
                resolution_minutes = 5;
            end
            
            if obj.ALT_deg<alt_limit
                val = NaN; % we cannot observe this object at the given time, so it returns NaN
                return; 
            end
            
            if obj.HA_deg<0
                val = round(abs(obj.HA_deg/15*60)/resolution_minutes)*resolution_minutes; 
            else
                
                e = head.Ephemeris(obj); % copy the object! 
                
                val = 0;
                
                for ii = 1:1e3
                    
                    e.time = e.time + minutes(resolution_minutes); 
                    
                    if e.ALT_deg>alt_limit
                        val = val + resolution_minutes;
                    else
                        break;
                    end
                    
                end
                
            end
            
        end
        
        function val = getMoonDistance(obj) % get distance between moon and target RA/Dec
        % Usage: val = getMoonDistance(obj) 
        % Calculate the distance of the moon (based on the latest call to 
        % updateMoon()) and return that value, in degrees. 
        % 
        % reference: https://en.wikipedia.org/wiki/Haversine_formula
            
            if isempty(obj.moon) || isempty(obj.RA_deg) || isempty(obj.DEC_deg)
                val = [];
            else
                
                havTheta = sind((obj.DEC_deg-obj.moon.Dec)/2).^2 + cosd(obj.DEC_deg).*cosd(obj.moon.Dec).*sind((obj.RA_deg-obj.moon.RA)/2).^2;
                
                havTheta(havTheta>1) = 1;
                havTheta(havTheta<-1) = -1;
                
                val = 2.*asind(sqrt(havTheta)); 
                
            end
            
        end
        
        function random(obj, varargin) % find a random point in the observable sky
        % Usage: random(obj, varargin)
        % Choose new values for RA/Dec somewhere in the sky that is observable. 
        % Will choose random values until one of them is above the alt limit, 
        % and above the south limit in declination (default is -20). 
        % 
        % OPTIONAL ARGUMENTS:
        %   -alt_limit: choose the minimal altitude that is considered 
        %               observable at the given time. Default 25 degrees. 
        %   -south_limit: choose how far south in declination you wish to 
        %                 use as valid target. Default -20 degrees. 
        %   -side: choose "east" or "west" or "both". 
        %   -meridian_distance: number of degrees from meridian (on either 
        %                       side) that we want to avoid. Default 5 degrees. 
        
            input = util.text.InputVars;
            input.input_var('alt_limit', obj.alt_limit); 
            input.input_var('south', -20, 'south_limit');
            input.input_var('side', 'both', 'hemisphere');
            input.input_var('meridian', 5, 'meridian_distance', 'meridien_distance'); 
            input.scan_vars(varargin{:}); 
            
            if util.text.cs(input.side, 'both')
                min_RA = 6; 
                range_RA = 12;
            elseif util.text.cs(input.side, 'east')
                min_RA = 6; 
                range_RA = 6;
            elseif util.text.cs(input.side, 'west')
                min_RA = 0;
                range_RA = 6;
            else
                error('Unknown side/hemisphere option "%s". Choose "east", "west" or "both"', input.side); 
            end
            
            for ii = 1:100
                
                obj.RA = obj.LST_deg/15 + rand.*range_RA-min_RA; % random between min_RA and min_RA+range_RA (hours)
                
                obj.Dec = rand.*(90-input.south)+input.south; % between south_limit and 90
                
                if obj.Alt_deg > input.alt_limit && abs(obj.HA_deg)>input.meridian
                    obj.updateSecondaryCoords;
                    return;
                end
                
            end
            
            error('Could not find a random pointing inside observational bounds!'); 
            
        end
        
        function val = better_than(obj, other, varargin) % find if this object has better observable conditions than that one
        % Usage: val = better_than(obj, other, varargin)
        % Check which target is better to observe right now. 
        % Returns true when this object is better. 
        % Must take two different Ephemeris objects. Assumes objects are 
        % updated to the same time (this is up to the user). 
        % 
        % In the simplest case, the target with lower airmass will be better. 
        %
        % OPTIONAL ARGUMENTS:
        %   -alt_limit: minimal altitude (angle) above horizon for a target
        %               to be viable/observable. Default is 25 degrees.
        %   -duration: minimal duration (in hours) that a target must be
        %              observable, until meridian or horizon, to be valid. 
        %              Default is 0.5 hours. 
        %   -side: can only choose viable targets on the given side, that is
        %          "East", "West", or "both" (default, same as empty). 
        %          Note that this overwrides the "wind_speed" input. 
        %   -wind_speed: if wind is above some limit (20km/h) then this 
        %                function will prefer the Eastern target regardless
        %                of airmass. If both are on the same side, this is 
        %                ignored. If not given (or empty) will skip this test. 
        % 
        % NOTE: if the other object is not observable, then this object is 
        %       accepted as default (val==1), even if it isn't observable. 
        %       If other is observable and this object isn't, only then is 
        %       the other object preferred (val==0). 
            
            import util.text.cs;
        
            if nargin==1, help('head.Ephemeris.better_than'); return; end
        
            input = util.text.InputVars;
            input.input_var('alt_limit', obj.alt_limit); 
            input.input_var('duration', obj.min_duration); 
            input.input_var('side', []); 
            input.input_var('wind_speed', []); 
            input.scan_vars(varargin{:}); 
        
            if isempty(other) || ~isa(other, 'head.Ephemeris') || obj==other
                error('Must supply another Ephemeris object to this function!'); 
            end
            
            if isempty(other.AIRMASS) || other.ALT_deg<input.alt_limit || other.calcObsTimeMinutes(input.alt_limit)<input.duration*60
                val = 1;
                return; 
            end
            
            if isempty(obj.AIRMASS) || obj.ALT_deg<input.alt_limit || obj.calcObsTimeMinutes(input.alt_limit)<input.duration*60
                val = 0;
                return;
            end
            
            val = obj.AIRMASS<=other.AIRMASS; % if this object has lower airmass, it is better
            
            %%%%%%%%% additional checks besides airmass / alt_limit %%%%%%%%%%
            
            % check if the wind is too high, prefer Eastern targets
            if ~isempty(input.wind_speed)
                if input.wind_speed>20 % only prefer Eastern targets when the wind is strong
                    if obj.HA_deg<0 && other.HA_deg>=0 % this is Eastern while other is Western
                        val = 1;
                    elseif obj.HA_deg>=0 && other.HA_deg<0 % this is Western while other is Eastern
                        val = 0;
                    else
                        % pass: if both targets are on the same side, do not use wind as a deciding factor
                    end
                end
            end
            
            % if given a side to prefer (to prevent meridian flips)
            if ~isempty(input.side) && ~cs(input.side, 'both')
                
                if cs(input.side, 'east')
                    
                    if obj.HA_deg<0 && other.HA_deg>=0 % this is Eastern while other is Western
                        val = 1;
                    elseif obj.HA_deg>=0 && other.HA_deg<0 % this is Western while other is Eastern
                        val = 0;
                    else
                        % pass: if both targets are on the same side, do not use wind as a deciding factor
                    end
                    
                elseif cs(input.side, 'west')
                    
                    if obj.HA_deg>=0 && other.HA_deg<0 % this is Western while other is Eastern
                        val = 1;
                    elseif obj.HA_deg<0 && other.HA_deg>=0 % this is Eastern while other is Western
                        val = 0;
                    else
                        % pass: if both targets are on the same side, do not use wind as a deciding factor
                    end
                    
                else
                    error('Unknown side "%s". Choose "East" or "West" or "both". ', input.side); 
                end
                
            end
            
        end
        
    end
    
    methods % utilities
        
        function struct2obj(obj, s) % must first construct an object and then use the struct to update its properties
            
            list = properties(obj);
            
            for ii = 1:length(list)
               
                if isfield(s, list{ii})
                    
                    if isobject(obj.(list{ii}))
                        
                        if isprop(obj.(list{ii}), 'struct2obj')
                            struct2obj(obj.(list{ii}), s.list{ii}); % call the sub-object method if it exists
                        else
                            
                            list2 = properties(obj.(list{ii}));
                            
                            for jj = 1:length(list2)
                            
                                if isfield(s.(list{ii}), list2{jj})
                                    
                                    try 
                                        obj.(list{ii}).(list2{jj}) = s.(list{ii}).(list2{jj}); 
                                    end
                                    
                                end
                                
                            end
                            
                        end
                        
                    else % regular data types
                        try % we can run into many problems trying to set unsetable properties
                            obj.(list{ii}) = s.(list{ii}); 
                        end
                    end
                    
                end
                
            end
            
        end
        
    end
    
    methods % produce and load default fields
        
        function gotoDefaultField(obj, type, varargin) % jump to one of the default fields by type and number
        % Usage: gotoDefaultField(obj, type, number=[])
        % Go to one of the default fields, from the given type. 
        % Choose one of the types: KBOs, ...
        % If the number is given, will go to that specific field. 
        % (use numberDefaultFields(type) to know how many there are). 
        % If no number is given, simply choose the best field. 
        % 
        % OPTIONAL ARGUMENTS are passed along to the better_than() function. 
        
            import util.text.cs;
            
            if nargin==1, help('head.Ephemeris.gotoDefaultField'); return; end

            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('number', []);
            input.scan_vars(varargin{:}); 
            
            if ischar(input.number)
                if cs(input.number, 'best')
                    input.number=[];
                else
                    error('Wrong input for number: "%s". Use "best" or a numeric value.', input.number);
                end
            end
            
            if isempty(obj.default_fields)
                obj.makeDefaultFields;
            end
            
            if ~isempty(input.number) && (input.number<1 || input.number>obj.numberDefaultFields(type))
                error('Requested field number %d is out of bounds for type %s (limited to %d fields).', input.number, type, obj.numberDefaultFields(type)); 
            end
            
            if cs(type, 'kbos', 'kuiper belt objects', 'ecliptic')
                s = obj.default_fields(strcmp({obj.default_fields.type}, 'ecliptic')); % get all the structs that fit this type of field
            elseif cs(type, 'galactic')
                s = obj.default_fields(strcmp({obj.default_fields.type}, 'galactic')); % get all the structs that fit this type of field
            % add other default field types...
            else
                error('Unknown default field type "%s". Try "ecliptic"... ', type); 
            end
            
            
            if isempty(input.number) % automatically choose the best field
                
                input.number = 1; % assume it is 1 and adjust it if a better field is found
                
                obj.RA_deg = s(1).RA;
                obj.Dec_deg = s(1).Dec;
                obj.star_density = s(1).stars;
                obj.updateSecondaryCoords;
                
                e = head.Ephemeris; % make a test object 
                e.time = obj.time;
                
                for ii = 2:length(s)
                    
                    e.RA_deg = s(ii).RA;
                    e.Dec_deg = s(ii).Dec;
                    e.star_density = s(ii).stars;
                    e.updateSecondaryCoords;
                    
                    if better_than(e, obj, varargin{:}) % found a better field, update to that one
                        input.number = ii;
                        obj.RA_deg = s(ii).RA;
                        obj.Dec_deg = s(ii).Dec;
                        obj.star_density = s(ii).stars;
                        obj.updateSecondaryCoords; 
                    end
                    
                end
                
            else
                s = s(input.number);
                obj.RA_deg = s.RA;
                obj.Dec_deg = s.Dec;
                obj.star_density = s.stars;
                obj.updateSecondaryCoords;
            end
            
            
            
        end
        
        function val = numberDefaultFields(obj, type) % check how many default fields are available for this type
        % Usage: val = numberDefaultFields(obj, type)
        % ...
        
            import util.text.cs;
            
            if nargin==1, help('head.Ephemeris.numberDefaultFields'); return; end
            
            if isempty(obj.default_fields)
                obj.makeDefaultFields;
            end
            
            if cs(type, 'kbos', 'kuiper belt objects', 'ecliptic')
                val = nnz(strcmp({obj.default_fields.type}, 'ecliptic'));
            elseif cs(type, 'galactic')
                val = nnz(strcmp({obj.default_fields.type}, 'galactic'));
            % add other default field types...
            else
                error('Unknown default field type "%s". Try "ecliptic"... ', type); 
            end
            
        end
        
        function makeDefaultFields(obj) % make a few fields that we know are useful for certain types of survey
  
            obj.default_fields = [];
            
            s = struct;
            
            % I chose these using the scripts/scan_ecliptic.m script and the sky map (down to 14th mag)
            
            %%%%%%%%%%%%%%%%%%%%%% Ecliptic fields %%%%%%%%%%%%%%%%%%%%%%%%
            
            s.type = 'ecliptic'; s.RA = 10.5; s.Dec = 6.0; s.ECL_LAT = 1.364; s.GAL_LAT = -56.793; s.stars = 1000;
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'ecliptic'; s.RA = 29; s.Dec = 17.0; s.ECL_LAT = 4.642; s.GAL_LAT = -43.01; s.stars = 1000; 
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'ecliptic'; s.RA = 58.0; s.Dec = 22.0; s.ECL_LAT = 1.772; s.GAL_LAT = -24.28; s.stars = 2000;
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'ecliptic'; s.RA = 98.5; s.Dec = 21.5; s.ECL_LAT = -1.706; s.GAL_LAT = 5.98; s.stars = 15000;
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'ecliptic'; s.RA = 120.0; s.Dec = 19.0; s.ECL_LAT = -1.548; s.GAL_LAT = 23.357; s.stars = 4000;
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'ecliptic'; s.RA = 156.5; s.Dec = 8.0; s.ECL_LAT = -1.684; s.GAL_LAT = 50.643;s.stars = 1000;
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'ecliptic'; s.RA = 209.5; s.Dec = -10.5; s.ECL_LAT = 1.455; s.GAL_LAT = 49.055; s.stars = 1700;
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'ecliptic'; s.RA = 235.5; s.Dec = -19.5; s.ECL_LAT = 0.158; s.GAL_LAT = 27.706; s.stars = 4000;
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'ecliptic'; s.RA = 255.5; s.Dec = -24.5; s.ECL_LAT = -1.72; s.GAL_LAT = 10.55; s.stars = 9400;
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'ecliptic'; s.RA = 264.0; s.Dec = -18.0; s.ECL_LAT = 5.320; s.GAL_LAT = 7.625; s.stars = 7000;
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'ecliptic'; s.RA = 275; s.Dec = -25.0; s.ECL_LAT = -1.639; s.GAL_LAT = -4.727; s.stars = 20000;
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'ecliptic'; s.RA = 305.5; s.Dec = -17.5; s.ECL_LAT = 1.889; s.GAL_LAT = -27.636; s.stars = 5000;
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'ecliptic'; s.RA = 331.0; s.Dec = -10.5; s.ECL_LAT = 1.285; s.GAL_LAT = -47.294; s.stars = 2000;
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'ecliptic'; s.RA = 331.0; s.Dec = -10.5; s.ECL_LAT = 1.285; s.GAL_LAT = -47.294; s.stars = 2000;
            obj.default_fields = [obj.default_fields; s];
            
            %%%%%%%%%%%%%%%%%%%%%% galactic fields %%%%%%%%%%%%%%%%%%%%%%%%
            
            s.type = 'galactic'; s.RA = 34.5; s.Dec = 53; s.ECL_LAT = 36.666; s.GAL_LAT = -7.683; s.stars = 9000; 
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'galactic'; s.RA = 89.5; s.Dec = 28.5; s.ECL_LAT = 5.061; s.GAL_LAT = 2.112; s.stars = 14000;
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'galactic'; s.RA = 119.0; s.Dec = -25.0; s.ECL_LAT = -44.672; s.GAL_LAT = 1.806; s.stars = 31000;
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'galactic'; s.RA = 150.0; s.Dec = -25.0; s.ECL_LAT = -34.611; s.GAL_LAT = 23.488; s.stars = 4000;
            obj.default_fields = [obj.default_fields; s];

            s.type = 'galactic'; s.RA= 209.5; s.Dec = -25.0; s.ECL_LAT = -12.135; s.GAL_LAT = 35.459; s.stars = 2600;
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'galactic'; s.RA = 281.5; s.Dec = -7.5; s.ECL_LAT = 15.468; s.GAL_LAT = -2.211; s.stars = 37500;
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'galactic'; s.RA = 300.0; s.Dec = 38.0; s.ECL_LAT= 56.753; s.GAL_LAT = 4.246; s.stars = 20500;
            obj.default_fields = [obj.default_fields; s];
            
            s.type = 'galactic'; s.RA = 333.5; s.Dec = 51.5; s.ECL_LAT = 55.947; s.GAL_LAT = -4.101; s.stars = 20000;
            obj.default_fields = [obj.default_fields; s];
            
            % add more fields! 
            
            
        end
        
        function printDefaultFields(obj, type) % print all fields of a specific type
            
            import util.text.cs;
            
            if cs(type, 'ecliptic', 'kbos', 'kuiper belt objects')
                type = 'ecliptic';
            elseif cs(type, 'galactic')
                type = 'galactic';
            else
                error('Unknown default field type "%s". Try "ecliptic" or "galactic".', type);
            end
            
            s = obj.default_fields;
            
            c = 1;
            
            for ii = 1:length(s)
                
                if cs(s(ii).type, type)
                    
                    fprintf('id= %02d | RA= %5.1f | Dec= %5.1f | ECL_LAT= %7.3f | GAL_LAT= %7.3f | stars= %d \n', ...
                        c, s(ii).RA, s(ii).Dec, s(ii).ECL_LAT, s(ii).GAL_LAT, s(ii).stars); 
                    
                    c = c + 1;
                    
                end
                    
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end
    
    methods (Static=true) % conversion tools and some additional calculations I stole from Eran's code
        
        function str = deg2hour(number) % convert the degrees into hours then into HH:MM:SS string
        % Usage: str = deg2hour(number)
        % Convert numeric degrees to sexagesimal hour string
        % Example: 180 --> '12:00:00.0'
        
            if nargin==0, help('head.Ephemeris.deg2hour'); return; end
            
            str = head.Ephemeris.numbers2hour(number/15); 
            
        end
        
        function str = deg2sex(number) % convert degrees into DD:MM:SS string
        % Usage: str = deg2sex(number)
        % Convert numeric degrees to sexagesimal degree string
        % Example: -45.5 --> '-45:30:00.0'
        
            if nargin==0, help('head.Ephemeris.deg2sex'); return; end
            
            str = head.Ephemeris.numbers2sex(number);
            
        end
        
        function val = hour2deg(str) % convert an HH:MM:DD hour string to degrees
        % Usage: val = hour2deg(str)
        % Convert sexagesimal hour string numeric degrees
        % Example: '18:26:30.0' --> 276.375 deg
        
            if nargin==0, help('head.Ephemeris.hour2deg'); return; end
            
            num = util.text.extract_numbers(str);
            
            if isempty(num)
                val = [];
            else
            
                hours = num{1}(1);
                degrees = hours*15;
                
                s = 1;
                if degrees<0, s=-1; end
                
                degrees = abs(degrees);
                
                if length(num{1})>1
                    minutes = num{1}(2)*15;
                else
                    minutes = 0;
                end
                
                if length(num{1})>2
                    seconds = num{1}(3)*15;
                else
                    seconds = 0;
                end
                
                val = degrees+minutes/60+seconds/3600;
                val = val.*s;
                
            end
            
        end
        
        function val = sex2deg(str) % convert a DD:MM:SS string (sexagesimal) to degrees
        % Usage: val = sex2deg(str)
        % Convert sexagesimal degree string to numeric degrees
        % Example: '22:33:44.5' --> 22.562375 deg
        
            if nargin==0, help('head.Ephemeris.sex2deg'); return; end 
            
            str = strip(str);
            
            num = util.text.extract_numbers(str);
            
            if isempty(num)
                val = [];
            else
            
                degrees = num{1}(1);
                s = 1;
                if degrees<0, s=-1; end
                if degrees==0 && str(1)=='-', s=-1; end
                
                degrees = abs(degrees);
                
                if length(num{1})>1
                    minutes = num{1}(2);
                else
                    minutes = 0;
                end
                
                if length(num{1})>2
                    seconds = num{1}(3);
                else
                    seconds = 0;
                end
                
                val = degrees+minutes/60+seconds/3600;
                val = val.*s;
                
            end
            
        end
        
        function str = numbers2hour(hours, minutes, seconds) % take hours or hours/minutes/seconds and produce sexagesimal hour string
        % Usage: str = numbers2hour(hours, minutes, seconds) 
        % Convert a number of hours (or hours, minutes, seconds) into a 
        % string of sexagesimal hours (e.g., for RA). 
        % Can input one number for the hours, or a vector with three components
        % for [hour minute second] or just input them in three inputs. 
        %
        % The output is a string with HH:MM:SS.S format. 
        
            if nargin==0, help('head.Ephemeris.numbers2hour'); return; end
        
            if isempty(hours)
                str = '';
                return;
            end
            
            if iscell(hours)
                hours = cell2mat(hours);
            end
            
            if nargin<2 || isempty(minutes)
                minutes = 0;
            end
            
            if nargin<3 || isempty(seconds)
                seconds = 0;
            end
            
            if size(hours,2)==3
                minutes = hours(:,2);
                seconds = hours(:,3);
                hours = hours(:,1);
            end
            
            s = sign(hours);
            hours = abs(hours); 
            
            total_secs = hours*3600 + minutes*60 + seconds;
            total_secs = mod(total_secs, 24*3600); 
            
            hours   = floor(total_secs/3600);
            minutes = floor(mod(total_secs,3600)/60);
            seconds = mod(total_secs,60);
            
            str = sprintf('%02d:%02d:%04.1f', hours, minutes, seconds);
            
            if s<0
                str = ['-' str];
            end
            
        end
        
        function str = numbers2sex(degrees, minutes, seconds) 
        % Usage: str = numbers2sex(degrees, minutes, seconds) 
        % Convert a number of degrees (or degrees, arcminutes, arcseconds) 
        % into a string of sexidecimal degrees. 
        % Can input one parameter with the degrees, or a vector with three
        % elements for [deg, arcmin, arcsec], or given these three inputs
        % as individual input values. 
        %
        % Output is degrees in a sexagesimal string: +DD:MM:SS.S with + or -
        
            if nargin==0, help('head.Ephemeris.numbers2sex'); return; end
            
            if isempty(degrees)
                str = '';
                return;
            end
            
            if iscell(degrees)
                degrees = cell2mat(degrees);
            end
            
            if nargin<2 || isempty(minutes)
                minutes = 0;
            end
            
            if nargin<3 || isempty(seconds)
                seconds = 0;
            end
            
            if size(degrees,2)==3
                minutes = degrees(:,2);
                seconds = degrees(:,3);
                degrees = degrees(:,1);
            end
            
            if minutes<0 || seconds<0
                error('Why did we get negative minutes/seconds??');
            end
            
            s = sign(degrees);
            degrees = abs(degrees);
            total_secs = degrees*3600 + minutes*60 + seconds;
            
            degrees = fix(total_secs/3600);
            minutes = fix(mod(total_secs,3600)/60);
            seconds = mod(total_secs,60);
            
            if s>=0
                sign_str = '+';
            else
                sign_str = '-';
            end
                
            str = sprintf('%s%02d:%02d:%04.1f', sign_str, degrees, minutes, seconds);
            
        end
        
        function Obl=obliquity(JulianDay,Type)
        % stolen from Eran's obliquity function 
            
            if nargin<2 || isempty(Type)
                Type = 'L';
            end
        
            switch Type
                case 'L'
                    T   = (JulianDay - 2451545.0)./36525.0;
                    Obl = 23.439291 - 0.0130042.*T - 0.00000016.*T.*T + 0.000000504.*T.*T.*T;

                case 'H'
                    T   = (JulianDay - 2451545.0)./36525.0;
                    U   = T./100;
                    Obl = 23.44484666666667 ...
                        +   (-4680.93.*U ...
                        - 1.55.*U.^2 ...
                        + 1999.25.*U.^3 ...
                        - 51.38.*U.^4 ...
                        - 249.67.*U.^5 ...
                        - 39.05.*U.^6 ...
                        + 7.12.*U.^7 ...
                        + 27.87.*U.^8 ...
                        + 5.79.*U.^9 ...
                        + 2.45.*U.^10)./3600;

                otherwise
                    error('Unknown calculation type in obliquity.m');
            end
            
        end
        
        function [nut_long, nut_obli] = nutation(JD)
        % stolen from Eran's nutation function 
        
            T   = (JD - 2451545.0)./36525.0;
            
            % Mean elongation of the Moon from the Sun:
            D = (297.85036 + 445267.111480.*T - 0.0019142.*T.*T + T.*T.*T./189474);
            
            % Mean anomaly of the Sun (Earth):
            M = (357.52772 + 35999.050340.*T - 0.0001603.*T.*T - T.*T.*T./300000);
            
            % Mean anomaly of the Moon:
            Mt = (134.96298 + 477198.867398.*T + 0.0086972.*T.*T + T.*T.*T./56250);
            
            % Moon's argument of latitude
            F = (93.27191 + 483202.017538.*T - 0.0036825.*T.*T + T.*T.*T./327270);
            
            % Longitude of ascending node of the Moon's mean orbit
            Om = (125.04452 - 1934.136261.*T + 0.0020708.*T.*T + T.*T.*T./450000);
            
            % nutation in longitude [ 0."0001]
            DLon = (-171996 - 174.2.*T).*sind(Om) + ...,
                ( -13187 +   1.6.*T).*sind(2.*(-D + F + Om)) + ...,
                (  -2274 -   0.2.*T).*sind(2.*(F + Om)) + ...,
                (   2062 +   0.2.*T).*sind(2.*Om) + ...,
                (   1426 -   3.4.*T).*sind(M) + ...,
                (    712 +   0.1.*T).*sind(Mt) + ...,
                (   -517 +   1.2.*T).*sind(2.*(-D + F + Om) + M) + ...,
                (   -386 -   0.4.*T).*sind(2.*F + Om) + ...,
                (   -301           ).*sind(Mt + 2.*(F + Om)) + ...,
                (    217 -   0.5.*T).*sind(-M + 2.*(-D + F + Om)) + ...,
                (   -158           ).*sind(-2.*D + Mt) + ...,
                (    129 +   0.1.*T).*sind(2.*(-D + F) + Om) + ...,
                (    123           ).*sind(-Mt + 2.*(F + Om)) + ...,
                (     63           ).*sind(2.*D) + ...,
                (     63 +   0.1.*T).*sind(Mt + Om) + ...,
                (    -59           ).*sind(-Mt + 2.*(D + F + Om)) + ...,
                (    -58 -   0.1.*T).*sind(-Mt + Om) + ...,
                (    -51           ).*sind(Mt + 2.*F + Om) + ...,
                (     48           ).*sind(2.*(-D + Mt)) + ...,
                (     46           ).*sind(Om + 2.*(-Mt + F)) + ...,
                (    -38           ).*sind(2.*(D + F + Om)) + ...,
                (    -31           ).*sind(2.*(Mt + F + Om)) + ...,
                (     29           ).*sind(2.*Mt) + ...,
                (     29           ).*sind(Mt + 2.*(-D + F + Om)) + ...,
                (     26           ).*sind(2.*F) + ...,
                (    -22           ).*sind(2.*(-D + F)) + ...,
                (     21           ).*sind(-Mt + Om + 2.*F) + ...,
                (     17 -   0.1.*T).*sind(2.*M) + ...,
                (     16           ).*sind(2.*D - Mt + Om) + ...,
                (    -16 +   0.1.*T).*sind(2.*(-D + M + F + Om)) + ...,
                (    -15           ).*sind(M + Om) + ...,
                (    -13           ).*sind(-2.*D + Mt + Om) + ...,
                (    -12           ).*sind(-M + Om) + ...,
                (     11           ).*sind(2.*(Mt - F)) + ...,
                (    -10           ).*sind(2.*(D + F) - Mt + Om) + ...,
                (     -8           ).*sind(2.*(D + F + Om) + Mt) + ...,
                (      7           ).*sind(2.*(F + Om) + M) + ...,
                (     -7           ).*sind(-2.*D + M + Mt) + ...,
                (     -7           ).*sind(-M + 2.*(F + Om)) + ...,
                (     -7           ).*sind(2.*(D + F) + Om) + ...,
                (      6           ).*sind(2.*D + Mt) + ...,
                (      6           ).*sind(2.*(-D + Mt + F + Om)) + ...,
                (      6           ).*sind(2.*(-D + F) + Mt + Om) + ...,
                (     -6           ).*sind(2.*(D - Mt) + Om) + ...,
                (     -6           ).*sind(2.*D + Om) + ...,
                (      5           ).*sind(-M + Mt) + ...,
                (     -5           ).*sind(2.*(F - D) + Om - M) + ...,
                (     -5           ).*sind(Om - 2.*D) + ...,
                (     -5           ).*sind(2.*(Mt + F) + Om) + ...,
                (      4           ).*sind(2.*(Mt - D) + Om) + ...,
                (      4           ).*sind(2.*(F - D) + M + Om) + ...,
                (      4           ).*sind(Mt - 2.*F) + ...,
                (     -4           ).*sind(Mt - D) + ...,
                (     -4           ).*sind(M -2.*D) + ...,
                (     -4           ).*sind(D) + ...,
                (      3           ).*sind(Mt + 2.*F) + ...,
                (     -3           ).*sind(2.*(F + Om - Mt)) + ...,
                (     -3           ).*sind(Mt - D - M) + ...,
                (     -3           ).*sind(M + Mt) + ...,
                (     -3           ).*sind(Mt - M + 2.*(F - Om)) + ...,
                (     -3           ).*sind(2.*(D + F + Om) - M - Mt) + ...,
                (     -3           ).*sind(3.*Mt + 2.*(F + Om)) + ...,
                (     -3           ).*sind(2.*(D + F + Om) - M);
            
            nut_long = DLon./10000/3600; % convert to degrees...
            
            if nargout>1
                % nutation in obliquity [ 0."0001]
                DObl = (  92025 +   8.9.*T).*cosd(Om) + ...,
                    (   5736 -   3.1.*T).*cosd(2.*(-D + F + Om)) + ...,
                    (    977 -   0.5.*T).*cosd(2.*(F + Om)) + ...,
                    (   -895 +   0.5.*T).*cosd(2.*Om) + ...,
                    (     54 -   0.1.*T).*cosd(M) + ...,
                    (     -7           ).*cosd(Mt) + ...,
                    (    224 -   0.6.*T).*cosd(2.*(-D + F + Om) + M) + ...,
                    (    200           ).*cosd(2.*F + Om) + ...,
                    (    129 -   0.1.*T).*cosd(Mt + 2.*(F + Om)) + ...,
                    (    -95 +   0.3.*T).*cosd(-M + 2.*(-D + F + Om)) + ...,
                    (    -70           ).*cosd(2.*(-D + F) + Om) + ...,
                    (    -53           ).*cosd(-Mt + 2.*(F + Om)) + ...,
                    (    -33           ).*cosd(Mt + Om) + ...,
                    (     26           ).*cosd(-Mt + 2.*(D + F + Om)) + ...,
                    (     32           ).*cosd(-Mt + Om) + ...,
                    (     27           ).*cosd(Mt + 2.*F + Om) + ...,
                    (    -24           ).*cosd(Om + 2.*(-Mt + F)) + ...,
                    (     16           ).*cosd(2.*(D + F + Om)) + ...,
                    (     13           ).*cosd(2.*(Mt + F + Om)) + ...,
                    (    -12           ).*cosd(Mt + 2.*(-D + F + Om)) + ...,
                    (    -10           ).*cosd(-Mt + Om + 2.*F) + ...,
                    (     -8           ).*cosd(2.*D - Mt + Om) + ...,
                    (      7           ).*cosd(2.*(-D + M + F + Om)) + ...,
                    (      9           ).*cosd(M + Om) + ...,
                    (      7           ).*cosd(-2.*D + Mt + Om) + ...,
                    (      6           ).*cosd(-M + Om) + ...,
                    (      5           ).*cosd(2.*(D + F) - Mt + Om) + ...,
                    (      3           ).*cosd(2.*(D + F + Om) + Mt) + ...,
                    (     -3           ).*cosd(2.*(F + Om) + M) + ...,
                    (      3           ).*cosd(-M + 2.*(F + Om)) + ...,
                    (      3           ).*cosd(2.*(D + F) + Om) + ...,
                    (     -3           ).*cosd(2.*(-D + Mt + F + Om)) + ...,
                    (     -3           ).*cosd(2.*(-D + F) + Mt + Om) + ...,
                    (      3           ).*cosd(2.*(D - Mt) + Om) + ...,
                    (      3           ).*cosd(2.*D + Om) + ...,
                    (      3           ).*cosd(2.*(F - D) + Om - M) + ...,
                    (      3           ).*cosd(Om - 2.*D) + ...,
                    (      3           ).*cosd(2.*(Mt + F) + Om);

                nut_obli = DObl./10000/3600; % convert to degrees...
            end
        end
        
    end
    
    properties (Transient=true, Hidden=true, Dependent=true) % shortcuts to old property names, kept for backward compatibility
        
        DE;
        DE_deg;
        
        ECL_LAMBDA;
        ECL_BETA;
        
    end
    
    methods % shortcuts to old property names, kept for backward compatibility
        
        function val = get.DE(obj)
            
            val = obj.DEC;
            
        end
        
        function set.DE(obj, val)
            
            obj.DEC = val;
            
        end
        
        function val = get.DE_deg(obj)
            
            val = obj.DEC_deg;
            
        end
        
        function set.DE_deg(obj, val)
            
            obj.DEC_deg = val;
            
        end
        
        function val = get.ECL_LAMBDA(obj)
            
            val = obj.ECL_long;
            
        end
        
        function set.ECL_LAMBDA(obj, val)
            
            obj.ECL_long = val;
            
        end
        
        function val = get.ECL_BETA(obj)
            
            val = obj.ECL_lat;
            
        end
        
        function set.ECL_BETA(obj, val)
            
            obj.ECL_lat = val;
            
        end
        
    end
    
end

