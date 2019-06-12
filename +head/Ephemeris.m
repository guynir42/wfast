classdef (CaseInsensitiveProperties, TruncatedProperties) Ephemeris < handle

    properties % objects
        
        time@datetime;
        
    end
    
    properties % inputs/outputs
        
        RA_deg;
        DEC_deg;
        
    end
    
    properties % switches/controls
        
        debug_bit = 0;
        
    end
    
    properties(Dependent=true)
        
        RA;
        DEC;
        
        STARTTIME;
        JD;
        MJD;

        LST;
        LST_deg;
        
        HA;
        HA_deg;
        
        ALT_deg;
        
        AZ_deg;
        AIRMASS;
        
    end
    
    properties(Hidden=true)
        
        sidereal_time_type = 'mean'; % can also choose "apparent"
        
        latitude = 30.59678; % from my phone inside the dome (using whatsmylocation.net)
        longitude = 34.76202; % from my phone inside the dome (using whatsmylocation.net)
        
%         latitude = 30.5968322; % from my phone (outside the dome)
%         longitude = 34.7619663; % from my phone (outside the dome)
        
%         latitude = 30.59583333333; % WISE observatory (from wikipedia)
%         longitude = 34.76333333333; % WISE observatory (from wikipedia)
        
%         latitude = 30.597296; % updated from google maps
%         longitude = -34.761897; % updated from google maps
        
%         latitude = 31.907867; % of the observatory
%         longitude = 34.811363; % of the observatory
                
        version = 1.02;
        
    end
    
    methods % constructor
        
        function obj = Ephemeris(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'head.Ephemeris')
                if obj.debug_bit, fprintf('Ephemeris copy-constructor v%4.2f\n', obj.version); end
                util.oop.full_copy(varargin{1});
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
        
    end
    
    methods % getters
        
        function val = get.RA(obj)
            
            val = obj.deg2hour(obj.RA_deg);
            
        end
        
        function val = get.DEC(obj)
            
            val = obj.deg2sex(obj.DEC_deg);
            
        end
        
        function val = get.STARTTIME(obj)
            
            val = util.text.time2str(obj.time);
            
        end
        
        function val = get.JD(obj)

            if isempty(obj.STARTTIME)
                val = [];
            else
                val = juliandate(obj.STARTTIME, 'yyyy-mm-ddThh:MM:ss');
            end
             
        end
        
        function [ha, l] = getHAandLST(obj) % do we need this anymore?

            l = obj.local_sidereal_time(obj.JD, abs(obj.longitude)); % in fractions of day
            
            ha = 360*l - obj.RA; % in degrees
            
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
        
        function val = get.LST(obj) % in text format (HH:MM:SS)
            
            val = obj.deg2hour(obj.LST_deg);
            
        end
        
        function val = get.HA_deg(obj)
            
            val = obj.LST_deg - obj.RA_deg;
%             val = obj.getHAandLST;
            
        end
        
        function val = get.HA(obj)
            
            val = obj.deg2hour(obj.HA_deg);
            
        end
        
        function val = get.ALT_deg(obj)
            
%             val = obj.ha2alt(obj.HA_deg, obj.DEC_deg, obj.latitude);

            val = asind(sind(obj.DEC_deg).*sind(obj.latitude) + cosd(obj.DEC_deg).*cosd(obj.latitude).*cosd(obj.HA_deg));
            
        end
        
        function val = get.AZ_deg(obj)
            
%             val = obj.ha2az(obj.HA_deg, obj.DEC_deg, deg2rad(obj.latitude));
            
            SinAlt = sind(obj.DEC_deg).*sind(obj.latitude) + cosd(obj.DEC_deg).*cosd(obj.HA_deg).*cosd(obj.latitude);
            CosAlt = sqrt(1-SinAlt.*SinAlt);

            SinAz  = (-cosd(obj.DEC_deg).*sind(obj.HA_deg))./CosAlt;
            CosAz  = (sind(obj.DEC_deg).*cosd(obj.latitude) - cosd(obj.DEC_deg).*cosd(obj.HA_deg).*sind(obj.latitude))./CosAlt;

            val = atan2d(SinAz, CosAz);
            
        end
        
        function val = get.AIRMASS(obj)
           
%             [~, val] = obj.ha2alt(obj.HA_deg, obj.DE_deg, obj.latitude);
            if obj.ALT_deg>87 % unreliable values near the horizon
                val = []; 
            else
                
                SecZ = 1./cosd(90-obj.ALT_deg);

                val = SecZ - 0.0018167.*(SecZ - 1) - 0.002875.*(SecZ - 1).*(SecZ - 1)...
                    - 0.0008083.*(SecZ - 1).*(SecZ - 1).*(SecZ - 1); % Hardie's polynomial formula

            end
            
        end
        
    end
    
    methods % setters
        
        function set.RA(obj, val)

            if isempty(val)
                obj.RA_deg = [];
            elseif isnumeric(val) && isscalar(val)
                obj.RA_deg = val*15;
            elseif isnumeric(val) && length(val)==3
                obj.RA_deg = 0;
            elseif ischar(val)
                obj.RA_deg = obj.hour2deg(val);
            else
                error('Input a string RA in HH:MM:SS format!');
            end
            
%             if isnumeric(val) && isvector(val) && length(val)==3
%                 obj.RA = sprintf('%02d:%02d:%05.2f', val(1), val(2), val(3));
%             elseif isnumeric(val) && isscalar(val)
%                 obj.RA = obj.deg2hour(val);
%             elseif ischar(val)
%                 obj.RA = val;
%             else
%                 error('Unknown format for RA');
%             end
            
        end 
        
        function set.DEC(obj, val)
            
            if isempty(val)
                obj.DEC_deg = [];
            elseif ischar(val)
                obj.DEC_deg = obj.sex2deg(val);
            else
                error('Input a string DEC in DD:MM:SS format!');
            end
            
%             if isnumeric(val) && isvector(val) && length(val)==3
%                 obj.DEC = sprintf('%+03d:%02d:%04.1f', val(1), val(2), val(3));
%             elseif isnumeric(val) && isscalar(val)
%                 obj.DEC = obj.deg2sex(val);
%                 obj.DEC_deg = val;
%             elseif ischar(val)
%                 obj.DEC = val;
%                 obj.DEC_deg = obj.hour2deg(val);
%             else
%                 error('Unknown format for DEC');
%             end
            
        end
        
    end
    
    methods % calculations
        
        function input(obj, RA, DEC, time)
            
            if nargin<3
                disp('Usage: input(RA,DEC,time)');
            end
            
            obj.RA = RA;
            obj.DEC = DEC;
            
            if nargin>3 && ~isempty(time)
                
                if isa(time, 'datetime')
                    obj.time = time;
                elseif ischar(time) && util.text.cs(time, 'now', 'update')
                    obj.update;
                elseif ischar(time)
                    obj.time = util.text.str2time(time);
                elseif isnumeric(time) % juldate??
                    obj.time = datetime(time, 'ConvertFrom', 'juliandate', 'TimeZone', 'UTC');
                end
                
            end
            
        end
        
        function update(obj)
            
            obj.time = datetime('now', 'timezone', 'UTC');
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end
    
    methods (Static=true)
        
        function str = deg2hour(number) % convert the degrees into hours then into HH:MM:SS string
            
            str = head.Ephemeris.numbers2hour(number/15); 
            
        end
        
        function str = deg2sex(number) % convert degrees into DD:MM:SS string
            
            str = head.Ephemeris.numbers2sex(number);
            
        end
        
        function val = hour2deg(str) % convert an HH:MM:DD hour string to degrees
            
            num = util.text.extract_numbers(str);
            
            if isempty(num)
                val = [];
            else
            
                hours = num{1}(1);
                degrees = hours*15;
                
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
                
            end
            
        end
        
        function val = sex2deg(str) % convert a DD:MM:SS string (sexagesimal) to degrees
            
            num = util.text.extract_numbers(str);
            
            if isempty(num)
                val = [];
            else
            
                degrees = num{1}(1);
                s = sign(degrees);
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
        
        function str = numbers2hour(hours, minutes, seconds)
        % convert a number of hours (or hours, minutes, seconds) into a 
        % string of right ascention. 
        
            if nargin==0, help('head.Ephemeris.num2ra'); return; end
        
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
            
            total_secs = hours*3600 + minutes*60 + seconds;
            
            hours   = floor(total_secs/3600);
            minutes = floor(mod(total_secs,3600)/60);
            seconds = mod(total_secs,60);
            
            str = sprintf('%02d:%02d:%04.1f', hours, minutes, seconds);
            
        end
        
        function str = numbers2sex(degrees, minutes, seconds) 
        % convert a number of degrees (or degrees, arcminutes, arcseconds) into a 
        % string of sexidecimal degrees
        
            if nargin==0, help('head.Ephemeris.num2dec'); return; end
            
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
            
            if s>0
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
    
    properties (Transient=true, Hidden=true, Dependent=true)
        
        DE;
        DE_deg;
        
    end
    
    methods
        
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
        
    end
    
end

