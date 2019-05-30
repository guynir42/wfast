classdef Ephemeris < handle

    properties % objects
        
%         pars@head.Parameters;        
        time@datetime;
        
    end
    
    properties % inputs/outputs
        
        RA = '';
        DE = '';
                
    end
    
    properties % switches/controls
        
        debug_bit = 0;
        
    end
    
    properties(Dependent=true)
        
        juldate;
        HA;
        LST;
        ALT;
        AZ;
        airmass;
        
    end
    
    properties(Hidden=true)
        
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
                
        version = 1.00;
        
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
                    
%                     if isa(varargin{ii}, 'head.Parameters')
%                         obj.pars = varargin{ii};
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
        
        function val = file_start_datestr(obj)
            
            val = util.text.time2str(obj.time);
            
        end
        
        function val = RA_rad(obj)
           
            if ischar(obj.RA)
                val = obj.ra2rad(obj.RA);
            else
                val = obj.RA;
            end
            
        end
        
        function val = DE_rad(obj)
           
            if ischar(obj.DE)
                val = obj.dec2rad(obj.DE);
            else
                val = obj.DE;
            end
            
        end
        
        function val = RA_deg(obj)
           
            val = 360/2/pi*obj.RA_rad;
            
        end
        
        function val = DE_deg(obj)
           
            val = 360/2/pi*obj.DE_rad;
            
        end
        
        function jd = get.juldate(obj)

            if isempty(obj.file_start_datestr)
                jd = [];
            else
                jd = juliandate(obj.file_start_datestr, 'yyyy-mm-ddThh:MM:ss');
            end
             
        end
        
        function [ha, l] = getHAandLST(obj)

            l = obj.lst(obj.juldate, deg2rad(abs(obj.longitude))); % in fractions of day
            
            ha = 2*pi*l - obj.RA_rad; % in radians
            
        end
        
        function val = HA_rad(obj)
            
            val = obj.getHAandLST;
            
        end
        
        function val = HA_deg(obj)
            
            val = 360/2/pi*obj.HA_rad;
            
        end
        
        function ha = get.HA(obj)
            
            ha = obj.rad2ra(obj.HA_rad);
            
        end
        
        function val = LST_rad(obj)
            
            [~, val] = obj.getHAandLST;
            
        end
        
        function val = LST_deg(obj)
            
            val = 360/2/pi*obj.LST_rad;
            
        end
        
        function lst = get.LST(obj)
            
            lst = obj.rad2ra(2*pi*obj.LST_rad); % in text format (HH MM SS)
            
        end
                
        function val = get.ALT(obj)
            
            val = obj.ha2alt(obj.HA_rad, obj.DE_rad, deg2rad(obj.latitude));
            
        end
        
        function val = ALT_deg(obj)
            
            val = 360/2/pi.*obj.ALT;
            
        end
        
        function val = get.AZ(obj)
            
            val = obj.ha2az(obj.HA_rad, obj.DE_rad, deg2rad(obj.latitude));
            
        end
        
        function a = get.airmass(obj)
           
            [~, a] = obj.ha2alt(obj.HA_rad, obj.DE_rad, degtorad(obj.latitude));
            
            if isnan(a)
                a = [];
            end
            
        end
        
    end
    
    methods % setters
        
        function set.RA(obj, val)
           
            if isnumeric(val) && isvector(val) && length(val)==3
                obj.RA = sprintf('%02d %02d %05.2f', val(1), val(2), val(3));
            elseif isnumeric(val) && isscalar(val)
                obj.RA = obj.rad2ra(val);
            else
                obj.RA = val;
            end
            
        end 
        
        function set.DE(obj, val)
           
            if isnumeric(val) && isvector(val) && length(val)==3
                obj.DE = sprintf('%+03d %02d %04.1f', val(1), val(2), val(3));
            elseif isnumeric(val) && isscalar(val)
                obj.DE = obj.rad2dec(val);
            else
                obj.DE = val;
            end
            
        end
    end
    
    methods % calculations
        
        function input(obj, RA,DE,time)
            
            if nargin<3
                disp('Usage: input(RA,DE,time)');
            end
            
            if ischar(RA)
                obj.RA = RA;
            elseif isnumeric(RA)
                obj.RA_rad = RA;
            end
            
            if ischar(DE)
                obj.DE = DE;
            elseif isnumeric(DE)
                obj.DE_rad = DE;
            end
            
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
        
        function [rads, arcsecs] = dec2rad(Dec_string)
            % convert declination string to radians (or arcsecs)
            % usage: [rads, arcsecs] = dec2rad(Dec_string)
            
            if nargin==0
                help('head.Ephemeris.dec2rad');
                return;
            end
            
            if isempty(Dec_string)
                rads = [];
                arcsecs = [];
                return;
            end
            
            if ischar(Dec_string)
                Dec_string = {Dec_string};
            end
            
            factors = [3600 60 1];
            
            numbers = util.text.extract_numbers(Dec_string);
            
            arcsecs = zeros(length(numbers),1);
            
            for ii = 1:length(numbers)
                
                sign = 1;
                
                for jj = 1:length(Dec_string{ii}) % pick up minus sign before first numeral if it is zero...
                    
                    if Dec_string{ii}(jj)=='0' && jj>1
                        if strcmp(Dec_string{ii}(jj-1), '-')
                            sign = -1;
                            break;
                        end
                    end
                end
                
                for jj = 1:min([3 length(numbers{ii})])
                    arcsecs(ii) = arcsecs(ii) + numbers{ii}(jj)*factors(jj);
                end
                
                arcsecs(ii) = arcsecs(ii)*sign;
                
            end
            
            rads = arcsecs./head.Ephemeris.rad2arcsec;
            
        end
        
        function [rads, arcsecs] = ra2rad(RA_string)
            % convert right ascention string to radians (or arcsecs)
            % usage: [rads, arcsecs] = ra2rad(RA_string)
            
            if nargin==0
                help('head.Ephemeris.ra2rad');
                return;
            end
            
            if isempty(RA_string)
                rads = [];
                arcsecs = [];
                return;
            end
            
            if ischar(RA_string)
                RA_string = {RA_string};
            end
            
            factors = 15.*[3600 60 1];
            
            numbers = util.text.extract_numbers(RA_string);
            
            arcsecs = zeros(length(numbers),1);
            
            for ii = 1:length(numbers)
                
                sign = 1;
                
                for jj = 1:length(RA_string{ii}) % pick up minus sign before first numeral if it is zero...
                    
                    if RA_string{ii}(jj)=='0' && jj>1
                        if strcmp(RA_string{ii}(jj-1), '-')
                            sign = -1;
                            break;
                        end
                    end
                end
                
                for jj = 1:min([3 length(numbers{ii})])
                    arcsecs(ii) = arcsecs(ii) + numbers{ii}(jj)*factors(jj);
                end
                
                arcsecs(ii) = arcsecs(ii)*sign;
                
            end
            
            rads = arcsecs./head.Ephemeris.rad2arcsec;
            
        end
        
        function val = rad2arcsec(rad)
            
            if nargin<1
                rad = 1;
            end
            
            val = rad*360/(2*pi)*3600;
            
        end
        
        function str = num2ra(hours, minutes, seconds)
        % convert a number of hours (or hours, minutes, seconds) into a 
        % string of right ascention. 
        
            if nargin==0, help('head.Ephemeris.num2ra'); return; end
        
            if isempty(hours)
                str = [];
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
            
            str = sprintf('%02d %02d %04.1f', hours, minutes, seconds);
            
        end
        
        function str = num2dec(degrees, minutes, seconds)
        % convert a number of degrees (or degrees, arcminutes, arcseconds) into a 
        % string of declination
        
            if nargin==0, help('head.Ephemeris.num2dec'); return; end
            
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
                
            str = sprintf('%s%02d %02d %04.1f', sign_str, degrees, minutes, seconds);
            
        end
        
        function str = rad2dec(radians)
            
            degrees = radians*360/(2*pi);
            
            str = head.Ephemeris.num2dec(degrees);
            
        end
        
        function str = rad2ra(radians)
            
            hours = radians*24/(2*pi);
            
            str = head.Ephemeris.num2ra(hours);
            
        end

        function LST = lst(JD,EastLong,STType)
        %--------------------------------------------------------------------------
        % lst function                                                       ephem
        % Description: Local Sidereal Time, (mean or apparent), for vector of
        %              JDs and a given East Longitude.
        % Input  : - Vector of JD [days], in UT1 time scale.
        %          - East Longitude in radians.
        %          - Sidereal Time Type,
        %            'm' - Mean (default).
        %            'a' - apparent.
        % Output : - vector of LST in fraction of day.
        % Tested : Matlab 5.3
        %     By : Eran O. Ofek                    Aug 1999
        %    URL : http://weizmann.ac.il/home/eofek/matlab/
        % Example: LST=lst(2451545+[0:1:5]',0);  % LST at Greenwhich 0 UT1
        % Reliable: 1
        %--------------------------------------------------------------------------

            RAD = 180./pi;

            if nargin==0, help('head.Ephemeris.lst'); return; end
            
            if (nargin==2)
               STType = 'm';
            elseif (nargin==3)
               % do nothing
            else
               error('Illigal number of input arguments');
            end

            % convert JD to integer day + fraction of day
            TJD = floor(JD - 0.5) + 0.5;
            DayFrac = JD - TJD;

            T = (TJD - 2451545.0)./36525.0;

            GMST0UT = 24110.54841 + 8640184.812866.*T + 0.093104.*T.*T - 6.2e-6.*T.*T.*T;

            % convert to fraction of day in range [0 1)
            GMST0UT = GMST0UT./86400.0;

            GMST0UT = GMST0UT - floor(GMST0UT);
            LST = GMST0UT + 1.0027379093.*DayFrac + EastLong./(2.*pi);
            LST = LST - floor(LST);


            switch STType
             case {'m'}
                % do nothing
             case {'a'}
                % calculate nutation
                NutMat = ephem.nutation(JD);
                Obl    = ephem.obliquity(JD);
                EquationOfEquinox = (RAD.*3600).*NutMat(:,1).*cos(Obl)./15;
                LST = LST + EquationOfEquinox./86400;    
             otherwise
                error('Unknown sidereal time type');
            end
            
        end
        
        function [Alt, AM] = ha2alt(HA,Dec,Lat)
        %--------------------------------------------------------------------------
        % ha2alt function                                                    ephem
        % Description: Given Hour Angle as measured from the meridian, the source
        %              declination and the observer Geodetic latitude, calculate
        %              the source altitude above the horizon and its airmass.
        % Input  : - Hour Angle [radians].
        %          - Declination [radians].
        %          - Latitude [radians].
        % Output : - Altitude [radians].
        %          - Airmass.
        % See also: horiz_coo.m, ha2az.m
        % Tested : Matlab 7.10
        %     By : Eran O. Ofek                    Aug 2010
        %    URL : http://weizmann.ac.il/home/eofek/matlab/
        % Reliable: 1
        %--------------------------------------------------------------------------

            if nargin==0, help('head.Ephemeris.ha2alt'); return; end
        
            Alt = asin(sin(Dec).*sin(Lat) + cos(Dec).*cos(Lat).*cos(HA));
            AM  = head.Ephemeris.hardie(pi./2-Alt);

        end
        
        function [Az,Alt,AM]=ha2az(HA,Dec,Lat)
        %--------------------------------------------------------------------------
        % ha2az function                                                     ephem
        % Description: Given Hour Angle as measured from the meridian, the source
        %              declination and the observer Geodetic latitude, calculate
        %              the horizonal source azimuth
        % Input  : - Hour Angle [radians].
        %          - Declination [radians].
        %          - Latitude [radians].
        % Output : - Azimuth [radians].
        %          - Altitude [radians].
        %          - Airmass.
        % See also: horiz_coo.m, ha2alt.m
        % Tested : Matlab 7.10
        %     By : Eran O. Ofek                    Aug 2010
        %    URL : http://weizmann.ac.il/home/eofek/matlab/
        % Example: [Az,Alt,AM]=ha2az(1,1,1)
        % Reliable: 1
        %--------------------------------------------------------------------------
            
            if nargin==0, help('head.Ephemeris.ha2az'); return; end

            SinAlt = sin(Dec).*sin(Lat) + cos(Dec).*cos(HA).*cos(Lat);
            CosAlt = sqrt(1-SinAlt.*SinAlt);

            SinAz  = (-cos(Dec).*sin(HA))./CosAlt;
            CosAz  = (sin(Dec).*cos(Lat) - cos(Dec).*cos(HA).*sin(Lat))./CosAlt;

            Az     = atan2(SinAz, CosAz);
            if (nargin>1)
                Alt = asin(sin(Dec).*sin(Lat) + cos(Dec).*cos(Lat).*cos(HA));
                if (nargin>2)
                    AM  = head.Ephemeris.hardie(pi./2-Alt);
                end
            end

        end
        
        function AM = hardie(X,Algo)
        %--------------------------------------------------------------------------
        % hardie function                                                    ephem
        % Description: Calculate airmass using the Hardie formula.
        % Input  : - Matrix of zenith distances [radians].
        %          - Algorith: {'hardie','csc'}. Default is 'hardie'.
        % Output : - Air mass.
        % Tested : Matlab 4.2
        %     By : Eran O. Ofek                    Jan 1994
        %    URL : http://weizmann.ac.il/home/eofek/matlab/
        % See also: airmass.m; hardie_inv.m
        % Example: AM=hardie([1;1.1]);
        % Reliable: 1
        %--------------------------------------------------------------------------

            RAD = 180./pi;

            Def.Algo = 'hardie'; 
            if nargin==1
               Algo = Def.Algo;
            elseif nargin==2
               % do nothing
            else
               error('Illegal number of input arguments');
            end

            SecZ = 1./cos(X);
            switch lower(Algo)
                case 'hardie'
                    AM   = SecZ - 0.0018167.*(SecZ - 1) - 0.002875.*(SecZ - 1).*(SecZ - 1) - 0.0008083.*(SecZ - 1).*(SecZ - 1).*(SecZ - 1);
                case 'csc'
                    AM = SecZ;
                otherwise
                    error('Unknown airmass algorith option');
            end

            FlagI = X>87./RAD;
            AM(FlagI) = NaN;

        end
        
    end
    
end

