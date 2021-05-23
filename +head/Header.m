classdef (CaseInsensitiveProperties, TruncatedProperties) Header < dynamicprops
% This class contains all the useful information to be accessed by user
% or saved as metadata for a specific observation. 
% 
% The properties of this object are often saved as keywords to files, 
% either as HDF5 or as text files. 
%
% A handle to this object should be given to any pipeline objects that 
% need access to metadata (e.g., exposure time, plate scale or sky coordinates). 
% Objects that calculate additional info can save it back into the header. 
%
% The header also contains a few objects (e.g., Ephemeris) that contain some
% of the information or do some calculations. See the documentation for those
% classes for more info. 
% 
% The properties can be called by shorter expressions: e.g., call PROJ instead
% of PROJECT. Also you can use lower or upper case. 
% Some keywords (property names) are accessible with older names for backward
% compatibility (see the last two blocks of properties/methods).
%


    properties(Transient=true) 
        
        gui@head.gui.HeaderGUI; % GUI objects are never saved, and generated on demand
        
    end

    properties % objects 
        
        filter_obj; % must be a head.Filter (enforced in the setter...)
        ephem@head.Ephemeris; % all sky coordinates are managed through this object
        stars@head.Star; % for keeping track of position/magnitude of a small number of stars
        
        WCS@head.WorldCoordinates; % transformation between sky and pixel coordinates. Reorganized version of Eran's WCS object
        
        % NOTE: do not share the handles to objects in this class. They are
        % often replaced to other handles when loading Header from file. 
        
    end
    
    properties 
        
%         OBJECT = 'star1'; % the name of the target field or object
        TYPE; % type of observation... choose science/dark/flat etc... 
        
        COMMENT = ''; % general comment
        
        PROJECT = 'WFAST'; % project can be Kraar or WFAST
        INST = 'Zyla_5.5'; % The camera model used. In new headers this must be updated by hardware to Balor
        
        STARTTIME; % start time of current batch of images (in YYYY-MM-DDThh:mm:ss.sss format) in UTC!
        ENDTIME; % end time of current batch of images (in YYYY-MM-DDThh:mm:ss.sss format) in UTC!
        RUNSTART;  % start time of the observation run (in YYYY-MM-DDThh:mm:ss.sss format) in UTC!
        
        END_STAMP; % the camera timestamp correspinding to ENDTIME
        
        TEL_APER = 57; % telescope aperture in cm 
        FOCLEN = 108; % telescope focal length in cm
        F_RATIO; % focal ratio is FOCLEN/TEL_APER (updated automatically)
        
        PIXSIZE = 6.5; % size of pixels (microns) 
        SCALE; % plate scale in arcseconds per pixel
        
        EXPTIME; % exposure time (seconds)
        FRAMERATE; % frame rate (Hz). If NaN or empty, camera tried to take images as fast as possible with current EXPTIME
        ACT_FRAMERATE; % actual frame rate measured while taking images
        
        THRESH_DETECTION; % detection threshold for findStars
        THRESH_STACK; % detection threshold for stars on the stack (for calculating limiting magnitude)
        THRESH_INDIVIDUAL; % threshold for using this star in the event finding algorithm
        
        LIMMAG_DETECTION; % the GAIA_BP magnitude of the faintest star from findStars
        LIMMAG_STACK; % the limiting magnitude for a stack image 
        LIMMAG_INDIVIDUAL; % the limiting magnitude for individual frames at 25 Hz
        ZEROPOINT; % the conversion between flux and magnitude zero
        BACKGROUND; % the estimated number of counts per pixel from the background (per frame, not stacked!)
        
        FILTER = 'F505W'; % optical filter installed on telescope
        
        NAXIS; % number of axes for the full frame images
        NAXIS1; % size of the y axis of the images
        NAXIS2; % size of the x axis of the images
        NAXIS3; % number of images per batch
        NAXIS4; % not used right now
        
        BINX = 1; % binning factor in x. Not currently used
        BINY = 1; % binning factor in y. Not currently used
        ROI; % region of interest: a 4-vector with [top left height width] in pixels
        CCDSEC; % which CCD sector was used (not in use right now)
        
        GAIN; % the camera gain in photons per ADU (0.6 for Zyla and 0.8 for Balor)
        READNOISE; % estimated read noise per pixel (counts per pixel per frame)
        DARKCUR; % estimated dark current (counts per pixel per second)
        
        IS_DARK = 0; % set this to 1 for dark exposures
        IS_FLAT = 0; % set this to 1 for flat exposures
        IS_SIMULATED = 0; % set this to 1 for simulation images
       
        FOCUS_POS; % position of focuser (in mm, the average of all actuators)
        FOCUS_TIP; % value of focus spider tip
        FOCUS_TILT; % value of focus spider tilt
         
        SEEING; % arcsec FWHM of stars
        TEMP_DET; % detector temperature (celsius)
        TEMP_OUT; % outside temperature (celsius)
        
        WIND_DIR; % in km/h
        WIND_SPEED; % in degrees
        HUMID_IN; % humidity inside dome (percent)
        HUMID_OUT; % humidity outside dome (percent)
        PRESSURE; % pressure in mbar
        LIGHT; % light value in arbitrary units (from Boltwood, max is 1024, dark is around 100)

        SENSOR_TEMP; % temperature in C of the sensor
        
        % coordinates of the object as given by user/scheduler
        OBJRA; 
        OBJDEC;
        OBJRA_DEG;
        OBJDEC_DEG;
        OBJMAG; % keep track of the magnitude of the object, if known 
        
        % hardware (mount) coordinates
        TELRA; % telescope RA reported from hardware in sexagesimal hour string
        TELDEC; % telescope Dec reported from hardware in sexagesimal degree string
        TELRA_DEG; % telescope RA reported from hardware in numeric degrees
        TELDEC_DEG; % telescope Dec reported from hardware in numeric degrees
        
    end
    
    properties(Dependent=true)
        
        % ephemeris: time & coordinates
        OBJECT; 
        
        % measured from observations (the data for these is in the WCS object)
        OBSRA; % center of field RA as measured by fitting to GAIA, in sexagesimal hour string
        OBSDEC; % center of field Dec as measured by fitting to GAIA, in sexagesimal degree string
        OBSRA_DEG; % center of field RA as measured by fitting to GAIA, in numeric degrees
        OBSDEC_DEG; % center of field Dec as measured by fitting to GAIA, in numeric degrees
        FIELDROT; % rotation of field vs. celestial north (in degrees)
        
        % the following properties are stored in the ephem object
        RA; % this is the RA given to telescope by scheduler/user, in sexagesimal hour string
        DEC; % this is the Dec given to telescope by scheduler/user, in sexagesimal degree string
        RA_DEG; % this is the RA given to telescope by scheduler/user, in numeric degrees
        DEC_DEG; % this is the Dec given to telescope by scheduler/user, in numeric degrees
        FIELD_ID; % numeric identifier for fields in a target list/bank
        
        HA; % hour angle, in sexagesimal hour string
        HA_DEG; % hour angle, in numeric degrees
        LST; % local sidereal time, in sexagesimal hour string
        ALT; % altitude in numeric degrees
        AZ; % azimuth in numeric degrees
        AIRMASS; % length of atmosphere in the direction of center of field
        
        MOONAZ; % azimuth of moon to observer, in numeric degrees
        MOONALT; % altitude of moon above horizon, in numeric degrees
        MOONILL; % illumination fraction of the moon
        MOONDIST; % distance between center of field and moon, in numeric degrees

        SUNAZ; % azimuth of sun to observer, in numeric degrees
        SUNALT; % altitude of sum above horizon, in numeric degrees
        
        OBSLONG; % observatory longitude on the Earth, in numeric degrees
        OBSLAT; % observatory latitude on the Earth, in numeric degrees
        OBSEL; % observatory elevation above sea level, in meters
         
        JD; % julian day of start of observation
        MJD; % modified julian day
        
        % filter properties from the filter_obj
        FILT_WAVE; % central wavelength of the filter, in nano-meters
        FILT_WIDTH; % bandwidth of the filter, in nano-meters

    end
    
    properties(Hidden=true) 
       
        run_identifier = '';  
        
        DEADTIME = 0; % time lost after each frame, in seconds (needs to be updated someday)
        
        QE = 1; % quantum efficiency (needs to be updated someday)
        
        % these are defaults used for recovering values when using the GUI
        default_APERTURE; 
        default_F_RATIO;
        default_FILTER;
        default_OBJECT;
        
        filter_name_full; % for backward compatibility with older versions.
        
        CAMS_VER; % camera firmware version
        TELS_VER; % telescope firmware version
        
        PHOT_PARS; % a structure with the parameters used for photometry
        
        debug_bit = 1; 
        version = 4.00;
          
    end
    
    methods % constructor
       
        function obj = Header(other)
            
            if nargin>0 && isa(other, 'head.Parameters')
                
                if obj.debug_bit>1, fprintf('Header conversion-constructor v%4.2f\n', obj.version); end
                
                obj = cast(util.oop.full_copy(other));
                
            elseif nargin>0 && isa(other, 'head.Header')
                
                if obj.debug_bit>1, fprintf('Header copy-constructor v%4.2f\n', obj.version); end
                
                obj = util.oop.full_copy(other);
                    
            else
            
                if obj.debug_bit>1, fprintf('Header constructor v%4.2f\n', obj.version); end
                
                obj.filter_obj = head.Filter(obj.FILTER);
                obj.ephem = head.Ephemeris;
                obj.WCS = head.WorldCoordinates;
%                 obj.cat = head.Catalog(obj);
                
                obj.FOCLEN = obj.FOCLEN; % this should fill the values for SCALE and F_RATIO
                
                util.oop.save_defaults(obj);
                
                obj.update; % set the time to now. 
                
            end
            
        end
                
    end
    
    methods % reset methods
        
        function reset(obj)
            
            obj.STARTTIME = '';
            obj.ENDTIME = '';
            obj.END_STAMP = []; 
            
            obj.SEEING = [];
            obj.TEMP_DET = [];
            obj.TEMP_OUT = [];
            obj.WIND_DIR = [];
            obj.WIND_SPEED = [];
            obj.HUMID_IN = [];
            obj.HUMID_OUT = [];
            obj.PRESSURE = [];
            obj.LIGHT = [];
            obj.TELRA = [];
            obj.TELDEC = [];
            obj.TELRA_DEG = [];
            obj.TELDEC_DEG = [];
            
        end
        
        function clearStars(obj) % remove all stars from the stars vector
           
            obj.stars = head.Star.empty;
            
        end
        
        function clearObsConditions(obj) % remove all seeing and weather measurements
           
            obj.SEEING = [];
            obj.TEMP_DET = [];
            obj.TEMP_OUT = [];
            obj.WIND_DIR = [];
            obj.WIND_SPEED = [];
            obj.HUMID_IN = [];
            obj.HUMID_OUT = [];
            obj.PRESSURE = [];
            obj.LIGHT = [];
        
        end
                
        function resetTarget(obj) % reset the target name, the ephem object and all coordinates of the target
            
            obj.OBJECT = obj.default_OBJECT;
            obj.ephem.reset;
            
        end
        
        function resetDrifts(obj) % reset drifts in the star objects
           
            util.oop.setprop(obj.stars, 'drift_x', []);
            util.oop.setprop(obj.stars, 'drift_y', []);
            
        end
        
    end
    
    methods % getters
        
        function val = get.filter_obj(obj)
            
            if ~isempty(obj.filter_obj)
                val = obj.filter_obj;
            else

                try 
                    obj.filter_obj = head.Filter(obj.FILTER);
                catch 
                    obj.filter_obj = head.Filter.empty;
                end

            end
            
        end
        
        function val = filter_range(obj)
            
            val = obj.FILT_WAVE + [-0.5 0.5].*obj.FILT_WIDTH;
            
        end
        
        function val = diff_limit(obj)
            
            if isempty(obj.aperture) || isempty(obj.wavelength)
                val = [];
            else
                val = obj.wavelength*1e-7./obj.APERTURE*360*3600/2/pi;
            end
            
        end
        
        function val = diff_limit_pix(obj)
           
            if isempty(obj.SCALE) || isempty(obj.diff_limit)
                val = [];
            else            
                val = obj.diff_limit./obj.SCALE;
            end
            
        end
        
        function val = aperture_area(obj)
            
%             if isempty(obj.slit_size)
                val = pi*(obj.APERTURE/2).^2;
%             else
%                 val = obj.aperture*obj.slit_size;
%             end
            
        end
        
        function val = image_sampling(obj)
            
            if isempty(obj.diff_limit) || isempty(obj.plate_scale)
                val = [];
            else
                val = obj.diff_limit/obj.plate_scale;
            end
            
        end
                
        % camera and filter 
        function str = cam_name(obj)
           
            c = strsplit(obj.INST, {' ','_'});
            
            str = c{1};
            
        end
        
        function val = get.FILTER(obj)
        
            if ~isempty(obj.FILTER) && util.text.cs('[1x1 head.Filter]', obj.FILTER) % patch for backward compatibility with older headers where "filter" was the object
                obj.FILTER = obj.filter_name_full;
            end
                
            val = obj.FILTER;

        end
     
        function val = get.FILT_WAVE(obj)
            
            if isempty(obj.filter_obj)
                val = [];
            else
                val = obj.filter_obj.wavelength;
            end
            
        end
        
        function val = get.FILT_WIDTH(obj)
            
            if isempty(obj.filter_obj)
                val = [];
            else
                val = obj.filter_obj.bandwidth;
            end
            
        end
        
        % ephemeris
        function val = get.OBJECT(obj)
            
            if isempty(obj.ephem)
                val = '';
            else
                val = obj.ephem.name;
            end
            
        end
        
        function val = get.JD(obj)
            
            val = obj.ephem.JD;
            
        end
        
        function val = get.MJD(obj)
            
            val = obj.ephem.MJD;
            
        end
        
        function val = get.OBSRA(obj)
            
            val = head.Ephemeris.deg2hour(obj.OBSRA_DEG);
            
        end
        
        function val = get.OBSDEC(obj)
            
            val = head.Ephemeris.deg2sex(obj.OBSDEC_DEG);
            
        end
        
        function val = get.OBSRA_DEG(obj)
            
            if isempty(obj.WCS)
                val = [];
            else
                val = obj.WCS.RA_deg_center;
            end
            
        end
        
        function val = get.OBSDEC_DEG(obj)
            
            if isempty(obj.WCS)
                val = [];
            else
                val = obj.WCS.DE_deg_center;
            end
            
        end
        
        function val = get.FIELDROT(obj)
            
            
            if isempty(obj.WCS)
                val = [];
            else
                val = obj.WCS.rotation;
            end
            
        end
        
        function val = get.RA(obj)
            
            val = obj.ephem.RA;
            
        end
        
        function val = get.DEC(obj)
            
            val = obj.ephem.DEC;
            
        end
        
        function val = get.RA_DEG(obj)
            
            val = obj.ephem.RA_deg;
            
        end
        
        function val = get.DEC_DEG(obj)
            
            val = obj.ephem.DEC_deg;
            
        end
        
        function val = get.FIELD_ID(obj)
            
            val = obj.ephem.field_id; 
            
        end
        
        function val = get.HA(obj)
            
            val = obj.ephem.HA;
            
        end
        
        function val = get.HA_DEG(obj)
            
            val = obj.ephem.HA_deg;
            
        end
        
        function val = get.LST(obj)
            
            val = obj.ephem.LST;
            
        end
        
        function val = get.ALT(obj)
            
            val = obj.ephem.ALT_deg;
            
        end
        
        function val = get.AZ(obj)
            
            val = obj.ephem.AZ_deg;
            
        end
        
        function val = get.AIRMASS(obj)
            
            val = obj.ephem.airmass;
            
        end        
        
        function val = get.MOONAZ(obj)
            
            if isempty(obj.ephem) || isempty(obj.ephem.moon)
                val = [];
            else
                val = obj.ephem.moon.Az;
            end
            
        end
        
        function val = get.MOONALT(obj)
            
            if isempty(obj.ephem) || isempty(obj.ephem.moon)
                val = [];
            else
                val = obj.ephem.moon.Alt;
            end
            
        end
        
        function val = get.MOONILL(obj)
            
            if isempty(obj.ephem) || isempty(obj.ephem.moon)
                val = [];
            else
                val = obj.ephem.moon.IllF;
            end
            
        end
        
        function val = get.MOONDIST(obj)
            
            if isempty(obj.ephem) || isempty(obj.ephem.moon) || isempty(obj.ephem.moon.Dist)
                val = [];
            else
                val = obj.ephem.moon.Dist;
            end
            
        end
        
        function val = get.SUNAZ(obj)
            
            if isempty(obj.ephem) || isempty(obj.ephem.sun)
                val = [];
            else
                val = obj.ephem.sun.Az;
            end
            
        end
        
        function val = get.SUNALT(obj)
            
            if isempty(obj.ephem) || isempty(obj.ephem.sun)
                val = [];
            else
                val = obj.ephem.sun.Alt;
            end
            
        end
        
        % stars
        function val = magnitude(obj)
           
            if isempty(obj.stars)
                val = [];
                return;
            end
            
            val = [obj.stars.mag];
                        
        end
        
        function val = separation(obj)
                        
            if isempty(obj.stars)
                val = [];
                return;
            end
                  
            val = zeros(size(obj.stars));
            
            for ii = 1:length(val)
                temp = obj.stars(ii).getSeparationAndAngle(obj.stars(1), 'arcsec');
                if isempty(temp)
                    val(ii) = NaN;
                else
                    val(ii) = temp;
                end
            end
            
        end
        
        function val = pos_angle(obj)
                        
            if isempty(obj.stars)
                val = [];
                return;
            end
                  
            val = zeros(size(obj.stars));
            
            for ii = 1:length(val)
                [~, temp] = obj.stars(ii).getSeparationAndAngle(obj.stars(1), '', 'degrees');
                if isempty(temp)
                    val(ii) = NaN;
                else
                    val(ii) = temp;
                end
            end
            
        end
                
        function val = star_x(obj)
                        
            if isempty(obj.stars)
                val = [];
                return;
            end
                  
            val = zeros(size(obj.stars));
            
            for ii = 1:length(val)
                val(ii) = obj.stars(ii).getFinalX('pixels');
            end
            
            
        end
                
        function val = star_y(obj)
                        
            if isempty(obj.stars)
                val = [];
                return;
            end
                  
            val = zeros(size(obj.stars));
            
            for ii = 1:length(val)
                val(ii) = obj.stars(ii).getFinalX('pixels');
            end
            
            
        end
        
        function str = getApertureString(obj)
           
%            if isempty(obj.slit_size)
               str = sprintf('%dcm', round(obj.aperture));
%            else
%                str = sprintf('%dx%dcm', round(obj.aperture), round(obj.slit_size));
%            end
            
        end
        
        function val = getCount(obj) % add options later??
            
            % add varargin overrides to internal properties
            
            val = obj.filter.getCount(obj.magnitude, obj.aperture_area, obj.expT);
            
        end
        
        function val = get.OBSLONG(obj)
            
            val = obj.ephem.longitude;
            
        end
        
        function val = get.OBSLAT(obj)
            
            val = obj.ephem.latitude;
            
        end
        
        function val = get.OBSEL(obj)
            
            val = obj.ephem.elevation;
            
        end
        
        function val = is_galactic_center(obj)
            
            if isempty(obj.RA_DEG) || isempty(obj.DEC_DEG)
                val = 0;
            elseif abs(obj.RA_DEG-270)<25 || abs(obj.DEC_DEG+20)<25
                val = 1; 
            else
                val = 0;
            end
            
        end
        
        function val = get_juldates(obj, timestamps)
            
            if isempty(obj.END_STAMP) || isempty(obj.ENDTIME) 
                val = [];
            else
                val = juliandate(util.text.str2time(obj.ENDTIME) + seconds(timestamps - obj.END_STAMP));
            end
            
        end
        
        function val = get_datetimes(obj, timestamps)
            
            if isempty(obj.END_STAMP) || isempty(obj.ENDTIME) 
                val = [];
            else
                val = util.text.str2time(obj.ENDTIME) + seconds(timestamps - obj.END_STAMP);
            end
            
        end
        
    end
           
    methods % setters
        
        % filter/telescope parameters
        function set.FOCLEN(obj, val)
           
            obj.FOCLEN = val;
            if ~isempty(obj.FOCLEN) && ~isempty(obj.TEL_APER)
                ratio = obj.FOCLEN./obj.TEL_APER; 
                if ~isequal(ratio, obj.F_RATIO)
                    obj.F_RATIO = ratio;
                end
            elseif ~isempty(obj.FOCLEN) && ~isempty(obj.F_RATIO)
                aperture = obj.FOCLEN./obj.F_RATIO; 
                if ~isequal(aperture, obj.TEL_APER)
                    obj.TEL_APER = aperture;
                end
            end
            
            if ~isempty(obj.FOCLEN) && ~isempty(obj.PIXSIZE)
                scale = obj.PIXSIZE.*1e-4./obj.FOCLEN./2./pi.*360.*3600;
                if ~isequal(scale, obj.SCALE)
                    obj.SCALE = scale;
                end
            elseif ~isempty(obj.FOCLEN) && ~isempty(obj.SCALE)
                pix = obj.SCALE.*obj.FOCLEN.*2.*pi./360./3600.*1e-4;
                if ~isequal(pix, obj.PIXSIZE)
                    obj.PIXSIZE = pix;
                end
            end
            
        end
        
        function set.TEL_APER(obj, val)
            
            obj.TEL_APER = val;
            if ~isempty(val) && ~isempty(obj.FOCLEN)
                ratio = obj.FOCLEN./obj.TEL_APER; 
                if ~isequal(ratio, obj.F_RATIO)
                    obj.F_RATIO = ratio;
                end
            elseif ~isempty(val) && ~isempty(obj.F_RATIO)
                focal = obj.F_RATIO.*obj.TEL_APER;
                if ~isequal(focal, obj.FOCLEN)
                    obj.FOCLEN = focal;
                end
            end

        end
        
        function set.F_RATIO(obj, val)
            
            obj.F_RAT = val;
            if ~isempty(obj.F_RATIO) && ~isempty(obj.TEL_APER)
                focal = obj.TEL_APER.*obj.F_RATIO; % 
                if ~isequal(focal, obj.FOCLEN)
                    obj.FOCLEN = focal;
                end
            elseif ~isempty(obj.F_RATIO) && ~isempty(obj.FOCLEN)
                aperture = obj.FOCLEN./obj.F_RATIO;
                if ~isequal(aperture, obj.TEL_APER)
                    obj.TEL_APER = aperture;
                end
            end
            
        end
        
        function set.SCALE(obj, val)
            
            obj.SCALE = val;
            if ~isempty(obj.SCALE) && ~isempty(obj.FOCLEN) 
                pix = 2.*pi./360./3600.*1.e4.*obj.FOCLEN.*obj.SCALE;
                if ~isequal(pix, obj.PIXSIZE)
                    obj.PIXSIZE = pix;
                end
            elseif ~isempty(obj.SCALE) && ~isempty(obj.PIXSIZE)
                focal = obj.PIXSIZE./(2.*pi./360./3600.*1.e4.*obj.SCALE);
                if ~isequal(focal, obj.FOCLEN)
                    obj.FOCLEN = focal;
                end
            end
            
        end
        
        function set.PIXSIZE(obj, val)
            
            obj.PIXSIZE = val;
            if ~isempty(obj.PIXSIZE) && ~isempty(obj.FOCLEN)
                scale = obj.PIXSIZE./(2.*pi./360./3600.*1.e4.*obj.FOCLEN);
                if ~isequal(scale, obj.SCALE)
                    obj.SCALE = scale;
                end
            elseif ~isempty(obj.PIXSIZE) && ~isempty(obj.SCALE)
                focal = obj.PIXSIZE./(2.*pi./360./3600.*1.e4.*obj.SCALE);
                if ~isequal(focal, obj.FOCLEN)
                    obj.FOCLEN = focal;
                end
            end
            
        end
        
        function set.FILT_WAVE(obj, val)
            
            if ~isempty(obj.filter_obj)
                obj.filter_obj.wavelength = val; 
            end
            
        end
        
        function set.FILT_WIDTH(obj, val)
            
            if ~isempty(obj.filter_obj)
                obj.filter_obj.bandwidth= val; 
            end
            
        end

        % ephemeris
        function set.OBJECT(obj, val)
            
            obj.ephem.name = val;
            
        end
        
        function set.RA(obj, val)
            
            if iscell(val) && isscalar(val)
                obj.ephem.RA = val{1};
            else
                obj.ephem.RA = val;
            end
            
        end 
        
        function set.RA_DEG(obj, val)
            
            obj.ephem.RA = val/15;
            
        end
        
        function set.DEC_DEG(obj, val)
            
            obj.ephem.Dec = val;
            
        end
        
        function set.DEC(obj, val)
            
            if iscell(val) && isscalar(val)
                obj.ephem.Dec = val{1};
            else
                obj.ephem.Dec = val;
            end
            
        end
        
        function set.FIELD_ID(obj, val)
            
            obj.ephem.field_id = val; 
            
        end
        
        function set.STARTTIME(obj, val)
            
            obj.STARTTIME = val;
            obj.ephem.start_time = val;
            
        end
        
        function set.ENDTIME(obj, val)
            
            obj.ENDTIME = val;
            obj.ephem.end_time = val;
            
        end
        
        function set.OBSLONG(obj, val)
            
            obj.ephem.longitude = val;
            
        end
        
        function set.OBSLAT(obj, val)
            
            obj.ephem.latitude = val;
            
        end
        
        function addStar(obj, varargin)
                        
            obj.stars(end+1) = head.Star('pars', obj, varargin{:});
        
        end

        function set.stars(obj, val)
            
            if ~isempty(val)            
                util.oop.setprop(val, 'head', obj); % make sure all stars share the same "pars" object
                if ~isempty(obj.stars) && ~isempty(obj.stars(1).gui)
                    util.oop.setprop(val, 'gui', obj.stars(1).gui); % make sure all stars share the same "pars" object                
                end
            end
            
            obj.stars = val;
            
        end
        
        % target properties
        function name_plus_plus(obj)
            
            obj.TARGET = util.text.strPlusPlus(obj.TARGET);
            
        end
        
    end
    
    methods % save/load 
       
        function save(obj, filename, varargin) % not yet implemented! 
        % Usage: save(filename, varargin)
        % Saves the Header object to file named "filename". 
        % Simply calls the util.oop.save() function and passes it the 
        % object and filename and any additional arguments. 
        % The varargin is pre-appended with {'name', 'header'} which sets
        % the object name and HDF5 group name (location) to "header", but
        % this can be overwritten by specifying the 'name' property in the
        % optional arguments. 
        % 
        % See the help for util.oop.save for more details. 

            if nargin<2
                help('head.Header.save');
                help('util.oop.save'); 
                return;
            end
            
            % add the default name=header keyword-value pair (can be overwritten)
            varargin = [{'name', 'header'}, varargin];            
            
            util.oop.save(obj, filename, varargin{:}); 
            
        end
        
        function load(obj, filename, locations) % not yet implemented! 
        % Usage: load(filename, locations={header, head, pars})
        % Uses the util.oop.load() function to get the header information 
        % from a file named "filename". 
        % The function can understand from the file extension if it is a 
        % text file, mat-file, or HDF5 file. 
        % This function will try to get a header from a few different locations
        % in the file: header, head, pars and maybe a few more. 
        % User can overwrite this list by specifying a location string or 
        % multiple locations as a cell of strings. 
        % 
        % The output from the loading will be cast into a Header object if 
        % it is loaded as the older version Parameters object. 
        % 
        
            if nargin<2
                help('head.Header.load');
                return;
            end
            
            if nargin<3 || isempty(locations)
                locations = {'header', 'head', 'pars'}; 
            end
            
            if ischar(locations)
                locations = {locations};
            end
            
            loaded_header = [];
            
            for ii = 1:length(locations)
                
                try
                    
                    try
                        loaded_header = util.oop.load(filename, 'location', location{ii}); 
                    catch
                        if location{ii}(1)=='/'
                            loaded_header = util.oop.load(filename, 'location', location{ii}(2:end)); % try without the / at the beginning
                        else
                            loaded_header = util.oop.load(filename, 'location', ['/' location{ii}]); % try after adding a / at the beginning
                        end
                    end
                    
                end
                
            end
            
            if isempty(loaded_header)
                error('Could not find a header object in file: %s', filename);
            else
                
                if isa(loaded_header, 'head.Parameters')
                    loaded_header = cast(loaded_header); 
                end
                
                util.oop.copy_props(obj, loaded_header); 
                
            end
            
        end
                
        function writeFITS(obj, filename, timestamp, num_sum) % open a fits file and write keywords into it
            
            file_ptr = matlab.io.fits.openFile(filename,'readwrite');
            cleanup = onCleanup(@() matlab.io.fits.closeFile(file_ptr)); % make sure file closes at the end
            
            if nargin<3 || isempty(timestamp)
                timestamp = [];
            end
            
            if nargin<4 || isempty(num_sum)
                num_sum = 1;
            end
            
            obj.writeFitsHeader(file_ptr, timestamp, num_sum);
            
        end
        
        function writeFitsHeader(obj, file_ptr, timestamp, num_sum) % write keywords into an open FITS file 
            
            if nargin<3 || isempty(timestamp)
                start_time_str = obj.STARTTIME;
            else
                start_time_str = util.text.time2str(obj.observation_time(timestamp)); 
            end
            
            if nargin<4 || isempty(num_sum)
                num_sum = 1;
            end
            
            if ~isempty(start_time_str)
                obj.ephem.time = util.text.str2time(start_time_str);
                matlab.io.fits.writeKey(file_ptr, 'DATE-OBS', start_time_str); 
            end
            
            % this format was adopted from the work I did with microFUN, this is their standard header
            if ~isempty(obj.EXPTIME), matlab.io.fits.writeKey(file_ptr, 'EXPTIME',obj.EXPTIME.*num_sum, 'seconds'); end
            if ~isempty(obj.PIXSIZE), matlab.io.fits.writeKey(file_ptr, 'XPIXSZ',obj.PIXSIZE, 'microns'); end
            if ~isempty(obj.PIXSIZE), matlab.io.fits.writeKey(file_ptr, 'YPIXSZ',obj.PIXSIZE, 'microns'); end
            if ~isempty(obj.BINX), matlab.io.fits.writeKey(file_ptr, 'XBINNING',obj.BINX); end
            if ~isempty(obj.BINY), matlab.io.fits.writeKey(file_ptr, 'YBINNING',obj.BINY); end
            if ~isempty(obj.FILTER), matlab.io.fits.writeKey(file_ptr, 'FILTER',obj.FILTER); end            
            if ~isempty(obj.FILTER), matlab.io.fits.writeKey(file_ptr, 'FILT_WAVE',obj.FILT_WAVE); end
            if ~isempty(obj.FILTER), matlab.io.fits.writeKey(file_ptr, 'FILT_WIDTH',obj.FILT_WIDTH); end
            if ~isempty(obj.TYPE), matlab.io.fits.writeKey(file_ptr, 'IMAGETYP', obj.TYPE); end
            if ~isempty(obj.FOCLEN), matlab.io.fits.writeKey(file_ptr, 'FOCALLEN',obj.FOCLEN*10, 'mm'); end
            if ~isempty(obj.APERTURE), matlab.io.fits.writeKey(file_ptr, 'APTDIA', obj.APERTURE*10, 'mm'); end
            if ~isempty(obj.APERTURE), matlab.io.fits.writeKey(file_ptr, 'APTAREA',obj.aperture_area*100, 'mm^2'); end
            if ~isempty(obj.GAIN), matlab.io.fits.writeKey(file_ptr, 'EGAIN', obj.GAIN, 'e-/ADU'); end
            if ~isempty(obj.TARGET), matlab.io.fits.writeKey(file_ptr, 'OBJECT',obj.TARGET); end
            if ~isempty(obj.RA), matlab.io.fits.writeKey(file_ptr, 'OBJCTRA',obj.RA_DEG, 'degree'); end
            if ~isempty(obj.DE), matlab.io.fits.writeKey(file_ptr, 'OBJCTDEC',obj.DEC_DEG, 'degree'); end
            if ~isempty(obj.HA), matlab.io.fits.writeKey(file_ptr, 'OBJCTHA', obj.HA_DEG, 'degree'); end
            if ~isempty(obj.OBSLONG), matlab.io.fits.writeKey(file_ptr, 'SITELONG', obj.OBSLONG, 'degrees'); end
            if ~isempty(obj.OBSLAT), matlab.io.fits.writeKey(file_ptr, 'SITELAT', obj.OBSLAT, 'degrees'); end
            if ~isempty(obj.JD), matlab.io.fits.writeKey(file_ptr, 'JD', obj.JD); end
            if ~isempty(obj.AIRMASS), matlab.io.fits.writeKey(file_ptr, 'AIRMASS', obj.AIRMASS); end
            if ~isempty(obj.INST), matlab.io.fits.writeKey(file_ptr, 'INSTRUME', obj.INST); end
            matlab.io.fits.writeKey(file_ptr, 'INPUTFMT', 'FITS');
            
        end
        
        function readFitsHeader(obj, filename)
            
            in = fitsinfo(filename); 
            
            keys = in.PrimaryData.Keywords;
            
            obj.COMMENT = ''; 
            
            for ii = 1:size(keys,1)
                
                success = 1; 
                
                if strcmp(keys{ii,1}, 'COMMENT') && ~isempty(keys{ii,2})
                    obj.comment = [obj.comment ' ' keys{ii,2}]; % append to comment
                elseif isprop(obj, keys{ii,1})
                    try
                        obj.(keys{ii,1}) = keys{ii,2}; 
                    end
                elseif strcmp(keys{ii,1}, 'DATE-OBS')
                    obj.STARTTIME = keys{ii,2}; 
                    obj.ephem.time = keys{ii,2};
                elseif strcmp(keys{ii,1}, 'XPIXSZ') || strcmp(keys{ii,1}, 'YPIXSZ')
                    obj.PIXSIZE = keys{ii,2}; 
                elseif strcmp(keys{ii,1}, 'XBINNING')
                    obj.BINX = keys{ii,2}; 
                elseif strcmp(keys{ii,1}, 'YBINNING')
                    obj.BINY = keys{ii,2}; 
                elseif strcmp(keys{ii,1}, 'YBINNING')
                    obj.YBIN = keys{ii,2}; 
                elseif strcmp(keys{ii,1}, 'APTDIA')
                    if isempty(keys{ii,3}) || strcmpi(keys{ii,3}, 'cm')
                        obj.TEL_APER = keys{ii,2}; 
                    elseif strcmpi(keys{ii,3}, 'mm')
                        obj.TEL_APER = keys{ii,2}/10; 
                    elseif strcmpi(keys{ii,3}, 'm')
                        obj.TEL_APER = keys{ii,2}*100; 
                    end
                elseif strcmp(keys{ii,1}, 'FOCALLEN')
                    if isempty(keys{ii,3}) || strcmpi(keys{ii,3}, 'cm')
                        obj.FOCLEN = keys{ii,2}; 
                    elseif strcmpi(keys{ii,3}, 'mm')
                        obj.FOCLEN = keys{ii,2}/10; 
                    elseif strcmpi(keys{ii,3}, 'm')
                        obj.FOCLEN = keys{ii,2}*100; 
                    end
                elseif strcmp(keys{ii,1}, 'OBJCTRA')
                    obj.OBJRA = keys{ii,2}; 
                elseif strcmp(keys{ii,1}, 'OBJCTDEC')
                    obj.OBJDEC = keys{ii,2}; 
                elseif strcmp(keys{ii,1}, 'OBJCTHA')
                    obj.OBJRA = keys{ii,2}; 
                elseif strcmp(keys{ii,1}, 'INSTRUME')
                    obj.INST= keys{ii,2}; 
                elseif strcmp(keys{ii,1}, 'OBJCTHA')
                    obj.FOCLEN = keys{ii,2}; 
                elseif strcmp(keys{ii,1}, 'IMAGETYP')
                    obj.TYPE = keys{ii,2}; 
                else
                    success = 0; 
                end
                
                if obj.debug_bit>1
                    if success
                        fprintf('Reading keyword* "%s" with value "%s" of type "%s".\n', keys{ii,1}, util.text.print_value(keys{ii,2}), class(keys{ii,2}));
                    else
                        fprintf('No property matches the key "%s"...\n', keys{ii,1}); 
                    end
                end
                
            end
                        
        end
        
        function s = obj2struct(obj) % turn this object, and all sub objects, into structs
            
            warning('off', 'MATLAB:structOnObject');
            
            s = struct(obj); 
            
            list = properties(obj);
            
            for ii = 1:length(list)
               
                if isobject(obj.(list{ii}))
                    
                    % if any sub-object has embedded objects, must treat them individually
                    if isa(obj.(list{ii}), 'head.Star') && ~isempty(obj.(list{ii}))
                        s.(list{ii}) = []; % just until we figure this out
                         % need to figure out how to save a vector of
                         % structs and then add the primary_ref saved as
                         % numeric and not as object... 
                    elseif isa(obj.(list{ii}), 'head.gui.ParsGUI')
                        s.(list{ii}) = [];
                    else
                        s.(list{ii}) = struct(obj.(list{ii})); % turn the sub-objects into structs as well 
                    end
                    
                end
                
            end
            
        end
        
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
    
    methods % update and utilities
        
        function update(obj, varargin) % update the ephem object and any other internal calculations
           
            obj.ephem.update;
            obj.STARTTIME = util.text.time2str(obj.ephem.time); 
            
        end
        
        function val = observation_time(obj, timestamp)
            
            time = util.text.str2time(obj.ENDTIME); 
            
            val = time + seconds(timestamp - obj.END_STAMP); 
            
        end
        
        function val = stamp2str(obj, time_delay) % add some seconds to STARTTIME and return it as a string or cell array of strings
            
            time = util.text.str2time(obj.STARTTIME);
            
            time = time + seconds(time_delay);
            
            val = util.text.time2str(time);
            
        end
        
%         function val = run_identifier(obj) % get the identifier for the current run (date+object name)
%             
%             val = '';
%             
%             if length(obj.STARTTIME)>=10
%                 val = obj.STARTTIME(1:10);
%             end
%             
%             if ~isempty(obj.OBJECT)
%                 
%                 if ~isempty(val)
%                     val = [val '_'];
%                 end
%                 
%                 val = [val obj.OBJECT];
%                 
%             end
%             
%         end
        
    end
    
    methods (Static=true)
        
        function list = makeSyncList % list of parameters to give from dome-PC manager to the PcSync object (so it would be passed to the camera PC)
            
            list = {'OBJECT', 'OBJRA', 'OBJDEC', 'OBJRA_DEG', 'OBJDEC_DEG', 'TELRA', 'TELDEC', 'TELRA_DEG', 'TELDEC_DEG',...
                'FIELD_ID', 'TEMP_IN', 'TEMP_OUT', 'WIND_DIR', 'WIND_SPEED', 'HUMID_IN', 'HUMID_OUT', 'PRESSURE', 'LIGHT'};
                
        end
        
    end
    
    methods % plotting / GUI
       
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = head.gui.HeaderGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
    end
    
    % what follows is some boilerplate needed to keep this class backward
    % compatible after adhering to the ancient traditions of FITS headers. 
    % don't bother reading it, it just allows my old reasonable property
    % names to be used while displaying and saving only the FITS keywords. 
    
    properties(Transient=true, Hidden=true, Dependent=true)
        
        aperture;
        im_size;
        f_number;
        focal_length;
        frame_rate;
        frame_rate_measured;
        batch_size;
        target_name;
        instrument;
        t_start;
        t_end;
        run_start_datestr;
        notes;
        DE;
        focus;
        wavelength;
        bandwidth;
        pixel_size;
        plate_scale;
        longitude;
        latitude;
        SITELONG;
        SITELAT;
        
    end
    
    methods 
        
        function val = get.aperture(obj)
            
            val = obj.TEL_APER;
            
        end
        
        function set.aperture(obj, val)
            
            obj.TEL_APER = val;
            
        end
        
        function val = get.im_size(obj)
            
            val = [obj.NAXIS1 obj.NAXIS2];
            if obj.NAXIS>2
                val = [val obj.NAXIS3];
            end
            
            if obj.NAXIS>3
                val = [val obj.NAXIS4];
            end
            
        end
        
        function set.im_size(obj, val)
            
            obj.NAXIS = length(val);
            
            for ii = 1:obj.NAXIS
                obj.(['NAXIS' num2str(ii)]) = val(ii);
            end
            
        end
        
        function val = get.f_number(obj)
            
            val = obj.F_RATIO;
            
        end
        
        function set.f_number(obj, val)
            
            obj.F_RATIO = val;
            
        end
        
        function val = get.focal_length(obj)
            
            val = obj.FOCLEN;
            
        end
        
        function set.focal_length(obj, val)
            
            obj.FOCLEN = val;
            
        end
        
        function val = get.frame_rate(obj)
           
            val = obj.FRAMERATE;
            
        end
        
        function set.frame_rate(obj, val)
            
            obj.FRAMERATE = val;
            
        end
        
        function val = get.frame_rate_measured(obj)
           
            val = obj.ACT_FRAMERATE;
            
        end
        
        function set.frame_rate_measured(obj, val)
            
            obj.ACT_FRAMERATE = val;
            
        end
        
        function val = get.batch_size(obj)
            
            val = obj.NAXIS3;
            
        end
        
        function set.batch_size(obj, val)
            
            obj.NAXIS3 = val;
            
        end
        
        function val = get.target_name(obj)
           
            val = obj.OBJECT;
            
        end
        
        function set.target_name(obj, val)
            
            obj.OBJECT = val;
            
        end
        
        function val = get.instrument(obj)
            
            val = obj.INST;
            
        end
        
        function set.instrument(obj, val)
            
            obj.INST = val;
            
        end
        
        function val = get.notes(obj)
            
            val = obj.COMMENT;
            
        end
        
        function set.notes(obj, val)
            
            obj.COMMENT = val;
            
        end
        
        function val = get.DE(obj)
            
            val = obj.DEC;
            
        end
        
        function set.DE(obj, val)
            
            obj.DEC = val;
            
        end
        
        function val = get.focus(obj)
            
            val = obj.FOCUS_POS;
            
        end
        
        function set.focus(obj, val)
            
            obj.FOCUS_POS = val;
            
        end
        
        function val = get.t_start(obj)
            
            val = obj.STARTTIME;
            
        end
        
        function set.t_start(obj, val)
            
            obj.STARTTIME = val;
            
        end
        
        function val = get.t_end(obj)
            
            val = obj.ENDTIME;
            
        end
        
        function set.t_end(obj, val)
            
            obj.ENDTIME = val;
            
        end
        
        function val = get.run_start_datestr(obj)
            
            val = obj.RUNSTART;
            
        end
        
        function set.run_start_datestr(obj, val)
            
            obj.RUNSTART = val;
            
        end
        
        function val = get.wavelength(obj)
            
            val = obj.FILT_WAVE;
            
        end
        
        function set.wavelength(obj, val)
            
            obj.FILT_WAVE = val;
            
        end
        
        function val = get.bandwidth(obj)
            
            val = obj.FILT_WIDTH;
            
        end
        
        function set.bandwidth(obj, val)
            
            obj.FILT_WIDTH = val;
            
        end
        
        function val = get.plate_scale(obj)
            
            val = obj.SCALE;
            
        end
        
        function set.plate_scale(obj, val)
            
            obj.SCALE = val;
            
        end
        
        function val = get.pixel_size(obj)
            
            val = obj.PIXSIZE;
            
        end
        
        function set.pixel_size(obj, val)
            
            obj.PIXSIZE = val;
            
        end
        
        function val = get.longitude(obj)
            
            val = obj.OBSLONG;
            
        end
        
        function set.longitude(obj, val)
            
            obj.OBSLONG = val;
            
        end
        
        function val = get.latitude(obj)
            
            val = obj.OBSLAT;
            
        end
        
        function set.latitude(obj, val)
            
            obj.OBSLAT = val;
            
        end
        
        function val = get.SITELONG(obj)
            
            val = obj.OBSLONG;
            
        end
        
        function set.SITELONG(obj, val)
            
            obj.OBSLONG = val;
            
        end
        
        function val = get.SITELAT(obj)
            
            val = obj.OBSLAT;
            
        end
        
        function set.SITELAT(obj, val)
            
            obj.OBSLAT = val;
            
        end
        
    end
    
end
