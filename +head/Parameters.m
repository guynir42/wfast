classdef (CaseInsensitiveProperties, TruncatedProperties) Parameters < dynamicprops
% This class contains all the useful information to be accessed by user
% or saved as metadata for a specific observation. 
% also lets you save or read metadata from file

    properties(Transient=true)
        
        gui@head.gui.ParsGUI;
        
        cat@head.Catalog; % this is transient right now because I don't know how to save tables (yet)
        
    end

    properties % objects 
        
        filter_obj; % must be a head.Filter (enforced in the setter...)
        ephem@head.Ephemeris;
%         run_start_datetime;
        stars@head.Star;
        
        WCS@head.WorldCoordinates; % to be expanded later
        
    end
    
    properties 
        
%         target_name = 'star1';
        OBJECT = 'star1';
%         type = '';
        TYPE;
        
        COMMENT = '';
        
        PROJECT = 'WFAST';
        INST = 'Zyla_5.5';
        
%         t_start; % most recent file (when it started filming)
%         t_end; % most recent file (when it finished filming)
        STARTTIME;
        ENDTIME;
        RUNSTART; 
        
%         aperture = 57;
        TEL_APER = 57;
%         f_number = 1.8947; 
        FOCLEN = 108;
        F_RATIO;
        
        PIXSIZE = 6.5; % microns
        SCALE; 
        
%         expT;
        EXPTIME;
%         frame_rate;
        FRAMERATE;
%         frame_rate_measured;
        ACT_FRAMERATE; 
        
        FILTER = 'F505W'; 
        
%         batch_size;
        
%         im_size;
        NAXIS; 
        NAXIS1;
        NAXIS2;
        NAXIS3;
        NAXIS4;
        
        BINX = 1;
        BINY = 1;
        ROI; 
        CCDSEC;
        
%         is_dark = 0;
%         is_flat = 0;
%         is_sim = 0;
        
%         gain;
        GAIN; 
        READNOISE;
        DARKCUR;
        
        IS_DARK = 0;
        IS_FLAT = 0;
        IS_SIMULATED = 0;

%         notes = '';
               
        FOCUS_POS;
        FOCUS_TIP;
        FOCUS_TILT;
         
        SEEING; % arcsec
%         temperature; % celsius
        TEMP_DET; % detector temperature (celsius)
        TEMP_OUT; % outside temperature (celsius)
        
        WIND_DIR;
        WIND_SPEED;
        HUMID_IN;
        HUMID_OUT;
        PRESSURE;
        LIGHT;

        % hardware (mount) coordinates
        TELRA;
        TELDEC; 
        TELRA_DEG;
        TELDEC_DEG;
        
    end
    
    properties(Dependent=true)
        
        % camera/telescope properties
%         focal_length;        
%         plate_scale; % arcseconds per pixel
        
%         diff_limit; % arcseconds
        
        % stars (specific for each cutout - row vectors or matrices)
%         magnitude; % vector of magnitudes of all stars
%         separation; % vector separations from 1st star (arcsec)
%         pos_angle; % vector of position angles from 1st star (degrees)
%         star_x; % in pixels, relative to frame size (final position)
%         star_y; % in pixels, relative to frame size (final position)
        
        % ephemeris: time & coordinates
        
        % measured from observations
        OBSRA; 
        OBSDEC;
        OBSRA_DEG;
        OBSDEC_DEG;
        
        % given by scheduler
        RA; % for the center of the image
        DEC;% for the center of the image
        RA_DEG;
        DEC_DEG;
        
        
        HA;
        HA_DEG;
        LST;
        ALT;
        AZ;        
        AIRMASS;
        
        MOONAZ;
        MOONALT;
        MOONILL;
        MOONDIST;

        SUNAZ;
        SUNALT;
        
%         longitude;
%         latitude;
        OBSLONG;
        OBSLAT;
         
%         run_start_datestr;
        JD;
        MJD;
%         MIDJD;
        
        % filter properties
%         wavelength;
%         bandwidth;
        FILT_WAVE;
        FILT_WIDTH;

    end
    
    properties(Hidden=true) 
       
%         SLIT_SIZE; % = 2.2; % in cm (also use aperture - the long edge of the asymetric scope)
        
        DEADTIME = 0;
        
%         xbin = 1;
%         ybin = 1;
        
%         CCDID
%         AMPID 
           
%         pixel_size = 6.5; % microns
        
        QE = 1;
        
        default_APERTURE;
        default_F_RATIO;
        default_FILTER;
        default_OBJECT;
        
        filter_name_full; % for backward compatibility with older versions.
        
%         datapath;
%         dark_name = 'dark';
%         flat_name = 'flat';
%                 
%         use_folder_name = 1;
%         use_folder_date = 1;
%         
%         AOI_height;
%         AOI_top;
%         AOI_width;
%         AOI_left;
        
%         latitude = 31.907867; % of the observatory
%         longitude = 34.811363; % of the observatory
        
        CAMS_VER;
        TELS_VER;
        
        debug_bit = 0;
        version = 4.00;
          
    end
    
    methods % constructor
       
        function obj = Parameters(other)
                                 
            if nargin>0 && isa(other, 'head.Parameters')
                
                if obj.debug_bit, fprintf('Parameters copy-constructor v%4.2f\n', obj.version); end
                
                obj = util.oop.full_copy(other);
                
            else
            
                if obj.debug_bit, fprintf('Parameters constructor v%4.2f\n', obj.version); end
                
%                 obj.datapath = getenv('DATA');
%                 if isempty(obj.datapath)
%                     obj.datapath = pwd;
%                 end
                
                obj.filter_obj = head.Filter(obj.FILTER);
                obj.ephem = head.Ephemeris;
                
                obj.FOCLEN = obj.FOCLEN; % this should fill the values for SCALE and F_RATIO
                
                util.oop.save_defaults(obj);
                
                obj.update; % set the time to now. 
                
            end
            
        end
                
    end
    
    methods % reset methods
        
        function clearStars(obj)
           
            obj.stars = head.Star.empty;
            
        end
        
        function clearObsConditions(obj)
           
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
                
        function resetTarget(obj)
            
            obj.OBJECT = obj.default_OBJECT;
            obj.ephem.reset;
            
        end
        
        function resetDrifts(obj)
           
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
           
            c = strsplit(obj.INSTR, {' ','_'});
            
            str = c{1};
            
        end
        
        function val = get.FILTER(obj)
        
            if ~isempty(obj.FILTER) && util.text.cs('[1x1 head.Filter]', obj.FILTER) % patch for backward compatibility with older headers where "filter" was the object
                obj.FILTER = obj.filter_name_full;
            end
                
            val = obj.FILTER;

        end
%         
%         function val = get.filter_name_full(obj)
%             
%             if isempty(obj.filter_name_full) && ~isempty(obj.FILTER)
%                 val = obj.FILTER;
%             else
%                 val = obj.filter_name_full;
%             end
%             
%         end
%         
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
            
            if isempty(obj.WCS) || isempty(obj.WCS.CRVAL)
                val = [];
            else
                val = obj.WCS.CRVAL(1);
            end
            
        end
        
        function val = get.OBSDEC_DEG(obj)
            
            if isempty(obj.WCS) || isempty(obj.WCS.CRVAL)
                val = [];
            else
                val = obj.WCS.CRVAL(2);
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
            
            val = obj.ephem.ALT;
            
        end
        
        function val = get.AZ(obj)
            
            val = obj.ephem.AZ;
            
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
%         
%         function set.FILTER(obj, val)
%             
%             if ~util.text.cs('[1x1 head.Filter]', val) % patch for backward compatibility with older headers where "filter" was the object
%             
%             end
%                 
%             if ~strcmp(obj.FILTER, val)
%                 
%                 obj.FILTER = val;
%                 
%                 try
%                     obj.filter_obj = head.Filter(val);
%                 catch ME
%                     disp(['Unknown filter name "' val '", cannot generate filter object']);
%                 end
%                 
%             end
%             
%         end
%         
        % ephemeris
        function set.RA(obj, val)
            
            if iscell(val) && isscalar(val)
                obj.ephem.RA = val{1};
            else
                obj.ephem.RA = val;
            end
            
        end 
        
        function set.DEC(obj, val)
            
            if iscell(val) && isscalar(val)
                obj.ephem.Dec = val{1};
            else
                obj.ephem.Dec = val;
            end
            
        end
        
        function addStar(obj, varargin)
                        
            obj.stars(end+1) = head.Star('pars', obj, varargin{:});
        
        end

        function set.stars(obj, val)
            
            if ~isempty(val)            
                util.oop.setprop(val, 'pars', obj); % make sure all stars share the same "pars" object
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
       
        function save(obj, filename, varargin)
        % saves the Parameters object to file named "filename". 
        % usage: pars.save(filename, varargin)
        % OPTIONAL PARAMETERS
        %   -type: HDF5 (default), fits, text, mat.
        %   -location: (HDF5 only) where to save inside the file (default is '/').
        %   -name: (HDF5 only) what to call the dataset (default is name of object, e.g. "pars").
        %   -if_exist: overwrite (default), error, warning, nothing

            if nargin<2
                help('head.Parameters.save');
                return;
            end

            input = util.text.InputVars;
            input.input_var('type', 'hdf5');
            input.input_var('location', '/');
            input.input_var('name', inputname(1), 'dataset name');
            input.input_var('if_exists', 'overwrite', 'exists');
            input.scan_vars(varargin);

            warning('head.Parameters.save is not yet implemented!');

        end
        
        function load(obj, filename)
        % usage: pars.load(filename). 
        
            if nargin<2
                help('head.Parameters.load');
                return;
            end
            
            warning('head.Parameters.load is not yet implemented!');
            
        end
                
        function writeFITS(obj, filename, time_delay, num_sum) % write keywords into a FITS file
            
            file_ptr = matlab.io.fits.openFile(filename,'readwrite');
            cleanup = onCleanup(@() matlab.io.fits.closeFile(file_ptr));
            
            if nargin<3 || isempty(time_delay)
%                 timestamp = obj.obs_time_days;
%                 if isempty(timestamp)
%                     timestamp = datenum(clock);
%                 end
                time_delay = [];
            end
            
            if nargin<4 || isempty(num_sum)
                num_sum = 1;
            end
            
            try
                
                obj.writeFitsHeader(file_ptr, time_delay, num_sum);
                
            catch ME
                rethrow(ME);
            end
            
        end
        
        function writeFitsHeader(obj, file_ptr, time_delay, num_sum)
            
            if nargin<3 || isempty(time_delay)
                start_time_str = obj.t_start;
            else
                start_time_str = obj.stamp2str(time_delay);
            end
            
            if nargin<4 || isempty(num_sum)
                num_sum = 1;
            end
            
            if ~isempty(start_time_str)
                obj.ephem.time = util.text.str2time(start_time_str);
                matlab.io.fits.writeKey(file_ptr, 'DATE-OBS', start_time_str); 
            end
            
            if ~isempty(obj.EXPTIME), matlab.io.fits.writeKey(file_ptr, 'EXPTIME',obj.EXPTIME.*num_sum, 'seconds'); end
            if ~isempty(obj.PIXSIZE), matlab.io.fits.writeKey(file_ptr, 'XPIXSZ',obj.PIXSIZE, 'microns'); end
            if ~isempty(obj.PIXSIZE), matlab.io.fits.writeKey(file_ptr, 'YPIXSZ',obj.PIXSIZE, 'microns'); end
            if ~isempty(obj.BINX), matlab.io.fits.writeKey(file_ptr, 'XBINNING',obj.BINX); end
            if ~isempty(obj.BINY), matlab.io.fits.writeKey(file_ptr, 'YBINNING',obj.BINY); end
            if ~isempty(obj.FILTER), matlab.io.fits.writeKey(file_ptr, 'FILTER',obj.FILTER); end
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
            if ~isempty(obj.INSTR), matlab.io.fits.writeKey(file_ptr, 'INSTRUME', obj.INSTR); end
            matlab.io.fits.writeKey(file_ptr, 'INPUTFMT', 'FITS');
            
        end
        
    end
    
    methods % update and utilities
        
        function update(obj, varargin)
           
            obj.ephem.update;
            
        end
        
        function val = stamp2str(obj, time_delay)
            
            time = util.text.str2time(obj.t_start);
            
            time = time + seconds(time_delay);
            
            val = util.text.time2str(time);
            
        end
        
    end
    
    methods (Static=true)
        
        function list = makeSyncList
            
            list = {'OBJECT', 'RA', 'DEC', 'RA_DEG', 'DEC_DEG', 'TELRA', 'TELDEC', 'TELRA_DEG', 'TELDEC_DEG',...
                'TEMP_IN', 'TEMP_OUT', 'WIND_DIR', 'WIND_SPEED', 'HUMID_IN', 'HUMID_OUT', 'PRESSURE', 'LIGHT'};
                
        end
        
    end
    
    methods % plotting / GUI
       
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = head.gui.ParsGUI(obj);
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
