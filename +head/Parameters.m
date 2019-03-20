classdef Parameters < dynamicprops
% this class contains all the useful information to be accessed by user
% or saved as metadata for a specific observation. 
% also lets you save or read metadata from file

    properties(Transient=true)
        
        gui@head.gui.ParsGUI;
        
    end

    properties % objects 
        
        filter; % must be a head.Filter (enforced in the setter...)
        ephem@head.Ephemeris;
        run_start_datetime@datetime;
        stars@head.Star;
        
    end
    
    properties 
        
        im_size;
        
        aperture = 57;
        f_number = 1.8947; 
                     
        % needs to be updated from camera/simulation or read from file
        frame_rate = 30;
        expT = 0.025;
        
%         clockFreq;
%         shutter;
        serial = 0;
        
        gain;
        batch_size;
        
        psf_sampling;
        
        focus;
        
        t_start; % most recent file (when it started filming)
        t_end; % most recent file (when it finished filming)
        
        target_name = 'star1';
        type = '';
        
        project = 'WFAST';
        instrument = 'Zyla_5.5';
        
        is_dark = 0;
        is_flat = 0;
        is_sim = 0;
        
        notes = '';
                
        seeing; % arcsec
        temperature; % celsius
        
    end
    
    properties(Dependent=true)
        
        % camera/telescope properties
        focal_length;        
        plate_scale; % arcseconds per pixel
        aperture_area;
        diff_limit; % arcseconds
        
        % filter properties
        filter_name;
        wavelength;
        bandwidth;
        
        % stars (specific for each cutout - row vectors or matrices)
        magnitude; % vector of magnitudes of all stars
        separation; % vector separations from 1st star (arcsec)
        pos_angle; % vector of position angles from 1st star (degrees)
        star_x; % in pixels, relative to frame size (final position)
        star_y; % in pixels, relative to frame size (final position)
        
        % ephemeris: time & coordinates        
        RA = ''; % for the center of the image
        DE = '';% for the center of the image
        
        HA;
        LST;
        ALT;
        AZ;        
        airmass;
        
        run_start_datestr;
        night_start_datestr;
        juldate;
        
    end
       
    properties(Hidden=true) 
       
        slit_size; % = 2.2; % in cm (also use aperture - the long edge of the asymetric scope)
        
        dead_time = 0;
        
        xbin = 1;
        ybin = 1;
                
        pixel_size = 6.5; % microns
        
        QE = 1;
        
        default_aperture;
        default_f_number;
        default_filter_name = 'clear';
        default_im_size = [2160 2560];
        
        datapath;
        dark_name = 'dark';
        flat_name = 'flat';
                
        use_folder_name = 1;
        use_folder_date = 1;
        
        ccd_id;
        amp_id;

%         psf_size;
%         psf_binning;
        
        AOI_height;
        AOI_top;
        AOI_width;
        AOI_left;
        
        debug_bit = 0;
        
        latitude = 31.907867; % of the observatory
        longitude = 34.811363; % of the observatory
        
        default_project;
        default_instrument;
        
        numImagesPerFile;
        
        version = 3.00;
          
    end
    
    methods % constructor
       
        function obj = Parameters(other)
                                 
            if nargin>0 && isa(other, 'head.Parameters')
                
                if obj.debug_bit, fprintf('Parameters copy-constructor v%4.2f\n', obj.version); end
                
                obj = util.oop.full_copy(other);
                
            else
            
                if obj.debug_bit, fprintf('Parameters constructor v%4.2f\n', obj.version); end
                
                obj.datapath = getenv('DATA');
                if isempty(obj.datapath)
                    obj.datapath = pwd;
                end
                
                obj.filter = head.Filter;
                obj.ephem = head.Ephemeris(obj);
                
                util.oop.save_defaults(obj);
                
                obj.update; % set the time to now. 
                
            end
            
        end
                
    end
    
    methods % reset methods
       
        function reset(obj)

            obj.resetTarget;
            obj.resetInstrument;
            obj.resetStars;
                        
        end
        
        function resetTarget(obj)
            
            obj.RA = '';
            obj.DE = '';
            
            if regexp(obj.target_name, 'star\d')==1
                obj.name_plus_plus;
            else
                obj.target_name = 'star1';
            end
            
            obj.type = '';
%             obj.spec_type = '';
            obj.notes = '';
            
        end
        
        function resetInstrument(obj)
            
            obj.aperture = obj.default_aperture;
            obj.f_number = obj.default_f_number;
            obj.filter_name = obj.default_filter_name;
            
        end        
        
        function resetStars(obj)
           
            obj.stars = head.Star.empty;
            
        end
        
        function clearObsConditions(obj)
           
            obj.seeing = [];
            obj.temperature = [];
            
        end
                
        function resetDrifts(obj)
           
            util.oop.setprop(obj.stars, 'drift_x', []);
            util.oop.setprop(obj.stars, 'drift_y', []);
            
        end
        
    end
    
    methods % getters
        
        % telescope parameters
        function val = get.focal_length(obj)
           
            val = obj.aperture.*obj.f_number;
            
        end
        
        function val = get.plate_scale(obj)
           
            val = obj.ephem.rad2arcsec(obj.pixel_size.*1e-4./obj.focal_length);
            
        end
        
        function val = get.diff_limit(obj)
            
            if isempty(obj.aperture) || isempty(obj.wavelength)
                val = [];
            else
                val = obj.wavelength*1e-7./obj.aperture*360*3600/2/pi;
            end
            
        end
        
        function val = diff_limit_pix(obj)
           
            if isempty(obj.diff_limit) || isempty(obj.diff_limit)
                val = [];
            else            
                val = obj.diff_limit./obj.plate_scale;
            end
            
        end
        
        function val = get.aperture_area(obj)
            
            if isempty(obj.slit_size)
                val = pi*(obj.aperture/2).^2;
            else
                val = obj.aperture*obj.slit_size;
            end
            
        end
        
        function val = image_sampling(obj)
            
            if isempty(obj.diff_limit) || isempty(obj.plate_scale)
                val = [];
            else
                val = obj.diff_limit/obj.plate_scale;
            end
            
        end
                
        % camera and filter 
        function val = get.im_size(obj)
            
            if isempty(obj.im_size)
                val = obj.default_im_size;
            else
                val = obj.im_size;
            end
            
            val = util.vec.imsize(val);
            
        end
        
        function str = cam_name(obj)
           
            c = strsplit(obj.instrument, {' ','_'});
            
            str = c{1};
            
        end
       
        function val = get.filter_name(obj)
            
            if isempty(obj.filter)
                val = '';
            else
                val = obj.filter.name;
            end
            
        end
        
        function val = get.wavelength(obj)
            
            if isempty(obj.filter)
                val = [];
            else
                val = obj.filter.wavelength;
            end
            
        end
        
        function val = get.bandwidth(obj)
            
            if isempty(obj.filter)
                val = [];
            else
                val = obj.filter.bandwidth;
            end
            
        end
         
        % ephemeris
        function name = folder_name(obj)
            
            import util.text.cs;
            
            if cs(obj.type, 'dark') && ~isempty(obj.dark_name)
                name = obj.dark_name;
            elseif cs(obj.type, 'flat') && ~isempty(obj.flat_name)
                name = obj.flat_name;
            else 
                name = obj.target_name;
            end

        end
        
        function f = folder(obj) % to be depricated!
            
            import util.text.cs;
            
            f = obj.date_folder;
            
            if obj.use_folder_name 
                
                name = '';
                if cs(obj.type, 'dark') && ~isempty(obj.dark_name)
                    name = obj.dark_name;
                elseif cs(obj.type, 'flat') && ~isempty(obj.flat_name)
                    name = obj.flat_name;
                elseif ~isempty(obj.target_name)
                    name = obj.target_name;
                end
                
                if ~isempty(name)
                    f = fullfile(f,  name);
                end
                
            end
            
        end
           
        function p = path(obj) % to be depricated!
            
            p = obj.folder;
            
        end
        
        function f = date_folder(obj) % to be depricated!
            
            f = obj.datapath;
            
            if obj.use_folder_date && ~builtin('isempty', obj.run_start_datetime)
            
                f = fullfile(f, datestr(obj.run_start_datetime, 'yyyy-mm-dd/'));
                                
            elseif obj.use_folder_date && ~isempty(obj.t_end)
                
                f = [fullfile(f, obj.t_end(1:10)) '/'];
                
            end
            
        end
        
        function val = get.run_start_datestr(obj)
            
            val = util.text.time2str(obj.run_start_datetime);
            
        end

        function val = get.night_start_datestr(obj)
            
            val = '';
            
        end
        
        function val = get.juldate(obj)
            
            val = obj.ephem.juldate;
            
        end
        
        function val = get.RA(obj)
            
            val = obj.ephem.RA;
            
        end
        
        function val = get.DE(obj)
            
            val = obj.ephem.DE;
            
        end
        
        function val = get.HA(obj)
            
            val = obj.ephem.HA;
            
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
        
        function val = get.airmass(obj)
            
            val = obj.ephem.airmass;
            
        end        
        
        % stars
        function val = get.magnitude(obj)
           
            if isempty(obj.stars)
                val = [];
                return;
            end
            
            val = [obj.stars.mag];
                        
        end
        
        function val = get.separation(obj)
                        
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
        
        function val = get.pos_angle(obj)
                        
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
                
        function val = get.star_x(obj)
                        
            if isempty(obj.stars)
                val = [];
                return;
            end
                  
            val = zeros(size(obj.stars));
            
            for ii = 1:length(val)
                val(ii) = obj.stars(ii).getFinalX('pixels');
            end
            
            
        end
                
        function val = get.star_y(obj)
                        
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
            
           if isempty(obj.slit_size)
               str = sprintf('%dcm', round(obj.aperture));
           else
               str = sprintf('%dx%dcm', round(obj.aperture), round(obj.slit_size));
           end
            
        end
        
        function val = getCount(obj) % add options later??
            
            % add varargin overrides to internal properties
            
            val = obj.filter.getCount(obj.magnitude, obj.aperture_area, obj.T);
            
        end
        
    end
           
    methods % setters

        % objects with loop-back
        function set.ephem(obj, val)
            
            obj.ephem = val;
            obj.ephem.pars = obj;
            
        end
        
        % filter/telescope parameters
        function set.focal_length(obj, val)
           
            obj.f_number = val./obj.aperture;
            
        end
        
        function set.plate_scale(obj, val)
           
            obj.pixel_size = util.units.arcsec2rad(val).*1.e4.*obj.focal_length;
            
        end
        
        function set.filter(obj, val)
            
            if ischar(val)
                obj.filter_name = val;
                return;
            elseif isa(val, 'head.Filter')
                obj.filter = val;
            else
                error('Must enter a head.Filter object, or a name of a known filter (e.g. R,V,B...)');
            end
            
        end
        
        function set.filter_name(obj, val)
            
            if isempty(obj.filter) || ~strcmpi(obj.filter_name, val)
                obj.filter = head.Filter(val);
            end
            
        end
        
        function set.run_start_datestr(obj, val)
            
            obj.run_start_datetime = astro.text.str2time(val);
            
        end

        % ephemeris
        function set.RA(obj, val)
            
            obj.ephem.RA = val;
            
        end 
        
        function set.DE(obj, val)
            
            obj.ephem.DE = val;
            
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
            
            obj.target_name = util.text.strPlusPlus(obj.target_name);
            
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
        
        input = util.text.InputVars;
        input.input_var('type', 'hdf5');
        input.input_var('location', '/');
        input.input_var('name', inputname(1), 'dataset name');
        input.input_var('if_exists', 'overwrite', 'exists');
        
        if nargin<2
            help('head.Parameters.save');
            return;
        end
        
        
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
                
        function writeFITS(obj, filename, timestamp)
            
            file_ptr = matlab.io.fits.openFile(filename,'readwrite');
            cleanup = onCleanup(@() matlab.io.fits.closeFile(file_ptr));
            
            if nargin<3 || isempty(timestamp)
%                 timestamp = obj.obs_time_days;
%                 if isempty(timestamp)
%                     timestamp = datenum(clock);
%                 end
                timestamp = [];
            end
            
            try
                
                obj.writeFitsHeader(file_ptr, timestamp);
                
            catch ME
                rethrow(ME);
            end
            
        end
        
        function writeFitsHeader(obj, file_ptr, time_delay)
            
            if nargin<2 || isempty(time_delay)
                start_time_str = obj.stamp2str(time_delay);
            else
                start_time_str = obj.t_start;
            end
            
            if ~isempty(start_time_str)
                obj.ephem.time = util.text.str2time(start_time_str);
                matlab.io.fits.writeKey(file_ptr, 'DATE-OBS', start_time_str); 
            end
            
            if ~isempty(obj.expT), matlab.io.fits.writeKey(file_ptr, 'EXPTIME',obj.expT, 'seconds'); end
            if ~isempty(obj.pixel_size), matlab.io.fits.writeKey(file_ptr, 'XPIXSZ',obj.pixel_size, 'microns'); end
            if ~isempty(obj.pixel_size), matlab.io.fits.writeKey(file_ptr, 'YPIXSZ',obj.pixel_size, 'microns'); end
            if ~isempty(obj.xbin), matlab.io.fits.writeKey(file_ptr, 'XBINNING',obj.xbin); end
            if ~isempty(obj.ybin), matlab.io.fits.writeKey(file_ptr, 'YBINNING',obj.ybin); end
            if ~isempty(obj.filter_name), matlab.io.fits.writeKey(file_ptr, 'FILTER',obj.filter_name); end
            if ~isempty(obj.type), matlab.io.fits.writeKey(file_ptr, 'IMAGETYP', obj.type); end
            if ~isempty(obj.focal_length), matlab.io.fits.writeKey(file_ptr, 'FOCALLEN',obj.focal_length*10, 'mm'); end
            if ~isempty(obj.aperture), matlab.io.fits.writeKey(file_ptr, 'APTDIA', obj.aperture*10, 'mm'); end
            if ~isempty(obj.aperture), matlab.io.fits.writeKey(file_ptr, 'APTAREA',(obj.aperture*10/2).^2*pi, 'mm^2'); end
            if ~isempty(obj.gain), matlab.io.fits.writeKey(file_ptr, 'EGAIN', obj.gain, 'e-/ADU'); end
            if ~isempty(obj.target_name), matlab.io.fits.writeKey(file_ptr, 'OBJECT',obj.target_name); end
            if ~isempty(obj.RA), matlab.io.fits.writeKey(file_ptr, 'OBJCTRA',obj.RA, 'hours'); end
            if ~isempty(obj.DE), matlab.io.fits.writeKey(file_ptr, 'OBJCTDEC',obj.DE, 'degree'); end
            if ~isempty(obj.HA), matlab.io.fits.writeKey(file_ptr, 'OBJCTHA', head.Ephemeris.ra2rad(obj.HA)/2/pi*24, 'hours'); end
            if ~isempty(obj.latitude), matlab.io.fits.writeKey(file_ptr, 'SITELAT', obj.latitude, 'degrees'); end
            if ~isempty(obj.longitude), matlab.io.fits.writeKey(file_ptr, 'SITELONG', obj.longitude, 'degrees'); end
            if ~isempty(obj.juldate), matlab.io.fits.writeKey(file_ptr, 'JD', obj.juldate); end
            if ~isempty(obj.airmass), matlab.io.fits.writeKey(file_ptr, 'AIRMASS', obj.airmass); end
            if ~isempty(obj.instrument), matlab.io.fits.writeKey(file_ptr, 'INSTRUME', obj.instrument); end
            matlab.io.fits.writeKey(file_ptr, 'INPUTFMT', 'FITS');
            
        end
        
    end
    
    methods % update and utilities
        
        function update(obj, varargin)
           
            obj.run_start_datetime = datetime('now', 'timezone', 'utc');
            
        end
        
        function updateObsTimeOffset(obj, time_offset_seconds) % to be depricated?
            
            obj.obs_time = obj.obs_time + time_offset_seconds/24/3600;
                        
        end
        
        function val = stamp2str(obj, time_delay)
            
            time = util.text.str2time(obj.t_start);
            
            time.Second = time.Second + time_delay;
            
            val = util.text.time2str(time);
            
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
    
end