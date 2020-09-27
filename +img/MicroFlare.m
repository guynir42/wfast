classdef MicroFlare < handle

    properties(Transient=true)
        
        cat@head.Catalog;
        cal@img.Calibration;
        
        all_cutouts; % a transient array holding the calibrated cutouts 
        all_positions; % might as well keep the positions as well
        
    end
    
    properties % inputs/outputs
        
        filename; % where the flare was detected
        folder; % what folder the flare was in
        serial; % index of the total number of flares detected in this run
        file_index; % index in the run
        frame_index; % what frame did the peak appear in
        cut_index; % which cutout (star) from the file cutouts is this flare coming from
        pos; % x and y position of the peak pixel in the full-frame image
        peak; % value of the peak pixel
        
        peak_mag; % magnitude of peak matched to other stars
        peak_time; % time of peak (datetime object)
        peak_region; % frames around peak that are also counted
        
        pixel_var; % variance of underlying pixel (from calibration file)
        cutouts; % cutout images around the peak and all frames before/after it in that batch
        
        % from photometry2 or from rough estimates on the whole cutout
        timestamps;
        flux;
        error;
        background;
        area;
        offset;        
        width;
        bad_pixels;
        
        mean; % mean flux calculated outside the peak region
        std; % RMS flux calculated outside the peak region
        
        back_mean; % mean of images
        back_std; % RMS of images
        
        magnitudes; % by comparing the other stars (vs. GAIA_BP)
        mag_error; % translate the relative error to mag error
        
        image; % single cutout at the frame of the peak
        num_peaks; % how many distinct peaks can we find in the one image
        num_frames; % how many continuous frames we have around the peak
        num_pixels; % how many pixels above the threshold do we have 
        num_blinks; % how many distinct times it flashes in one batch
        
        fwhm; % calculted from the image, using dedicated FWHM function (e.g., using filters)
        
        % calcExtras results:
        seeing; % calculate the seeing based on all the stars in the cutouts 
        offset_shift; % maximum motion of the offsets during the peak region
        velocity; % the speed inside a group/cluster (in pixels/sec)
        is_streaked; % if the flare can be see to streak across the cutout (LEO satellites)
        is_bad_pos; % if cutout is in a bad place like edge of sensor
        
        type = ''; % can be cosmic ray, flare, satellite, pixel... 
        
    end
    
    properties % switches/controls
        
        pixel_thresh = 256; % threshold of the brightest pixel, used for triggering 
        region_thresh = 5; % each pixel is counted if it surpasses this many multiples of  image noise 
        frac_thresh = 0.1; % threshold of pixels around the brightest pixel, as fraction of peak, for finding number of pixels
        flux_thresh = 7.5; % threshold for peak intensity (of the flux, not of single pixel)
        side_thresh = 3; % threshold for flux in frames adjacent to the peak
        interval_peak = 2; % how many frames in either direction from peak should be excluded for calculating the mean/std of the flux
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        corrected_flux;
        corrected_mean; 
        corrected_std; 
        
        first_peak_time; % in case of multiple blinks
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = MicroFlare(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'img.MicroFlare')
                if obj.debug_bit>1, fprintf('MicroFlare copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('MicroFlare constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
    end
    
    methods % getters
        
        function str = printSummary(obj)
            
            c = {}; 
            
            if ~isempty(obj.seeing)
                c{end+1} = sprintf('seeing= %4.1f"', obj.seeing); 
            end
            
            if ~isempty(obj.offset_shift)
                c{end+1} = sprintf('shift= %4.2f', obj.offset_shift); 
            end
            
            if ~isempty(obj.is_streaked)
                c{end+1} = sprintf('streaked= %d', obj.is_streaked); 
            end
            
            if ~isempty(obj.is_bad_pos)
                c{end+1} = sprintf('bad pos= %d', obj.is_bad_pos); 
            end
            
            if ~isempty(obj.num_pixels)
                c{end+1} = sprintf('pixels= %d', obj.num_pixels); 
            end
            
            if ~isempty(obj.peak_region)
                c{end+1} = sprintf('frames= %d', length(obj.peak_region)); 
            elseif ~isempty(obj.num_frames)
                c{end+1} = sprintf('frames= %d', obj.num_frames); 
            end
            
            if ~isempty(obj.num_blinks)
                c{end+1} = sprintf('blinks= %d', obj.num_blinks); 
            end
            
            if ~isempty(obj.velocity)
                c{end+1} = sprintf('velocity= %4.1f', obj.velocity); 
            end
            
            str = strjoin(c, ' | '); 
            
        end
        
    end
    
    methods % setters
        
        function set.filename(obj, val)
            
            [a,b,c] = fileparts(val);
            
            obj.filename = [b,c];
            
            if ~isempty(a) && isempty(obj.folder)
                obj.folder = a;
            end
            
        end
        
    end
    
    methods(Static=true) % utilities
        
        function obj = readStruct( st)
            
            if ~isstruct(st)
                error('Must supply a struct array to this function! Instead got a "%s" type object.', class(st)); 
            end
            
            list = properties(img.MicroFlare); 
            
            for ii = 1:length(st)
                obj(ii) = img.MicroFlare; 
                for jj = 1:length(list)
                    if isfield(st(ii), list{jj})
                        obj(ii).(list{jj}) = st(ii).(list{jj}); 
                    end
                end
            end
            
        end
        
    end
    
    methods % calculations
        
        function calculate(obj, varargin)
            
            if numel(obj)>1
                
                for ii = 1:numel(obj)
                    obj(ii).calculate(varargin{:}); % recursively call this for individual objects in the array
                end
                
                return;
                
            end
            
            input = util.text.InputVars;
            input.input_var('pixel_thresh', [], 'pixel_threshold'); 
            input.input_var('area_thresh', [], 'area_threshold'); 
            input.input_var('flux_thresh', [], 'flux_threshold'); 
            input.input_var('side_thresh', [], 'side_threshold'); 
            input.input_var('interval_peak', []); 
            input.input_var('debug_bit', []); 
            input.scan_vars(varargin{:}); 
            
            for ii = 1:length(input.list_added_properties)
                if ~isempty(input.(input.list_added_properties{ii})) && isprop(obj, input.list_added_properties{ii})
                    obj.(input.list_added_properties{ii}) = input.(input.list_added_properties{ii}); 
                end
            end
            
            obj.image = obj.cutouts(:,:,obj.frame_index); 
            
            flux_norm = obj.getNormFlux; 

            % how many frames above the threshold
            counter = 1;
            
            for ii = 1:length(flux_norm)
                
                idx = obj.frame_index - ii;
                
                if idx<1, break; end
                if flux_norm(idx)<obj.side_thresh, break; end
                
                counter = counter + 1;
                
            end
            
            for ii = 1:length(flux_norm)
               
                idx = obj.frame_index + ii;
                
                if idx>length(flux_norm), break; end
                if flux_norm(idx)<obj.side_thresh, break; end
                
                counter = counter + 1;
                
            end
            
            obj.num_frames = counter;
            
            obj.back_mean = nanmedian(util.stat.mean2(obj.cutouts)); 
            obj.back_std = nanmedian(util.stat.std2(obj.cutouts)); 
            
            % how many pixels in the image are above the threshold
            BW = (obj.image-obj.back_mean)./obj.back_std>obj.region_thresh; % black/white image of pixels above the threshold
            N = nnz(BW); 
            
            if N<=1
                obj.num_pixels = 1;
                obj.num_peaks = 1;
            else
                
                C = bwconncomp(BW); 
                [X,Y] = meshgrid((1:size(obj.image,2))-floor(size(obj.image,2)/2)-1, (1:size(obj.image,1))-floor(size(obj.image,1)/2)-1);
                
                obj.num_peaks = C.NumObjects;
                
                for ii = 1:C.NumObjects
                    
                    r(ii) = min(sqrt(X(C.PixelIdxList{ii}).^2 + Y(C.PixelIdxList{ii}).^2)); 
                    
                end
                
                [~,idx] = min(r); % find the connected component that is nearest to the center of the cutout
                
                N = numel(C.PixelIdxList{idx}); 
                
            end
            
            obj.num_pixels = N;
            
            obj.fwhm = util.img.fwhm(obj.image, 'method', 'filters'); 
            
            if obj.num_pixels==1
                obj.type = 'bad pixel';
            elseif obj.num_frames==1
                obj.type = 'cosmic ray';
            else
                obj.type = 'flare';
            end
            
        end
        
        function flux_norm = getNormFlux(obj)
            
            indices = obj.frame_index;
            for ii = 1:obj.interval_peak
                
                idx = obj.frame_index - ii;
                if idx<1, break; end
                
                indices = [indices idx]; 
                    
            end
            
            for ii = 1:obj.interval_peak
                
                idx = obj.frame_index + ii;
                if idx>length(obj.flux)
                    break;
                end
                
                indices = [indices idx]; 
                    
            end
            
            indices = sort(indices);
            flux_bg = obj.flux;
            flux_bg(indices) = NaN;
            
            obj.mean = nanmean(flux_bg);
            obj.std = nanstd(flux_bg); 
            
            flux_norm = (obj.flux-obj.mean)./obj.std;
            
        end
        
        function calcExtras(obj, cal, cat)
            
            if length(obj)>1 % for vector inputs
                
                for ii = 1:length(obj)
                    obj(ii).calcExtras(obj(1).cal, obj(1).cat); % make sure to use the same calibration and catalog for all flares in this run
                end
                
                return;
                
            end
            
            %%%%%%% get the calibration and catalog objects %%%%%%%%
                        
            if nargin<2 || isempty(cal)
                if isempty(obj.cal)
                    d = util.text.run_id(obj.folder); 
                    obj.cal = img.Calibration;
                    obj.cal.loadByDate(d(1:10), 'Balor'); 
                end
            else
                obj.cal = cal;
            end
            
            d = util.sys.WorkingDirectory;
            
            d.find_run(obj.folder); 
            
            if nargin<3 || isempty(cat)
                if isempty(obj.cat)
                    obj.cat = head.Catalog;
                    obj.cat.loadMAT(fullfile(d.pwd, 'catalog.mat')); 
                end
            else
                obj.cat = cat; 
            end
            
            if isempty(obj.cat) || obj.cat.success==0
                disp('Cannot load catalog.'); 
                return; 
            end
            
            %%%%%%%%%%% get the new path of the file %%%%%%%%%%%%%%%
            
            filename = fullfile(d.pwd, obj.filename); 
            if ~exist(filename, 'file')
                filename = fullfile(d.pwd, [obj.filename, 'z']); % try with the defalted extension
            end
            
            if ~exist(filename, 'file')
                fprintf('Cannot find file: %s \n', filename); 
                return;
            end
            
            try % get the other stars to calculate magnitudes

                C = h5read(filename, '/cutouts'); 
                P = h5read(filename, '/positions'); 

                obj.all_cutouts = obj.cal.input(C, 'pos', P); 
                obj.all_positions = P; 
                
                s = util.img.photometry2(obj.all_cutouts, 'aperture', [3 8], 'use_aperture', 1, 'use_gaussian', 1); % the big aperture is used to zero-in on the source even on the edges of the cutout
                
                f = s.apertures_photometry.flux - s.apertures_photometry.area .* s.apertures_photometry.background;
                x = s.apertures_photometry.offset_x(:,obj.cut_index); 
                y = s.apertures_photometry.offset_y(:,obj.cut_index); 
                
                ff = f(:,obj.cut_index); % flux for this specific micro flare

                fr = util.fit.polyfit((1:length(ff))', ff, 'order', 1, 'sigma', 2); 

                f2 = ff-fr.ym; % subtract the linear model
                s2 = nanstd(f2(~fr.bad_idx)); % noise, excluding outliers (peaks)

                obj.corrected_flux = f2;
                obj.corrected_mean = nanmean(fr.ym(~fr.bad_idx)); 
                obj.corrected_std = s2; 
                
                m = obj.cat.data.Mag_BP'; % set the magnitudes into a row vector

                S = nanmean(f(:,1:length(m))); % signal -> mean flux for each star
                N = nanstd(f(:,1:length(m))); % noise -> standard deviation

    %             dm = m + 2.5*log10(F); % difference between GAIA mag and the instrumental flux (converted to mag)            
                dm = m - util.units.flux2lup(S, 5); % difference between GAIA mag and the instrumental flux (converted to mag)            

                zp = util.vec.weighted_average(dm, sqrt(abs(S./N)),2);
                zp_err = util.stat.rstd(dm, 2); 
%                 zp_err = sqrt(util.vec.weighted_average((dm-nanmean(dm)).^2, sqrt(abs(S./N)),2))./sqrt(nansum((S./N))); % zero point error: https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Standard_error
                
                obj.magnitudes = zp + util.units.flux2lup(f2, 5);

                obj.peak_mag = obj.magnitudes(obj.frame_index); 
                obj.mag_error = sqrt(zp_err.^2 + (obj.corrected_std./max(obj.corrected_flux)).^2);
                

            catch ME
                rethrow(ME);
%                 warning(ME.getReport); 
            end
            
            try % get a better estimate for the seeing
                
                N = min(50, size(obj.all_cutouts,4)); % 50 best stars (or all stars if there are not enough)
                
                CC = util.img.centering(obj.all_cutouts(:,:,:,1:N)); 
                
                w = util.img.fwhm(nansum(CC,3), 'method', 'filters'); 
                
                obj.seeing = nanmedian(w(:)).*obj.cat.head.SCALE; 
                
            catch ME
                rethrow(ME); 
            end
            
            try % find the number of blinks 
            
                obj.num_blinks = 0; 
                state = 0; 

                first_blink_frame = 0;

                for ii = 1:length(f2)

                    if f2(ii)>obj.flux_thresh*s2 % turn on 

                        if state==0 % if it was off, this counts as a blink
                            obj.num_blinks = obj.num_blinks + 1; 
                        end

                        if first_blink_frame==0
                            first_blink_frame = ii; % log when was the first blink (not necessarily the same as the peak)
                        end

                        state = 1; 

                    end

                    if f2(ii)<obj.side_thresh*s2 % turn off
                        state = 0; 
                    end

                end
                
            catch ME
                rethrow(ME);
%                 warning(ME.getReport); 
            end
            
            try % figure out if the peak pixel is moving
                
                obj.peak_region = obj.frame_index; 
                
                for ii = 1:length(f2) % go back
                    
                    idx = obj.frame_index - ii;
                    
                    if idx<1 || f2(idx)<obj.side_thresh*s2
                        break;
                    end
                    
                    obj.peak_region(end+1) = idx; 
                    
                end
                
                for ii = 1:length(f2) % go forward
                    
                    idx = obj.frame_index + ii;
                    
                    if idx>length(f2) || f2(idx)<obj.side_thresh*s2
                        break;
                    end
                    
                    obj.peak_region(end+1) = idx; 
                    
                end
                
                obj.peak_region = sort(obj.peak_region); 
                
                counter = 1;
                
                for ii = 1:length(obj.peak_region)
                    for jj = 1:length(obj.peak_region)
                        idx1 = obj.peak_region(ii); 
                        idx2 = obj.peak_region(jj); 
                        dr(counter) = sqrt( (x(idx2)-x(idx1)).^2 + (y(idx2)-y(idx1)).^2 );
                        counter = counter + 1; 
                    end
                end
                
                obj.offset_shift = max(dr); 
                
            catch ME
                rethrow(ME);
%                 warning(ME.getReport); 
            end
            
            try % find the time when this occured 
            
                t_end = h5readatt(filename, '/timestamps', 't_end'); 
                t_stamp = h5readatt(filename, '/timestamps', 't_end_stamp');

                obj.peak_time = util.text.str2time(t_end) + seconds(obj.timestamps(obj.frame_index) - t_stamp); 

                if first_blink_frame
                    obj.first_peak_time = util.text.str2time(t_end) + seconds(obj.timestamps(first_blink_frame) - t_stamp); 
                else
                    obj.first_peak_time = obj.peak_time; 
                end
            
            catch ME
                warning(ME.getReport); 
            end
            
            try % check if the flare is streaked
                obj.is_streaked = obj.checkStreaked; 
            catch ME
                rethrow(ME); 
            end
            
            obj.is_bad_pos = nnz(isnan(obj.image))>20; 
            
        end
        
        function [is_streaked, ps_max, s_max] = checkStreaked(obj)
            
            % first, zoom in on the middle of the streak (we don't have the
            % tools for short streak detection!)
%             [~, idx] = util.stat.max2(obj.image); % point source peak intensity 
%             ps_err = sqrt(ps_max); % very rough estimate of the error
            
%             c = floor(size(obj.image)./2)+1; % center of the image (2-element vector!)
%             I = circshift(obj.image, c-idx); 

            I = util.img.centering(obj.image); 
            I(isnan(I)) = util.stat.median2(obj.image); 
            I = filter2(util.img.gaussian2(1, 'norm', 2), I); % match filter with a plausible sized PSF
            I = util.img.crop2size(I, 7); % get the middle only
            
            ps_max = util.stat.max2(I); 
            
            R = radon(I)./sqrt(radon(ones(size(I)))); % the normalized Radon image (we can replace this with FRT at some point)
            s_max = util.stat.max2(R);
            
            is_streaked = s_max > ps_max*sqrt(2); % the radon image has a brighter peak than the regular image (with error margin!)
            
        end
        
        function out_vec = filterGeosats(obj_vec)
            
            out_vec = obj_vec([obj_vec.is_streaked]==0 & [obj_vec.is_bad_pos]==0); 
            
        end
            
        function groups = clusterFlares(obj_vec, speed)
            
            if nargin<2 || isempty(speed)
                speed = 7; % pixels/sec, should not be more than 15 "/sec (6.55 assuming Balor's 2.29" pixels)
            end
            
            t = [obj_vec.first_peak_time]; 
            [~,idx] = sort(t); 
            obj_vec = obj_vec(idx); % sort the objects in ascending time order
            
            groups = {}; 
            
            for ii = 1:length(obj_vec)
                
                if isempty(obj_vec(ii).pos) || isempty(obj_vec(ii).peak_time)
                    error('Flare number %d has no position/time...', ii); 
                end
                
                found = 0; 
                
                for jj = 1:length(groups) % scan all existing groups
                    
                    if checkSpeed(obj_vec(ii), groups{jj}(end), speed) % check if the last member of each group is close enough to this object
                        groups{jj}(end+1) = obj_vec(ii); 
                        found = 1;
                        break; % once added to a group, we can skip to next object
                    end
                    
                end
                
                if found==0 % if we did not hit a "break" we must start a new group
                    groups{end+1} = obj_vec(ii); 
                end
                
            end
            
            for ii = 1:length(groups)
                
                s = groups{ii}; 
                
                % relative position and time
                dt = seconds([s.first_peak_time] - s(1).first_peak_time); 
                xy = [s.pos]-s(1).pos; 
                dr = sqrt(sum(xy.^2)); 
                
                % remove duplicate flares
%                 duplicate_indices = 1 + find(diff(dt)<0.01 & diff(dr)<0.01);
                duplicate_indices = 1 + find(diff(dt)<0.2);
                
                dr(duplicate_indices) = [];
                dt(duplicate_indices) = [];
                s(duplicate_indices) = [];
                
                groups{ii} = s; 
                
                fr = util.fit.polyfit(dt,dr,'order',1); 
                
                for jj = 1:length(s)
                    s(jj).velocity = fr.coeffs(2);
                end
                
            end
            
            
        end
        
        function groups = multiVelocityClustering(obj_vec, velocities)
            
            if nargin<2 || isempty(velocities)
                velocities = [7 20 100]; 
            end
            
            groups = {}; 
            
            new_vec = obj_vec; 
            
            for ii = 1:length(velocities)
                
                if ~isempty(new_vec)
                    
                    g1 = clusterFlares(new_vec, velocities(ii)); % cluster on a high speed

                    N = cellfun(@length, g1); % size of the groups

                    idx = N > 1; % indices of groups with more than one glint

                    if any(idx)
                        groups = horzcat(groups, g1(idx)); 
                    end

                    g1 = g1(~idx); 

                    new_vec = horzcat(g1{:}); 
                
                end
                
            end
            
        end
        
        function val = checkSpeed(obj, other, speed)
            
            dx = obj.pos(1) - other.pos(1); 
            dy = obj.pos(2) - other.pos(2); 
            dt = seconds(obj.first_peak_time - other.first_peak_time); 
            
            v = abs(sqrt(dx.^2+dy.^2)./(dt+2)); % add a fudge factor of 2 seconds to the time delay 
            
            val = abs(dt)<600 && v<speed; % obj is close enough to be associated with "other"
            
            
        end
        
        function val = checkDuplicatePeaks(obj)
            
            if isempty(obj)
                
                val = [];
                
            elseif isscalar(obj)
            
                [mx1, idx] = nanmax(obj.flux);

                f = obj.flux; 
                f(idx) = NaN;

                mx2 = nanmax(f); 

                val = mx1==mx2; 

            else
                
                for ii = 1:length(obj)
                    val(ii) = obj(ii).checkDuplicatePeaks; 
                end
                
                val = reshape(val, size(obj)); 
                
            end
            
        end
        
    end
    
    methods % other utilities (like save)
        
        function save(obj, filename)
            
            if nargin<2 || isempty(filename)
                error('Must supply a file name!');
            end
            
            event = obj;
            
            save(filename, 'event'); 
            
        end
        
        function saveDialog(obj, filename)
            
            if nargin<2 || isempty(filename)
                
                filename = 'flare'; 

                run_name = ''; 

                if ~isempty(obj.folder)
                    [a, b] = fileparts(obj.folder);
                    [~, c] = fileparts(a);
                    run_name = [c '_' b];
                end
                
                if ~isempty(run_name)
                    filename = [filename '_' run_name];
                end
                
                filename = sprintf('%s_id%03d', filename, obj.serial);
                
            end
            
            [filepath,filename,ext] = fileparts(filename);
            
            if ~isempty(ext)
                filename = [filename ext];
            else
                filename = [filename '.mat']; 
            end
            
            if isempty(filepath)
                filepath = fullfile(getenv('DATA'), 'WFAST/saved/flares/'); 
            end
            
            [filename, filepath] = uiputfile(fullfile(filepath, filename)); 
            
            if ~isequal(filename, 0)
                obj.save(fullfile(filepath, filename)); 
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('index', []); 
            input.input_var('cutouts', 9, 'number'); 
            input.input_var('parent', []); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.parent)
                input.parent = gcf;
            end
            
            if ~isempty(input.index) % user input overrides all
                input.parent.UserData = input.index;
            elseif ~isempty(input.parent.UserData) % no user input, recover latest index
                input.index = input.parent.UserData;
                if input.index>length(obj)
                    input.index = 1;
                end
            else % the default is to start with the first index
                input.parent.UserData = 1;
                input.index = 1;
            end
            
            this = obj(input.index);
            
            clf(input.parent); 
            
            panel_filename = uipanel(input.parent, 'Units', 'Normalized', 'Position', [0 0.9 1 0.1]);
            
            uicontrol(panel_filename, 'Style', 'pushbutton', 'String', this.filename, ...
                'Units', 'Normalized', 'Position', [0 0 0.5 1], 'FontSize', 10); 
            
            uicontrol(panel_filename, 'Style', 'pushbutton', 'String', this.printSummary, ...
                'Units', 'Normalized', 'Position', [0.5 0 0.5 1], 'FontSize', 10); 
            
            panel_flux = uipanel(input.parent, 'Units', 'Normalized', 'Position', [0 0 0.5 0.9]); 
            
            uicontrol(panel_flux, 'Style', 'pushbutton', 'Units', 'Normalized', 'Position', [0.3 0.05 0.1 0.1], ...
                'String', '-', 'Callback', @obj.prevShow, 'FontSize', 16); 
            
            uicontrol(panel_flux, 'Style', 'edit', 'Units', 'Normalized', 'Position', [0.4 0.05 0.3 0.1], ...
                'String', sprintf('%d / %d ', input.index, length(obj)), 'Callback', @obj.chooseShow, 'FontSize', 16); 
            
            uicontrol(panel_flux, 'Style', 'pushbutton', 'Units', 'Normalized', 'Position', [0.7 0.05 0.1 0.1], ...
                'String', '+', 'Callback', @obj.nextShow, 'FontSize', 16); 
            
            uicontrol(panel_flux, 'Style', 'pushbutton', 'Units', 'Normalized', 'Position', [0.05 0.05 0.2 0.1], ...
                'String', 'save', 'Callback', @obj.saveCallback, 'FontSize', 16); 
            
            ax = axes('Parent', panel_flux, 'Position', [0.2 0.35 0.7 0.55]);
            
            plot(ax, this.flux, 'LineWidth', 2); 
            
            xlabel(ax, 'frame number'); 
            ylabel(ax, 'flux'); 
            title(ax, sprintf('peak = %4.2f | pos= %d %d | pix= %d | rms= %4.2f', this.peak, this.pos(1), this.pos(2), obj(input.index).num_pixels, sqrt(this.pixel_var))); 
            ax.FontSize = 14;
            
            panel_cutouts = uipanel(input.parent, 'Units', 'Normalized', 'Position', [0.5 0 0.5 0.9]); 
            
            ax = util.plot.show_cutouts(this.cutouts, 'parent', panel_cutouts, 'frame', this.frame_index, 'number', input.cutouts); 
            
            for ii = 1:length(ax)
                if ~isempty(ax{ii}.UserData) && ax{ii}.UserData==this.frame_index
                    
                    hold(ax{ii}, 'on'); 
                    x = size(this.image,2)/2 +0.5; % + this.offset(this.frame_index,1);
                    y = size(this.image,1)/2 +0.5; % + this.offset(this.frame_index,2);
                    plot(ax{ii}, x, y, 'ro', 'MarkerSize', 25); 
                    hold(ax{ii}, 'off'); 
                    
                end
            end
            
        end
        
        function nextShow(obj, hndl, ~)
            
            idx = hndl.Parent.Parent.UserData;
            
            idx = idx + 1;
            if idx>length(obj)
                idx = 1;
            end
            
            obj.show('parent', hndl.Parent.Parent, 'index', idx); 
            
        end
        
        function prevShow(obj, hndl, ~)
            
            idx = hndl.Parent.Parent.UserData;
            
            idx = idx - 1;
            if idx<1
                idx = length(obj);
            end
            
            obj.show('parent', hndl.Parent.Parent, 'index', idx); 
            
        end
        
        function chooseShow(obj, hndl, ~)
            
            idx = util.text.extract_numbers(hndl.String);
            
            if ~isempty(idx)
                idx = idx{1}; 
            end
            
            if ~isempty(idx)
                idx = idx(1);
            end
            
            if isnumeric(idx) && ~isempty(idx)
                obj.show('parent', hndl.Parent.Parent, 'index', idx); 
            end
            
        end
        
        function saveCallback(obj, hndl, ~)
            
            idx = hndl.Parent.Parent.UserData;
            
            obj(idx).saveDialog;
            
        end
        
        function plotPositions(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('scale', []); 
            input.input_var('axes', [], 'axis'); 
            input.input_var('font_size', 14); 
            input.input_var('marker', 10);
            input.input_var('serial', false); 
            input.scan_vars(varargin{:}); 
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            if isempty(input.scale)
                input.scale = 1; 
                units = 'pixel/s';
            else
                units = '"/s'; 
            end
            
            prev_state = input.axes.NextPlot;
            
            for ii = 1:length(obj)
                
                if input.serial
                    serial_str = sprintf('(%d) ', obj(ii).serial); 
                else
                    serial_str = ''; 
                end
                
                plot(input.axes, obj(ii).pos(1), obj(ii).pos(2), 'rx', 'MarkerSize', input.marker); 
                
                if ii==1 || ii==length(obj)
                    text(input.axes, obj(ii).pos(1), obj(ii).pos(2), sprintf('  %s%s', serial_str, datestr(obj(ii).peak_time, 'HH:MM:SS')), 'FontSize', input.font_size); 
                elseif ii==floor(length(obj)/2)+1
                    text(input.axes, obj(ii).pos(1), obj(ii).pos(2), sprintf('  %svelocity= %4.2f%s', serial_str, obj(1).velocity.*input.scale, units), 'FontSize', input.font_size); 
                else
                    text(input.axes, obj(ii).pos(1), obj(ii).pos(2), sprintf('  %s', serial_str), 'FontSize', input.font_size); 
                end
                
                input.axes.NextPlot = 'add'; 
                
            end
            
            input.axes.NextPlot = prev_state; 
            
            input.axes.XLim = [1 4104]; 
            input.axes.YLim = [1 4128]; 
            
%             title(sprintf('velocity= %4.2f', obj(1).velocity)); 
            
        end
        
    end    
    
end

