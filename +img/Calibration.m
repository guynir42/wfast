classdef Calibration < handle
% Makes and applies the dark+flat used for basic calibration of any images. 
% This class has 3 main functions:
% (1) Make "super-dark" and "super-flat" images. This is done by calling
%     the "reader_dark" and "reader_flat" file.Reader objects to loop on 
%     the dark and flat files and make an average. 
%     For darks, we keep track of mean and variance of each pixel, and use
%     them to make a bad pixel mask. We assume the dark images are the same
%     exposure times as the science images, so that dark-current is well
%     calibrated. The products are "dark_mean", "dark_var", and "dark_mask".
%     For flats, we are mainly interested in the average flux, normalized
%     so the "flat_field" is the relative flux of each pixel. 
%     To produce a "super-dark" use "makeDark" and for a "super-flat" use
%     "makeFlat". Each time this is completed a "calibration.mat" file is
%     saved in the both the dark and flat folders of the readers. This file
%     can be loaded to save time (use "load" or "load_browse").
%
% (2) Using the dark/flat results to calibrate images. This is done useing
%     the "input" method, which applies dark and flat (if available) to the
%     image. This includes removing the bias and dark current, removing 
%     known bad pixels, and normalizing the flat (if available). 
%     Optional calibrations include removing any residual bias with 
%     "use_subtract_median" and getting rid of any non-constant noise with 
%     "use_remove_noise" (based on util.img.maskBadPixels). 
%     Optional parameters to the "input" function include "num_sum" that
%     tells the Calibration object the image it has gotten is a sum of many
%     images (applies a multiplied "dark_mean"), and "cut_pos" which is a
%     handle to a img.Clipper object with cutout positions. This is used to
%     make sure the calibration cutouts are in the right location, as
%     explained below.
%
% (3) Applying calibration of cutout images using dark/flat cutouts. This
%     is done by keeping a copy of an img.Clipper object inside this
%     object. The clipper generates "dark_mean_cut" and "dark_mask_cut" and
%     "flat_field_cut" that are 4D matrices (with 3rd singleton dimension)
%     so they can easily be applied to the image cutouts. 
%     If the input image cutouts are in the same position as the last
%     batch, there is no need to change anything. But if the size or number
%     or position of the cutouts has changed (which is tested for by
%     "updateCutouts") then new versions of the calibration cutouts are
%     generated (this is what happens when we let the cutouts drift to
%     follow telescope position drift). 
%
% NOTE 1: when making dark/flat the object just calls "loop" with the
%         "mode" parameter set to 'dark' or 'flat'. This is set internally 
%         and does not need user input. 
% 
% NOTE 2: the dark_mask is lazy-loaded from dark_mean and dark_var, by
%         looking for pixels that have mean or variance to far from average. 
%         This cut is done using "dark_mask_sigma" that kills pixels with
%         mean N times above noise level, and "dark_mask_var_ratio" that
%         kills pixels in the highest variance percentile. 
%         Changing these parameters (or "dark_mean" or "dark_var") will
%         reset the "dark_mask" and "dark_mask_cut" and require them to be
%         re-calculated on demand.
% 
% A Graphic User Interface (GUI) for this class is also included in the 
% sub-package +gui. It is invoked using obj.makeGUI. 

    properties(Transient=true)
       
        gui@img.gui.CalGUI;
        
        audio@util.sys.AudioControl; % sounds an alert when starting/finishing making dark/flat
        
        prog@util.sys.ProgressBar; % prints how much time is left to make dark/flat
        
        dark_mean_transformed;
        dark_var_transformed;
        dark_mask_transformed;
        
        flat_field_transformed;
        flat_var_transformed;
        flat_mean_transformed;
        
        dark_mean_cut; % a 4D matrix with the cutouts of dark_mean (same positions as the data)
        dark_var_cut; % a 4D matrix with the cutouts of dark_var (same positions as the data)
        dark_mask_cut; % a 4D matrix with the cutouts of dark_mask (same positions as the data)
        flat_field_cut; % a 4D matrix with the cutouts of flat_field (same positions as the data)
        
    end
    
    properties % objects
       
        clip@img.Clipper; % this is used to get cutout dark/flat
        reader_dark@file.Reader; % reads the dark images to generate the dark frame
        reader_flat@file.Reader; % reads the flat images to generate the flat frame
        
    end
    
    properties % outputs
        
        num_darks = 0; % how many dark frames were summed
        dark_mean; % mean of all darks
        dark_var; % variance of all darks
        dark_mask; % mask containing the location of bad pixels
        date_dark = '';
        
        num_flats = 0; % how many flat frames were summed
        flat_field; % normalized flat field (relative intensities)
        flat_mean; % the mean of the flats (used for calculating gain, etc.)
        flat_var; % variance of the flats (used for calculating gain, etc.)
        date_flat = '';
        
        flat_pixel_mean; % used for calculating gain
        flat_pixel_var; % used for calculating gain
        pixel_gain; % used for calculating gain
        
        gain; % when gain calculation is performed
        
        num_pixels_removed = 0; % how many pixels are removed by "maskBadPixels"
        
    end
    
    properties (Dependent=true)

        height; % of the image
        width; % of the image
        
    end
    
    properties % switches and controls
        
        % switches for generating dark/flat
        use_median = 0;
        
        % switches for calibrating images
        use_flat = 1; % when calibrating images, choose if to divide by flat
        use_single = 1; % output everything in single precision... 
        
        use_noise_removal = 0; % when calibrating images, choose if to use maskBadPixels
        use_subtract_median = 0; % after subtracting dark, subtract residual median of each image
        use_gain = 0; % normalize by the gain
        
        use_interp_mask = 0; % interpolate over bad pixels
        use_conv_mask = 0; % replace bad pixels by mean of surrounding pixels (using convolution)
        replace_value = NaN;
        
        dark_mask_sigma = 5; % how hot a pixel must be (relative to dark noise) to be considered a "bad pixel"
        dark_mask_var_thresh = 120; % pixels with variance above this value are bad pixels. 
        dark_mask_var_sigma = 100; % how many "sigmas" above the mean variance value to cut (to be depricated)
        dark_mask_var_ratio = 0.99; % what fraction of pixel variance is considered "bad pixels" (to be depricated)
        
        use_remove_empty_frames = 1;
        
        % switches for making darks/flats (can we get rid of this??)
        use_calc_gain = 0;
        use_mv_points = 0; % choose if to collect mean and variance data when making flats (used for calculating gain)
        
        mode = ''; % keeps track if we are looping (calculating) dark or flat
        
        % transformations: 
        use_roi = 0;
        ROI = []; % in pixels, defined as [top, left, height, width]
        
        use_downsample = 0;
        downsampling = 2;
        
        % switches: 
        autosave = 1; % once dark or flat is calculated, save it as a .mat file
        brake_bit = 1; % stop in mid calculation loop
        debug_bit = 1;
        
    end
    
    properties (Hidden=true)
    
        default_dark_mask_sigma;
        default_dark_mask_var_thresh;
        default_dark_mask_var_sigma; % we're not using this anymore...
        default_dark_mask_var_ratio; % we're not using this anymore...
        
        version = 1.06;
        
    end
    
    methods % constructor
        
        function obj = Calibration(other)
            
            if nargin>0 && ~isempty(other) && isa(other, 'img.Calibration') % copy-constructor
                
                if obj.debug_bit, fprintf('Calibration copy-constructor v%4.2f\n', obj.version); end
                
                obj = util.oop.full_copy(other);
                
            else
                
                if obj.debug_bit, fprintf('Calibration constructor v%4.2f\n', obj.version); end
                
                util.oop.save_defaults(obj);
                
                obj.clip = img.Clipper;
                obj.reader_dark = file.Reader;
                obj.reader_flat = file.Reader;

                obj.audio = util.sys.AudioControl;
                obj.prog = util.sys.ProgressBar;

            end
            
        end
        
    end
    
    methods % reset and copy utilities
             
%         function reset(obj)
% 
%             obj.resetDark;
%             obj.resetFlat;
%             
%         end
        
        function resetDark(obj)
            
            obj.num_darks = 0;
            obj.dark_mean = [];
            obj.dark_var = [];
            obj.dark_mask = [];
            
        end
        
        function resetReaderDark(obj) % not sure we use this anymore...
            
            obj.reader_dark.reset;
            
        end
        
        function resetFlat(obj)
        
            obj.num_flats = 0;
            obj.flat_field = [];
            obj.flat_mean = [];
            obj.flat_var = [];
            
        end
        
        function resetPixelMeanVar(obj)
            
            obj.flat_pixel_mean = [];
            obj.flat_pixel_var = [];
            obj.pixel_gain = [];
            
        end
        
        function resetReaderFlat(obj)
            
            obj.reader_flat.reset;
            
        end
        
        function clear(obj)
            
            obj.num_pixels_removed = 0;
            
        end
        
    end
    
    methods % getters
        
        function val = get.height(obj)
            
            if isempty(obj.dark_mean)
                val = [];
            else            
                val = size(obj.dark_mean, 1);
            end
            
        end
        
        function val = get.width(obj)
            
            if isempty(obj.dark_mean)
                val = [];
            else            
                val = size(obj.dark_mean, 2);
            end
            
        end
        
        function M = get.dark_mask(obj)
            
            import util.stat.mean2;
            
            if isempty(obj.dark_mean) || isempty(obj.dark_var) % without dark mean/var cannot produce a mask
                M = [];
                return;
            end
            
            if isempty(obj.dark_mask)
                
                obj.calcDarkMask;

            end
            
            M = obj.dark_mask;
            
        end
        
        function M = get.dark_mask_cut(obj)
                        
            if isempty(obj.dark_mean) || isempty(obj.dark_var) || isempty(obj.clip) % without dark mean/var (or clipper) cannot make a dark_mask_cut
                M = [];
                return;
            end
                        
            if isempty(obj.dark_mask_cut) && ~isempty(obj.clip.positions)
                obj.dark_mask_cut = obj.clip.cutMatrix(obj.dark_mask); % lazy load the cutouts of the dark mask
            end
            
            M = obj.dark_mask_cut;
            
        end
            
        function res = checkDark(obj) % did we already calculate/load a dark
            
            res = ~isempty(obj.dark_mean);
            
        end
        
        function res = checkFlat(obj) % did we already calculate/load a dark
           
            res = ~isempty(obj.flat_mean);
            
        end
        
        function val = get.dark_mean_transformed(obj) % lazy loading! 
            
            if isempty(obj.dark_mean_transformed)
                
                obj.dark_mean_transformed = obj.transform(obj.dark_mean);
                
            end
            
            val = obj.dark_mean_transformed;
            
        end
        
        function val = get.dark_mask_transformed(obj) % lazy loading! 
            
            if isempty(obj.dark_mask_transformed)
                
                obj.dark_mask_transformed = obj.transform(obj.dark_mask);

            end
            
            val = obj.dark_mask_transformed;
            
        end
        
        function val = get.dark_var_transformed(obj) % lazy loading! 
            
            if isempty(obj.dark_var_transformed)
                
                obj.dark_var_transformed = obj.transform(obj.dark_var);
                
            end
            
            val = obj.dark_var_transformed;
            
        end
        
        function val = get.flat_mean_transformed(obj) % lazy loading! 
            
            if isempty(obj.flat_mean_transformed)
                
                obj.flat_mean_transformed = obj.transform(obj.flat_mean);

            end
            
            val = obj.flat_mean_transformed;
            
        end
        
        function val = get.flat_var_transformed(obj) % lazy loading! 
            
            if isempty(obj.flat_var_transformed)
                
                obj.flat_var_transformed = obj.transform(obj.flat_var);
                
            end
            
            val = obj.flat_var_transformed;
            
        end
        
        function val = get.flat_field_transformed(obj) % lazy loading! 
            
            if isempty(obj.flat_field_transformed)
                
                obj.flat_field_transformed = obj.transform(obj.flat_field);
                
            end
            
            val = obj.flat_field_transformed;
            
        end
        
        function M_trans = transform(obj, M)
            
            M_trans = M;
            
            if isempty(M), return; end
            
            if isa(M_trans, 'double') && obj.use_single
                M_trans = single(M_trans);
            end
            
            if obj.use_roi
                M_trans = M_trans(obj.roi_y1:obj.roi_y2, obj.roi_x1:obj.roi_x2);
            end
            
            if obj.use_downsample
                M_trans = util.img.downsample(M_trans, obj.downsampling);
            end
            
        end
        
        function val = roi_x1(obj)
            
            val = obj.ROI(2);
            
        end
        
        function val = roi_x2(obj)
            
            val = obj.ROI(2)+obj.ROI(4);
            
        end
        
        function val = roi_y1(obj)
            
            val = obj.ROI(1);
            
        end
        
        function val = roi_y2(obj)
            
            val = obj.ROI(1)+obj.ROI(3);
            
        end
        
    end
 
    methods % setters
                
        function set.dark_mean(obj, val)
           
            obj.dark_mean = val;
            obj.dark_mask = [];
            obj.dark_mask_cut = [];
            obj.dark_mean_transformed = [];
            
        end
        
        function set.dark_mean_transformed(obj, val)
            
            obj.dark_mean_transformed = val;
            obj.dark_mean_cut = [];
            
        end
        
        function set.dark_var(obj, val)
           
            obj.dark_var = val;
            obj.dark_mask = [];
            obj.dark_mask_cut = [];
            obj.dark_var_transformed = [];
            
        end
        
        function set.dark_var_transformed(obj, val)
            
            obj.dark_var_transformed = val;
            obj.dark_var_cut = [];
            
        end
            
        function set.dark_mask(obj, val)
            
            obj.dark_mask = val;
            obj.dark_mask_transformed = [];
            
        end
        
        function set.dark_mask_transformed(obj, val)
            
            obj.dark_mask_transformed = val;
            obj.dark_mask_cut = [];
            
        end
        
        function set.dark_mask_sigma(obj, sigma)
           
            obj.dark_mask_sigma = sigma;
            obj.dark_mask = [];
            obj.dark_mask_cut = [];
            
        end
        
        function set.dark_mask_var_thresh(obj, val)
           
            obj.dark_mask_var_thresh = val;
            obj.dark_mask = [];
            obj.dark_mask_cut = [];
            
        end
        
        function set.dark_mask_var_sigma(obj, sigma)
           
            obj.dark_mask_var_sigma = sigma;
            obj.dark_mask = [];
            obj.dark_mask_cut = [];
            
        end
        
        function set.dark_mask_var_ratio(obj, ratio)
            
            obj.dark_mask_var_ratio = ratio;
            obj.dark_mask = [];
            obj.dark_mask_cut = [];
           
        end
 
        function set.use_interp_mask(obj, val)
            
            obj.use_interp_mask = val;
            
            if val
                obj.use_conv_mask = 0;
            end
            
        end
        
        function set.use_conv_mask(obj, val)
            
            obj.use_conv_mask = val;
            
            if val
                obj.use_interp_mask = 0;
            end
            
        end
        
        function set.flat_mean(obj, val)
            
            obj.flat_mean = val;
            obj.flat_mean_transformed = [];
            
        end
        
        function set.flat_var(obj, val)
            
            obj.flat_var = val;
            obj.flat_var_transformed = [];
            
        end
        
        function set.flat_field(obj, val)
            
            obj.flat_field = val;
            obj.flat_field_transformed = [];
            
        end
        
        function set.flat_field_transformed(obj, val)
            
            obj.flat_field_transformed = val;
            obj.flat_field_cut = [];
            
        end
        
        function set.use_single(obj, val)
            
            if val==obj.use_single
                % pass
            else
                obj.use_single = val;
                
                obj.dark_mean_transformed = [];
                obj.dark_var_transformed = [];
                obj.flat_mean_transformed = [];
                obj.flat_var_transformed = [];
                obj.flat_field_transformed = [];
                
            end
            
        end
        
        function set.use_roi(obj, val)
            
            if val==obj.use_roi
                % pass
            else
                obj.use_roi = val;
                obj.dark_mean_transformed = [];
                obj.dark_var_transformed = [];
                obj.dark_mask_transformed = [];
                obj.flat_mean_transformed = [];
                obj.flat_var_transformed = [];
                obj.flat_field_transformed = [];
            end
            
        end
        
        function set.ROI(obj, val)
            
            if isequal(val, obj.ROI)
                % pass
            else
                obj.ROI = val;
                obj.dark_mean_transformed = [];
                obj.dark_var_transformed = [];
                obj.dark_mask_transformed = [];
                obj.flat_mean_transformed = [];
                obj.flat_var_transformed = [];
                obj.flat_field_transformed = [];
            end
            
        end
        
        function set.use_downsample(obj, val)
            
            if val==obj.use_downsample
                % pass
            else
                obj.use_downsample = val;
                obj.dark_mean_transformed = [];
                obj.dark_var_transformed = [];
                obj.dark_mask_transformed = [];
                obj.flat_mean_transformed = [];
                obj.flat_var_transformed = [];
                obj.flat_field_transformed = [];
            end
            
        end
        
        function set.downsampling(obj, val)
            
            if val==obj.downsampling
                % pass
            else
                obj.downsampling = val;
                obj.dark_mean_transformed = [];
                obj.dark_var_transformed = [];
                obj.dark_mask_transformed = [];
                obj.flat_mean_transformed = [];
                obj.flat_var_transformed = [];
                obj.flat_field_transformed = [];
            end
            
        end
        
    end
    
    methods % actions! 
        
        function browse(obj) % choose the date directory (will try to locate dark and flat dirs automatically)
        
            dir = obj.reader_dark.browseDir; 
            
            if ~isempty(dir) && ~isnumeric(dir) % successfully chose a directory
                
                obj.reader_flat.dir.cd(obj.reader_dark.dir); % move reader_flat to same directory
                
                if strfind(obj.reader_dark.dir.cwd, 'dark') % if we already specified a dark dir
                    obj.reader_flat.dir.up; % assume the level above has a "flat" dir, too
                else % need to go into a dark dir, if it exists
                    if ~obj.reader_dark.dir.smart_cd('dark')
                        warning('Cannot locate any "dark" folders inside %s', obj.reader_dark.cwd);
                    end
                end
                
                if ~obj.reader_flat.dir.smart_cd('flat')
                    warning('Cannot locate any "flat" folders inside %s', obj.reader_flat.cwd);
                end
                
            end
            
        end
            
        function dir = browseDark(obj) % set the directory for dark files
            
            dir = obj.reader_dark.browseDir;
            
        end
        
        function dir = browseFlat(obj) % set the directory for flat files
            
            dir = obj.reader_flat.browseDir;
            
        end
        
        function addDark(obj, images) % add a batch of dark files
        
            import util.stat.runningMean;
            import util.stat.runningSum;

            if nargin<2 || isempty(images)
%                 error('add some images to "addDark"!');
                return;
            end
                  
            if (~isempty(obj.height) && size(images,1)~=obj.height)
                error(['size mismatch. size(images)= ' num2str(size(images)) ' | height= ' num2str(obj.height) ' | width= ' num2str(obj.width)])
            end
            
            if (~isempty(obj.width) && size(images,2)~=obj.width)
                error(['size mismatch. size(images)= ' num2str(size(images)) ' | height= ' num2str(obj.height) ' | width= ' num2str(obj.width)])
            end
            
            current_frames = size(images,3);
            
            % dark mean
            if obj.use_median
                d_sum = current_frames.*double(median(images,3));
            else
                d_sum = sum(images,3);
            end
            
            obj.dark_mean = runningMean(obj.dark_mean, obj.num_darks, d_sum, current_frames);
            
            % dark variance
            d_var = var(double(images), [], 3)*current_frames;
            obj.dark_var = runningMean(obj.dark_var, obj.num_darks, d_var, current_frames);
                        
            % number of darks
            obj.num_darks = runningSum(obj.num_darks, current_frames);
                        
        end
        
        function addFlat(obj, images) % add a batch of flat files

            import util.stat.runningMean;            
            import util.stat.runningSum;
            
            if nargin<2 || isempty(images)
%                 warning('no data given to "addFlat"');
                return;
            end
            
            images = double(images);
            
            if ~isempty(obj.dark_mean) % subtract darks if available
                images = bsxfun(@minus, images, obj.dark_mean);
            end
                        
            current_frames = size(images,3);
            
            mask = obj.dark_mask;
            mask_3d = repmat(mask, [1 1 current_frames]); % make mask 3D matrix

            images(mask_3d) = NaN;

            norms = util.stat.mean2(images);
            
            % flat mean
            f_sum = sum(images,3);
            f_sum(mask) = 0; % maybe leave it as NaNs? or put in the mean values?
            obj.flat_mean = runningMean(obj.flat_mean, obj.num_flats, f_sum, current_frames);            
                        
            if obj.use_calc_gain % NOT FINISHED:
                
                % flat var
                m = util.stat.mean2(images); % 3D vector of the mean of each frame (to de-trend before calculating variance across frames)
                m = m./mean(m);
                f_var = var(images./m, [], 3)*current_frames;
                f_var(mask) = 0; % maybe leave it as NaNs? or put in the mean values?          
                obj.flat_var = runningMean(obj.flat_var, obj.num_flats, f_var, current_frames);

                obj.flat_pixel_mean = cat(3, obj.flat_pixel_mean, obj.flat_mean);
                obj.flat_pixel_var = cat(3, obj.flat_pixel_var, obj.flat_var);
                
            end
                        
            % the flat field is the flat mean, normalized
            f_field = sum(bsxfun(@rdivide, images, norms),3); % normalize to the mean of each image
            f_field(mask) = current_frames;
            
            obj.flat_field = runningMean(obj.flat_field, obj.num_flats, f_field, current_frames);
            
            % number of flat frames 
            obj.num_flats = runningSum(obj.num_flats, current_frames);
            
        end
                
        function loop(obj) % loops over all files in the reader (dark/flat based on "mode")
            
            import util.text.cs;
            
            if isempty(obj.mode)
                return;
            elseif cs(obj.mode, {'dark'})
                reader = obj.reader_dark;
            elseif cs(obj.mode, {'flat'})
                reader = obj.reader_flat;
            else
                error(['unknown calibration mode: ' obj.mode]);
            end
            
            assert(~isempty(reader), 'make sure to give a useful looper before running "loop"');            
            
            obj.prog.unpause;
            
            for ii = 1:1000000
                
                try 
                    
                    if obj.brake_bit
                        obj.prog.pause;
                        return;
                    end
                    
                    if reader.is_finished
                        break;
                    end

                    reader.batch;

                    if cs(obj.mode, 'dark')
                        obj.addDark(reader.images);
                        obj.date_dark = reader.t_start(1:10);
                    elseif cs(obj.mode, 'flat')
                        obj.addFlat(reader.images);
                        obj.date_flat = reader.t_start(1:10);
                    else
                        error(['unknown calibration mode: ' obj.mode]);
                    end

                    try 
                        obj.prog.showif(reader.counter_batches);
                    catch ME
                        warning(ME.getReport);
                    end
                    
                    if ~isempty(obj.gui)
                        obj.show;
                        obj.gui.update;
                    end
                catch ME
                    obj.brake_bit = 1;
                    obj.audio.playError;
                    rethrow(ME);
                end
                
            end
            
            if obj.brake_bit==0
                
                obj.brake_bit = 1;
                obj.mode = '';
                
                if obj.autosave
                    obj.save;
                end

                try
                    obj.prog.show(reader.getNumBatches);
                    obj.prog.pause;
                    obj.audio.playShowsOver;
                catch ME
                    warning(ME.getReport);
                end
            
            end
            
        end
        
        function makeDark(obj) % sets "mode=dark" and runs the loop
           
            if isempty(obj.reader_dark)
                error('cannot make dark without a reader_dark!');
            end
            
            obj.reader_dark.reset;
            
            if obj.reader_dark.is_finished
                error('no images found for making darks...');
            end
            
            obj.mode = 'dark';
            
            obj.brake_bit = 0;
            
            try
                obj.audio.playTakeForever;
                obj.prog.start(obj.reader_dark.getNumBatches);
            catch ME
                warning(ME.getReport)
            end
            
            obj.loop;            
            
        end
        
        function makeFlat(obj) % sets "mode=flat" and runs the loop
            
            if isempty(obj.reader_flat)
                error('cannot make flat without a reader_flat!');
            end
            
            obj.reader_flat.reset;
            
            if obj.reader_flat.is_finished
                error('no images found for making flats...');
            end
                        
            obj.mode = 'flat';
            
            obj.brake_bit = 0;
            
            if obj.use_mv_points
                obj.resetPixelMeanVar;
            end
                
            try
                obj.audio.playTakeForever;
                obj.prog.start(obj.reader_flat.getNumBatches);
            catch ME
                warning(ME.getReport)
            end
            
            obj.loop; 
            
        end
        
        function calcDarkMask(obj)
            
            % kill pixels with mean above/below the average dark_mean value
%             pix = bsxfun(@rdivide, bsxfun(@minus, obj.dark_mean, mean2(obj.dark_mean)), sqrt(mean2(obj.dark_var)))./obj.dark_mask_sigma;
%             pix(isnan(pix)) = 1;
%             M = logical(floor(abs(pix))); % if there are pixels whose means are over Xsigma from the norm
% 
%             % kill pixels with the highest variance
%             V = obj.dark_var(:);
%             V = sort(V);
%             index = ceil(obj.dark_mask_var_ratio*length(V)); 
%             V_thresh = V(index);
%             M(obj.dark_var>=V_thresh) = 1;

            if obj.debug_bit>1, disp('making dark mask'); end

            M = false(size(obj.dark_mean)); % start with an empty mask

            % find the pixels with unusual mean values
            m = obj.dark_mean(:);
            [mu,sig] = util.stat.sigma_clipping(m, 'dist', 'gauss', 'iterations', 5);
            
            idx = abs(m - mu)>sig*obj.dark_mask_sigma;
            M(idx) = true;
            
            if obj.debug_bit>2, fprintf('dark_mean stats: mu= %4.2f | sig= %4.2f | num_pix= %d | frac= %g\n', mu, sig, nnz(idx), nnz(idx)/numel(idx)); end
            
            % variance values
            v = obj.dark_var(:);
            M(v>obj.dark_mask_var_thresh) = true; % replace old method with simple cut value
            
            % find the pixels with unusual variance values (old method)
%             [mu,sig] = util.stat.sigma_clipping(v, 'dist', 'weibul', 'iterations', 5);
%             idx = abs(v - mu)>sig*obj.dark_mask_var_sigma;
%             M(idx) = true;
            
            if obj.debug_bit>2, fprintf('dark_mean stats: mu= %4.2f | sig= %4.2f | num_pix= %d | frac= %g\n', mu, sig, nnz(idx), nnz(idx)/numel(idx)); end
            
            obj.dark_mask = M;  
            
        end
        
        function updateCutouts(obj, clipper, cut_size) % check if we need to adjust the dark_mean_cut etc. 
            
            if nargin<3 || isempty(cut_size)
                cut_size = obj.clip.cut_size;
            end
            
            if isempty(clipper) || (isa(clipper, 'img.Clipper') && is_empty(clipper))
                return;
            end
                        
            if isempty(obj.dark_mean_cut) || isempty(obj.dark_var_cut) || ...
                    isempty(obj.dark_mask_cut) || isempty(obj.flat_field_cut) || isempty(obj.clip) ...
                    || ~obj.clip.is_equal(clipper) || cut_size~=obj.clip.cut_size
                
                if obj.debug_bit>1, disp('updating cutouts for calibration!'); end            
                
                if isa(clipper, 'img.Clipper')
                    obj.clip.positions = clipper.positions;
                    obj.clip.cut_size = clipper.cut_size;
                else
                    obj.clip.positions = clipper; % assumes clipper is just numeric matrix of cut_pos values! 
                    obj.clip.cut_size = cut_size; % if we only change the cut_size for some bizzare reason... 
                end

                % actually update the cutouts... 
                if obj.checkDark
                    obj.dark_mean_cut = obj.clip.cutMatrix(obj.dark_mean_transformed);
                    obj.dark_var_cut = obj.clip.cutMatrix(obj.dark_var_transformed);
                    obj.dark_mask_cut = obj.clip.cutMatrix(obj.dark_mask_transformed);
                end
                
                if obj.checkFlat
                    obj.flat_field_cut = obj.clip.cutMatrix(obj.flat_field_transformed);
                end
                
            end
            
        end
                
        function I = input(obj, varargin) % input images to be calibrated
        % usage: I = obj.input(varargin)
        % Can enter images as the first argument, or use input('images', I). 
        % If cutout images are given, should supply the Clipper object.
        % using the keyword "clipper". 
        % If summed image is given, should supply the number of images in
        % the sum using keyword "num_sum".
        % OTHER OPTIONS (these can be controller through object pars)
        % -flat: use_flat
        % -median: use_subtract_median
        % -gain: use_gain
        % -conv: use_conv_mask
        % -interp: use_interp_mask
        % -pixels: use_noise_removal (using maskBadPixels)
        % 
        
            import util.text.cs;
            
            if nargin==1, help('img.Calibration.input'); return; end
            
            assert(obj.checkDark, 'cannot do calibration without darks loaded...');
            
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('images', []);
            input.input_var('clipper', [], 'cut_pos');
            input.input_var('num_sum', 1, 'sum');
            input.input_var('flat', obj.use_flat);
            input.input_var('median', obj.use_subtract_median, 'subtract_median');
            input.input_var('gain', obj.use_gain);
            input.input_var('conv', obj.use_conv_mask, 'conv_mask');
            input.input_var('interp', obj.use_interp_mask, 'interp_mask');
            input.input_var('pixels', obj.use_noise_removal, 'noise_removal', 'bad_pixels');
            input.input_var('replace_value', obj.replace_value);
            input.scan_vars(varargin{:});
            
            assert(~isempty(input.images), 'cannot do calibration without any images!');
            if obj.use_single
                I = single(input.images);
            else
                I = double(input.images); 
            end
            
            %%%%%%%%%%%%%%%%% Take care of zero frames %%%%%%%%%%%%%%%%%%%%
            
            if obj.use_remove_empty_frames
                
                idx = squeeze(sum(util.stat.sum2(abs(I)),4)==0); % vector of true if the whole frame (dim 3) is zero
                I(:,:,idx,:) = NaN;
                
%                 I(I==0) = NaN;
                
            end
            
            if ~isempty(input.clipper) % this means we are working with cutouts
                
                if ~isa(input.clipper, 'img.Clipper') && ~(isnumeric(input.clipper) && size(input.clipper, 2)==2)
                    error('Must input the "clipper" object as an img.Clipper type object, or an Nx2 matrix. What is given is a %s "%s"',...
                        class(input.clipper), util.text.print_vec(size(input.clipper), 'x'));
                end

                obj.updateCutouts(input.clipper, size(input.images,1));
                D = obj.dark_mean_cut;
                DM = obj.dark_mask_cut;
                F = obj.flat_field_cut;
                
            else % full frame images
                
                D = obj.dark_mean_transformed;
                DM = obj.dark_mask_transformed;
                F = obj.flat_field_transformed;
                
            end
            
            DM = repmat(DM, [1,1,size(I,3),1]);
            
            if input.num_sum>1
                D = D*input.num_sum; % if we need to subtract multiple darks from a summed image...
                I = sum(I,3); % just making sure the input images are summed 
            end
            
            if size(I,1)~=size(D,1) || size(I,2)~=size(D,2)
%                 error(['size mismatch! ' util.text.print_size(I,D)]);
                error('size mismatch!');
            end
            
            %%%%%%%%%%%%%%%%% SUBTRACT DARK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            I = bsxfun(@minus, I, D); % subtract the mean dark image
            
            %%%%%%%%%%%%%%%%% DIVIDE BY FLAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if input.flat && obj.checkFlat
                I = bsxfun(@rdivide, I, F);
            end
                        
            %%%%%%%%%%%%%%%%% SUBTRACT MEDIAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
            if input.median
                I = bsxfun(@minus, I, util.stat.median2(I));
            end
            
            %%%%%%%%%%%%%%%%% DIVIDE BY GAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if input.gain && ~isempty(obj.gain)
                I = I./obj.gain;
            end
            
            %%%%%%%%%%%%%%%%% PIXEL MASK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if ~isempty(DM)
                                
                if obj.use_interp_mask % use regionfill for each bad pixel (the slow method!)
                    
                    I2 = zeros(size(I));
                    
                    for ii = 1:size(I,3)
                        for jj = 1:size(I,4)
                            I2(:,:,ii,jj) = regionfill(I(:,:,ii,jj), DM(:,:,1,jj));
                        end
                    end
                    
                    I = I2;
                    
                elseif obj.use_conv_mask % convolve with nearest neighbor average value
                                    
                    I2 = zeros(size(I)); 
                    k = ones(3);
                    k(2,2) = 0;
                    k = k./8;                    
                    
                    for ii = 1:size(I,3)
                        for jj = 1:size(I,4)
                            I2(:,:,ii,jj) = filter2(k, I(:,:,ii,jj));
                        end
                    end
                    
                    I = bsxfun(@times, ~logical(DM), I)+bsxfun(@times, logical(DM), I2); % replace the bad pixels with the interpolated values. 
                    
                else % just use zeros 
%                     I = bsxfun(@times, ~logical(DM), I);
                    I(DM) = input.replace_value;
                end
                                
            end
            
            %%%%%%%%%%%%%%%%% HANDLE BAD PIXELS %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if obj.use_noise_removal
                [I, pix_mask] = util.img.maskBadPixels(I);
                obj.num_pixels_removed = sum(pix_mask(:)); % number of noisy pixels detected
            end
            
        end
        
    end
    
    methods % utilities and calculations
               
        function calcGain(obj)
            
            import util.vec.tocolumn;
            
            if isempty(obj.flat_pixel_mean) || isempty(obj.flat_pixel_var)
                error('cannot calculate gain without pixel mean/var !');
            end

            % add progress bar
            obj.prog.start(size(obj.flat_pixel_mean,1));
            
            obj.pixel_gain = zeros(size(obj.flat_pixel_mean,1), size(obj.flat_pixel_mean,2)); % gain for each pixel... 
            
            for ii = 1:size(obj.flat_pixel_mean,1)
                
                obj.prog.showif(ii);
                
                for jj = 1:size(obj.flat_pixel_mean,2)
                
                    if obj.dark_mask(ii,jj), continue; end
                    
                    M = tocolumn(obj.flat_pixel_mean(ii,jj,:));
                    V = tocolumn(obj.flat_pixel_var(ii,jj,:));
                    fr = fit(M,V, 'poly1');
                    obj.pixel_gain(ii,jj) = fr.p1;
                    
                end
                
                plot(fr, M,V,'.');

            end

            obj.prog.finish; % should I input ii as well??
                        
            obj.gain = util.stat.median2(obj.pixel_gain);
            
        end
        
        function new = makeROIcalibration(obj, ROI_vec, offset) % create a new calibration object, for a smaller region of interest (to be depricated when using transformations)
            
            new = img.Calibration(obj);
            
            if ~obj.checkDark && ~obj.checkFlat                
                return;
            end
            
            if nargin<2 || isempty(ROI_vec)
                ROI_vec = 512;
            end
            
            if nargin<3 || isempty(offset)
                offset = 1;
            end
            
            if isscalar(ROI_vec)
                S = [1 1]*ROI_vec;            
            end
            
            left = floor((size(obj.dark_mean,2)-S(2))/2)+offset;
            right = left+S(2)-1;
            top = floor((size(obj.dark_mean,1)-S(1))/2)+offset;
            bottom = top+S(1)-1;

            if obj.checkDark
                new.dark_mean = obj.dark_mean(top:bottom, left:right);
                new.dark_var = obj.dark_var(top:bottom, left:right);
            end 
            
            if obj.checkFlat
                new.flat_mean = obj.flat_mean(top:bottom, left:right);
                new.flat_var = obj.flat_var(top:bottom, left:right);
                new.flat_field = obj.flat_field(top:bottom, left:right);
                
                if ~isempty(obj.flat_pixel_mean)
                    new.flat_pixel_mean = obj.flat_pixel_mean(top:bottom, left:right, :);
                end
                
                if ~isempty(obj.flat_pixel_var)
                    new.flat_pixel_var = obj.flat_pixel_var(top:bottom, left:right,:);
                end
                
            end
        end
        
    end
    
    methods % save and load
        
        function save(obj, filename, directory)
            
            import util.text.sa;
            
            if nargin<2 || isempty(filename)
                filename = 'calibration.mat';
            end
            
            if nargin<3 || isempty(directory)
                directory = [];
            end
            
            % don't need to save the loaded dark/flat files
            obj.reader_dark.clear;
            obj.reader_flat.clear;
            obj.clip.clear;
            deflate = 1;
            
            if isempty(directory)
                if isempty(obj.reader_dark) && isempty(obj.reader_flat)
                    util.oop.save(obj, filename, 'name', 'cal', 'deflate', deflate);
                else
                    if ~isempty(obj.reader_dark)
                        if obj.debug_bit, disp(['saving calibration data to: ' sa(obj.reader_dark.dir.pwd, filename)]); end
                        util.oop.save(obj, sa(obj.reader_dark.dir.pwd, filename), 'name', 'cal', 'deflate', deflate);
                    end
                    
                    if ~isempty(obj.reader_flat)
                        if obj.debug_bit, disp(['saving calibration data to: ' sa(obj.reader_flat.dir.pwd, filename)]); end
                        util.oop.save(obj, sa(obj.reader_flat.dir.pwd, filename), 'name', 'cal', 'deflate', deflate);
                    end
                    
                end
                
            else
                util.oop.save(obj, sa(directory, filename), 'name', 'cal', 'deflate', deflate);
            end
            
        end
           
        function save_browse(obj)
            
            import util.text.sa;
            
            if ~isempty(obj.reader_flat)
                [filename, directory] = uiputfile('', 'Save calibration file', sa(obj.reader_flat.dir.pwd, 'calibration.h5'));
            elseif ~isempty(obj.reader_dark)
                [filename, directory] = uiputfile('', 'Save calibration file', sa(obj.reader_dark.dir.pwd, 'calibration.h5'));
            else                
                [filename, directory] = uiputfile('', 'Save calibration file', 'calibration.h5');
            end
            
            if ischar(filename) && ischar(directory)
                [~, filename, ext] = fileparts(filename);
                if isempty(ext)
                    filename = [filename '.h5'];
                end
                obj.save(filename, directory);
            end
           
        end
        
        function load(obj, filename, directory)
                    
            import util.text.sa;
            
            if nargin<2 || isempty(filename)
                filename = 'calibration.mat';
            end
            
            if nargin<3 || isempty(directory)
                directory = [];
            end            
            
            fullname = '';
            
            if ~isempty(directory)
                fullname = util.text.sa(directory, filename);             
            else
                if isempty(fullname) && ~isempty(obj.reader_flat)
                    f_temp = sa(obj.reader_flat.dir.pwd, filename);
                    if ~isempty(f_temp) && exist(f_temp, 'file')
                        fullname = f_temp;
                    end
                end
                
                if isempty(fullname) && ~isempty(obj.reader_dark)
                    f_temp = sa(obj.reader_dark.dir.pwd, filename);
                    if ~isempty(f_temp) && exist(f_temp, 'file')
                        fullname = f_temp;
                    end
                end
                
            end
            
            if isempty(fullname) || ~exist(fullname, 'file')
                disp('cannot find calibration file.');
            else
                if obj.debug_bit, disp(['Loading calibration from ' fullname]); end
%                 util.oop.copy_props(obj, util.oop.load(f_temp, 'classname', 'img.Calibration', 'name', 'cal'),1);
                temp = util.oop.load(fullname, 'classname', 'img.Calibration', 'name', 'cal');
                
                % copy the outputs
                obj.num_darks = temp.num_darks;
                obj.dark_mean = temp.dark_mean;
                obj.dark_var = temp.dark_var;
                
                obj.num_flats = temp.num_flats;
                obj.flat_field = temp.flat_field;
                obj.flat_mean = temp.flat_mean;
                obj.flat_var = temp.flat_var;

                obj.flat_pixel_mean = temp.flat_pixel_mean;
                obj.flat_pixel_var = temp.flat_pixel_var;
                obj.pixel_gain = temp.pixel_gain;

                obj.gain = temp.gain;

                obj.num_pixels_removed = temp.num_pixels_removed;

                if obj.gui.check
                    obj.show;
                end
                
            end
            
        end
        
        function load_browse(obj)
            
            import util.text.sa;
            
            if ~isempty(obj.reader_flat)
                [filename, directory] = uigetfile('', 'Load calibration file', sa(obj.reader_flat.dir.pwd, 'calibration.h5'));
            elseif ~isempty(obj.reader_dark)
                [filename, directory] = uigetfile('', 'Load calibration file', sa(obj.reader_dark.dir.pwd, 'calibration.h5'));
            else                
                [filename, directory] = uigetfile('', 'Load calibration file', 'calibration.h5');
            end
            
            if ischar(filename) && ischar(directory)
                [~, filename, ext] = fileparts(filename);
                if isempty(ext)
                    filename = [filename '.h5'];
                end
                obj.load(filename, directory);
            end
           
            if isempty(obj.audio)
                obj.audio = util.sys.AudioControl;
            end
            
            if isempty(obj.prog)
                obj.prog = util.sys.ProgressBar;
            end
            
        end
        
    end
    
    methods % plotting
                
        function show(obj, varargin)
           
            import util.plot.show;
            import util.plot.inner_title;            
            import util.stat.mean2;
            
            if ~obj.checkDark
                disp('no dark has been calculated yet...');
                return
            end
            
            if isempty(obj.gui) || ~obj.gui.check
                parent = gcf;
            else
                parent = obj.gui.panel_image;
            end
            
            delete(parent.Children);
            
            ax = axes('Parent', parent, 'Position', [0.0 0.5 0.5 0.5]);
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.axes_dark = ax;
            end
            
            if ~isempty(obj.dark_mean) 
                show(obj.dark_mean, 'fancy',0, 'autodyn', 'on', 'axes', ax, varargin{:});
                colorbar(ax, 'off');
                inner_title(ax, 'dark mean', 'color', 'red');
                inner_title(ax, ['mean= ' num2str(mean2(obj.dark_mean))], 'Position', 'bottom', 'Color', 'red');
            else
                axis(ax,'image');
            end
  
            ax = axes('Parent', parent, 'Position', [0.5 0.5 0.5 0.5]);
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.axes_var = ax;
            end
            
            if ~isempty(obj.dark_var)
                show(obj.dark_var(:,:,1,1), 'fancy',0, 'autodyn', 'on', 'axes', ax, varargin{:});
                colorbar(ax, 'off');
                inner_title(ax, 'dark variance', 'color', 'red');
                inner_title(ax, ['mean= ' num2str(mean2(obj.dark_var))], 'Position', 'bottom', 'Color', 'red');
            else
                axis(ax,'image');
            end
            
            ax = axes('Parent', parent, 'Position', [0.0 0.0 0.5 0.5]);
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.axes_mask = ax;
            end
                        
            if ~isempty(obj.dark_mask)
                show(obj.dark_mask, 'fancy', 0, 'ax', ax, varargin{:});
                inner_title(ax, 'dark mask', 'color', 'red');
                inner_title(ax, ['fraction= ' num2str(mean2(obj.dark_mask))], 'Position', 'bottom', 'Color', 'red');
            end
            
            axis(ax,'image');

            
            ax = axes('Parent', parent, 'Position', [0.5 0.0 0.5 0.5]);
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.axes_flat = ax;
            end
            
            if ~isempty(obj.flat_field)
                show(obj.flat_field, 'fancy',0, 'autodyn', 'on', 'axes', ax, varargin{:});
                colorbar(ax, 'off');
                inner_title(ax, 'flat field', 'color', 'red');
                inner_title(ax, ['mean= ' num2str(mean2(obj.flat_mean))], 'ax', ax, 'Position', 'bottom', 'Color', 'red');
            else
                axis(ax,'image');
            end
            
        end
        
        function showComparison(obj, images, varargin)
            
            import util.text.cs;
            import util.plot.show;
            import util.plot.inner_title;
            
            parent = [];
            
            for ii = 1:2:length(varargin)
                if cs(varargin{ii}, 'parent')
                    parent = varargin{ii+1};
                end
            end
            
            if isempty(parent)
                parent = gcf;
            end
            
            delete(parent.Children);
            
            ax1 = axes('Parent', parent, 'Position', [0.05 0.1 0.4 0.8]);
            show(images, 'ax', ax1, varargin{:});
            inner_title('no calibration');
            
            ax2 = axes('Parent', parent, 'Position', [0.55 0.1 0.4 0.8]);
            show(obj.input(images), 'ax', ax2, varargin{:});
            inner_title('with calibration');
            
        end
        
        function makeGUI(obj)
           
            if isempty(obj.gui)
                obj.gui = img.gui.CalGUI(obj);
            end
            
            obj.gui.makeGUI;
            obj.show;
            
        end
        
        function makeLooperDarkGUI(obj)
           
            if isempty(obj.reader_dark)
                obj.reader_dark = img.FileLooper;
            end
            
            obj.reader_dark.makeGUI;
            
        end
        
        function makeLooperFlatGUI(obj)
           
            if isempty(obj.reader_flat)
                obj.reader_flat = img.FileLooper;
            end
            
            obj.reader_flat.makeGUI;
            
        end
        
    end
    
end