classdef Analysis < file.AstroData

    properties(Transient=true)
        
        gui;
        audio@util.sys.AudioControl;
        
        aux_figure;
        
    end
    
    properties % objects
        
        pars@head.Parameters;
        reader@file.Reader;
        cal@img.Calibration;
        clip@img.Clipper;
        clip_bg@img.Clipper;
        back@img.Background;
        phot@img.Photometry;
        
%         light_original@img.Lightcurves;
%         light_basic@img.Lightcurves;
        lightcurves@img.Lightcurves;
%         light_gauss@img.Lightcurves;
%         light_cosqrt@img.Lightcurves;
        % light_fit@img.Lightcurves;
        
        model_psf@img.ModelPSF;
        
        finder@trig.Finder;
        
        prog@util.sys.ProgressBar;
        
        func; % any function that takes first argument this object and runs custom analysis
        
    end
    
    properties % inputs/outputs
        
        cutouts_proc;
        cutouts_sub;
        cutouts_bg_proc;
        
        stack_cutouts; 
        stack_cutouts_sub;
        stack_cutouts_bg;
        stack_proc;
        deep_stack;
        deep_stack_aligned;
        subtract_stack;
        
        prev_stack;
        
        FWHM; % latest measured full width half maximum
        
        batch_counter = 0;
        
    end
    
    properties % switches/controls
        
        use_background_stack = 1; % subtract b/g from the full-frame stack
        use_background_cutouts = 1; % subtract b/g from the cutouts (and stack cutouts!)
        use_refine_bg = 0; % need to figure out exactly how to do this
        
        use_save_fits = 0;
        
        use_audio = 0;
        
        display_num_rect_stars = 30;
        
        brake_bit = 1;
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        num_batches;
        
    end
    
    properties(Hidden=true)
       
        % these are used for backward compatibility with older versions of 
        % img.Clipper that made a slightly different cutout based on the 
        % same positions. This affects ONLY the cutting of CALIBRATION! 
        use_cutout_adjustment = 0; % turn adjustments on/off
        cutout_adjustment_pixels = 0; % how many pixels to push (back or forward) relative to today's positions
        use_cutout_adjustment_floor = 0; % use floor of positions before (instead of) using round(). 
        
        num_batches_limit;
        version = 1.01;
        
    end
    
    methods % constructor
        
        function obj = Analysis(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'img.Analysis')
                if obj.debug_bit, fprintf('Analysis copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Analysis constructor v%4.2f\n', obj.version); end
            
                obj.reader = file.Reader;
                obj.cal = img.Calibration;
                obj.cal.load;
                obj.clip = img.Clipper;
                obj.clip.use_adjust = 0; % this should be disabled and depricated!
                obj.clip_bg = img.Clipper;
                obj.clip_bg.use_adjust = 0; % this should be disabled and depricated!
                obj.back = img.Background;
                obj.phot = img.Photometry;
                obj.phot.use_basic = 0;
                obj.phot.use_aperture = 1;
                obj.phot.use_gaussian = 0;
%                 obj.light_original = img.Lightcurves; 
%                 obj.light_basic = img.Lightcurves; obj.light_basic.signal_method = 'square'; obj.light_basic.background_method = 'corners';
                obj.lightcurves = img.Lightcurves; obj.lightcurves.signal_method = 'aperture'; obj.lightcurves.background_method = 'annulus';
%                 obj.light_gauss = img.Lightcurves; obj.light_gauss.signal_method = 'gauss'; obj.light_gauss.background_method = 'annulus';
                
                obj.model_psf = img.ModelPSF;
                
                obj.finder = trig.Finder;
                obj.finder.setupKernels;
                
                obj.prog = util.sys.ProgressBar;
                obj.audio = util.sys.AudioControl;
                
                obj.pars = head.Parameters; % this also gives "pars" to all sub-objects
                
                util.oop.save_defaults(obj); % make sure each default_XXX property is updated with the current XXX property value. 
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            list = properties(obj);
            
            for ii = 1:length(list)
                
                if isobject(obj.(list{ii})) && ~isempty(obj.(list{ii})) && ismethod(obj.(list{ii}), 'reset') 
                    obj.(list{ii}).reset;
                end
                
            end
            
            obj.batch_counter = 0;

            obj.prev_stack = [];
            
        end
        
        function clear(obj)
            
            list = properties(obj);
            
            for ii = 1:length(list)
                
                if isobject(obj.(list{ii})) && ~isempty(obj.(list{ii})) && ismethod(obj.(list{ii}), 'clear') 
                    obj.(list{ii}).clear;
                end
                
            end
            
        end
        
    end
    
    methods % getters
        
        function val = get.num_batches(obj)
            
            val = min([obj.num_batches_limit, obj.reader.getNumBatches]);
            
        end
        
        function val = seeing(obj)
            
            if isempty(obj.pars)
                val = obj.FWHM.*1.24;
            else
                val = obj.FWHM.*obj.pars.SCALE;
            end
            
        end
        
    end
    
    methods % setters
        
        function set.pars(obj,val)
            
            obj.pars = val;
            
            list = properties(obj);
            
            for ii = 1:length(list)
                
                if isobject(obj.(list{ii})) && ~isempty(obj.(list{ii})) && isprop(obj.(list{ii}), 'pars') 
                    obj.(list{ii}).pars = val;
                end
                
            end
            
        end
        
        function set.brake_bit(obj, val)
            
            obj.brake_bit = val;
            obj.reader.brake_bit = val;
            
        end
    
        function set.num_batches(obj, val)
            
            obj.num_batches_limit = val;
            
        end
        
    end
    
    methods % calculations
        
        function chooseDir(obj, dirname)
            
            if nargin<2 || isempty(dirname)
                dirname = '';
            end
            
            if isempty(dirname)
                obj.reader.browseDir;
            else
                if ~obj.reader.dir.cd(dirname)
                    error('cannot find directory %s', dirname);
                end
            end
            
            obj.cal.reader_dark.dir.cd(obj.reader.dir.pwd);
            obj.cal.reader_dark.dir.cd('..');
            if obj.cal.reader_dark.dir.smart_cd('dark')
                obj.cal.reader_dark.loadFiles;
            end
            
            obj.cal.reader_flat.dir.cd(obj.reader.dir.pwd);
            obj.cal.reader_flat.dir.cd('..');
            if obj.cal.reader_flat.dir.smart_cd('flat')
                obj.cal.reader_flat.loadFiles;
            end
            
            obj.cal.load; 
            
        end
        
        function startup(obj)
            
            obj.brake_bit = 0;
            
            if obj.use_audio
                try obj.audio.playTakeForever; catch ME, warning(ME.getReport); end
            end
            
            obj.prog.start(obj.num_batches);
            
        end
        
        function finishup(obj)
            
            obj.prog.finish;
            
            obj.finder.finishup;
            
            obj.brake_bit = 1;
            
            if obj.use_audio
                try
                    obj.audio.playShowsOver;
                catch ME
                    warning(ME.getReport);
                end
            end
            
            if obj.debug_bit, disp(['Finished run with ' num2str(obj.batch_counter) ' batches.']); end
            
        end
        
        function run(obj)
            
            if ~obj.cal.checkDark
                error('Cannot start a new run without loading darks into calibration object!');
            end
            
            cleanup = onCleanup(@() obj.finishup);
            obj.startup;
            
            for ii = obj.batch_counter+1:obj.num_batches
                
                if obj.brake_bit
                    return;
                end
                
                obj.batch;
                
                if ~isempty(obj.func)
                    feval(obj.func, obj);
                end
                
                obj.prog.showif(ii);
                
                if ~isempty(obj.gui) && obj.gui.check
                    obj.gui.update; 
                end

                drawnow;
                
                obj.batch_counter = obj.batch_counter + 1;
                
            end
            
        end
        
        function batch(obj)
            
            %%%%%%%%%%%%%%%%%%%%% GET DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.reader.batch;
            obj.copyFrom(obj.reader); 
            obj.clip.positions = obj.positions;
            obj.clip.cut_size = size(obj.cutouts,1);
            
            if obj.use_cutout_adjustment
                
                if obj.use_cutout_adjustment_floor
                    obj.clip.positions = floor(obj.positions) + obj.cutout_adjustment_pixels;
                else
                    obj.clip.positions = obj.positions + obj.cutout_adjustment_pixels;
                end
                
            end
            
            obj.clip_bg.positions = obj.positions_bg;
            obj.clip_bg.cut_size = size(obj.cutouts_bg,1);
            
            if isempty(obj.positions_bg)
                obj.clip_bg.num_stars = 50;
                obj.clip_bg.cut_size = 20;
                obj.clip_bg.arbitraryPositions;
            end
            
            if isempty(obj.num_sum)
                obj.num_sum = size(obj.cutouts,3);
            end
            
            %%%%%%%%%%%%%%%%%%%%% STACK ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % calibrate the stack if needed
            if nnz(isnan(obj.stack)) % stack is already calibrated (has NaN values...)
                obj.stack_proc = obj.stack;
            else
                obj.stack_proc = obj.cal.input(obj.stack, 'sum', obj.num_sum);
            end
            
            % cutouts of the stack
            obj.stack_cutouts = obj.clip.input(obj.stack_proc); 
            obj.stack_cutouts_bg = obj.clip_bg.input(obj.stack_proc); 
            
            % make background model based on stack cutouts
            if obj.use_background_cutouts || obj.use_background_stack
                obj.back.input(obj.stack_cutouts_bg, obj.clip_bg.positions);
            end
            
            % subtract background from stack (full frame)
            if obj.use_background_stack
                B = obj.back.getImage(size(obj.stack));
                obj.stack_proc = obj.stack_proc - B;
            end
            
            % subtract background from cutouts stack
            if obj.use_background_cutouts % if we already used background subtraction on all pixels, why not just cutout from stack_proc??
                BC = obj.back.getPoints(obj.clip.positions);
                BC = permute(BC, [4,3,2,1]); % turn the column vector into a 4D vector
                obj.stack_cutouts_sub = obj.stack_cutouts - BC; 
            else
                obj.stack_cutouts_sub = obj.stack_cutouts;
            end
            
            %%%%%%%%%%%%%%%%%%%%% CUTOUT ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.cutouts_proc = obj.cal.input(obj.cutouts, 'clip', obj.clip);
            
            if ~isempty(obj.cutouts_bg)
                obj.cutouts_bg_proc = obj.cal.input(obj.cutouts_bg, 'clip', obj.clip_bg);
            end
            
            B = obj.back.getPoints(obj.clip.positions); 
            B = permute(B, [4,3,2,1]); % turn the column vector into a 4D vector
            % can also get variance from background object...
            
            if obj.use_refine_bg
                % use bg_cutouts to calculate overall differences between
                % frames to correct for regional results from background
                % object (based on the stack). 
            end
            
            obj.cutouts_sub = obj.cutouts_proc - B./obj.num_sum;
            
            %%%%%%%%%%%%%%%%%%%%% PHOTOMETRY ANALYSIS %%%%%%%%%%%%%%%%%%%%%
            
            obj.phot.input('images', obj.cutouts_sub, 'timestamps', obj.timestamps, ...
                'positions', obj.positions, 'variance', single(2.5)); % need to add the sky background too
            
            obj.lightcurves.getData(obj.phot);
            if obj.lightcurves.gui.check, obj.lightcurves.gui.update; end
            
            %%%%%%%%%%%%%%%%%%%%% PSF modeling %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.model_psf.input(obj.cutouts_sub, obj.phot.offsets_x, obj.phot.offsets_y);
            
            obj.FWHM = util.img.fwhm(obj.model_psf.stack);
            
            %%%%%%%%%%%%%%%%%%%%% Event finding %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            f = obj.phot.fluxes;
            e = obj.phot.errors;
            a = obj.phot.areas;
            b = obj.phot.backgrounds;
            v = obj.phot.variances;
            x = obj.phot.offsets_x;
            y = obj.phot.offsets_y;
            w = obj.phot.widths;
            p = obj.phot.bad_pixels;
            phot_pars = obj.phot.pars_struct; % maybe also give this to model_psf??
            
            r = [];
            g = [];
            
            if obj.phot.use_gaussian
                g = obj.phot.gauss_sigma;
            elseif obj.phot.use_aperture
                
                r = obj.phot.aperture;
            end
            
            obj.finder.input(f, e, a, b, v, x, y, w, p, r, g, ...
                obj.timestamps, obj.cutouts_proc, obj.positions, obj.stack_proc, ...
                obj.batch_counter+1, 'filename', obj.reader.this_filename, ...
                't_end', obj.t_end, 't_end_stamp', obj.t_end_stamp,...
                'used_background', obj.phot.use_backgrounds, 'pars', phot_pars);
            
            %%%%%%%%%%%%%%%%%%%% save FITS files of stacks %%%%%%%%%%%%%%%%
            
            if obj.use_save_fits
                
                [d, f] = fileparts(obj.reader.this_filename);
                
                d = strrep(d, ' (Weizmann Institute)', '');
                
                d = fullfile(d, 'FITS/');
                
                if ~exist(d, 'dir')
                    mkdir(d);
                end
                
                fullname = fullfile(d,[f,'.fits']);
                fprintf('Saving "stack_proc" in FITS file: %s\n', fullname);
                
                fitswrite(double(obj.stack_proc), fullname); 
                obj.pars.writeFITS(fullname, [], obj.num_sum);
                
            end
            
            %%%%%%%%%%%%%%%%%%%% Update GUI and show stuff %%%%%%%%%%%%%%%%
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.show('ax', obj.gui.axes_image);
            end
            
            if ~isempty(obj.aux_figure) && isvalid(obj.aux_figure)
                obj.showLastEvents('parent', obj.aux_figure);
            end
            
            drawnow;
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.scan_vars(varargin{:});
            
            if isempty(input.ax)
                input.ax = gca;
            end
            
            util.plot.setImage(obj.stack_proc, input.ax);

            obj.clip.showRectangles('color', 'black', 'ax', input.ax, 'delete', 1, 'text', 1, 'num', obj.display_num_rect_stars);
            obj.clip_bg.showRectangles('color', 'red', 'ax', input.ax, 'delete', 0, 'text', 0);
            
        end
        
        function showLastEvents(obj, varargin)
            
            if isempty(obj.finder.last_events)
                return;
            end
            
            input = util.text.InputVars;
            input.input_var('parent', []);
            input.scan_vars(varargin{:});
            
            if isempty(input.parent)
                input.parent = gcf;
            end
            
            for ii = 1:length(obj.finder.last_events)
                
                if ii>1, pause(2); end
                
                obj.finder.last_events(ii).show('parent', input.parent);
                
            end
            
        end
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = img.gui.AnalysisGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
    end    
    
end

