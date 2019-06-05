classdef Analysis < file.AstroData

    properties(Transient=true)
        
        gui;
        audio@util.sys.AudioControl;
        
    end
    
    properties % objects
        
        pars@head.Parameters;
        reader@file.Reader;
        cal@img.Calibration;
        clip@img.Clipper;
        clip_bg@img.Clipper;
        back@img.Background;
        phot@img.Photometry;
        
        light_original@img.Lightcurves;
        light_basic@img.Lightcurves;
        light_ap@img.Lightcurves;
        light_gauss@img.Lightcurves;
%         light_cosqrt@img.Lightcurves;
        % light_fit@img.Lightcurves;
        
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
        prev_stack;
        
        batch_counter = 0;
        
    end
    
    properties % switches/controls
        
        use_background_stack = 1; % subtract b/g from the full-frame stack
        use_background_cutouts = 1; % subtract b/g from the cutouts (and stack cutouts!)
        use_refine_bg = 0; % need to figure out exactly how to do this
        
        show_cutouts_cal = 0; % show calibrated cutouts
        
        use_audio = 0;
        
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
        version = 1.00;
        
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
                obj.light_original = img.Lightcurves;
                obj.light_basic = img.Lightcurves;
                obj.light_ap = img.Lightcurves;
                obj.light_gauss = img.Lightcurves;
                
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
                
                obj.prog.show(ii);
                
                obj.show;

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
            
            %%%%%%%%%%%%%%%%%%%%% LIGHTCURVE ANALYSIS %%%%%%%%%%%%%%%%%%%%%
            
            obj.phot.input('images', obj.cutouts_sub, 'timestamps', obj.timestamps, 'positions', obj.positions); 
            
            obj.light_original.input('fluxes', obj.fluxes, 'timestamps', obj.timestamps);
            obj.light_basic.getData(obj.phot, 'basic');
            if obj.light_basic.gui.check, obj.light_basic.gui.update; end
            
            obj.light_ap.getData(obj.phot, 'ap');
            if obj.light_ap.gui.check, obj.light_ap.gui.update; end
            
            obj.light_gauss.getData(obj.phot, 'gauss');
            if obj.light_gauss.gui.check, obj.light_gauss.gui.update; end
            
            drawnow;
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
            
            
        end
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = img.gui.AnalysisGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
    end    
    
end

