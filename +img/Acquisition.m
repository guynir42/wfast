classdef Acquisition < file.AstroData

    properties(Transient=true)
        
        gui@img.gui.PipeGUI;
        audio@util.sys.AudioControl;
        
    end
    
    properties % objects
        
        % general objects
        pars@head.Parameters;
        
        % input objects
        cam@obs.cam.CameraControl;
        reader@file.Reader;
        sim; % later add the class for the simulator        
        src; % can be camera, reader or simulator
        
        % image processing
        cal@img.Calibration;
        back@img.Background;
        clip@img.Clipper;
        clip_bg@img.Clipper;
        fit_psf@img.FitPSF;
        photo@img.PhotoSimple;
        flux_buf@util.vec.CircularBuffer;
        % output to file
        buf@file.BufferWheel;
        
        % utilities
        prog@util.time.ProgressBar;
        
        % various observational parameter inputs
        input_recording@util.text.InputVars;
        input_preview@util.text.InputVars;
        input_live@util.text.InputVars;
        input_focus@util.text.InputVars;
        
    end
    
    properties % inputs/outputs
        
        cutouts_sub;
        stack_cutouts; 
        stack_cutouts_sub;
        stack_cutouts_bg;
        adjust_pos; % latest adjustment to x and y
        average_width; % of all stars in the stack
        stack_rem; % do we need this???
        
        
    end
    
    properties % switches/controls
        
        use_adjust_cutouts = 1; % use adjustments in software (not by moving the mount)
        use_background_refine = 0;
        
        % display parameters
        display_what = 'image'; % can choose "image", "raw" or "stack"
        use_flip = 0; % flip view by 180 degrees (for meridien flip)
        
        % size of samples used for autofocus
        focus_num_batches = 30;
        focus_batch_size = 100;
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        % these are for recording only
        num_batches;
        batch_size;
        num_stars;
        expT;
        
    end
    
    properties(Hidden=true)
       
        brake_bit = 1; % when this is set to 1 (using the GUI, for example), the run stops. 
        
        default_num_batches = 2;
        default_batch_size = 100;
        default_num_stars = 10;
        default_expT = 0.025;
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Acquisition(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'Acquisition')
                if obj.debug_bit, fprintf('Acquisition copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Acquisition constructor v%4.2f\n', obj.version); end
                
                obj.reader = file.Reader;
                obj.sim; % fill this when we have a simulator
        
                obj.cal = img.Calibration;
                obj.back = img.Background;
                obj.clip = img.Clipper;
                obj.clip.use_adjust = 0; % adjust through photometry object?
                obj.clip_bg = img.Clipper;
                obj.clip_bg.use_adjust = 0;
                obj.fit_psf = img.FitPSF;
                obj.photo = img.PhotoSimple;
                obj.flux_buf = util.vec.CircularBuffer;
                
                obj.buf = file.BufferWheel;
                
                % initialize the paramter parser for different work modes
                obj.input_recording = obj.makeInputVars('save', 1);
                obj.input_preview = obj.makeInputVars('num_batches', 1, 'batch_size', 1, 'expT', 1, 'num_stars', 0);
                obj.input_live = obj.makeInputVars('num_batches', 1e6, 'batch_size', 1, 'num_stars', 0);
                obj.input_focus = obj.makeInputVars('num_batches', obj.focus_num_batches, 'batch_size', obj.focus_batch_size, 'expT', 1, 'num_stars', 0);
                
                obj.pars = head.Parameters; % this also gives "pars" to all sub-objects
                
                obj.src = obj.reader;
                
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
            
            val = obj.input_recording.num_batches;
            
        end
        
        function val = get.batch_size(obj)
            
            val = obj.input_recording.batch_size;
            
        end
        
        function val = get.num_stars(obj)
            
            val = obj.input_recording.num_stars;
            
        end
        
        function val = get.expT(obj)
            
            val = obj.input_recording.expT;
            
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
        
        function set.num_batches(obj, val)
            
            obj.input_recording.num_batches = val;
            
        end
        
        function set.batch_size(obj, val)
            
            obj.input_recording.batch_size = val;
            
        end
        
        function set.num_stars(obj, val)
            
            obj.input_recording.num_stars = val;
            
        end
        
        function set.expT(obj, val)
            
            obj.input_recording.expT = val;
            
        end
        
        function set.brake_bit(obj, val)
            
            obj.brake_bit = val;
            
            if isprop(obj.src, 'brake_bit')
                obj.src.brake_bit = val;
            end
            
        end
        
    end
    
    methods % utilities
        
        function input = makeInputVars(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('name', 'test_run', 'run_name');
            input.input_var('num_batches', obj.default_num_batches, 'Nbatches');
            input.input_var('batch_size', obj.default_batch_size, 'frames');
            input.input_var('num_stars', obj.default_num_stars, 'Nstars');
            input.input_var('expT', obj.default_expT, 'T', 'exposure time');
            input.input_var('save', 0, 'use_save');
            input.input_var('show', 0, 'use_show');
            input.input_var('use_random_pos', 0);
            input.input_var('use_mextractor', 0);
            input.input_var('use_sounds', 1);
            input.input_var('use_rough_focus', 1);
%             input.input_var('use
            input.scan_vars(varargin{:});
            
        end
        
        function chooseSource(obj, source)
            
            import util.text.cs;
            
            if nargin<2 || isempty(source)
                source = questdlg('Which source?', 'choose source', 'Reader', 'Camera', 'Simulator', 'Reader');
                if isempty(source)
                    return;
                end
                
                if cs(source, 'file reader', 'reader')
                    obj.chooseDir;
                end
                
            end
            
            if ischar(source)
                if cs(source, 'Zyla')
                    
                    if isempty(obj.cam) || ~isa(obj.cam.cam, 'obs.cam.ZylaControl')
                        obj.cam = obs.cam.CameraControl('zyla');
                        obj.cam.pars = obj.pars;
                    end
                    
                    obj.src = obj.cam;
                    
                elseif cs(source, 'Dhyana')
                    
                    if isempty(obj.cam) || ~isa(obj.cam.cam, 'obs.DhyanaControl')
                        obj.cam = obs.CameraControl('dhyana');
                        obj.cam.pars = obj.pars;
                    end
                    
                    obj.src = obj.cam;
                    
                elseif cs(source, 'simcamera')
                    
                    if isempty(obj.cam) || ~isa(obj.cam.cam, 'obs.cam.SimCamera')
                        obj.cam = obs.cam.CameraControl('sim');
                        obj.cam.pars = obj.pars;
                    end
                    obj.src = obj.cam;
                    
                    
                elseif cs(source, 'simulator')
                    
                    if isempty(obj.sim)
                        error('Simulator is not yet implemented!');
                        obj.sim = sim.Simulator;
                        obj.sim.pars = obj.pars;
                    end
                    
                    obj.src = obj.sim;

                elseif cs(source, {'file reader', 'reader'})
                    
                    if isempty(obj.reader)
                        obj.reader = file.Reader;
                        obj.cam.pars = obj.pars;
                    end
                    
                    obj.src = obj.reader;
                    
                elseif cs(source, 'camera')
                    if isempty(obj.cam)
                        obj.cam = obs.cam.CameraControl;
                        obj.cam.pars = obj.pars;
                    end
                    
                    obj.src = obj.cam;
                else
                    warning('unknown source "%s"', source);
                end
            else
                if isa(source, 'obs.cam.CameraControl') || isa(source, 'file.Reader') % add simulator check, too
                    obj.src = source;
                    obj.src.pars = obj.pars;
                else
                    warning(['unknown source class: ' class(source)]);
                end
            end
           
        end
        
        function chooseDir(obj, dirname)
            
            if isempty(obj.src)
                error('Cannot chooseDir with an empty source');
            end
            
            if nargin<2 || isempty(dirname)
                dirname = '';
            end
            
            if ~isa(obj.src, 'file.Reader')
                error('Cannot chooseDir for a source of type %s', class(obj.src));
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
        
    end
    
    methods % calculations
        
        function run(obj, varargin)
            
            if isempty(varargin) || ~isa(varargin{1}, 'util.text.InputVars')
                input = util.oop.full_copy(obj.input_recording);
                input.scan_vars(varargin{:});
            elseif isa(varargin{1}, 'util.text.InputVars')
                input = varargin{1};
                input.scan_vars(varargin{2:end});
            else
                input = obj.makeInputVars(varargin{:});
            end
            
            cleanup = onCleanup(@() obj.finishup(input));
            obj.startup(input);
            
            for ii = 1:obj.num_batches

                if obj.brake_bit
                    return;
                end
                
                obj.batch(input);
                
            end
            
        end
        
        function update(obj)
            
            
            
        end
        
        function startup(obj, input)
            
            obj.brake_bit = 0;
            
            if input.use_mextractor
                % add the code for mextractor+astrometry here
            end
            
            obj.update;
            
            if input.save
                filename = obj.buf.getReadmeFilename;
                util.oop.save(obj, filename, 'name', 'camera'); 
            end
            
            if input.use_audio
                try obj.audio.playTakeForever; catch ME, warning(ME.getReport); end
            end
            
            obj.src.startup(input.num_batches);
            
            obj.prog.start(input.num_batches);
            
        end
        
        function finishup(obj, input)
            
            obj.prog.finish;
            
            obj.src.finishup;
            
            obj.brake_bit = 1;
            
        end
        
        function batch(obj, input)
            
            if obj.src.is_finished
                obj.brake_bit = 1;
                return;
            end
            
            obj.src.batch; % produce the data (from camera, file, or simulator)
            
            obj.gui.update; drawnow; % make sure commands to the GUI and other callbacks are noticed... 
            
            obj.clear; % only clear data once a new batch is done successfully. 
            
            obj.copyFrom(obj.src); % get the data into this object
            
            obj.calcStack;
            obj.calcCutouts;
            obj.calcLightcurves;
            obj.calcTrigger;
            
            if input.save
                obj.buf.copyFrom(obj);
                obj.buf.save;
                
                % check triggering then call the camera save mode
                
            end
            
            obj.batch_counter = obj.batch_counter + 1;
            
            if input.show
                obj.show;
            end
            
        end
        
        function calcStack(obj)
            
            obj.num_sum = size(obj.images,3);
            obj.stack = obj.cal.input(sum(obj.images,3), 'sum', obj.num_sum);
                
            obj.stack_cutouts = obj.clip.input(obj.stack); % if no positions are known, it will call "findStars"
            
            if isempty(obj.bg_clip.positions) % only if we didn't already assign positions to the bg_cutouts
                obj.bg_clip.arbitraryPositions('im_size', size(obj.stack));
            end
            
            obj.stack_cutouts_bg = obj.clip_bg.input(obj.stack); % dim 1&2 are y&x, dim 3 is scalar, dim 4 is star number.
            
            obj.back.input(obj.stack_cutouts, obj.clip_bg.positions);
            
            B = obj.back.getPoints(obj.clip.positions);
            B = permute(B, [4,3,2,1]); % turn the column vector into a 4D vector
            % can also get variance from background object...
            
            obj.stack_cutouts_sub = obj.stack_cutouts - B; 
            
            obj.photo.input(obj.stack_cutouts_sub, 'moments', 1); % run photometry on the stack to verify flux and adjust positions
            
            obj.flux_buf.input(obj.photo.fluxes);
            
            % get the average width and offsets (weighted by the flux of each star...)
            obj.adjust_pos = [median(obj.photo.fluxes.*obj.photo.offsets_x), median(obj.photo.fluxes.*obj.photo.offsets_y)];
            obj.average_width = median(obj.photo.fluxes.*obj.photo.widths);
            
            if obj.use_adjust_cutouts
                obj.clip.positions = obj.clip.positions + obj.adjust_pos;
            else
                % must send the average adjustment back to mount controller
            end
            
        end
        
        function calcCutouts(obj)
            
            [obj.cutouts, obj.positions] = obj.clip.input(obj.images);
            
            obj.cutouts_proc = obj.cal.input(obj.cutouts, 'clip', obj.clip);
            
            [obj.bg_cutouts, obj.bg_positions] = obj.clip_bg.input(obj.image);
            
            obj.bg_cutouts_proc = obj.cal.input(obj.bg_cutouts, 'clip', obj.clip_bg);
            
            B = obj.back.getPoints(obj.clip.positions); 
            B = permute(B, [4,3,2,1]); % turn the column vector into a 4D vector
            % can also get variance from background object...
            
            if obj.use_background_refine
                % use bg_cutouts to calculate overall differences between
                % frames to correct for regional results from background
                % object (based on the stack). 
            end
            
            obj.cutouts_sub = obj.cutouts_proc - B;
            
        end
        
        function calcLightcurves(obj)
            
            obj.photo.input('images', obj.cutouts_sub, 'timestamps', obj.timestamps); % add variance input? 
            
            obj.lightcurves = obj.photo.fluxes;
            obj.widths = obj.photo.widths;
            
        end
        
        function calcTrigger(obj)
            % to be implemented!
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin)
            
            import util.text.cs;
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.scan_vars(varargin{:});
            
            if isempty(input.ax)
                if obj.gui.check
                    input.ax = obj.gui.axes_image;
                else
                    input.ax = gca;
                end
            end
            
            if cs(obj.display_what, 'images')
                I = obj.images_proc(:,:,end);
            elseif cs(obj.display_what, 'raw')
                I = obj.images(:,:,end);
            elseif cs(obj.display_what, 'stack')
                I = obj.stack;
            else
                error('Unknown option for "display_what". Use "images", "raw", or "stack"');
            end
            
            if obj.use_flip
                I = rot90(I,2);
            end
            
            util.plot.setImage(I(:,:,ii), ax);
            drawnow;
            
        end
        
    end    
    
end

