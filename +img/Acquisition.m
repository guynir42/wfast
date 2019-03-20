classdef Acquisition < file.AstroData

    properties(Transient=true)
        
        gui;
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
        photo@img.PhotoSimple;
        flux_buf@util.vec.CircularBuffer;
        af@obs.focus.AutoFocus;
        
        % output to file
        buf@file.BufferWheel;
        
        % utilities
        prog@util.sys.ProgressBar;
        
        runtime_buffer@util.vec.CircularBuffer;
        
    end
    
    properties % inputs/outputs
        
        frame_rate = NaN;
        sensor_temp = NaN;
        
        cutouts_proc;
        cutouts_sub;
        cutouts_bg_proc;
        
        stack_cutouts; 
        stack_cutouts_sub;
        stack_cutouts_bg;
        stack_sub;
        prev_stack;        
        ref_stack;
        ref_positions;
        
        full_lightcurves;
        
        adjust_pos; % latest adjustment to x and y
        average_width; % of all stars in the stack
        
%         focus_curve_pos;
%         focus_curve_width;
        
        batch_counter = 0;
        batch_index = 1;

    end
    
    properties % switches/controls
        
        run_name = 'run1';
        
        expT = 0.025;
        
        num_batches = 2;
        batch_size = 100;
        
        num_stars = 100;
        cut_size = 11;
        
        num_backgrounds = 20;
        cut_size_bg = 20;
        
        use_background = 1;
        use_refine_bg = 0;
        
        use_adjust_cutouts = 1; % use adjustments in software (not by moving the mount)
        
        use_mextractor = 0;
        use_arbitrary_pos = 0;
        
        use_autofocus = 1;
        use_rough_focus = 0;
        
        use_roi = 0;
        roi_x1 = 1;
        roi_x2 = 512;
        roi_y1 = 1;
        roi_y2 = 512;
        
        use_save = 0; % must change this when we are ready to really start
        use_triggered_save = 0;
        
        % display parameters
        use_show = 1;
        show_what = 'images'; % can choose "images" or "stack"
        num_rect_stars = 30;
        num_rect_bg = 30;
        
        use_flip = 0; % flip view by 180 degrees (for meridien flip)
        
        use_audio = 1;
        
        % size of samples used for autofocus
%         focus_batch_size = 10;
%         
%         focus_curve_range = 0.2;
%         focus_curve_step = 0.01;
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        buf_full; % camera's buffers are used for full-frame dump on triggers
        
    end
    
    properties(Hidden=true)
       
        brake_bit = 1; % when this is set to 1 (using the GUI, for example), the run stops. 
        
        default_run_name;
        
        default_num_batches;
        
        default_batch_size;
        default_num_stars;
        default_cut_size;
        default_num_backgrounds;
        default_cut_size_bg;
        default_expT;
        
        show_what_list = {'images', 'stack'};
        
        version = 1.01;
        
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
                obj.photo = img.PhotoSimple;
                obj.flux_buf = util.vec.CircularBuffer;
                obj.af = obs.focus.AutoFocus;
                
                obj.buf = file.BufferWheel;
                obj.buf.product_type = 'Cutouts';
                
                obj.runtime_buffer = util.vec.CircularBuffer;
                obj.runtime_buffer.titles = {'time', 'num_frames'};
                
                util.oop.save_defaults(obj); % make sure each default_XXX property is updated with the current XXX property value. 

                obj.pars = head.Parameters; % this also gives "pars" to all sub-objects
                
                obj.prog = util.sys.ProgressBar;
                obj.audio = util.sys.AudioControl;
                
                obj.src = obj.reader;
                
            end
            
            
        end
        
        function updateInputObjects(obj) % to be depricated
            
            obj.input_main_survey = obj.makeInputVars;
            obj.input_preview = obj.makeInputVars('num_batches', 1, 'batch_size', 1, 'expT', 1, 'num_stars', 0);
            obj.input_live = obj.makeInputVars('num_batches', 1e6, 'batch_size', 1, 'num_stars', 0);
            obj.input_focus = obj.makeInputVars('batch_size', obj.focus_batch_size, 'expT', 1, 'num_stars', 0);
            
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
            obj.batch_index = 1;

            obj.full_lightcurves = [];
            
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
%         
%         function val = get.num_stars(obj)
%             
%             if isempty(obj.clip)
%                 val = obj.num_stars;
%             else
%                 val = obj.clip.num_stars;
%             end
%             
%         end
%         
%         function val = get.cut_size(obj)
%             
%             if isempty(obj.clip)
%                 val = obj.cut_size;
%             else
%                 val = obj.clip.cut_size;
%             end
%             
%         end
%         
%         function val = get.num_backgrounds(obj)
%             
%             if isempty(obj.clip_bg)
%                 val = obj.num_backgrounds;
%             else
%                 val = obj.clip_bg.num_stars;
%             end
%             
%         end
%         
%         function val = get.cut_size_bg(obj)
%             
%             if isempty(obj.clip_bg)
%                 val = obj.cut_size_bg;
%             else
%                 val = obj.clip_bg.cut_size;
%             end
%             
%         end
%         
        function val = get.buf_full(obj)
            
            if isempty(obj.cam)
                val = [];
            else
                val = obj.cam.buffers;
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
        
        function set.batch_size(obj, val)
            
            obj.batch_size = val;
            
            % any other sub-objects? 
            
        end
        
%         function set.num_stars(obj, val)
%             
%             if ~isempty(obj.clip)
%                 obj.clip.num_stars = val;
%             end
%             
%             obj.num_stars = val;
%             
%         end
%         
%         function set.cut_size(obj, val)
%             
%             if ~isempty(obj.clip)
%                 obj.clip.cut_size = val;
%             end
%             
%             obj.cut_size = val;
%             
%         end
%         
%         function set.num_backgrounds(obj, val)
%             
%             if ~isempty(obj.clip_bg)
%                 obj.clip_bg.num_stars = val;
%             end
%             
%             obj.num_backgrounds = val;
%             
%         end
%         
%         function set.cut_size_bg(obj, val)
%             
%             if ~isempty(obj.clip_bg)
%                 obj.clip_bg.cut_size = val;
%             end
%             
%             obj.cut_size_bg = val;
%             
%         end
%         
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
            input.input_var('run_name', 'test_run', 'name');
            input.input_var('expT', [], 'T', 'exposure time');
            input.input_var('num_batches', [], 'Nbatches');
            input.input_var('use_cutouts', 1);
            input.input_var('batch_size', [], 'frames');
            input.input_var('num_stars', [], 'Nstars');
            input.input_var('cut_size', []);
            input.input_var('num_backgrounds', [], 'Nbackgrounds');
            input.input_var('cut_size_bg', []);
            input.input_var('use_background', []);
            input.input_var('use_refine_bg', []);
            input.input_var('use_adjust_cutouts', []);
            input.input_var('use_mextractor', []);
            input.input_var('use_arbitrary_pos', []);
            input.input_var('use_autofocus', []);
            input.input_var('use_rough_focus', []);
            input.input_var('use_save', [], 'save');
            input.input_var('use_trigger_save', []);
            input.input_var('use_show', [], 'show');
            input.input_var('use_audio', []);
            input.input_var('axes', [], 'axis');
            input.input_var('use_index', 1);
            input.input_var('debug_bit', []);
            
            input.scan_obj(obj);
            input.scan_vars(varargin{:});
            
        end
        
        function parseInput(obj, input)
            
            list = properties(input);
            
            for ii = 1:length(list)
            
                if isprop(obj, list{ii})
                    obj.(list{ii}) = input.(list{ii});
                end
                
            end
            
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
                        obj.sim = img.Simulator;
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
    
    methods % commands/calculations
        
        function new_run(obj, varargin)
            
            input = obj.makeInputVars(varargin);
            
            obj.parseInput(input);
            
            obj.reset;
            
        end
        
        function run(obj, varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'util.text.InputVars')
                input = varargin{1};
                input.scan_vars(varargin{2:end});
            else
                input = obj.makeInputVars(varargin{:});
            end
            
            if obj.debug_bit>1, disp(['Starting run "' input.run_name '" for ' num2str(input.num_batches) ' batches.']); end
            
            cleanup = onCleanup(@() obj.finishup(input));
            obj.startup(input);
            
            if input.use_index
                idx = obj.batch_index;
            else
                idx = 1;
            end
            
            for ii = idx:input.num_batches

                if obj.brake_bit
                    return;
                end
                
                obj.batch(input);
                
                obj.prog.showif(obj.batch_counter);
                
            end
            
        end
        
        function update(obj) % do we need this??
            
            
            
        end
        
        function startup(obj, input)
            
            if ~obj.cal.checkDark
                error('Cannot start a new run without loading darks into calibration object!');
            end
                
            obj.brake_bit = 0;
            
            obj.update;
            
            if input.use_save
                filename = obj.buf.getReadmeFilename;
                util.oop.save(obj, filename, 'name', 'camera'); 
            end
            
            if input.use_audio
                try obj.audio.playTakeForever; catch ME, warning(ME.getReport); end
            end
            
            if isa(obj.src, 'file.Reader')
                obj.src.num_frames_per_batch = input.batch_size;
                obj.src.num_files_per_batch = 1;
                % what if batch_size is bigger than 100??
            else
                obj.src.batch_size = input.batch_size;
            end
            
            if input.use_cutouts
                
                obj.clip.num_stars = input.num_stars;
                obj.clip.cut_size = input.cut_size;

                obj.clip_bg.num_stars = input.num_backgrounds;
                obj.clip_bg.cut_size = input.cut_size_bg;

            end
            
            obj.frame_rate = NaN;
            obj.sensor_temp = NaN;
            
            obj.src.startup(input.num_batches);
            
            obj.prog.start(input.num_batches);
            
        end
        
        function finishup(obj, input)
            
            obj.prog.finish;
            
            obj.src.finishup;
            
            if obj.debug_bit, disp(['Finished run "' input.run_name '" with ' num2str(obj.batch_counter) ' batches.']); end
            
            if input.use_audio
                try
                    obj.audio.playShowsOver;
                catch ME
                    warning(ME.getReport);
                end
            end
            
            obj.brake_bit = 1;
            
        end
        
        function batch(obj, input)
            
            if obj.src.is_finished
                obj.brake_bit = 1;
                return;
            end
            
            obj.prev_stack = obj.stack_sub; % keep one stack from last batch
            obj.clear;
            
            t = tic;
            
            obj.src.batch; % produce the data (from camera, file, or simulator)
            
            if obj.debug_bit>1, fprintf('Starting batch %d. Loaded %d images from "%s" source.\n', obj.batch_counter+1, size(obj.images,3), class(obj.src)); end
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.update; 
            end
            
            drawnow; % make sure commands to the GUI and other callbacks are noticed... 
            
            obj.copyFrom(obj.src); % get the data into this object
            
            if obj.use_roi
                
                S = size(obj.images);
                
                x1 = obj.roi_x1;
                if x1<1, x1 = 1; end
                if x1>S(2), x1 = S(2)-1; end
                
                x2 = obj.roi_x2;
                if x2<x1, x2 = x1+1; end
                if x2>S(2), x2 = S(2); end
                
                y1 = obj.roi_y1;
                if y1<1, y1 = 1; end
                if y1>S(1), y1 = S(1)-1; end
                
                y2 = obj.roi_y2;
                if y2<y1, y2 = y1+1; end
                if y2>S(1), y2 = S(1); end
                
                obj.images = obj.images(y1:y2, x1:x2,:);
                
                obj.cal.use_roi = 1;
                obj.cal.roi_x1 = x1;
                obj.cal.roi_x2 = x2;
                obj.cal.roi_y1 = y1;
                obj.cal.roi_y2 = y2;
                
            else
                obj.cal.use_roi = 0;
            end
            
            obj.calcStack(input);
            
            if input.use_cutouts
                obj.calcCutouts(input);
                obj.calcLightcurves(input);
                obj.calcTrigger(input);
            end
            
            if input.use_save
                obj.buf.input(obj);
                obj.buf.clearImages;
                obj.buf.save;
                obj.buf.nextBuffer;
                % check triggering then call the camera save mode
                
            end
            
            if ismethod(obj.src, 'getTemperature')
                obj.sensor_temp = obj.src.getTemperature;
            end
            
            obj.runtime_buffer.input([toc(t), size(obj.images,3)]);
            
            T = sum(obj.runtime_buffer.data(:,1));
            N = sum(obj.runtime_buffer.data(:,2));
            
            obj.frame_rate = N./T;
            
            if ismethod(obj.src, 'next')
                obj.src.next;
            end
            
            obj.batch_counter = obj.batch_counter + 1;
            
            if input.use_index
                obj.batch_index = obj.batch_index + 1;
            end
            
            if input.use_show
                obj.show;
            end
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.update; 
            end
            
            drawnow;
            
        end
        
        function single(obj, input)
            
            input = util.oop.full_copy(input);
            input.use_audio = 0;
            input.num_batches = 1;
            input.batch_size = 1;
            
            
            cleanup = onCleanup(@() obj.finishup(input));
            obj.startup(input);
            
            obj.batch;
            
        end
        
        function calcStack(obj, input)
            
            if nargin<2 || isempty(input)
                input = obj.input_main_survey;
            end
            
            % make the basic, calibrated stack image
            obj.num_sum = size(obj.images,3);
            obj.stack = obj.cal.input(sum(obj.images,3), 'sum', obj.num_sum);
            
            % make the background cutouts of the stack 
            if isempty(obj.clip_bg.positions) % only if we didn't already assign positions to the bg_cutouts
                obj.clip_bg.arbitraryPositions('im_size', size(obj.stack));
            end
            
            obj.stack_cutouts_bg = obj.clip_bg.input(obj.stack); % dim 1&2 are y&x, dim 3 is scalar, dim 4 is star number.

            if input.use_background
            
                obj.back.input(obj.stack_cutouts_bg, obj.clip_bg.positions);            
                % should also get variance from background object...
                
                B = obj.back.getImage(size(obj.stack));
                obj.stack_sub = obj.stack - B;

            else
                obj.stack_sub = obj.stack;
            end
            
            
            if obj.batch_counter==0
                obj.findStars(input); % this replaces the clipper findStars with somethinig better, or just runs it explicitely... 
            end
            
            obj.stack_cutouts = obj.clip.input(obj.stack); % if no positions are known, it will call "findStars"
            
            if input.use_background
                BC = obj.back.getPoints(obj.clip.positions);
                BC = permute(BC, [4,3,2,1]); % turn the column vector into a 4D vector
                obj.stack_cutouts_sub = obj.stack_cutouts - BC; 
            else
                obj.stack_cutouts_sub = obj.stack_cutouts;
            end
            
            obj.photo.input(obj.stack_cutouts_sub, 'moments', 1); % run photometry on the stack to verify flux and adjust positions
            
            obj.checkRealign(input);
            
            % store the latest fluxes from the stack cutouts
            obj.flux_buf.input(obj.photo.fluxes);
            
            % get the average width and offsets (weighted by the flux of each star...)
            M = mean(obj.photo.fluxes);
            obj.adjust_pos = [median(obj.photo.fluxes./M.*obj.photo.offsets_x, 'omitnan'), median(obj.photo.fluxes./M.*obj.photo.offsets_y, 'omitnan')];
            obj.adjust_pos(isnan(obj.adjust_pos)) = 0;
            
            obj.average_width = median(obj.photo.fluxes./M.*obj.photo.widths, 'omitnan'); % maybe find the average width of each image and not the stack??
            
            if obj.use_adjust_cutouts
                obj.clip.positions = obj.clip.positions + obj.adjust_pos;
            else
                % must send the average adjustment back to mount controller
            end
            
        end
        
        function findStars(obj, input)
            
            if input.use_arbitrary_pos
                obj.clip.arbitraryPositions; % maybe add some input parameters?
            elseif input.use_mextractor
                % add the code for mextractor+astrometry here
            else
                obj.clip.findStars(obj.stack_sub);                
            end
            
            obj.ref_stack = obj.stack_sub;
            obj.ref_positions = obj.clip.positions;
            
        end
        
        function checkRealign(obj, input)
            
            if ~is_empty(obj.flux_buf) % check that stars are still aligned properly... 
                
                mean_fluxes = obj.flux_buf.mean;
                mean_fluxes(mean_fluxes<=0) = NaN;

                new_fluxes = obj.photo.fluxes;
                new_fluxes(isnan(mean_fluxes)) = [];
                mean_fluxes(isnan(mean_fluxes)) = [];

                % maybe 0.5 is arbitrary and should be turned into parameters?
                if sum(new_fluxes<0.5*mean_fluxes)>0.5*numel(mean_fluxes) % lost half the flux in more than half the stars...

                    disp('Lost star positions, using quick_align');
                    
                    [~,shift] = util.img.quick_align(obj.stack_sub, obj.ref_stack);
                    obj.clip.positions = obj.ref_positions + flip(shift);

                    % this shift should also be reported back to mount controller? 

                    obj.stack_cutouts = obj.clip.input(obj.stack);

                    if input.use_background
                        BC = obj.back.getPoints(obj.clip.positions);
                        BC = permute(BC, [4,3,2,1]); % turn the column vector into a 4D vector
                        % should also get variance from background object...

                        obj.stack_cutouts_sub = obj.stack_cutouts - BC;
                    else
                        obj.stack_cutouts_sub = obj.stack_cutouts;
                    end

                    obj.photo.input(obj.stack_cutouts_sub, 'moments', 1); % run photometry on the stack to verify flux and adjust positions

                end

                % add second test and maybe quit the run if it fails...
                
            end
            
        end
        
        function calcCutouts(obj, input)
            
            if nargin<2 || isempty(input)
                input = obj.input_main_survey;
            end
            
            obj.cutouts = obj.clip.input(obj.images);
            obj.positions = obj.clip.positions;
            
            obj.cutouts_proc = obj.cal.input(obj.cutouts, 'clip', obj.clip);
            
            obj.cutouts_bg = obj.clip_bg.input(obj.images);
            obj.positions_bg = obj.clip_bg.positions;
            
            obj.cutouts_bg_proc = obj.cal.input(obj.cutouts_bg, 'clip', obj.clip_bg);
            
            B = obj.back.getPoints(obj.clip.positions); 
            B = permute(B, [4,3,2,1]); % turn the column vector into a 4D vector
            % can also get variance from background object...
            
            if obj.use_refine_bg
                % use bg_cutouts to calculate overall differences between
                % frames to correct for regional results from background
                % object (based on the stack). 
            end
            
            obj.cutouts_sub = obj.cutouts_proc - B;
            
        end
        
        function calcLightcurves(obj, input)
            
            if nargin<2 || isempty(input)
                input = obj.input_main_survey;
            end
            
            obj.photo.input('images', obj.cutouts_sub, 'timestamps', obj.timestamps); % add variance input? 
            
            obj.lightcurves = obj.photo.fluxes;
            obj.full_lightcurves = cat(1, obj.full_lightcurves, obj.lightcurves); % consider replacing this with some storage object...
%             obj.widths = obj.photo.widths;
            
        end
        
        function calcTrigger(obj, input)
            
            if nargin<2 || isempty(input)
                input = obj.input_main_survey;
            end
            
            % to be implemented!
        end
        
        function roi(obj, position)
            
            import util.text.cs;
            
            if cs(position, 'none')
                obj.use_roi = 0;
            elseif cs(position, 'center')
                obj.use_roi = 1;
                obj.roi_x1 = 2160-256; % must change this to width/height of source
                obj.roi_x2 = 2160+255;
                obj.roi_y1 = 2560-256;
                obj.roi_y2 = 2560+255;
            elseif cs(position, 'northwest')
                obj.use_roi = 1;
                obj.roi_x1 = 1;
                obj.roi_x2 = 512;
                obj.roi_y1 = 1;
                obj.roi_y2 = 512;
            elseif cs(position, 'northeast')
                obj.use_roi = 1;
                obj.roi_x1 = 2160-512+1;
                obj.roi_x2 = 2160;
                obj.roi_y1 = 1;
                obj.roi_y2 = 512;
            elseif cs(position, 'southeast')
                obj.use_roi = 1;
                obj.roi_x1 = 2160-512+1;
                obj.roi_x2 = 2160;
                obj.roi_y1 = 2560-512+1;
                obj.roi_y2 = 2560;
            elseif cs(position, 'southwest')
               obj.use_roi = 1;
               obj.roi_x1 = 1;
                obj.roi_x2 = 512;
                obj.roi_y1 = 2560-512+1;
                obj.roi_y2 = 2560;
            else
                error('Unknown ROI position "%s". Use "center", "none", "NorthEast" etc...', position);
            end
            
        end
        
    end
    
    methods % optional run commands
        
        function runPreview(obj, varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'util.text.InputVars')
                input = varargin{1};
                input.scan_vars(varargin{2:end});
            else
                input = obj.makeInputVars('num_batches', 1, 'batch_size', 1, 'expT', 1, 'num_stars', 0, ...
                    'use_save', 0, 'use_audio', 0, varargin{:});
            end
            
            obj.run(input); % take the same run-loop but with different input parameters
            
        end
        
        function runLive(obj, varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'util.text.InputVars')
                input = varargin{1};
                input.scan_vars(varargin{2:end});
            else
                input = obj.makeInputVars('num_batches', 1e6, 'batch_size', 1, 'num_stars', 0, ...
                    'use_save', 0, 'use_audio', 0, 'use_show', 1, varargin{:});
            end
            
            obj.run(input); % take the same run-loop but with different input parameters
            
        end
        
        function runFocus(obj, varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'util.text.InputVars')
                input = varargin{1};
                input.scan_vars(varargin{2:end});
            else
                input = obj.makeInputVars('batch_size', 1, 'expT', 3, 'num_stars', 5, 'cut_size', 20, 'use audio', 0, varargin{:});
            end
            
            if isempty(obj.cam) || isempty(obj.cam.focuser)
                error('must be connected to camera and focuser!');
            end
            
            obj.single(input);
            obj.findStarsFocus(input);
            
            obj.af.pos = obj.cam.focuser.pos + (-obj.af.range:obj.af.step:obj.af.range)';
            obj.af.width = NaN(length(obj.af.pos), input.num_stars);
            obj.af.weight = NaN(length(obj.af.pos), input.num_stars);
            obj.af.xy_pos = obj.clip.positions;
            
            input.num_batches = length(obj.focus_curve_pos);
            
            obj.cam.focuser.pos = obj.af.pos(1);
            pause(0.1);
            
            obj.startup(input);
            cleanup = onCleanup(@() obj.finishup(input));
                        
            for ii = 1:input.num_batches

                if obj.brake_bit
                    return;
                end
                
                try 
                    obj.cam.focuser.pos = obj.focus_curve_pos(ii);
                catch ME
                    warning(ME.getReport);
                end
                
                obj.batch(input);
                
                %%% calculate the widths/fluxes somehow
                
%                 s = SIM;
%                 s.Im = obj.stack;
%                 s = s.mextractor;
%                 [~,H] = s.curve_growth_psf;
%                 W = H.*2;
                
                obj.af.pos(ii) = obj.cam.focuser.pos;
                obj.af.width(ii,:) = obj.photo.widths;
                obj.af.weights(ii,:) = obj.photo.fluxes;
                
                obj.af.plot;
                
                drawnow;
                                
            end
            
            obj.af.calculate;
            
%             fr = fit(obj.focus_curve_pos, obj.focus_curve_width, 'poly2')
%             best_pos = -fr.p2./2./fr.p1
%             
%             if best_pos<obj.focus_curve_pos(1) || best_pos>obj.focus_curve_pos(2)
%                 error('focus best pos is out of range');
%             end
            
            fprintf('FOCUSER RESULTS: pos= %f | tip= %f | tilt= %f\n', obj.af.found_pos, obj.af.found_tip, obj.af.found_tilt);
            
            obj.cam.focuser.pos = obj.af.found_pos;
            
            if isprop(obj.cam.focuse, 'tip')
                obj.cam.focuser.tip = obj.af.found_tip;
            end
            
            if isprop(obj.cam.focuse, 'tilt')
                obj.cam.focuser.tilt = obj.af.found_tilt;
            end
            
            obj.af.plot;
            
            obj.clip.reset; % don't save these star positions! 
            
        end
        
        function findStarsFocus(obj, input) % find stars in order in five locations around the sensor
            
            I = obj.stack;
            S = size(I);
            C = input.cut_size;
            
            I = util.img.maskBadPixels(I);
            I = conv2(I, util.img.gaussian2(2), 'same'); % smoothing filter
            
            markers = round(S.*[1/3 2/3]); % divide the sensor to 1/3rds 
            
            mask{1} = false(S);
            mask{1}(markers(1,1):markers(1,2), markers(2,1):markers(2,2)) = 1; % only the central part is unmasked
            
            mask{2} = false(S);
            mask{2}(C:markers(1,1), C:markers(2,1)) = 1; % upper left corner
            
            mask{3} = false(S);
            mask{3}(markers(1,2):end-C+1, C:markers(2,1)) = 1; % lower left corner
            
            mask{4} = false(S);
            mask{4}(C:markers(1,1), markers(2,2):end-C+1) = 1; % upper right corner
            
            mask{5} = false(S);
            mask{5}(markers(1,2):end-C+1, markers(2,2):end-C+1) = 1; % lower right corner
            
            pos = zeros(input.num_stars, 2);
            
            for ii = 1:input.num_stars
                
                [~,idx] = util.stat.max2(I.*mask{mod(ii,5)}); % find the maximum in each masked area
                
                pos(ii,:) = flip(idx); % x then y!
                
                I(idx(1)-floor(C/2):idx(1)+floor(C/2)+1, idx(2)-floor(C/2):idx(2)+floor(C/2)+1) = NaN; % remove found stars
                
            end
            
            obj.clip.positions = pos;
            
            if obj.gui.check
                obj.show;
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function makeGUI(obj)
            
            if isempty(obj.gui)
                obj.gui = img.gui.AcqGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
        function show(obj, varargin)
            
            import util.text.cs;
            
            try
            
            input = util.text.InputVars;
            input.input_var('ax', [], 'axes', 'axis');
            input.scan_vars(varargin{:});
            
            if isempty(input.ax)
                if ~isempty(obj.gui) && obj.gui.check
                    input.ax = obj.gui.axes_image;
                else
                    input.ax = gca;
                end
            end
            
            if cs(obj.show_what, 'images')
                I = obj.images(:,:,end);
            elseif cs(obj.show_what, 'stack')
                I = obj.stack_sub;
            else
                error('Unknown option for "display_what". Use "images", "raw", or "stack"');
            end
            
            if obj.use_flip
                I = rot90(I,2);
            end
            
            util.plot.setImage(I, input.ax);
            
            obj.clip.showRectangles('num', obj.num_rect_stars, 'color', 'black', 'ax', input.ax, 'flip', obj.use_flip, 'delete', 1);
            obj.clip_bg.showRectangles('num', obj.num_rect_bg, 'color', 'red', 'ax', input.ax, 'flip', obj.use_flip, 'delete', 0);
            
            catch ME
                warning(ME.getReport);
            end
            
        end
        
    end    
    
end

