classdef Simulator < file.AstroData

    properties(Transient=true)
        
    end
    
    properties % objects
        
        head;
        src;
        focus; 
        
    end
    
    properties % inputs/outputs
        
    end
    
    properties % switches/controls
        
        batch_size = 100;
        im_size; 
        
        flux = 1e4;
        psf_sigma = 2;
        
        use_external_source = 1;
        use_focus_simulator = 1;
        use_source_noise = 0;
        use_background = 0;
        use_read_noise = 0;
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Simulator(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'img.Simulator')
                if obj.debug_bit>1, fprintf('Simulator copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit>1, fprintf('Simulator constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
    end
    
    methods % getters
        
        function val = is_finished(obj)
            
            if obj.use_external_source && ~isempty(obj.src)
                val = obj.src.is_finished;
            else
                val = 0;
            end
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function startup(obj, varargin)
            
            if obj.use_external_source && ~isempty(obj.src)
                obj.src.startup(varargin{:});
            end
            
        end
        
        function finishup(obj, varargin)
            
            if obj.use_external_source && ~isempty(obj.src)
                obj.src.finishup(varargin{:});
            end
            
        end
        
        function batch(obj)
            
            if obj.use_external_source && ~isempty(obj.src)
                obj.src.batch;
                obj.images = double(obj.src.images);
                obj.im_size = util.vec.imsize(obj.images);
            end
            
            if isempty(obj.im_size)
                error('Must specify an image size!');
            end
            
            if isempty(obj.images)
                obj.images = obj.makeSources;
            else
                obj.images = obj.images + obj.makeSources;
            end
            
            
            if obj.use_background
            
            end
                
            if obj.use_read_noise
                
            end
            
            obj.images = uint16(round(obj.images));
            
        end
        
        function next(obj)
            
            if obj.use_external_source && ~isempty(obj.src) && ismethod(obj.src, 'next')
                obj.src.next;
            end

        end
        
        function I = makeSources(obj)
            
            I = zeros(obj.im_size);
            center = floor(obj.im_size/2)+1;
            
            I(center(1),center(2)) = obj.flux; % replace flux and coordinates with something smarter
            
            obj.psfs = util.img.gaussian2(obj.psf_sigma, 'norm', 1);
            
            if obj.use_focus_simulator && ~isempty(obj.focus) && isa(obj.focus, 'obs.focus.Simulator')
                obj.psfs = filter2(obj.focus.makePSF, obj.psfs);
            end
            
            I = filter2(obj.psfs,I);
            
            I = repmat(I, [1,1,obj.batch_size]);
            
            if obj.use_source_noise
                
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

