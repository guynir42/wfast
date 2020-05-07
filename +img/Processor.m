classdef Processor < file.AstroData

    properties(Transient=true)
        
        gui;
        
    end
    
    properties % objects
        
        head@head.Header;
        cat@head.Catalog;
        
        reader@file.Reader;
        cal@img.Calibration;
        
        phot@img.Photometry;
        lightcurves@img.Lightcurves;
        
        % ... add streak detection maybe?
        
    end
    
    properties % inputs/outputs
        
    end
    
    properties % switches/controls
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Processor(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'img.Processor')
                if obj.debug_bit, fprintf('Processor copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Processor constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

