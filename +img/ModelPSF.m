classdef ModelPSF < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
        offsets_x;
        offsets_y;
        cutouts;
        cutouts_shifted;
        stack; 
        mask;
        
        
    end
    
    properties % switches/controls
        
        use_mex = 0;
        radius = 5;
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = ModelPSF(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'img.ModelPSF')
                if obj.debug_bit, fprintf('ModelPSF copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('ModelPSF constructor v%4.2f\n', obj.version); end
            
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
        
        function input(obj, varargin)
            
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('cutouts', [], 'images');
            input.input_var('offsets_x', [], 'dx', 9);
            input.input_var('offsets_y', [], 'dx', 9);
            input.input_var('radius', [], 'aperture'); 
            input.scan_vars(varargin{:});
            
            if isempty(input.cutouts)
                return;
            end
            
            obj.cutouts = input.cutouts;
            obj.offsets_x = input.offsets_x;
            obj.offsets_y = input.offsets_y;
            
            if ~isempty(input.radius)
                obj.radius = input.radius;
            end
            dx = obj.offsets_x;
            dx(isnan(dx)) = 0;
            dy = obj.offsets_y;
            dy(isnan(dy)) = 0;
                
            if obj.use_mex
                obj.cutouts_shifted = util.img.shift(obj.cutouts, -dx, -dy);
            else
                
                obj.cutouts_shifted = zeros(size(obj.cutouts), 'like', obj.cutouts); % preallocate
                
                for ii = 1:size(obj.cutouts,3)
                    
                    for jj = 1:size(obj.cutouts,4)
                        obj.cutouts_shifted(:,:,ii,jj) = util.img.imshift(obj.cutouts(:,:,ii,jj), -dy(ii,jj), -dx(ii,jj));
                    end
                    
                end
            end
            
            obj.stack = sum(sum(obj.cutouts_shifted,3, 'omitnan'),4, 'omitnan');
            
            obj.mask = util.img.ellipse('radius', obj.radius, 'size', [size(obj.cutouts,1), size(obj.cutouts,2)]);
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

