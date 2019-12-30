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
        
        fwhm;
        
        m2x;
        m2y;
        mxy;
        maj_axis;
        min_axis;
        angle;
        
    end
    
    properties % switches/controls
        
        use_mex = 1;
        radius = 5;
        
        use_gaussian = 0;
        gauss_sigma = 5;
        
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
        
        function reset(obj)
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            obj.offsets_x = [];
            obj.offsets_y = [];
            obj.cutouts = [];
            obj.cutouts_shifted = [];
            obj.stack = []; 
            obj.mask = [];

%             obj.fwhm = [];

%             obj.m2x = [];
%             obj.m2y = [];
%             obj.mxy = [];
%             obj.maj_axis = [];
%             obj.min_axis = [];
%             obj.angle = [];
    
        end
        
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
            
            S = util.vec.imsize(obj.cutouts);
            
            if ~isempty(input.radius)
                obj.radius = input.radius;
            end
            
            dx = obj.offsets_x;
            dy = obj.offsets_y;
            
            % average offsets (I would prefer the flux-weighted average)
            Adx = nanmean(dx);
            Ady = nanmean(dy); 
            
            dx(isnan(dx)) = Adx;
            dy(isnan(dy)) = Ady;
            
            dx(abs(dx)>S(2)/2) = Adx;
            dy(abs(dy)>S(1)/2) = Ady;
            
            if isempty(dx) || isempty(dy)
                obj.cutouts_shifted = obj.cutouts;
            end
            
            if obj.use_mex
                obj.cutouts_shifted = util.img.shift(obj.cutouts, -dx, -dy);
            else
                
                obj.cutouts_shifted = zeros(size(obj.cutouts), 'like', obj.cutouts); % preallocate
                
                for ii = 1:size(obj.cutouts,3)
                    
                    for jj = 1:size(obj.cutouts,4)
                        if isnan(dx(ii,jj)) || isnan(dy(ii,jj))
                            obj.cutouts_shifted(:,:,ii,jj) = NaN(size(cutouts,1),size(cutouts,2), 'like', cutouts);
                        else
                            obj.cutouts_shifted(:,:,ii,jj) = util.img.imshift(obj.cutouts(:,:,ii,jj), -dy(ii,jj), -dx(ii,jj));
                        end
                    end
                    
                end
            end
            
            obj.calcStack;
            
        end
        
        function calcStack(obj)
            
            import util.stat.sum2;
            
            obj.stack = sum(sum(obj.cutouts_shifted,3, 'omitnan'),4, 'omitnan');
            
            if obj.use_gaussian
                obj.stack = obj.stack.*util.img.gaussian2(obj.gauss_sigma, 'size', [size(obj.cutouts,1),size(obj.cutouts,2)]); 
            end
            
            obj.mask = util.img.ellipse('radius', obj.radius, 'size', [size(obj.cutouts,1), size(obj.cutouts,2)]);
            
            c = size(obj.cutouts); c = c(1:2);
            [X, Y] = meshgrid((1:c(2))-floor(c(2)/2)-1, (1:c(1))-floor(c(1)/2)-1); 
            
            I = obj.stack.*obj.mask;
            
            S = sum2(I);
            
            obj.m2x = sum2(X.^2.*I)./S;
            obj.m2y = sum2(Y.^2.*I)./S;
            obj.mxy = sum2(X.*Y.*I)./S;
            
            M = [obj.m2x obj.mxy; obj.mxy obj.m2y];
            
            [R,E] = eig(M);
            
            obj.angle = asind(R(4));
            obj.maj_axis = max(sqrt(diag(E)));
            obj.min_axis = min(sqrt(diag(E)));
            
%             obj.fwhm = sqrt(mean([obj.m2x,obj.m2y])).*2.355;
            obj.fwhm = util.img.fwhm(I);
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

