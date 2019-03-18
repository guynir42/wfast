classdef PhotoSimple < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        aperture@img.Aperture;
        
    end
    
    properties % inputs/outputs
        
        % inputs
        cutouts;
        
        % outputs
        fluxes;
        weights;
        offsets_x;
        offsets_y
        widths;
        
    end
    
    properties % switches/controls
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        ap_size; % how many pixels is the aperture?
        
    end
    
    properties(Hidden=true)
       
        default_ap_size = 3;
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = PhotoSimple(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'CLASS')
                if obj.debug_bit, fprintf('PhotoSimple copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('PhotoSimple constructor v%4.2f\n', obj.version); end
                
                obj.aperture = img.Aperture;
                obj.aperture.plateau_size = obj.default_ap_size;
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function clear(obj)
            
            obj.cutouts = [];
            obj.fluxes = [];
            obj.weights = [];
            obj.offsets_x = [];
            obj.offsets_y = [];
            obj.widths = [];
            
        end
        
    end
    
    methods % getters
        
        function val = get.ap_size(obj)
            
            val = obj.aperture.plateau_size;
            
        end
        
    end
    
    methods % setters
        
        function set.ap_size(obj, val)
            
            obj.aperture.plateau_size = val;
            
        end
        
    end
    
    methods % calculations
        
        function input(obj, varargin) 
            
            import util.stat.sum2;

            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('cutouts', [], 'images');
            input.input_var('moments', 1);
            input.scan_vars(varargin{:});
            
            obj.cutouts = input.cutouts;
            
            cut_size = size(input.cutouts);
            cut_size = cut_size(1:2);
            
            if input.moments
                
                [X,Y] = meshgrid((1:cut_size(2))-floor(cut_size(2)/2)-1, (1:cut_size(1))-floor(cut_size(1)/2)-1); 

                I = input.cutouts;
                I = I - util.stat.median2(I);
%                 I(I<0) = NaN;
                S = sum2(I);
                I2 = I;
                I2(I2<3.*util.stat.std2(I2)) = 0;
                S2 = sum2(I2);
                
                m1x = sum2(I.*X)./S;
                m1y = sum2(I.*Y)./S;
                m2x = sum2(I2.*(X-m1x).^2)./S2;
                m2y = sum2(I2.*(Y-m1y).^2)./S2;
                mxy = sum2(I2.*(X-m1x).*(Y-m1y))./S2;
                
                % make sure there are no S==0 elements
                m1x(S==0) = 0;
                m1y(S==0) = 0;
                m2x(S==0) = NaN;
                m2y(S==0) = NaN;
                mxy(S==0) = NaN;
                
                % make sure there are no negative second moments
                m2x(m2x<0) = NaN;
                m2y(m2y<0) = NaN;
                
                % should we ever reach such values??
                m1x(abs(m1x)>cut_size(2)) = NaN;
                m1y(abs(m1y)>cut_size(1)) = NaN;
                
                obj.offsets_x = permute(m1x, [3,4,2,1]);
                obj.offsets_y = permute(m1y, [3,4,2,1]);

                obj.widths = permute(sqrt(m2x+m2y), [3,4,2,1]); % should we add mxy too?
                
            end
            
            % use the first moment to adjust positions inside the aperture?
            
            obj.aperture.tile_size = cut_size;
            
            obj.fluxes = permute(sum2(obj.aperture.mask.*input.cutouts), [3,4,2,1]); 
            obj.weights = repmat(obj.aperture.weight, size(obj.fluxes));
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end
    
end

