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
        
        function input(obj, cutouts, varargin) 
            
            import util.stat.sum2;
            import util.text.cs;
            import util.text.parse_bool;
            
            use_moments = 0;
            
            for ii = 1:2:varargin
                key = varargin{ii};
                val = varargin{ii+1};
                
                if cs(key, 'moments', 'width', 'offsets')
                    use_moments = parse_bool(val);
                end
                
            end
            
            obj.cutouts = cutouts;
            
            cut_size = size(cutouts);
            cut_size = cut_size(1:2);
            
            if use_moments
                
                [X,Y] = meshgrid((1:cut_size(2))-floor(cut_size(2)/2)-1, 1:cut_size(1)-floor(cut_size(1)/2)-1); 

                S = sum2(cutouts);
                m1x = sum2(cutouts.*X)./S;
                m1y = sum2(cutouts.*Y)./S;
                m2x = sum2(cutouts.*(X-m1x).^2)./S;
                m2y = sum2(cutouts.*(Y-m1y).^2)./S;

                obj.widths = permute(sqrt(m2x+m2y), [3,4,2,1]);
                obj.offsets_x = permute(m1x, [3,4,2,1]);
                obj.offsets_y = permute(m1y, [3,4,2,1]);

            end
            
            % use the first moment to adjust positions inside the aperture?
            
            obj.aperture.tile_size = cut_size;
            
            obj.fluxes = permute(sum2(obj.aperture.mask.*cutouts), [3,4,2,1]); 
            obj.weights = repmat(obj.aperture.weight, size(obj.fluxes));
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end
    
end

