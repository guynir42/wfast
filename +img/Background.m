classdef Background < handle

    properties % objects
        
        interpolant@scatteredInterpolant; % interpolant model
        var_interp@scatteredInterpolant; % model for the variance
        
    end
    
    properties % inputs/outputs
        
        % inputs 
        im_size; % must give this to calculate some b/g models... 
        
        cutouts; % cutouts in areas not containing stars! 
        positions; % x then y, for each cutout center
        
        num_frames; % do we need this???
        
        
        % intermediate calculation of sample points... 
        backgrounds; % the amount of light, per pixel, measured in those points
        variances; % the noise level, per pixel, measured in those points
        
        % model results
        median_background;
        median_variance;
        mean_background;
        mean_variance;
        
        poly_coeffs_background;
        poly_coeffs_variance;

        fourier_coeffs_background;
        fourier_coeffs_variance;

        % outputs
        output_image; % full res background image.
        output_var_image; % full res variance image.
        query_positions; % do we want the backgrounds in specific points?
        output_backgrounds; % the calculated backgrounds in the query points. 
        
    end
    
    properties % switches/controls
        
        pixel_mode = 'median'; % can choose "mean" or "sigma clipping" (future) 
        
        model_type = 'median'; % can choose "median", "mean", "polynomial", "fourier" (future) or "interp"
        
        poly_number = 2;
        
        interp_mode = 'linear';
        extrap_mode = 'linear';
        
        num_sigmas = 3;
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Background(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'img.Background')
                if obj.debug_bit, fprintf('Background copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Background constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function clear(obj)
           
            obj.cutouts = [];
            obj.positions = [];
            obj.backgrounds = [];
            obj.variances = [];
            
            obj.output_image = [];
            obj.output_var_image = [];
            obj.query_positions = [];
            obj.output_backgrounds = [];
            
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
            input.input_var('cutouts',[]);
            input.input_var('positions', []); 
            input.input_var('num_frames', []);
            input.scan_vars(varargin{:});
            
            if isempty(input.cutouts)
                error('Must supply cutouts for background estimation!');
            end
            
            if isempty(input.positions)
                error('Must supply positions for background estimation!');
            end
            
            obj.clear;
            
            obj.cutouts = input.cutouts;
            obj.positions = input.positions;
            
            if isempty(input.num_frames)
                input.num_frames = size(input.cutouts,3);
            end
            
            obj.num_frames = input.num_frames;
            
            obj.makeBackgrounds;
            obj.updateModel;
            
        end
        
        function makeBackgrounds(obj)
            
            import util.text.cs;
            
            if cs(obj.pixel_mode, 'mean')
                obj.backgrounds = permute(util.stat.mean2(obj.cutouts), [4,3,1,2]);
                obj.variances = permute(util.stat.var2(obj.cutouts), [4,3,1,2]);
            elseif cs(obj.pixel_mode, 'median')
                obj.backgrounds = permute(util.stat.median2(obj.cutouts), [4,3,1,2]);
                obj.variances = permute(util.stat.var2(obj.cutouts), [4,3,1,2]); % need something like median only for variance...
            elseif cs(obj.pixel_mode, 'sigma clipping')
                for ii = 1:size(obj.cutouts,4)
                    for jj = 1:size(obj.cutouts,3)
                        C = obj.cutouts(:,:,jj,ii);
                        C = C(:);
                        [mu, sig] = util.stat.sigma_clipping(C, 'nsigma', obj.num_sigmas);
                        obj.backgrounds(ii,jj) = mu;
                        obj.variances(ii,jj) = sig.^2;
                    end
                end
            else
                error('Unknown pixel_mode "%s". Use "mean" or "median", "sigma_clipping", etc...', obj.pixel_mode);
            end
            
        end
        
        function updateModel(obj)
            
            import util.text.cs;
            
            if cs(obj.model_type, 'median')
                obj.median_background = util.stat.median2(obj.backgrounds);
                obj.median_variance = util.stat.median2(obj.variances);
            elseif cs(obj.model_type, 'mean')
                obj.median_background = util.stat.mean2(obj.backgrounds);
                obj.median_variance = util.stat.mean2(obj.variances);
            elseif cs(obj.model_type, 'polynomial')
                obj.updatePolynomial;
            elseif cs(obj.model_type, 'interp')
                obj.updateInterpolant;
            end
            
        end
        
        function updatePolynomial(obj)
            
            m = size(obj.positions,1); % number of measurements
            n = 1+obj.poly_number*2; % number of parameters (piston+polynomial in each axis)
            
            B = obj.backgrounds; % measured backgrounds
            w = 1./obj.variances; % weights
            
            A = ones(m,1); % design matrix! 

            x = obj.positions(:,1);
            y = obj.positions(:,2);

            if obj.poly_number==1
                A = [A x y];
            elseif obj.poly_number==2
                A = [A x y x.^2 y.^2 x.*y]; % add single cross term
            elseif obj.poly_number==3
                A = [A x y x.^2 y.^2 x.*y x.^3 y.^3 x.^2.*y x.*y.^2]; % add multiple cross terms
            else
                error('no support for higher order polynomials...');
            end
            
%             for ii = 1:obj.poly_number
%                 
%                 A = [A obj.positions.^ii];
%                 
%             end
            
            obj.poly_coeffs_background = lscov(A,B,w);
            
        end
        
        function updateInterpolant(obj)
            
            obj.interpolant = scatteredInterpolant;
            obj.var_interp = scatteredInterpolant;
            
            cut_num = size(obj.positions,1);
            
            if obj.num_frames==1
                
                b = obj.backgrounds';
                obj.interpolant.Points = obj.positions(~isnan(b),:);
                obj.interpolant.Values = b(~isnan(b));
                
                v = obj.variances';
                obj.var_interp.Points = obj.positions(~isnan(v),:);
                obj.var_interp.Values = v(~isnan(v));
                
            else % need to make sure this works or get rid of it...
                
                Px = repmat(obj.positions(:,1), [obj.num_frames, 1]);
                Py = repmat(obj.positions(:,2), [obj.num_frames, 1]);
                Pz = 1:cut_num:(obj.num_frames*cut_num);
                Pz = Pz';
                Pz = repmat(Pz, [cut_num, 1]);

                obj.interpolant.Points = [Px, Py, Pz];
                obj.interpolant.Values = reshape(obj.backgrounds, [size(obj.backgrounds,1).*size(obj.backgrounds,2), 1]);
            
            end
            
            obj.interpolant.Method = obj.interp_mode;
            obj.interpolant.ExtrapolationMethod = obj.extrap_mode;
            
            obj.var_interp.Method = obj.interp_mode;
            obj.var_interp.ExtrapolationMethod = obj.extrap_mode;
            
        end
        
        function bg = getPoints(obj, positions)
            
            import util.text.cs;
            
            if nargin<2 || isempty(positions)
                positions = obj.query_positions;
            else
                obj.query_positions = positions;
            end
            
            if cs(obj.model_type, 'median')
                obj.output_backgrounds = repmat(obj.median_background, [size(positions,1),1]);
            elseif cs(obj.model_type, 'mean')
                obj.output_backgrounds = repmat(obj.mean_background, [size(positions,1),1]);
            elseif cs(obj.model_type, 'poly')
                
                A = ones(size(positions,1),1);

                x = positions(:,1);
                y = positions(:,2);
                
                if obj.poly_number==1
                    A = [A x y];
                elseif obj.poly_number==2
                    A = [A x y x.^2 y.^2 x.*y]; % add single cross term
                elseif obj.poly_number==3
                    A = [A x y x.^2 y.^2 x.*y x.^3 y.^3 x.^2.*y x.*y.^2]; % add multiple cross terms
                else
                    error('no support for higher order polynomials...');
                end
                
%                 m = (length(obj.poly_coeffs_background)-1)/2;
%                 for ii = 1:m
%                     A = [A positions.^ii];
%                 end
%                 
                obj.output_backgrounds = A*obj.poly_coeffs_background;
                
            elseif cs(obj.model_type, 'fourier')
                error('Fourier model is not yet implemented!');
            elseif cs(obj.model_type, 'interp')

                num_cuts = size(positions,1);

                if obj.num_frames==1

                    obj.output_backgrounds = obj.interpolant(positions);

                else

                    Px = positions(:,1);
                    Py = positions(:,2);
                    Pz = 1:obj.num_frames;

                    obj.output_backgrounds = obj.interpolant({Px, Py, Pz});

                end

                obj.output_backgrounds = reshape(obj.output_backgrounds, [obj.num_frames, num_cuts]);

            else
                error('Unknown model_type "%s". Use "median" or "mean" or "polynomial" or "interp".');
            end
            
            if nargout>0
                bg = obj.output_backgrounds;
            end
            
        end
        
        function I = getImage(obj, im_size)
            
            import util.text.cs;
            
            if nargin<2 || isempty(im_size)
                im_size = obj.im_size;
            end
            
            im_size = util.vec.imsize(im_size);
            
            if cs(obj.model_type, 'median')
                obj.output_image = ones(im_size).*obj.median_background;
            elseif cs(obj.model_type, 'mean')
                obj.output_image = ones(im_size).*obj.mean_background;
            elseif cs(obj.model_type, 'polynomial')
                
                [X,Y] = meshgrid(1:im_size(2), 1:im_size(1));
                m = (length(obj.poly_coeffs_background)-1)/2;
                c = obj.poly_coeffs_background;
                
                obj.output_image = c(1);
                
                if obj.poly_number>=1
                    obj.output_image = obj.output_image + c(2).*X + c(3).*Y;
                end
                
                if obj.poly_number>=2
                    obj.output_image = obj.output_image + c(4).*X.^2 + c(5).*Y.^2 + c(6).*X.*Y;
                end
                
                if obj.poly_number>=3
                    obj.output_image = obj.output_image + c(7).*X.^3 + c(8).*Y.^3 + c(9).*X.^2.*Y + c(10).*X.*Y.^2;
                end
                
                
%                 for ii = 1:m
%                     
%                     obj.output_image = obj.output_image + obj.poly_coeffs_background(2.*m).*X.^m + obj.poly_coeffs_background(2.*m+1).*Y.^m;
%                     
%                 end
                
            elseif cs(obj.model_type, 'fourier')
                error('Fourier model is not yet implemented!');
            elseif cs(obj.model_type, 'interp')

                x = 1:im_size(2);
                y = 1:im_size(1);
            
                if obj.num_frames==1
                    obj.output_image = obj.interpolant({x,y});
                    obj.output_image = obj.output_image';
                else

                end

            else
                error('Unknown model_type "%s". Use "median" or "mean" or "polynomial" or "interp".', obj.model_type);
            end
            
            if nargout>0
                I = obj.output_image;
            end
            
        end
        
        function I = getVarImage(obj, im_size)
            
            im_size = util.vec.imsize(im_size);
            
            x = 1:im_size(2);
            y = 1:im_size(1);
            
            if obj.num_frames==1
                obj.output_var_image = obj.var_interp({x,y});
                obj.output_var_image = obj.output_var_image';
            else
                
            end
            
            if nargout>0
                I = obj.output_var_image;
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

