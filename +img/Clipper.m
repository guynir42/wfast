classdef Clipper < handle
% Object used for cutting out stars, keeping track of the positions of
% stars, and for finding and adjusting the positions. 
% This object has four main responsibilities:
% (1) Cut stamps of a few pixels around each star (as many as requested by
%     the user). This is done using util.img.cutStars (or util.img.mexCutout)
%     that produces a 4D matrix with the cutouts themselves. 
%
% (2) Keep track of the cutout positions and size for many stars, so that
%     we do not waste time finding them for each image/batch, and so we can
%     keep track of where the stars were in the image (for matching
%     calibration data or for reconstructing their sky-positions). 
%
% (3) Find the positions of N stars in an image, by locating the brightest
%     points that are (hopefully) not bad pixels or cosmic rays. This is
%     done once at the beginning of the run, and again if the stars are
%     somehow displaced outside the cutouts or just not observed at all for
%     a while. 
%
% (4) Measure changes in cutout positions when there is telescope drift of 
%     the field or refraction drift of individual stars. This is done at 
%     the end of each batch. If there are drifts, they can be used to
%     adjust the cutout positions or as feedback to the telescope control. 
% 
% (5) Realign the position of the cutouts relative to an input image and 
%     find the new positions in case the stars went out of the cutouts due 
%     to telescope drift or jump. 

    properties(Transient=true)
        
        gui@img.gui.ClipperGUI;
        
    end
    
    properties % outputs
        
        cutouts; % a 4D matrix of the image cutouts
        mean_drift; % in x and y (the same as pos), dim 1 is the number of batches
        
    end
   
    properties % switches and controls
        
        positions; % dim 1 is all the different cutouts, dim 2 is x and then y
        cut_size = 32; % size of all cutouts (square, uniform for all stars)
        
        start_pos; % starting positions for each star (same as positions)
        start_stack; % the stack image from which we found the stars (for realign).
        start_cuts;  % the cutouts from when we found the stars (for realign).
        
        use_adjust = 0; % we shouldn't be using internal adjustments anymore... 
        use_lock_adjust = 1; % force adjustment of all cutouts together (e.g., telescope drift)
        use_mex = 1; % use util.img.mexCutout (this is about x10 faster)
        use_moments = 1;
        pad_value = 0; % when clipping outside the edges of the frame...
        use_padding_warning = 0;
        use_find_stars_for_realign = 0; % if you want to refind all stars when realign is called. (we need to cancel this option...)
        
        % for use in findStars:
        num_stars = 1; % how many stars to find?
        star_threshold = 7; % in units of noise STD?
        use_masking = 1;
        use_filtering = 1;
        filter_size = 10;
        filter_sigma = 3;
        avoid_edges; % how many pixels away from the edge you want to scan when using findStars
        
        subframe_size; % cut a subframe from the cutouts for adjustCuts only
        subframe_pos; % cut a specific position (default is center) from the cutouts for adjustCuts only
        
        number_cuts_display = 4;
        
        rem_value = 0; % = NaN; % what to put when removing stars... to be depricated
        
        debug_bit = 1;        
        
    end
    
    properties % timing data
       
        cut_runtime = 0;
        adjust_runtime = 0;
        find_runtime = 0;
        remove_runtime = 0;
        
    end
    
    properties(Hidden=true)
        
        default_cut_size;
        default_filter_size;
        default_filter_sigma;
        
        default_number_cuts_display;
        
        version = 1.06;
        
    end
    
    properties(Dependent=true)
       
        N;
        
    end
    
    methods % constructor
        
        function obj = Clipper(other)
           
            if nargin>0 && ~isempty(other) && isa(other, 'img.Clipper')
               
                if other.debug_bit, fprintf('img.Clipper copy-constructor v%4.2f\n', obj.version); end
                
                obj = util.oop.full_copy(other);
                                
            else
                
                if obj.debug_bit, fprintf('img.Clipper constructor v%4.2f\n', obj.version); end
                
                util.oop.save_defaults(obj);
                obj.avoid_edges = obj.cut_size;
                
            end
            
        end
        
    end
    
    methods % reset methods
       
        function startup(obj, val)
           
            obj.start_pos = val;
            obj.positions = val;
            obj.resetOutputs;
             
        end
        
        function reset(obj)
           
            obj.resetPositions;
            obj.resetOutputs;
            obj.resetRuntime;
            obj.clear;
            
        end
        
        function resetPositions(obj)
           
            obj.positions = []; 
            obj.start_pos = [];
            obj.start_stack = [];
            obj.start_cuts = [];
            
        end
        
        function resetOutputs(obj)
            
            obj.mean_drift = [];
            
        end
        
        function resetRuntime(obj)
            
            obj.cut_runtime = 0;
            obj.find_runtime = 0;
            obj.adjust_runtime = 0;
            obj.remove_runtime = 0;

        end
        
        function clear(obj)
            
            obj.cutouts = [];
            
        end
        
    end
    
    methods % getters
        
        function check = is_empty(obj)
            
            check = isempty(obj.positions);
            
        end
        
        function check = is_equal(obj, other)
           
            if isempty(other)
                check = 0;
            elseif isa(other, 'img.Clipper')
            
                check = obj.N==other.N && ...
                    all(obj.positions(:)==other.positions(:)) &&...
                    obj.cut_size == other.cut_size;
                
            elseif isnumeric(other) && size(other, 2)==2
                check = obj.N==size(other,1) && all(obj.positions(:)==other(:));
            end
            
        end
        
        function val = get.N(obj)
           
            val = size(obj.positions,1);
            
        end
                
        function ind = lower_corner(obj, center_ind, im_size)
            
            if nargin<2 || isempty(center_ind)
                center_ind = obj.positions;
            end
            
            if nargin<3 || isempty(im_size)
                im_size = obj.cut_size;
            end
            
            if isempty(center_ind)
                ind = [];
            else
                ind = center_ind - floor(im_size/2);
            end
            
        end
        
        function ind = upper_corner(obj, center_ind, im_size)
            
            if nargin<2 || isempty(center_ind)
                center_ind = obj.positions;
            end
            
            if nargin<3 || isempty(im_size)
                im_size = obj.cut_size;
            end
            
            ind = obj.lower_corner(center_ind, im_size) + im_size - 1;
            
        end
        
        function ind = upper_corner_old(obj, center_ind, im_size)
        
            if nargin<2 || isempty(center_ind)
                center_ind = obj.positions;
            end
            
            if nargin<3 || isempty(im_size)
                im_size = obj.cut_size;
            end
            
            if isempty(center_ind)
                ind = [];
            else
                ind = center_ind + ceil((im_size)/2);
            end
            
        end
                
    end
    
    methods % setters
        
        function set.positions(obj, val)
            
            if ~isempty(val) && size(val,2)~=2
                error(['size of input to positions should be a Nx2 matrix. size(val)= ' num2str(size(val))]);
            end
            
            if isempty(obj.start_pos)
                obj.start_pos = val;
            end
            if all(size(obj.positions)==size(val)) && all(obj.positions(:)==val(:))
                return;
            end
            
            obj.positions = val;
                        
%             obj.subframe_pos = [];
            
        end
        
        function set.cut_size(obj, val)
            
            if obj.cut_size==val
                return;
            end
            
            obj.cut_size = val;
            obj.subframe_pos = [];
            
        end
        
        function set.number_cuts_display(obj, val)
            
            if obj.number_cuts_display==val
                return;
            end
            
            obj.number_cuts_display = val;
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.gui.makeAxes;
            end
            
        end
        
        function chooseSubframe(obj)
            
            if isempty(obj.gui) || ~obj.gui.check
                return;
            end
            
            [X,Y] = getpts(obj.gui.axes_plot{1});
            
            obj.subframe_pos = [X,Y];
            
            obj.show;
            
        end
        
    end
    
    methods % actions
        
        function M_out = cutMatrix(obj, M, positions, cut_size, pad_value)
            
            % two modes, regular or mex file
                        
            if isempty(M)
                M_out = [];
                return;
            end
            
            if nargin<3 || isempty(positions)
                positions = obj.positions;
            end
            
            if nargin<4 || isempty(cut_size)
                cut_size = obj.cut_size;
            end
            
            if nargin<5 || isempty(pad_value)
                pad_value = obj.pad_value;
            end
            
            if nnz(isnan(positions))
                error('Some (or all) the positions are NaN!');
            end
            
            cut_t = tic;
            
            if obj.use_mex==0
                M_out = zeros(cut_size, 'like', M);
                M_out = repmat(M_out, [1 1 size(M,3) size(positions, 1)]);
                N = size(M,3);

                x_max = size(M,2);
                y_max = size(M,1);

                for ii = 1:size(positions, 1)
                                        
                    l_x = round(obj.lower_corner(positions(ii,1), cut_size));
                    u_x = round(obj.upper_corner(positions(ii,1), cut_size));
                    l_y = round(obj.lower_corner(positions(ii,2), cut_size));
                    u_y = round(obj.upper_corner(positions(ii,2), cut_size));

                    x1 = max([1, l_x]);
                    x2 = min([x_max, u_x]);

                    y1 = max([1, l_y]);
                    y2 = min([y_max, u_y]);
                    
                    if l_x>x_max || u_x<1 || l_y>y_max || u_y<1
                        M_temp = zeros(cut_size, cut_size, N);
                    else

                        M_temp = M(y1:y2, x1:x2, :);

                        if l_x<1
                            if obj.use_padding_warning, warning(['adding left padding with ' num2str(1-l_x) ' columns']); end
                            M_temp = horzcat(ones(size(M_temp,1), 1-l_x, N)*pad_value, M_temp); % left padding
                        end

                        if u_x>x_max
                            if obj.use_padding_warning, warning(['adding right padding with ' num2str(u_x-x_max) ' columns']); end
                            M_temp = horzcat(M_temp, ones(size(M_temp,1), u_x-x_max, N)*pad_value); % right padding
                        end

                        if l_y<1
                            if obj.use_padding_warning, warning(['adding top padding with ' num2str(1-l_y) ' rows']); end
                            M_temp = vertcat(ones(1-l_y, size(M_temp,2), N)*pad_value, M_temp); % top padding
                        end

                        if u_y>y_max
                            if obj.use_padding_warning, warning(['adding bottom padding with ' num2str(u_y-y_max) ' rows']); end
                            M_temp = vertcat(M_temp, ones(u_y-y_max, size(M_temp,2), N)*pad_value); % bottom padding
                        end

                    end
                    
                    M_out(:,:,:,ii) = M_temp;

                end
                
            else
                M_out = util.img.mexCutout(M, positions, cut_size, pad_value);
            end
            
            obj.cut_runtime = obj.cut_runtime + toc(cut_t);
            
        end
        
        function M = input(obj, images, varargin) 
            
            import util.text.cs;
            import util.text.parse_bool;
            
            adjust = [];
            
            for ii = 1:2:length(varargin)
                key = varargin{ii};
                val = varargin{ii+1};
                
                if cs(key, 'adjust', 'need adjust')
                    adjust = parse_bool(val);
                end
                
            end
            
            if isempty(adjust)
                adjust = obj.use_adjust;
            end
            
            if isempty(images)
                M = [];
                return;
            end
            
            if ~isempty(obj.cut_size) && obj.cut_size==size(images,1) && obj.cut_size==size(images,2) % if the cut is the same size as the images!
                obj.cutouts = images;
            else % otherwise go on to cut
            
                if isempty(obj.positions) % if you need to find stars

                    if ~isempty(obj.start_pos)
                        obj.positions = obj.start_pos;
                    else
                        if obj.debug_bit, disp('no starting positions given, using findStars on input image'); end

                        if obj.findStars(sum(images, 3))
                            adjust = 0;
                            if obj.debug_bit, disp(['successfully found stars at ' num2str(obj.positions(1,:))]); end
                        else
                            error('failed to find stars!');
                        end
                    end

                end
            
                obj.cutouts = obj.cutMatrix(images);
                
                if adjust
                    obj.adjustCuts(obj.cutouts);
                end
            
            end
            
            if nargin>0
                M = obj.cutouts;
            end
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.show;
            end
            
        end
        
        function shift = adjustCuts(obj, cutouts) % to be depricated!
            
            error('Please do not use this function any more! Instead, adjust the positions from external class...');
            
            import util.stat.sum2;
            import util.stat.max2;
            import util.img.gaussian2;
            
            csize = obj.cut_size; % cut size
            if isscalar(csize)
                csize = [1 1]*csize;
            end
            
            if ~isempty(obj.subframe_size) % using a smaller "sub frame"
                
                S = obj.subframe_size;
                if isscalar(S)
                    S = [1 1]*S;
                end
                
                csize = S;
                
                if isempty(obj.subframe_pos)
                    obj.subframe_pos = floor([size(cutouts,2) size(cutouts,1)]/2)+1; % x then y!
                end
                
                if size(obj.subframe_pos,1)==size(cutouts,4)
                    pos = obj.subframe_pos;
                elseif size(obj.subframe_pos,1)==1
                    pos = repmat(obj.subframe_pos, [size(cutouts,4) 1]);
                else
                    error(['size mismatch: size(subframe_pos)= ' num2str(size(obj.subframe_pos)) ' | size(image_sum_cut)= ' num2str(size(cutouts))]);
                end
                
                new_matrix = zeros(S(1), S(2), size(cutouts,3), size(cutouts,4));
                
                for ii = 1:size(cutouts, 4)
                    
                    new_matrix(:,:,:,ii) = obj.cutMatrix(cutouts(:,:,:,ii), pos(ii,:), S);
                    
                end
                
                cutouts = new_matrix;
                
            end
            
            adjust_t = tic;
            
            all_drifts = [];
            weight = ones(obj.N,1);
            
            cutouts = sum(cutouts,3); % verify we are adjusting based on the summed cutouts
            
            if obj.use_moments
                
                M1x= zeros(size(obj.positions,1), 1);
                M1y= zeros(size(obj.positions,1), 1);
                
%                 [x,y] = meshgrid((1:csize(2))-csize(2)/2, (1:csize(1))-csize(1)/2);
                [x,y] = util.vec.centerGrid(zeros(csize));
                
                for ii = 1:obj.N % go over all stars/cutouts
                    
                    I = abs(cutouts(:,:,:,ii)); % this should be a 2D matrix
                    I = util.img.maskBadPixels(I);
                    I = I - util.stat.corner_median(I);
                    
%                     I = conv2(I, gaussian2(obj.filter_sigma, [], [], obj.filter_size), 'same'); % smoothing filter
                    weight(ii) = sum2(I);
                    
                    % calculate the first moment of the smoothed image
                    M1x(ii) = sum2(I.*x)/weight(ii);
                    M1y(ii) = sum2(I.*y)/weight(ii);
                    
                end
                
                all_drifts = [M1x, M1y];
            
            else
                
                mxy = zeros(size(obj.positions,1),2);                
                weight = zeros(size(obj.positions,1),1);
                
                for ii = 1:obj.N % go over all stars/cutouts
                    
                    I = cutouts(:,:,:,ii); % this should be a 2D matrix
                    I = util.img.maskBadPixels(I);
                    I = conv2(I, gaussian2(obj.filter_sigma, [], [], obj.filter_size), 'same'); % smoothing filter

                    [mval, ind] = max2(I);
                    
                    mxy(ii,:) = (flip(ind) - csize/2);
                    weight(ii) = mval;
                    
                end
                
                
%                 disp(['mxy= ' num2str(mxy) ' | weight= ' num2str(weight)]);
                
                all_drifts = mxy;
                
            end
            
            if mean(weight)>0, weight = weight./mean(weight); end
            
            shift = median(weight.*all_drifts, 1, 'omitnan');
            
            if obj.use_lock_adjust % move all cuts together by the weighted average drift
                obj.positions = obj.positions + round(shift);
            else
                obj.positions = obj.positions - round(all_drifts);
            end
            
            obj.adjust_runtime = obj.adjust_runtime + toc(adjust_t);
            
            obj.mean_drift(end+1, :) = shift;
            
        end
        
        function shift = realignCuts(obj, stack) % can work on calibrated or non calibrated images...

%             disp('realigning stars!');
            
            if obj.use_find_stars_for_realign
                obj.findStars(stack);
            else
            
                [~,shift,confidence] = util.img.quick_align(obj.start_stack, stack, 'static',1, 'flip', 1);
            
                obj.positions = obj.start_pos - shift;
                
            end
            
        end
        
        function [I_points, idx] = makePointImage(obj, im_size)
            
            im_size = util.vec.imsize(im_size);
            
            x = obj.positions(:,2);
            y = obj.positions(:,1);
            
            % remove points outside the image size
            y(x<1 | x>im_size(2)) = [];
            x(x<1 | x>im_size(2)) = [];
            x(y<1 | y>im_size(1)) = [];
            y(y<1 | y>im_size(1)) = [];
            
            idx = sub2ind(im_size, y, x);
            
            I_points = zeros(im_size);
            
            I_points(idx) = 1; 
            
        end
        
        function [I_squares, x, y] = makeSquaresImage(obj, im_size)
           
            im_size = util.vec.imsize(im_size);
            
            x0 = obj.positions(:,2);
            y0 = obj.positions(:,1);
            
            spread = -floor(obj.cut_size/2):floor((obj.cut_size-1)/2);
            spread_x = repmat(spread, [1, obj.cut_size]);
            spread_y = reshape(repmat(spread, [obj.cut_size,1]), [1, obj.cut_size.^2]);
            
            x = x0 + spread_x;
            y = y0 + spread_y;
            
            idx = sub2ind(im_size, x, y);
            
            I_squares = zeros(im_size);
            
            I_squares(idx) = 1; 
            
        end
        
        function I_cutouts = makeCutoutsImage(obj, im_size, cutouts, positions)
            
            if nargin<3 || isempty(cutouts)
                cutouts = obj.cutouts;
            end
            
            if nargin<4 || isempty(positions)
                positions = obj.positions;
            end
            
            if isempty(cutouts)
                error('Cannot do "makeCutoutsImage" without cutouts!');
            end
            
            im_size = util.vec.imsize(im_size);
            
            I_cutouts = zeros();
            
            C = sum(cutouts,3);

            for ii = 1:size(C,4)
                
                c0 = floor(obj.cut_size/2)+1;
                x0 = positions(ii,1);
                x_low = x0 - floor(obj.cut_size/2);
                x_high = x0 + floor((obj.cut_size-1)/2);
                if x_low>im_size(2) || x_high<1, continue; end
                if x_low<1, x_low = 1; end
                if x_high>im_size(2), x_high = im_size(2); end
                
                y0 = positions(ii,2);
                y_low = y0 - floor(obj.cut_size/2);
                y_high = y0 + floor((obj.cut_size-1)/2);
                if y_low>im_size(1) || y_high<1, continue; end
                if y_low<1, y_low = 1; end
                if y_high>im_size(1), y_high = im_size(1); end
                
                I_cutouts(y_low:y_high,x_low:x_high) = C(y_low-y0+c0:y_high-y0+c0,x_low-x0+c0:x_high-x0+c0,1,ii);
                
            end
            
        end
        
        function success = findStars(obj, images_full) % can work on calibrated or non calibrated images...
            
            import util.img.maskBadPixels;
            import util.img.gaussian2;
            import util.stat.max2;
            import util.stat.var2;
                        
            success = 0; % will only be considered successfull if one star passes the "star_threshold" 
            
            obj.reset; % is this the best place to put this...?
            
            I = images_full;
                                    
            find_t = tic;
            
            if obj.use_masking
                I = maskBadPixels(I); % in case we don't get a clean image... 
            end
            
%             I(isnan(I)) = util.stat.median2(I);
            
            if obj.use_filtering
                I(isnan(I)) = 0;
                I = conv2(I, gaussian2(obj.filter_sigma, 'size', obj.filter_size), 'same'); % smoothing filter
            end
            
            noise_std = sqrt(var2(I));
            
            height = size(I,1);
            width = size(I,2);
            
            if obj.avoid_edges
                y1 = obj.avoid_edges;
                y2 = height-obj.avoid_edges;
                x1 = obj.avoid_edges;
                x2 = width-obj.avoid_edges;
            end
            
            for ii = 1:obj.num_stars
                
                if obj.avoid_edges
                    [mx, ind] = max2(I, y1:y2, x1:x2);
                else % use the whole frame... 
                    [mx, ind] = max2(I);
                end
                
                if mx>=noise_std*obj.star_threshold
                    success = 1;
                end
                
                ind = fliplr(ind); % this turns the y then x index into x then y...
            
                if mx>obj.star_threshold
                   
                    obj.start_pos(ii,:) = ind; % this stores the positions found into the clip object
                    obj.positions(ii,:) = ind;
                    
                    low = obj.lower_corner(ind);
                    up = obj.upper_corner(ind);
                    
                    if low(2)<1
                        warning(['star number ' num2str(ii) ' out of bounds. corners: ' num2str([low up])]);
                        low(2) = 1;
                    end
                    
                    if low(1)<1
                        warning(['star number ' num2str(ii) ' out of bounds. corners: ' num2str([low up])]);
                        low(1) = 1;
                    end
                    
                    if up(2)>size(I,1)
                        warning(['star number ' num2str(ii) ' out of bounds. corners: ' num2str([low up])]);
                        up(2) = height;
                    end
                    
                    if up(1)>size(I,2)
                        warning(['star number ' num2str(ii) ' out of bounds. corners: ' num2str([low up])]);
                        up(1) = width;
                    end
                    
                    I(low(2):up(2), low(1):up(1)) = 0; % get rid of this star to look for the others...
                    
                end
                
            end
            
            obj.start_stack = images_full;
            obj.start_cuts = obj.cutouts;
            
            obj.find_runtime = obj.find_runtime + toc(find_t);
            
        end
        
        function M_out = removeStars(obj, M)
            
            if isempty(obj.positions) % if you need to find stars
                
                if ~isempty(obj.start_pos)
                    obj.positions = obj.start_pos;
                else
                    if obj.debug_bit, disp('no starting positions given, using findStars on input image'); end
                    
                    if obj.findStars(sum(M, 3))
                        if obj.debug_bit, disp(['successfully found stars at ' num2str(obj.positions(1,:))]); end
                    else
                        error('failed to find stars!');
                    end
                end
                
            end
            
            remove_t = tic;
            
            M_out = M;
            
            x_max = size(M,2);
            y_max = size(M,1);
            
            for ii = 1:size(obj.positions, 1)
                
                l_x = round(obj.lower_corner(obj.positions(ii,1), obj.cut_size));
                u_x = round(obj.upper_corner(obj.positions(ii,1), obj.cut_size));
                l_y = round(obj.lower_corner(obj.positions(ii,2), obj.cut_size));
                u_y = round(obj.upper_corner(obj.positions(ii,2), obj.cut_size));
                
                x1 = max([1, l_x]);
                x2 = min([x_max, u_x]);
                
                y1 = max([1, l_y]);
                y2 = min([y_max, u_y]);
                
                M_out(y1:y2,x1:x2,:) = obj.rem_value;
                
            end
            
            obj.remove_runtime = obj.remove_runtime + toc(remove_t);
         
        end
        
        function arbitraryPositions(obj, varargin)
            
            import util.text.cs;
            
            input = util.text.InputVars;
            input.input_var('num_cuts', obj.num_stars, 'number_cutouts');
            input.input_var('cut_size', obj.cut_size, 'size');
            input.input_var('image_size', 2048, 'im_size');
            input.input_var('clipper', []);
            input.input_var('check_mode', 'overlap'); % can also choose "radius"
            input.input_var('choose_mode', 'random'); % can also choose "grid"
            input.scan_vars(varargin{:});
            
            obj.reset;
            
            S = util.vec.imsize(input.image_size);
            
            obj.cut_size = input.cut_size;
            
            if ~isempty(input.clipper) && ~isa(input.clipper, 'img.Clipper')
                error('Must input an "img.Clipper" object to "clipper" field. Instead got "%s"', class(input.clipper));
            end
            
            if ~isempty(input.clipper) && (isempty(input.clipper.positions) || isempty(input.clipper.cut_size))
                input.clipper = img.Clipper.empty; 
            end
            
            if cs(input.choose_mode, 'grid')
                
                num_on_side = ceil(sqrt(input.num_cuts)); % how many points on each side of the square
                
                dx = S(2)/num_on_side;
                dy = S(1)/num_on_side;
                
                x = round(-dx/2);
                y = round(dy/2);
                
            end
            
            for ii = 1:input.num_cuts
                
                for jj = 1:100 % try several locations 
                    
                    check = 1;
                    if cs(input.choose_mode, 'random')
                        x = randi([ceil(input.cut_size/2) S(2)-ceil(input.cut_size/2)]);
                        y = randi([ceil(input.cut_size/2) S(1)-ceil(input.cut_size/2)]);
                    elseif cs(input.choose_mode, 'grid')
                        x = round(x + dx);
                        if x>S(2)
                            x = round(dx/2);
                            y = round(y + dy);
                        end
                    end
                    
                    if ~isempty(input.clipper) % if another clipper is given, check that there's no overlap
                        check = ~input.clipper.checkOverlap([x,y],input.cut_size, input.check_mode);
                    end
                    
                    check = ~obj.checkOverlap([x,y],input.cut_size, input.check_mode); % check this coordinate has no overlap with this object
                    
                    if check
                        obj.positions(end+1,:) = [x y]; 
                        break;
                    end
                    
                end
                
            end
            
            obj.start_pos = obj.positions;
            obj.num_stars = size(obj.positions,1);
            
        end
        
        function result = checkOverlap(obj, position_xy, cut_size, mode)
            
            if nargin<4 || isempty(mode)
                mode = 'overlap';
            end
            
            if isnumeric(mode)
                % pass
            elseif util.text.cs(mode, 'overlap')
                mode = 1;
            elseif util.text.cs(mode, 'radius')
                mode = 2;
            else
                error('No such mode "%s". Use "overlap" or "radius"', mode);
            end
            
            X = position_xy(1); % X coordinate of the test position
            Y = position_xy(2); % Y coordinate of the test position
            S = cut_size;       % size of the test position
            s = obj.cut_size;   % size of this object's cutouts
            R = ceil((S+s)/2);  % distance needed between any of these positions
                
            result = 0;
            
            for ii = 1:size(obj.positions, 1)
                
                x = obj.positions(ii,1);
                y = obj.positions(ii,2);
                
                if mode==1 % overlap (square)
                    if abs(x-X)<R && abs(y-Y)<R
                        result = 1;
                        return;
                    end
                elseif mode==2 % radius (circle)
                    if sqrt((x-X).^2+(y-Y).^2)<R
                        result = 1;
                        return;
                    end
                end
                
            end
            
        end
        
        function makeGUI(obj)
           
            if isempty(obj.gui)
                obj.gui = img.gui.ClipperGUI(obj);
            end
            
            obj.gui.makeGUI;
            
        end
        
    end
    
    methods % plotting tools
        
        function show(obj, varargin)
            
            import util.plot.show;
            
            if isempty(obj.cutouts)
                return;
            end
            
            N = min(size(obj.cutouts, 4), obj.number_cuts_display);
            pos = obj.subframe_pos;
            
            if size(pos,1)==1 && size(pos,2)==2
                pos = repmat(pos, [N,1]);
            end
            
            for ii = 1:N
                
                idx = ceil(size(obj.cutouts,3)/2);
                show(obj.cutouts(:,:,idx,ii), 'ax', obj.gui.axes_plot{ii}, 'autodyn', 'on', 'fancy', 'off');
                
                if ~isempty(obj.subframe_size) && ~isempty(obj.subframe_pos)
                    
                    delete(findobj(obj.gui.axes_plot{ii}, 'type', 'rectangle'));
                    
                    rectangle('Position', [obj.lower_corner(pos(ii,:), obj.subframe_size), obj.subframe_size*[1 1]], 'Parent', obj.gui.axes_plot{ii}, 'EdgeColor', 'Green');
                    
                end
                
            end
            
            drawnow;
            
        end
        
        function showRectangles(obj, varargin)
            
            import util.text.cs;
            import util.text.parse_bool;
            
            if isempty(obj.positions)
                return;
            end
            
            ax = [];
            flip = [];
            num = [];
            color = [];
            del = 1;
            use_text = 1;
            
            for ii = 1:2:length(varargin)
               if cs(varargin{ii}, {'axes', 'axis'})
                   ax = varargin{ii+1};
               elseif cs(varargin{ii}, 'flip')
                   flip = parse_bool(varargin{ii+1});
               elseif cs(varargin{ii}, 'number')
                   num = varargin{ii+1};
               elseif cs(varargin{ii}, 'color')
                   color = varargin{ii+1};
               elseif cs(varargin{ii}, 'delete')
                   del = parse_bool(varargin{ii+1});
               elseif cs(varargin{ii}, 'text', 'use_text')
                   use_text = parse_bool(varargin{ii+1});
               end
            end

            if isempty(ax)
                ax = gca;
            end
            
            h = findobj(ax, 'Type', 'Image');
            S = size(h.CData);
            
            C = double(obj.positions);
            if flip
                C = fliplr(S) - C;
            end
            
            if del
                delete(findobj(ax, 'type', 'rectangle'));
                if use_text, delete(findobj(ax, 'type', 'text')); end
            end
            
            if ~isempty(num) && num<obj.N
                if use_text 
                    text(0.05*S(2),0.95*S(1), sprintf('showing %d/%d cutouts!', num, obj.N), 'FontSize', 16, 'Parent', ax);
                end
            else 
                num = obj.N;
            end
            
            for ii = 1:num
                try
                    
                    if use_text, text(C(ii,1), C(ii,2), ['clip ' num2str(ii)],'FontSize', 16, 'Parent', ax); end
                    rectangle('Position', [obj.lower_corner(C(ii,:))-0.5 obj.cut_size obj.cut_size], 'Parent', ax, 'EdgeColor', color);
                    
                catch ME
                    warning(ME.getReport);
                end
            end
            
        end
        
    end
    
end