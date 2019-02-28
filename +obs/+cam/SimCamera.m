classdef SimCamera < handle
   
    properties % objects
        
        star@head.Star;
        hndl = 0;
        
    end
    
    properties % outputs
        
        images;
        timestamps;
        
    end
    
    properties % switches and inputs
        
        mode = '';
        
        name = 'SimCam';
        
        use_bad_pixels = 1;
        
        slow_mode = 0;
        
        product_type = 'RAW';
        
        debug_bit = 1;
        
    end
    
    properties(Hidden=false) 
        
        cooler_state = 0;
        temperature = 25;
        read_noise = 1.5;
        bias_level = 102;
        
        T;
        f;

        start_time = clock;
        
        width = 2048;
        height = 2048;
        top = 1;
        left = 1;
        
        star_vx;
        star_vy;
        
        bad_pixels;
        
        gain = 1;
        
        version = 1.00;
        
    end
    
    properties(Dependent=true)
        
    end
    
    methods % constructor
        
        function obj = SimCamera
           
            if obj.debug_bit, fprintf('SimCamera constructor v%4.2f\n', obj.version); end
             
            obj.initialize;
                       
        end
        
        function initialize(obj)
            
            obj.resetBadPixels;
            obj.setupStar;
            
        end
        
        function setupStar(obj)
           
            if isempty(obj.star), obj.star = head.Star; end
            
%             if isempty(obj.star.x), obj.star.x = rand.*(obj.width-1)+1; end
%             if isempty(obj.star.y), obj.star.y = rand.*(obj.height-1)+1; end
            
            if isempty(obj.star.anchor_x), obj.star.anchor_x = rand.*100-50 + obj.width/2 + obj.left; end
            if isempty(obj.star.anchor_y), obj.star.anchor_y = rand.*100-50 + obj.height/2 + obj.top; end
                        
            if isempty(obj.star_vx), obj.star_vx = 10*rand; end
            if isempty(obj.star_vy), obj.star_vy = 10*rand; end
            
            if isempty(obj.star.mag), obj.star.mag = 7; end
            
            obj.star.im_size = [obj.height, obj.width];
            
        end
        
    end
    
    methods % reset methods
        
        function clearOutputs(obj)
           
           obj.images = [];
           obj.timestamps = [];
           
        end
        
        function resetStar(obj)
           
            obj.star = head.Star;
            obj.star.mag = 5;
            
            obj.star_vx = normrnd(0, 1);
            obj.star_vy = normrnd(0, 1);
            
        end
        
        function resetBadPixels(obj)
           
            obj.bad_pixels = zeros(obj.maxHeight, obj.maxWidth);
            
            obj.bad_pixels(800, 900) = 300;
            obj.bad_pixels(830, 890) = 500;
            obj.bad_pixels(1250, 870) = 1000;
            obj.bad_pixels(870, 1000) = 6000;
            obj.bad_pixels(1100, 1200) = 1200;
            obj.bad_pixels(1000, 850) = 1200;
            obj.bad_pixels(1025, 900) = 300;
            obj.bad_pixels(750, 990) = 500;
            obj.bad_pixels(900, 1250) = 1100;
            obj.bad_pixels(780, 1020) = 3000;
            
        end
        
    end
    
    methods % getters
        
        function val = getTimestamp(obj)
            
            val = etime(clock, obj.start_time);
            
        end
        
        function val = getTemperature(obj)
           
            val = obj.temperature;
            
        end     
        
        function val = getCoolerState(obj)
           
            val = obj.cooler_state;
            
        end
        
        function val = getExpTime(obj)
            
            val = obj.T;
            
        end
        
        function val = getExpTimeMax(obj)
            
            val = 30;
            
        end
        
        function val = getExpTimeMin(obj)
            
            val = 0.00001;
            
        end
        
        function val = getFrameRate(obj)
           
            val = obj.f;
            
        end 
        
        function val = getFrameRateMax(obj)
           
            val = 1000;
            
        end  
        
        function val = getFrameRateMin(obj)
           
            val = 1/30;
            
        end
        
        function val = getHeight(obj)
           
            val = obj.height;
            
        end
        
        function val = getWidth(obj)
           
            val = obj.width;
            
        end
                
        function val = getLeft(obj)
           
            val = obj.left;
            
        end
        
        function val = getTop(obj)
           
            val = obj.top;
            
        end
        
        function val = maxHeight(obj)
           
            val = 2048;
            
        end
        
        function val = maxWidth(obj)
           
            val = 2048;
            
        end
        
    end
    
    methods % setters
        
        function setTemperature(obj, temp)
           
            obj.temperature = temp;
            
        end
        
        function setCoolerState(obj, state)
           
            obj.cooler_state = state;
            
        end
        
        function setExpTime(obj, time)
            
            obj.T = time;
            
        end
        
        function setFrameRate(obj, rate)
            
            obj.f = rate;
            
        end
          
        function setWidth(obj, val)
            
            obj.width = val;
            
        end
        
        function setHeight(obj, val)
            
            obj.height = val;
            
        end
        
        function setTop(obj, val)
            
            obj.top = val;
            
        end
        
        function setLeft(obj, val)
            
            obj.left = val;
            
        end
                        
        function setCenterX(obj, val)
            
            width = obj.getWidth;
            left = val - floor(width/2);
            
            assert(left>=1, ['cannot set center_x to ' num2str(val) ' the width ' num2str(width) ' is too big!']);
            
            obj.setLeft(left);
            
        end
        
        function setCenterY(obj, val)
            
            height = obj.getHeight;
            top = val - floor(height/2);
            
            assert(top>=1, ['cannot set center_y to ' num2str(val) ' the height ' num2str(height) ' is too big!']);
            
            obj.setTop(top);
            
        end
                
        function zoom(obj, width, height, center_x, center_y)
           
            if nargin<2 || isempty(width)
                width = 512;
            end
            
            if nargin<3 || isempty(height)
                height = width;
            end
            
            if nargin<4 || isempty(center_x)
                center_x = floor((obj.maxWidth)/2)+1;
            end
            
            if nargin<5 || isempty(center_y)
                center_y = floor((obj.maxHeight)/2)+1;
            end
            
            obj.setWidth(width);
            obj.setHeight(height);
            obj.setCenterX(center_x);
            obj.setCenterY(center_y);
        
        end
            
        function full_frame(obj)
           
            obj.width = obj.maxWidth;
            obj.height = obj.maxHeight;
            obj.left = 1;
            obj.top = 1;
            
        end
        
    end
    
    methods % actions
        
        function startup(obj, varargin)
            
            if obj.debug_bit>2, disp('starting up exposure with SimCamera'); end
            
            obj.start_time = clock;
            
%             obj.setupStar;
            
        end
        
        function finishup(obj, varargin)
            
            if obj.debug_bit, disp('finishing exposure with SimCamera'); end
            
        end
        
        function batch(obj, buffer, top_level_obj)
            
            import util.text.cs;
            
            if nargin<3 || isempty(top_level_obj)
                top_level_obj = [];
            end
            
            dark_mode = cs(obj.mode, 'dark');
            flat_mode = cs(obj.mode, 'flat');
            
            obj.clearOutputs;
            
            Nframes = top_level_obj.batch_size;
            
            temp_images = zeros(top_level_obj.height, top_level_obj.width, Nframes);
            temp_timestamps = zeros(Nframes, 1);
            
            for ii = 1:Nframes

                noise = zeros(obj.height, obj.width);

                if obj.use_bad_pixels
                    p = obj.bad_pixels(obj.top:obj.height+obj.top-1, obj.left:obj.width+obj.left-1);
                    noise = noise + p;
                end
                
                noise = normrnd(noise, obj.read_noise)+obj.bias_level;
%                 noise = rand(size(noise))*obj.read_noise+obj.bias_level;

                if ~dark_mode && ~flat_mode
                    
                    psf = util.img.gaussian2(3);
            
                    psf = poissrnd(psf);
                    
                    x = obj.star.final_x-obj.left;
                    y = obj.star.final_y-obj.top;
                
                    temp_images(:,:,ii) = obj.add_tile(noise, psf, [y, x]);
                                
                    border = 10;
                
                    if x<=border || x>=obj.width-border
                        obj.star_vx = -obj.star_vx;
                    elseif y<=border || y>=obj.height-border
                        obj.star_vy = -obj.star_vy;
                    end
                
                    obj.star.drift_x = obj.star.drift_x + obj.star_vx;
                    obj.star.drift_y = obj.star.drift_y + obj.star_vy;
                    
                elseif dark_mode
                    temp_images(:,:,ii) = noise;
                elseif flat_mode
                    error('flat mode is not implemented yet...');
                end
                
                temp_timestamps(ii) = etime(clock, obj.start_time);
                
                if ~isempty(top_level_obj)
                    top_level_obj.gui.update;
%                     top_level_obj.live_counter = top_level_obj.live_counter + 1;
%                     top_level_obj.show(temp_images(:,:,ii)); % also calls draw now after putting the image in place...
                   
                    if top_level_obj.brake_bit || ~top_level_obj.gui.check
                        break;
                    end
                    
                end
                
                if obj.slow_mode
                    drawnow;
                end
                
            end
            
            temp_images(:,:,ii+1:end) = [];
            temp_timestamps(ii+1:end) = [];
            
            obj.images = uint16(floor(temp_images));
            obj.timestamps = temp_timestamps;
            
        end
        
        function M_out = add_tile(obj, M_in, tile, position_yx)
            
            x1 = round(position_yx(2))-floor(size(tile,2)/2);
            x2 = round(position_yx(2))+ceil(size(tile,2)/2)-1;
            y1 = round(position_yx(1))-floor(size(tile,1)/2);
            y2 = round(position_yx(1))+ceil(size(tile,1)/2)-1;
            
            if x1<1 % if we are too far to the left
                tile = tile(:,2-x1:end);
                x1 = 1;
            end
            
            if x2>size(M_in, 2) % if we are too far to the right
                tile = tile(:, 1:end-(x2-size(M_in,2)));
                x2 = size(M_in,2);
            end
            
            if y1<1 % if we are too far to the top
                tile = tile(2-y1:end,:);
                y1 = 1;
            end
            
            if y2>size(M_in, 1) % if we are too far to the bottom
                tile = tile(1:end-(y2-size(M_in,1)),:);
                y2 = size(M_in,1);
            end
            
            M_out = M_in;
            M_out(y1:y2, x1:x2) = M_out(y1:y2, x1:x2) + tile;
            
        end
        
    end
    
end